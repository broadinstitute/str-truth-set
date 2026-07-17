"""This pipeline runs ATaRVa (https://github.com/SowpatiLab/ATaRVa) on long read data.

ATaRVa genotypes tandem repeats from long reads. It takes a BAM/CRAM, a reference fasta, and a bgzipped + tabix-indexed
regions bed (chrom, start_0based, end_1based, motif, motif_length) and writes a single-sample VCF whose FORMAT column
(GT:AL:CN:LPM:AR:SD:DP:SN:SQ:MA:MR:DS:MV) gives, per allele, the allele length in base pairs (AL) and the motif copy
number (CN).

This pipeline feeds ATaRVa the truth set's single-motif loci bed (the same {sample}.bed.gz catalog that inquiSTR uses),
runs `atarva genotype` in whole-genome mode, then converts the VCF to the ExpansionHunter .json format and combines it
into the .variants.tsv.gz / .alleles.tsv.gz tables used by the tool comparison scripts.
"""

import hailtop.fs as hfs
import os
import re

from step_pipeline import pipeline, Backend, Localize

DOCKER_IMAGE = "weisburd/atarva@sha256:39f19cc6d8fc5c195b9b167c580f6d6a65771889a22ff37c7f2ad779ed997cb5"

# str-analysis commit baked-in check: the converter (convert_atarva_vcf_to_expansion_hunter_json) is installed at
# runtime from this pinned commit so the converter output stays reproducible across runs (matches the vamos pipeline).
STR_ANALYSIS_COMMIT = "3d5e3dc37161d41a1bc6f92afdcf3fc99e81b30a"

REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/atarva"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", required=True, help="Path of the long-read BAM/CRAM to genotype")
    parser.add_argument("--input-bai", required=True, help="Path of the .bai/.crai index for --input-bam")
    parser.add_argument("--regions-bed", required=True,
                        help="Path of the bgzipped + tabix-indexed regions bed (chrom, start, end, motif, motif_length)")
    parser.add_argument("--male-or-female", default="female", choices=["male", "female"])
    parser.add_argument("--amplicon", action="store_true", help="Run ATaRVa in targeted (amplicon) mode")
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("--cpu", type=int, default=4)
    args = bp.parse_known_args()

    regions_bed_paths = [x.path for x in hfs.ls(args.regions_bed)]
    if len(regions_bed_paths) == 0:
        raise ValueError(f"No files found matching {args.regions_bed}")

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: ATaRVa: {bam_path_ending}")

    create_atarva_step(
        bp,
        reference_fasta=args.reference_fasta,
        reference_fasta_fai=args.reference_fasta_fai,
        input_bam=args.input_bam,
        input_bai=args.input_bai,
        regions_bed_path=args.regions_bed,
        output_dir=args.output_dir,
        output_prefix=re.sub("(.bam|.cram)$", "", os.path.basename(args.input_bam)),
        male_or_female=args.male_or_female,
        amplicon=args.amplicon,
        cpu=args.cpu)
    bp.run()


def create_atarva_step(bp, *, reference_fasta, input_bam, input_bai, regions_bed_path, output_dir, output_prefix,
                       reference_fasta_fai=None, male_or_female="female", amplicon=False, cpu=4, catalog_step=None):
    # catalog_step: optional upstream Step that produces the input catalog bed; the genotyping step depends on it so it
    # never runs before the catalog exists.
    # Size the disk to the localized bam/cram plus 30Gi of headroom (reference fasta + vcf/json outputs), rather than a
    # fixed 200Gi -- ATaRVa reads the alignment in place and writes only a small vcf, so it never needs much beyond the
    # input file itself.
    hfs_ls_results = hfs.ls(input_bam)
    if len(hfs_ls_results) == 0:
        raise ValueError(f"No files found matching {input_bam}")
    storage_gb = int(hfs_ls_results[0].size / 10**9) + 30

    s1 = bp.new_step(f"Run ATaRVa on {os.path.basename(input_bam)}  ({os.path.basename(regions_bed_path)})",
                     arg_suffix="run-atarva-step",
                     step_number=1,
                     image=DOCKER_IMAGE,
                     cpu=cpu,
                     # ATaRVa fans out one process per thread (-t {cpu}), each with its own read/SNV buffers, so memory
                     # scales ~linearly with cpu. Per the paper it peaks near ~1GB/thread (PacBio HiFi/ONT Duplex) and
                     # ~3GB/thread (ONT Simplex), so cpu=4 fits comfortably in standard memory (4GB/core).
                     memory="standard",
                     storage=f"{storage_gb}Gi",
                     output_dir=output_dir)
    if catalog_step is not None:
        s1.depends_on(catalog_step)
    s1.command("set -ex")

    # install str-analysis at a pinned commit so the ATaRVa->ExpansionHunter json converter output stays reproducible
    # across runs (matches the vamos pipeline). Ideally bake a pinned str-analysis into the ATaRVa image and drop this.
    s1.command("python3 -m pip install --no-cache-dir "
               f"'git+https://github.com/broadinstitute/str-analysis@{STR_ANALYSIS_COMMIT}'")

    local_fasta = s1.input(reference_fasta, localize_by=Localize.COPY)
    if reference_fasta_fai:
        s1.input(reference_fasta_fai, localize_by=Localize.COPY)
    else:
        s1.input(f"{reference_fasta}.fai", localize_by=Localize.COPY)

    local_bam = s1.input(input_bam, localize_by=Localize.COPY)
    if input_bai:
        s1.input(input_bai, localize_by=Localize.COPY)

    # ATaRVa reads the regions from a bgzipped + tabix-indexed bed, so localize the .tbi alongside the bed.
    local_regions_bed = s1.input(regions_bed_path, localize_by=Localize.COPY)
    s1.input(f"{regions_bed_path}.tbi", localize_by=Localize.COPY)

    # The bw2 ATaRVa fork passes the reference to pysam when opening the alignment file, so CRAM input decodes reliably
    # and is genotyped directly (no BAM conversion). --format tells ATaRVa which reader to use.
    atarva_format = "cram" if input_bam.endswith(".cram") else "bam"

    s1.command("df -kh")

    # ATaRVa needs the sample's karyotype to genotype the sex chromosomes as haploid in male samples (chrX/chrY outside
    # the PAR). Without it every locus is genotyped as diploid.
    karyotype = "XY" if male_or_female == "male" else "XX"
    amplicon_arg = "--amplicon " if amplicon else ""

    # Lower --min-reads from ATaRVa's default of 10 to 2, matching LongTR's setting in this benchmark, so the low-depth
    # comparison is coverage-matched. At the default 10, ATaRVa no-calls most loci at low coverage (8x) and on male
    # haploid chrX/chrY (half depth), which counts as discordant and unfairly sinks its "% green" versus tools run at a
    # lower threshold -- even though it genotypes those loci at 95-100% accuracy when it does call them.
    s1.command(f"/usr/bin/time --verbose atarva genotype "
               f"-f {local_fasta} "
               f"-b {local_bam} "
               f"-r {local_regions_bed} "
               f"--format {atarva_format} "
               f"--karyotype {karyotype} "
               f"--min-reads 2 "
               f"-t {cpu} "
               f"{amplicon_arg}"
               f"-o {output_prefix}.vcf")
    s1.command(f"bgzip {output_prefix}.vcf")
    s1.command("ls -lhrt")

    s1.command(f"python3 -m str_analysis.convert_atarva_vcf_to_expansion_hunter_json --discard-hom-ref "
               f"--sample-id {output_prefix.split('.')[0]} "
               f"{output_prefix}.vcf.gz")

    s1.command("ls -lhrt")
    s1.output(f"{output_prefix}.vcf.gz")
    s1.output(f"{output_prefix}.json")

    # step2: combine the json into the variants.tsv.gz / alleles.tsv.gz tables
    s2 = bp.new_step(name=f"Combine ATaRVa outputs for {os.path.basename(input_bam)}",
                     step_number=2,
                     arg_suffix="combine-atarva-step",
                     image=DOCKER_IMAGE,
                     cpu=2,
                     memory="highmem",
                     storage="20Gi",
                     output_dir=output_dir)

    s2.depends_on(s1)

    s2.command("python3 -m pip install --no-cache-dir "
               f"'git+https://github.com/broadinstitute/str-analysis@{STR_ANALYSIS_COMMIT}'")

    s2.command("set -ex")
    s2.command("mkdir /io/run_dir; cd /io/run_dir")
    local_json = s2.input(os.path.join(output_dir, f"{output_prefix}.json"), localize_by=Localize.COPY)
    s2.command(f"ln -s {local_json}")

    s2.command(f"python3 -m str_analysis.combine_str_json_to_tsv --output-prefix {output_prefix}")

    s2.command(f"mv {output_prefix}.1_json_files.bed {output_prefix}.bed")
    s2.command(f"mv {output_prefix}.1_json_files.variants.tsv.gz {output_prefix}.variants.tsv.gz")
    s2.command(f"mv {output_prefix}.1_json_files.alleles.tsv.gz {output_prefix}.alleles.tsv.gz")

    s2.command(f"bgzip {output_prefix}.bed")
    s2.command(f"tabix {output_prefix}.bed.gz")
    s2.command("ls -lhrt")
    s2.output(f"{output_prefix}.variants.tsv.gz")
    s2.output(f"{output_prefix}.alleles.tsv.gz")
    s2.output(f"{output_prefix}.bed.gz")
    s2.output(f"{output_prefix}.bed.gz.tbi")

    return s2


if __name__ == "__main__":
    main()
