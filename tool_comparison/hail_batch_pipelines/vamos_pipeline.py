"""This pipeline runs Vamos (https://github.com/ChaissonLab/vamos) on long read data.

vamos annotates the motif composition of tandem repeat loci from long reads. It takes a BAM/CRAM and a vamos catalog
(a headerless TSV: chrom, start, end, motifs, version, STR_or_VNTR, ...) and writes a single-sample diploid VCF whose
INFO fields give, per haplotype, the motif annotation (ALTANNO_Hx) and total motif count (LEN_Hx).

This pipeline derives the vamos catalog from the truth set's single-motif ExpansionHunter catalog (so each locus has
one motif), runs vamos in its default 1-based mode, then converts the VCF to the ExpansionHunter .json format and
combines it into the .variants.tsv.gz / .alleles.tsv.gz tables used by the tool comparison scripts.
"""

import hailtop.fs as hfs
import os
import re

from step_pipeline import pipeline, Backend, Localize

DOCKER_IMAGE = "weisburd/vamos@sha256:05bac98f8da9587ca00207d10eb619ee51a4864c3716ed9de3739d92a315b629"

REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/vamos"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_CRAI_PATH)
    parser.add_argument("--expansion-hunter-catalog", required=True,
                        help="Path of the ExpansionHunter catalog JSON to derive the vamos catalog from")
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("--cpu", type=int, default=16)
    args = bp.parse_known_args()

    expansion_hunter_catalog_paths = [x.path for x in hfs.ls(args.expansion_hunter_catalog)]
    if len(expansion_hunter_catalog_paths) == 0:
        raise ValueError(f"No files found matching {args.expansion_hunter_catalog}")

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: VAMOS: {bam_path_ending}")

    output_prefix = re.sub(".json(.gz)?$", "", os.path.basename(expansion_hunter_catalog_paths[0]))
    create_vamos_step(
        bp,
        reference_fasta=args.reference_fasta,
        reference_fasta_fai=args.reference_fasta_fai,
        input_bam=args.input_bam,
        input_bai=args.input_bai,
        expansion_hunter_catalog_paths=expansion_hunter_catalog_paths,
        output_dir=args.output_dir,
        output_prefix=output_prefix,
        cpu=args.cpu)
    bp.run()


def create_vamos_step(bp, *, reference_fasta, input_bam, input_bai, expansion_hunter_catalog_paths, output_dir,
                      output_prefix, reference_fasta_fai=None, cpu=16):
    if len(expansion_hunter_catalog_paths) > 1:
        raise ValueError("Only one ExpansionHunter catalog file is currently supported")
    if len(expansion_hunter_catalog_paths) == 0:
        raise ValueError("No ExpansionHunter catalog file provided")
    expansion_hunter_catalog_path = expansion_hunter_catalog_paths[0]

    vamos_catalog = f"{output_prefix}.vamos_catalog.tsv"

    s1 = bp.new_step(f"Run Vamos on {os.path.basename(input_bam)}  {os.path.basename(expansion_hunter_catalog_path)}",
                     arg_suffix="run-vamos-step",
                     step_number=1,
                     image=DOCKER_IMAGE,
                     cpu=cpu,
                     storage="200Gi",
                     output_dir=output_dir)
    s1.command("set -ex")

    # the vamos image doesn't bundle str-analysis, so install it (with its dependencies) before running the catalog
    # and vcf converters. Pin to a commit (not @main) so the converter output stays reproducible across runs; ideally
    # bake a pinned str-analysis into the vamos image as docker_with_inquistr does and drop this runtime install.
    s1.command("python3 -m pip install --no-cache-dir "
               "'git+https://github.com/broadinstitute/str-analysis@62608ee5b68dec6e66b55a07b3e00389d7d0e31c'")

    local_bam = s1.input(input_bam, localize_by=Localize.COPY)
    if input_bai:
        s1.input(input_bai, localize_by=Localize.COPY)
    local_eh_catalog = s1.input(expansion_hunter_catalog_path)

    # vamos has no reference option and decodes alignments via htslib, so reference-compressed CRAM input must be
    # decoded against the reference first. For CRAM input, localize the reference FASTA and convert to BAM with
    # samtools (bundled in the vamos image); BAM input is used as-is.
    if input_bam.endswith(".cram"):
        local_fasta = s1.input(reference_fasta, localize_by=Localize.COPY)
        if reference_fasta_fai:
            s1.input(reference_fasta_fai, localize_by=Localize.COPY)
        else:
            s1.input(f"{reference_fasta}.fai", localize_by=Localize.COPY)
        s1.command(f"samtools view -@ {cpu} -b -T {local_fasta} -o {output_prefix}.input.bam {local_bam}")
        s1.command(f"samtools index -@ {cpu} {output_prefix}.input.bam")
        vamos_input_bam = f"{output_prefix}.input.bam"
    else:
        vamos_input_bam = local_bam

    s1.command("df -kh")

    # build the vamos catalog (chrom, start, end, motifs, ...) from the single-motif truth set ExpansionHunter catalog.
    # vamos aborts (exit 1) on unsorted or overlapping loci, so sort by position and greedily drop any locus that
    # overlaps the previously kept one. The catalog uses 1-based fully-closed start/end (start=start_0based+1,
    # end=end_1based), so two loci overlap when the next start is <= the previous end (prev_end >= $2, not > $2).
    s1.command(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_vamos_catalog "
               f"-o {output_prefix}.vamos_catalog.unsorted.tsv {local_eh_catalog}")
    s1.command(f"sort -k1,1 -k2,2n {output_prefix}.vamos_catalog.unsorted.tsv | "
               f"""awk 'BEGIN{{OFS="\\t"; prev_chrom=""; prev_end=0}} """
               f"""{{ if ($1==prev_chrom && prev_end>=$2) next; print; prev_chrom=$1; prev_end=$3 }}' """
               f"> {vamos_catalog}")
    s1.command(f"echo Genotyping $(wc -l < {vamos_catalog}) loci in {local_bam.filename}")

    # Example: vamos --read -b ../example/demo.aln.bam -r vamos.effMotifs-0.1.GRCh38.tsv -s NA24385_CCS_h1 -o reads.vcf -t 8
    s1.command(f"/usr/bin/time --verbose vamos --read -b {vamos_input_bam} -r {vamos_catalog} "
               f"-s {output_prefix} -o {output_prefix}.vcf -t {cpu}")
    s1.command(f"bgzip {output_prefix}.vcf")
    s1.command("ls -lhrt")

    s1.command(f"python3 -m str_analysis.convert_vamos_vcf_to_expansion_hunter_json --discard-hom-ref "
               f"{output_prefix}.vcf.gz")

    s1.command("ls -lhrt")
    s1.output(f"{vamos_catalog}")
    s1.output(f"{output_prefix}.vcf.gz")
    s1.output(f"{output_prefix}.json")

    # step2: combine the json into the variants.tsv.gz / alleles.tsv.gz tables
    s2 = bp.new_step(name=f"Combine vamos outputs for {os.path.basename(input_bam)}",
                     step_number=2,
                     arg_suffix="combine-vamos-step",
                     image=DOCKER_IMAGE,
                     cpu=2,
                     memory="highmem",
                     storage="20Gi",
                     output_dir=output_dir)

    s2.depends_on(s1)

    s2.command("python3 -m pip install --no-cache-dir "
               "'git+https://github.com/broadinstitute/str-analysis@62608ee5b68dec6e66b55a07b3e00389d7d0e31c'")

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
