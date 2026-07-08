"""This pipeline runs inquiSTR (https://github.com/wdecoster/inquiSTR) on long read data.

inquiSTR genotypes tandem repeats from long reads by measuring, per haplotype, the median length difference from the
reference genome (in base pairs) across all supporting reads. The `inquiSTR call` subcommand takes a BAM/CRAM and a
region BED file (columns: chromosome, begin, end, info) and writes a tab-separated table to stdout with columns:

    chromosome  begin  end  info  <sample>_H1  <sample>_H2

This pipeline generates the region BED from the truth set positive_loci.bed.gz catalog, putting the locus id
("{chrom}-{start_0based}-{end}-{motif}") in the info column so that the downstream converter can recover the repeat
motif and reference coordinates. The inquiSTR output is then converted to the ExpansionHunter .json format and combined
into the .variants.tsv.gz / .alleles.tsv.gz tables used by the tool comparison scripts.
"""

import hailtop.fs as hfs
import os
import re

from step_pipeline import pipeline, Backend, Localize

# Bakes in str-analysis pinned in docker_with_inquistr/Dockerfile (includes the inquiSTR -> ExpansionHunter json
# converter), so no str-analysis is installed at runtime. Rebuild that image and update this digest when str-analysis changes.
DOCKER_IMAGE = "weisburd/str-analysis-with-inquistr@sha256:7aa72cfb8388b1e4e57c75e75a07fd23601d4687bdb2699b439e59dafe55cc7b"

REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

CHM1_CHM13_BAM_PATH = "gs://bw2-delete-after-60-days/long-reads/CHM1_CHM13.subreads.bam"
CHM1_CHM13_BAI_PATH = "gs://bw2-delete-after-60-days/long-reads/CHM1_CHM13.subreads.bam.bai"

POSITIVE_LOCI_BED = "gs://str-truth-set/hg38/positive_loci.bed.gz"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/inquistr"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_BAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_BAI_PATH)
    parser.add_argument("--positive-loci-bed", default=POSITIVE_LOCI_BED,
                        help="Path of the truth set positive_loci.bed.gz catalog (columns: chrom, start0, end, motif)")
    parser.add_argument("--male-or-female", default="female", choices=["male", "female"])
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("--cpu", type=int, default=16)
    args = bp.parse_known_args()

    positive_loci_bed_paths = [x.path for x in hfs.ls(args.positive_loci_bed)]
    if len(positive_loci_bed_paths) == 0:
        raise ValueError(f"No files found matching {args.positive_loci_bed}")

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: inquiSTR: {bam_path_ending}")

    output_prefix = re.sub(".bed(.gz)?$", "", os.path.basename(positive_loci_bed_paths[0]))
    create_inquistr_steps(
        bp,
        reference_fasta=args.reference_fasta,
        reference_fasta_fai=args.reference_fasta_fai,
        input_bam=args.input_bam,
        input_bai=args.input_bai,
        inquistr_catalog_bed_paths=positive_loci_bed_paths,
        output_dir=args.output_dir,
        output_prefix=output_prefix,
        male_or_female=args.male_or_female,
        cpu=args.cpu)
    bp.run()


def create_inquistr_steps(bp, *, reference_fasta, input_bam, input_bai, inquistr_catalog_bed_paths, output_dir,
                          output_prefix, reference_fasta_fai=None, male_or_female="female", cpu=16, unphased=True,
                          catalog_step=None):
    # catalog_step: optional upstream Step that produces the input catalog; the genotyping step(s) depend on it so they
    # never run before the catalog exists.
    if len(inquistr_catalog_bed_paths) > 1:
        raise ValueError("Only one inquiSTR catalog bed file is currently supported")
    if len(inquistr_catalog_bed_paths) == 0:
        raise ValueError("No inquiSTR catalog bed file provided")
    inquistr_catalog_bed_path = inquistr_catalog_bed_paths[0]

    s1 = bp.new_step(f"Run inquiSTR on {os.path.basename(input_bam)}  {os.path.basename(inquistr_catalog_bed_path)}",
                     arg_suffix="run-inquistr-step",
                     step_number=1,
                     image=DOCKER_IMAGE,
                     cpu=cpu,
                     storage="200Gi",
                     output_dir=output_dir)
    if catalog_step is not None:
        s1.depends_on(catalog_step)
    s1.command("set -ex")

    local_fasta = s1.input(reference_fasta, localize_by=Localize.COPY)
    if reference_fasta_fai:
        s1.input(reference_fasta_fai, localize_by=Localize.COPY)
    else:
        s1.input(f"{reference_fasta}.fai", localize_by=Localize.COPY)
    local_bam = s1.input(input_bam, localize_by=Localize.COPY)
    if input_bai:
        s1.input(input_bai, localize_by=Localize.COPY)
    local_catalog = s1.input(inquistr_catalog_bed_path)

    s1.command("df -kh")

    # build the inquiSTR region bed (chrom, start0, end, locus_id) from the positive_loci.bed.gz catalog
    # (chrom, start0, end, motif, ...). The locus id "{chrom}-{start0}-{end}-{motif}" is stored in the info column so
    # the converter can recover the motif and reference coordinates from the inquiSTR output.
    decompress = "gunzip -c" if str(local_catalog).endswith(".gz") else "cat"
    inquistr_region_bed = f"{output_prefix}.inquistr_regions.bed"
    s1.command(f"""{decompress} {local_catalog} | awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, $1"-"$2"-"$3"-"$4}}' """
               f"> {inquistr_region_bed}")
    s1.command(f"echo Genotyping $(wc -l < {inquistr_region_bed}) loci in {local_bam.filename}")

    haploid_arg = "--haploid chrX,chrY" if male_or_female == "male" else ""
    # The truth-set bams are not haplotagged (no HP tags), so inquiSTR is run in --unphased mode, which splits the
    # reads at each locus into two size-based alleles instead of using read haplotype assignments.
    unphased_arg = "--unphased" if unphased else ""
    s1.command(f"""/usr/bin/time --verbose inquiSTR call {local_bam} \
                                     -R {inquistr_region_bed} \
                                     --reference {local_fasta} \
                                     --threads {cpu} \
                                     {unphased_arg} \
                                     {haploid_arg} \
                                     > {output_prefix}.inq
    """)
    s1.command("ls -lhrt")

    s1.command(f"python3 -m str_analysis.convert_inquistr_calls_to_expansion_hunter_json --discard-hom-ref "
               f"{output_prefix}.inq")

    s1.command("ls -lhrt")
    s1.output(f"{output_prefix}.inq")
    s1.output(f"{output_prefix}.json")

    # step2: combine the json into the variants.tsv.gz / alleles.tsv.gz tables
    s2 = bp.new_step(name=f"Combine inquiSTR outputs for {os.path.basename(input_bam)}",
                     step_number=2,
                     arg_suffix="combine-inquistr-step",
                     image=DOCKER_IMAGE,
                     cpu=2,
                     memory="highmem",
                     storage="20Gi",
                     output_dir=output_dir)

    s2.depends_on(s1)

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
