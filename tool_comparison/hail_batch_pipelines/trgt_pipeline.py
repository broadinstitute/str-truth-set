"""This pipeline runs TRGT (https://github.com/PacificBiosciences/trgt)

Usage:
    trgt [OPTIONS] --genome <FASTA> --repeats <REPEATS> --vcf <VCF> --spanning-reads <SPANNING_READS> --repeat-id <REPEAT_ID> --image <IMAGE>

Options:
      --genome <FASTA>
          Path to reference genome FASTA
      --repeats <REPEATS>
          BED file with repeat coordinates
      --vcf <VCF>
          VCF file generated by TRGT
      --spanning-reads <SPANNING_READS>
          BAM file with spanning reads generated by TRGT
      --repeat-id <REPEAT_ID>
          ID of the repeat to plot
      --image <IMAGE>
          Output image path
      --plot-type <PLOT_TYPE>
          Type of plot to generate [default: allele] [possible values: allele, waterfall]
      --show <SHOW>
          What to show in the plot [default: motifs] [possible values: motifs, meth]
      --flank-len <FLANK_LEN>
          Length of flanking regions [default: 50]
      -v, --verbose...
"""

import hail as hl
import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/trgt@sha256:031ed186025494b715fdf3decd2906b0fff82d9ce99f52f788ca8d89689e47dc"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

TRGT_CATALOG_BED_POSITIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/trgt/positive_loci.TRGT_repeat_catalog.bed"
TRGT_CATALOG_BED_NEGATIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/trgt/negative_loci.TRGT_repeat_catalog.bed"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/trgt"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser_group = parser.add_mutually_exclusive_group(required=True)
    parser_group.add_argument("--positive-loci", action="store_true", help="Genotype truth set loci")
    parser_group.add_argument("--negative-loci", action="store_true", help="Genotype negative (hom-ref) loci")
    parser_group.add_argument("--trgt-catalog-bed", action="append", help="Path of TRGT catalog bed file(s) to process")

    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_CRAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    args = bp.parse_known_args()

    if args.positive_loci:
        positive_or_negative_loci = "positive_loci"
        trgt_catalog_bed_file_stats_list = hl.hadoop_ls(TRGT_CATALOG_BED_POSITIVE_LOCI)
        if len(trgt_catalog_bed_file_stats_list) == 0:
            raise ValueError(f"No files found matching {TRGT_CATALOG_BED_POSITIVE_LOCI}")
    elif args.negative_loci:
        positive_or_negative_loci = "negative_loci"
        trgt_catalog_bed_file_stats_list = hl.hadoop_ls(TRGT_CATALOG_BED_NEGATIVE_LOCI)
        if len(trgt_catalog_bed_file_stats_list) == 0:
            raise ValueError(f"No files found matching {TRGT_CATALOG_BED_NEGATIVE_LOCI}")
    elif args.trgt_catalog_bed:
        positive_or_negative_loci = os.path.basename(args.trgt_catalog_bed[0]).replace(".bed", "").replace(".gz", "")
        trgt_catalog_bed_file_stats_list = [{"path": path} for path in args.trgt_catalog_bed]
    else:
        parser.error("Must specify either --positive-loci or --negative-loci")

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: TRGT: {positive_or_negative_loci}: {bam_path_ending}")
    output_dir = os.path.join(args.output_dir, positive_or_negative_loci)
    if not args.force:
        log_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.log.gz"))
        logging.info(f"Precached {len(log_paths)} log files")
        vcf_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.vcf.gz"))
        logging.info(f"Precached {len(vcf_paths)} vcf files")
        bam_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.bam"))
        logging.info(f"Precached {len(bam_paths)} log files")

    step1s = []
    for trgt_catalog_i, trgt_catalog_bed_file_stats in enumerate(trgt_catalog_bed_file_stats_list):
        trgt_catalog_bed_path = trgt_catalog_bed_file_stats["path"]

        if args.n and trgt_catalog_i >= args.n:
            break

        cpu = 16
        s1 = bp.new_step(f"Run TRGT #{trgt_catalog_i}",
                         arg_suffix=f"trgt",
                         step_number=1,
                         image=DOCKER_IMAGE,
                         cpu=cpu,
                         storage="200Gi",
                         output_dir=output_dir)
        step1s.append(s1)

        s1.command("set -ex")
        local_fasta = s1.input(args.reference_fasta, localize_by=Localize.COPY)
        if args.reference_fasta_fai:
            s1.input(args.reference_fasta_fai, localize_by=Localize.COPY)

        local_bam = s1.input(args.input_bam, localize_by=Localize.COPY)
        if args.input_bai:
            s1.input(args.input_bai, localize_by=Localize.COPY)

        local_trgt_catalog_bed = s1.input(trgt_catalog_bed_path)
        s1.command("df -kh")
        output_prefix = re.sub(".bed(.gz)?$", "", local_trgt_catalog_bed.filename)
        s1.command(f"echo Genotyping $(cat {local_trgt_catalog_bed} | wc -l) loci")
        s1.command(f"""/usr/bin/time --verbose trgt \
                                     --genome {local_fasta} \
                                     --reads {local_bam} \
                                     --repeats {local_trgt_catalog_bed} \
                                     --output-prefix {output_prefix} \
                                     --threads {cpu} \
                                     --verbose
        """)

        s1.command("ls -lhrt")

        #s1.command(f"python3.9 -m str_analysis.convert_hipstr_vcf_to_expansion_hunter_json {output_prefix}.vcf.gz")
        s1.output(f"{output_prefix}.vcf.gz")
        s1.output(f"{output_prefix}.spanning.bam")

    bp.run()
    return

    # step2: combine json files
    s2 = bp.new_step(name="Combine TRGT outputs",
                     step_number=2,
                     image=DOCKER_IMAGE,
                     cpu=2,
                     memory="highmem",
                     storage="20Gi",
                     output_dir=output_dir)
    for step1 in step1s:
        s2.depends_on(step1)

    s2.command("mkdir /io/run_dir; cd /io/run_dir")
    for json_path in step1_output_json_paths:
        local_path = s2.input(json_path, localize_by=Localize.COPY)
        s2.command(f"ln -s {local_path}")

    output_prefix = f"combined.{positive_or_negative_loci}"
    s2.command("set -x")
    s2.command(f"python3.9 -m str_analysis.combine_str_json_to_tsv "
               f"--output-prefix {output_prefix}")
    s2.command(f"bgzip {output_prefix}.{len(step1_output_json_paths)}_json_files.bed")
    s2.command(f"tabix {output_prefix}.{len(step1_output_json_paths)}_json_files.bed.gz")
    s2.command("gzip *.tsv")
    s2.command("ls -lhrt")
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.variants.tsv.gz")
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.alleles.tsv.gz")
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.bed.gz")
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.bed.gz.tbi")
    bp.run()


if __name__ == "__main__":
    main()


