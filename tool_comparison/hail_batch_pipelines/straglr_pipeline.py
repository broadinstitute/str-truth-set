"""This pipeline runs Straglr (https://github.com/bcgsc/straglr)

Usage:
    python straglr.py <mm2.bam> <reference_fasta> <output_prefix> [--loci loci.bed] [--exclude skip_regions.bed] [--chroms chr] [--regions regions.bed] [--min_support N] [--min_ins_size N] [--min_str_len N] [--max_str_len N] [--nprocs N] [--genotype_in_size] [--max_num_clusters N] [--min_cluster_size N] [--working_dir] [--tmpdir] [--debug]

Options:
    --loci: a BED file containing loci to be genotyped. 4 column BED format: chromosome start end repeat
    --exclude: a BED file containing regions to be skipped in genome-scan (e.g. long segmental duplications or pericentromeric regions)
    --chroms: space-separated list of specific chromosomes for genome-scan
    --regions: a BED file containing regions to be used only in genome-scan
    --include_alt_chroms: include ALT chromosomes (chromosomes with "_" in names) in genome scan (Default: NOT included)
    --use_unpaired_clips: include examination of unpaired clipped alignments in genome scan to detect expansion beyond read size (Default:NOT used)
    --min_support: minimum number of suppport reads for an expansion to be captured in genome-scan (Default:2)
    --min_ins_size: minimum increase in size (relative to the reference genome) for an expansion to be captured in genome-scan (Default:100)
    --min_str_len: minimum length of repeat-motif for an expansion to be captured in genome-scan (Default:2)
    --max_str_len: maximum length of repeat-motif for an expansion to be captured in genome-scan (Default:50)
    --nprocs: number of processes to use in Python's multiprocessing (Default:1)
    --genotype_in_size: report genotype (column 5 of TSV output) in terms of allele sizes instead of copy numbers
    --max_num_clusters: maximum number of clusters to be tried in Gausssian Mixture Model (GMM) clustering (Default:2)
    --min_cluster_size: minimum number of reads required to constitute a cluster (allele) in GMM clustering (Default:2)
    --trf_args: TRF arguments (Default:2 5 5 80 10 10 500)
    --tmpdir: user-specified directory for holding temporary files
"""

import hailtop.fs as hfs
import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/straglr@sha256:b95c51f44b55655fad6d48a1b694ed2f8115b2633b528f09673057bccd2c1719"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

STRAGLR_CATALOG_BED_POSITIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/straglr/positive_loci.straglr_catalog.bed"
STRAGLR_CATALOG_BED_NEGATIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/straglr/negative_loci.straglr_catalog.bed"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/straglr"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser_group = parser.add_mutually_exclusive_group(required=True)
    parser_group.add_argument("--positive-loci", action="store_true", help="Genotype truth set loci")
    parser_group.add_argument("--negative-loci", action="store_true", help="Genotype negative (hom-ref) loci")
    parser_group.add_argument("--straglr-catalog-bed", action="append", help="Path of straglr catalog bed file(s) to process")

    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_CRAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    args = bp.parse_known_args()

    if args.positive_loci:
        positive_or_negative_loci = "positive_loci"
        straglr_catalog_bed_paths = [x.path for x in hfs.ls(STRAGLR_CATALOG_BED_POSITIVE_LOCI)]
        if len(straglr_catalog_bed_paths) == 0:
            raise ValueError(f"No files found matching {STRAGLR_CATALOG_BED_POSITIVE_LOCI}")
    elif args.negative_loci:
        positive_or_negative_loci = "negative_loci"
        straglr_catalog_bed_paths = [x.path for x in hfs.ls(STRAGLR_CATALOG_BED_NEGATIVE_LOCI)]
        if len(straglr_catalog_bed_paths) == 0:
            raise ValueError(f"No files found matching {STRAGLR_CATALOG_BED_NEGATIVE_LOCI}")
    elif args.straglr_catalog_bed:
        positive_or_negative_loci = os.path.basename(args.straglr_catalog_bed[0]).replace(".bed", "").replace(".gz", "")
        straglr_catalog_bed_paths = [x.path for x in hfs.ls(args.straglr_catalog_bed)]
    else:
        parser.error("Must specify either --positive-loci or --negative-loci")

    if args.n:
        straglr_catalog_bed_paths = straglr_catalog_bed_paths[:args.n]

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: straglr: {positive_or_negative_loci}: {bam_path_ending}")
    output_dir = os.path.join(args.output_dir, positive_or_negative_loci)
    if not args.force:
        #tsv_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.tsv.gz"))
        #logging.info(f"Precached {len(tsv_paths)} tsv files")
        vcf_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.bed.gz*"))
        logging.info(f"Precached {len(vcf_paths)} bed and tbi files")

    create_straglr_steps(
        bp,
        reference_fasta=args.reference_fasta,
        input_bam=args.input_bam,
        input_bai=args.input_bai,
        straglr_catalog_bed_paths=straglr_catalog_bed_paths,
        output_dir=output_dir,
        output_prefix=f"combined.{positive_or_negative_loci}",
        reference_fasta_fai=args.reference_fasta_fai)
    bp.run()


def create_straglr_steps(bp, *, reference_fasta, input_bam, input_bai, straglr_catalog_bed_paths, output_dir, output_prefix, reference_fasta_fai=None):
    step1s = []
    step1_output_paths = []
    hfs_ls_results = hfs.ls(input_bam)
    if len(hfs_ls_results) == 0:
        raise ValueError(f"No files found matching {input_bam}")
    input_bam_file_stats = hfs_ls_results[0]

    for straglr_catalog_i, straglr_catalog_bed_path in enumerate(straglr_catalog_bed_paths):
        s1 = bp.new_step(f"Run straglr #{straglr_catalog_i}",
                         arg_suffix=f"straglr",
                         step_number=1,
                         image=DOCKER_IMAGE,
                         cpu=16,
                         localize_by=Localize.COPY,
                         storage=f"{int(input_bam_file_stats.size/10**9) + 25}Gi",
                         output_dir=output_dir)
        step1s.append(s1)

        local_fasta = s1.input(reference_fasta)
        if reference_fasta_fai:
            s1.input(reference_fasta_fai)

        local_bam = s1.input(input_bam)
        if input_bai:
            s1.input(input_bai)

        local_straglr_catalog_bed = s1.input(straglr_catalog_bed_path)

        output_prefix = re.sub(".bed(.gz)?$", "", local_straglr_catalog_bed.filename)
        s1.command(f"echo Genotyping $(cat {local_straglr_catalog_bed} | wc -l) loci")
        s1.command("set -ex")
        s1.command(f"""/usr/bin/time --verbose python3 /usr/local/bin/straglr.py {local_bam} {local_fasta} {output_prefix} \
                --loci {local_straglr_catalog_bed} \
                --nprocs {cpu//2}""")

        s1.command("ls -lhrt")

        #s1.command(f"gzip {output_prefix}.tsv")
        s1.command(f"python3.9 -m str_analysis.convert_straglr_bed_to_expansion_hunter_json {output_prefix}.bed")
        s1.command("ls -lhtr")
        s1.command(f"gzip {output_prefix}.json")
        s1.command(f"bgzip {output_prefix}.bed")
        s1.command(f"tabix {output_prefix}.bed.gz")
        s1.output(f"{output_prefix}.bed.gz")
        s1.output(f"{output_prefix}.bed.gz.tbi")
        s1.output(f"{output_prefix}.json.gz")
        step1_output_paths.append(os.path.join(output_dir, f"{output_prefix}.json.gz"))

    # step2: combine json files
    s2 = bp.new_step(name=f"Combine straglr outputs for {os.path.basename(input_bam)}",
                     step_number=2,
                     image=DOCKER_IMAGE,
                     cpu=2,
                     memory="highmem",
                     storage="20Gi",
                     output_dir=output_dir)
    for step1 in step1s:
        s2.depends_on(step1)
    s2.command("mkdir /io/run_dir; cd /io/run_dir")
    for json_path in step1_output_paths:
        local_path = s2.input(json_path, localize_by=Localize.COPY)
        s2.command(f"ln -s {local_path}")

    s2.command("set -x")
    s2.command(f"python3.9 -m str_analysis.combine_str_json_to_tsv "
               f"--output-prefix {output_prefix}")
    s2.command(f"mv {output_prefix}.{len(step1_output_json_paths)}_json_files.bed {output_prefix}.bed")
    s2.command(f"mv {output_prefix}.{len(step1_output_json_paths)}_json_files.variants.tsv.gz {output_prefix}.variants.tsv.gz")
    s2.command(f"mv {output_prefix}.{len(step1_output_json_paths)}_json_files.alleles.tsv.gz {output_prefix}.alleles.tsv.gz")

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


