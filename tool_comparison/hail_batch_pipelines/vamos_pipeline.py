"""This pipeline runs Vamos (https://github.com/PacificBiosciences/vamos)"""

import hailtop.fs as hfs
import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/vamos@sha256:05bac98f8da9587ca00207d10eb619ee51a4864c3716ed9de3739d92a315b629"

REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

VAMOS_CATALOG_POSITIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/vamos/positive_loci.vamos_catalog.tsv"
VAMOS_CATALOG_NEGATIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/vamos/negative_loci.vamos_catalog.tsv"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/vamos"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser_group = parser.add_mutually_exclusive_group(required=True)
    parser_group.add_argument("--positive-loci", action="store_true", help="Genotype truth set loci")
    parser_group.add_argument("--negative-loci", action="store_true", help="Genotype negative (hom-ref) loci")
    parser_group.add_argument("--vamos-catalog", help="Path of Vamos catalog file(s) to process")

    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_CRAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("--cpu", type=int, default=16)
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    args = bp.parse_known_args()

    if args.positive_loci:
        positive_or_negative_loci = "positive_loci"
        vamos_catalog_paths = [x.path for x in hfs.ls(VAMOS_CATALOG_POSITIVE_LOCI)]
        if len(vamos_catalog_paths) == 0:
            raise ValueError(f"No files found matching {VAMOS_CATALOG_POSITIVE_LOCI}")
    elif args.negative_loci:
        positive_or_negative_loci = "negative_loci"
        vamos_catalog_paths = [x.path for x in hfs.ls(VAMOS_CATALOG_NEGATIVE_LOCI)]
        if len(vamos_catalog_paths) == 0:
            raise ValueError(f"No files found matching {VAMOS_CATALOG_NEGATIVE_LOCI}")
    elif args.vamos_catalog:
        positive_or_negative_loci = os.path.basename(args.vamos_catalog).replace(".tsv", "").replace(".gz", "")
        vamos_catalog_paths = [x.path for x in hfs.ls(args.vamos_catalog)]
    else:
        parser.error("Must specify either --positive-loci or --negative-loci")

    if args.n:
        vamos_catalog_paths = vamos_catalog_paths[:args.n]

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: VAMOS: {positive_or_negative_loci}: {bam_path_ending}")
    output_dir = os.path.join(args.output_dir, positive_or_negative_loci)
    if not args.force:
        log_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.log.gz"))
        logging.info(f"Precached {len(log_paths)} log files")
        vcf_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.vcf.gz"))
        logging.info(f"Precached {len(vcf_paths)} vcf files")
        bam_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.bam"))
        logging.info(f"Precached {len(bam_paths)} log files")

    for vamos_catalog_i, vamos_catalog_path in enumerate(vamos_catalog_paths):
        output_prefix = re.sub(".tsv(.gz)?$", "", os.path.basename(vamos_catalog_path))
        s1 = create_vamos_step(
            bp,
            reference_fasta=args.reference_fasta,
            input_bam=args.input_bam,
            input_bai=args.input_bai,
            vamos_catalog_paths=vamos_catalog_paths,
            output_dir=output_dir,
            output_prefix=output_prefix,
            reference_fasta_fai=args.reference_fasta_fai,
            cpu=args.cpu)
    bp.run()


def create_vamos_step(bp, *, reference_fasta, input_bam, input_bai, vamos_catalog_paths, output_dir, output_prefix,
                     reference_fasta_fai=None, male_or_female="female", parse_reference_region_from_locus_id=False, cpu=16):
    if len(vamos_catalog_paths) > 1:
        raise ValueError("Can only process one Vamos catalog file at a time")
    if len(vamos_catalog_paths) == 0:
        raise ValueError("No Vamos catalog file provided")
    vamos_catalog_path = vamos_catalog_paths[0]

    s1_output_json_paths = []

    s1 = bp.new_step(f"Run Vamos on {os.path.basename(input_bam)}  {os.path.basename(vamos_catalog_path)}",
                     arg_suffix=f"run-vamos-step",
                     step_number=1,
                     image=DOCKER_IMAGE,
                     cpu=cpu,
                     storage="200Gi",
                     output_dir=output_dir)
    s1.command("set -ex")
    #local_fasta = s1.input(reference_fasta, localize_by=Localize.COPY)
    #if reference_fasta_fai:
    #    s1.input(reference_fasta_fai, localize_by=Localize.COPY)
    #else:
    #    s1.input(f"{reference_fasta}.fai", localize_by=Localize.COPY)
    local_bam = s1.input(input_bam, localize_by=Localize.COPY)
    if input_bai:
        s1.input(input_bai, localize_by=Localize.COPY)
    local_vamos_catalog = s1.input(vamos_catalog_path)
    s1.command("df -kh")
    s1.command(f"echo Genotyping $(cat {local_vamos_catalog} | wc -l) loci in {local_bam.filename}")
    # Example: vamos --read -b ../example/demo.aln.bam -r vamos.effMotifs-0.1.GRCh38.tsv -s NA24385_CCS_h1 -o reads.vcf -t 8

    s1.command(f"/usr/bin/time --verbose vamos --read -b {local_bam} -r {local_vamos_catalog} -o {output_prefix}.vcf -t {cpu}")
    s1.command(f"bgzip {output_prefix}.vcf")
    s1.command("ls -lhrt")

    #s1.command(f"python3 -m str_analysis.convert_vamos_vcf_to_expansion_hunter_json --discard-hom-ref {output_prefix}.vcf.gz")

    s1.output(f"{output_prefix}.vcf.gz")

    s1_output_json_paths.append(os.path.join(output_dir, f"{output_prefix}.vcf"))
    
    # step2: combine vcf files
    #s2 = bp.new_step(name=f"Combine TRGT outputs for {os.path.basename(input_bam)}", 
    #                 step_number=2,
    #                 arg_suffix=f"combine-trgt-step",
    #                 image=DOCKER_IMAGE,
    #                 cpu=2,
    #                 memory="highmem",
    #                 storage="20Gi",
    #                 output_dir=output_dir)
    
    #s2.depends_on(s1)

    #s2.command("mkdir /io/run_dir; cd /io/run_dir")
    #for json_path in s1_output_json_paths:
    #    local_path = s2.input(json_path, localize_by=Localize.COPY)
    #    s2.command(f"ln -s {local_path}")

    #s2.command("set -x")
    #s2.command(f"python3 -m str_analysis.combine_str_json_to_tsv --include-extra-trgt-fields "
    #           f"--output-prefix {output_prefix}")

    #s2.command(f"mv {output_prefix}.{len(s1_output_json_paths)}_json_files.bed {output_prefix}.bed")
    #s2.command(f"mv {output_prefix}.{len(s1_output_json_paths)}_json_files.variants.tsv.gz {output_prefix}.variants.tsv.gz")
    #s2.command(f"mv {output_prefix}.{len(s1_output_json_paths)}_json_files.alleles.tsv.gz {output_prefix}.alleles.tsv.gz")

    #s2.command(f"bgzip {output_prefix}.bed")
    #s2.command(f"tabix {output_prefix}.bed.gz")
    #s2.command("ls -lhrt")
    #s2.output(f"{output_prefix}.variants.tsv.gz")
    #s2.output(f"{output_prefix}.alleles.tsv.gz")
    #s2.output(f"{output_prefix}.bed.gz")
    #s2.output(f"{output_prefix}.bed.gz.tbi")

    return s1

if __name__ == "__main__":
    main()


