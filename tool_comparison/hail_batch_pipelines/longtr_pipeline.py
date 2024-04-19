"""This pipeline runs LongTR (https://github.com/gymrek-lab/LongTR.git)

Usage: LongTR --bams <list_of_bams> --fasta <genome.fa> --regions <region_file.bed> --tr-vcf <tr_gts.vcf.gz> [OPTIONS]

Required parameters:
	--bams          <list_of_bams>        	Comma separated list of BAM/CRAM files. Either --bams or --bam-files must be specified
	--fasta         <genome.fa>           	FASTA file containing all of the relevant sequences for your organism
	                                      	  When analyzing CRAMs, this FASTA file must match the file used for compression
	--regions       <region_file.bed>     	BED file containing coordinates for each TR region
	--tr-vcf       <tr_gts.vcf.gz>      	Bgzipped VCF file to which TR genotypes will be written

Optional input parameters:
	--bam-files  <bam_files.txt>          	File containing BAM/CRAM files to analyze, one per line
	--ref-vcf    <tr_ref_panel.vcf.gz>   	Bgzipped input VCF file of a reference panel of TR genotypes. VCF alleles will be
	                                      	 used as candidate variants instead of finding candidates in the BAMs/CRAMs (Default)
	--snp-vcf    <phased_snps.vcf.gz>     	Bgzipped input VCF file containing phased SNP genotypes for the samples
	                                      	 to be genotyped. These SNPs will be used to physically phase TRs
	--skip-assembly                       	Skip assembly for genotyping with long reads
	--min-mean-qual       <threshold>     	Minimum average quality threshold for sequencing read data (Default = 30)
	--min-mapq            <threshold>     	Minimum MAPQ per read (Default = 20)
	--stutter-align-len   <threshold>     	Use stutter alignment for repeats with length less than threshold (Default = 0)
	--phased-bam	                      	Use phasing information from haplotagged sequencing data.
	--indel-flank-len     <max_bp>        	Include InDels in max_bp base pair around repeat as part of the repeath (Default = 5)
Optional output parameters:
	--log           <log.txt>             	Output the log information to the provided file (Default = Standard error)
Optional VCF formatting parameters:
	--output-gls                          	Write genotype likelihoods to the VCF (Default = False)
	--output-pls                          	Write phred-scaled genotype likelihoods to the VCF (Default = False)
	--output-phased-gls                   	Write phased genotype likelihoods to the VCF (Default = False)
	--output-filters                      	Write why individual calls were filtered to the VCF (Default = False)

Optional BAM/CRAM tweaking parameters:
	--bam-samps     <list_of_samples>     	Comma separated list of read groups in same order as BAM/CRAM files.
	                                      	  Assign each read the read group corresponding to its file. By default,
	                                      	  each read must have an RG tag and the sample is determined from the SM field
	--bam-libs      <list_of_libraries>   	Comma separated list of libraries in same order as BAM/CRAM files.
	                                      	  Assign each read the library corresponding to its file. By default,
	                                      	  each read must have an RG tag and the library is determined from the LB field
	--lib-from-samp                       	 Assign each read the library corresponding to its sample name. By default,
	                                      	  each read must have an RG tag and the library is determined from the LB field

Optional haplotype filtering parameters:
	--max-haps <max_haplotypes>           	Maximum allowable candidate haplotypes for an TR (Default = 1000)
	                                      	 Loci with more candidate haplotypes will not be genotyped
Other optional parameters:
	--help                                	Print this help message and exit
	--version                             	Print LongTR version and exit
	--quiet                               	Only output terse logging messages (Default = output all messages)
	--silent                              	Don't output any logging messages  (Default = output all messages)
	--chrom              <chrom>          	Only consider TRs on this chromosome
	--haploid-chrs       <list_of_chroms> 	Comma separated list of chromosomes to treat as haploid (Default = all diploid)
	--hap-chr-file       <hap_chroms.txt> 	File containing chromosomes to treat as haploid, one per line
	--min-reads          <num_reads>      	Minimum total reads required to genotype a locus (Default = 10)
	--max-reads          <num_reads>      	Skip a locus if it has more than NUM_READS reads (Default = 1000000)
	--max-tr-len         <max_bp>         	Only genotype TRs in the provided BED file with length < MAX_BP (Default = 1000)
"""

import hail as hl
import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/longtr@sha256:1b4ab7bba5c4642cb938ee7c60fb388b3ebdc99b5de09e0fba15f44e4750b4c1"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

REGIONS_BED_POSITIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/hipstr/positive_loci.*_of_015.bed"
REGIONS_BED_NEGATIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/hipstr/negative_loci.*_of_015.bed"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/longtr"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser_group = parser.add_mutually_exclusive_group(required=True)
    parser_group.add_argument("--positive-loci", action="store_true", help="Genotype truth set loci")
    parser_group.add_argument("--negative-loci", action="store_true", help="Genotype negative (hom-ref) loci")
    parser_group.add_argument("--regions-bed", action="append", help="Path of regions bed file(s) to process")

    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_CRAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    args = bp.parse_known_args()

    if args.positive_loci:
        positive_or_negative_loci = "positive_loci"
        regions_bed_file_stats_list = hl.hadoop_ls(REGIONS_BED_POSITIVE_LOCI)
        if len(regions_bed_file_stats_list) == 0:
            raise ValueError(f"No files found matching {REGIONS_BED_POSITIVE_LOCI}")
    elif args.negative_loci:
        positive_or_negative_loci = "negative_loci"
        regions_bed_file_stats_list = hl.hadoop_ls(REGIONS_BED_NEGATIVE_LOCI)
        if len(regions_bed_file_stats_list) == 0:
            raise ValueError(f"No files found matching {REGIONS_BED_NEGATIVE_LOCI}")
    elif args.regions_bed:
        positive_or_negative_loci = os.path.basename(args.regions_bed[0]).replace(".bed", "").replace(".gz", "")
        regions_bed_file_stats_list = [{"path": path} for path in args.regions_bed]
    else:
        parser.error("Must specify either --positive-loci or --negative-loci")

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: LongTR: {positive_or_negative_loci}: {bam_path_ending}")
    output_dir = os.path.join(args.output_dir, positive_or_negative_loci)
    if not args.force:
        json_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.json"))
        logging.info(f"Precached {len(json_paths)} json files")

    step1s = []
    step1_output_json_paths = []
    for repeat_spec_i, regions_bed_file_stats in enumerate(regions_bed_file_stats_list):
        regions_bed_path = regions_bed_file_stats["path"]

        if args.n and repeat_spec_i >= args.n:
            break

        s1 = bp.new_step(f"Run LongTR #{repeat_spec_i}",
                         arg_suffix=f"longtr",
                         step_number=1,
                         image=DOCKER_IMAGE,
                         cpu=1,
                         storage="200Gi")
        step1s.append(s1)

        local_fasta = s1.input(args.reference_fasta, localize_by=Localize.COPY)
        if args.reference_fasta_fai:
            s1.input(args.reference_fasta_fai, localize_by=Localize.COPY)

        local_bam = s1.input(args.input_bam, localize_by=Localize.COPY)
        if args.input_bai:
            s1.input(args.input_bai, localize_by=Localize.COPY)

        local_regions_bed = s1.input(regions_bed_path)

        output_prefix = re.sub(".bed(.gz)?$", "", local_regions_bed.filename)
        s1.command(f"echo Genotyping $(cat {local_regions_bed} | wc -l) loci")
        s1.command("set -ex")
        s1.command(f"""/usr/bin/time --verbose LongTR \
                --min-reads 2 \
                --bams {local_bam} \
                --fasta {local_fasta} \
                --regions {local_regions_bed} \
                --log {output_prefix}.log \
                --tr-vcf {output_prefix}.vcf.gz""")

        s1.command("ls -lhrt")

        s1.command(f"python3.9 -m str_analysis.convert_hipstr_vcf_to_expansion_hunter_json {output_prefix}.vcf.gz")
        s1.command(f"gzip {output_prefix}.log")

        s1.output(f"{output_prefix}.vcf.gz", output_dir=os.path.join(output_dir, f"vcf"))
        s1.output(f"{output_prefix}.log.gz", output_dir=os.path.join(output_dir, f"log"))
        s1.output(f"{output_prefix}.json", output_dir=os.path.join(output_dir, f"json"))

        step1_output_json_paths.append(os.path.join(output_dir, f"json", f"{output_prefix}.json"))

    # step2: combine json files
    s2 = bp.new_step(name="Combine LongTR outputs",
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
    s2.command(f"python3.9 -m str_analysis.combine_str_json_to_tsv --include-extra-longtr-fields "
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


