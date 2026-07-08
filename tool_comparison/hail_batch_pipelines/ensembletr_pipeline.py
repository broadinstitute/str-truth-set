import hailtop.fs as hfs
import os

from step_pipeline import pipeline, Backend, Localize

# str-analysis-with-ensembletr image (built by docker_with_ensembletr/ in broadinstitute/str-analysis via
# build_docker_images.yml): htslib/samtools/bcftools + bw2/EnsembleTR (fork of gymrek-lab/EnsembleTR with a fix
# for a crash on hemizygous male chrX/chrY genotype calls, see docker_with_ensembletr/Dockerfile) + str_analysis
# (for the convert_expansion_hunter_json_to_vcf / convert_ensembletr_vcf_to_expansion_hunter_json /
# combine_str_json_to_tsv steps below). Pinned by @sha256 digest (from docker_with_ensembletr/sha256_dockerhub.txt)
# like the sibling pipelines.
DOCKER_IMAGE = "weisburd/str-analysis-with-ensembletr@sha256:390a0f3952422060b18fa23165ac28396fbf3a9ca5ccdabcc7258a6a2cce26e5"

REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"


def create_ensembletr_steps(bp, *, reference_fasta, eh_json_paths, hipstr_vcf_paths, variant_catalog_path, output_dir,
                            output_prefix, sample_id, gangstr_vcf_paths=None, reference_fasta_fai=None,
                            male_or_female="female", catalog_step=None):
    """Build the steps that merge per-caller genotypes with EnsembleTR and convert the consensus to the truth-set format.

    EnsembleTR (https://github.com/gymrek-lab/EnsembleTR) merges the already-computed ExpansionHunter, HipSTR and
    (optionally) GangSTR genotypes for one sample into a single consensus VCF. Downstream str-truth-set comparison
    consumes ExpansionHunter-format json, so this runs, per (sample, data_type, coverage):

      step1 (EnsembleTR merge + convert):
        1. reconstruct one ExpansionHunter-style VCF from the EH json shard(s) with
           str_analysis.convert_expansion_hunter_json_to_vcf, then bgzip + tabix it;
        2. concatenate + sort each of the HipSTR and (optional) GangSTR native VCF shards into one bgzipped, tabixed VCF
           per caller, forcing the sample column to {sample_id} so all input VCFs share one sample for EnsembleTR;
        3. run `EnsembleTR --vcfs <eh>,<hipstr>[,<gangstr>] --ref <fasta> --out <prefix>`;
        4. reconcile the consensus VCF back to truth-set LocusIds and write EH-format json with
           str_analysis.convert_ensembletr_vcf_to_expansion_hunter_json --variant-catalog <catalog>.

      step2 (combine): str_analysis.combine_str_json_to_tsv --include-extra-ensembletr-fields on that json, producing the
        {output_prefix}.1_json_files.variants.tsv.gz / .alleles.tsv.gz the downstream add-columns step reads.

    Args:
        bp: the step_pipeline pipeline object.
        reference_fasta: gs:// path of the reference fasta (EnsembleTR reads reference allele sequences from it).
        eh_json_paths: list of gs:// ExpansionHunter output json (or json.gz) paths for this sample (one or more shards).
        hipstr_vcf_paths: list of gs:// HipSTR native VCF (.vcf.gz) paths for this sample (one or more shards).
        variant_catalog_path: gs:// path of the unsharded ExpansionHunter (EHv5) variant catalog json the upstream
            callers were run against, used to reconcile EnsembleTR records back to truth-set LocusIds.
        output_dir: directory the outputs are written under (json/ and vcf/ subdirs plus the combined tables).
        output_prefix: prefix for the combined step-2 outputs (e.g. "{sample}.EnsembleTR-EH+HipSTR").
        sample_id: sample id, written into every input VCF's sample column and the output json.
        gangstr_vcf_paths: list of gs:// GangSTR native VCF (.vcf.gz) paths, or None for the EH+HipSTR (2-caller) mode.
        reference_fasta_fai: gs:// path of the reference fasta .fai index (defaults to reference_fasta + ".fai").
        male_or_female: retained for a uniform create_*_steps signature; EnsembleTR has no sex option, and hemizygous
            (male chrX/chrY) calls are made diploid inside convert_expansion_hunter_json_to_vcf.
        catalog_step: optional upstream Step that produces the input catalog; the genotyping step(s) depend on it so
            they never run before the catalog exists.

    Returns:
        The combine step (step2), whose .variants.tsv.gz output the caller passes to add_tool_comparison_columns_step.
    """
    if reference_fasta_fai is None:
        reference_fasta_fai = reference_fasta + ".fai"
    if not eh_json_paths:
        raise ValueError(f"No ExpansionHunter json paths given for EnsembleTR on {sample_id}")
    if not hipstr_vcf_paths:
        raise ValueError(f"No HipSTR vcf paths given for EnsembleTR on {sample_id}")

    num_callers = 3 if gangstr_vcf_paths else 2

    # step1: reconstruct the EH VCF, normalize the HipSTR/GangSTR VCFs, run EnsembleTR, convert the consensus to json.
    # cpu=1/standard right-sizes this (measured via /usr/bin/time --verbose below: EnsembleTR itself is
    # single-threaded, ~10min wall/user time, peak RSS ~125MB on the 1-sample HG002 31x test run) instead of the
    # 2/highmem (13GB) it was over-provisioned with before that measurement existed.
    s1 = bp.new_step(
        f"Run EnsembleTR ({num_callers} callers) on {sample_id}",
        arg_suffix="run-ensembletr-step",
        step_number=1,
        image=DOCKER_IMAGE,
        cpu=1,
        memory="standard",
        storage="20Gi",
        localize_by=Localize.GSUTIL_COPY,
        output_dir=output_dir)

    if catalog_step is not None:
        s1.depends_on(catalog_step)

    local_fasta = s1.input(reference_fasta)
    s1.input(reference_fasta_fai)
    local_catalog = s1.input(variant_catalog_path)
    local_eh_jsons = [s1.input(p) for p in eh_json_paths]
    local_hipstr_vcfs = [s1.input(p) for p in hipstr_vcf_paths]
    local_gangstr_vcfs = [s1.input(p) for p in gangstr_vcf_paths] if gangstr_vcf_paths else []

    s1.command("set -exuo pipefail")
    s1.command("mkdir -p /io/run_dir && cd /io/run_dir")
    s1.command(f"echo {sample_id} > sample_name.txt")

    # 1. reconstruct a single ExpansionHunter-style VCF from the (possibly sharded) EH json; the converter tags it with
    #    the ##ALT=<ID=STR..> + INFO VARID/RU/RL + FORMAT REPCN fields EnsembleTR's ExpansionHunter harmonizer needs.
    eh_json_args = " ".join(str(x) for x in local_eh_jsons)
    s1.command(f"python3.9 -m str_analysis.convert_expansion_hunter_json_to_vcf "
               f"--sample-id {sample_id} -o expansion_hunter.vcf {eh_json_args}")
    s1.command("bgzip -f expansion_hunter.vcf && tabix -f -p vcf expansion_hunter.vcf.gz")

    # 2. per caller: normalize each shard to a bgzipped + indexed VCF (the native HipSTR/GangSTR shards can be plain-gzip
    #    and unindexed), concat (-a tolerates overlapping/unordered shards), sort, and force the sample column to
    #    {sample_id} so all three caller VCFs share one sample for EnsembleTR to merge on.
    for caller_name, local_vcfs in [("hipstr", local_hipstr_vcfs), ("gangstr", local_gangstr_vcfs)]:
        if not local_vcfs:
            continue
        normalized_shards = []
        for i, local_vcf in enumerate(local_vcfs):
            s1.command(f"bcftools view -Oz -o {caller_name}.shard{i}.vcf.gz {local_vcf}")
            s1.command(f"tabix -f -p vcf {caller_name}.shard{i}.vcf.gz")
            normalized_shards.append(f"{caller_name}.shard{i}.vcf.gz")
        s1.command(f"bcftools concat -a -Oz -o {caller_name}.concat.vcf.gz {' '.join(normalized_shards)}")
        s1.command(f"bcftools sort -Oz -o {caller_name}.sorted.vcf.gz {caller_name}.concat.vcf.gz")
        s1.command(f"bcftools reheader -s sample_name.txt -o {caller_name}.vcf.gz {caller_name}.sorted.vcf.gz")
        s1.command(f"tabix -f -p vcf {caller_name}.vcf.gz")

    # 3. run EnsembleTR over the reconstructed EH VCF + the normalized HipSTR (+ GangSTR) VCFs
    vcf_list = "expansion_hunter.vcf.gz,hipstr.vcf.gz" + (",gangstr.vcf.gz" if local_gangstr_vcfs else "")
    # EnsembleTR requires --out to end in "vcf" and writes to that path verbatim (main.py: endswith("vcf") guard +
    # Writer(args.out)), so the output file is named ensembletr_output.vcf.
    s1.command(f"/usr/bin/time --verbose EnsembleTR --vcfs {vcf_list} --ref {local_fasta} --out ensembletr_output.vcf")
    s1.command("ls -lhrt")

    # 4. reconcile the consensus VCF (ensembletr_output.vcf) back to truth-set LocusIds and write EH-format json
    json_filename = f"{output_prefix}.json"
    s1.command(f"python3.9 -m str_analysis.convert_ensembletr_vcf_to_expansion_hunter_json "
               f"--variant-catalog {local_catalog} --sample-id {sample_id} "
               f"-o {json_filename} ensembletr_output.vcf")
    s1.command("gzip -f ensembletr_output.vcf")
    s1.command("ls -lhrt")
    s1.output("ensembletr_output.vcf.gz", output_dir=os.path.join(output_dir, "vcf"))
    s1.output(json_filename, output_dir=os.path.join(output_dir, "json"))
    step1_json_output_path = os.path.join(output_dir, "json", json_filename)

    # step2: combine the single json into the .variants.tsv.gz / .alleles.tsv.gz tables the downstream comparison
    # reads. cpu=2/standard right-sizes this (the 1-sample HG002 31x test run combined a 259MB json of 568,979 loci
    # in ~2.4min single-threaded pandas work) instead of the 4/highmem (26GB) it was over-provisioned with before.
    s2 = bp.new_step(
        name=f"Combine EnsembleTR output for {sample_id}",
        step_number=2,
        arg_suffix="combine-ensembletr-step",
        image=DOCKER_IMAGE,
        cpu=2,
        memory="standard",
        storage="20Gi",
        localize_by=Localize.GSUTIL_COPY,
        output_dir=output_dir)
    s2.depends_on(s1)

    s2.command("mkdir /io/run_dir; cd /io/run_dir")
    local_json = s2.input(step1_json_output_path)
    s2.command(f"ln -s {local_json}")
    s2.command("set -x")
    s2.command(f"python3.9 -m str_analysis.combine_str_json_to_tsv --include-extra-ensembletr-fields "
               f"--output-prefix {output_prefix}")
    s2.command(f"bgzip {output_prefix}.1_json_files.bed")
    s2.command(f"tabix {output_prefix}.1_json_files.bed.gz")
    s2.command("ls -lhrt")
    s2.output(f"{output_prefix}.1_json_files.variants.tsv.gz")
    s2.output(f"{output_prefix}.1_json_files.alleles.tsv.gz")
    s2.output(f"{output_prefix}.1_json_files.bed.gz")
    s2.output(f"{output_prefix}.1_json_files.bed.gz.tbi")

    return s2


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--male-or-female", choices=["male", "female"], default="female")
    parser.add_argument("--eh-json", action="append", required=True,
                        help="gs:// path (or glob) of the ExpansionHunter output json for this sample. Repeatable.")
    parser.add_argument("--hipstr-vcf", action="append", required=True,
                        help="gs:// path (or glob) of the HipSTR native VCF(s) for this sample. Repeatable.")
    parser.add_argument("--gangstr-vcf", action="append",
                        help="gs:// path (or glob) of the GangSTR native VCF(s). Omit for the EH+HipSTR (2-caller) mode.")
    parser.add_argument("--variant-catalog", required=True,
                        help="gs:// path of the unsharded ExpansionHunter (EHv5) variant catalog json.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-prefix", required=True)
    args = bp.parse_known_args()

    def _resolve(globs):
        paths = []
        for glob in globs or []:
            paths.extend(x.path for x in hfs.ls(glob))
        return sorted(paths)

    bp.set_name(f"STR Truth Set: EnsembleTR: {args.sample_id}")
    create_ensembletr_steps(
        bp,
        reference_fasta=args.reference_fasta,
        reference_fasta_fai=args.reference_fasta_fai,
        eh_json_paths=_resolve(args.eh_json),
        hipstr_vcf_paths=_resolve(args.hipstr_vcf),
        gangstr_vcf_paths=_resolve(args.gangstr_vcf) or None,
        variant_catalog_path=args.variant_catalog,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        sample_id=args.sample_id,
        male_or_female=args.male_or_female)
    bp.run()


if __name__ == "__main__":
    main()
