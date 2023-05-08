import argparse
import gzip
import os
import re

p = argparse.ArgumentParser()
p.add_argument("-i", "--input-bed", default="./STR_truth_set.v1.variants.bed.gz", help="Truth set bed file")
args = p.parse_args()

if not os.path.isfile(args.input_bed):
    p.error(f"{args.input_bed} not found")

output_bed = re.sub(".bed(.gz|.bgz)?", "", args.input_bed) + ".with_reviewer_image_urls.bed"
fopen = gzip.open if args.input_bed.endswith("gz") else open
with fopen(args.input_bed, "rt") as f, open(output_bed, "wt") as f2:
    for i, line in enumerate(f):
        fields = line.strip().split("\t")
        chrom = fields[0]
        start_0based = fields[1]
        end = fields[2]
        name = fields[3]
        if ":" in name:
            name_fields = name.split(":")
            motif = name_fields[1]
        else:
            motif = name
        url = "https://storage.googleapis.com/str-truth-set/hg38/tool_results/expansion_hunter/positive_loci/svg/CHM1_CHM13_WGS2."
        url += chrom.replace("chr", "") + f"-{start_0based}-{end}-{motif}"
        url += ".svg"

        name = f"Label={name}; REViewer=<a href={url}>image</a>;"
        f2.write("\t".join([chrom, start_0based, end, name, fields[4]]) + "\n")

os.system(f"bgzip -f {output_bed}")
os.system(f"tabix {output_bed}.gz")
print(f"Wrote {i+1:,d} rows to {output_bed}.gz")
#%%
