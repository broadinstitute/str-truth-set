**STRuth v1.0**  

A genome-wide short tandem repeat (STR) truth set based on the 
[Synthetic Diploid Benchmark](https://github.com/lh3/CHM-eval) [[Li et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6341484/)]

---

### Downloads

Truth set vcf and tsv files are available on the [release page](https://github.com/broadinstitute/str-truth-set/releases/tag/v1).

The hg38-aligned genome and exome sequencing files used in this truth set project are publicly available in Google Storage buckets maintained by the Broad Institute:

**Genome:** [gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram](https://console.cloud.google.com/storage/browser/broad-public-datasets/CHM1_CHM13_WGS2)  
**Exome:** [gs://broad-public-datasets/CHM1_CHM13_WES/CHMI_CHMI3_Nex1.cram](https://console.cloud.google.com/storage/browser/broad-public-datasets/CHM1_CHM13_WES)   

The original Illumina short-read genome sequencing data used in [<a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6341484">Li 2018</a>] for the CHM1-CHM13 sample is available from EBI @ https://www.ebi.ac.uk/ena/browser/view/ERR1341796

---
### Creating the Truth Set

Besides the truth set itself, this repo also contains the code used to generate it.  

To recreate the truth set files:

```
# clone this repo
git clone git@github.com:broadinstitute/str-truth-set.git

cd str-truth-set

# download reference data and create the truth set
./create_STR_truth_set.sh  
```

