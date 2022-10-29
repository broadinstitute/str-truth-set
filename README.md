**STRuth v1.0**  

A genome-wide short tandem repeat (STR) truth set based on the 
[Synthetic Diploid Benchmark](https://github.com/lh3/CHM-eval) [[Li et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6341484/)]


---
### Creating the Truth Set

This repo contains the code used to generate the STR truth set.  
To recreate the truth set files:

```
# clone this repo
git clone git@github.com:broadinstitute/str-truth-set.git

cd str-truth-set

# download reference data and create the truth set
./create_STR_truth_set.sh  
```




