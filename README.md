# Predicting Gene Expression with Methylation-GEX Relationships Across Promoter Regions with an HMM

## Data 

Data was obtained from TCGA-COAD. 

## Project Structure

```
methyl-hmm
├── README.md
├── download_encode.py
├── exports
│   ├── rna_seq_manifest.csv
│   ├── rna_seq_merged.csv
│   ├── wgbs_manifest.csv
│   └── wgbs_merged.csv
├── fit_hmm_per_gene.R
├── fit_lasso_per_gene.R
├── hmm_analyses.Rmd
├── methyl_hmm.Rproj
├── old
│   ├── fit_hmm.R
│   ├── old_hmm_analysis.R
│   └── per_sample_table.R
└── sample_input.csv
```

`hmm_analyses.Rmd`: 

`fit_lasso_per_gene.R`: File containing a function to fit a LASSO model for a gene. Takes in the data frame and a target gene ID, and returns the fitted model and a summary data frame including the provided gene ID, number of CpGs and samples found for that gene, and the $R^2$ of the learned model after it was run on the test set, or `NULL` if the function failed at any point. This was used to evaluate the performance of the LASSO model for each gene.

`fit_hmm_per_gene.R`: File containing a function to fit an HMM for a gene. 

`methyl_hmm.Rproj`: 

`sample_input.csv`: Small sample input file

`download_encode.py`: do we want to talk about this?

`exports/`: Folder containing the processed data files from `download_encode.py`. should we talk about each individual file? (remove if we don't want to include the old results)

`old/`: (can remove this and all of the below files)

`old/fit_hmm.R`:

`old/old_hmm_analysis.R`:

`old/per_sample_table.R`:

## Reproducing Results

R version 4.5.2 was used with the following packages installed, as well as their dependencies: 
- `tidyverse` version 2.0.0
- `GenomicRanges` version 1.62.1
- `TxDb.Hsapiens.UCSC.hg38.knownGene` version 3.22.0
- `AnnotationDbi` version 1.72.0
- `org.Hs.eg.db` version 3.22.0
- `glmnet` version 4.1-10
- `depmixS4` version 1.5-1
- `here` version 1.0.2

To run 




requirements (DELETE):
• All source code/scripts. This includes any code used to generate display items used in your report, including supplementary materials.
• Sample, small input and output files.
• If your project generally operates on “large” data files that are available publicly, list the names of these files and where to get them.
• A ReadMe file that tells exactly what is each file, how to compile (if relevant) and test-run the project on sample inputs to get the sample output, how to really run the project on large files, if relevant, what are the parameters, ~~what are the system requirements (e.g. are you using a particular version of matlab/Python/R/Ruby/Java/C++? Do you run on a particularly powerful machine/cloud instance? Which standard or add-on libraries would you need to have been previously installed?)~~