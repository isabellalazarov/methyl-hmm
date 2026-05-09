# Predicting Gene Expression with Methylation-GEX Relationships Across Promoter Regions with an HMM

some description here?? idk

## Data 

Data was obtained from TCGA-COAD. The raw and preprocessed data files were not included in this repository due to their size, but they can be obtained through the following steps:

1. **Download raw data**
   - Files will be downloaded from the following respositories: [methylation](https://xenabrowser.net/datapages/?dataset=TCGA-COAD.methylation450.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443), [gene expression](https://xenabrowser.net/datapages/?dataset=TCGA-COAD.star_fpkm-uq.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443), and [phenotype](https://xenabrowser.net/datapages/?dataset=TCGA-COAD.clinical.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443).
   - From the methylation repository, download the methylation assay data file [TCGA-COAD.methylation450.tsv.gz](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.methylation450.tsv.gz) (next to "download" on the website) and rename it `methylation_data.tsv.gz`. Then, download the probe mapping [HM450.hg38.manifest.gencode.v36.probeMap](https://gdc-hub.s3.us-east-1.amazonaws.com/download/HM450.hg38.manifest.gencode.v36.probeMap) (next to "ID/Gene Mapping" on the website) and rename it `probe_map`.
   - From the gene expression repository, download the RNA-seq data file [TCGA-COAD.star_fpkm-uq.tsv.gz](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.star_fpkm-uq.tsv.gz) (next to "download" on the website) and call it `gene_expression_data.tsv.gz`.
   - From the phenotype repository, download the clinical data file [TCGA-COAD.clinical.tsv.gz](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.clinical.tsv.gz) (next to "download" on the website) and call it `phenotype_data.tsv.gz`.
2. **Organize directory structure**
   - Create a `data/` folder with two subfolders, `methylation/` and `gene_expression/`. 
   - Inside `methylation/`, place `methylation_data.tsv.gz` and `probe_map`.
   - Inside `gene_expression/`, place `gene_expression_data.tsv.gz`.
   - In `data/` but outside the two subfolders, place `phenotype_data.tsv.gz`.
   - The final directory structure should match the following:
   ```

   ```
3. **Run preprocessing steps**
   - test
4. 

## Project Structure

```
methyl-hmm
├── README.md
├── data  (not in repository)
│   ├── methylation
│   │   ├── methylation_data.tsv.gz
│   │   └── probe_map
│   ├── gene_expression
│   │   └── gene_expression_data.tsv.gz
│   └── phenotype_data.tsv.gz
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

└── sample_input.csv
```

`hmm_analyses.Rmd`: 

`fit_lasso_per_gene.R`: File containing a function to fit a LASSO model for a gene. Takes in the data frame and a target gene ID, and returns the fitted model and a summary data frame including the provided gene ID, number of CpGs and samples found for that gene, and the $R^2$ of the learned model after it was run on the test set, or `NULL` if the function failed at any point. This was used to evaluate the performance of the LASSO model for each gene.

`fit_hmm_per_gene.R`: File containing a function to fit an HMM for a gene. 

`methyl_hmm.Rproj`: 

`sample_input.csv`: Small sample input file

*Old files:*

`download_encode.py`: Initial Python script used to 

`exports/`: Folder containing the processed data files from `download_encode.py`. should we talk about each individual file?

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



To produce the old results, Python 
Simply run `download_encode.py` to produce the output files in the `exports/` folder.




requirements (DELETE):
• All source code/scripts. This includes any code used to generate display items used in your report, including supplementary materials.
• Sample, small input and output files.
• If your project generally operates on “large” data files that are available publicly, list the names of these files and where to get them.
• A ReadMe file that tells exactly what is each file, how to compile (if relevant) and test-run the project on sample inputs to get the sample output, how to really run the project on large files, if relevant, what are the parameters, ~~what are the system requirements (e.g. are you using a particular version of matlab/Python/R/Ruby/Java/C++? Do you run on a particularly powerful machine/cloud instance? Which standard or add-on libraries would you need to have been previously installed?)~~