library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationDbi)
library(org.Hs.eg.db)

select <- dplyr::select

rna_seq <- read_csv("exports/rna_seq_merged.csv")
wgbs <- read_tsv("exports/wgbs_k562_cpg.bed.gz",
                  col_names = c("chrom", "start", "end", "name", "score",
                               "strand", "thickStart", "thickEnd", "rgb",
                               "coverage", "methylation_percent")) %>%
  select(chrom, start, end, strand, methylation_percent) %>%
  mutate(accession = "ENCSR765JPC")

hg38_txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

id_type_counts <- rna_seq %>%
  mutate(
    id_raw = as.character(gene_id),
    id_no_ver = sub("\\..*$", "", id_raw),
    id_type = case_when(
      str_detect(id_no_ver, "^ENSG[0-9]+$") ~ "ensembl_gene",
      str_detect(id_no_ver, "^[0-9]+$") ~ "entrez_numeric",
      TRUE ~ "other"
    )
  ) %>%
  count(id_type, sort = TRUE)

id_type_counts

rna_std <- rna_seq %>%
  mutate(
    gene_id_raw = as.character(gene_id),
    gene_id_no_ver = sub("\\..*$", "", gene_id_raw),
    id_type = case_when(
      str_detect(gene_id_no_ver, "^ENSG[0-9]+$") ~ "ensembl_gene",
      str_detect(gene_id_no_ver, "^[0-9]+$") ~ "entrez_numeric",
      str_detect(gene_id_no_ver, "^ENST[0-9]+$") ~ "ensembl_transcript",
      TRUE ~ "other"
    )
  )

ensg_keys <- unique(rna_std$gene_id_no_ver[rna_std$id_type == "ensembl_gene"])
ensg_to_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ensg_keys,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

rna_std <- rna_std %>%
  mutate(
    gene_id_entrez = case_when(
      id_type == "entrez_numeric" ~ gene_id_no_ver,
      id_type == "ensembl_gene" ~ as.character(ensg_to_entrez[gene_id_no_ver]),
      TRUE ~ NA_character_
    )
  )

# Extract promoter regions

wgbs_gr <- GRanges(seqnames = wgbs$chrom,
        ranges = IRanges(start = wgbs$start, end = wgbs$end),
        strand = wgbs$strand,
        accession = wgbs$accession,
        methylation_percent = wgbs$methylation_percent)

promoter_regions <- promoters(genes(hg38_txdb, single.strand.genes.only=FALSE), upstream=2000, downstream=400)

wgbs$in_promoter <- overlapsAny(wgbs_gr, promoter_regions, ignore.strand = TRUE)

wgbs_filtered <- wgbs %>%
    filter(in_promoter == TRUE) %>%
    select(-c(in_promoter))

hg38_genes <- genes(hg38_txdb)
rna_seq_ranges <- hg38_genes[hg38_genes$gene_id %in%
                unique(na.omit(rna_std$gene_id_entrez))]

rna_seq_promoters <- promoters(rna_seq_ranges, upstream=2000, downstream=400)

wgbs_filtered_gr <- GRanges(seqnames = wgbs_filtered$chrom,
        ranges = IRanges(start = wgbs_filtered$start, end = wgbs_filtered$end),
        strand = wgbs_filtered$strand,
        accession = wgbs_filtered$accession,
        methylation_percent = wgbs_filtered$methylation_percent)

hits <- findOverlaps(rna_seq_promoters, wgbs_filtered_gr, ignore.strand = TRUE)

wgbs_df <- data.frame(
    wgbs_accession = mcols(wgbs_filtered_gr)$accession[subjectHits(hits)],
    gene_id_entrez = mcols(rna_seq_promoters)$gene_id[queryHits(hits)],
    chrom = seqnames(wgbs_filtered_gr)[subjectHits(hits)],
    start = start(wgbs_filtered_gr)[subjectHits(hits)],
    end = end(wgbs_filtered_gr)[subjectHits(hits)],
    strand = strand(wgbs_filtered_gr)[subjectHits(hits)],
    methylation_percent = 
        mcols(wgbs_filtered_gr)$methylation_percent[subjectHits(hits)]
)

rna_std <- rna_std %>%
    filter(!is.na(gene_id_entrez))

merged_df <- inner_join(wgbs_df, rna_std, by = c("gene_id_entrez"))

merged_df <- merged_df %>%
    select(-c("gene_id", "gene_id_raw", "gene_id_no_ver", "id_type")) %>%
    group_by(gene_id_entrez, chrom, start, end) %>%
    summarize(
        mean_methylation = mean(methylation_percent),
        mean_expression = mean(expected_count)
    )

# LASSO analysis

# for glmnet format: one row per gene, one column per CpG site
merged_wide <- pivot_wider(
  ungroup(
    mutate(
      group_by(merged_df, gene_id_entrez),
      cpg_index = paste0("cpg_", row_number())
    )
  ),
  id_cols = c(gene_id_entrez, mean_expression),
  names_from = cpg_index,
  values_from = mean_methylation,
  values_fill = 0
)

X <- as.matrix(select(merged_wide, starts_with("cpg_")))
Y <- merged_wide$mean_expression

set.seed(123)
k <- 10
foldid <- sample(rep(1:k, length.out = nrow(merged_wide)))

cv_fit <- cv.glmnet(
  x = X,
  y = Y,
  family = "gaussian",
  alpha = 1,
  foldid = foldid,
  standardize = TRUE
)

best_lambda = cv_fit$lambda.min

plot(cv_fit)

predicted <- predict(cv_fit, newx = X, s = "lambda.min")
lasso_r_squared <- cor(Y, predicted)^2
paste("r^2 value = ", lasso_r_squared)

# Correlation analysis

# average methylation per gene
gene_summary <- summarize(
  group_by(merged_df, gene_id_entrez),
  avg_methylation = mean(mean_methylation),
  avg_expression = mean(mean_expression)
  )

# correlation using all genes including those w/o expression 
# some may be silent due to reasons other than methylation
cor_all <- cor(gene_summary$avg_methylation, gene_summary$avg_expression)
paste("correlation (all genes) = ", cor_all)

# correlation using only genes with expression > 0
gene_summary_expressed <- filter(gene_summary, avg_expression > 0)
cor_expressed <- cor(gene_summary_expressed$avg_methylation, 
     gene_summary_expressed$avg_expression)
paste("correlation (expressed genes only) = ", cor_expressed)

ggplot(gene_summary_expressed,
       aes(x = avg_methylation, y = log1p(avg_expression))) + 
  geom_point(alpha = 0.3, size = 0.5) +
  labs(x = "Average Promoter Methylation (%)",
       y = "Log Gene Expression/Log(Expected Count + 1)",
       title = "Promoter Methylation vs Gene Expression (Expressed Genes Only)")
