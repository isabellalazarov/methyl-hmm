library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationDbi)
library(org.Hs.eg.db)

select <- dplyr::select

rna_seq <- read_csv("exports/rna_seq_merged.csv")
wgbs <- read_csv("exports/wgbs_merged.csv")

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

# Save to file (for easy loading)

# write_csv(merged_df, "exports/merged_final.csv")

# average methylation per gene
gene_summary <- summarize(
    group_by(merged_df, gene_id_entrez),
    avg_methylation = mean(beta, na.rm = TRUE),
    avg_expression = mean(expression, na.rm = TRUE)
)

# correlation using all genes including those w/o expression 
# some may be silent due to reasons other than methylation
cor_all <- cor(
    gene_summary$avg_methylation, gene_summary$avg_expression, use = "complete.obs")
paste("correlation (all genes) = ", round(cor_all, 4))

# correlation using only genes with expression > 0
gene_summary_expressed <- filter(gene_summary, avg_expression > 0)
cor_expressed <- cor(gene_summary_expressed$avg_methylation,
                     gene_summary_expressed$avg_expression,
                     use = "complete.obs")
paste("correlation (expressed genes only) = ", round(cor_expressed, 4))

ggplot(
    gene_summary_expressed,
    aes(x = avg_methylation, y = log1p(avg_expression))) +
    geom_point(alpha = 0.3, size = 0.5) +
    labs(x = "Average Promoter Methylation (%)",
         y = "Log Gene Expression/Log(Expected Count + 1)",
         title = "Promoter Methylation vs Gene Expression (Expressed Genes Only)")


set.seed(123)

cpgs_per_gene <- summarize(group_by(merged_df, gene_id_entrez), n_cpgs = n_distinct(probe_id))

# get 1000 genes with most CpG sites 
target_genes <- arrange(cpgs_per_gene, desc(n_cpgs)) 
target_genes <- slice_head(target_genes, n = 1000)
target_genes <- pull(target_genes, gene_id_entrez)

lasso_result <- NULL
lasso_result <- data.frame()
all_fits <- list()
for (gene in target_genes) {
    result <- fit_lasso_per_gene(gene, merged_df)
    if (!is.null(result)) {
        lasso_result <- bind_rows(lasso_result, result$summary)
        all_fits[[as.character(gene)]] <- result$cv_fit
    }
}

lasso_result

# better R^2 results than Zhong et al.,
# but we run lasso for 1000 genes as opposed to their 4000+
paste0("Percentage of genes with R^2 > 0.3 = ", mean(lasso_result$r_squared > 0.3, na.rm = TRUE)*100)
paste0("Percentage of genes with R^2 > 0.5 = ", mean(lasso_result$r_squared > 0.5, na.rm = TRUE)*100)
paste0("Percentage of genes with R^2 > 0.8 = ", mean(lasso_result$r_squared > 0.8, na.rm = TRUE)*100)

# signal visualization of worst (lowest R^2)... 
lasso_gene_rankings_asc <- arrange(lasso_result, r_squared)
worst_gene_lasso <- slice(lasso_gene_rankings_asc, 1)
worst_gene_lasso_id <- pull(worst_gene_lasso, gene_id_entrez)
worst_gene_r_squared <- pull(worst_gene_lasso, r_squared)

# ... and best (highest R^2) genes
lasso_gene_rankings_desc <- arrange(lasso_result, desc(r_squared))
best_gene_lasso <- slice(lasso_gene_rankings_desc, 1)
best_gene_lasso_id <- pull(best_gene_lasso, gene_id_entrez)
best_gene_lasso_r_squared <- pull(best_gene_lasso, r_squared)

best_fit <- all_fits[[as.character(best_gene_lasso_id)]]
plot(best_fit)
title("LASSO fit: best gene", line = 2.5)
paste0("best gene r^2 = ", best_gene_lasso_r_squared)

worst_fit <- all_fits[[as.character(worst_gene_lasso_id)]]
plot(worst_fit)
title("LASSO fit: worst gene", line = 2.5)
paste0("worst gene r^2 = ", worst_gene_r_squared)

