fit_lasso_per_gene <- function(target_gene, merged_df) {
  
  gene_data <- filter(merged_df, gene_id_entrez == target_gene)
  
  n_cpgs <- length(unique(gene_data$probe_id))
  if (n_cpgs < 2) {
    return(NULL)
  } 
  
  gene_data <- arrange(gene_data, start)
  unique_cpgs <- unique(gene_data$probe_id)
  gene_data$cpg_index <- paste0("CpG_", match(gene_data$probe_id, unique_cpgs))
  
  wide <- pivot_wider(select(gene_data, sample_id, cpg_index, beta, expression),
                      names_from = cpg_index,
                      values_from = beta)
  wide <- drop_na(wide)
  n_samples <- nrow(wide)
  if (n_samples < 100) {
    return(NULL)
  }
  
  X <- as.matrix(select(wide, starts_with("CpG_")))
  Y <- log1p(wide$expression)
  
  set.seed(123)
  train_ids <- sample(nrow(wide), floor(0.8 * nrow(wide)))
  X_train <- X[train_ids, ]
  X_test <- X[-train_ids, ]
  Y_train <- Y[train_ids]
  Y_test <- Y[-train_ids]
  
  
  if (var(Y_train) == 0) {
    return(NULL)
  }
  
  cv_fit <- tryCatch(
    cv.glmnet(x = X_train,
              y= Y_train,
              family = "gaussian",
              alpha = 1,
              standardize = TRUE,
              nfolds = 10),
    error = function(e) NULL)
  
  if (is.null(cv_fit)) {
    return(NULL)
  }
  
  predicted_expression <- predict(cv_fit, newx = X_test, s = "lambda.min")
  r_squared <- as.numeric(cor(Y_test, predicted_expression)^2)
  
  list(
    summary = data.frame(gene_id_entrez = target_gene,
                         n_cpgs = n_cpgs,
                         n_samples = n_samples,
                         r_squared = r_squared),
    cv_fit = cv_fit
  )
  
  
  
}