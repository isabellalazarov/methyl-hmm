fit_hmm_per_gene <- function(gene, hmm_data) {
  
  gene_data <- filter(hmm_data, gene_id_entrez == gene)
  
  num_cpgs <- length(unique(gene_data$probe_id))
  # gene needs to have 10+ CpG sites for meaningful HMM traversal
  if (num_cpgs < 10) {
    return(NULL)
  } 
  
  #transforming (0, 1) dist to (-inf, inf)
  gene_data$beta_logit <- qlogis(pmin(pmax(gene_data$beta, 1e-4), 1 - 1e-4))
  gene_data <- gene_data[complete.cases(gene_data[ , c("beta_logit", "log_expression")]), ]
  
  all_samples <- unique(gene_data$sample_id)
  n_samples <- length(all_samples)
  if (n_samples < 100) {
    return(NULL)
  }
  
  set.seed(123)
  
  # 80% train, 20% test
  train_sample_ids <- sample(all_samples, floor(0.8 * n_samples))
  train_data <- arrange(gene_data[gene_data$sample_id %in% train_sample_ids, ], sample_id, start)
  test_sample_ids <- setdiff(all_samples, train_sample_ids)
  test_data <- arrange(gene_data[gene_data$sample_id %in% test_sample_ids, ], sample_id, start)
  
  train_sequences <- pull(summarize(group_by(train_data, sample_id), n = n(), .groups = "drop"), n)
  test_sequences <- pull(summarize(group_by(test_data, sample_id), n = n(), .groups = "drop"), n)
  
  if (length(train_sequences) < 2 || length(test_sequences) < 2) {
    return(NULL)
  }
  
  train_hmm <- tryCatch({
    hmm_model <- depmix(
      list(beta_logit ~ 1, log_expression ~ beta_logit),
      data = train_data,
      nstates = 3,
      family = list(gaussian(), gaussian()),
      ntimes = train_sequences) 
    suppressWarnings(fit(hmm_model, verbose = FALSE))
  }, error = function(e) NULL)
  
  # convergence failure: if methyl-gex relationship not identified by states
  if (is.null(train_hmm)) {
    return(NULL)
  }
  
  learned_params <- getpars(train_hmm)
  train_data$state <- posterior(train_hmm, type = "viterbi")$state
  
  test_hmm <- tryCatch({
    hmm_model <- depmix(
      list(beta_logit ~ 1, log_expression ~ beta_logit),
      data = test_data,
      nstates = 3,
      family = list(gaussian(), gaussian()),
      ntimes = test_sequences) 
    setpars(hmm_model, learned_params) # use parameters learned from train data
  }, error = function(e) NULL)
  
  if (is.null(test_hmm)) {
    return(NULL)
  }
  
  test_data$state <- posterior(test_hmm, type = "viterbi")$state
  
  # proportion of CpG sites per sample in each state
  calculate_proportions <- function(df) {
    summarize(group_by(df, sample_id), log_expression = dplyr::first(log_expression), 
              prop_in_state_1 = mean (state == 1),
              prop_in_state_2 = mean (state == 2),
              prop_in_state_3 = mean (state == 3),
              .groups = "drop")
  }
  
  train_features <- calculate_proportions(train_data)
  test_features <- calculate_proportions(test_data)
  
  lm_fit <- tryCatch(
    lm(log_expression ~ prop_in_state_1 + prop_in_state_2 + prop_in_state_3,
       data = train_features), error = function(e) NULL
  )
  
  if (is.null(lm_fit)) {
    return(NULL)
  }
  
  predicted_expression <- tryCatch(predict(lm_fit, newdata = test_features), error = function(e) NULL)
  
  if (is.null(predicted_expression)) {
    return(NULL)
  } 
  
  valid_values <- !is.na(predicted_expression)
  
 if (sum(valid_values) < 5) {
    return(NULL)
  }
  
  r_squared <- as.numeric(cor(test_features$log_expression[valid_values], 
                              predicted_expression[valid_values])^2)
  
  list(summary = data.frame(
    gene_id_entrez = gene,
    num_cpgs = num_cpgs,
    n_samples = n_samples,
    r_squared = r_squared),
    lm_fit = lm_fit,
    hmm_fit = train_hmm)
  
}