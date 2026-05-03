fit_hmm_per_gene <- function(gene, hmm_data, num_states, sample_ids = NULL) {
  num_states <- as.integer(num_states)
  if (length(num_states) != 1 || is.na(num_states) || num_states < 2) {
    stop("nstates must be a single integer >= 2")
  }

  gene_data <- hmm_data[hmm_data$gene_id_entrez == gene, ]
  #gene_data <- gene_data[order(gene_data$start), ]
  
  if (!is.null(sample_ids)) {
    gene_data <- gene_data[gene_data$sample_id %in% sample_ids, ]
  }
  
  gene_data <- gene_data[order(gene_data$sample_id, gene_data$start), ]
  
  # at least 10 CpG sites for better model fit
  if (nrow(gene_data) < 10)
    return(NULL)
  
  # hmm set-up - each state emits methylation and a state-specific methylation to expression regression
  gene_data$beta_logit <- qlogis(pmin(pmax(gene_data$beta, 1e-4), 1 - 1e-4))  # logit transform betas to gaussian, avoid issues with 0 and 1
  gene_data <- gene_data[complete.cases(gene_data[, c("beta_logit", "log_expression")]), ]  # remove missing values

  if (nrow(gene_data) < 10) {
    return(NULL)
  }
  
  sample_sequences <- pull(summarize(group_by(gene_data, sample_id), n = n(), .groups = "drop"), n)
  
  if (length(sample_sequences) < 10) {
    return(NULL)
  }

  hmm_model <- depmix(
    list(beta_logit ~ 1, log_expression ~ beta_logit),
    data = gene_data,
    nstates = num_states,
    family = list(gaussian(), gaussian()),
    ntimes = sample_sequnces)
  
  # baum-welch
  fit_model <- tryCatch(
    {capture.output(result <- suppressWarnings(fit(hmm_model, verbose = FALSE))); result},
    error = function(e) NULL)
  if (is.null(fit_model)) 
    return(NULL)
  
  # viterbi
  post <- posterior(fit_model, type = "viterbi")

  state_params <- tryCatch(
    bind_rows(lapply(seq_len(num_states), function(state_id) {
      methyl_resp <- getmodel(fit_model, which = "response", state = state_id, number = 1)
      expr_resp <- getmodel(fit_model, which = "response", state = state_id, number = 2)

      methyl_pars <- getpars(methyl_resp)
      expr_pars <- getpars(expr_resp)

      data.frame(
        state = state_id,
        state_prop = mean(post$state == state_id),
        methyl_intercept = unname(methyl_pars[1]),
        methyl_sd = unname(methyl_pars[length(methyl_pars)]),
        expr_intercept = unname(expr_pars[1]),
        expr_slope = unname(expr_pars[2]),
        expr_sd = unname(expr_pars[length(expr_pars)]),
        stringsAsFactors = FALSE
      )
    })),
    error = function(e) NULL
  )

  if (is.null(state_params))
    return(NULL)

  # sort/relabel states by methyl_intercept so state 1 = lowest methylation, state 2 = medium, state 3 = highest
  state_order_map <- state_params %>%
    distinct(state, methyl_intercept) %>%
    arrange(methyl_intercept) %>%
    mutate(ordered_state = row_number()) %>%
    select(state, ordered_state)

  state_params <- state_params %>%
    left_join(state_order_map, by = "state") %>%
    mutate(state = ordered_state) %>%
    select(-ordered_state)
  
  gene_data$state <- post$state
  gene_data <- select(mutate(left_join(gene_data, state_order_map, by = "state"), state = ordered_state), -ordered_state)
  state_params$gene_id_entrez <- gene
  
  list(
    gene_data = gene_data,
    state_params = state_params
  )
}