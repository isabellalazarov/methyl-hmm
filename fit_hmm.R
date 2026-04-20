fit_hmm <- function(gene, hmm_data, nstates) {
  gene_data <- hmm_data[hmm_data$gene_id_entrez == gene, ]
  
  # at least 10 CpG sites for better model fit
  if (nrow(gene_data) < 10)
    return(NULL)
  
  #hmm set-up: each state emits methylation value according to own gaussian
  hmm_model <- depmix(
    mean_methylation ~ 1,
    data = gene_data,
    nstates = 3,
    family = gaussian(),
    ntimes = nrow(gene_data))
  
  # baum-welch
  fit_model <- tryCatch(
    {capture.output(result <- suppressWarnings(fit(hmm_model, verbose = FALSE))); result},
    error = function(e) NULL)
  if (is.null(fit_model)) 
    return(NULL)
  
  # viterbi
  # prop_# = proportion of CpG sites that belong to each state
  post <- posterior(fit_model, type = "viterbi")
  data.frame(
    gene_id_entrez = gene,
    prop_state1 = mean(post$state == 1),
    prop_state2 = mean(post$state == 2),
    prop_state3 = mean(post$state == 3),
    log_expression = gene_data$log_expression[1])
}