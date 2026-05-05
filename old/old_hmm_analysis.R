# # train hmm
# state_exp <- data.frame()
# source(here("fit_hmm.R"))
# i <- 0
# for (gene in train_genes) {
#   i <- i + 1
#   if (gene %in% state_exp$gene_id_entrez) {
#       print(paste0("[", i, "/", length(train_genes), "] ", "Skipping ", gene))
#       next
#   }
#   print(paste0("[", i, "/", length(train_genes), "] ", "Processing ", gene))
#   result <- fit_hmm(gene, hmm_data, num_states = num_states)
#   if(!is.null(result)) 
#     state_exp <- rbind(state_exp, result)
# }
# 
# # summarize state-specific regression coefficients to gene level
# train_state_summary <- state_exp %>%
#   group_by(gene_id_entrez) %>%
#   summarize(
#     log_expression = dplyr::first(log_expression),
#     weighted_methyl_intercept = sum(state_prop * methyl_intercept, na.rm = TRUE),
#     weighted_expr_intercept = sum(state_prop * expr_intercept, na.rm = TRUE),
#     weighted_expr_slope = sum(state_prop * expr_slope, na.rm = TRUE),
#     dominant_state = state[which.max(state_prop)][1],
#     .groups = "drop"
#   )
# 
# train_lm <- lm(log_expression ~ weighted_expr_intercept + weighted_expr_slope,
#                data = train_state_summary)
# summary(train_lm)
# 
# # predict test gene expression
# test_pred <- data.frame()
# test_state_exp <- data.frame()
# i <- 0
# for (gene in test_genes) {
#   i <- i + 1
#   if (gene %in% test_state_exp$gene_id_entrez) {
#       print(paste0("[", i, "/", length(test_genes), "] ", "Skipping ", gene))
#       next
#   }
#   print(paste0("[", i, "/", length(test_genes), "] ", "Processing ", gene))
# 
#   result <- fit_hmm(gene, hmm_data, num_states = num_states)
#   if (is.null(result))
#     next
#   
#   test_state_exp <- rbind(test_state_exp, result)
# }
# 
# test_state_summary <- test_state_exp %>%
#   group_by(gene_id_entrez) %>%
#   summarize(
#     log_expression = dplyr::first(log_expression),
#     weighted_methyl_intercept = sum(state_prop * methyl_intercept, na.rm = TRUE),
#     weighted_expr_intercept = sum(state_prop * expr_intercept, na.rm = TRUE),
#     weighted_expr_slope = sum(state_prop * expr_slope, na.rm = TRUE),
#     dominant_state = state[which.max(state_prop)][1],
#     .groups = "drop"
#   )
# 
# pred_exp <- predict(train_lm, newdata = test_state_summary)
# test_pred <- data.frame(actual = test_state_summary$log_expression,
#                         pred = pred_exp)
#   
# hmm_r_squared <- (cor(test_pred$actual, test_pred$pred))^2

