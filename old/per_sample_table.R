# table: one row per sample, columns are log_expression, proportion of sample
# that is in state 1, proportion of sample in state 2, proportion of sample in 
# state 3
per_sample_table <- function(gene_data) {
  
  summarize(group_by(gene_data, sample_id), log_expression = first(log_expression),
            prop_in_state_1 = mean(state == 1),
            prop_in_state_2 = mean(state == 2),
            prop_in_state_3 = mean(state == 3),
            .groups = "drop"
  )
  
}