library(shazam)
library(alakazam)
library(stringr)
library(dplyr)

# Including some functions
source("/home/ls/hadasn/new_pipeline/wrk_env/80-Elul-poptree/unite_program/treez/1_3/rscripts/functions.R")

# Defining file names
args <- commandArgs(TRUE)
#args <- 'C:\\Users\\user\\Documents\\biolab\\selection\\pat3_ex_1000.tab'

if (length(args) < 2) {
  stop("No input. Should be: selection.r <IgTreez_output1> <IgTreez_output2> ...", call. = FALSE)
}

# Collecting and updating the input databases
df_list <-  lapply(args, FUN = update_col_names)

# Testing the sample names
sample_names <- lapply(df_list, {function (x) x$sample})
if (length(sample_names) != length(unique(sample_names)))
  stop("The sample names of the input databases are not unique.", call. = FALSE)

# Binding the dataframes by rows
df <- do.call("rbind", df_list)

df$clone_id <- seq.int(nrow(df))

# Converting the int to numeric
df$mu_count_cdr_r <- as.numeric(df$mu_count_cdr_r)
df$mu_count_cdr_s <- as.numeric(df$mu_count_cdr_s)
df$mu_count_fwr_r <- as.numeric(df$mu_count_fwr_r)
df$mu_count_fwr_s <- as.numeric(df$mu_count_fwr_s)


# Creaing the output prefix
plot_prefix <- paste0(unique(as.character(df$sample)), collapse = "_")
if (nchar(plot_prefix) > 50){
  plot_prefix <- paste0(substr(plot_prefix, 1, 40),sample(1:900, 1))
  cat('***Output saved as ',plot_prefix, '\n')
  
}

# Count expected mutations and append to the output
expected <- expectedMutations(df,
                              sequenceColumn = "clonal_germline",
                              germlineColumn = "clonal_germline",
                              targetingModel = HH_S5F,
                              regionDefinition = IMGT_V)

# write.csv(expected, 'expected_no_cdr3.csv', quote = F)

if ('cdr3_end' %in% names(df)){
  # Calculate the expected row by row using the CDR3 length
  for (r in 1:nrow(df)){
    expected[r,] <- calc_expected_on_one_row(df[r,], sequenceColumn = "clonal_germline") 
  }
  
  # write.csv(expected, 'expected_with_cdr3.csv', quote = F)
  p_sum_name <- paste0(plot_prefix, '_selection_sum_with_cdr3.pdf', collapse = "_")
  p_den_name <- paste0(plot_prefix, '_selection_density_with_cdr3.pdf', collapse = "_")
  file_name <- paste0(plot_prefix, '_selection_pvals_with_cdr3.txt', collapse = "_")
  csv_name <- paste0(plot_prefix, '_selection_sigma_with_cdr3.csv', collapse = "_")
  
  
  
} else {
  p_sum_name <- paste0(plot_prefix, '_selection_sum_no_cdr3.pdf', collapse = "_")
  p_den_name <- paste0(plot_prefix, '_selection_density_no_cdr3.pdf', collapse = "_")
  file_name <- paste0(plot_prefix, '_selection_pvals_no_cdr3.txt', collapse = "_")
  csv_name <- paste0(plot_prefix, '_selection_sigma_with_cdr3.csv', collapse = "_")
}



# Calculate selection scores
baseline <- calcBaseline(expected, testStatistic = "focused",
                         sequenceColumn = "clonal_germline",
                         germlineColumn = "clonal_germline",
                         regionDefinition = IMGT_V)

# Combine selection scores by time-point
grouped <- groupBaseline(baseline, groupBy="sample")



# Plot mean and confidence interval by time-point
p <- plotBaselineSummary(grouped, "sample") + theme(text = element_text(size=20), axis.text.x = element_text(size = 12))
ggsave(p_sum_name, p)

# Plot selection PDFs for a subset of the data
p <- plotBaselineDensity(grouped, "sample", sigmaLimits=c(-1, 1)) + theme(text = element_text(size=20), legend.text = element_text(size = 10))
ggsave(p_den_name, p)

# Testing the difference in selection PDFs between groups
con <- file(file_name)
sink(con)
testBaseline(grouped, groupBy="sample")
write.csv(grouped@stats, csv_name)






