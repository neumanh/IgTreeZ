library(shazam)
library(alakazam)
library(stringr)
library(dplyr)

# Getting the source script path (from https://stackoverflow.com/a/1815743)
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(file_arg_name, "", initial_options[grep(file_arg_name, initial_options)])
script_basename <- dirname(script_name)
source_name <- paste0(script_basename, "/functions.r")
# Including some functions
source(source_name)

# Defining file names
args <- commandArgs(TRUE)

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


# Creating the output prefix
plot_prefix <- paste0(unique(as.character(df$sample)), collapse = "_")
if (nchar(plot_prefix) > 50){
  plot_prefix <- paste0(substr(plot_prefix, 1, 40),"_", sample(1:900, 1))
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

# Replacing "_" in space
sample_order = sapply(sample_names, levels)
sample_order_spaces = gsub('_', ' ', sample_order)
grouped@db$sample <- factor(gsub('_', ' ', grouped@db$sample), levels = sample_order_spaces)

# Plot selection PDFs for a subset of the data
p <- plotBaselineDensity(grouped, "sample", sigmaLimits=c(-1, 1), silent = T) + 
  theme(text = element_text(size=20), legend.text = element_text(size = 20), 
        legend.title = element_blank(), legend.key.size = unit(0.9, "cm"), 
        legend.title.align = 0.25) +
  xlab("Selection score")
  
ggsave(p_den_name, p, width=9)


# Plot mean and confidence interval by time-point

# Replacing the stat-sample names with " " and wrapping
#grouped@stats$sample <- gsub('_', ' ', grouped@stats$sample)
#grouped@stats$sample <- str_wrap(grouped@stats$sample, width = 10)


p <- plotBaselineSummary(grouped@stats, "sample", size = 1.15, silent = T) + 
  ylab("Selection score") +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 15, color = "black",
                                   angle = 0, hjust = 0.5))

ggsave(p_sum_name, p)

##### Adding numbers to the sample names
new_stats <- grouped@stats

i = 1
for (sn in unique(grouped@stats$sample)) {

  new_stats$sample <- gsub(paste0('\\b',sn,'\\b'), paste0(i, '_', sn), new_stats$sample) # The \b is used to replace the whole phrase
  i = i+1
}

p_sum_name <- paste0(plot_prefix, '_selection_sum_ordered_with_cdr3.pdf', collapse = "_")
p <- plotBaselineSummary(new_stats, "sample", size = 1.15, silent = T) + 
  ylab("Selection score") +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 15, color = "black",
                                   angle = 0, hjust = 0.5))
ggsave(p_sum_name, p)
############


# Testing the difference in selection PDFs between groups
con <- file(file_name)
sink(con)
testBaseline(grouped, groupBy="sample")
write.csv(grouped@stats, csv_name)


# Ordering the sample names by the input order
#grouped@stats <- grouped@stats[order(match(grouped@stats$sample, sample_order)),]



