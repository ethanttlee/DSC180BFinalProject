library(coloc)

# Set the disease information
disease <- "schizophrenia"
disease_file <- "gwas_data/schizophrenia.tsv"
sample_size <- 282871
case_control_prop <- 0.0019
eQTL_files <- list.files(path="eQTL_subsets", pattern='ENSG.*', full.names = TRUE)

error_file_path <- "coloc_error_log.txt"
result_matrix <- list()

# Read the GWAS data for the first disease
# curr_gwas <- na.omit(read.csv(paste("gwas_data/", disease, ".tsv", sep = ""), 
curr_gwas <- na.omit(read.csv(paste("GCST010774_buildGRCh37.tsv", sep = ""), 
                              sep = '\t', header = TRUE,
                              colClasses = c("character", "double", "NULL", "double", "NULL",
                                             "NULL", "double", "NULL", "NULL", "NULL",
                                             "NULL", "NULL", "NULL", "NULL", "NULL",
                                             "NULL", "NULL", "NULL", "NULL")))

# Remove duplicate SNPs
curr_gwas <- curr_gwas[!duplicated(curr_gwas[, "variant_id"]), ]

# Define the genome-wide significance threshold
# p_value_threshold <- 5e-2
p_value_threshold <- 5e-2

# Filter the dataframe to include only SNPs with p-values below the threshold
filtered_gwas <- curr_gwas[curr_gwas$p_value < p_value_threshold, ]

# Prepare the GWAS list
gwas_list <- list()
gwas_list$MAF <- filtered_gwas$MAF_calculated_from_dosage_data
gwas_list$snp <- filtered_gwas$variant_id
gwas_list$position <- filtered_gwas$base_pair_location
gwas_list$N <- sample_size
gwas_list$pvalues <- filtered_gwas$p_value
gwas_list$type <- "cc"
gwas_list$s <- case_control_prop

# Initialize an empty data frame to store results with the appropriate columns
results_df <- data.frame(
  gene_name = character(),
  nsnps = numeric(),
  PP_H0_abf = numeric(),
  PP_H1_abf = numeric(),
  PP_H2_abf = numeric(),
  PP_H3_abf = numeric(),
  PP_H4_abf = numeric(),
  stringsAsFactors = FALSE
)

batch_start <- 1
batch_end <- 21

# Check the GWAS dataset
gwas_fail <- tryCatch({
    check_dataset(gwas_list)
}, error = function(cond) {
    error_message <- paste("Error: disease file ", disease, " has error: \n ", cond, sep = "")
    cat(error_message, file = error_file_path, append = TRUE)
})

if (inherits(gwas_fail, "error")) {
    print("Error in GWAS dataset. Skipping further processing.")
} else {
#     print(paste("Processed disease file ", disease_file))
    
    # Initialize a list to store results
    curr_results <- list()
    
    # Iterate over eQTL files
    for (i in batch_start:batch_end) {
        if(i > length(eQTL_files)) {
            break 
          }
        
        eQTL_file <- eQTL_files[i]
        gene_name <- gsub(".csv$", "", basename(eQTL_file))  # Extract gene name from file name
        eQTL_data <- read.csv(eQTL_file)

        eQTL_list <- list()
        eQTL_list$beta <- eQTL_data$beta
        eQTL_list$varbeta <- eQTL_data$varbeta
        eQTL_list$snp <- eQTL_data$snp
        eQTL_list$position <- eQTL_data$pos
        eQTL_list$type <- eQTL_data[, 'type'][1]
        eQTL_list$N <- eQTL_data$N[0]
        eQTL_list$MAF <- eQTL_data$MAF

        # Check the eQTL dataset
        eQTL_fail <- tryCatch({
            suppressWarnings(check_dataset(eQTL_list))
        }, error = function(cond) {
            error_message <- paste("Error: eQTL file ", eQTL_file, " has error: \n ", cond$message, sep = "")
            cat(error_message, file = error_file_path, append = TRUE)
            return(NULL)  # Return NULL to indicate failure
        })

        if (inherits(eQTL_fail, "error")) {
            # Skip further processing if there was an issue with the dataset
            next
        }
        
        # Check for missing values in the "type" column
        if (any(is.na(eQTL_data$type))) {
            cat("Skipping processing for", eQTL_file, "due to missing values in 'type' column.\n")
            next  # Move to the next iteration
        }
        
        #THIS IS WHERE I FILTER THE DATASETS
        common_snps <- intersect(gwas_list$snp, eQTL_list$snp)
#         print(length(common_snps))
#         print(length(eQTL_list$snp))

        # Filter both gwas_list and eQTL_list to only include these common SNPs
        # And ensure they are in the same order

        # For gwas_list - Assuming it's already filtered globally, just reorder based on common SNPs
        gwas_list_filtered <- list(
          MAF = gwas_list$MAF[match(common_snps, gwas_list$snp)],
          snp = common_snps,  # This ensures the order matches
          position = gwas_list$position[match(common_snps, gwas_list$snp)],
          N = gwas_list$N,
          pvalues = gwas_list$pvalues[match(common_snps, gwas_list$snp)],
          type = gwas_list$type,
          s = gwas_list$s
        )

        # For eQTL_list - this part goes inside your loop where you prepare eQTL_list
        eQTL_list_filtered <- list(
          beta = eQTL_list$beta[match(common_snps, eQTL_list$snp)],
          varbeta = eQTL_list$varbeta[match(common_snps, eQTL_list$snp)],
          snp = common_snps,  # This ensures the order matches
          position = eQTL_list$position[match(common_snps, eQTL_list$snp)],
          type = eQTL_list$type,
          N = eQTL_list$N,
          MAF = eQTL_list$MAF[match(common_snps, eQTL_list$snp)]
        )
        
        # Run coloc.abf inside a tryCatch block
        coloc_results <- tryCatch({
#             coloc.abf(gwas_list, eQTL_list)
            coloc.abf(gwas_list_filtered, eQTL_list_filtered)
        }, error = function(err) {
            # Print an error message
            cat("Error in coloc.abf:", conditionMessage(err), "\n")
            # Return NULL to indicate failure
            return(NULL)
        })

        # Check if coloc_results is NULL (indicating an error occurred)
        if (is.null(coloc_results)) {
            # Skip further processing if there was an issue with coloc.abf
            next
        }
        
        temp_results <- data.frame(
          gene_name = gene_name,
          nsnps = coloc_results$summary["nsnps"],
          PP_H0_abf = coloc_results$summary["PP.H0.abf"],
          PP_H1_abf = coloc_results$summary["PP.H1.abf"],
          PP_H2_abf = coloc_results$summary["PP.H2.abf"],
          PP_H3_abf = coloc_results$summary["PP.H3.abf"],
          PP_H4_abf = coloc_results$summary["PP.H4.abf"],
          stringsAsFactors = FALSE
        )
        
        # Append temp_results to results_df
        results_df <- rbind(results_df, temp_results)
        
    }
}

# Ensure the dplyr package is installed and loaded
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)

# Use distinct() to keep only unique rows based on gene_name
results_df_unique <- results_df %>% distinct(gene_name, .keep_all = TRUE)

# Specify the file path
file_path <- "results_df_schizophrenia.csv"

# Check if the file already exists to determine if the header should be included
if(file.exists(file_path)) {
  # File exists, append without header
  # Note: write.table is used here with sep="," to produce CSV formatted output
  write.table(results_df_unique, file = file_path, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE, quote = TRUE)
} else {
  # File does not exist, include header and do not append (implicitly creates a new file)
  write.table(results_df_unique, file = file_path, sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE, quote = TRUE)
}

