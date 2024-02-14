#install.packages("coloc")
library(coloc)

error_file_path <- "coloc_error_log.txt"

disease_list <- c("essential_hypertension", "abdominal_hernia", "hyperlipidemia", "osteoarthrosis", "cardiac_dysrhythmias", "asthma", "cataract", "coronary_atherosclerosis", "type_2_diabetes", "parkinsons_disease", "alzheimers_disease", "schizophrenia")
# from https://www.ebi.ac.uk/gwas/publications/32589924
disease_files <- c("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010774/GCST010774_buildGRCh37.tsv",
                   "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010773/GCST010773_buildGRCh37.tsv",
                   "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010772/GCST010772_buildGRCh37.tsv",
                   "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010771/GCST010771_buildGRCh37.tsv",
                   "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010770/GCST010770_buildGRCh37.tsv",
                   "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010769/GCST010769_buildGRCh37.tsv",
                   "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010768/GCST010768_buildGRCh37.tsv",
                   "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010767/GCST010767_buildGRCh37.tsv",
                   "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010766/GCST010766_buildGRCh37.tsv",
                   "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010765/GCST010765_buildGRCh37.tsv",
                   "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010764/GCST010764_buildGRCh37.tsv",
                   "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010763/GCST010763_buildGRCh37.tsv")
sample_sizes <- rep(282871, 12)
case_control_props <- c(0.271, 0.0154, 0.1259, 0.1028, 0.0904, 0.0085, 0.0800, 0.0674, 0.0656, 0.0048, 0.0023, 0.0019)

# from 1000G gene expression data
eQTL_files <- list.files(path="eQTL_subsets", pattern='ENSG.*', full.names = TRUE)

result_matrix <- list()

for (i in 1:length(disease_files)) {
    print(disease_list[i])
    
    curr_results <- list()
    
#    #download current gwas file
#    system2(command = "curl", 
#        args    = c(disease_files[i]), 
#        stdout  = paste("gwas_data/",disease_list[i],".tsv", sep="", collapse=""))
    
    curr_gwas <- na.omit(read.csv(paste("gwas_data/",disease_list[i],".tsv", sep="", collapse=""), 
                                  sep = '\t', header = TRUE,
                                  colClasses = c("character", "double", "NULL", "double", "NULL",
                                            "NULL", "double", "NULL", "NULL", "NULL",
                                            "NULL", "NULL", "NULL", "NULL", "NULL",
                                        "NULL", "NULL", "NULL", "NULL")))
    #remove duplicate snps
    curr_gwas <- curr_gwas[!duplicated(curr_gwas[,c("variant_id")]),]
    
    gwas_list <- list()
    gwas_list$MAF <- curr_gwas$MAF_calculated_from_dosage_data
    gwas_list$snp <- curr_gwas$variant_id
    gwas_list$position <- curr_gwas$base_pair_location
    gwas_list$N <- sample_sizes[i]
    gwas_list$pvalues <- curr_gwas$p_value
    gwas_list$type <- "cc"
    gwas_list$s <- case_control_props[i]
    
    tryCatch( 
	{
		check_dataset(gwas_list)
	},
	error = function(cond)	{
        	error_message <- paste("Error: disease file ", disease_list[i], " has error: \n ", cond, sep="", collapse="")
        	cat(error_message, file = error_file_path, append = TRUE)
		next
    	}
	)

    print(paste("Processed disease file ", disease_files[i], sep="", collapse=""))
    
    
    num_files <- length(eQTL_files)
    counter <- 0
    
    for (eQTL_file in eQTL_files) {
        eQTL_data <- read.csv(eQTL_file)
        
        eQTL_list <- list()
        eQTL_list$beta <- eQTL_data$beta
        eQTL_list$varbeta <- eQTL_data$varbeta
        eQTL_list$snp <- eQTL_data$snp
        eQTL_list$position <- eQTL_data$pos
        eQTL_list$type <- eQTL_data[, 'type'][1]
        eQTL_list$N <- eQTL_data$N[0]
        eQTL_list$MAF <- eQTL_data$MAF

	tryCatch( 
	{
		suppressWarnings(check_dataset(eQTL_list))
	},
	error = function(cond) {
		error_message <- paste("Error: eQTL file ", eQTL_file, " has error: \n ", check_dataset(gwas_list), sep="", collapse="")
                cat(error_message, file = error_file_path, append = TRUE)
                curr_results <- append(curr_results, list(NULL))
	},
	
	)

        
        coloc_results <- coloc.abf(gwas_list, eQTL_list)
    
        curr_results <- append(curr_results, coloc_results$summary[[6]])

        counter <- counter + 1
        
        print(paste("Processed ", counter/num_files, " eQTL files", sep="", collapse=""))
    
    }
    
    result_matrix[[disease_list[[i]]]] <- curr_results
    
    #remove file for storage space
    #system2(command = "rm", 
    #    args    = c(disease_list[i]+".tsv"), 
    #    stdout  = paste("gwas_data/",disease_list[i],".tsv", sep="", collapse=""))
    break 
}

result_frame <- t(data.frame(unlist(result_matrix)))
colnames(result_frame) <- list.files(path="/Users/johndriscoll/Downloads/180B/DSC180BFinalProject/eQTL_subsets", pattern='ENSG.*', full.names = FALSE)
rownames(result_frame) <- disease_list[1]

write.csv(result_frame, "coloc_matrix.csv")

