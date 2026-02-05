library(tidyverse)

# Arguments parsing
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 2){
  stop("Call 'qPCR_analyzer' with the following arguments:
# arg1: directory containing qPCR xlsx result files
# arg2: name of the housekeeping gene"
       , call.=FALSE)
}

input_dir = args[1]
housekeeping_gene = args[2]

print(paste("Using input files from:", input_dir))
print(paste("Using housekeeping gene:", housekeeping_gene))

run_qPCR_analyzer <- function(input_dir, housekeeping_gene){
  
  input_files <- list.files(input_dir, full.names = T)
  
  # Declare housekeeping gene
  hk <- housekeeping_gene
  
  # Read input files
  file_list <- list()
  
  for(i in 1:length(input_files)){
    file_extension <- tools::file_ext(input_files[i])
    
    if(file_extension == "xls"){
      file_list[[i]] <- readxl::read_xls(path = input_files[i], sheet = "Results")
    }else if(file_extension == "xlsx"){
      file_list[[i]] <- readxl::read_xlsx(path = input_files[i], sheet = "Results")
    }
  }

  # Remove headers and endings
  removeHeader <- function(input_file){
    starting_line <- which(input_file[,1] == "Well")
    input_file <- input_file[starting_line:nrow(input_file),]
    colnames(input_file) <- input_file[1,]
    input_file <- input_file[2:nrow(input_file),]
    return(input_file)
  }
  processed_files <- lapply(X = file_list, FUN = removeHeader)
  
  removeEnd <- function(input_file){
    end_line <- which(input_file[,"Well"] == 96)
    input_file <- input_file[1:end_line,]
    return(input_file)
  }
  processed_files <- lapply(X = processed_files, FUN = removeEnd)
  
  # Rbind all dataframes and keep useful columns
  master_df <- bind_rows(processed_files)
  master_df <- master_df %>% 
    dplyr::select(Well, `Well Position`, `Sample Name`, `Target Name`, CT, `Ct SD`) %>% 
    dplyr::mutate_at(c("CT", "Ct SD"), as.numeric)
  
  # Store warning messages for samples with a CT SD > 0.5
  warning_list <- c()
  for(i in 1:nrow(master_df)){
    if(!is.na(master_df$`Ct SD`[i]) & master_df$`Ct SD`[i] > 0.5){
      print(i)
      warning_list <- c(warning_list, paste("Warning, SD higher than 0.5 detected for sample", master_df$`Sample Name`[i], "in gene", master_df$`Target Name`[i], "\n"))
      warning_list <- unique(warning_list)
    }
  }
  
  # Print warnings
  for(i in 1:length(warning_list)){
    writeLines(warning_list[i])
  }
  
  # Create df with average CTs
  avg_cts <- master_df %>% 
    dplyr::group_by(`Sample Name`, `Target Name`) %>% 
    summarise(avg_cts = mean(CT, na.rm = TRUE)) %>% 
    pivot_wider(names_from = `Target Name`, values_from = avg_cts)
  
  # Create df with delta CTs 
  delta_cts <- master_df %>% 
    dplyr::group_by(`Sample Name`, `Target Name`) %>% 
    summarise(avg_cts = mean(CT, na.rm = TRUE))
  
  delta_cts$delta_cts <- rep(NA, nrow(delta_cts))
  
  for(i in 1:nrow(delta_cts)){
    
    sample_name <- delta_cts$`Sample Name`[i]
    
    delta_cts$delta_cts[i] <- delta_cts$avg_cts[i] - (
      delta_cts %>% dplyr::filter(`Target Name` == hk & `Sample Name` == sample_name) %>% pull(avg_cts)
    )
  }
  
  delta_cts <- delta_cts %>%
    dplyr::select(-avg_cts) %>% 
    pivot_wider(names_from = `Target Name`, values_from = delta_cts) %>% 
    dplyr::select(-!!sym(hk)) 
  
  # Create df with: 2^-(delta ct)
  
  levels_df <- master_df %>% 
    dplyr::group_by(`Sample Name`, `Target Name`) %>% 
    summarise(avg_cts = mean(CT, na.rm = TRUE))
  
  levels_df$delta_cts <- rep(NA, nrow(levels_df))
  
  for(i in 1:nrow(levels_df)){
    
    sample_name <- levels_df$`Sample Name`[i]
    
    levels_df$delta_cts[i] <- levels_df$avg_cts[i] - (
      levels_df %>% dplyr::filter(`Target Name` == hk & `Sample Name` == sample_name) %>% pull(avg_cts)
    )
  }
  
  levels_df <- levels_df %>%
    dplyr::mutate(level = 2^(-delta_cts)) %>% 
    dplyr::select(-c("avg_cts", "delta_cts")) %>% 
    pivot_wider(names_from = `Target Name`, values_from = level) %>% 
    dplyr::select(-!!sym(hk))
  
  final_results <- list(master_df, avg_cts, delta_cts, levels_df)
  names(final_results) <- c("CTs", "Avg_CTs", "Delta_CTs", "Levels")
  
  return(final_results)
}

res <- run_qPCR_analyzer(input_dir = input_dir, housekeeping_gene = housekeeping_gene)
writexl::write_xlsx(res, gsub(pattern = "/", replacement = "", x = paste0(input_dir, ".xlsx")))
