###### HEADER ######################################################
# Hyphy - absrel
# Author: Maggie
# Date: 20250611
###### HEADER ######################################################

# 1. Setup ====

# https://hendrikvanb.gitlab.io/2018/07/nested_data-json_to_tibble/
# install.packages("jsonlite")
library(jsonlite)
library(tidyverse)
library(parallel)
library(treeio)

options(warn=-1) # turn off warning globally, turn back on: options(warn=0)
args <- commandArgs(trailingOnly = TRUE)

#args 1 is where the dir with absrel json output
#args 2 is file pattern to find
#args 3 is number of cores

project_path = args[1]
# project_path = "/Users/meng-ching.ko/ownCloud/Baldwin/Projects/Nectar/Seq_Data/globus/findgenes/json_nexus_nointernalnodes/ABSREL/"
# project_path = "/Users/meng-ching.ko/ownCloud/Baldwin/Projects/Nectar/Seq_Data/globus/test/"

pattern = args[2]
# pattern = ".json$"

core = args[3]
# core = 1

# 2. find json files ====
target_working_files = list.files(path = paste0(project_path), 
                                  pattern = pattern, 
                                  full.names = T, recursive = TRUE); target_working_files
# target_working_files[grepl("CDHR2", target_working_files, ignore.case = T)]

# 3. function ====
get_absrel_table = function(s){
  # s = target_working_files[2]
  print(s)
  temp_f_path = s;
  tmp = fromJSON(txt = temp_f_path) ; # glimpse(tmp)
  
  current_gene = gsub("\\.json", "", basename(temp_f_path))
  
  omega_highcutoff = 10
  
  tmp_aic = tmp$fits$`Full adaptive model`$`AIC-c`; 
  
  ### get tree
  tmp_tree = tmp$input$trees$`0`; # glimpse(tmp_tree)
  tree <- read.tree(text = paste0(tmp_tree, ";")) 
  nTips = length(tree$tip.label)
  
  tree = tree %>% as_tibble() %>% 
    bind_cols("gene" = current_gene, "nTips" = nTips, "AIC-c" = tmp_aic) 
  # tree_p = ggtree(tree) +
  #   geom_tiplab(offset = 0.2, size = FontSize) +
  #   geom_nodelab() +
  #   # geom_text(aes(label=node), hjust=-.3) +
  #   xlim(0, 20) +
  #   # xlim(0, max(tree[[1]])) +
  #   theme_tree() ; tree_p
  # tree_p = tree_p %>% rotate(51); tree_p
  
  ### get branch data (e.g. LRT)
  tmp_b = tmp$`branch attributes`$`0`; # glimpse(tmp_b)
  
  df = lapply(1:length(tmp_b), function(i){
    # tmp_b[[i]]
    tmp_df = tmp_b[[i]] %>% 
      map_if(is.matrix, function(s){data.frame(matrix(s, nrow = 1))}) %>% 
      map_if(is.data.frame, list) %>%
      as_tibble() %>% unnest(cols = c(`Rate Distributions`))
    return(tmp_df)
  })
  df_2 = df %>% bind_rows() %>% add_column("Name" = names(tmp_b)) %>% 
    select(-`original name`) 
  
  # handling different num. of omega classes 
  ind = grep("X\\d+", names(df_2))
  trouble = df_2[, ind] %>% as.matrix()
  trouble_sorted = lapply(1:nrow(trouble), function(row_i){
    x = trouble[row_i, ] 
    step1 = na.omit(x)
    n = length(step1)
    half = n/2
    tmp_omega = matrix(step1[1:half], nrow = 1) %>% as_tibble()
    names(tmp_omega) = paste0("Omega", 1:half)
    tmp_proportion = matrix(step1[(half+1):n], nrow = 1) %>% as_tibble()
    names(tmp_proportion) = paste0("Proportion", 1:half)
    high_omega = step1[half] %>% as_tibble()
    names(high_omega) = "HighestOmega"
    high_proportion = step1[n] %>% as_tibble()
    names(high_proportion) = "HigheshProportion"
    out = cbind("Name" = df_2$Name[row_i], tmp_omega, tmp_proportion, 
                high_omega, high_proportion) %>% as_tibble()
  }) %>% bind_rows()
  
  ## modify df_2 for more data variations
  df_3 = df_2 %>% select(-all_of(ind)) %>% left_join(trouble_sorted) %>%
    # adding Uncorrected_P.value<0.05 column
    mutate(`Uncorrected_P.value<0.05` = ifelse(`Uncorrected P-value` < 0.05, 
                                               `Uncorrected P-value`, NA )) %>%
    # make a high cutoff on highest omega
    mutate(HighestOmega_cutoff = ifelse(HighestOmega >= omega_highcutoff, NA, HighestOmega)) 
  
  df_4 = df_3 %>% left_join(tree, by = c("Name" = "label"))
  
  return(df_4)
}

get_absrel_aic = function(s){
  # s = target_working_files[2]
  print(s)
  temp_f_path = s;
  tmp = fromJSON(txt = temp_f_path) ; # glimpse(tmp)
  
  current_gene = gsub("\\.json", "", basename(temp_f_path))
  
  omega_highcutoff = 10
  
  tmp_aic = tmp$fits$`Full adaptive model`$`AIC-c`; 
  
  tmp_out = tibble("gene" = current_gene, "AIC-c" = tmp_aic) 
  
  return(tmp_out)
}

# output = lapply(target_working_files, get_absrel_table) %>% bind_rows()
output = mclapply(target_working_files,  get_absrel_table, mc.cores= core, mc.preschedule = FALSE) %>%
  bind_rows()

out_dir = paste0(project_path, "/../combined/")
dir.create(out_dir)

path = paste0(out_dir, "/combined_", pattern)
write_tsv(out, path, na = "")

print(paste0(path, " done!"))
