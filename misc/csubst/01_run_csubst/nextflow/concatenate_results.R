library(tidyverse)
library(parallel)
options(warn=-1) # turn off warning globally, turn back on: options(warn=0)

args <- commandArgs(trailingOnly = TRUE)
#args 1 is where the dir with csubst output
#args 2 is file pattern to find
#args 3 is number of cores

project_path = args[1]
# project_path = "/home/mpg08/mko/Nectar/analysis/csubst/01_run_csubst/nextflow/results/csubst_out/" 
# project_path = "/home/mpg08/mko/Nectar/analysis/csubst/01_run_csubst_2024_proto/nf_chunk/results/csubst_out"

pattern = args[2]
# pattern = "csubst_cb_stats.tsv"

core = args[3]
# core = 1

concat_func = function(path){
  # tmp_path = target_files[1]
  tmp_path = path 
  tmp_dir = dirname(tmp_path)
  tmp = unlist(strsplit(gsub(project_path, "", tmp_dir), "/"))
  tmp_transcript = tmp[2]
  tb = read_tsv(tmp_path, show_col_types = FALSE) 
  tb2 = bind_cols(transcript = tmp_transcript) %>% bind_cols(tb)
}

# 2. csubst_cb_2.tsv ====
target_files = list.files(project_path, pattern, full.names = TRUE, recursive = TRUE); # target_files

out = mclapply(target_files,  concat_func, mc.cores= core, mc.preschedule = FALSE) %>%
  bind_rows()

out_dir = paste0(project_path, "/../combined/")
dir.create(out_dir)

path = paste0(out_dir, "/combined_", pattern)
write_tsv(out, path, na = "")

print(paste0(path, " done!"))