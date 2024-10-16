###### HEADER ######################################################
# Hyphy - absrel
# Author: Maggie
# Date: 20210328
# original file name: 5_1_Hyphy_parse_0706.R > Hyphy_parse_absrel_json_1208_glenn.R

# major goal for this script:
# 1. read result from json files
# 2. plot

# log:
# Significantly changed compared to (2020)1208 version
# 1. handling cases when there are more than 2 omega classes
# 2. ditch phylobase for treeio for easier appending data to tree
# 3. omega color ver5
# https://github.com/veg/hyphy/issues/1198  "aBSREL baseline == free-ratio models, except MG94xREVxCF3x4 vs GY94 (in PAML)"

###### HEADER ######################################################

# 1. Setup ====

# https://hendrikvanb.gitlab.io/2018/07/nested_data-json_to_tibble/
# install.packages("jsonlite")
library(jsonlite)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(ggtree)
library(treeio)
library(readxl)

FontSize = 2.5
AxisTxFontSizeSize = 5
AxisTxFontSizeSize_s = 4
AxisTitleFontSizeSize = 6
Factor_mmtoin = 0.0393701
Width_HalfCol = 98*Factor_mmtoin 

openfile <- function(filepath){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  type = tolower(os);
  switch(os, Linux = system(paste0("xdg-open ", filepath)),
         Windows = system(paste0("open \"", filepath, "\"")),
         osx = system(paste0("open \"", filepath, "\""))
  )}

theme_m = theme(panel.background = element_blank(),
                plot.background = element_blank(), # defualt is white
                plot.title = element_text(hjust = 0.5, size = AxisTitleFontSizeSize),
                panel.border = element_blank(), 
                axis.line = element_line(color = "black", linewidth = 0.1),
                panel.grid = element_blank(), 
                axis.ticks = element_line(colour = "black", linewidth = 0.05),
                legend.position = "none",
                # legend.position = c(0.5, 0.9),
                legend.background = element_blank())

# 2. Load files ====
project_path = paste0(getwd(), "/PhyloDiet/fig3/absrel/20221018/")

path = paste0(getwd(), "/FromMaude/20221012/AbsrelForMaggieOct11_3.xlsx") 
genelist = read_excel(path, sheet = 1) %>% filter(!is.na(Gene))

# 3. aBSREL ####
target_ori_files = "\\.json"
# target_ori_files = "_mod.json"

target_working_files = list.files(path = paste0(project_path), 
                                  pattern = target_ori_files, 
                                  full.names = T, recursive = TRUE); target_working_files
# target_working_files[grepl("CDHR2", target_working_files, ignore.case = T)]

# 3.1 plot function ====

plot_absrel_variHeatMapVal = function(s, withBranchLen, tmp_heatmap_val, x_limit, tipLabelSize){
  # s = target_working_files[7]
  print(s)
  temp_f_path = s;
  tmp = fromJSON(txt = temp_f_path) ; # glimpse(tmp)
  
  omega_highcutoff = 10
  
  ### get tree
  tmp_tree = tmp$input$trees$`0`; # glimpse(tmp_tree)
  tree <- read.tree(text = paste0(tmp_tree, ";"));
  # tree_p = ggtree(tree) +
  #   geom_tiplab(offset = 1, size = FontSize) +
  #   geom_nodelab(color = 'red', size = 2) +
  #   geom_text(aes(label=node), hjust=-.3) +
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
  
  print(paste0("Max Rate class: ", max(df_3$`Rate classes`)))
  
  ## append data using treeio
  tree_phy = tree %>% as_tibble() %>% left_join(df_3, by = c("label" = "Name")) %>%
    as.treedata()
  
  # # modify tip labels
  # tip_lab = tree_phy@phylo$tip.label
  # tip_lab = gsub("_R[1|3]", "", tip_lab)
  # 
  # tip_lab_new = sapply(tip_lab, function(t){
  #   if(sum(grepl(t, speciesNAME$ID, ignore.case = TRUE)) == 1) {
  #     ind = grep(t, speciesNAME$ID, ignore.case = TRUE)
  #     # speciesNAME$`Common Name`[ind]
  #     targetName[ind]
  #   }else{t}
  # })
  # # rbind(tip_lab_new, tip_lab)
  # tree_phy@phylo$tip.label = tip_lab_new
  
  ###  plot
  # tmp_heatmap_val = c("LRT",
  #                     "Uncorrected P-value",
  #                     "Uncorrected_P.value<0.05", 
  #                     "Corrected P-value",
  #                     "HighestOmega", "HighestOmega_cutoff) # Omega
  # tmp_heatmap_val = "Uncorrected_P.value<0.05"
  # tmp_heatmap_val = "HighestOmega_cutoff"
  # withBranchLen = 'n'
  # tt = tmp_heatmap_val
  
  p = lapply(tmp_heatmap_val, function(tt){
    
    ## add branch length
    if (withBranchLen == 'y'){
      tmp_branchlen = df_3$`Full adaptive model (synonymous subs/site)`
      tree_phy@phylo$edge.length = tmp_branchlen
      max_branchlength = max(tree_phy@phylo$edge.length , na.rm = T)
      # max_branchlength = sum(tree_phy@phylo$edge.length , na.rm = T)
      x_limit = max_branchlength * 3
      time_scale_x = max_branchlength * 3.1
      rootedge = max_branchlength
    } else {
      # max_branchlength = max(edgeLen, na.rm = T)
      # x_limit = max_branchlength*30
      rootedge = 2
    }
    
    
    ## theme and panel title
    title_ana = unlist(strsplit(basename(temp_f_path), "\\."))[1]
    # set = basename(dirname(temp_f_path)) ## only for chris, webserver
    # title = paste0(title_ana, " ", set, "\n", tt) # only for chris, webserver
    title = paste0(title_ana, "\n", tt)
    theme_tmp = theme_m + theme(legend.position = c(0.9, 0.5), # c(0.9, 0.2), 
                                legend.title = element_text(size = AxisTxFontSizeSize),
                                legend.text = element_text(size = AxisTxFontSizeSize_s),
                                legend.key.size = unit(6, 'pt'),
                                axis.text = element_text(size = AxisTxFontSizeSize),
                                plot.title = element_text(hjust = 0.5, size = AxisTitleFontSizeSize))
    
    ## color scheme depending on wanted value types
    if (grepl("value", tt)){
      if( grepl("Uncorrected_P.value<0.05", tt)){
        title = gsub("Uncorrected_P.value<0.05", "Uncorrected_P.value < 0.05", title)}
      legend_title = "p value"
      color_tmp = 'red'
      color_highp_tmp = 'rosybrown1'
      
    }else if(grepl("LRT", tt)){
      legend_title = "LRT"
      color_tmp = 'steelblue'
      color_lowval_tmp = 'skyblue1'
      
    }else if(grepl("Omega", tt)){
      legend_title = expression(omega)
      # ver5
      colpal = brewer.pal(name = 'Greens', n=9)
      color_highval_tmp = colpal[8] # "#00aa44" # 'steelblue'
      color_lowval_tmp = "lightgrey" #   # 'palegreen'
      color_midval_tmp = "khaki1" # colpal[5] 
      missingcol = "grey10" 
      
    }else{print('which tmp_heatmap_val ?')}
    
    # tipLabelSize = 1.5
    # tipLabelSize = 2
    linewidth = 0.5  
    offset_tip = x_limit/50
    
    if(grepl("value", tt) == FALSE){
      if (grepl("LRT", tt) ){
        
        #LRT 
        tree_p1 = ggtree(tree_phy, size = linewidth, aes_(col = as.name(tt))) + 
          geom_tiplab(offset = offset_tip, size = tipLabelSize, color = "grey10") +
          xlim_tree(x_limit) +
          scale_color_continuous(low = color_lowval_tmp, high = color_tmp, na.value = "grey50") +
          guides(color = guide_colorbar(title = legend_title)) +
          ggtitle( title ) +
          theme_tmp ;  tree_p1
        # tree_p1 = tree_p1 %>% rotate(25) ; tree_p1
      } else {
        
        #omega
        col = colorRampPalette(c(color_lowval_tmp, color_midval_tmp, color_highval_tmp))(10)
        tree_p1 = ggtree(tree_phy, size = linewidth, aes_(col = as.name(tt))) + 
          geom_tiplab(offset = offset_tip, size = tipLabelSize, color = "grey10") +
          xlim_tree(x_limit) +
          scale_color_gradient2(low = color_lowval_tmp, mid = color_midval_tmp, high = color_highval_tmp, 
                                midpoint = 1,
                                limits = c(0, omega_highcutoff),
                                na.value = missingcol,
                                breaks = c(0, 1, 5, 10)) +
          guides(color = guide_colorbar(title = legend_title)) +
          ggtitle( title ) +
          theme_tmp ; tree_p1
        # tree_p1 = tree_p1 %>% rotate(25) ; tree_p1
      }
    }else{
      # p value
      if(sum(is.na(df_3$`Uncorrected_P.value<0.05`)) == nrow(df_3)){
        
        # if all rows are na (no color for the whole tree)
        tree_p1 = ggtree(tree_phy, size = linewidth, aes_(col = as.name(tt)), color = "grey50") + 
          geom_tiplab(offset = offset_tip, size = tipLabelSize, color = "grey10") +
          xlim_tree(x_limit) +
          guides(color = guide_colorbar(title = legend_title, reverse = TRUE)) +
          ggtitle(title ) +
          theme_tmp ; tree_p1
        # tree_p1 = tree_p1 %>% rotate(25) ; tree_p1
        
      }else{
        
        # normal cases:
        tree_p1 = ggtree(tree_phy, size = linewidth, aes_(col = as.name(tt))) +
          geom_tiplab(offset = offset_tip, size = tipLabelSize, color = "grey10") +
          xlim_tree(x_limit) +
          scale_color_continuous(low = color_tmp, high = color_highp_tmp, na.value = "grey50",
                                 limits = c(0, 0.05)) +
          guides(color = guide_colorbar(title = legend_title, reverse = TRUE)) +
          ggtitle(title ) +
          theme_tmp ; tree_p1
        # tree_p1 = tree_p1 %>% rotate(25) ; tree_p1
      }
    }
    if (withBranchLen == 'y'){
      tree_p1 = tree_p1 + geom_treescale(x = time_scale_x, y = 0.1, fontsize = tipLabelSize) } 
    if (is.rooted(tree_phy) ){
      tree_p1 = tree_p1 + geom_rootedge(rootedge = rootedge, color = "grey50", size = linewidth)
    }
    
    return(tree_p1)
  })
  pp = arrangeGrob(grobs = p, nrow = 1)
  # dev.off()
  # grid.draw(pp)
  return(pp)
}


get_absrel_table = function(s){
  # s = target_working_files[2]
  print(s)
  temp_f_path = s;
  tmp = fromJSON(txt = temp_f_path) ; # glimpse(tmp)
  
  class = gsub(project_path, "", dirname(temp_f_path))
  class = gsub("/", "", class)
  current_gene = gsub("\\.json", "", basename(temp_f_path))
  
  omega_highcutoff = 10
  
  ### get tree
  tmp_tree = tmp$input$trees$`0`; # glimpse(tmp_tree)
  tree <- read.tree(text = paste0(tmp_tree, ";")) 
  nTips = length(tree$tip.label)
  
  tree = tree %>% as_tibble() %>% 
    bind_cols("class" = class, "gene" = current_gene, "nTips" = nTips)
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


# plot tree for checking tree topology
plot_tree_withNodelabel = function(path){
  s = path
  print(s)
  
  tmp_name = gsub("\\.json", "", basename(s))
  
  temp_f_path = s;
  tmp = fromJSON(txt = temp_f_path) ; # glimpse(tmp)
  
  tmp_tree = tmp$input$trees$`0`; # glimpse(tmp_tree)
  tree <- read.tree(text = paste0(tmp_tree, ";"))
  x_limit = max(tree$edge) * 0.8
  rootedge = min(tree$edge, na.rm = TRUE) * 0.002
  
  # change name
  # ind = sapply(tree$tip.label, function(s){grep(s, name_exchange$X1)})
  # tree$tip.label = name_exchange$X2[ind]
  
  tree_tb = tree %>% as_tibble()
  internalnodes_nodeID = sapply(internalnodes, function(s){
    grep(s, tree_tb$label)
  })
  
  tree_p = ggtree(tree, size = 0.25) + 
    geom_tiplab(offset = 2, size = FontSize*0.75) +
    geom_point2(aes(subset=(node %in% internalnodes_nodeID)), shape=19, size=1,
                color='gold') +
    geom_nodelab(aes(subset=(node %in% internalnodes_nodeID)),
                 size = FontSize*0.75, hjust = 1.1, # vjust = 1.5, nudge_x = -1,
                 nudge_y = -0.4 ) +
    xlim_tree(x_limit) +
    ggtitle(tmp_name) +
    geom_text(aes(label=node), size = FontSize*0.75, hjust=-.3, color = 'red') +
    theme_tree() + theme(legend.position = 'none',
                         plot.background = element_blank()); tree_p
}

# 3.2 plot ====
# 3.2.3 cladogram, uncorrected pvalues less than 0.05  ====
# tmp_heatmap_val = c("LRT","Uncorrected.P.value", "Uncorrected_P.value.0.05",  "Corrected.P.value", "Omega1", "Omega2")

# test
# pl = plot_absrel_variHeatMapVal(s =target_working_files[1], withBranchLen = 'n', tmp_heatmap_val =  "Uncorrected_P.value.0.05")
# pl = plot_absrel_variHeatMapVal(s =target_working_files[1], withBranchLen = 'n', tmp_heatmap_val =  "Omega2")
# dev.off()
# grid.draw(pl)
run_tests = c("Uncorrected_P.value<0.05", "Corrected P-value") # "HighestOmega_cutoff"
rep_n = length(run_tests)

pl = mapply(plot_absrel_variHeatMapVal,
            s = rep(target_working_files, each = rep_n),
            withBranchLen = rep('n', length(target_working_files)*rep_n),
            tmp_heatmap_val = rep(run_tests, length(target_working_files)),
            x_limit = rep(c(16, 16), each = rep_n), 
            tipLabelSize = rep(2, length(target_working_files)*2), 
            SIMPLIFY = F)

p_out = arrangeGrob(grobs = pl, ncol = rep_n)
# dev.off()
# grid.draw(p_out)

graph_path = paste0(project_path, "absrel_batch1_cladogram_", Sys.Date(), ".pdf")       # all data
pdf(graph_path, width = Width_HalfCol*0.4*rep_n*1.5, height = Width_HalfCol*0.5*length(target_working_files),
    pointsize = AxisTxFontSizeSize, onefile = TRUE)
# png(graph_path, width = Width_HalfCol*1, height = Width_HalfCol*1, unit = "in", res = 600)
grid.draw(p_out)
dev.off()
openfile(graph_path)

# 3.2.3.5 cladogram [common name] [baseline omega] ====
path = paste0(getwd(), "/FromMaude/20221130/HighOmega_ForMaggie.xlsx")
wanted = read_excel(path, col_names = F)
wanted_g = na.omit(wanted$...1)

wanted_files = sapply(wanted_g, function(s){target_working_files[grepl(s, target_working_files)]})

run_tests = c("Uncorrected_P.value<0.05", "Corrected P-value", "Baseline MG94xREV omega ratio")
rep_n = length(run_tests)

pl = mapply(plot_absrel_variHeatMapVal,
            s = rep(target_working_files, each = rep_n),
            withBranchLen = rep('n', length(target_working_files)*rep_n),
            tmp_heatmap_val = rep(run_tests, length(target_working_files)),
            x_limit = rep(30, length(target_working_files)), 
            tipLabelSize = rep(2, length(target_working_files)*rep_n), 
            targetName = list(speciesNAME$`Common Name`),
            branchlen_col = rep("Full adaptive model (synonymous subs/site)", 
                                length(target_working_files)*rep_n),
            SIMPLIFY = F)

p_out = arrangeGrob(grobs = pl, ncol = rep_n)
# dev.off()
# grid.draw(p_out)

graph_path = paste0(project_path, "absrel_T1Rs_cladogram_commonName_baselineOmega_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*0.6*rep_n, height = Width_HalfCol*0.7*length(target_working_files),
    pointsize = AxisTxFontSizeSize, onefile = TRUE)
# png(graph_path, width = Width_HalfCol*1, height = Width_HalfCol*1, unit = "in", res = 600)
grid.draw(p_out)
dev.off()
openfile(graph_path)

# 5. get result ====
path = paste0(getwd(), "/PhyloDiet/fig3/speciesNames_10_14.txt")
name_exchange = read_tsv(file = path, col_names = F)

# The nodes we care about most are 12,18, and 19. But if it's possible to plot also 4,5,8,10,, then 12,18,19,20,21, and all manakin tips (from neopelma=nech onward), then 13,15; CEOR, PICH, EMTR,that would be awesome, then we can see what we show in the end.

# wanted = c("Node12", "Node18", "Node19", 
#            "Node20", "Node21",
#            "NECH", "COAL", "LECO", "MAVI", "PIFI", 
#            "Node4", "Node5", "Node8", "Node10",
#            "Node13", "Node15", "CEOR", "PICH", "EMTR")
# internalnodes = wanted[grepl("Node", wanted)]

wanted = c("Node12", "Node18", "Node19", "Node21")
internalnodes = wanted[grepl("Node", wanted)]


# # 5.1 get example/representative tree [1: from absrel] ====
s = target_working_files[1] # 20221018/1_Muscle/CLCN1.json
print(s)
temp_f_path = s;
tmp = fromJSON(txt = temp_f_path) ; # glimpse(tmp)

tmp_tree = tmp$input$trees$`0`; # glimpse(tmp_tree)
tree <- read.tree(text = paste0(tmp_tree, ";"))
x_limit = max(tree$edge) * 1
rootedge = min(tree$edge, na.rm = TRUE) * 0.002

# change name
ind = sapply(tree$tip.label, function(s){grep(s, name_exchange$X1)})
tree$tip.label = name_exchange$X2[ind]

tree_tb = tree %>% as_tibble()
internalnodes_nodeID = sapply(internalnodes, function(s){
  grep(s, tree_tb$label)
})

tree_p = ggtree(tree, size = 0.25) + 
  geom_tiplab(offset = 0.1, size = FontSize*0.75) +
  geom_point2(aes(subset=(node %in% internalnodes_nodeID)), shape=19, size=1,
              color='gold') +
  geom_nodelab(aes(subset=(node %in% internalnodes_nodeID)),
               size = FontSize*0.75, hjust = 1.1, 
               # vjust = 1.5
               # nudge_x = -1,
               nudge_y = -0.4
  )+
  xlim_tree(x_limit) +
  # geom_text(aes(label=node), hjust=-.3) +
  # geom_rootedge(rootedge, color = 'black', size = 0.25) +
  # geom_tippoint(aes(color = group)) +
  # scale_color_manual(values=c("black"), na.value =  "black") + 
  # coord_flip() +
  theme_tree() + theme(legend.position = 'none',
                       plot.background = element_blank()); tree_p
tree_p = tree_p %>% rotate(15); tree_p

# graph_path = paste0(project_path, "/tree_", Sys.Date(), ".pdf")
# pdf(graph_path, width = Width_HalfCol*0.5, height = Width_HalfCol*0.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
# tree_p
# dev.off()
# openfile(graph_path)

### new but stop working 20221222
s = target_working_files[1] # 20221018/1_Muscle/CLCN1.json

input_path = s
rep_n = length(input_path)

pl = mapply(plot_tree_withNodelabel,
            path = input_path,
            SIMPLIFY = F)

p_out = arrangeGrob(grobs = pl, nrow = 1)
# dev.off()
# grid.draw(p_out)

# graph_path = paste0(project_path, "/tree_CLCN1_2_", Sys.Date(), ".pdf")
# pdf(graph_path, width = Width_HalfCol*0.5, height = Width_HalfCol*0.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# tree_p
# dev.off()
# openfile(graph_path)
 
# 5.1.2 example tree maude sent 20221027
tree = read.nexus(file = paste0(getwd(), "/FromMaude/20221012/SPECIES_TREE_2.treNoLength.nex"))
plot(tree)

# graph_path = paste0(project_path, "/Maude_sp_trees_", Sys.Date(), ".pdf")
# pdf(graph_path, width = Width_HalfCol*0.5*1, height = Width_HalfCol*0.5*1, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# plot(tree)
# dev.off()
# openfile(graph_path)


# 5.2 get absrel results ====
output = lapply(target_working_files, get_absrel_table) %>% bind_rows()
unique(output$Name)

output_m = output %>% 
  mutate(Name = gsub("Aquila", "AQCH", Name, ignore.case = TRUE),
         Name = gsub("Cephalopterus", "CEOR", Name, ignore.case = TRUE),
         Name = gsub("Corapipo", "COAL", Name, ignore.case = TRUE),
         Name = gsub("Corvus", "COCO", Name, ignore.case = TRUE),
         Name = gsub("Empidonax", "EMTR", Name, ignore.case = TRUE),
         Name = gsub("Ficedula", "FIAL", Name, ignore.case = TRUE),
         Name = gsub("Lepidothrix", "LECO", Name, ignore.case = TRUE),
         Name = gsub("Manacus", "MAVI", Name, ignore.case = TRUE),
         Name = gsub("Melopsittacus", "MEUN", Name, ignore.case = TRUE),
         Name = gsub("Neopelma", "NECH", Name, ignore.case = TRUE),
         Name = gsub("Pipra", "PIFI", Name, ignore.case = TRUE),
         Name = gsub("Piprites", "PICH", Name, ignore.case = TRUE),
         Name = gsub("Rhegmatorhina", "RHHO", Name, ignore.case = TRUE),
         Name = gsub("Smithornis", "SMCA", Name, ignore.case = TRUE))
unique(output_m$Name)
# name_exchange

output_diff_tree = output_m %>% select(gene, nTips) %>% distinct()
output_diff_treenot14 = output_diff_tree %>% filter(nTips != 14)

# 5.3. node 12, 18, 19, 21; all genes ====

# 5.3.1 check for 14 tips's tree topology  ====
output_14 = output_m %>% filter(nTips == 14)
length(unique(output_14$gene)) # 91 > 94 > 100 > 104

# ## comment start 20221221
# 
# output_14_check = output_14 %>% filter(Name %in% c(wanted)) %>% 
#   dplyr::select(Name, parent, node) %>% distinct() %>% 
#   arrange(Name)
# # Name   parent  node
# # <chr>   <int> <int>
# # 1 Node12     19    20 # major
# # 2 Node12     18    19
# # 3 Node12     16    18
# # 4 Node18     20    23 # major
# # 5 Node18     20    21
# # 6 Node18     19    22
# # 7 Node18     20    22
# # 8 Node18     19    20
# # 9 Node18     18    19
# # 10 Node19     23    24 # major
# # 11 Node19     21    22
# # 12 Node19     22    23
# # 13 Node19     20    21
# # 14 Node21     25    26 # major
# # 15 Node21     24    25
# # 16 Node21     22    23
# # 17 Node21     23    24
# # > seems to have some variation
# 
# ## correct/major
# output_14_look = output_14 %>% filter(Name == "Node12", node == 20)
# m = unique(output_14_look$gene)
# length(m) # 89/95 genes are major tree 
# 
# # plot these trees 
# input_path = sapply(m, function(s){ target_working_files[grepl(s, target_working_files)]})
# input_path = unique(unlist(input_path))
# rep_n = length(input_path)
# Nrow = ceiling(rep_n/4)
# 
# pl = mapply(plot_tree_withNodelabel, path = input_path, SIMPLIFY = F)
# p_out = arrangeGrob(grobs = pl, ncol = 4)
# # dev.off()
# # grid.draw(p_out)
# 
# graph_path = paste0(project_path, "/N14trees_", Sys.Date(), ".pdf")
# pdf(graph_path, width = Width_HalfCol*0.5*4, height = Width_HalfCol*0.5*Nrow, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# grid.draw(p_out)
# dev.off()
# openfile(graph_path)
# 
# # after plotting out the 87 tree, they aren't always the same variations
# output_14_look2 = output_14 %>% filter(gene %in% m) %>% 
#   # filter(Name == "PIFI", node == 10)
#   filter(Name == "Node18", node == 23)
# x = unique(output_14_look2$gene); x # 
# length(x) # 80 > plot these tree to see if they are all major
# 
# output_14_look3 = output_14 %>% filter(gene %in% x) %>% 
#   filter(Name == "Node19", node == 24) 
# x2 = unique(output_14_look3$gene); #x2  
# length(x2) # 80
# 
# output_14_look4 = output_14 %>% filter(gene %in% x2) %>% 
#   filter(Name == "Node21", node == 26) 
# x3 = unique(output_14_look4$gene); #x2  
# length(x3) # 65 >> plot this
# # plot these trees 
# 
# input_path = sapply(x3, function(s){ target_working_files[grepl(paste0(s, ".json"), target_working_files)]})
# input_path = unique(unlist(input_path))
# rep_n = length(input_path); #65
# Nrow = ceiling(rep_n/4)
# 
# pl = mapply(plot_tree_withNodelabel, path = input_path, SIMPLIFY = F)
# p_out = arrangeGrob(grobs = pl, ncol = 4)
# # dev.off()
# # grid.draw(p_out)
# 
# graph_path = paste0(project_path, "/N14trees_node12eq20_node18eq23_node21eq26_", Sys.Date(), ".pdf")
# pdf(graph_path, width = Width_HalfCol*0.5*4, height = Width_HalfCol*0.5*Nrow, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# grid.draw(p_out)
# dev.off()
# openfile(graph_path)
# 
# 
# x4 = setdiff(x2, x3)
# length(x4) # 15  >> plot this
# 
# input_path = sapply(x4, function(s){ target_working_files[grepl(paste0(s, ".json"), target_working_files)]})
# input_path = unique(unlist(input_path))
# rep_n = length(input_path); #15
# Nrow = ceiling(rep_n/4)
# 
# pl = mapply(plot_tree_withNodelabel, path = input_path, SIMPLIFY = F)
# p_out = arrangeGrob(grobs = pl, ncol = 4)
# # dev.off()
# # grid.draw(p_out)
# 
# graph_path = paste0(project_path, "/N14trees_node12eq20_node18eq23_node21Noteq26_", Sys.Date(), ".pdf")
# pdf(graph_path, width = Width_HalfCol*0.5*4, height = Width_HalfCol*0.5*Nrow, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# grid.draw(p_out)
# dev.off()
# openfile(graph_path)
# 
# 
# output_14_look = output_14 %>% filter(Name == "Node12", node != 20)
# m2 = unique(output_14_look$gene)
# length(m2) # 5
# 
# input_path = sapply(m2, function(s){ target_working_files[grepl(paste0(s, ".json"), target_working_files)]})
# input_path = unique(unlist(input_path))
# rep_n = length(input_path); 
# Nrow = ceiling(rep_n/5)
# 
# pl = mapply(plot_tree_withNodelabel, path = input_path, SIMPLIFY = F)
# p_out = arrangeGrob(grobs = pl, ncol = 5)
# # dev.off()
# # grid.draw(p_out)
# 
# graph_path = paste0(project_path, "/N14trees_node12Noteq20_", Sys.Date(), ".pdf")
# pdf(graph_path, width = Width_HalfCol*0.5*5, height = Width_HalfCol*0.5*Nrow, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# grid.draw(p_out)
# dev.off()
# openfile(graph_path)
# 
# # 
# output_14_look5 = output_14 %>% filter(gene %in% m) %>% 
#   filter(Name == "Node18", node != 23)
# y = unique(output_14_look5$gene); y
# length(y) # 9
# 
# input_path = sapply(y, function(s){ target_working_files[grepl(paste0(s, ".json"), target_working_files)]})
# input_path = unique(unlist(input_path))
# rep_n = length(input_path); 
# Nrow = ceiling(rep_n/5)
# 
# pl = mapply(plot_tree_withNodelabel, path = input_path, SIMPLIFY = F)
# p_out = arrangeGrob(grobs = pl, ncol = 5)
# # dev.off()
# # grid.draw(p_out)
# 
# graph_path = paste0(project_path, "/N14trees_node12eq20_node18Noteq23_2_", Sys.Date(), ".pdf")
# pdf(graph_path, width = Width_HalfCol*0.5*5, height = Width_HalfCol*0.5*Nrow, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# grid.draw(p_out)
# dev.off()
# openfile(graph_path)
# 
# xxx = output_14 %>% filter(gene %in% y) %>% 
#   dplyr::select(class, gene, Name, node, `Uncorrected P-value`, `Corrected P-value`) %>% 
#   filter(Name %in% wanted) #%>% 
#   # group_by(Name) %>% tally()
# 
# 
# # 5.3.2 check for 13 tips's tree topology  ====
output_13 = output_m %>% filter(nTips != 14)
length(unique(output_13$gene)) # 4 >5

# m13 = unique(output_13$gene)
# length(m13) # 4
# 
# output_13_check = output_13 %>% filter(Name %in% c(wanted)) %>%
#   dplyr::select(Name, parent, node) %>% distinct() %>%
#   arrange(Name)
# # Name   parent  node
# # <chr>   <int> <int>
# # 1 Node12     18    19
# # 2 Node18     22    23
# # 3 Node19     23    24
# # > after checking the trees of these 4
# 
# input_path = sapply(m13, function(s){ target_working_files[grepl(paste0(s, ".json"), target_working_files)]})
# input_path = unique(unlist(input_path))
# rep_n = length(input_path); #15
# Nrow = ceiling(rep_n/4)
# 
# pl = mapply(plot_tree_withNodelabel, path = input_path, SIMPLIFY = F)
# p_out = arrangeGrob(grobs = pl, ncol = 4)
# dev.off()
# grid.draw(p_out)
# 
# graph_path = paste0(project_path, "/N13trees_", Sys.Date(), ".pdf")
# pdf(graph_path, width = Width_HalfCol*0.5*4, height = Width_HalfCol*0.5*Nrow, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# grid.draw(p_out)
# dev.off()
# openfile(graph_path)
# 
# ## comment end 20221221

# 5.3.3 mod atypical trees ====
# summarize here:
# /Users/maggie/ownCloud/Baldwin/Projects/Manakin/PhyloDiet/fig3/absrel/20221018/check_tree.pptx

output_14_selCol = output_14 %>% select(class, gene, Name, node, `Uncorrected P-value`, `Corrected P-value`)
output_13_selCol = output_13 %>% select(class, gene, Name, node, `Uncorrected P-value`, `Corrected P-value`) 

# ## 9+ different topologies for nTip=14 tree
# # x3 # 65
# x5 = c("HYAL2", "GUSB", "SIAE")
# x6 = "SLC9A3"
# x7 = setdiff(x4, c(x5, x6)) # 11
# x8 = "CPLX4"
# x9 = "SWS1"
# x10 = "DAB1"
# x11 = "CBX3"
# x12 = "TERB1"
# # y # 9
# 
# # x3: 65
# output_14_sel1 = output_14_selCol %>%
#   filter(gene %in% x3) %>% 
#   filter(Name %in% wanted)
# length(unique(output_14_sel1$gene)) # 65
# 
# # x5: 3
# output_14_sel2 = output_14_selCol %>%
#   filter(gene %in% x5) %>% 
#   filter(node %in% c(20, 23, 24, 25))
# length(unique(output_14_sel2$gene)) # 3
# 
# # x6: 1
# output_14_sel3 = output_14_selCol %>%
#   filter(gene %in% x6) %>% 
#   filter(node %in% c(20, 23, 24))
# length(unique(output_14_sel3$gene)) # 1
# 
# # x7: 11
# output_14_sel4 = output_14_selCol %>%
#   filter(gene %in% x7) %>% 
#   filter(node %in% c(20, 23, 24))
# length(unique(output_14_sel4$gene)) # 11
# 
# # x8: 1
# output_14_sel5 = output_14_selCol %>%
#   filter(gene %in% x8) %>% 
#   filter(node %in% c(19, 22, 23, 25))
# length(unique(output_14_sel5$gene)) # 1
# 
# # x9: 1
# output_14_sel6 = output_14_selCol %>%
#   filter(gene %in% x9) %>%
#   filter(node %in% c(19, 22, 25))
# length(unique(output_14_sel6$gene)) # 1
# 
# # x10: 1
# output_14_sel7 = output_14_selCol %>%
#   filter(gene %in% x10) %>%
#   filter(node %in% c(20, 21))
# length(unique(output_14_sel7$gene)) # 1
# 
# # x11: 1
# output_14_sel8 = output_14_selCol %>%
#   filter(gene %in% x11) %>%
#   filter(node %in% c(18:19))
# length(unique(output_14_sel8$gene)) # 1
# 
# # x12: 1
# output_14_sel9 = output_14_selCol %>%
#   filter(gene %in% x12) %>%
#   filter(node %in% c(19, 20, 21, 23))
# length(unique(output_14_sel9$gene)) # 1
# 
# output_14_selFinal = output_14_sel1 %>% 
#   bind_rows(output_14_sel2) %>% 
#   bind_rows(output_14_sel3) %>% 
#   bind_rows(output_14_sel4) %>% 
#   bind_rows(output_14_sel5) %>% 
#   bind_rows(output_14_sel6) %>% 
#   bind_rows(output_14_sel7) %>% 
#   bind_rows(output_14_sel8) %>% 
#   bind_rows(output_14_sel9)
# unique(output_14_selFinal$Name)
# length(unique(output_14_selFinal$gene)) # 85
# # conclusion, the names are pretty consistent, thus I don't really need to check so deligently

output_14_selFinal = output_14_selCol %>% 
  filter(Name %in% wanted)
unique(output_14_selFinal$Name)
length(unique(output_14_selFinal$gene)) # 94 > 100 > 105

# N13 tips tree
x13 = c("LMOD1", "SLC26A5")
x14 = "SLC52A3"
x15 = "LDHB"
x16 = "COL4A3"

output_13_sel1 = output_13_selCol %>%
  filter(gene %in% x13) %>%
  filter(node %in% c(19, 21, 22)) %>% 
  mutate(Name = gsub("Node16", "Node18", Name),
         Name = gsub("Node17", "Node19", Name))
length(unique(output_13_sel1$gene)) # 2

output_13_sel2 = output_13_selCol %>%
  filter(gene %in% x14) %>%
  filter(node %in% c(19, 21, 22, 24)) %>% 
  mutate(Name = gsub("Node16", "Node18", Name),
         Name = gsub("Node17", "Node19", Name),
         Name = gsub("Node19", "Node21", Name))
length(unique(output_13_sel2$gene)) # 1

output_13_sel3 = output_13_selCol %>%
  filter(gene %in% x15) %>%
  filter(node %in% c(18, 21, 22)) %>% 
  mutate(Name = gsub("Node10", "Node12", Name),
         Name = gsub("Node16", "Node18", Name),
         Name = gsub("Node17", "Node19", Name))
length(unique(output_13_sel3$gene)) # 1

output_13_sel4 = output_13_selCol %>%
  filter(gene %in% x16) %>%
  filter(node %in% c(17, 18)) %>% 
  mutate(Name = gsub("Node14", "Node18", Name),
         Name = gsub("Node17", "Node21", Name))
length(unique(output_13_sel4$gene)) # 1

# combine
output_14_sel2 = output_14_selFinal %>% 
  bind_rows(output_13_sel1) %>% bind_rows(output_13_sel2) %>% 
  bind_rows(output_13_sel3) %>% bind_rows(output_13_sel4) %>% 
  filter(grepl("^\\d", class)) %>%
  mutate(class = gsub("^\\d\\w{0,1}_", "", class)) %>%
  mutate(corrected5 = ifelse(`Corrected P-value` <=0.05, TRUE, NA)) 
unique(output_14_sel2$class)
dim(output_14_sel2) #  394   7 >414   7

# make dummy data for not-existing nodes and get 'class'
gene_runthrough = unique(output_14_sel2$gene)
output_14_sel2_count = output_14_sel2 %>%
  group_by(gene, class) %>%
  mutate(n = n())
tmp_out = lapply(gene_runthrough, function(s){
  # s = "SPO11"
  tmp_sub = output_14_sel2_count %>% filter(grepl(s, gene))
  tmp_n = unique(tmp_sub$n)
  tmp_class = unique(tmp_sub$class)
  tmp_nodes = unique(tmp_sub$Name)
  if(tmp_n<4){
    toadd = setdiff(wanted, tmp_nodes)
    # names(tmp_sub)
    df_out = tibble("Name" =  toadd, "gene" = s, 
                    "Uncorrected P-value" = 1,
                    "Corrected P-value" = 1,  
                    "class" = tmp_class, "corrected5" = NA )
  }}) %>% bind_rows()

output_mod = output_14_sel2 %>% bind_rows(tmp_out)
# output_mod$gene = sapply(output_mod$gene, function(s){unlist(strsplit(s, "_"))[1]})
dim(output_mod) #  432   7 > 440   7

# # checking start ====
# output_14_sel3 = output_mod %>%
#   group_by(gene, class) %>%
#   summarise(n = n())
# 
# # target_working_files[grepl("CBX3", target_working_files)]
# x = output_14_sel3 %>% filter(n<4)
# y = unique(x$gene)
# x = unlist(sapply(y, function(s){target_working_files[grep(s, target_working_files)]}))
# 
# run_tests = c("Corrected P-value") # "HighestOmega_cutoff"
# rep_n = length(run_tests)
# 
# pl = mapply(plot_absrel_variHeatMapVal_nodelab,
#             s = rep(x, each = rep_n),
#             withBranchLen = rep('n', length(x)*rep_n),
#             tmp_heatmap_val = rep(run_tests, length(x)),
#             x_limit = rep(c(10), each = rep_n), 
#             tipLabelSize = rep(2, length(x)), 
#             SIMPLIFY = F)
# 
# p_out = arrangeGrob(grobs = pl, ncol = rep_n)
# # dev.off()
# # grid.draw(p_out)
# 
# graph_path = paste0(project_path, "absrel_checking_cladogram_", Sys.Date(), ".pdf")       # all data
# pdf(graph_path, width = Width_HalfCol*0.4*rep_n*1.5, height = Width_HalfCol*0.5*length(x),
#     pointsize = AxisTxFontSizeSize, onefile = TRUE)
# # png(graph_path, width = Width_HalfCol*1, height = Width_HalfCol*1, unit = "in", res = 600)
# grid.draw(p_out)
# dev.off()
# openfile(graph_path)
# 
# # checking end

# # change tip name
# ind_tips = which(!grepl("^Node", output_mod$Name) & !grepl("FUFI", output_mod$Name))
# ind = sapply(output_mod$Name[ind_tips], function(s){grep(s, name_exchange$X1, ignore.case = T)})
# output_mod$Name[ind_tips] = name_exchange$X2[ind]
# 
# ind_tips = which(!grepl("^Node", wanted) & !grepl("FUFI", wanted))
# ind = sapply(wanted[ind_tips], function(s){grep(s, name_exchange$X1, ignore.case = T)})
# wanted[ind_tips] = name_exchange$X2[ind]

unique(output_mod$Name)
output_mod$Name = factor(output_mod$Name, levels = rev(wanted))
unique(output_mod$class)
class_order = c("Diet_Carbohydrates", "Diet_Protein", "Diet_Fat", "Diet",
                "Vision", "Audition", "Telomere", "Muscle", "Sperm" )
output_mod$class = factor(output_mod$class, levels = class_order)

length(unique(output_mod$gene)) #110
length(unique(genelist$Gene)) # 110
setdiff(genelist$Gene, output_mod$gene) # 0 > good

output_mod$gene = factor(output_mod$gene, levels = genelist$Gene)
# output_mod$corrected5 = factor(output_mod$corrected5, levels = c(TRUE))

theme_m2 = theme(legend.position = "right",
                 axis.text.y = element_text(size = AxisTxFontSizeSize_s),
                 axis.title.y = element_blank(),
                 axis.title.x = element_blank(), 
                 strip.text = element_text(size = AxisTxFontSizeSize_s, 
                                           margin = margin(2,0,2,0, "pt")),
                 legend.key.size = unit(8, 'points'), 
                 legend.background = element_blank(),
                 legend.key = element_blank(),
                 legend.title = element_text(size = AxisTxFontSizeSize_s),
                 legend.text = element_text(size = AxisTxFontSizeSize_s))

gene_for_mainFig = genelist %>% filter(is.na(`Keep for SI not main figure`))
output_mod_main = output_mod %>% filter(gene %in% gene_for_mainFig$Gene)

### color 2 
color_highp_tmp2 = 'white'
color_low_tmp2 = 'red3'

color_highp_tmp = 'chocolate2'
color_midp_tmp = 'khaki1' 
color_tmp = 'red3'

theme_m2 = theme(legend.position = "right",
                 axis.text.y = element_text(size = AxisTxFontSizeSize_s),
                 axis.title.y = element_blank(),
                 axis.title.x = element_blank(), 
                 strip.text = element_text(size = AxisTxFontSizeSize_s, 
                                           margin = margin(2,0,2,0, "pt")),
                 legend.key.size = unit(8, 'points'), 
                 legend.background = element_blank(),
                 legend.key = element_blank(),
                 legend.title = element_text(size = AxisTxFontSizeSize_s),
                 legend.text = element_text(size = AxisTxFontSizeSize_s))


p = ggplot(output_mod_main, aes(y = Name, x = gene, fill = `Uncorrected P-value`)) +
  geom_tile(color = "grey90", size = 0.1) +
  # geom_tile(aes(color = corrected5, size = `Corrected P-value`)) +
  geom_tile(aes(color = corrected5), size = 0.2, height = 0.9, width = 0.9) +
  scale_fill_gradientn(colors = c(color_tmp, 'rosybrown1', color_midp_tmp, color_highp_tmp),
                       na.value = "transparent",
                       limits = c(0, 0.1), breaks = c(0, 0.05, 0.1)) +
  # scale_fill_gradient(low = color_low_tmp2, high = color_highp_tmp2, 
  # na.value = "transparent",
  # limits = c(0, 0.05)) +
  # scale_color_manual(values = c("black", "transparent"), na.value = "transparent") +
  scale_color_manual(values = c("black"), na.value = "transparent") +
  # scale_size_continuous(range = c(0.1, 0.5)) +
  scale_y_discrete(position = "right", expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(.~class, scales = 'free', space = 'free', switch = "x") +
  guides(colour = "none") +
  # facet_wrap(.~class, scales = 'free', strip.position = "bottom", nrow = 2, shrink = TRUE) +
  theme_m + theme_m2 + 
  theme(axis.text.x = element_text(size = AxisTxFontSizeSize_s, 
                                   angle = 40, hjust = 1),
        legend.key.size = unit(5, 'pt')); p


# color 0.5 to 0.1 gray
output_mod2 = output_mod_main %>%
  mutate(unco = ifelse(`Uncorrected P-value` <= 0.1, `Uncorrected P-value`, NA))

p = ggplot(output_mod2, aes(y = Name, x = gene, fill = unco)) +
  geom_tile(color = "grey90", size = 0.1) +
  # geom_tile(aes(color = corrected5, size = `Corrected P-value`)) +
  geom_tile(aes(color = corrected5), size = 0.2, height = 0.9, width = 0.9) +
  # scale_fill_gradientn(colors = c(color_tmp, 'rosybrown1'),
                       # na.value = "grey90",
                       # limits = c(0, 0.1), breaks = c(0, 0.05, 0.1)) +
  scale_fill_gradient2(low = color_low_tmp2, mid = color_highp_tmp2, high = 'grey80',
                       midpoint = 0.05,
                      na.value = "transparent",
                      limits = c(0, 0.1), breaks = c(0, 0.05, 0.1)) +
  # scale_color_manual(values = c("black", "transparent"), na.value = "transparent") +
  scale_color_manual(values = c("black"), na.value = "transparent") +
  # scale_size_continuous(range = c(0.1, 0.5)) +
  scale_y_discrete(position = "right", expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(.~class, scales = 'free', space = 'free', switch = "x") +
  # facet_wrap(.~class, scales = 'free', strip.position = "bottom", nrow = 2, shrink = TRUE) +
  guides(colour = "none") +
  # facet_wrap(.~class, scales = 'free', strip.position = "bottom", nrow = 2, shrink = TRUE) +
  theme_m + theme_m2 +
  theme(axis.text.x = element_text(size = AxisTxFontSizeSize_s,
                                   angle = 40, hjust = 1),
        legend.key.size = unit(5, 'pt'),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()
        # panel.border = element_rect(color = "black", size = 0.1, fill = NA)
  ); p

# color 0.5 to 0.1 gray -ver3 [using] ====
color_highp_tmp2 = 'mistyrose'
output_mod2 = output_mod_main %>%
  mutate(unco = ifelse(`Uncorrected P-value` < 0.1, `Uncorrected P-value`, NA))
x = output_mod2 %>% filter(`Uncorrected P-value` < 0.1 & `Uncorrected P-value` >0.05 )
nrow(x) # 9

p = ggplot(output_mod2, aes(y = Name, x = gene, fill = unco)) +
  geom_tile(color = "grey90", size = 0.1) +
  # geom_tile(aes(color = corrected5, size = `Corrected P-value`)) +
  geom_tile(aes(color = corrected5), size = 0.2, height = 0.9, width = 0.9) +
  # scale_fill_gradientn(colors = c(color_tmp, 'rosybrown1'),
  # na.value = "grey90",
  # limits = c(0, 0.1), breaks = c(0, 0.05, 0.1)) +
  scale_fill_gradient2(low = color_low_tmp2, mid = color_highp_tmp2, high = 'grey20',
                       midpoint = 0.05,
                       na.value = "transparent",
                       limits = c(0, 0.1), breaks = c(0, 0.05, 0.1)) +
  # scale_color_manual(values = c("black", "transparent"), na.value = "transparent") +
  scale_color_manual(values = c("black"), na.value = "transparent") +
  # scale_size_continuous(range = c(0.1, 0.5)) +
  scale_y_discrete(position = "right", expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(.~class, scales = 'free', space = 'free', switch = "x") +
  # facet_wrap(.~class, scales = 'free', strip.position = "bottom", nrow = 2, shrink = TRUE) +
  guides(colour = "none") +
  # facet_wrap(.~class, scales = 'free', strip.position = "bottom", nrow = 2, shrink = TRUE) +
  theme_m + theme_m2 +
  theme(axis.text.x = element_text(size = AxisTxFontSizeSize_s,
                                   angle = 40, hjust = 1),
        legend.key.size = unit(5, 'pt'),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()
        # panel.border = element_rect(color = "black", size = 0.1, fill = NA)
  ); p

# x = output_mod %>% filter(grepl("SI", gene))

n_row = 10
n_empty_row = 0
n_col = 8
n_col_tree = 1
lay = matrix(c( rep(c(rep( 1, n_col_tree), rep(2, (n_col-n_col_tree) )), (n_row-n_empty_row) ),
                rep(c(rep(NA, n_col_tree), rep(2, (n_col-n_col_tree) )), n_empty_row)), 
             nrow = n_row, byrow = TRUE); lay
p_out = arrangeGrob(grobs = list(tree_p, p), layout_matrix = lay)
# p_out = arrangeGrob(grobs = list(tree_p, p), nrow = 1, widths = c(2,3))
# dev.off()
# grid.draw(p_out)

graph_path = paste0(project_path, "/Summary_4nodes_mainFig_", Sys.Date(), "_c2g.pdf")
pdf(graph_path, width = Width_HalfCol*2.5, height = Width_HalfCol*0.24, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
grid.draw(p_out)
dev.off()
openfile(graph_path)

x  = genelist %>% filter(is.na(`Keep for SI not main figure`))

# 2. all node; all genes ====
### get tree
s = target_working_files[6]
print(s)
temp_f_path = s;
tmp = fromJSON(txt = temp_f_path) ; # glimpse(tmp)

tmp_tree = tmp$input$trees$`0`; # glimpse(tmp_tree)
tree <- read.tree(text = paste0(tmp_tree, ";"))
x_limit = max(tree$edge) * 1
rootedge = min(tree$edge, na.rm = TRUE) * 0.002

# change name
ind = sapply(tree$tip.label, function(s){grep(s, name_exchange$X1)})
tree$tip.label = name_exchange$X2[ind]

tree_tb = tree %>% as_tibble()
internalnodes_nodeID = sapply(internalnodes, function(s){
  grep(s, tree_tb$label)
})

tree_p = ggtree(tree, size = 0.25) + 
  geom_tiplab(offset = 0.1, size = FontSize*0.5) +
  geom_point2(aes(subset=(node %in% internalnodes_nodeID)), shape=19, size=1,
              color='gold') +
  geom_nodelab(#aes(subset=(node %in% internalnodes_nodeID)),
    size = FontSize*0.5, hjust = 1.1, 
    # vjust = 1.5
    # nudge_x = -1,
    nudge_y = -0.4
  )+
  xlim_tree(x_limit) +
  # geom_text(aes(label=node), hjust=-.3) +
  # geom_rootedge(rootedge, color = 'black', size = 0.25) +
  # geom_tippoint(aes(color = group)) +
  # scale_color_manual(values=c("black"), na.value =  "black") + 
  # coord_flip() +
  theme_tree() + theme(legend.position = 'none',
                       plot.background = element_blank()); tree_p
tree_p = tree_p %>% rotate(15); tree_p

## prep data

wanted = c("Node12", "Node18", "Node19",
           "Node20", "Node21",
           "NECH", "COAL", "LECO", "MAVI", "PIFI",
           "Node4", "Node5", "Node8", "Node10",
           "Node13", "Node15", "CEOR", "PICH", "EMTR")
internalnodes = wanted[grepl("Node", wanted)]
wanted_tips = wanted[!grepl("Node", wanted)]

# wanted = c("Node12", "Node18", "Node19", "Node21")
# internalnodes = wanted[grepl("Node", wanted)]

output_14_selCol = output_14 %>% select(class, gene, Name, node, `Uncorrected P-value`, `Corrected P-value`)
output_13_selCol = output_13 %>% select(class, gene, Name, node, `Uncorrected P-value`, `Corrected P-value`) 

output_14_selFinal = output_14_selCol %>% 
  filter(Name %in% wanted)
unique(output_14_selFinal$Name)
length(unique(output_14_selFinal$gene)) # 94

# N13 tips tree
x13 = c("LMOD1", "SLC26A5")
x14 = "SLC52A3"
x15 = "LDHB"

output_13_sel1_tips = output_13_selCol %>%
  filter(gene %in% x13) %>%
  filter(Name %in% wanted_tips) 
output_13_sel1_internal = output_13_selCol %>%
  filter(gene %in% x13) %>%
  filter(node %in% c(23, 22, 21, 19, 20, 17, 16, 15)) %>% 
  arrange(gene, node) %>% # pause here to check what to rename
  mutate(Name = gsub("Node18", "Node20", Name), # 23
         Name = gsub("Node17", "Node19", Name), # 22
         Name = gsub("Node16", "Node18", Name)  # 21
         ) 
output_13_sel1 = output_13_sel1_tips %>% bind_rows(output_13_sel1_internal)
length(unique(output_13_sel1$gene)) # 2

output_13_sel2_tips = output_13_selCol %>%
  filter(gene %in% x14) %>%
  filter(Name %in% wanted_tips) 
output_13_sel2_internal = output_13_selCol %>%
  filter(gene %in% x14) %>%
  filter(node %in% c(24, 23, 22, 21, 19, 20, 18, 17, 16, 15)) %>% 
  arrange(gene, node) %>% # pause here to check what to rename
  mutate(Name = gsub("Node19", "Node21", Name), # 24
         Name = gsub("Node18", "Node20", Name), # 23
         Name = gsub("Node17", "Node19", Name), # 22
         Name = gsub("Node16", "Node18", Name)  # 21
  )
output_13_sel2 = output_13_sel2_tips %>% bind_rows(output_13_sel2_internal)
length(unique(output_13_sel2$gene)) # 1

output_13_sel3_tips = output_13_selCol %>%
  filter(gene %in% x15) %>%
  filter(Name %in% wanted_tips) 
output_13_sel3_internal = output_13_selCol %>%
  filter(gene %in% x15) %>%
  filter(node %in% c(22, 21, 19, 20, 18, 17, 16, 15)) %>% 
  arrange(gene, node) %>% # pause here to check what to rename
  mutate(Name = gsub("Node17", "Node19", Name), # 22 m
         Name = gsub("Node16", "Node18", Name), # 21 m
         Name = gsub("Node13", "Node15", Name), # 20 m
         Name = gsub("Node11", "Node13", Name), # 19 m
         Name = gsub("Node10", "Node12", Name), # 18 m
         Name = gsub("Node8", "Node10", Name) # 17 m
  )
output_13_sel3 = output_13_sel3_tips %>% bind_rows(output_13_sel3_internal)
length(unique(output_13_sel3$gene)) # 1

output_13_sel4_tips = output_13_selCol %>%
  filter(gene %in% x16) %>%
  filter(Name %in% wanted_tips) 
output_13_sel4_internal = output_13_selCol %>%
  filter(gene %in% x16) %>%
  filter(node %in% c(18, 17, 16, 15, 14)) %>% 
  arrange(gene, node) %>% # pause here to check what to rename
  mutate(#Name = gsub("Node17", "Node4", Name), # 14
         #Name = gsub("Node16", "Node5", Name), # 15
         #Name = gsub("Node13", "Node8", Name), # 16
         Name = gsub("Node14", "Node18", Name), # 17
         Name = gsub("Node17", "Node21", Name)  # 18
  )
output_13_sel4 = output_13_sel4_tips %>% bind_rows(output_13_sel4_internal)
length(unique(output_13_sel3$gene)) # 1

# combine
output_14_sel2 = output_14_selFinal %>% 
  bind_rows(output_13_sel1) %>% bind_rows(output_13_sel2) %>% 
  bind_rows(output_13_sel3) %>%  bind_rows(output_13_sel4) %>% 
  filter(grepl("^\\d", class)) %>%
  mutate(class = gsub("^\\d\\w{0,1}_", "", class)) %>%
  mutate(corrected5 = ifelse(`Corrected P-value` <=0.05, TRUE, NA)) 
unique(output_14_sel2$class)
dim(output_14_sel2) #  1913   7 > 1998    7

# make dummy data for not-existing nodes and get 'class'
gene_runthrough = unique(output_14_sel2$gene)
output_14_sel2_count = output_14_sel2 %>%
  group_by(gene, class) %>%
  mutate(n = n())
tmp_out = lapply(gene_runthrough, function(s){
  # s = "SI"
  print(s)
  tmp_sub = output_14_sel2_count %>% filter( gene == s)
  tmp_n = unique(tmp_sub$n)
  tmp_class = unique(tmp_sub$class)
  tmp_nodes = unique(tmp_sub$Name)
  if( tmp_n < length(internalnodes)){
    toadd = setdiff(wanted, tmp_nodes)
    # names(tmp_sub)
    df_out = tibble("Name" =  toadd, "gene" = s, 
                    "Uncorrected P-value" = 1,
                    "Corrected P-value" = 1,  
                    "class" = tmp_class, "corrected5" = NA )
  }}) %>% bind_rows()

output_mod = output_14_sel2 %>% bind_rows(tmp_out)
# output_mod$gene = sapply(output_mod$gene, function(s){unlist(strsplit(s, "_"))[1]})
dim(output_mod) #  1913 >1998   7
length(unique(output_mod$gene)) == length(target_working_files)


# change tip name
ind_tips = which(!grepl("^Node", output_mod$Name) & !grepl("FUFI", output_mod$Name))
ind = sapply(output_mod$Name[ind_tips], function(s){grep(s, name_exchange$X1, ignore.case = T)})
output_mod$Name[ind_tips] = name_exchange$X2[ind]

ind_tips = which(!grepl("^Node", wanted) & !grepl("FUFI", wanted))
ind = sapply(wanted[ind_tips], function(s){grep(s, name_exchange$X1, ignore.case = T)})
wanted[ind_tips] = name_exchange$X2[ind]

unique(output_mod$Name)
output_mod$Name = factor(output_mod$Name, levels = rev(wanted))
unique(output_mod$class)
class_order = c("Diet_Carbohydrates", "Diet_Protein", "Diet_Fat", "Diet", "Muscle", 
                "Vision", "Audition", "Telomere", "Sperm" )
output_mod$class = factor(output_mod$class, levels = class_order)
# output_mod$corrected5 = factor(output_mod$corrected5, levels = c(TRUE))

# groups
manakins = c("Neopelma chrysocephalum", "Corapipo altera", "Lepidothrix coronata", 
             "Manacus vitellinus", "Pipra filicauda")
tc = c("Cephalopterus ornatus", "Piprites chloris", "Empidonax traillii")
output_mod = output_mod %>% 
  mutate(group = case_when(Name %in% paste0("Node", c(12,18:21)) ~ "Manakin ancestors",
                           Name %in% manakins ~ "Manakin",
         Name %in% paste0("Node", 4:5) ~ "Passerine and oscine ancestors",
         Name %in% paste0("Node", c(8,10)) ~ "Suboscine ancestors",
         Name %in% paste0("Node", c(13,15)) ~ "Tyrant + cotinga ancestors",
         Name %in% tc ~ "Tyrant + cotinga"
         )) %>%
  mutate(group = fct_relevel(group, c("Manakin ancestors", "Manakin", 
                                      "Suboscine ancestors", "Tyrant + cotinga ancestors", 
                                      "Tyrant + cotinga",
                        "Passerine and oscine ancestors")))
unique(output_mod$group)


### color 2 
color_highp_tmp2 = 'white'
color_low_tmp2 = 'red3'

color_highp_tmp = 'chocolate2'
# color_highp_tmp = 'grey50'
color_midp_tmp = 'khaki1' 
color_tmp = 'red3'

theme_m2 = theme(legend.position = "right",
                 axis.text.y = element_text(size = AxisTxFontSizeSize_s),
                 axis.title.y = element_blank(),
                 axis.title.x = element_blank(), 
                 strip.text = element_text(size = AxisTxFontSizeSize_s, 
                                           margin = margin(2,0,2,0, "pt")),
                 legend.key.size = unit(8, 'points'), 
                 legend.background = element_blank(),
                 legend.key = element_blank(),
                 legend.title = element_text(size = AxisTxFontSizeSize_s),
                 legend.text = element_text(size = AxisTxFontSizeSize_s))

length(unique(output_mod$gene)) # 98

p = ggplot(output_mod, aes(y = Name, x = gene, fill = `Uncorrected P-value`)) +
  geom_tile(color = "grey90", size = 0.1) +
  # geom_tile(aes(color = corrected5, size = `Corrected P-value`)) +
  geom_tile(aes(color = corrected5), size = 0.2, height = 0.9, width = 0.9) +
  scale_fill_gradientn(colors = c(color_tmp, 'rosybrown1', color_midp_tmp, color_highp_tmp),
                       na.value = "transparent",
                       limits = c(0, 0.1), breaks = c(0, 0.05, 0.1)) +
  # scale_fill_gradient(low = color_low_tmp2, high = color_highp_tmp2, 
  #                     na.value = "transparent",
  #                     limits = c(0, 0.05)) +
  # scale_color_manual(values = c("black", "transparent"), na.value = "transparent") +
  scale_color_manual(values = c("black"), na.value = "transparent") +
  # scale_size_continuous(range = c(0.1, 0.5)) +
  scale_y_discrete(position = "right", expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(group~class, scales = 'free', space = 'free', switch = "x") +
  # facet_wrap(.~class, scales = 'free', strip.position = "bottom", nrow = 2, shrink = TRUE) +
  guides(colour = "none") +
  # facet_wrap(.~class, scales = 'free', strip.position = "bottom", nrow = 2, shrink = TRUE) +
  theme_m + theme_m2 + 
  theme(axis.text.x = element_text(size = AxisTxFontSizeSize_s, 
                                   angle = 40, hjust = 1),
        legend.key.size = unit(5, 'pt'),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()
        # panel.border = element_rect(color = "black", size = 0.1, fill = NA)
  ); p


# grey 0.05 to 0.1
output_mod2 = output_mod %>%
  mutate(unco = ifelse(`Uncorrected P-value` <= 0.1, `Uncorrected P-value`, NA))

p = ggplot(output_mod2, aes(y = Name, x = gene, fill = unco)) +
  geom_tile(color = "grey90", size = 0.1) +
  
  # geom_tile(aes(color = corrected5, size = `Corrected P-value`)) +
  geom_tile(aes(color = corrected5), size = 0.2, height = 0.9, width = 0.9) +
  # scale_fill_gradientn(colors = c(color_tmp, 'rosybrown1'),
  #                      na.value = "grey90",
  #                      limits = c(0, 0.05), breaks = c(0, 0.05, 0.1)) +
  scale_fill_gradient2(low = color_low_tmp2, mid = color_highp_tmp2, high = 'grey80',
                       midpoint = 0.05,
                       na.value = "transparent",
                       limits = c(0, 0.1), breaks = c(0, 0.05, 0.1)) +
  # scale_color_manual(values = c("black", "transparent"), na.value = "transparent") +
  scale_color_manual(values = c("black"), na.value = "transparent") +
  # scale_size_continuous(range = c(0.1, 0.5)) +
  scale_y_discrete(position = "right", expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(group~class, scales = 'free', space = 'free', switch = "x") +
  # facet_wrap(.~class, scales = 'free', strip.position = "bottom", nrow = 2, shrink = TRUE) +
  guides(colour = "none") +
  # facet_wrap(.~class, scales = 'free', strip.position = "bottom", nrow = 2, shrink = TRUE) +
  theme_m + theme_m2 +
  theme(axis.text.x = element_text(size = AxisTxFontSizeSize_s,
                                   angle = 40, hjust = 1),
        legend.key.size = unit(5, 'pt'),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()
        # panel.border = element_rect(color = "black", size = 0.1, fill = NA)
  ); p

# color 0.5 to 0.1 gray -ver3 [using] ====
color_highp_tmp2 = 'mistyrose'
output_mod2 = output_mod %>%
  mutate(unco = ifelse(`Uncorrected P-value` < 0.1, `Uncorrected P-value`, NA)) #%>% 
  # tidyr::complete(nesting(gene, Name, class, group))
output_mod2$gene = factor(output_mod2$gene, levels = genelist$Gene)

x = output_mod2 %>% filter(`Uncorrected P-value` < 0.1 & `Uncorrected P-value` >0.05 ) %>% 
  arragne(class)
nrow(x) # 74
  
p = ggplot(output_mod2, aes(y = Name, x = gene, fill = unco)) +
  geom_tile(color = "grey90", fill = 'transparent', size = 0.1) +
  # geom_tile(aes(color = corrected5, size = `Corrected P-value`)) +
  geom_tile(aes(color = corrected5), size = 0.2, height = 0.9, width = 0.9) +
  # scale_fill_gradientn(colors = c(color_tmp, 'rosybrown1'),
  # na.value = "grey90",
  # limits = c(0, 0.1), breaks = c(0, 0.05, 0.1)) +
  scale_fill_gradient2(low = color_low_tmp2, mid = color_highp_tmp2, high = 'grey20',
                       midpoint = 0.05,
                       na.value = "transparent",
                       limits = c(0, 0.1), breaks = c(0, 0.05, 0.1)) +
  # scale_color_manual(values = c("black", "transparent"), na.value = "transparent") +
  scale_color_manual(values = c("black"), na.value = "transparent") +
  # scale_size_continuous(range = c(0.1, 0.5)) +
  scale_y_discrete(position = "right", expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  facet_grid(group~class, scales = 'free', space = 'free', switch = "x") +
  # facet_wrap(.~class, scales = 'free', strip.position = "bottom", nrow = 2, shrink = TRUE) +
  guides(colour = "none") +
  # facet_wrap(.~class, scales = 'free', strip.position = "bottom", nrow = 2, shrink = TRUE) +
  theme_m + theme_m2 +
  theme(axis.text.x = element_text(size = AxisTxFontSizeSize_s,
                                   angle = 40, hjust = 1),
        legend.key.size = unit(5, 'pt'),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()
        # panel.border = element_rect(color = "black", size = 0.1, fill = NA)
  ); p
# x = output_mod %>% filter(grepl("SI", gene))

n_row = 10
n_empty_row = 0
n_col = 8
n_col_tree = 1
lay = matrix(c( rep(c(rep( 1, n_col_tree), rep(2, (n_col-n_col_tree) )), (n_row-n_empty_row) ),
                rep(c(rep(NA, n_col_tree), rep(2, (n_col-n_col_tree) )), n_empty_row)), 
             nrow = n_row, byrow = TRUE); lay
p_out = arrangeGrob(grobs = list(tree_p, p), layout_matrix = lay)
# p_out = arrangeGrob(grobs = list(tree_p, p), nrow = 1, widths = c(2,3))
# dev.off()
# grid.draw(p_out)

graph_path = paste0(project_path, "/Summary_14sp_allnodes_14tips_", Sys.Date(), "_c2g.pdf")
pdf(graph_path, width = Width_HalfCol*2.95, height = Width_HalfCol*0.6, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
grid.draw(p_out)
dev.off()
openfile(graph_path)

