library(tidyverse)
library(tidygraph)
library(pracma)


#all_nodes_file <- "/mnt/dv/wid/projects7/Roy-Aspergillus/Results/RnaSeq/NcaResults/Output_20211122111849/Lambda_0100/Merlinp_inputs/net1_nodes.txt"
#edge_list_file <- "/mnt/dv/wid/projects7/Roy-Aspergillus/Results/RnaSeq/MerlinpResults/Afum_RnaSeq_results/Condor_results/PostBatchCorr_I02/Lambda_0100/output_net_0_8.txt"
#module2gene_file <- "/mnt/dv/wid/projects7/Roy-Aspergillus/Results/RnaSeq/MerlinpResults/Afum_RnaSeq_results/Condor_results/PostBatchCorr_I02/Lambda_0100/consensus_module_0_3_geneset.txt"
#module_file <- "/mnt/dv/wid/projects7/Roy-Aspergillus/Results/RnaSeq/MerlinpResults/Afum_RnaSeq_results/Condor_results/PostBatchCorr_I02/Lambda_0100/consensus_module_0_3_geneset_enrichAnalyzer.txt"
#go_file = "/mnt/dv/wid/projects7/Roy-Aspergillus/Data/GeneOntology/GO_enrichAnalyzer_idx/afumgotermap.txt"
#regulator_enrich_file <- "/mnt/dv/wid/projects7/Roy-Aspergillus/Results/RnaSeq/MerlinpResults/Afum_RnaSeq_results/Condor_results/PostBatchCorr_I02/Lambda_0100/Enrichments/merlin.0_8.0_3_details.txt"
#go_enrich_file <- "/mnt/dv/projects7/Roy-Aspergillus/Results/RnaSeq/MerlinpResults/Afum_RnaSeq_results/Condor_results/PostBatchCorr_I02/Lambda_0100/Enrichments/go.0_3_details.txt"


################### Make R Data Files *Only need to run once###################
makePostProcessDataStruct <- function (all_nodes_file, edge_list_file,
                          module2gene_file, go_file, module_file, 
                          regulator_enrich_file, go_enrich_file)
{
###### Generate labeled Node Set
genes2modules <- read_tsv(module2gene_file, c("feature", "module"))
genes2modules <- genes2modules %>% 
  group_by(module) %>% 
  summarize(count = n()) %>% 
  right_join(genes2modules, by = "module") %>%
  ungroup() %>% 
  filter(count > 4) %>%
  select(feature, module)

go <- read_tsv(go_file, col_names = TRUE) %>% 
  rename("feature"= "GeneName") %>%
  rename("go" = "GOTerm") %>%
  select(feature, go) %>%
  group_by(feature) %>%
  chop(go) %>%
  ungroup()

nodes <- read_tsv(file = all_nodes_file, col_names = "feature") %>% 
          rowid_to_column("id") %>%
          left_join(genes2modules) %>%
          left_join(go)


  
### Generate Edge Set
routes <- read_tsv(edge_list_file, c("source", "target", "weight"))
edges <- routes %>% 
  left_join(nodes, by = c("source" = "feature")) %>% 
  rename(from = id)

edges <- edges %>% 
  left_join(nodes, by = c("target" = "feature")) %>% 
  rename(to = id)

edges <- edges %>%
  select('from', 'to', 'weight')

### Make Tidy Graph Structure 
Net <- tbl_graph(nodes= nodes, edges = edges)
Net <- Net %>%
  convert(to_undirected) %>%
  mutate(neighbors = map_local(order = 1, .f = function(neighborhood, node, ...) {
         as_tibble(neighborhood, active = 'nodes')$feature
      }))
degree_v <-  Net %N>% as_tibble() %>% rowwise() %>% summarize(length(neighbors)) -1
Net <- Net %>% mutate(degree = degree_v$`length(neighbors)`)

### Generate Module Structure 
Module <- read_tsv(module_file, c("module", "gene_list"))  
Module$gene_list <- str_split(Module$gene_list, '#')
go_enrich <- read_tsv(go_enrich_file, col_names = c("module", "go", 
                      "go_pvalue", "go_correct_pvalue", "go_num_tot_genes", 
                      "num_tot_go_genes", "num_module_genes", 
                      "num_module_go_genes", "go_enrichment", "go_genes"))
go_enrich$go_genes <- str_split(go_enrich$go_genes, ';')
go_enrich <- nest(go_enrich, GO=c("go", "go_pvalue", "go_correct_pvalue", 
                                  "go_num_tot_genes", "num_tot_go_genes", 
                                  "num_module_genes", "num_module_go_genes", 
                                  "go_enrichment", "go_genes"))

merlin_enrich <- read_tsv(regulator_enrich_file, col_names = c("module", 
                          "regulator", "reg_pvalue", "reg_correct_pvalue", 
                          "reg_num_tot_genes", "num_tot_regulator_genes", 
                          "reg_num_module_genes", "num_module_regulator_genes", 
                          "regulator_enrichment", "reg_target_genes"))
merlin_enrich$reg_target_genes <- str_split(merlin_enrich$reg_target_genes, ';')
merlin_enrich <- nest(merlin_enrich, "regulators"= c("regulator", "reg_pvalue", 
                          "reg_correct_pvalue", "reg_num_tot_genes", 
                          "num_tot_regulator_genes", "reg_num_module_genes", 
                          "num_module_regulator_genes", "regulator_enrichment", 
                          "reg_target_genes") )

Module <- Module %>% 
  left_join(go_enrich, by = "module") %>%
  left_join(merlin_enrich, by = "module")

Module$module <- strtoi(str_replace_all(Module$module, "Cluster", ""))

module_ids <- Module %>% pull(module)
enrich_2_module <- Module %>% select(module, GO) %>% unnest(GO) %>% select(go, module) %>% group_by(go) %>% chop(module)

enriched_go_terms <- enrich_2_module %>%
pull(go)

 
genes <- unique(Net %N>% as_tibble() %>% pull(feature))

save(list = c("Net", "Module", "enriched_go_terms", "module_ids","enrich_2_module", "genes"), file = "net_data.Rdata")


return(list(Net, Module, enriched_go_terms, module_ids, enrich_2_module, genes)) 
}


makeLaplacian <- function(Net){
  degree_vec <- Net %N>% pull(degree)
  D <- diag(degree_vec)
  adj = matrix(0, length(degree_vec), length(degree_vec))
  Edges <- Net %E>% as_tibble()
  from <- Edges %>% pull(from)
  to <- Edges %>% pull(to)
  for(i in 1:length(from)){
    adj[from[i], to[i]]= 1
    adj[to[i], from[i]]= 1
  }
  L <- D - adj;
  save(list = c("Net", "Module", "enriched_go_terms", "module_ids","enrich_2_module", "genes", "L"), file = "net_data.Rdata")
  return(L)
}


MakeKernel <- function(L, lambda){
  num_nodes <- size(L)[1]
  I <- eye(num_nodes)
  inside <- I + lambda*L
  kernel <- inv(inside)
  return(kernel)
}

###############################################################################
##################### Search Functions ########################################
###############################################################################
searchForModule<- function(Module, moduleID, celltype){
  modIdx<-which(Module$module == moduleID & Module$Cell == celltype)
  gene_in_module <- Module$gene_list[[modIdx]]
  regulators <- Module$regulators[[modIdx]]$regulator
  mod_genes<- c(gene_in_module, regulators)
return(mod_genes)
}

searchForGeneset<- function(TGS, gene_set,celltype){
  modIdx<-which(TGS$Transitioning == gene_set & TGS$Cell == celltype)
  gene_in_module <- TGS$gene_list[[modIdx]]
return(gene_in_module)
}

searchForGene <-function(Net,Module, gene){
  moduleIDs <- Net %N>%
    filter(feature %in% gene) %>%
    pull(module)
  moduleIDs <- unique(moduleIDs[!is.na(moduleIDs)])
  genes <-list()
  for (ID in moduleIDs) {
    new_genes <- searchForModule(Module, ID )
    genes <- append(unlist(genes), unlist(new_genes))
  }
  
  neighbors <- Net %N>%
    filter(feature %in% gene) %>%
    pull(neighbors)
  
  genes <- append(unlist(genes), unlist(unlist(neighbors)))
  print(unique(genes))
  return(unique(genes))
}

searchForGeneList <-function(Net,Module, gene_list){
  moduleIDs <- Net %N>%
    filter(str_detect(feature, paste(unlist(gene_list), collapse="|"))) %>%
    pull(module)
  moduleIDs <- unique(moduleIDs[!is.na(moduleIDs)])
  genes <-list()
  for (ID in moduleIDs) {
    new_genes <- searchForModule(Module, ID )
    genes <- append(unlist(genes), unlist(new_genes))
  }

  neighbors <- Net %N>%
    filter(str_detect(feature, paste(unlist(gene_list), collapse="|"))) %>%
    pull(neighbors)

  genes <- append(unlist(genes), unlist(unlist(neighbors)))
  print(unique(genes))
return(unique(genes))
}

computeEnrichment <- function(Module, gl, num_genes){
    Module <- Module %>% rowwise() %>% 
      mutate(m_size = length(gene_list)) %>%
      mutate(intersect_size = length(intersect(gl, gene_list))) %>%
      mutate(enrich_pval = phyper(intersect_size, length(gl), num_genes, m_size, lower.tail=FALSE)) %>%
      mutate(corrected_pval = ifelse(enrich_pval * dim(Module)[1] < 1, enrich_pval * dim(Module)[1], 1))
    return(Module)
}

########################## Node Diffusion #####################################

generateScoreVector <- function(Net, gene_list){
  nodes <- Net %N>% as_tibble() %>% select(feature)
  scores <- left_join(nodes, Net %N>% as_tibble() %>% 
    filter(str_detect(feature, paste(unlist(gene_list), collapse="|"))) %>%
    select(feature) %>% 
    mutate(score = 100)) %>% 
    mutate (score = replace_na(score, 0))
  return(scores$score)
}


loadScoreVector <- function(Net, score_data) {
  nodes <- Net %N>% as_tibble() %>% mutate(score = 0)
  for( i in 1:length(score_data[[1]])){ 
    nodes <- nodes %>%
      mutate(score = replace(score, 
             str_detect(feature, score_data %>% slice(i) %>% pull(feature)),
             score_data %>% slice(i) %>% pull(score)))
  }
  return(nodes$score)
}

computeDiffusionScore <- function (Net, score_data, kernel){
  score <- loadScoreVector(Net, score_data)
  diff_score <- kernel %*% score
  Net <- Net %N>% mutate("score" = as.vector(diff_score))
  return(Net)
}



###############################################################################
##################### Subgraph Functions ######################################
###############################################################################
moduleSubgraph <- function(Net, Module, module_id,celltype){
  mod_genes <- searchForModule(Module, module_id,celltype)
  sub_graph <-induceSubraph(Net, mod_genes)
  return(sub_graph)
}

genesetSubgraph <- function(Net, TGS, gene_set,celltype){
  mod_genes <- searchForGeneset(TGS, gene_set,celltype)
  sub_graph <-induceSubraph(Net, mod_genes)
  return(sub_graph)
}

geneSubgraph <- function(Net, Module, gene){
  list_genes <- searchForGene(Net, Module, gene)
  sub_graph <- induceSubraph(Net, list_genes)
  return(sub_graph)
}

geneListSubgraph <- function(Net, Module, gene_list){
  list_genes <- searchForGeneList(Net, Module, gene_list)
  sub_graph <- induceSubraph(Net, list_genes)
  return(sub_graph)
}

goSubgraph <- function(Net, Module, enrich_2_module, go_term){
  modules_list<-unlist(enrich_2_module %>%
                         filter(go == go_term) %>%
                         pull(module))
  mod_genes = list()
  for(i in 1:length(modules_list)){
    mod_genes<- append(unlist(mod_genes), unlist(searchForModule(Module, modules_list[i])))
  }
  sub_graph <- induceSubraph(Net, mod_genes)
  return(sub_graph)
}

diffScoreSubgraph <- function(Net, percentile){
  val <- quantile(Net %N>% pull(score), probs = percentile)
  genes <- Net %N>% filter(score >= val) %>% pull(feature)
  subgraph <- induceSubraph(Net, genes)
}


induceSubraph <- function(Net, list){
  sub_graph <- Net %N>%
    convert(to_subgraph, feature %in% list)
  return(sub_graph)
}


graph2NodeEdgeTables <- function(Net,celltype){
  graph_nodes <- Net %N>%
	mutate(id = row_number()-1)
  graph_nodes <- graph_nodes %N>%
    filter(Cell == celltype)

  graph_edges <- graph_nodes %E>%
    filter(Cell == celltype) %>%
    as_tibble() %>%
    mutate(from =from - 1, to = to - 1)

  graph_nodes <- graph_nodes %N>%
     mutate(id = id -min(id)) %>%
     as_tibble()

  return( list(graph_nodes, graph_edges))
}


###############################################################################
##################### Stiener Tree Construction ###############################
###############################################################################
getDistMatrix <- function (Net, gene_list){
gene_id <- Net %>% 
  filter(feature %in% gene_list) %>%
  pull(id)

gene_name <-  Net %>% 
  filter(feature %in% gene_list) %>%
  pull(feature)

dist_matrix <- Net %N>%
  as_tibble() %>%
  select(feature)

for(idx in 1:length(gene_id)){
  root <- gene_id[idx]
  name <- gene_name[idx]
  dist_matrix <- Net %N>% 
  mutate(dist = bfs_dist(root)) %>%
  as_tibble() %>%
  select(feature, dist) %>%
  rename(!!name := dist) %>%
  right_join(dist_matrix, by="feature")  
}

dist_matrix <- dist_matrix %>% replace(. ==0 , length(dist_matrix$feature)+1)
return(dist_matrix)
}


buildStienerTrees <- function(Net, gene_list){
  dist_matrix <- getDistMatrix(Net, gene_list)
  gene_names <- colnames(dist_matrix)[2:length(dist_matrix)]
  dist2graph <- tibble(gene_names)
  dist2graph <- dist2graph %>%
    mutate(Dist = Inf) %>%
    mutate(Closest = "")
  
  
## Find Closest Nodes
  stiener_tree=tbl_graph()
  restrict_dist_matrix <- dist_matrix %>%
    filter(feature %in% gene_names) 
  
  for(i in 2:length(restrict_dist_matrix)){
    gene <- colnames(restrict_dist_matrix[i])
    idx <-which(dist2graph$gene_names == gene)
    if(all(is.na(restrict_dist_matrix[i]))){
      dist2graph$Closest[idx[1]] <- gene
      dist2graph$Dist[idx[1]] <- length(dist_matrix$feature)+1
    }else{
      dist <- min(restrict_dist_matrix[i], na.rm = TRUE)
      match_idx <- which(restrict_dist_matrix[i] == min(restrict_dist_matrix[i], na.rm = TRUE))
      closest <- restrict_dist_matrix$feature[match_idx[1]]
      dist2graph$Dist[idx[1]] <- dist
      dist2graph$Closest[idx[1]] <- closest
    }
  }
  
## Select minimum length path   
  path_2_add<-dist2graph %>% 
    arrange(Dist) %>%
    slice(1) %>%
    select("gene_names", "Closest")
  path_2_add<-c(path_2_add$gene_names,path_2_add$Closest)
  
## Find id of node_2_add  
  nodes_2_add<- Net %N>% 
    filter(feature %in% path_2_add) %>%
    pull(id)
    

## Find path and initialize Stiener tree. 
  stiener_tree <- Net %N>%
    convert(to_shortest_path, nodes_2_add[1], nodes_2_add[2])
  
  nodes <- stiener_tree %N>%
    pull(feature)
  
## Prune search matrix & dist2graph    
  dist_matrix <- dist_matrix %>% 
    select(-nodes)
  dist2graph <- dist2graph %>%
    filter(gene_names %in% setdiff(gene_names, nodes))
  
    
## Now we to find distance to nodes in graph.
  while(length(intersect(nodes, gene_names)) < length(gene_names)){
    restrict_dist_matrix <- dist_matrix %>%
      filter(feature %in% nodes)

    for(i in 2:length(restrict_dist_matrix)){
      gene <- colnames(restrict_dist_matrix[i])
      idx <-which(dist2graph$gene_names == gene)
      if(all(is.na(restrict_dist_matrix[i]))){
        dist2graph$Closest[idx[1]] <- gene
        dist2graph$Dist[idx[1]] <- length(dist_matrix$feature)+1
      }else{
        dist <- min(restrict_dist_matrix[i], na.rm = TRUE)
        match_idx <- which(restrict_dist_matrix[i] == min(restrict_dist_matrix[i], na.rm = TRUE))
        closest <- restrict_dist_matrix$feature[match_idx[1]]
        dist2graph$Dist[idx[1]] <- dist
        dist2graph$Closest[idx[1]] <- closest
      }
    }
  
    path_2_add<-dist2graph %>% 
      arrange(Dist) %>%
      slice(1) %>%
      select("gene_names", "Closest")
    path_2_add<-c(path_2_add$gene_names,path_2_add$Closest)
  
    nodes_2_add<- Net %N>% 
      filter(feature %in% path_2_add) %>%
      pull(id)
    
    if( length(nodes_2_add) == 1){
      node_info<-Net %N>%
        filter(id %in% nodes_2_add) %>%
        as_tibble()
      stiener_tree<- stiener_tree %>%
        bind_nodes(node_info)
    }else{
      Path <- Net %N>%
        convert(to_shortest_path, nodes_2_add[1], nodes_2_add[2])
       
        
    
      stiener_tree <- stiener_tree %>% 
        graph_join(Path)
    }
    
    nodes <- stiener_tree %N>%
      pull(feature)
    
    ## Prune search matrix & dist2graph    
    dist_matrix <- dist_matrix %>% 
      select(setdiff(colnames(dist_matrix), nodes))
    dist2graph <- dist2graph %>%
      filter(gene_names %in% setdiff(gene_names, nodes))
  }
return(stiener_tree)
}

###############################################################################
######################### Print Functions #####################################
###############################################################################

printSNPInfo <- function(SNPs, node_name,cell) {
if(length(node_name)==0){
    return(" ")
  }else{
      module_info <- SNPs %>%
      filter(Gene == node_name & Cell == cell)
      if(nrow(module_info)>0){
      	text <- sprintf("Gene Name:%s<br/>SNPs:%s", node_name,paste(unlist(module_info$SNP), collapse = ','))
      }else{
         text <- sprintf("Gene Name:%s<br/>SNPs: No SNPs", node_name)
      }
    return(text)
  }
}

printNodeInfo <- function(Net, node_name){
  if(is.na(node_name)){
    return(" ")
  }else{
  node_data <- Net %N>%
    filter(feature == node_name) %>%
    as_tibble()
  node_id <- node_data %>% pull(id)
  
  neighbors <- setdiff(Net %>%
                          convert(to_local_neighborhood, node_id) %N>%
                          as_tibble() %>%
                          pull(feature), node_name)
  
  text <- sprintf("Node Name: %s<br/>Module: %d<br/>GO Terms:  %s<br/>Neighbors: %s", 
                  node_name, node_data$module, 
                  paste(unlist(node_data$go), collapse = ', '), 
                  paste(unlist(neighbors), collapse = ', '))  
  return(text)
  }
}

printModuleInfo <- function(Module, module_id, gene_list, genes){
  if(is.na(module_id)){
    return(" ")
  }else{
    
    if(!is_empty(gene_list)){
      Module<-computeEnrichment(Module, gene_list, length(genes))
    }
    
    module_info <- Module %>%
      filter(module == module_id)
    
    module_regulators <- module_info %>%
      select(regulators) %>%
      unnest(regulators)
    if(nrow(module_regulators) == 0){
      regulator_text <- "No enriched regulators"
    }else{
      regulator<- module_regulators$regulator[1]
      p_value <- module_regulators$reg_correct_pvalue[1]
      target <- module_regulators$reg_target_genes[[1]]
      regulator_text <-sprintf ("Regulator Name: %s &emsp;P-Value: %.02e<br/> Target Genes: %s<br/>", regulator, p_value, paste(unlist(target), collapse = ', '))
      if(length(module_regulators$regulator) > 1){
        for(i in 2:length(module_regulators$regulator)){
          regulator<- module_regulators$regulator[i]
          p_value <- module_regulators$reg_correct_pvalue[i]
          target <- module_regulators$reg_target_genes[[i]]
          regulator_text <-sprintf ("%sRegulator Name: %s &emsp;P-Value: %.02e<br/> Target Genes: %s<br/>", regulator_text, regulator, p_value, paste(unlist(target), collapse = ', '))
        }
      }
    }
    
    module_go <- module_info %>% 
      select(GO) %>%
      unnest(GO) 
    
    if(nrow(module_go) == 0){
      go_text <- "No enriched GO terms"
    }else{
      go <- module_go$go[1]
      p_value <- module_go$go_correct_pvalue[1]
      go_genes <- module_go$go_genes[[1]]
      go_text <-sprintf ("Go Term: %s &emsp;P-Value: %.02e<br/> GO Genes: %s<br/>", go, p_value, paste(unlist(go_genes), collapse = ', '))
      if(length(module_go$go) >1 ){
        for(i in 2:length(module_go$go)){
          go <- module_go$go[i]
          p_value <- module_go$go_correct_pvalue[i]
          go_genes <- module_go$go_genes[[i]]
          go_text <-sprintf ("%sGo Term: %s &emsp;P-Value: %.02e<br/> GO Genes: %s<br/>", go_text, go, p_value, paste(unlist(go_genes), collapse = ', '))
        }
      }
    }
    
    if(is_empty(gene_list)){
      text <- sprintf("Module Name: %s<br/>Module Genes: %s<br/><br/>Enriched Regulators:<br/>%s<br/><br/>Enriched GO:<br/>%s <br/><br/>", 
                      module_id, paste(unlist(module_info$gene_list), collapse = ', '), 
                      regulator_text, 
                      go_text)  
      
    }else{
      text <- sprintf("Module Name: %s<br/>module enrichment p-value: %.02e<br/>module corrected p-value: %.02e<br/>Genes from gene list: %s<br/>Module Genes: %s<br/><br/>Enriched Regulators:<br/>%s<br/><br/>Enriched GO: %s <br/><br/>", 
                      module_id, 
                      module_info$enrich_pval, 
                      module_info$corrected_pval, 
                      paste(ifelse( length(intersect(module_info$gene_list[[1]], gene_list)) > 0, unlist(intersect(module_info$gene_list[[1]], gene_list)), ""), collapse = ', '),
                      paste(unlist(module_info$gene_list), collapse = ', '),
                      regulator_text, 
                      go_text)  
    }
    return(text)
  }
}

printAllModuleInfo <- function(SubNet, Module, gene_list, genes){
 if(!is_empty(gene_list)){
   text_info<-""
   unique_modules <-unique(SubNet %N>% 
                             as_tibble() %>%
                             pull(module))
   
   for(id in unique_modules){
     text_info<- sprintf('%s %s', text_info, printModuleInfo(Module, id, gene_list, genes))
   }
 }else{
  text_info<-""
  unique_modules <-unique(SubNet %N>% 
   as_tibble() %>%
   pull(module))
  for(id in unique_modules){
   text_info<- sprintf('%s %s', text_info, printModuleInfo(Module, id, list()))
  }
 }
 return(text_info)
}

getModuleID <- function(Net, node_name){
  module_id <- Net %N>%
    filter(feature == node_name) %>%
    as_tibble() %>% 
    pull(module)
  print(module_id)
  return(module_id)
}
