library(tidyverse)
library(tidygraph)
library(pracma)

all_nodes_file <-"/mnt/dv/wid/projects3/Roy-enhancer-promoter/RMEC/Graph_Diffusion/Roadmap_RegulatoryVariation/RShiny/igraphexample4/cd56_primary_cells_k7_L1.txt"
edge_list_file <- "/mnt/dv/wid/projects3/Roy-enhancer-promoter/RMEC/Graph_Diffusion/Roadmap_RegulatoryVariation/RShiny/igraphexample4/subnetwork_cd56_primary_cells.txt"

Module <- read_tsv(all_nodes_file, col_names=TRUE)
Module <- Module %>% 
  group_by(ClusterID) %>%
  summarise(Gene = list(Gene)) %>%
  rename("module" = "ClusterID") %>%
  rename("gene_list"="Gene")

module_ids=0:6

TGS <- read_tsv(all_nodes_file, col_names=TRUE)
transitioning_sets = unique(TGS$Transitioning)
TGS <- TGS %>%
  group_by(Transitioning) %>%
  summarise(Gene = list(Gene)) %>%
  rename("gene_list"="Gene")


nodes <- read_tsv(all_nodes_file, col_names = TRUE) %>%
	rename("feature"= "Gene") %>%
	rowid_to_column("id")

routes <- read_tsv(edge_list_file, col_names = TRUE)

### Generate Edge Set
edges <- routes %>%
  left_join(nodes, by = c("Source" = "feature")) %>%
  rename(from = id)

edges <- edges %>%
  left_join(nodes, by = c("Target" = "feature")) %>%
  rename(to = id)

edges <- edges %>%
  select('from', 'to')

### Make Tidy Graph Structure 
Net <- tbl_graph(nodes= nodes, edges = edges)
Net <- Net %>%
  convert(to_undirected) %>%
  mutate(neighbors = map_local(order = 1, .f = function(neighborhood, node, ...) {
         as_tibble(neighborhood, active = 'nodes')$feature
      }))
degree_v <-  Net %N>% as_tibble() %>% rowwise() %>% summarize(length(neighbors)) -1
Net <- Net %>% mutate(degree = degree_v$`length(neighbors)`)

genes <- unique(Net %N>% as_tibble() %>% pull(feature))

save(list = c("Net","nodes","TGS","transitioning_sets","Module","module_ids","genes"), file = "netdata.Rdata")
