#!/usr/bin/env Rscript
# 03_module_discovery.R
# Detects functional modules within the rewired network using Louvain Clustering.

library(igraph)
library(dplyr)

# --- USER CONFIGURATION ---
REWIRED_LIST <- "candidates_rewired_proteins.csv"      # Output from Step 02
INTERACTOME_FILE <- "path/to/reference_interactome.csv" # Original PPI network

# --- 1. Construct the Rewired Sub-Network ---
cat("Constructing rewired sub-network...\n")
rewired_prots <- read.csv(REWIRED_LIST)$Protein
full_ppi <- read.csv(INTERACTOME_FILE)

# Keep edges only if BOTH nodes are in the "Rewired" list
# (Strict modularity: looking for dense cores of rewired proteins)
sub_network_edges <- full_ppi %>%
    filter(Protein_A %in% rewired_prots & Protein_B %in% rewired_prots)

g <- graph_from_data_frame(sub_network_edges, directed = FALSE)

# --- 2. Remove Disconnected Components ---
# We focus on the Largest Connected Component (LCC) or significant subgraphs
components <- components(g)
annotated_nodes <- names(components$membership)
to_keep <- annotated_nodes[components$membership == which.max(components$csize)]
g_lcc <- induced_subgraph(g, to_keep)

cat(sprintf("Rewired Network Stats:\n  Nodes: %d\n  Edges: %d\n", vcount(g_lcc), ecount(g_lcc)))

# --- 3. Louvain Community Detection ---
cat("Running Louvain Modularity...\n")
set.seed(123)
cluster_res <- cluster_louvain(g_lcc)

# Assign module membership
membership_df <- data.frame(
    Protein = names(membership(cluster_res)),
    Module = as.numeric(membership(cluster_res))
)

# --- 4. Module Summary ---
module_stats <- membership_df %>%
    group_by(Module) %>%
    summarise(Size = n()) %>%
    arrange(desc(Size)) %>%
    filter(Size >= 5) # Filter small artifacts

print(head(module_stats))

# --- 5. Export ---
write.csv(membership_df, "rewired_modules.csv", row.names = FALSE)
cat("Saved module assignments. Ready for Functional Enrichment analysis.\n")
