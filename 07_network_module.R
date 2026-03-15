#!/usr/bin/env Rscript
# 07_network_module.R - Network generation + Louvain modules + Cytoscape export
#
# Louvain community detection: igraph::cluster_louvain() with default resolution
# parameter (gamma = 1). Seed set via set.seed(42) in 00_setup.R for reproducibility.
#
# Module enrichment uses the same hypergeometric test as 06_enrichment.R:
# phyper(observed - 1, m, n, k, lower.tail = FALSE) computes P(X >= observed),
# since phyper(q, ..., lower.tail = FALSE) = P(X > q).
local({
    args <- commandArgs(trailingOnly = FALSE)
    f <- grep("--file=", args, value = TRUE)
    d <- if (length(f)) dirname(normalizePath(sub("--file=", "", f))) else getwd()
    source(file.path(d, "00_setup.R"))
})

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(igraph)
    library(ggraph)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(tidyr)
    library(svglite)
    library(openxlsx)
    library(stringr)
})

# --- Parameters ---
MIN_MODULE_SIZE      <- 3
TOP_MODULES          <- 20
TOP_HUBS_PER_MODULE  <- 5
TOP_NETWORK_HUBS     <- 10
module_colors <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3"),
                   brewer.pal(8, "Paired"))

# --- Helper functions ---

get_gene_description <- function(gene_id, aliases_df) {
    matches <- aliases_df %>%
        filter(alias == gene_id) %>%
        pull(string_id) %>%
        unique()
    if (length(matches) == 0) return("Unknown")
    desc <- aliases_df %>%
        filter(string_id %in% matches,
               source == "UniProt_DE_SubName_Full") %>%
        pull(alias) %>%
        first()
    if (is.na(desc)) {
        desc <- aliases_df %>%
            filter(string_id %in% matches,
                   source == "KEGG_PRODUCT") %>%
            pull(alias) %>%
            first()
    }
    if (is.na(desc)) return("Uncharacterized")
    desc
}

create_network <- function(interactions_df, rewired_proteins, name,
                           top_hubs = TOP_NETWORK_HUBS) {
    if (nrow(interactions_df) == 0) return(NULL)

    edges <- interactions_df %>%
        select(from = ProteinA, to = ProteinB, q_value, State) %>%
        mutate(
            interaction_type = case_when(
                q_value <= 0.1 ~ "Repressed",
                q_value >= 0.9 ~ "Activated",
                TRUE           ~ "Neutral"
            ),
            edge_color = case_when(
                interaction_type == "Repressed" ~ COLOR_REPRESSED,
                interaction_type == "Activated" ~ COLOR_ACTIVATED,
                TRUE                            ~ "#CCCCCC"
            )
        ) %>%
        filter(interaction_type != "Neutral")

    if (nrow(edges) == 0) return(NULL)

    g <- graph_from_data_frame(edges, directed = FALSE)

    V(g)$is_rewired <- V(g)$name %in% rewired_proteins
    V(g)$degree     <- degree(g)

    rewired_degrees <- data.frame(
        name       = V(g)$name,
        degree     = V(g)$degree,
        is_rewired = V(g)$is_rewired,
        stringsAsFactors = FALSE
    ) %>%
        filter(is_rewired) %>%
        arrange(desc(degree)) %>%
        head(top_hubs)

    hub_proteins   <- rewired_degrees$name
    V(g)$is_hub    <- V(g)$name %in% hub_proteins
    V(g)$node_color <- ifelse(V(g)$is_hub, COLOR_HUB,
                       ifelse(V(g)$is_rewired, "#4DAF4A", COLOR_REGULAR))
    V(g)$node_size  <- ifelse(V(g)$is_hub, 15,
                       ifelse(V(g)$is_rewired, 8, 4))
    V(g)$label      <- ifelse(V(g)$is_hub, V(g)$name, "")

    E(g)$color            <- edges$edge_color
    E(g)$interaction_type <- edges$interaction_type

    list(graph = g, edges = edges, hubs = rewired_degrees, name = name)
}

export_to_cytoscape_format <- function(network_obj, aliases_df) {
    if (is.null(network_obj)) return(invisible(NULL))
    g    <- network_obj$graph
    name <- network_obj$name

    # GraphML
    write_graph(g, file.path(NETWORK_DIR, paste0(name, "_network.graphml")),
                format = "graphml")

    # Node table
    node_table <- data.frame(
        id         = V(g)$name,
        is_rewired = V(g)$is_rewired,
        is_hub     = V(g)$is_hub,
        degree     = V(g)$degree,
        node_color = V(g)$node_color,
        node_size  = V(g)$node_size,
        stringsAsFactors = FALSE
    )
    write_csv(node_table, file.path(NETWORK_DIR, paste0(name, "_nodes.csv")))

    # Edge table
    write_csv(network_obj$edges, file.path(NETWORK_DIR, paste0(name, "_edges.csv")))

    # Hub table with descriptions
    hub_table <- network_obj$hubs
    hub_table$description <- sapply(hub_table$name, get_gene_description,
                                    aliases_df = aliases_df)
    write_csv(hub_table, file.path(NETWORK_DIR, paste0(name, "_hubs.csv")))

    # Cytoscape style XML
    style_xml <- sprintf(
'<?xml version="1.0" encoding="UTF-8"?>
<vizmap documentVersion="3.0" id="%s_style">
  <visualStyle name="%s">
    <network>
      <visualProperty name="NETWORK_BACKGROUND_PAINT" default="#FFFFFF"/>
    </network>
    <node>
      <visualProperty name="NODE_FILL_COLOR" default="#999999">
        <passthroughMapping attributeName="node_color" attributeType="string"/>
      </visualProperty>
      <visualProperty name="NODE_SIZE" default="30">
        <passthroughMapping attributeName="node_size" attributeType="float"/>
      </visualProperty>
      <visualProperty name="NODE_LABEL">
        <passthroughMapping attributeName="label" attributeType="string"/>
      </visualProperty>
      <visualProperty name="NODE_LABEL_FONT_SIZE" default="12"/>
    </node>
    <edge>
      <visualProperty name="EDGE_STROKE_UNSELECTED_PAINT" default="#CCCCCC">
        <discreteMapping attributeName="interaction_type" attributeType="string">
          <discreteMappingEntry value="%s" attributeValue="Activated"/>
          <discreteMappingEntry value="%s" attributeValue="Repressed"/>
        </discreteMapping>
      </visualProperty>
      <visualProperty name="EDGE_WIDTH" default="2"/>
    </edge>
  </visualStyle>
</vizmap>', name, name, COLOR_ACTIVATED, COLOR_REPRESSED)
    writeLines(style_xml, file.path(NETWORK_DIR, paste0(name, "_style.xml")))

    invisible(file.path(NETWORK_DIR, paste0(name, "_network.graphml")))
}

create_network_visualization <- function(network_obj, layout_type = "fr") {
    if (is.null(network_obj)) return(invisible(NULL))
    g    <- network_obj$graph
    name <- network_obj$name

    if (vcount(g) > 500) {
        rewired_nodes <- which(V(g)$is_rewired)
        neighbors_of_rewired <- unique(unlist(ego(g, order = 1,
                                                  nodes = rewired_nodes)))
        g <- induced_subgraph(g, neighbors_of_rewired)
    }

    set.seed(42)
    p <- ggraph(g, layout = layout_type) +
        geom_edge_link(aes(color = interaction_type), alpha = 0.5, width = 0.5) +
        geom_node_point(aes(size = node_size, color = I(node_color))) +
        geom_node_text(aes(label = label), repel = TRUE, size = 3,
                       fontface = "bold", max.overlaps = 20) +
        scale_edge_color_manual(
            values = c("Activated" = COLOR_ACTIVATED,
                       "Repressed" = COLOR_REPRESSED),
            name = "Interaction"
        ) +
        scale_size_identity() +
        theme_graph(base_family = "sans") +
        theme(legend.position = "right")

    ggsave(file.path(FIG_DIR, paste0("network_", name, ".png")), p,
           width = 14, height = 12, dpi = 300)
    ggsave(file.path(FIG_DIR, paste0("network_", name, ".svg")), p,
           width = 14, height = 12)
    invisible(p)
}

run_module_enrichment <- function(module_genes, module_id, fgr_ko, ko_pathway) {
    if (length(module_genes) < 5) return(NULL)
    query_ko <- fgr_ko %>% filter(gene %in% module_genes) %>%
        pull(ko) %>% unique()
    if (length(query_ko) < 3) return(NULL)

    query_pathways <- ko_pathway %>% filter(ko %in% query_ko)
    if (nrow(query_pathways) == 0) return(NULL)

    pathway_counts <- query_pathways %>%
        group_by(pathway) %>%
        summarise(observed = n_distinct(ko), .groups = "drop")

    bg_counts <- ko_pathway %>%
        filter(ko %in% unique(fgr_ko$ko)) %>%
        group_by(pathway) %>%
        summarise(expected = n_distinct(ko), .groups = "drop")

    results <- pathway_counts %>%
        inner_join(bg_counts, by = "pathway") %>%
        mutate(
            M = length(unique(fgr_ko$ko)),
            N = length(query_ko),
            # phyper(k-1, m, n, k_draw, lower.tail=FALSE) = P(X >= k)
            p_value = phyper(observed - 1, expected, M - expected, N,
                             lower.tail = FALSE)
        ) %>%
        mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
        filter(p_adj < 0.05) %>%
        arrange(p_adj)

    if (nrow(results) == 0) return(NULL)
    results$module_id <- module_id
    results
}

run_module_go_enrichment <- function(module_genes, module_id, gene_go_df,
                                      go_names_df) {
    if (length(module_genes) < 5 || nrow(gene_go_df) == 0) return(NULL)

    query_go <- gene_go_df %>% filter(gene %in% module_genes)
    if (nrow(query_go) < 3) return(NULL)

    query_counts <- query_go %>%
        group_by(go_id) %>%
        summarise(observed = n_distinct(gene), .groups = "drop")

    bg_counts <- gene_go_df %>%
        group_by(go_id) %>%
        summarise(background = n_distinct(gene), .groups = "drop")

    q_genes <- n_distinct(query_go$gene)
    b_genes <- n_distinct(gene_go_df$gene)

    results <- query_counts %>%
        inner_join(bg_counts, by = "go_id") %>%
        filter(observed >= 3) %>%
        mutate(
            total_bg        = b_genes,
            query_size      = q_genes,
            # phyper(k-1, m, n, k_draw, lower.tail=FALSE) = P(X >= k)
            p_value         = phyper(observed - 1, background,
                                     total_bg - background, query_size,
                                     lower.tail = FALSE),
            fold_enrichment = (observed / query_size) / (background / total_bg)
        ) %>%
        mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
        filter(p_adj < 0.05) %>%
        arrange(p_adj)

    if (nrow(results) == 0) return(NULL)

    results <- results %>%
        left_join(go_names_df, by = "go_id") %>%
        mutate(go_name = ifelse(is.na(go_name), go_id, go_name))

    results$module_id <- module_id
    results
}

# --- PART A -- Network Generation ---

# --- Load rewired proteins ---
rewired_all       <- read_csv(file.path(REWIRED_DIR, "all_rewired_proteins.csv"),
                              show_col_types = FALSE)
fungal_rewired    <- read_csv(file.path(REWIRED_DIR, "fungal_rewired_proteins.csv"),
                              show_col_types = FALSE)
bacterial_rewired <- read_csv(file.path(REWIRED_DIR, "bacterial_rewired_proteins.csv"),
                              show_col_types = FALSE)

# --- Load STRING aliases ---
fgr_aliases <- read_tsv(file.path(STRING_DIR, "229533.protein.aliases.v12.0.txt"),
                         col_names = c("string_id", "alias", "source"),
                         col_types = "ccc", comment = "#", show_col_types = FALSE)

# --- Load all interactions per study ---
all_interactions <- lapply(1:8, function(sid) {
    f <- file.path(DI_DIR, paste0("Study_", sid), "all_interactions.csv")
    if (!file.exists(f)) return(NULL)
    df <- read_csv(f, show_col_types = FALSE)
    df$Study    <- sid
    df$BCA_Type <- ifelse(sid %in% 1:4, "Fungal", "Bacterial")
    df
}) %>% bind_rows()

# --- Host metadata ---
host_map <- data.frame(
    Study = 1:8,
    Host  = STUDY_INFO$organism,
    stringsAsFactors = FALSE
)
rewired_all       <- rewired_all       %>% left_join(host_map, by = "Study")
fungal_rewired    <- fungal_rewired    %>% left_join(host_map, by = "Study")
bacterial_rewired <- bacterial_rewired %>% left_join(host_map, by = "Study")

# Fungal aggregated network
fungal_interactions <- all_interactions %>% filter(BCA_Type == "Fungal")
if (nrow(fungal_interactions) > 0) {
    fungal_net <- create_network(fungal_interactions,
                                 unique(fungal_rewired$Protein),
                                 "fungal_rewired")
    if (!is.null(fungal_net)) {
        export_to_cytoscape_format(fungal_net, fgr_aliases)
        create_network_visualization(fungal_net)
        save_fig(last_plot(), "network_fungal_rewired", width = 14, height = 12,
                 title = "Fungal BCA Rewired Protein Network",
                 description = paste(
                     "PPI network of rewired proteins under fungal BCAs.",
                     "Hub proteins (orange) are the top-10 rewired by degree.",
                     "Red = activated, blue = repressed interactions."))
    }
}

# Bacterial aggregated network
bacterial_interactions <- all_interactions %>% filter(BCA_Type == "Bacterial")
if (nrow(bacterial_interactions) > 0) {
    bacterial_net <- create_network(bacterial_interactions,
                                    unique(bacterial_rewired$Protein),
                                    "bacterial_rewired")
    if (!is.null(bacterial_net)) {
        export_to_cytoscape_format(bacterial_net, fgr_aliases)
        create_network_visualization(bacterial_net)
        save_fig(last_plot(), "network_bacterial_rewired", width = 14, height = 12,
                 title = "Bacterial BCA Rewired Protein Network",
                 description = paste(
                     "PPI network of rewired proteins under bacterial BCAs.",
                     "Hub proteins (orange) are the top-10 rewired by degree.",
                     "Red = activated, blue = repressed interactions."))
    }
}

# Shared rewired network
shared_file <- file.path(REWIRED_DIR, "rewired_shared.csv")
if (file.exists(shared_file)) {
    shared_proteins <- read_csv(shared_file, show_col_types = FALSE)$Protein
    if (length(shared_proteins) > 0) {
        shared_interactions <- all_interactions %>%
            filter(ProteinA %in% shared_proteins | ProteinB %in% shared_proteins)
        if (nrow(shared_interactions) > 0) {
            shared_net <- create_network(shared_interactions,
                                         shared_proteins,
                                         "shared_rewired")
            if (!is.null(shared_net)) {
                export_to_cytoscape_format(shared_net, fgr_aliases)
                create_network_visualization(shared_net)
                save_fig(last_plot(), "network_shared_rewired",
                         width = 14, height = 12,
                         title = "Shared Rewired Protein Network",
                         description = paste(
                             "PPI network of proteins rewired by both fungal",
                             "and bacterial BCAs, representing core Fusarium",
                             "biocontrol responses."))
            }
        }
    }
}

# Per-study networks
for (sid in 1:8) {
    int_file <- file.path(DI_DIR, paste0("Study_", sid), "all_interactions.csv")
    if (!file.exists(int_file)) next

    study_interactions <- read_csv(int_file, show_col_types = FALSE)
    study_rewired <- rewired_all %>%
        filter(Study == sid) %>%
        pull(Protein) %>%
        unique()
    if (length(study_rewired) == 0) next

    study_net <- create_network(study_interactions,
                                study_rewired,
                                paste0("study_", sid))
    if (!is.null(study_net)) {
        export_to_cytoscape_format(study_net, fgr_aliases)
        if (vcount(study_net$graph) <= 1000) {
            create_network_visualization(study_net)
        }
    }
}

# Network legend
legend_data <- data.frame(
    Element     = c("Hub Rewired Protein", "Rewired Protein", "Other Protein",
                    "Activated Interaction", "Repressed Interaction"),
    Color       = c(COLOR_HUB, "#4DAF4A", COLOR_REGULAR,
                    COLOR_ACTIVATED, COLOR_REPRESSED),
    Type        = c("Node", "Node", "Node", "Edge", "Edge"),
    Description = c("Top connected rewired protein (labeled)",
                    "Protein with high connectivity change but stable expression",
                    "Interaction partner of rewired protein",
                    "Interaction increased in treatment vs control",
                    "Interaction decreased in treatment vs control"),
    stringsAsFactors = FALSE
)
write_csv(legend_data, file.path(NETWORK_DIR, "network_legend.csv"))

legend_plot <- ggplot() +
    geom_point(aes(x = 1, y = 3), color = COLOR_HUB, size = 12) +
    geom_text(aes(x = 1.5, y = 3, label = "Hub Rewired Protein"),
              hjust = 0, size = 4) +
    geom_point(aes(x = 1, y = 2), color = "#4DAF4A", size = 8) +
    geom_text(aes(x = 1.5, y = 2, label = "Rewired Protein"),
              hjust = 0, size = 4) +
    geom_point(aes(x = 1, y = 1), color = COLOR_REGULAR, size = 4) +
    geom_text(aes(x = 1.5, y = 1, label = "Other Protein"),
              hjust = 0, size = 4) +
    geom_segment(aes(x = 4, xend = 4.8, y = 3, yend = 3),
                 color = COLOR_ACTIVATED, linewidth = 2) +
    geom_text(aes(x = 5, y = 3, label = "Activated Interaction"),
              hjust = 0, size = 4) +
    geom_segment(aes(x = 4, xend = 4.8, y = 2, yend = 2),
                 color = COLOR_REPRESSED, linewidth = 2) +
    geom_text(aes(x = 5, y = 2, label = "Repressed Interaction"),
              hjust = 0, size = 4) +
    xlim(0.5, 8) + ylim(0.5, 3.5) +
    theme_void()

ggsave(file.path(FIG_DIR, "network_legend.png"), legend_plot,
       width = 8, height = 4, dpi = 300)
ggsave(file.path(FIG_DIR, "network_legend.svg"), legend_plot,
       width = 8, height = 4)

# --- PART B -- Module Analysis (Louvain) ---

# --- 1. Aggregated network ---
edge_summary <- all_interactions %>%
    group_by(ProteinA, ProteinB) %>%
    summarise(
        n_studies    = n(),
        n_repressed  = sum(q_value <= 0.1),
        n_activated  = sum(q_value >= 0.9),
        mean_q       = mean(q_value),
        n_fungal     = sum(BCA_Type == "Fungal"),
        n_bacterial  = sum(BCA_Type == "Bacterial"),
        .groups      = "drop"
    ) %>%
    mutate(
        interaction_type = case_when(
            n_repressed > n_activated ~ "Repressed",
            n_activated > n_repressed ~ "Activated",
            TRUE                      ~ "Mixed"
        )
    )

g <- graph_from_data_frame(edge_summary, directed = FALSE)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

V(g)$is_rewired <- V(g)$name %in% rewired_all$Protein
V(g)$degree     <- degree(g)

# --- 2. Louvain clustering ---
# Seed is set globally in 00_setup.R (set.seed(42)) for reproducibility.
# Default resolution parameter (gamma = 1) is used.
set.seed(42)
louvain_result <- cluster_louvain(g)
V(g)$module    <- membership(louvain_result)
n_modules      <- max(V(g)$module)

module_sizes <- sizes(louvain_result)
module_stats <- data.frame(
    module_id = as.integer(names(module_sizes)),
    size      = as.integer(module_sizes)
) %>% arrange(desc(size))

significant_modules <- module_stats %>% filter(size >= MIN_MODULE_SIZE)

# --- 3. Module-trait analysis ---
fg_proteins <- rewired_all$Protein[rewired_all$Host == "F. graminearum"]
fo_proteins <- rewired_all$Protein[rewired_all$Host == "F. oxysporum"]

module_traits <- data.frame(module_id = significant_modules$module_id)

module_traits$prop_rewired <- sapply(module_traits$module_id, function(m) {
    members <- V(g)$name[V(g)$module == m]
    sum(members %in% rewired_all$Protein) / length(members)
})

module_traits$prop_fungal <- sapply(module_traits$module_id, function(m) {
    members <- V(g)$name[V(g)$module == m]
    nf <- sum(members %in% fungal_rewired$Protein)
    nb <- sum(members %in% bacterial_rewired$Protein)
    if ((nf + nb) == 0) return(0.5)
    nf / (nf + nb)
})

module_traits$prop_fg <- sapply(module_traits$module_id, function(m) {
    members <- V(g)$name[V(g)$module == m]
    nfg <- sum(members %in% fg_proteins)
    nfo <- sum(members %in% fo_proteins)
    if ((nfg + nfo) == 0) return(0.5)
    nfg / (nfg + nfo)
})

module_traits$mean_degree <- sapply(module_traits$module_id, function(m) {
    mean(degree(g, v = which(V(g)$module == m)))
})

module_traits$module_size <- significant_modules$size

# --- 4. Module classification ---
module_traits <- module_traits %>%
    mutate(
        classification = case_when(
            prop_fg > 0.8 ~ "F. graminearum-specific",
            prop_fg < 0.2 ~ "F. oxysporum-specific",
            TRUE          ~ "Shared/Mixed"
        )
    )

# module_traits saved after enrichment adds functional labels (section 7b)

# --- 5. Module-trait heatmap ---
# (Created after enrichment + functional labeling in section 7c below)

# --- 6. Hub gene identification per module ---
module_hubs_list <- list()

n_to_process <- min(TOP_MODULES, nrow(significant_modules))
for (m_id in significant_modules$module_id[seq_len(n_to_process)]) {
    member_idx <- which(V(g)$module == m_id)
    subg <- induced_subgraph(g, member_idx)
    intra_deg <- degree(subg)
    names(intra_deg) <- V(subg)$name

    hub_df <- data.frame(
        protein      = names(intra_deg),
        intra_degree = as.numeric(intra_deg),
        is_rewired   = names(intra_deg) %in% rewired_all$Protein,
        stringsAsFactors = FALSE
    ) %>%
        arrange(desc(intra_degree)) %>%
        head(TOP_HUBS_PER_MODULE)

    hub_df$module_id <- m_id
    module_hubs_list[[as.character(m_id)]] <- hub_df
}

all_hubs <- bind_rows(module_hubs_list)
write_csv(all_hubs, file.path(MODULE_DIR, "module_hub_genes.csv"))

# --- 7. Module KEGG + GO enrichment ---
fgr_ko <- tryCatch({
    read_tsv(file.path(STRING_DIR, "fgr_ko.txt"),
             col_names = c("gene", "ko"),
             col_types = "cc", show_col_types = FALSE) %>%
        mutate(gene = gsub("^fgr:", "", gene),
               ko   = gsub("^ko:", "", ko))
}, error = function(e) {
    data.frame(gene = character(), ko = character(),
               stringsAsFactors = FALSE)
})

ko_pathway <- tryCatch({
    lines <- readLines("https://rest.kegg.jp/link/pathway/ko")
    data.frame(
        ko      = gsub("ko:", "", gsub("\t.*", "", lines)),
        pathway = gsub(".*\tpath:", "", lines),
        stringsAsFactors = FALSE
    ) %>% filter(grepl("^map", pathway))
}, error = function(e) {
    data.frame(ko = character(), pathway = character(),
               stringsAsFactors = FALSE)
})

# Pathway names for labeling
pathway_names <- tryCatch({
    raw <- readLines("https://rest.kegg.jp/list/pathway/fgr", warn = FALSE)
    data.frame(
        pathway_id   = gsub("path:fgr", "map", gsub("\t.*", "", raw)),
        pathway_name = gsub(".*\t", "", raw) %>%
            gsub(" - Fusarium graminearum.*", "", .),
        stringsAsFactors = FALSE
    )
}, error = function(e) {
    data.frame(pathway_id = character(), pathway_name = character())
})

# Load GO annotations from UniProt file
go_file <- file.path(STRING_DIR, "fgr_uniprot_go.tsv")
gene_go_mod <- data.frame(gene = character(), go_id = character(),
                           stringsAsFactors = FALSE)
go_names_mod <- data.frame(go_id = character(), go_name = character(),
                            stringsAsFactors = FALSE)

if (file.exists(go_file)) {
    go_raw <- read_tsv(go_file, show_col_types = FALSE)
    names(go_raw) <- c("entry", "gene_names", "go_ids", "go_terms")

    # Build FGRAMPH1 -> FGSG mapping from STRING aliases
    fgsg_al  <- fgr_aliases %>%
        filter(grepl("^FGSG_\\d+$", alias)) %>%
        select(string_id, fgsg = alias) %>%
        distinct(string_id, .keep_all = TRUE)
    fgram_al <- fgr_aliases %>%
        filter(grepl("^FGRAMPH1_", alias)) %>%
        select(string_id, fgramph1 = alias) %>%
        distinct(string_id, .keep_all = TRUE)
    fgram_to_fgsg <- inner_join(fgram_al, fgsg_al, by = "string_id") %>%
        select(fgramph1, fgsg)

    go_gene_list <- list()
    go_name_list <- list()

    for (i in seq_len(nrow(go_raw))) {
        gn <- go_raw$gene_names[i]
        go_ids_str <- go_raw$go_ids[i]
        go_terms_str <- go_raw$go_terms[i]
        if (is.na(gn) || is.na(go_ids_str)) next

        fgsg_ids <- str_extract_all(gn, "FGSG_\\d+")[[1]]
        if (length(fgsg_ids) == 0) {
            fgram_ids <- str_extract_all(gn, "FGRAMPH1_\\d+T\\d+")[[1]]
            if (length(fgram_ids) > 0) {
                mapped <- fgram_to_fgsg %>% filter(fgramph1 %in% fgram_ids)
                fgsg_ids <- mapped$fgsg
            }
        }
        if (length(fgsg_ids) == 0) next

        go_ids <- trimws(unlist(strsplit(go_ids_str, ";")))
        go_ids <- go_ids[grepl("^GO:", go_ids)]
        if (length(go_ids) == 0) next

        for (gg in unique(fgsg_ids)) {
            go_gene_list[[length(go_gene_list) + 1]] <- data.frame(
                gene = gg, go_id = go_ids, stringsAsFactors = FALSE)
        }

        go_desc_parts <- trimws(unlist(strsplit(go_terms_str, ";")))
        for (part in go_desc_parts) {
            m <- regmatches(part, regexec("^(.+?)\\s*\\[(GO:\\d+)\\]$", part))[[1]]
            if (length(m) == 3) {
                go_name_list[[length(go_name_list) + 1]] <- data.frame(
                    go_id = m[3], go_name = trimws(m[2]),
                    stringsAsFactors = FALSE)
            }
        }
    }

    gene_go_mod <- bind_rows(go_gene_list) %>% distinct()
    go_names_mod <- bind_rows(go_name_list) %>% distinct(go_id, .keep_all = TRUE)
}

# --- Run KEGG + GO enrichment per module ---
module_enrichment_kegg <- list()
module_enrichment_go   <- list()

for (m_id in significant_modules$module_id[seq_len(n_to_process)]) {
    members <- V(g)$name[V(g)$module == m_id]

    # KEGG
    kegg_res <- run_module_enrichment(members, m_id, fgr_ko, ko_pathway)
    if (!is.null(kegg_res)) {
        kegg_res <- kegg_res %>%
            left_join(pathway_names, by = c("pathway" = "pathway_id")) %>%
            mutate(pathway_name = ifelse(is.na(pathway_name), pathway,
                                          pathway_name))
        module_enrichment_kegg[[as.character(m_id)]] <- kegg_res
    }

    # GO
    go_res <- run_module_go_enrichment(members, m_id, gene_go_mod, go_names_mod)
    if (!is.null(go_res)) {
        module_enrichment_go[[as.character(m_id)]] <- go_res
    }
}

# Combine and save KEGG
all_enrichment <- data.frame()
if (length(module_enrichment_kegg) > 0) {
    all_enrichment <- bind_rows(module_enrichment_kegg)
    write_csv(all_enrichment, file.path(MODULE_DIR, "module_enrichment_results.csv"))
}

# Combine and save GO
all_go_enrichment <- data.frame()
if (length(module_enrichment_go) > 0) {
    all_go_enrichment <- bind_rows(module_enrichment_go)
    write_csv(all_go_enrichment,
              file.path(MODULE_DIR, "module_go_enrichment_results.csv"))
}

# Excel workbook with both KEGG and GO sheets
wb <- createWorkbook()
for (m_id in unique(c(names(module_enrichment_kegg),
                       names(module_enrichment_go)))) {
    sheet_name <- paste0("Module_", m_id)
    addWorksheet(wb, sheet_name)
    kegg_part <- module_enrichment_kegg[[m_id]]
    go_part   <- module_enrichment_go[[m_id]]
    combined <- bind_rows(
        if (!is.null(kegg_part))
            kegg_part %>% mutate(source = "KEGG", term = pathway_name) %>%
            select(source, term, observed, p_adj),
        if (!is.null(go_part))
            go_part %>% mutate(source = "GO", term = go_name) %>%
            select(source, term, observed, p_adj)
    )
    if (nrow(combined) > 0) writeData(wb, sheet_name, combined)
}
saveWorkbook(wb, file.path(MODULE_DIR, "module_enrichment.xlsx"),
             overwrite = TRUE)

# --- 7b. Functional labeling of modules ---

# For each module, pick the best functional label from KEGG + GO results.
# Priority: most specific enriched term (highest fold enrichment among top hits)

assign_module_label <- function(m_id) {
    terms <- c()

    # Top KEGG pathway
    if (as.character(m_id) %in% names(module_enrichment_kegg)) {
        kegg <- module_enrichment_kegg[[as.character(m_id)]] %>%
            arrange(p_adj) %>% head(3)
        # Skip generic pathways like "Metabolic pathways", "Biosynthesis of
        # secondary metabolites" if more specific ones exist
        generic <- c("Metabolic pathways", "Biosynthesis of secondary metabolites",
                      "Biosynthesis of amino acids", "Carbon metabolism",
                      "Microbial metabolism in diverse environments")
        specific <- kegg %>% filter(!pathway_name %in% generic)
        if (nrow(specific) > 0) {
            terms <- c(terms, specific$pathway_name[1])
        } else if (nrow(kegg) > 0) {
            terms <- c(terms, kegg$pathway_name[1])
        }
    }

    # Top GO biological process (prefer BP terms for labeling)
    if (as.character(m_id) %in% names(module_enrichment_go)) {
        go <- module_enrichment_go[[as.character(m_id)]] %>%
            arrange(p_adj) %>% head(5)
        # Try to find a BP term (not CC or MF-like terms)
        # Simple heuristic: skip terms that are clearly CC/MF
        cc_mf_patterns <- c("^ribosome$", "^cytoplasm$", "^nucleus$",
                             "^membrane$", "binding$", "activity$",
                             "^mitochondri", "^cytosol")
        bp_terms <- go %>%
            filter(!grepl(paste(cc_mf_patterns, collapse = "|"),
                          go_name, ignore.case = TRUE))
        if (nrow(bp_terms) > 0) {
            terms <- c(terms, bp_terms$go_name[1])
        } else if (nrow(go) > 0) {
            terms <- c(terms, go$go_name[1])
        }
    }

    if (length(terms) == 0) return("Uncharacterized")

    # Use the most informative term (shortest unique descriptor)
    # If KEGG and GO give the same info, use KEGG. Otherwise combine.
    label <- terms[1]
    if (length(terms) > 1 && !grepl(tolower(terms[2]),
                                      tolower(terms[1]), fixed = TRUE)) {
        label <- paste0(terms[1], " / ", terms[2])
    }
    # Truncate if too long
    if (nchar(label) > 80) label <- paste0(substr(label, 1, 77), "...")
    label
}

module_traits$functional_label <- sapply(module_traits$module_id,
                                          assign_module_label)

# Add top KEGG terms (comma-separated, up to 3)
module_traits$top_kegg <- sapply(module_traits$module_id, function(m_id) {
    if (!as.character(m_id) %in% names(module_enrichment_kegg)) return("")
    kegg <- module_enrichment_kegg[[as.character(m_id)]] %>%
        arrange(p_adj) %>% head(3)
    paste(kegg$pathway_name, collapse = "; ")
})

# Add top GO terms (comma-separated, up to 3)
module_traits$top_go <- sapply(module_traits$module_id, function(m_id) {
    if (!as.character(m_id) %in% names(module_enrichment_go)) return("")
    go <- module_enrichment_go[[as.character(m_id)]] %>%
        arrange(p_adj) %>% head(3)
    paste(go$go_name, collapse = "; ")
})

# Re-save module_traits with functional labels
write_csv(module_traits, file.path(MODULE_DIR, "module_traits.csv"))

# Save a detailed module characterization table
module_char <- module_traits %>%
    select(module_id, module_size, prop_rewired, prop_fungal, prop_fg,
           mean_degree, classification, functional_label, top_kegg, top_go) %>%
    arrange(desc(module_size))
write_csv(module_char,
          file.path(MODULE_DIR, "module_functional_characterization.csv"))

# Top terms summary
top_terms_summary <- bind_rows(
    if (nrow(all_enrichment) > 0)
        all_enrichment %>%
        group_by(module_id) %>%
        slice_min(p_adj, n = 3) %>%
        mutate(source = "KEGG",
               term = ifelse("pathway_name" %in% names(.), pathway_name,
                              pathway)) %>%
        select(module_id, source, term, observed, p_adj) %>%
        ungroup(),
    if (nrow(all_go_enrichment) > 0)
        all_go_enrichment %>%
        group_by(module_id) %>%
        slice_min(p_adj, n = 3) %>%
        mutate(source = "GO", term = go_name) %>%
        select(module_id, source, term, observed, p_adj) %>%
        ungroup()
)
if (nrow(top_terms_summary) > 0) {
    write_csv(top_terms_summary,
              file.path(MODULE_DIR, "module_top_terms_summary.csv"))
}

# --- 7c. Module-trait heatmap (with functional labels) ---
n_show <- min(15, nrow(module_traits))
heatmap_data <- module_traits %>%
    select(module_id, prop_rewired, prop_fungal, prop_fg, mean_degree) %>%
    head(n_show)

heatmap_mat <- as.matrix(heatmap_data[, -1])
# Row labels include module number + functional label
row_labels <- sapply(seq_len(n_show), function(i) {
    m_id <- module_traits$module_id[i]
    label <- module_traits$functional_label[i]
    if (nchar(label) > 45) label <- paste0(substr(label, 1, 42), "...")
    sprintf("M%d: %s", m_id, label)
})
rownames(heatmap_mat) <- row_labels
colnames(heatmap_mat) <- c("Rewired\nProp", "Fungal\nBias",
                            "F. graminearum\nBias", "Mean\nDegree")

heatmap_mat_scaled <- scale(heatmap_mat)

row_anno <- data.frame(
    Classification = as.character(module_traits$classification[seq_len(n_show)]),
    row.names = rownames(heatmap_mat)
)

anno_colors <- list(Classification = c(
    "F. graminearum-specific" = "#388E3C",
    "F. oxysporum-specific"   = "#D32F2F",
    "Shared/Mixed"            = "#FFA000"
))

tryCatch({
    png(file.path(FIG_DIR, "heatmap_module_traits.png"),
        width = 14, height = 10, units = "in", res = 300)
    pheatmap(heatmap_mat_scaled,
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(50),
             cluster_rows = TRUE, cluster_cols = FALSE,
             annotation_row = row_anno, annotation_colors = anno_colors,
             display_numbers = TRUE, number_format = "%.2f",
             fontsize = 11, fontsize_number = 9, fontsize_row = 9,
             main = "Module Traits with Functional Annotations")
    dev.off()

    svglite(file.path(FIG_DIR, "heatmap_module_traits.svg"),
            width = 14, height = 10)
    pheatmap(heatmap_mat_scaled,
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(50),
             cluster_rows = TRUE, cluster_cols = FALSE,
             annotation_row = row_anno, annotation_colors = anno_colors,
             display_numbers = TRUE, number_format = "%.2f",
             fontsize = 11, fontsize_number = 9, fontsize_row = 9,
             main = "Module Traits with Functional Annotations")
    dev.off()
}, error = function(e) {
    if (dev.cur() > 1) dev.off()
})

# --- 8. Module network plots (top 5 modules) ---
n_plot <- min(5, nrow(significant_modules))
for (m_id in significant_modules$module_id[seq_len(n_plot)]) {
    members_idx <- which(V(g)$module == m_id)
    if (length(members_idx) < 3) next

    subg <- induced_subgraph(g, members_idx)
    V(subg)$node_color <- ifelse(V(subg)$is_rewired, "#4DAF4A", "#999999")
    V(subg)$node_size  <- ifelse(V(subg)$is_rewired, 8, 4)
    V(subg)$label <- ifelse(
        V(subg)$name %in% all_hubs$protein[all_hubs$module_id == m_id],
        V(subg)$name, ""
    )

    # Get functional label for subtitle
    func_label <- ""
    trait_row <- module_traits %>% filter(module_id == m_id)
    if (nrow(trait_row) > 0 && trait_row$functional_label[1] != "Uncharacterized") {
        func_label <- paste0("\n", trait_row$functional_label[1])
    }

    set.seed(42)
    p_mod <- ggraph(subg, layout = "fr") +
        geom_edge_link(alpha = 0.3, width = 0.5) +
        geom_node_point(aes(size = node_size, color = I(node_color))) +
        geom_node_text(aes(label = label), repel = TRUE, size = 3,
                       fontface = "bold") +
        scale_size_identity() +
        theme_graph(base_family = "sans") +
        labs(subtitle = sprintf("Module %d (%d genes)%s",
                                m_id, vcount(subg), func_label))

    ggsave(file.path(FIG_DIR, sprintf("network_module_%d.png", m_id)),
           p_mod, width = 10, height = 10, dpi = 300)
    ggsave(file.path(FIG_DIR, sprintf("network_module_%d.svg", m_id)),
           p_mod, width = 10, height = 10)
}

# --- 9. Export module GraphML files ---
for (m_id in significant_modules$module_id[seq_len(n_to_process)]) {
    members_idx <- which(V(g)$module == m_id)
    if (length(members_idx) < MIN_MODULE_SIZE) next
    subg <- induced_subgraph(g, members_idx)
    write_graph(subg,
                file.path(NETWORK_DIR, sprintf("module_%d.graphml", m_id)),
                format = "graphml")
}

# --- 10. BCA comparison bar plot ---
rewired_modules <- data.frame(
    protein  = rewired_all$Protein,
    study    = rewired_all$Study,
    bca_type = rewired_all$BCA_Type,
    stringsAsFactors = FALSE
)
rewired_modules$module <- V(g)$module[match(rewired_modules$protein, V(g)$name)]
rewired_modules <- rewired_modules %>% filter(!is.na(module))

module_bca_counts <- rewired_modules %>%
    group_by(module, bca_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = bca_type, values_from = n, values_fill = 0)
write_csv(module_bca_counts, file.path(MODULE_DIR, "module_bca_distribution.csv"))

# Ensure both columns exist even if all rewired belong to one type
if (!"Fungal"    %in% names(module_bca_counts)) module_bca_counts$Fungal    <- 0L
if (!"Bacterial" %in% names(module_bca_counts)) module_bca_counts$Bacterial <- 0L

module_bca_long <- module_bca_counts %>%
    filter(module %in% significant_modules$module_id[seq_len(n_to_process)]) %>%
    pivot_longer(cols = c(Fungal, Bacterial),
                 names_to = "BCA_Type", values_to = "Count")

p_compare <- ggplot(module_bca_long,
                    aes(x = factor(module), y = Count, fill = BCA_Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = COLORS_BCA) +
    labs(x = "Module", y = "Number of Rewired Proteins", fill = "BCA Type") +
    theme_bw(base_size = 14) +
    theme(legend.position = "right")

ggsave(file.path(FIG_DIR, "barplot_module_bca_comparison.png"), p_compare,
       width = 12, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "barplot_module_bca_comparison.svg"), p_compare,
       width = 12, height = 6)

save_fig(p_compare, "barplot_module_bca_comparison", width = 12, height = 6,
         title = "Rewired Protein Distribution Across Modules by BCA Type",
         description = paste(
             "Grouped bar plot comparing fungal vs bacterial BCA rewired",
             "protein counts per Louvain module."))

# --- 11. Module membership table ---
membership_table <- data.frame(
    protein    = V(g)$name,
    module     = as.integer(V(g)$module),
    degree     = V(g)$degree,
    is_rewired = V(g)$is_rewired,
    stringsAsFactors = FALSE
) %>%
    left_join(module_traits %>% select(module_id, classification, functional_label),
              by = c("module" = "module_id")) %>%
    arrange(module, desc(degree))

write_csv(membership_table, file.path(MODULE_DIR, "full_module_membership.csv"))

# --- Cytoscape Integration via RCy3 ---
tryCatch({
    if (requireNamespace("RCy3", quietly = TRUE)) {
        library(RCy3)
        tryCatch({
            cytoscapePing()

            for (m_id in significant_modules$module_id) {
                members_idx <- which(V(g)$module == m_id)
                if (length(members_idx) < MIN_MODULE_SIZE) next
                subg <- induced_subgraph(g, members_idx)

                # Add functional label as network title
                func_lbl <- module_traits %>%
                    filter(module_id == m_id) %>%
                    pull(functional_label)
                func_lbl <- if (length(func_lbl) > 0 &&
                                 func_lbl[1] != "Uncharacterized")
                    paste0(" - ", func_lbl[1]) else ""
                net_title <- paste0("Module_", m_id, func_lbl)

                tryCatch({
                    createNetworkFromIgraph(subg, title = net_title,
                                            collection = "Louvain_Modules")
                    Sys.sleep(1)  # allow Cytoscape to render

                    exportImage(file.path(MODULE_DIR,
                                          sprintf("cytoscape_module_%d.png",
                                                  m_id)),
                                type = "PNG", resolution = 300)
                    exportImage(file.path(MODULE_DIR,
                                          sprintf("cytoscape_module_%d.svg",
                                                  m_id)),
                                type = "SVG")
                }, error = function(e2) {
                })
            }
        }, error = function(e) {
            # Cytoscape not running -- skip silently
        })
    }
}, error = function(e) {
    # RCy3 not available -- skip silently
})
