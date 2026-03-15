#!/usr/bin/env Rscript
# 06_enrichment.R - KEGG + GO enrichment for rewired proteins and DEGs
#
# Enrichment test: hypergeometric (Fisher's exact) via phyper().
# The call phyper(observed - 1, ..., lower.tail = FALSE) computes
# P(X >= observed), because phyper(q, ..., lower.tail = FALSE) = P(X > q),
# so P(X >= k) = P(X > k-1) = phyper(k-1, ..., lower.tail = FALSE).
#
# KEGG pathway names are fetched live from rest.kegg.jp for F. graminearum.
# If the KEGG REST API is unavailable, pathway IDs are used as fallback names.
local({
    args <- commandArgs(trailingOnly = FALSE)
    f <- grep("--file=", args, value = TRUE)
    d <- if (length(f)) dirname(normalizePath(sub("--file=", "", f))) else getwd()
    source(file.path(d, "00_setup.R"))
})

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(tidyr)
    library(svglite)
    library(openxlsx)
    library(httr)
    library(stringr)
})

# --- Load KEGG and GO annotation data ---
# KEGG KO assignments for F. graminearum
fgr_ko <- read_tsv(file.path(STRING_DIR, "fgr_ko.txt"),
                    col_names = c("gene", "ko"),
                    col_types = "cc", show_col_types = FALSE)
fgr_ko$gene <- gsub("^fgr:", "", fgr_ko$gene)
fgr_ko$ko   <- gsub("^ko:", "", fgr_ko$ko)

# KEGG pathway names (fetched live from rest.kegg.jp)
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

# KO to pathway mapping
ko_pathway <- tryCatch({
    raw <- readLines("https://rest.kegg.jp/link/pathway/ko", warn = FALSE)
    data.frame(
        ko      = gsub("ko:", "", gsub("\t.*", "", raw)),
        pathway = gsub(".*\tpath:", "", raw),
        stringsAsFactors = FALSE
    ) %>% filter(grepl("^map", pathway))
}, error = function(e) {
    data.frame(ko = character(), pathway = character())
})

# --- GO annotations from UniProt GO file + STRING aliases for FGSG mapping ---
go_file <- file.path(STRING_DIR, "fgr_uniprot_go.tsv")
if (file.exists(go_file)) {
    go_raw <- read_tsv(go_file, show_col_types = FALSE)
    names(go_raw) <- c("entry", "gene_names", "go_ids", "go_terms")

    # Build FGRAMPH1 -> FGSG mapping from STRING aliases
    alias_file <- file.path(STRING_DIR, "229533.protein.aliases.v12.0.txt")
    if (file.exists(alias_file)) {
        al <- read_tsv(alias_file, show_col_types = FALSE)
        names(al) <- c("string_id", "alias", "source")
        fgsg_al  <- al %>% filter(grepl("^FGSG_\\d+$", alias)) %>%
            select(string_id, fgsg = alias) %>% distinct(string_id, .keep_all = TRUE)
        fgram_al <- al %>% filter(grepl("^FGRAMPH1_", alias)) %>%
            select(string_id, fgramph1 = alias) %>% distinct(string_id, .keep_all = TRUE)
        fgram_to_fgsg <- inner_join(fgram_al, fgsg_al, by = "string_id") %>%
            select(fgramph1, fgsg)
    } else {
        fgram_to_fgsg <- data.frame(fgramph1 = character(), fgsg = character())
    }

    # Parse each row: extract FGSG IDs from gene_names, or map FGRAMPH1 -> FGSG
    gene_go_list <- list()
    go_names_list <- list()

    for (i in seq_len(nrow(go_raw))) {
        gn <- go_raw$gene_names[i]
        go_ids_str <- go_raw$go_ids[i]
        go_terms_str <- go_raw$go_terms[i]
        if (is.na(gn) || is.na(go_ids_str)) next

        # Extract FGSG IDs directly
        fgsg_ids <- str_extract_all(gn, "FGSG_\\d+")[[1]]

        # If no FGSG, try mapping via FGRAMPH1
        if (length(fgsg_ids) == 0) {
            fgram_ids <- str_extract_all(gn, "FGRAMPH1_\\d+T\\d+")[[1]]
            if (length(fgram_ids) > 0) {
                mapped <- fgram_to_fgsg %>% filter(fgramph1 %in% fgram_ids)
                fgsg_ids <- mapped$fgsg
            }
        }
        if (length(fgsg_ids) == 0) next

        # Parse GO IDs
        go_ids <- trimws(unlist(strsplit(go_ids_str, ";")))
        go_ids <- go_ids[grepl("^GO:", go_ids)]
        if (length(go_ids) == 0) next

        # Build gene-GO pairs
        for (g in unique(fgsg_ids)) {
            gene_go_list[[length(gene_go_list) + 1]] <- data.frame(
                gene = g, go_id = go_ids, stringsAsFactors = FALSE
            )
        }

        # Parse GO term names from the descriptions column
        go_desc_parts <- trimws(unlist(strsplit(go_terms_str, ";")))
        for (part in go_desc_parts) {
            m <- regmatches(part, regexec("^(.+?)\\s*\\[(GO:\\d+)\\]$", part))[[1]]
            if (length(m) == 3) {
                go_names_list[[length(go_names_list) + 1]] <- data.frame(
                    go_id = m[3], go_name = trimws(m[2]), stringsAsFactors = FALSE
                )
            }
        }
    }

    gene_go <- bind_rows(gene_go_list) %>% distinct()
    go_term_names <- bind_rows(go_names_list) %>% distinct(go_id, .keep_all = TRUE)
} else {
    gene_go <- data.frame(gene = character(), go_id = character())
    go_term_names <- data.frame(go_id = character(), go_name = character())
}

# --- Enrichment functions ---

# --- KEGG pathway enrichment (hypergeometric test, FDR < 0.05) ---
kegg_enrichment <- function(gene_list, background_genes = NULL, name = "query") {
    if (length(gene_list) < 3) return(NULL)

    query_ko <- fgr_ko %>%
        filter(gene %in% gene_list) %>%
        pull(ko) %>% unique()
    if (length(query_ko) < 3) return(NULL)

    if (is.null(background_genes)) {
        background_ko <- unique(fgr_ko$ko)
    } else {
        background_ko <- fgr_ko %>%
            filter(gene %in% background_genes) %>%
            pull(ko) %>% unique()
    }

    query_pathways <- ko_pathway %>% filter(ko %in% query_ko)
    if (nrow(query_pathways) == 0) return(NULL)

    pathway_counts <- query_pathways %>%
        group_by(pathway) %>%
        summarise(observed = n_distinct(ko), .groups = "drop")

    bg_counts <- ko_pathway %>%
        filter(ko %in% background_ko) %>%
        group_by(pathway) %>%
        summarise(background = n_distinct(ko), .groups = "drop")

    results <- pathway_counts %>%
        inner_join(bg_counts, by = "pathway") %>%
        mutate(
            total_bg        = length(background_ko),
            query_size      = length(query_ko),
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
        left_join(pathway_names, by = c("pathway" = "pathway_id")) %>%
        mutate(pathway_name = ifelse(is.na(pathway_name), pathway, pathway_name)) %>%
        select(pathway_id = pathway, pathway_name, observed, background,
               total_background = total_bg, query_size, fold_enrichment,
               p_value, p_adj)

    return(results)
}

# --- GO term enrichment (hypergeometric test, FDR < 0.05) ---
go_enrichment <- function(gene_list, background_genes = NULL, name = "query") {
    if (length(gene_list) < 3 || nrow(gene_go) == 0) return(NULL)

    query_go_df <- gene_go %>% filter(gene %in% gene_list)
    if (nrow(query_go_df) < 3) return(NULL)

    if (is.null(background_genes)) {
        background_go <- gene_go
    } else {
        background_go <- gene_go %>% filter(gene %in% background_genes)
    }

    query_counts <- query_go_df %>%
        group_by(go_id) %>%
        summarise(observed = n_distinct(gene), .groups = "drop")

    bg_counts <- background_go %>%
        group_by(go_id) %>%
        summarise(background = n_distinct(gene), .groups = "drop")

    query_genes_with_go <- n_distinct(query_go_df$gene)
    bg_genes_with_go    <- n_distinct(background_go$gene)

    results <- query_counts %>%
        inner_join(bg_counts, by = "go_id") %>%
        filter(observed >= 3) %>%
        mutate(
            total_bg        = bg_genes_with_go,
            query_size      = query_genes_with_go,
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
        left_join(go_term_names, by = "go_id") %>%
        mutate(go_name = ifelse(is.na(go_name), go_id, go_name)) %>%
        select(go_id, go_name, observed, background,
               total_background = total_bg, query_size,
               fold_enrichment, p_value, p_adj)

    return(results)
}

# --- Bar plot for enrichment results ---
create_enrichment_plot <- function(enrich_df, filename,
                                   term_col = "pathway_name", top_n = 20) {
    if (is.null(enrich_df) || nrow(enrich_df) == 0) return(invisible(NULL))

    plot_data <- enrich_df %>%
        arrange(p_adj) %>%
        head(top_n) %>%
        mutate(
            neglog10p    = -log10(p_adj),
            term         = .data[[term_col]],
            term_display = ifelse(nchar(term) > 50,
                                  paste0(substr(term, 1, 47), "..."), term),
            term_display = factor(term_display, levels = rev(term_display))
        )

    p <- ggplot(plot_data,
                aes(x = neglog10p, y = term_display, fill = fold_enrichment)) +
        geom_bar(stat = "identity") +
        scale_fill_gradient(low = "#FFC107", high = "#D32F2F",
                            name = "Fold\nEnrichment") +
        labs(x = "-log10(FDR-adjusted p-value)", y = "", title = "") +
        theme_bw(base_size = 11) +
        theme(axis.text.y   = element_text(size = 9),
              legend.position = "right",
              plot.margin    = margin(10, 10, 10, 10))

    ggsave(file.path(FIG_DIR, paste0(filename, ".png")), p,
           width = 12, height = 8, dpi = 300)
    ggsave(file.path(FIG_DIR, paste0(filename, ".svg")), p,
           width = 12, height = 8)

    return(invisible(p))
}

# --- Comparison dot plot (fungal vs bacterial) ---
create_comparison_dotplot <- function(results_list, ids_fungal, ids_bacterial,
                                      id_col, name_col, filename,
                                      colour_low = "#FFC107",
                                      colour_high = "#D32F2F") {

    fungal_df <- results_list[[ids_fungal]]
    bacterial_df <- results_list[[ids_bacterial]]
    if (is.null(fungal_df) || is.null(bacterial_df)) return(invisible(NULL))

    top_f <- fungal_df %>% arrange(p_adj) %>% head(10)
    top_b <- bacterial_df %>% arrange(p_adj) %>% head(10)

    all_ids <- unique(c(top_f[[id_col]], top_b[[id_col]]))

    comparison <- bind_rows(
        fungal_df    %>% filter(.data[[id_col]] %in% all_ids) %>%
            mutate(BCA_Type = "Fungal"),
        bacterial_df %>% filter(.data[[id_col]] %in% all_ids) %>%
            mutate(BCA_Type = "Bacterial")
    )
    if (nrow(comparison) == 0) return(invisible(NULL))

    comparison <- comparison %>%
        mutate(term_short = ifelse(nchar(.data[[name_col]]) > 40,
                                   paste0(substr(.data[[name_col]], 1, 37), "..."),
                                   .data[[name_col]]))

    p <- ggplot(comparison,
                aes(x = BCA_Type, y = reorder(term_short, -p_adj))) +
        geom_point(aes(size = -log10(p_adj), color = fold_enrichment)) +
        scale_color_gradient(low = colour_low, high = colour_high,
                             name = "Fold\nEnrichment") +
        scale_size_continuous(range = c(3, 10), name = "-log10\n(FDR p-value)") +
        labs(x = "BCA Type", y = "") +
        theme_bw(base_size = 12) +
        theme(axis.text.y = element_text(size = 9))

    ggsave(file.path(FIG_DIR, paste0(filename, ".png")), p,
           width = 14, height = 10, dpi = 300)
    ggsave(file.path(FIG_DIR, paste0(filename, ".svg")), p,
           width = 14, height = 10)

    return(invisible(p))
}

# --- Rewired protein enrichment ---
# Load rewired protein sets
rewired_all <- read_csv(file.path(REWIRED_DIR, "all_rewired_proteins.csv"),
                        show_col_types = FALSE)
fungal_rewired    <- read_csv(file.path(REWIRED_DIR, "fungal_rewired_proteins.csv"),
                              show_col_types = FALSE)
bacterial_rewired <- read_csv(file.path(REWIRED_DIR, "bacterial_rewired_proteins.csv"),
                              show_col_types = FALSE)

shared_file   <- file.path(REWIRED_DIR, "rewired_shared.csv")
fungal_sp_file <- file.path(REWIRED_DIR, "rewired_fungal_specific.csv")
bact_sp_file   <- file.path(REWIRED_DIR, "rewired_bacterial_specific.csv")

shared_proteins    <- if (file.exists(shared_file))
    read_csv(shared_file, show_col_types = FALSE)$Protein else character()
fungal_specific    <- if (file.exists(fungal_sp_file))
    read_csv(fungal_sp_file, show_col_types = FALSE)$Protein else character()
bacterial_specific <- if (file.exists(bact_sp_file))
    read_csv(bact_sp_file, show_col_types = FALSE)$Protein else character()

# Analysis groups
rewired_groups <- list(
    Fungal_All         = unique(fungal_rewired$Protein),
    Bacterial_All      = unique(bacterial_rewired$Protein),
    Fungal_Specific    = fungal_specific,
    Bacterial_Specific = bacterial_specific,
    Shared             = shared_proteins
)

all_kegg_rewired <- list()
all_go_rewired   <- list()

for (gname in names(rewired_groups)) {
    genes <- rewired_groups[[gname]]
    if (length(genes) < 3) next

    # KEGG
    kegg_res <- kegg_enrichment(genes, name = gname)
    if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
        all_kegg_rewired[[gname]] <- kegg_res
        write_csv(kegg_res,
                  file.path(ENRICH_DIR,
                            paste0("kegg_enrichment_", tolower(gname), ".csv")))
        create_enrichment_plot(kegg_res,
                               paste0("enrichment_kegg_", tolower(gname)),
                               term_col = "pathway_name")
    }

    # GO
    go_res <- go_enrichment(genes, name = gname)
    if (!is.null(go_res) && nrow(go_res) > 0) {
        all_go_rewired[[gname]] <- go_res
        write_csv(go_res,
                  file.path(ENRICH_DIR,
                            paste0("go_enrichment_", tolower(gname), ".csv")))
        create_enrichment_plot(go_res,
                               paste0("enrichment_go_", tolower(gname)),
                               term_col = "go_name")
    }
}

# Per-study enrichment (studies with >= 10 genes)
per_study_kegg_rew <- list()
per_study_go_rew   <- list()

for (sid in seq_len(nrow(STUDY_INFO))) {
    study_genes <- rewired_all %>%
        filter(Study == sid) %>%
        pull(Protein) %>% unique()

    if (length(study_genes) >= 10) {
        kegg_res <- kegg_enrichment(study_genes, name = paste0("Study_", sid))
        if (!is.null(kegg_res)) per_study_kegg_rew[[paste0("Study_", sid)]] <- kegg_res

        go_res <- go_enrichment(study_genes, name = paste0("Study_", sid))
        if (!is.null(go_res)) per_study_go_rew[[paste0("Study_", sid)]] <- go_res
    }
}

# Comparison dot plots
create_comparison_dotplot(all_kegg_rewired, "Fungal_All", "Bacterial_All",
                          "pathway_id", "pathway_name",
                          "enrichment_kegg_comparison")
create_comparison_dotplot(all_go_rewired, "Fungal_All", "Bacterial_All",
                          "go_id", "go_name",
                          "enrichment_go_comparison",
                          colour_low = "#4CAF50", colour_high = "#1976D2")

# --- DEG enrichment ---
# Load DEG data
deg_data <- list()
for (i in seq_len(nrow(STUDY_INFO))) {
    deg_file <- file.path(DE_DIR, paste0("Study_", i), "deseq2_full_results.csv")
    if (file.exists(deg_file)) {
        deg_data[[i]] <- read_csv(deg_file, show_col_types = FALSE) %>%
            mutate(Study = i)
    }
}

if (length(deg_data) > 0) {
    all_degs <- bind_rows(deg_data)

    # Build DEG gene lists
    deg_lists <- list()

    for (bca_type in c("Fungal", "Bacterial")) {
        studies <- STUDY_INFO$study_id[STUDY_INFO$bca_type == bca_type]

        # All DEGs
        deg_lists[[paste0(bca_type, "_All")]] <- all_degs %>%
            filter(Study %in% studies, !is.na(log2FoldChange), !is.na(padj),
                   abs(log2FoldChange) > 1, padj < 0.05) %>%
            pull(gene) %>% unique()

        # Upregulated
        deg_lists[[paste0(bca_type, "_Up")]] <- all_degs %>%
            filter(Study %in% studies, !is.na(log2FoldChange), !is.na(padj),
                   log2FoldChange > 1, padj < 0.05) %>%
            pull(gene) %>% unique()

        # Downregulated
        deg_lists[[paste0(bca_type, "_Down")]] <- all_degs %>%
            filter(Study %in% studies, !is.na(log2FoldChange), !is.na(padj),
                   log2FoldChange < -1, padj < 0.05) %>%
            pull(gene) %>% unique()
    }

    # Specific and shared
    deg_lists[["Fungal_Specific"]]    <- setdiff(deg_lists[["Fungal_All"]],
                                                  deg_lists[["Bacterial_All"]])
    deg_lists[["Bacterial_Specific"]] <- setdiff(deg_lists[["Bacterial_All"]],
                                                  deg_lists[["Fungal_All"]])
    deg_lists[["Shared"]]             <- intersect(deg_lists[["Fungal_All"]],
                                                    deg_lists[["Bacterial_All"]])

    # Per-study DEGs
    for (i in seq_len(nrow(STUDY_INFO))) {
        deg_lists[[paste0("Study_", i)]] <- all_degs %>%
            filter(Study == i, !is.na(log2FoldChange), !is.na(padj),
                   abs(log2FoldChange) > 1, padj < 0.05) %>%
            pull(gene) %>% unique()
    }

    background_deg <- unique(all_degs$gene)

    all_kegg_deg <- list()
    all_go_deg   <- list()

    for (cname in names(deg_lists)) {
        genes <- deg_lists[[cname]]
        if (length(genes) < 3) next

        # KEGG
        kegg_res <- kegg_enrichment(genes, background_deg, name = cname)
        if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
            all_kegg_deg[[cname]] <- kegg_res
            write_csv(kegg_res,
                      file.path(ENRICH_DIR,
                                paste0("deg_kegg_enrichment_",
                                       tolower(cname), ".csv")))
            create_enrichment_plot(kegg_res,
                                   paste0("deg_kegg_enrichment_",
                                          tolower(cname)),
                                   term_col = "pathway_name")
        }

        # GO
        go_res <- go_enrichment(genes, background_deg, name = cname)
        if (!is.null(go_res) && nrow(go_res) > 0) {
            all_go_deg[[cname]] <- go_res
            write_csv(go_res,
                      file.path(ENRICH_DIR,
                                paste0("deg_go_enrichment_",
                                       tolower(cname), ".csv")))
            create_enrichment_plot(go_res,
                                   paste0("deg_go_enrichment_",
                                          tolower(cname)),
                                   term_col = "go_name")
        }
    }
}

# --- Combined Excel workbooks ---
# --- Rewired KEGG ---
if (length(all_kegg_rewired) > 0 || length(per_study_kegg_rew) > 0) {
    wb_kegg <- createWorkbook()
    for (nm in names(all_kegg_rewired)) {
        addWorksheet(wb_kegg, nm)
        writeData(wb_kegg, nm, all_kegg_rewired[[nm]])
    }
    for (nm in names(per_study_kegg_rew)) {
        addWorksheet(wb_kegg, nm)
        writeData(wb_kegg, nm, per_study_kegg_rew[[nm]])
    }
    saveWorkbook(wb_kegg,
                 file.path(ENRICH_DIR, "kegg_enrichment_all_results.xlsx"),
                 overwrite = TRUE)
}

# --- Rewired GO ---
if (length(all_go_rewired) > 0 || length(per_study_go_rew) > 0) {
    wb_go <- createWorkbook()
    for (nm in names(all_go_rewired)) {
        addWorksheet(wb_go, nm)
        writeData(wb_go, nm, all_go_rewired[[nm]])
    }
    for (nm in names(per_study_go_rew)) {
        addWorksheet(wb_go, nm)
        writeData(wb_go, nm, per_study_go_rew[[nm]])
    }
    saveWorkbook(wb_go,
                 file.path(ENRICH_DIR, "go_enrichment_all_results.xlsx"),
                 overwrite = TRUE)
}

# --- DEG combined ---
if (exists("all_kegg_deg") &&
    (length(all_kegg_deg) > 0 || length(all_go_deg) > 0)) {
    wb_deg <- createWorkbook()
    for (nm in names(all_kegg_deg)) {
        sname <- substr(paste0("KEGG_", nm), 1, 31)
        addWorksheet(wb_deg, sname)
        writeData(wb_deg, sname, all_kegg_deg[[nm]])
    }
    for (nm in names(all_go_deg)) {
        sname <- substr(paste0("GO_", nm), 1, 31)
        addWorksheet(wb_deg, sname)
        writeData(wb_deg, sname, all_go_deg[[nm]])
    }
    saveWorkbook(wb_deg,
                 file.path(ENRICH_DIR, "deg_enrichment_all_results.xlsx"),
                 overwrite = TRUE)
}

# --- Summary tables ---
# --- Rewired enrichment summary ---
build_summary <- function(groups, kegg_results, go_results) {
    rows <- lapply(names(groups), function(nm) {
        n_kegg <- if (nm %in% names(kegg_results)) nrow(kegg_results[[nm]]) else 0L
        n_go   <- if (nm %in% names(go_results))   nrow(go_results[[nm]])   else 0L
        top_kegg <- if (n_kegg > 0) kegg_results[[nm]]$pathway_name[1] else "N/A"
        top_go   <- if (n_go   > 0) go_results[[nm]]$go_name[1]       else "N/A"
        data.frame(Analysis      = nm,
                   Gene_Count    = length(groups[[nm]]),
                   KEGG_Pathways = n_kegg,
                   GO_Terms      = n_go,
                   Top_KEGG      = substr(top_kegg, 1, 50),
                   Top_GO        = substr(top_go,   1, 50),
                   stringsAsFactors = FALSE)
    })
    bind_rows(rows)
}

rew_summary <- build_summary(rewired_groups, all_kegg_rewired, all_go_rewired)
write_csv(rew_summary, file.path(ENRICH_DIR, "enrichment_summary.csv"))

if (exists("deg_lists")) {
    deg_summary <- build_summary(deg_lists,
                                 if (exists("all_kegg_deg")) all_kegg_deg else list(),
                                 if (exists("all_go_deg"))   all_go_deg   else list())
    write_csv(deg_summary, file.path(ENRICH_DIR, "deg_enrichment_summary.csv"))
}
