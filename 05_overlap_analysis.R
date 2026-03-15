#!/usr/bin/env Rscript
# 05_overlap_analysis.R - Venn/UpSet diagrams for rewired, DEGs, DiffInt
# Optional visualization script — not required for core analysis pipeline.
# Generates Venn diagrams, UpSet plots, and overlap statistics for
# cross-study comparison of rewired proteins, DEGs, and DiffInt results.
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
    library(VennDiagram)
    library(UpSetR)
    library(RColorBrewer)
    library(gridExtra)
    library(svglite)
    library(pheatmap)
    library(readxl)
})

# Suppress VennDiagram log file creation
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# --- Counters for final summary ---
n_venn   <- 0L
n_upset  <- 0L
n_tables <- 0L

# --- Rewired Protein Overlaps ---
# Load rewired data
rewired_all       <- read_csv(file.path(REWIRED_DIR, "all_rewired_proteins.csv"),
                              show_col_types = FALSE)
fungal_rewired    <- read_csv(file.path(REWIRED_DIR, "fungal_rewired_proteins.csv"),
                              show_col_types = FALSE)
bacterial_rewired <- read_csv(file.path(REWIRED_DIR, "bacterial_rewired_proteins.csv"),
                              show_col_types = FALSE)

# Per-study rewired protein sets
rewired_by_study <- split(rewired_all$Protein, rewired_all$Study)

# --- 1a. 4-Way Venn: Fungal rewired proteins (Studies 1-4) ---
fungal_sets <- lapply(1:4, function(i) {
    key <- as.character(i)
    if (key %in% names(rewired_by_study)) unique(rewired_by_study[[key]]) else character()
})
names(fungal_sets) <- paste0("Study_", 1:4)

non_empty_fungal <- sapply(fungal_sets, length) > 0
if (sum(non_empty_fungal) >= 2) {
    fsets <- fungal_sets[non_empty_fungal]
    venn_obj <- venn.diagram(
        x        = fsets,
        filename = NULL,
        fill     = COLORS_FUNGAL[seq_along(fsets)],
        alpha    = 0.5,
        cat.cex  = 1.2,
        cex      = 1.5,
        main     = "",
        main.cex = 1.5
    )
    save_fig(venn_obj, "venn_4way_fungal_rewired",
             width = 10, height = 10, is_grid = TRUE,
             title = "Four-Way Venn Diagram of Fungal BCA Rewired Proteins",
             description = "Overlap of rewired proteins across fungal BCA studies (Studies 1-4: T. gamsii, T. atroviride, T. hamatum, A. pullulans). Rewired proteins have significant PPI changes but stable gene expression.")
    n_venn <- n_venn + 1L
}

# --- 1b. 4-Way Venn: Bacterial rewired proteins (Studies 5-8) ---
bacterial_sets <- lapply(5:8, function(i) {
    key <- as.character(i)
    if (key %in% names(rewired_by_study)) unique(rewired_by_study[[key]]) else character()
})
names(bacterial_sets) <- paste0("Study_", 5:8)

non_empty_bacterial <- sapply(bacterial_sets, length) > 0
if (sum(non_empty_bacterial) >= 2) {
    bsets <- bacterial_sets[non_empty_bacterial]
    venn_obj <- venn.diagram(
        x        = bsets,
        filename = NULL,
        fill     = COLORS_BACTERIAL[seq_along(bsets)],
        alpha    = 0.5,
        cat.cex  = 1.2,
        cex      = 1.5,
        main     = "",
        main.cex = 1.5
    )
    save_fig(venn_obj, "venn_4way_bacterial_rewired",
             width = 10, height = 10, is_grid = TRUE,
             title = "Four-Way Venn Diagram of Bacterial BCA Rewired Proteins",
             description = "Overlap of rewired proteins across bacterial BCA studies (Studies 5-8: P. aeruginosa, B. velezensis, B. subtilis, B. amyloliquefaciens). Rewired proteins have significant PPI changes but stable gene expression.")
    n_venn <- n_venn + 1L
}

# --- 1c. UpSet plot: all 8 studies rewired proteins ---
all_study_sets <- lapply(1:8, function(i) {
    key <- as.character(i)
    if (key %in% names(rewired_by_study)) unique(rewired_by_study[[key]]) else character()
})
names(all_study_sets) <- paste0("Study_", 1:8)

non_empty_all <- sapply(all_study_sets, length) > 0
if (sum(non_empty_all) >= 2) {
    sets_filt   <- all_study_sets[non_empty_all]
    all_proteins <- unique(unlist(sets_filt))
    upset_matrix <- data.frame(Protein = all_proteins, stringsAsFactors = FALSE)
    for (sn in names(sets_filt)) {
        upset_matrix[[sn]] <- as.integer(all_proteins %in% sets_filt[[sn]])
    }

    # UpSetR does not produce ggplot objects -- save directly
    for (ext in c("png", "svg")) {
        fpath <- file.path(FIG_DIR, paste0("upset_all_studies_rewired.", ext))
        if (ext == "png") png(fpath, width = 14, height = 8, units = "in", res = 300)
        else              svglite(fpath, width = 14, height = 8)
        print(upset(upset_matrix[, -1],
                    sets           = names(sets_filt),
                    order.by       = "freq",
                    decreasing     = TRUE,
                    mb.ratio       = c(0.6, 0.4),
                    number.angles  = 30,
                    text.scale     = c(1.5, 1.3, 1.3, 1.2, 1.5, 1.2),
                    point.size     = 3.5,
                    line.size      = 1.5,
                    mainbar.y.label = "Intersection Size",
                    sets.x.label   = "Rewired Proteins per Study"))
        dev.off()
    }
    writeLines(c(
        "Title: UpSet Plot of Rewired Proteins Across All Studies", "",
        "Description: Intersection analysis of rewired proteins across all 8 BCA studies. Studies 1-4 are fungal BCAs, Studies 5-8 are bacterial BCAs."
    ), file.path(FIG_DIR, "upset_all_studies_rewired_description.txt"))

    write_csv(upset_matrix, file.path(TABLE_DIR, "upset_matrix_rewired_proteins.csv"))
    n_upset  <- n_upset  + 1L
    n_tables <- n_tables + 1L
}

# --- 1d. Cross-study heatmap of differential degree ---
behavior_matrix_file <- file.path(REWIRED_DIR, "cross_study_behavior_matrix.csv")
if (file.exists(behavior_matrix_file)) {
    behavior_matrix <- read_csv(behavior_matrix_file, show_col_types = FALSE)
    degree_cols <- grep("_Degree$", names(behavior_matrix), value = TRUE)

    if (length(degree_cols) > 1) {
        heatmap_data <- behavior_matrix %>%
            select(Protein, all_of(degree_cols)) %>%
            filter(rowSums(!is.na(across(-Protein))) >= 2)

        if (nrow(heatmap_data) > 0) {
            heatmap_mat <- as.matrix(heatmap_data[, -1])
            rownames(heatmap_mat) <- heatmap_data$Protein
            heatmap_mat[is.na(heatmap_mat)] <- 0

            col_anno <- data.frame(
                BCA_Type = ifelse(grepl("Study_[1-4]", colnames(heatmap_mat)),
                                  "Fungal", "Bacterial"),
                row.names = colnames(heatmap_mat)
            )
            colnames(heatmap_mat) <- gsub("_Degree", "", colnames(heatmap_mat))
            rownames(col_anno)    <- gsub("_Degree", "", rownames(col_anno))

            for (ext in c("png", "svg")) {
                fpath <- file.path(FIG_DIR, paste0("heatmap_rewired_cross_study.", ext))
                if (ext == "png") png(fpath, width = 12, height = 10, units = "in", res = 300)
                else              svglite(fpath, width = 12, height = 10)
                pheatmap(heatmap_mat,
                         color            = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
                         cluster_rows     = TRUE,
                         cluster_cols     = FALSE,
                         show_rownames    = nrow(heatmap_mat) <= 50,
                         annotation_col   = col_anno,
                         annotation_colors = list(BCA_Type = COLORS_BCA),
                         fontsize         = 10,
                         fontsize_row     = 8,
                         main             = "")
                dev.off()
            }
            writeLines(c(
                "Title: Cross-Study Differential Degree Heatmap of Rewired Proteins", "",
                "Description: Heatmap of differential degree for rewired proteins present in at least 2 studies. Studies annotated by BCA type (Fungal/Bacterial). Rows hierarchically clustered."
            ), file.path(FIG_DIR, "heatmap_rewired_cross_study_description.txt"))
        }
    }
}

# --- 1e. Comparison bar plot by classification ---
comparison_data <- rewired_all %>%
    group_by(Study, BCA_Type, Classification) %>%
    summarise(Count = n(), .groups = "drop")

p_comparison <- ggplot(comparison_data,
                       aes(x = factor(Study), y = Count, fill = Classification)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~BCA_Type, scales = "free_x") +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Study", y = "Number of Rewired Proteins", fill = "Classification") +
    theme_bw(base_size = 14) +
    theme(legend.position  = "right",
          strip.background = element_rect(fill = "white"),
          strip.text       = element_text(face = "bold", size = 12))

save_fig(p_comparison, "barplot_rewired_by_study_classification",
         width = 12, height = 6,
         title = "Classification of Rewired Proteins by Study",
         description = "Stacked bar plot of rewired proteins per study classified by dominant interaction type (Mostly Repressed, Mostly Activated, or Balanced), grouped by BCA type.")

# --- 1f. Summary statistics and statistical tests ---
summary_stats <- rewired_all %>%
    group_by(BCA_Type) %>%
    summarise(
        Total_Rewired        = n(),
        Unique_Proteins      = n_distinct(Protein),
        Mean_Degree          = round(mean(Ntotal, na.rm = TRUE), 2),
        Mean_Repressed       = round(mean(Nrep,   na.rm = TRUE), 2),
        Mean_Activated       = round(mean(Nact,   na.rm = TRUE), 2),
        Mostly_Repressed_Pct = round(100 * sum(Classification == "Mostly_Repressed", na.rm = TRUE) / n(), 1),
        Mostly_Activated_Pct = round(100 * sum(Classification == "Mostly_Activated", na.rm = TRUE) / n(), 1),
        .groups = "drop"
    )
write_csv(summary_stats, file.path(TABLE_DIR, "rewired_summary_by_bca_type.csv"))
n_tables <- n_tables + 1L

if (nrow(fungal_rewired) > 0 && nrow(bacterial_rewired) > 0) {
    wilcox_degree <- wilcox.test(fungal_rewired$Ntotal, bacterial_rewired$Ntotal)

    class_table <- table(rewired_all$BCA_Type, rewired_all$Classification)
    chisq_class <- if (all(dim(class_table) > 1)) chisq.test(class_table)
                   else list(p.value = NA)

    stat_tests <- data.frame(
        Test        = c("Differential Degree (Wilcoxon)", "Classification Distribution (Chi-sq)"),
        P_Value     = c(wilcox_degree$p.value, chisq_class$p.value),
        Significant = c(wilcox_degree$p.value < 0.05, chisq_class$p.value < 0.05)
    )
    write_csv(stat_tests, file.path(TABLE_DIR, "statistical_comparisons_bca_type.csv"))
    n_tables <- n_tables + 1L
}

# --- DEG Overlaps ---
# Load DEG data from DE_DIR/Study_X/deseq2_full_results.csv
deg_data <- list()
for (i in 1:8) {
    deg_file <- file.path(DE_DIR, paste0("Study_", i), "deseq2_full_results.csv")
    if (file.exists(deg_file)) {
        deg_data[[i]] <- read_csv(deg_file, show_col_types = FALSE) %>%
            mutate(Study = i)
    }
}

if (length(deg_data) > 0) {
    all_degs_df <- bind_rows(deg_data)

    # Extract DEG sets per study
    get_degs <- function(study_id, direction = "all") {
        d <- all_degs_df %>%
            filter(Study == study_id, !is.na(log2FoldChange), !is.na(padj))
        if (direction == "up")   d <- d %>% filter(log2FoldChange >  1, padj < 0.05)
        if (direction == "down") d <- d %>% filter(log2FoldChange < -1, padj < 0.05)
        if (direction == "all")  d <- d %>% filter(abs(log2FoldChange) > 1, padj < 0.05)
        unique(d$gene)
    }

    deg_sets      <- lapply(1:8, get_degs, direction = "all")
    deg_up_sets   <- lapply(1:8, get_degs, direction = "up")
    deg_down_sets <- lapply(1:8, get_degs, direction = "down")

    # --- 2a. 4-Way Venn: Fungal all DEGs ---
    venn_obj <- venn.diagram(
        x = setNames(deg_sets[1:4], paste0("Study_", 1:4)),
        filename = NULL,
        fill     = COLORS_FUNGAL,
        alpha    = 0.5,
        cat.cex  = 1.2,
        cex      = 1.5,
        main     = ""
    )
    save_fig(venn_obj, "venn_4way_fungal_degs",
             width = 10, height = 10, is_grid = TRUE,
             title = "DEG Overlap: Fungal BCAs (All DEGs)",
             description = "4-way Venn diagram of all DEGs across fungal BCA studies (Studies 1-4).")
    n_venn <- n_venn + 1L

    # --- 2b. 4-Way Venn: Bacterial all DEGs ---
    venn_obj <- venn.diagram(
        x = setNames(deg_sets[5:8], paste0("Study_", 5:8)),
        filename = NULL,
        fill     = COLORS_BACTERIAL,
        alpha    = 0.5,
        cat.cex  = 1.2,
        cex      = 1.5,
        main     = ""
    )
    save_fig(venn_obj, "venn_4way_bacterial_degs",
             width = 10, height = 10, is_grid = TRUE,
             title = "DEG Overlap: Bacterial BCAs (All DEGs)",
             description = "4-way Venn diagram of all DEGs across bacterial BCA studies (Studies 5-8).")
    n_venn <- n_venn + 1L

    # --- 2c. 4-Way Venn: Fungal up-regulated ---
    venn_obj <- venn.diagram(
        x = setNames(deg_up_sets[1:4], paste0("Study_", 1:4)),
        filename = NULL,
        fill     = COLORS_FUNGAL,
        alpha    = 0.5,
        cat.cex  = 1.2,
        cex      = 1.5,
        main     = ""
    )
    save_fig(venn_obj, "venn_4way_fungal_degs_up",
             width = 10, height = 10, is_grid = TRUE,
             title = "DEG Overlap: Fungal BCAs (Upregulated)",
             description = "4-way Venn diagram of upregulated DEGs across fungal BCA studies.")
    n_venn <- n_venn + 1L

    # --- 2d. 4-Way Venn: Fungal down-regulated ---
    venn_obj <- venn.diagram(
        x = setNames(deg_down_sets[1:4], paste0("Study_", 1:4)),
        filename = NULL,
        fill     = COLORS_FUNGAL,
        alpha    = 0.5,
        cat.cex  = 1.2,
        cex      = 1.5,
        main     = ""
    )
    save_fig(venn_obj, "venn_4way_fungal_degs_down",
             width = 10, height = 10, is_grid = TRUE,
             title = "DEG Overlap: Fungal BCAs (Downregulated)",
             description = "4-way Venn diagram of downregulated DEGs across fungal BCA studies.")
    n_venn <- n_venn + 1L

    # --- 2e. 4-Way Venn: Bacterial up-regulated ---
    venn_obj <- venn.diagram(
        x = setNames(deg_up_sets[5:8], paste0("Study_", 5:8)),
        filename = NULL,
        fill     = COLORS_BACTERIAL,
        alpha    = 0.5,
        cat.cex  = 1.2,
        cex      = 1.5,
        main     = ""
    )
    save_fig(venn_obj, "venn_4way_bacterial_degs_up",
             width = 10, height = 10, is_grid = TRUE,
             title = "DEG Overlap: Bacterial BCAs (Upregulated)",
             description = "4-way Venn diagram of upregulated DEGs across bacterial BCA studies.")
    n_venn <- n_venn + 1L

    # --- 2f. 4-Way Venn: Bacterial down-regulated ---
    venn_obj <- venn.diagram(
        x = setNames(deg_down_sets[5:8], paste0("Study_", 5:8)),
        filename = NULL,
        fill     = COLORS_BACTERIAL,
        alpha    = 0.5,
        cat.cex  = 1.2,
        cex      = 1.5,
        main     = ""
    )
    save_fig(venn_obj, "venn_4way_bacterial_degs_down",
             width = 10, height = 10, is_grid = TRUE,
             title = "DEG Overlap: Bacterial BCAs (Downregulated)",
             description = "4-way Venn diagram of downregulated DEGs across bacterial BCA studies.")
    n_venn <- n_venn + 1L

    # --- 2g. UpSet plot: all studies DEGs ---
    deg_named <- setNames(deg_sets, paste0("Study_", 1:8))
    non_empty_deg <- sapply(deg_named, length) > 0
    if (sum(non_empty_deg) >= 2) {
        dsets     <- deg_named[non_empty_deg]
        all_genes <- unique(unlist(dsets))
        deg_matrix <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)
        for (sn in names(dsets)) {
            deg_matrix[[sn]] <- as.integer(all_genes %in% dsets[[sn]])
        }

        for (ext in c("png", "svg")) {
            fpath <- file.path(FIG_DIR, paste0("upset_all_studies_degs.", ext))
            if (ext == "png") png(fpath, width = 14, height = 8, units = "in", res = 300)
            else              svglite(fpath, width = 14, height = 8)
            print(upset(deg_matrix[, -1],
                        sets           = names(dsets),
                        order.by       = "freq",
                        decreasing     = TRUE,
                        mb.ratio       = c(0.6, 0.4),
                        number.angles  = 30,
                        text.scale     = c(1.5, 1.3, 1.3, 1.2, 1.5, 1.2),
                        point.size     = 3.5,
                        line.size      = 1.5,
                        mainbar.y.label = "Intersection Size",
                        sets.x.label   = "DEGs per Study"))
            dev.off()
        }
        writeLines(c(
            "Title: UpSet Plot of DEGs Across All Studies", "",
            "Description: Intersection analysis of differentially expressed genes across all 8 BCA studies."
        ), file.path(FIG_DIR, "upset_all_studies_degs_description.txt"))

        write_csv(deg_matrix, file.path(TABLE_DIR, "upset_matrix_degs.csv"))
        n_upset  <- n_upset  + 1L
        n_tables <- n_tables + 1L
    }

    # --- 2h. Pairwise overlap statistics (Jaccard) ---
    compute_pairwise <- function(sets, study_ids) {
        results <- data.frame(
            Comparison    = character(),
            Study_1_Count = integer(),
            Study_2_Count = integer(),
            Shared        = integer(),
            Jaccard       = numeric(),
            stringsAsFactors = FALSE
        )
        for (a in 1:(length(study_ids) - 1)) {
            for (b in (a + 1):length(study_ids)) {
                s1 <- study_ids[a]; s2 <- study_ids[b]
                shared     <- length(intersect(sets[[s1]], sets[[s2]]))
                union_size <- length(union(sets[[s1]], sets[[s2]]))
                results <- rbind(results, data.frame(
                    Comparison    = paste(s1, "vs", s2),
                    Study_1_Count = length(sets[[s1]]),
                    Study_2_Count = length(sets[[s2]]),
                    Shared        = shared,
                    Jaccard       = if (union_size > 0) shared / union_size else 0
                ))
            }
        }
        results
    }

    fungal_overlaps    <- compute_pairwise(deg_sets, 1:4)
    bacterial_overlaps <- compute_pairwise(deg_sets, 5:8)

    write_csv(fungal_overlaps,    file.path(TABLE_DIR, "venn_deg_fungal_overlaps.csv"))
    write_csv(bacterial_overlaps, file.path(TABLE_DIR, "venn_deg_bacterial_overlaps.csv"))
    n_tables <- n_tables + 2L

    # --- 2i. Shared / specific gene lists ---
    fungal_all_degs    <- Reduce(union, deg_sets[1:4])
    bacterial_all_degs <- Reduce(union, deg_sets[5:8])
    deg_shared              <- intersect(fungal_all_degs, bacterial_all_degs)
    deg_fungal_specific     <- setdiff(fungal_all_degs,    bacterial_all_degs)
    deg_bacterial_specific  <- setdiff(bacterial_all_degs, fungal_all_degs)

    write_csv(data.frame(Gene = deg_shared),             file.path(TABLE_DIR, "deg_shared_proteins.csv"))
    write_csv(data.frame(Gene = deg_fungal_specific),    file.path(TABLE_DIR, "deg_fungal_specific.csv"))
    write_csv(data.frame(Gene = deg_bacterial_specific), file.path(TABLE_DIR, "deg_bacterial_specific.csv"))
    n_tables <- n_tables + 3L
}

# --- DiffInt Overlaps ---
# Load protein_stats.csv from DI_DIR/Study_X/
activated_proteins  <- list()
repressed_proteins  <- list()

for (i in 1:8) {
    pstats_file <- file.path(DI_DIR, paste0("Study_", i), "protein_stats.csv")
    if (file.exists(pstats_file)) {
        stats <- read_csv(pstats_file, show_col_types = FALSE)
        activated_proteins[[i]]  <- stats %>% filter(Nact > Nrep) %>% pull(Protein)
        repressed_proteins[[i]]  <- stats %>% filter(Nrep > Nact) %>% pull(Protein)
    } else {
        # Fallback: try to derive from all_interactions.csv
        int_file <- file.path(DI_DIR, paste0("Study_", i), "all_interactions.csv")
        if (file.exists(int_file)) {
            ints <- read_csv(int_file, show_col_types = FALSE)
            act_int <- ints %>% filter(!is.na(q_value), q_value >= 0.9)
            rep_int <- ints %>% filter(!is.na(q_value), q_value <= 0.1)
            activated_proteins[[i]]  <- unique(c(act_int$ProteinA, act_int$ProteinB))
            repressed_proteins[[i]]  <- unique(c(rep_int$ProteinA, rep_int$ProteinB))
        } else {
            activated_proteins[[i]]  <- character()
            repressed_proteins[[i]]  <- character()
        }
    }
}

# --- 3a. 4-Way Venn: Fungal activated proteins ---
if (any(sapply(activated_proteins[1:4], length) > 0)) {
    venn_obj <- venn.diagram(
        x = setNames(activated_proteins[1:4], paste0("Study_", 1:4)),
        filename = NULL,
        fill     = COLORS_FUNGAL,
        alpha    = 0.5,
        cat.cex  = 1.2,
        cex      = 1.5,
        main     = ""
    )
    save_fig(venn_obj, "venn_4way_fungal_activated_proteins",
             width = 10, height = 10, is_grid = TRUE,
             title = "Activated Interaction Proteins: Fungal BCAs",
             description = "4-way Venn diagram of proteins with predominantly activated interactions across fungal BCA studies (Studies 1-4).")
    n_venn <- n_venn + 1L
}

# --- 3b. 4-Way Venn: Bacterial activated proteins ---
if (any(sapply(activated_proteins[5:8], length) > 0)) {
    venn_obj <- venn.diagram(
        x = setNames(activated_proteins[5:8], paste0("Study_", 5:8)),
        filename = NULL,
        fill     = COLORS_BACTERIAL,
        alpha    = 0.5,
        cat.cex  = 1.2,
        cex      = 1.5,
        main     = ""
    )
    save_fig(venn_obj, "venn_4way_bacterial_activated_proteins",
             width = 10, height = 10, is_grid = TRUE,
             title = "Activated Interaction Proteins: Bacterial BCAs",
             description = "4-way Venn diagram of proteins with predominantly activated interactions across bacterial BCA studies (Studies 5-8).")
    n_venn <- n_venn + 1L
}

# --- 3c. 4-Way Venn: Fungal repressed proteins ---
if (any(sapply(repressed_proteins[1:4], length) > 0)) {
    venn_obj <- venn.diagram(
        x = setNames(repressed_proteins[1:4], paste0("Study_", 1:4)),
        filename = NULL,
        fill     = COLORS_FUNGAL,
        alpha    = 0.5,
        cat.cex  = 1.2,
        cex      = 1.5,
        main     = ""
    )
    save_fig(venn_obj, "venn_4way_fungal_repressed_proteins",
             width = 10, height = 10, is_grid = TRUE,
             title = "Repressed Interaction Proteins: Fungal BCAs",
             description = "4-way Venn diagram of proteins with predominantly repressed interactions across fungal BCA studies (Studies 1-4).")
    n_venn <- n_venn + 1L
}

# --- 3d. 4-Way Venn: Bacterial repressed proteins ---
if (any(sapply(repressed_proteins[5:8], length) > 0)) {
    venn_obj <- venn.diagram(
        x = setNames(repressed_proteins[5:8], paste0("Study_", 5:8)),
        filename = NULL,
        fill     = COLORS_BACTERIAL,
        alpha    = 0.5,
        cat.cex  = 1.2,
        cex      = 1.5,
        main     = ""
    )
    save_fig(venn_obj, "venn_4way_bacterial_repressed_proteins",
             width = 10, height = 10, is_grid = TRUE,
             title = "Repressed Interaction Proteins: Bacterial BCAs",
             description = "4-way Venn diagram of proteins with predominantly repressed interactions across bacterial BCA studies (Studies 5-8).")
    n_venn <- n_venn + 1L
}

# --- 3e. Overlap statistics for DiffInt ---
activated_fungal_all    <- Reduce(union, activated_proteins[1:4])
activated_bacterial_all <- Reduce(union, activated_proteins[5:8])
repressed_fungal_all    <- Reduce(union, repressed_proteins[1:4])
repressed_bacterial_all <- Reduce(union, repressed_proteins[5:8])

activated_shared             <- intersect(activated_fungal_all, activated_bacterial_all)
activated_fungal_specific    <- setdiff(activated_fungal_all,    activated_bacterial_all)
activated_bacterial_specific <- setdiff(activated_bacterial_all, activated_fungal_all)

repressed_shared             <- intersect(repressed_fungal_all, repressed_bacterial_all)
repressed_fungal_specific    <- setdiff(repressed_fungal_all,    repressed_bacterial_all)
repressed_bacterial_specific <- setdiff(repressed_bacterial_all, repressed_fungal_all)

write_csv(data.frame(Protein = activated_shared),             file.path(TABLE_DIR, "diffint_activated_shared.csv"))
write_csv(data.frame(Protein = activated_fungal_specific),    file.path(TABLE_DIR, "diffint_activated_fungal_specific.csv"))
write_csv(data.frame(Protein = activated_bacterial_specific), file.path(TABLE_DIR, "diffint_activated_bacterial_specific.csv"))
write_csv(data.frame(Protein = repressed_shared),             file.path(TABLE_DIR, "diffint_repressed_shared.csv"))
write_csv(data.frame(Protein = repressed_fungal_specific),    file.path(TABLE_DIR, "diffint_repressed_fungal_specific.csv"))
write_csv(data.frame(Protein = repressed_bacterial_specific), file.path(TABLE_DIR, "diffint_repressed_bacterial_specific.csv"))

diffint_summary <- data.frame(
    Category = c("Activated_Fungal", "Activated_Bacterial", "Activated_Shared",
                 "Repressed_Fungal", "Repressed_Bacterial", "Repressed_Shared"),
    Count    = c(length(activated_fungal_all), length(activated_bacterial_all),
                 length(activated_shared),
                 length(repressed_fungal_all), length(repressed_bacterial_all),
                 length(repressed_shared)),
    stringsAsFactors = FALSE
)
write_csv(diffint_summary, file.path(TABLE_DIR, "diffint_statistics_summary.csv"))
n_tables <- n_tables + 7L
