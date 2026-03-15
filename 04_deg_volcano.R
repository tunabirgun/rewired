#!/usr/bin/env Rscript
# 04_deg_volcano.R - Volcano plots from existing DEG results
# Optional visualization script — not required for core analysis pipeline.
# Generates per-study and combined 8-panel volcano plots with optional
# rewired protein overlay.
local({
    args <- commandArgs(trailingOnly = FALSE)
    f <- grep("--file=", args, value = TRUE)
    d <- if (length(f)) dirname(normalizePath(sub("--file=", "", f))) else getwd()
    source(file.path(d, "00_setup.R"))
})

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(readr)
})

# Load orthology mappings
fox_to_fgr <- NULL
fpu_to_fgr <- NULL

fox_ortho_file <- file.path(RESOURCE_DIR, "fox_to_fgr_orthology.csv")
fpu_ortho_file <- file.path(RESOURCE_DIR, "fpu_to_fgr_orthology.csv")

if (file.exists(fox_ortho_file)) {
    fox_to_fgr <- read.csv(fox_ortho_file, stringsAsFactors = FALSE)
}

if (file.exists(fpu_ortho_file)) {
    fpu_to_fgr <- read.csv(fpu_ortho_file, stringsAsFactors = FALSE)
}

# Map gene IDs to FGSG (used for rewired protein matching)
map_to_fgsg <- function(gene_ids, gene_prefix) {
    if (gene_prefix == "FGSG") {
        return(gene_ids)
    } else if (gene_prefix == "FOXG" && !is.null(fox_to_fgr)) {
        mapping <- setNames(fox_to_fgr$FGSG_ID, fox_to_fgr$FOXG_ID)
        return(mapping[gene_ids])
    } else if (gene_prefix == "FPSE" && !is.null(fpu_to_fgr)) {
        mapping <- setNames(fpu_to_fgr$FGSG_ID, fpu_to_fgr$FPSE_ID)
        return(mapping[gene_ids])
    }
    return(rep(NA, length(gene_ids)))
}

# Prepare volcano data.frame from DESeq2 results
prepare_volcano_df <- function(res_df) {
    res_df %>%
        filter(!is.na(padj) & !is.infinite(-log10(padj))) %>%
        mutate(
            minusLog10Padj = -log10(padj),
            sig = case_when(
                padj < PADJ_THRESHOLD & log2FoldChange >  LOG2FC_THRESHOLD ~ "Up",
                padj < PADJ_THRESHOLD & log2FoldChange < -LOG2FC_THRESHOLD ~ "Down",
                TRUE ~ "NS"
            )
        )
}

# Create a plain volcano plot
create_plain_volcano <- function(volc_df, study_info) {

    p <- ggplot(volc_df, aes(x = log2FoldChange, y = minusLog10Padj)) +
        geom_point(aes(color = sig), alpha = 0.6, size = 1.2) +
        scale_color_manual(values = SIG_COLORS, name = "Regulation") +
        geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD),
                   linetype = "dashed", color = "grey40", linewidth = 0.4) +
        geom_hline(yintercept = -log10(PADJ_THRESHOLD),
                   linetype = "dashed", color = "grey40", linewidth = 0.4) +
        labs(
            x = expression(Log[2]~fold~change),
            y = expression(-Log[10]~adjusted~italic(p))
        ) +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 11) +
        theme(
            legend.position = "right",
            plot.title = element_blank(),
            axis.text = element_text(color = "black"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt")
        )

    n_up   <- sum(volc_df$sig == "Up",   na.rm = TRUE)
    n_down <- sum(volc_df$sig == "Down", na.rm = TRUE)
    p <- p + annotate("text", x = Inf, y = Inf,
                       label = sprintf("Up: %d\nDown: %d", n_up, n_down),
                       hjust = 1.1, vjust = 1.5, size = 3, color = "grey30")
    p
}

# Create a rewired-highlighted volcano plot
create_rewired_volcano <- function(volc_df, rewired_proteins, study_info,
                                    n_labels = 10) {

    gene_prefix <- study_info$gene_prefix
    volc_df$fgsg_id <- map_to_fgsg(volc_df$gene, gene_prefix)

    # Mark rewired proteins
    volc_df$is_rewired <- volc_df$fgsg_id %in% rewired_proteins$Protein

    # Join connectivity / classification
    if (nrow(rewired_proteins) > 0 && "Ntotal" %in% colnames(rewired_proteins)) {
        rewired_scores <- rewired_proteins %>%
            select(Protein, Ntotal, Classification) %>%
            rename(fgsg_id = Protein)
        volc_df <- volc_df %>%
            left_join(rewired_scores, by = "fgsg_id")
    } else {
        volc_df$Ntotal <- NA
        volc_df$Classification <- NA
    }

    # Stats
    n_rewired        <- sum(volc_df$is_rewired, na.rm = TRUE)
    n_rewired_stable <- sum(volc_df$is_rewired &
                                abs(volc_df$log2FoldChange) < LOG2FC_THRESHOLD,
                            na.rm = TRUE)
    n_rewired_deg    <- sum(volc_df$is_rewired & volc_df$sig != "NS", na.rm = TRUE)

    # Top N rewired by connectivity for labels
    top_rewired <- volc_df %>%
        filter(is_rewired & !is.na(Ntotal)) %>%
        arrange(desc(Ntotal)) %>%
        head(n_labels) %>%
        mutate(label = fgsg_id)

    # Plot
    p <- ggplot(volc_df, aes(x = log2FoldChange, y = minusLog10Padj)) +
        # Background points
        geom_point(aes(color = sig), alpha = 0.6, size = 1.2) +
        scale_color_manual(values = SIG_COLORS, name = "Regulation") +
        # Threshold lines
        geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD),
                   linetype = "dashed", color = "grey40", linewidth = 0.4) +
        geom_hline(yintercept = -log10(PADJ_THRESHOLD),
                   linetype = "dashed", color = "grey40", linewidth = 0.4) +
        # Highlight rewired proteins with green circles
        geom_point(data = filter(volc_df, is_rewired),
                   aes(x = log2FoldChange, y = minusLog10Padj),
                   shape = 21, color = "#009E73", fill = NA,
                   size = 2.5, stroke = 1.2) +
        # Labels for top rewired proteins (improved repelling)
        geom_text_repel(
            data = top_rewired,
            aes(x = log2FoldChange, y = minusLog10Padj, label = label),
            color = "#009E73",
            bg.color = "white",
            bg.r = 0.15,
            size = 2.8,
            fontface = "bold",
            box.padding = 1.0,
            point.padding = 0.5,
            segment.color = "#009E73",
            segment.size = 0.3,
            segment.linetype = 1,
            segment.curvature = -0.1,
            segment.ncp = 3,
            min.segment.length = 0,
            max.overlaps = Inf,
            force = 15,
            force_pull = 0.3,
            max.iter = 10000,
            max.time = 2,
            direction = "both",
            seed = 42
        ) +
        labs(
            x = expression(Log[2]~fold~change),
            y = expression(-Log[10]~adjusted~italic(p))
        ) +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 11) +
        theme(
            legend.position = "right",
            plot.title = element_blank(),
            axis.text = element_text(color = "black"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt")
        ) +
        annotate("text", x = Inf, y = Inf,
                 label = sprintf("Rewired: %d\n(stable: %d)", n_rewired, n_rewired_stable),
                 hjust = 1.1, vjust = 1.5, size = 3, color = "#009E73")

    list(plot = p, data = volc_df)
}

# Per-study volcano plots (plain + rewired)
all_volcano_data <- list()
summary_rows     <- list()

for (study_id in 1:8) {

    study_info <- STUDY_INFO[STUDY_INFO$study_id == study_id, ]

    # Read existing DESeq2 full results
    results_file <- file.path(DE_DIR, paste0("Study_", study_id),
                              "deseq2_full_results.csv")
    if (!file.exists(results_file)) next

    res_df <- read.csv(results_file, stringsAsFactors = FALSE)

    # Prepare base volcano data
    volc_df <- prepare_volcano_df(res_df)

    # Plain volcano plot
    p_plain <- create_plain_volcano(volc_df, study_info)

    plain_fname <- paste0("volcano_plain_study_", study_id)
    save_fig(p_plain, plain_fname,
             width = 7.87, height = 6.30,
             title = sprintf("Study %d: %s vs %s - Volcano Plot",
                              study_id, study_info$organism, study_info$bca),
             description = sprintf(
                 "Volcano plot for Study %d (%s challenged by %s). Points colored by significance: Up (|log2FC| > %.1f & padj < %.2f), Down, or NS.",
                 study_id, study_info$organism, study_info$bca,
                 LOG2FC_THRESHOLD, PADJ_THRESHOLD))

    # Rewired-highlighted volcano plot
    rewired_file <- file.path(REWIRED_DIR, paste0("Study_", study_id),
                              "rewired_proteins.csv")
    if (file.exists(rewired_file)) {
        rewired_proteins <- read.csv(rewired_file, stringsAsFactors = FALSE)

        rw_result <- create_rewired_volcano(volc_df, rewired_proteins,
                                             study_info, n_labels = 10)
        p_rewired <- rw_result$plot
        volc_df   <- rw_result$data

        # Save to FIG_DIR
        rewired_fname <- paste0("volcano_rewired_study_", study_id)
        save_fig(p_rewired, rewired_fname,
                 width = 7.87, height = 6.30,
                 title = sprintf("Study %d: %s vs %s - Volcano with Rewired Proteins",
                                  study_id, study_info$organism, study_info$bca),
                 description = sprintf(
                     "Volcano plot for Study %d with rewired proteins highlighted (green circles). Top 10 rewired proteins labelled by Ntotal connectivity.",
                     study_id))

        # Also save into per-study DE_DIR
        study_de_dir <- file.path(DE_DIR, paste0("Study_", study_id))
        ggsave(file.path(study_de_dir, "volcano_rewired.png"),
               p_rewired, width = 200, height = 160, units = "mm", dpi = 300)
        ggsave(file.path(study_de_dir, "volcano_rewired.svg"),
               p_rewired, width = 200, height = 160, units = "mm")

    } else {
        # Ensure columns exist even without rewired data
        volc_df$fgsg_id       <- map_to_fgsg(volc_df$gene, study_info$gene_prefix)
        volc_df$is_rewired    <- FALSE
        volc_df$Ntotal        <- NA
        volc_df$Classification <- NA
    }

    # Accumulate data
    volc_df$study_id <- study_id
    all_volcano_data[[as.character(study_id)]] <- volc_df

    # Summary row
    summary_rows[[as.character(study_id)]] <- data.frame(
        Study          = study_id,
        Organism       = study_info$organism,
        BCA            = study_info$bca,
        Total_Genes    = nrow(volc_df),
        Up_DEGs        = sum(volc_df$sig == "Up",   na.rm = TRUE),
        Down_DEGs      = sum(volc_df$sig == "Down", na.rm = TRUE),
        Rewired_Total  = sum(volc_df$is_rewired,    na.rm = TRUE),
        Rewired_Stable = sum(volc_df$is_rewired &
                                 abs(volc_df$log2FoldChange) < LOG2FC_THRESHOLD,
                             na.rm = TRUE),
        stringsAsFactors = FALSE
    )
}

# Combined 8-panel plain volcano
if (length(all_volcano_data) > 0) {

    combined_df <- bind_rows(all_volcano_data)

    combined_df <- combined_df %>%
        left_join(STUDY_INFO %>% select(study_id, bca, bca_type),
                  by = "study_id") %>%
        mutate(study_label = paste0("Study ", study_id, ": ", bca))

    # --- Plain panel ---
    p_plain_panel <- ggplot(combined_df,
                            aes(x = log2FoldChange, y = minusLog10Padj)) +
        geom_point(aes(color = sig), alpha = 0.5, size = 0.8) +
        scale_color_manual(values = SIG_COLORS, name = "Regulation") +
        geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD),
                   linetype = "dashed", color = "grey40", linewidth = 0.3) +
        geom_hline(yintercept = -log10(PADJ_THRESHOLD),
                   linetype = "dashed", color = "grey40", linewidth = 0.3) +
        facet_wrap(~ study_label, ncol = 4, scales = "free") +
        labs(
            x = expression(Log[2]~fold~change),
            y = expression(-Log[10]~adjusted~italic(p))
        ) +
        theme_classic(base_size = 9) +
        theme(
            legend.position = "bottom",
            strip.background = element_rect(fill = "grey95", color = NA),
            strip.text = element_text(size = 8, face = "bold"),
            axis.text = element_text(color = "black", size = 7),
            panel.border = element_rect(color = "grey50", fill = NA,
                                        linewidth = 0.3)
        )

    save_fig(p_plain_panel, "volcano_panel_plain_all_studies",
             width = 13.39, height = 7.87,
             title = "Combined plain volcano plots for all 8 studies",
             description = "8-panel volcano plots showing DEG significance for each study without rewired protein overlay.")

    # Combined 8-panel rewired volcano
    # Select top 5 rewired per study for labeling in panel
    top_rewired_panel <- combined_df %>%
        filter(is_rewired & !is.na(Ntotal)) %>%
        group_by(study_id) %>%
        arrange(desc(Ntotal)) %>%
        slice_head(n = 5) %>%
        ungroup() %>%
        mutate(label = fgsg_id)

    p_rewired_panel <- ggplot(combined_df,
                              aes(x = log2FoldChange, y = minusLog10Padj)) +
        geom_point(aes(color = sig), alpha = 0.5, size = 0.8) +
        scale_color_manual(values = SIG_COLORS, name = "Regulation") +
        geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD),
                   linetype = "dashed", color = "grey40", linewidth = 0.3) +
        geom_hline(yintercept = -log10(PADJ_THRESHOLD),
                   linetype = "dashed", color = "grey40", linewidth = 0.3) +
        # Highlight rewired proteins
        geom_point(data = filter(combined_df, is_rewired),
                   shape = 21, color = "#009E73", fill = NA,
                   size = 1.5, stroke = 0.8) +
        # Labels for top 5 rewired per study (improved repelling)
        geom_text_repel(
            data = top_rewired_panel,
            aes(x = log2FoldChange, y = minusLog10Padj, label = label),
            color = "#009E73",
            bg.color = "white",
            bg.r = 0.12,
            size = 1.6,
            fontface = "bold",
            box.padding = 0.6,
            point.padding = 0.3,
            segment.color = "#009E73",
            segment.size = 0.2,
            min.segment.length = 0,
            max.overlaps = Inf,
            force = 12,
            force_pull = 0.2,
            max.iter = 8000,
            direction = "both",
            seed = 42
        ) +
        facet_wrap(~ study_label, ncol = 4, scales = "free") +
        labs(
            x = expression(Log[2]~fold~change),
            y = expression(-Log[10]~adjusted~italic(p))
        ) +
        theme_classic(base_size = 9) +
        theme(
            legend.position = "bottom",
            strip.background = element_rect(fill = "grey95", color = NA),
            strip.text = element_text(size = 8, face = "bold"),
            axis.text = element_text(color = "black", size = 7),
            panel.border = element_rect(color = "grey50", fill = NA,
                                        linewidth = 0.3)
        )

    save_fig(p_rewired_panel, "volcano_panel_all_studies",
             width = 13.39, height = 7.87,
             title = "Combined volcano plots with rewired protein overlay",
             description = "8-panel volcano plots showing DEGs with rewired proteins highlighted in green. Top 5 rewired proteins per study labelled by Ntotal connectivity.")

    # Summary outputs
    # Combined volcano data
    export_cols <- intersect(
        c("gene", "log2FoldChange", "padj", "sig", "fgsg_id",
          "is_rewired", "Ntotal", "Classification", "study_id"),
        colnames(combined_df)
    )
    write.csv(combined_df[, export_cols],
              file.path(DE_DIR, "all_studies_volcano_data.csv"),
              row.names = FALSE)

    # Summary statistics
    summary_stats <- bind_rows(summary_rows)
    write.csv(summary_stats,
              file.path(DE_DIR, "deg_summary_statistics.csv"),
              row.names = FALSE)

    # Print summary
    print(summary_stats)

}
