#!/usr/bin/env Rscript

# ============================================================================
# Pre-analysis of the single-cell reference for the deconvolution pipeline.
#
# Profiles the cell-type composition of the reference .rds and simulates how
# different `downsample_n_cells` caps would thin it, so you can pick optimal
# downsampling parameters before running deconv_pipeline.R.
#
# Outputs (in <output_dir>/PREANALYSIS/):
#   tables : *.csv   (composition + downsampling simulation)
#   plots  : preanalysis_<gse>.pdf  (multi-page)
#
# Standalone:
#   Rscript preanalysis.R --single_cell_rds ref.rds --output_dir OUT
# Or via Docker:
#   ./run_pipeline.sh preanalysis config.yaml
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

option_list = list(
  make_option("--single_cell_rds", type = "character"),
  make_option("--celltype_col", type = "character", default = "cell_type"),
  make_option("--patient_col", type = "character", default = "patient"),
  make_option("--output_dir", type = "character", default = "."),
  make_option("--gse_id", type = "character", default = ""),
  # Candidate per-group caps to simulate (comma-separated)
  make_option("--downsample_grid", type = "character",
              default = "20,40,60,100,150,200,300,500"),
  # The cap currently set in config (highlighted on plots); <=0 to disable
  make_option("--current_cap", type = "integer", default = 40L)
)
opt = parse_args(OptionParser(option_list = option_list))

if (is.null(opt$single_cell_rds)) stop("--single_cell_rds is required")
stopifnot(file.exists(opt$single_cell_rds))

grid = as.integer(trimws(unlist(strsplit(opt$downsample_grid, ","))))
grid = sort(unique(grid[!is.na(grid) & grid > 0]))

out_dir = file.path(opt$output_dir, "PREANALYSIS")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tag = if (nzchar(opt$gse_id)) paste0("_", opt$gse_id) else ""

write_csv = function(x, name) {
  f = file.path(out_dir, paste0(name, tag, ".csv"))
  write.csv(x, f, row.names = FALSE)
  message("  wrote ", f)
}

# ---- Load metadata only ----------------------------------------------------
message("Reading reference: ", opt$single_cell_rds)
sc = readRDS(opt$single_cell_rds)
meta = as.data.frame(sc@meta.data)

for (col in c(opt$celltype_col, opt$patient_col)) {
  if (!col %in% colnames(meta)) {
    stop("Column '", col, "' not found in meta.data. Available columns:\n  ",
         paste(colnames(meta), collapse = ", "))
  }
}

celltype = as.character(meta[[opt$celltype_col]])
patient  = as.character(meta[[opt$patient_col]])

n_na = sum(is.na(celltype) | is.na(patient))
if (n_na > 0) {
  message("Dropping ", n_na, " cells with NA cell_type/patient")
  keep = !(is.na(celltype) | is.na(patient))
  celltype = celltype[keep]; patient = patient[keep]
}

n_cells = length(celltype)
message("Total cells: ", n_cells,
        " | cell types: ", length(unique(celltype)),
        " | patients: ", length(unique(patient)))

# ---- Tables ----------------------------------------------------------------
# Cells per cell type
ct = as.data.frame(table(cell_type = celltype), stringsAsFactors = FALSE)
names(ct)[2] = "n_cells"
ct$pct = round(100 * ct$n_cells / n_cells, 2)
ct = ct[order(-ct$n_cells), ]
write_csv(ct, "cells_per_celltype")

# Cells per patient
pt = as.data.frame(table(patient = patient), stringsAsFactors = FALSE)
names(pt)[2] = "n_cells"
pt$pct = round(100 * pt$n_cells / n_cells, 2)
pt = pt[order(-pt$n_cells), ]
write_csv(pt, "cells_per_patient")

# Patient x cell type (long) -> these are the "patientxcelltype" groups
pc = as.data.frame(table(patient = patient, cell_type = celltype),
                   stringsAsFactors = FALSE)
names(pc)[3] = "n_cells"
pc_nz = pc[pc$n_cells > 0, ]              # only groups that actually exist
write_csv(pc_nz, "patient_x_celltype_long")

# Patient x cell type (wide matrix)
pc_wide = as.data.frame.matrix(table(celltype, patient))
pc_wide = cbind(cell_type = rownames(pc_wide), pc_wide)
write_csv(pc_wide, "patient_x_celltype_matrix")

# Per-celltype summary of group sizes across patients (key for downsampling)
grp = split(pc_nz$n_cells, pc_nz$cell_type)
ct_summary = data.frame(
  cell_type      = names(grp),
  n_patients     = vapply(grp, length, integer(1)),
  total_cells    = vapply(grp, sum, numeric(1)),
  min_per_group  = vapply(grp, min, numeric(1)),
  median_per_grp = vapply(grp, median, numeric(1)),
  mean_per_group = round(vapply(grp, mean, numeric(1)), 1),
  max_per_group  = vapply(grp, max, numeric(1)),
  stringsAsFactors = FALSE
)
ct_summary = ct_summary[order(-ct_summary$total_cells), ]
write_csv(ct_summary, "celltype_group_summary")

# ---- Downsampling simulation ----------------------------------------------
# For each candidate cap N:
#   patientxcelltype retained = sum over (patient,celltype) groups of min(size, N)
#   celltype          retained = sum over celltype groups        of min(size, N)
g_pc = pc_nz$n_cells                       # patient x celltype group sizes
g_ct = ct$n_cells                          # celltype group sizes
sim = do.call(rbind, lapply(grid, function(N) {
  data.frame(
    cap = N,
    retained_patientxcelltype = sum(pmin(g_pc, N)),
    retained_celltype         = sum(pmin(g_ct, N)),
    full_size                 = n_cells,
    groups_capped_pxc         = sum(g_pc > N),
    groups_total_pxc          = length(g_pc),
    stringsAsFactors = FALSE
  )
}))
sim$pct_of_full_pxc = round(100 * sim$retained_patientxcelltype / n_cells, 1)
write_csv(sim, "downsampling_simulation")

# Per-celltype retained cells at each cap (patientxcelltype) -> wide
ret_by_ct = sapply(grid, function(N)
  vapply(grp, function(s) sum(pmin(s, N)), numeric(1)))
colnames(ret_by_ct) = paste0("cap_", grid)
ret_by_ct = data.frame(cell_type = names(grp), total_cells = ct_summary$total_cells[
  match(names(grp), ct_summary$cell_type)], ret_by_ct, check.names = FALSE)
write_csv(ret_by_ct, "retained_per_celltype_by_cap")

# ---- Console summary -------------------------------------------------------
message("\n--- Cell types (top by abundance) ---")
print(utils::head(ct, 20), row.names = FALSE)
message("\n--- Group-size distribution (cells per patient x cell type) ---")
print(round(quantile(g_pc, c(0, .25, .5, .75, .9, 1))))
message("\n--- Downsampling simulation ---")
print(sim, row.names = FALSE)
cur = opt$current_cap
if (cur > 0) {
  message(sprintf(
    "\nAt the current cap (%d): %d/%d patientxcelltype groups are capped (%.0f%%); reference ~%d cells (%.0f%% of full).",
    cur, sum(g_pc > cur), length(g_pc), 100 * mean(g_pc > cur),
    sum(pmin(g_pc, cur)), 100 * sum(pmin(g_pc, cur)) / n_cells))
}

# ---- Plots (multi-page PDF) ------------------------------------------------
pdf_file = file.path(out_dir, paste0("preanalysis", tag, ".pdf"))
pdf(pdf_file, width = 10, height = 7)
on.exit(dev.off(), add = TRUE)

base_theme = theme_bw(base_size = 12)
cap_lines = if (cur > 0) sort(unique(c(grid, cur))) else grid

# 1. Cells per cell type
print(
  ggplot(ct, aes(x = reorder(cell_type, n_cells), y = n_cells)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = n_cells), hjust = -0.1, size = 3) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Cells per cell type", x = NULL, y = "n cells") +
    base_theme
)

# 2. Cells per patient
print(
  ggplot(pt, aes(x = reorder(patient, n_cells), y = n_cells)) +
    geom_col(fill = "darkorange") +
    coord_flip() +
    labs(title = "Cells per patient", x = NULL, y = "n cells") +
    base_theme +
    theme(axis.text.y = element_text(size = 7))
)

# 3. Heatmap patient x cell type (log10 counts)
print(
  ggplot(pc_nz, aes(x = patient, y = cell_type, fill = log10(n_cells))) +
    geom_tile(color = "grey90") +
    scale_fill_gradient(low = "#fff5eb", high = "#7f2704",
                        name = "log10(n)") +
    labs(title = "Cells per patient x cell type",
         subtitle = "blank = combination absent", x = "patient", y = NULL) +
    base_theme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
)

# 4. Group-size distribution per cell type, with candidate caps
print(
  ggplot(pc_nz, aes(x = reorder(cell_type, n_cells, FUN = median), y = n_cells)) +
    geom_boxplot(outlier.size = 0.6, fill = "grey95") +
    geom_jitter(width = 0.15, alpha = 0.4, size = 0.7, color = "steelblue") +
    geom_hline(yintercept = cap_lines, linetype = "dashed",
               color = "firebrick", alpha = 0.6) +
    coord_flip() +
    scale_y_log10() +
    labs(title = "Group sizes (cells per patient x cell type)",
         subtitle = paste("dashed = candidate caps:",
                          paste(cap_lines, collapse = ", ")),
         x = NULL, y = "n cells per group (log10)") +
    base_theme
)

# 5. Downsampling retention curve
curve_df = rbind(
  data.frame(cap = sim$cap, retained = sim$retained_patientxcelltype,
             grouping = "patientxcelltype"),
  data.frame(cap = sim$cap, retained = sim$retained_celltype,
             grouping = "celltype")
)
p5 = ggplot(curve_df, aes(cap, retained, color = grouping)) +
  geom_hline(yintercept = n_cells, linetype = "dotted") +
  annotate("text", x = min(grid), y = n_cells, vjust = -0.5, hjust = 0,
           label = paste("full =", n_cells), size = 3) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  labs(title = "Reference size vs downsampling cap",
       x = "downsample_n_cells (per group)", y = "cells retained",
       color = "grouping") +
  base_theme
if (cur > 0) p5 = p5 + geom_vline(xintercept = cur, linetype = "dashed",
                                  color = "grey40")
print(p5)

message("\nWrote plots: ", pdf_file)
message("Pre-analysis done.")
