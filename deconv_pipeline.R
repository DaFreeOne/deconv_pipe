#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyestimate)
  library(MCPcounter)
  library(omnideconv)
  library(org.Hs.eg.db)
  library(Seurat)
  library(jsonlite)
  library(data.table)
  library(Matrix)
})

# Used only when --ncores is not passed (i.e. ncores left NULL/blank in the YAML)
default_cores = parallel::detectCores(logical = TRUE)
if (is.na(default_cores)) default_cores = 2L
default_cores = max(1L, default_cores - 4L)

option_list = list(
  make_option("--bulk_counts", type = "character"),
  make_option("--bulk_TPM", type = "character"),
  make_option("--bulk_ID", type = "character", default = "NONAME"),
  make_option("--clinic_dir", type = "character"),
  make_option("--clinic_file", type = "character"),
  make_option("--single_cell_rds", type = "character"),
  make_option("--gse_id", type = "character"),
  make_option("--signatures_dir", type = "character"),
  make_option("--output_dir", type = "character"),
  make_option("--deconv_methods", type = "character", default = "DWLS"),
  make_option("--celltype_col", type = "character", default = "cell_type"),
  make_option("--patient_col", type = "character", default = "patient"),
  make_option("--downsample_n_cells", type = "integer", default = 1000L),
  make_option("--downsampling_method", type = "character", default = "patientxcelltype"),
  make_option("--ncores", type = "integer", default = default_cores),
  make_option("--seed", type = "integer", default = 1L),
  make_option("--dwls_method", type = "character", default = "mast_optimized")
)

opt = parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

parse_csv_arg = function(x) {
  if (is.null(x) || is.na(x) || x == "") return(character(0))
  out = trimws(unlist(strsplit(x, ",")))
  out[nzchar(out)]
}

write_csv_mkdir = function(x, file, ...) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file = file, ...)
  invisible(file)
}

save_rds_mkdir = function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file)
  invisible(file)
}

rm_if_exists = function(...) {
  objs = unlist(list(...))
  objs = objs[objs %in% ls(envir = .GlobalEnv)]
  if (length(objs) > 0) {
    rm(list = objs, envir = .GlobalEnv)
  }
  invisible(NULL)
}

message("Starting deconvolution pipeline")
message("Methods: ", opt$deconv_methods)

deconv_to_use = parse_csv_arg(opt$deconv_methods)
ncores = as.integer(opt$ncores)

stopifnot(file.exists(opt$bulk_counts) | dir.exists(opt$bulk_TPM))
stopifnot(dir.exists(opt$clinic_dir))
stopifnot(file.exists(opt$single_cell_rds))

dir.create(opt$signatures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# Save run config
jsonlite::write_json(
  list(
    bulk_counts = opt$bulk_counts,
    bulk_TPM = opt$bulk_TPM,
    clinic_dir = opt$clinic_dir,
    clinic_file = opt$clinic_file,
    single_cell_rds = opt$single_cell_rds,
    gse_id = opt$gse_id,
    signatures_dir = opt$signatures_dir,
    output_dir = opt$output_dir,
    deconv_methods = deconv_to_use,
    celltype_col = opt$celltype_col,
    patient_col = opt$patient_col,
    downsample_n_cells = opt$downsample_n_cells,
    downsampling_method = opt$downsampling_method,
    ncores = ncores,
    seed = opt$seed
  ),
  path = file.path(opt$output_dir, paste0("run_config_", opt$gse_id, ".json")),
  auto_unbox = TRUE,
  pretty = TRUE
)

# Inputs
RNAseq_TPM_path = opt$bulk_TPM
RNAseq_counts_path = opt$bulk_counts
clinic_annot_dir = opt$clinic_dir
clinic_annot_filename = opt$clinic_file
single_cell_object_path = opt$single_cell_rds
single_cell_GSE_ID = opt$gse_id
signatures_directory = opt$signatures_dir
output_dir = opt$output_dir

# Load bulk and clinic annotations
message("Reading bulk and clinic files")
if (any(c("AutogeneS", "Bseq-Sc", "CIBERSORTx", "DWLS", "MuSiC", "Rectangle", "Scaden", "SCDC") %in% unlist(deconv_to_use))) {
  RNAseq_TPM = read.csv(
    file.path(RNAseq_TPM_path),
    row.names = 1,
    check.names = FALSE
  )
}

if (any(c("CDSeq", "BayesPrism", "Bisque", "CPM", "MOMF") %in% unlist(deconv_to_use))) {
  RNAseq_counts = read.csv(
    file.path(RNAseq_counts_path),
    row.names = 1,
    check.names = FALSE
  )
}

# clinic_annot = read.csv(
#   file.path(clinic_annot_dir, clinic_annot_filename),
#   row.names = 1,
#   check.names = FALSE
# )

# # ensure patient ID column exists for joins
# if (!"ID_Patient" %in% colnames(clinic_annot)) {
#   clinic_annot$ID_Patient = rownames(clinic_annot)
# }



################# PREPROCESS DES DATA ###################
message("Processing single-cell dataset")
## Chargement des single-cell et annotations
sc_annotations = readRDS(single_cell_object_path)
sc_annotations$celltype = sc_annotations[[opt$celltype_col]]
sc_annotations$patient_id = sc_annotations[[opt$patient_col]]


## Downsampling patients/celltypes for a more representative signature (less skewed by dominant donor)
if (opt$downsampling_method == "patientxcelltype"){
  cells_keep = unlist(
    split(Cells(sc_annotations), interaction(sc_annotations$patient_id, sc_annotations$celltype, drop = TRUE)) |>
      lapply(function(cells) if (length(cells) > opt$downsample_n_cells) sample(cells, opt$downsample_n_cells) else cells)
  )
ref_ds = subset(sc_annotations, cells = cells_keep)
}

# Prepare count matrix and labels
DefaultAssay(ref_ds) = "RNA"
counts_mat = as.matrix(GetAssayData(ref_ds, assay = "RNA", layer = "counts"))
labels = ref_ds$celltype
batch_ids = ref_ds$patient_id


# Intersect on genes and keep common genes
if (any(c("AutogeneS", "Bseq-Sc", "CIBERSORTx", "DWLS", "MuSiC", "Rectangle", "Scaden", "SCDC") %in% unlist(deconv_to_use))) {
  bulk_tpm_mat = as.matrix(RNAseq_TPM)
  common_genes = intersect(rownames(counts_mat), rownames(bulk_tpm_mat))
  bulk_tpm_mat = bulk_tpm_mat[common_genes, , drop = FALSE]
}
if (any(c("CDSeq", "BayesPrism", "Bisque", "CPM", "MOMF") %in% unlist(deconv_to_use))) {
  bulk_counts_mat = as.matrix(RNAseq_counts)
  common_genes = intersect(rownames(counts_mat), rownames(bulk_counts_mat))
  bulk_counts_mat = bulk_counts_mat[common_genes, , drop = FALSE]
}

stopifnot(length(labels) == ncol(counts_mat))
stopifnot(length(batch_ids) == ncol(counts_mat))

message("Data preparation and downsampling : OK")



### DWLS ###
if ("DWLS" %in% deconv_to_use) {
  message("Running DWLS")

  sig_path = file.path(signatures_directory, "DWLS", paste0("signature_", single_cell_GSE_ID, "_downsample-", opt$downsample_n_cells))
  deconv_path = file.path(output_dir, "DWLS", paste0("deconv_DWLS_", opt$bulk_ID, "_with_", single_cell_GSE_ID, "_downsample-", opt$downsample_n_cells))

  if (!file.exists(paste0(sig_path, ".rds")) | !file.exists(paste0(sig_path, ".csv"))) {
    message(paste0("Computing ", single_cell_GSE_ID," DWLS signature"))
    signature_DWLS = omnideconv::build_model(
      single_cell_object = counts_mat,
      cell_type_annotations = labels,
      method = "DWLS",
      bulk_gene_expression = bulk_tpm_mat,
      dwls_method = opt$dwls_method,  #"seurat" or "mast_optimized"
      ncores = ncores,
      batch_ids = batch_ids
    )
    save_rds_mkdir(signature_DWLS, paste0(sig_path, ".rds"))
    write_csv_mkdir(signature_DWLS, paste0(sig_path, ".csv"))
  }else{
    message(paste0("Loading preexisting ", single_cell_GSE_ID, " DWLS signature"))
    signature_DWLS = readRDS(paste0(sig_path, ".rds"))
  }
  message(paste0("Running DWLS deconvolution using ", single_cell_GSE_ID))

  deconvolution_DWLS = omnideconv::deconvolute(
    bulk_gene_expression = bulk_tpm_mat,
    model = signature_DWLS,
    method = "DWLS",
    dwls_submethod = "DampenedWLS",
    batch_ids = batch_ids
  )

  save_rds_mkdir(deconvolution_DWLS, paste0(deconv_path, ".rds"))
  write_csv_mkdir(deconvolution_DWLS, paste0(deconv_path, ".csv"))

  rm_if_exists("signature_DWLS", "deconvolution_DWLS")
  gc()

  message(paste0("DWLS deconvolution : OK"))
}





### CDSeq ###
if ("CDSeq" %in% deconv_to_use) {
  message("Running CDSeq")
  deconv_path = file.path(output_dir, "CDSeq", paste0("deconv_CDSeq_", opt$bulk_ID, "_with_", single_cell_GSE_ID, "_downsample-", opt$downsample_n_cells))

  if (!file.exists(paste0(deconv_path, ".rds")) | !file.exists(paste0(deconv_path, ".csv"))) {
    message(paste0("Running CDSeq deconvolution using ", single_cell_GSE_ID))
    deconvolution_CDSeq = omnideconv::deconvolute(
      bulk_gene_expression = bulk_counts_mat,
      single_cell_object = counts_mat,
      cell_type_annotations = labels,
      batch_ids = batch_ids,
      method = "CDSeq"
    )
    save_rds_mkdir(deconvolution_CDSeq, paste0(deconv_path, ".rds"))
    write_csv_mkdir(deconvolution_CDSeq, paste0(deconv_path, ".csv"))

  }else{
    message(paste0("CDseq deconvolution results already exists, CDseq deconv aborted."))
    deconvolution_CDSeq = readRDS(paste0(deconv_path, ".rds"))
  }
  rm_if_exists("signature_CDSeq", "deconvolution_CDSeq")
  gc()

  message(paste0("CDSeq deconvolution : OK\n\n"))
}





### MuSiC ###
if ("MuSiC" %in% deconv_to_use) {
  message("Running MuSiC")

  # MuSiC does signature creation + deconvolution in one step (no build_model):
  # it weights genes by cross-subject consistency using the single-cell reference
  # directly. It needs single_cell_object + cell_type_annotations + batch_ids,
  # and at least two bulk samples.
  deconv_path = file.path(output_dir, "MuSiC", paste0("deconv_MuSiC_", opt$bulk_ID, "_with_", single_cell_GSE_ID, "_downsample-", opt$downsample_n_cells))

  message(paste0("Running MuSiC deconvolution using ", single_cell_GSE_ID))
  deconvolution_MuSiC = omnideconv::deconvolute(
    bulk_gene_expression  = bulk_tpm_mat,
    single_cell_object    = counts_mat,
    cell_type_annotations = labels,
    batch_ids             = batch_ids,
    method                = "MuSiC"
  )

  save_rds_mkdir(deconvolution_MuSiC, paste0(deconv_path, ".rds"))
  write_csv_mkdir(deconvolution_MuSiC, paste0(deconv_path, ".csv"))

  rm_if_exists("deconvolution_MuSiC")
  gc()

  message(paste0("MuSiC deconvolution : OK\n\n"))
}







### BayesPrism ###
# BayesPrism is one-step like CDSeq/MuSiC (no build_model): it uses the
# single-cell reference directly and takes raw counts as the bulk input.
if ("BayesPrism" %in% deconv_to_use){
  message("Running BayesPrism")
  deconv_path = file.path(output_dir, "BayesPrism", paste0("deconv_BayesPrism_", opt$bulk_ID, "_with_", single_cell_GSE_ID, "_downsample-", opt$downsample_n_cells))

  if (!file.exists(paste0(deconv_path, ".rds")) | !file.exists(paste0(deconv_path, ".csv"))) {
    message(paste0("Running BayesPrism deconvolution using ", single_cell_GSE_ID))
    deconvolution_BayesPrism = omnideconv::deconvolute(
      bulk_gene_expression  = bulk_counts_mat,
      single_cell_object    = counts_mat,
      cell_type_annotations = labels,
      batch_ids             = batch_ids,
      method                = "bayesprism",
      n_cores               = ncores
    )
    save_rds_mkdir(deconvolution_BayesPrism, paste0(deconv_path, ".rds"))
    write_csv_mkdir(deconvolution_BayesPrism, paste0(deconv_path, ".csv"))
  }else{
    message(paste0("BayesPrism deconvolution results already exist, BayesPrism deconv aborted."))
    deconvolution_BayesPrism = readRDS(paste0(deconv_path, ".rds"))
  }
  rm_if_exists("deconvolution_BayesPrism")
  gc()

  message(paste0("BayesPrism deconvolution : OK\n\n"))
}




# # Scaden
# if ("Scaden" %in% deconv_to_use){
#   sig_path = file.path(signatures_directory, "Scaden", paste0("signature_", single_cell_GSE_ID, "_downsample-", opt$downsample_n_cells))
#   deconv_path = file.path(output_dir, "Scaden", paste0("deconv_Scaden_", opt$bulk_ID, "_with_", single_cell_GSE_ID, "_downsample-", opt$downsample_n_cells))
  
#   if (!file.exists(paste0(sig_path,".rds")) | !file.exists(paste0(sig_path,".csv"))) {
#     signature_Scaden <- omnideconv::build_model(single_cell_object = counts_mat, 
#                                                 cell_type_annotations = labels,
#                                                 method = "scaden", 
#                                                 bulk_gene_expression = RNAseq_TPM,
#                                                 batch_ids = sc_annotations$patient_id)
#     saveRDS(signature_Scaden, paste0(sig_path,".rds"))
#   }else{
#     message(paste0("Loading preexisting ", single_cell_GSE_ID, " Scaden signature"))
#     signature_Scaden = readRDS(paste0(sig_path,".rds"))
#   }

#   deconvolution_Scaden <- omnideconv::deconvolute(bulk_gene_expression = RNAseq_TPM, 
#                                                   signature = signature_Scaden)
#   saveRDS(deconvolution_Scaden, paste0(deconv_path,".rds"))
# }
# rm(signature_Scaden); rm(deconvolution_Scaden)
# gc()




# ### ESTIMATE ###
# if ("ESTIMATE" %in% deconv_to_use) {
#   message("Running ESTIMATE")

#   df_vst = data.frame(gene = rownames(RNAseq_vst), RNAseq_vst, check.names = FALSE)
#   est_vst = estimate_score(df_vst, is_affymetrix = FALSE)

#   df_tpm = data.frame(gene = rownames(RNAseq_TPM), RNAseq_TPM, check.names = FALSE)
#   est_tpm = estimate_score(df_tpm, is_affymetrix = FALSE)

#   names(est_vst)[names(est_vst) == "sample"] = "ID_Patient"
#   names(est_tpm)[names(est_tpm) == "sample"] = "ID_Patient"

#   write_csv_mkdir(est_vst, file.path(output_dir, "ESTIMATE_vst.csv"), row.names = FALSE)
#   write_csv_mkdir(est_tpm, file.path(output_dir, "ESTIMATE_tpm.csv"), row.names = FALSE)

#   clinic_annot_out = dplyr::left_join(clinic_annot, est_vst, by = "ID_Patient")
#   write_csv_mkdir(
#     clinic_annot_out,
#     file.path(output_dir, "clinic_merged_cellularity.csv"),
#     row.names = FALSE
#   )

#   rm_if_exists("est_vst", "est_tpm", "df_vst", "df_tpm", "clinic_annot_out")
#   gc()
# }
