#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-run}"
CONFIG_YAML="${2:-config.yaml}"

if [ ! -f "$CONFIG_YAML" ]; then
  echo "Error: config file not found: $CONFIG_YAML"
  exit 1
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "Error: docker is not installed or not in PATH"
  exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
  echo "Error: python3 is required on the host to parse $CONFIG_YAML"
  exit 1
fi

if ! python3 -c "import yaml" >/dev/null 2>&1; then
  echo "Error: the PyYAML package is required on the host to parse $CONFIG_YAML"
  echo "Install it with: python3 -m pip install pyyaml"
  exit 1
fi

read_yaml() {
  local key="$1"
  python3 - "$CONFIG_YAML" "$key" <<'PY'
import sys, yaml
config_path = sys.argv[1]
key_path = sys.argv[2].split(".")
with open(config_path, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
val = cfg
for k in key_path:
    val = val[k]
if val is None:
    print("")
elif isinstance(val, bool):
    print(str(val).lower())
else:
    print(val)
PY
}

IMAGE_NAME="$(read_yaml image_name)"

HOST_BULK_COUNTS_DIR="$(read_yaml paths.bulk_counts)"
HOST_BULK_TPM_DIR="$(read_yaml paths.bulk_TPM)"
HOST_CLINIC_DIR="$(read_yaml paths.clinic_dir)"
HOST_SINGLE_CELL_RDS="$(read_yaml paths.single_cell_rds)"
HOST_SIGNATURES_DIR="$(read_yaml paths.signatures_dir)"
HOST_OUTPUT_DIR="$(read_yaml paths.output_dir)"

CLINIC_FILE="$(read_yaml inputs.clinic_file)"
GSE_ID="$(read_yaml inputs.gse_id)"
BULK_ID="$(read_yaml inputs.bulk_ID)"

DECONV_METHODS="$(read_yaml parameters.deconv_methods)"
DWLS_METHOD="$(read_yaml parameters.dwls_method)"
CELLTYPE_COL="$(read_yaml parameters.celltype_col)"
PATIENT_COL="$(read_yaml parameters.patient_col)"
DOWNSAMPLE_N_CELLS="$(read_yaml parameters.downsample_n_cells)"
DOWNSAMPLING_METHOD="$(read_yaml parameters.downsampling_method)"
NCORES="$(read_yaml parameters.ncores)"
SEED="$(read_yaml parameters.seed)"

CONTAINER_BULK_COUNTS_DIR="/data/bulk_counts"
CONTAINER_BULK_TPM_DIR="/data/bulk_TPM"
CONTAINER_CLINIC_DIR="/data/clinic"
CONTAINER_SC_DIR="/data/sc"
CONTAINER_SIGNATURES_DIR="/work/signatures"
CONTAINER_OUTPUT_DIR="/work/output"

check_paths() {
  [ -f "$HOST_BULK_COUNTS_DIR" ]      || { echo "Missing bulk_counts dir: $HOST_BULK_COUNTS_DIR"; exit 1; }
  [ -f "$HOST_BULK_TPM_DIR" ]      || { echo "Missing bulk_TPM dir: $HOST_BULK_TPM_DIR"; exit 1; }
  [ -d "$HOST_CLINIC_DIR" ]      || { echo "Missing clinic_dir: $HOST_CLINIC_DIR"; exit 1; }
  [ -f "$HOST_SINGLE_CELL_RDS" ] || { echo "Missing single_cell_rds: $HOST_SINGLE_CELL_RDS"; exit 1; }

  mkdir -p "$HOST_SIGNATURES_DIR"
  mkdir -p "$HOST_OUTPUT_DIR"

  [ -f "$HOST_CLINIC_DIR/$CLINIC_FILE" ] || {
    echo "Missing clinic file: $HOST_CLINIC_DIR/$CLINIC_FILE"
    exit 1
  }
}

build_image() {
  docker build --pull -t "$IMAGE_NAME" -f DockerFile .
}

run_preanalysis() {
  [ -f "$HOST_SINGLE_CELL_RDS" ] || { echo "Missing single_cell_rds: $HOST_SINGLE_CELL_RDS"; exit 1; }
  mkdir -p "$HOST_OUTPUT_DIR"

  # Mount the project dir so preanalysis.R runs from the existing image
  # (no rebuild needed); override the entrypoint to call Rscript directly.
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

  docker run --rm \
    -v "$(dirname "$HOST_SINGLE_CELL_RDS"):$CONTAINER_SC_DIR:ro" \
    -v "$HOST_OUTPUT_DIR:$CONTAINER_OUTPUT_DIR" \
    -v "$SCRIPT_DIR:/pipeline:ro" \
    --entrypoint Rscript \
    "$IMAGE_NAME" /pipeline/preanalysis.R \
    --single_cell_rds "$CONTAINER_SC_DIR/$(basename "$HOST_SINGLE_CELL_RDS")" \
    --celltype_col "$CELLTYPE_COL" \
    --patient_col "$PATIENT_COL" \
    --gse_id "$GSE_ID" \
    --bulk_ID "$BULK_ID" \
    --current_cap "$DOWNSAMPLE_N_CELLS" \
    --output_dir "$CONTAINER_OUTPUT_DIR"
}

run_container() {
  check_paths

  # ncores is optional: leave it NULL/blank in the YAML to let the container
  # auto-detect (all logical cores - 4). When set, pass it through explicitly.
  ncores_args=()
  if [ -n "$NCORES" ]; then
    ncores_args=(--ncores "$NCORES")
  fi

  docker run --rm -it \
    -v "$HOST_BULK_COUNTS_DIR:$CONTAINER_BULK_COUNTS_DIR:ro" \
    -v "$HOST_BULK_TPM_DIR:$CONTAINER_BULK_TPM_DIR:ro" \
    -v "$HOST_CLINIC_DIR:$CONTAINER_CLINIC_DIR:ro" \
    -v "$(dirname "$HOST_SINGLE_CELL_RDS"):$CONTAINER_SC_DIR:ro" \
    -v "$HOST_SIGNATURES_DIR:$CONTAINER_SIGNATURES_DIR" \
    -v "$HOST_OUTPUT_DIR:$CONTAINER_OUTPUT_DIR" \
    "$IMAGE_NAME" \
    --bulk_counts "$CONTAINER_BULK_COUNTS_DIR" \
    --bulk_TPM "$CONTAINER_BULK_TPM_DIR" \
    --clinic_dir "$CONTAINER_CLINIC_DIR" \
    --clinic_file "$CLINIC_FILE" \
    --single_cell_rds "$CONTAINER_SC_DIR/$(basename "$HOST_SINGLE_CELL_RDS")" \
    --gse_id "$GSE_ID" \
    --bulk_ID "$BULK_ID" \
    --signatures_dir "$CONTAINER_SIGNATURES_DIR" \
    --output_dir "$CONTAINER_OUTPUT_DIR" \
    --deconv_methods "$DECONV_METHODS" \
    --dwls_method "$DWLS_METHOD" \
    --celltype_col "$CELLTYPE_COL" \
    --patient_col "$PATIENT_COL" \
    --downsample_n_cells "$DOWNSAMPLE_N_CELLS" \
    --downsampling_method "$DOWNSAMPLING_METHOD" \
    --seed "$SEED" \
    "${ncores_args[@]}"
}

case "$MODE" in
  build)
    build_image
    ;;
  run)
    run_container
    ;;
  preanalysis)
    run_preanalysis
    ;;
  all)
    build_image
    run_container
    ;;
  *)
    echo "Usage: $0 [build|run|preanalysis|all] [config.yaml]"
    exit 1
    ;;
esac