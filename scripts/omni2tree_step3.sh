#!/usr/bin/env bash
set -euo pipefail

PROGNAME="$(basename "$0")"
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
MAIN_DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"   # .../omni2tree/scripts
REPO_ROOT="$(cd "$MAIN_DIR/.." && pwd)"              # .../omni2tree
UTILS_DIR="$REPO_ROOT/utils"
VIEW_DIR="$REPO_ROOT/view"

ROOT_DIR=""
OUT_DIR=""
THREADS=4
BOOTSTRAP=1000
SEQ_MODE="aa"
SEQ_LC="aa"
SEQ_UC="AA"
LABEL="Omni2tree_Tree"
METADATA_FILE=""
TEMP_DIR=""
DEBUG=false

EXCLUDE_PATTERN=""
FILTER_COLUMN=""
FILTER_VALUE=""
EXCLUDE_GAPS=false
MIN_SAMPLES=""
ADD_DOMAIN_FILE=""
GROUP_BY=()
FIVE_LETTER_FILE="five_letter_taxon.tsv"
OG_ENTROPY_TABLE="stats/entropy/OG_genes_entropy.csv"

if [ -t 1 ]; then
  RED="\033[1;31m"
  GREEN="\033[1;32m"
  YELLOW="\033[1;33m"
  NC="\033[0m"
else
  RED=""
  GREEN=""
  YELLOW=""
  NC=""
fi

log_info() {
  echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]${NC} $*"
}

log_warn() {
  echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]${NC} $*"
}

log_error() {
  echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR]${NC} $*" >&2
}

usage() {
  log_info "Usage: ${PROGNAME} --o2t_out <dir> -m <metadata.csv> [-l <label>] [options]"
}

show_help() {
  cat << EOF
$(usage)

Runs Omni2tree Step 3 end-to-end:
  1) read2tree --step 3combine
  2) iqtree tree inference
  3) visualization metadata/tree preparation + omni2treeview
  4) Shannon entropy pipeline (msa_to_position_table.py -> calculate_entropy.py -> plot_entropy.R)

Required:
  --o2t_out <dir>                               Base Omni2tree output directory that contains O2T_RESULTS
  -m, --metadata <file>                         Metadata CSV (header + first data row with column types)

General Optional:
  --seq_type <aa|dna>                           Selects concatenated alignment used by IQ-TREE and entropy input [default: aa]
  -T, --threads <int>                           Threads for iqtree [default: 4]
  --bootstrap <int>                             IQ-TREE ultrafast bootstrap replicates (0 or >=1000) [default: 1000]
  --temp_dir <dir>                              Optional temp dir. If relative, it will be relative to o2t_out.
  --debug                                       Keep temp dir and print extra messages
  -h, --help                                    Show this help

Visualization:
  -l, --label <text>                            Optional visualization label [default: Omni2tree_Tree]

Entropy Step 1 (msa_to_position_table.py) optional filters:
  --exclude_pattern <regex>                     Python regex to exclude sample IDs (case-sensitive)
  --filter_column <name>                        Metadata column to filter on (requires --filter_value)
  --filter_value <value>                        Metadata value to keep (requires --filter_column)

Entropy Step 2 (calculate_entropy.py) optional:
  --group_by <col1> [col2 ...]                  Metadata columns for grouped entropy (space-separated)
  --min_samples <int>                           Minimum samples per position (default: script default)
  --exclude_gaps                                Exclude '-' characters before entropy calculation

Entropy Step 3 (plot_entropy.R) optional:
  --add_domain <csv>                            Domain annotations CSV with columns: gene,domain,start,end
                                                Gene names are validated against stats/entropy/OG_genes_entropy.csv
                                                Note: currently validation-only; plot_entropy.R does not render domains yet.

Expected inputs from previous steps (inside --o2t_out):
  O2T_RESULTS/                                  read2tree output from step1/step2
  marker_genes/                                 marker genes directory from step1
  dna_ref.fa                                    read2tree DNA reference from step1
  five_letter_taxon.tsv                         mapping generated in step1
  stats/entropy/OG_genes_entropy.csv            generated in step1 from entropy_msa_og_gene_table.py

Examples:
  $PROGNAME --o2t_out virus2tree_rsv -m metadata.csv

  $PROGNAME --o2t_out virus2tree_rsv -l RSV_Run -m metadata.csv \\
    --seq_type aa --exclude_pattern s0 --group_by subgroup time_phase --exclude_gaps

  $PROGNAME --o2t_out virus2tree_rsv -l RSV_Run -m metadata.csv \\
    --seq_type dna --bootstrap 2000 --add_domain domains.csv

EOF
  exit 0
}

run_cmd() {
  local cmd_str=""
  printf -v cmd_str '%q ' "$@"
  cmd_str="${cmd_str% }"
  log_info "Executing command: $cmd_str"
  "$@"
}

clean_alnum() {
  local x="$1"
  echo "$x" | tr -cd '[:alnum:]'
}

sanitize_label_id() {
  local x="$1"
  x="${x%_1}"
  x="${x%_2}"
  x="$(echo "$x" | sed -E 's/[^[:alnum:]_]+/_/g; s/_+/_/g; s/^_+//; s/_+$//')"
  [[ -z "$x" ]] && x="NA"
  echo "$x"
}

require_file() {
  local path="$1"
  local msg="${2:-Required file not found or empty}"
  if [[ ! -f "$path" ]] || [[ ! -s "$path" ]]; then
    log_error "$msg: $path"
    exit 1
  fi
}

require_dir() {
  local path="$1"
  local msg="${2:-Required directory not found}"
  if [[ ! -d "$path" ]]; then
    log_error "$msg: $path"
    exit 1
  fi
}

check_r_packages() {
  if ! Rscript -e "quit(status=!(requireNamespace('tidyverse', quietly=TRUE) && requireNamespace('RColorBrewer', quietly=TRUE)))" >/dev/null 2>&1; then
    log_error "Missing R packages: tidyverse and/or RColorBrewer"
    log_error "Install with: Rscript -e \"install.packages(c('tidyverse','RColorBrewer'), repos='https://cran.rstudio.com/')\""
    exit 1
  fi
}

check_dependencies() {
  declare -A tools=(
    ["read2tree"]="read2tree"
    ["iqtree"]="IQ-TREE"
    ["python3"]="python 3"
    ["Rscript"]="Rscript"
    ["awk"]="AWK"
    ["grep"]="grep"
    ["sed"]="sed"
  )

  for cmd in "${!tools[@]}"; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
      log_error "Missing requirement: ${tools[$cmd]}"
      exit 1
    fi
  done

  for pkg in pandas numpy; do
    if ! python3 -c "import $pkg" >/dev/null 2>&1; then
      log_error "Missing python3 package: $pkg"
      exit 1
    fi
  done
  if ! python3 -c "import Bio" >/dev/null 2>&1; then
    log_error "Missing python3 package: biopython"
    exit 1
  fi

  check_r_packages

  require_file "$MAIN_DIR/msa_to_position_table.py" "Missing script"
  require_file "$MAIN_DIR/calculate_entropy.py" "Missing script"
  require_file "$MAIN_DIR/plot_entropy.R" "Missing script"
  require_file "$UTILS_DIR/validate_metadata.py" "Missing script"
  require_file "$UTILS_DIR/prepare_metadata_o2t_view.py" "Missing script"
  require_file "$VIEW_DIR/omni2treeview.py" "Missing script"
  require_file "$VIEW_DIR/template_v5.html" "Missing HTML template"
}

validate_domain_csv() {
  local domain_csv="$1"
  local og_table_csv="$2"
  python3 - "$domain_csv" "$og_table_csv" <<'PY'
import sys
import pandas as pd

domain_path, og_table_path = sys.argv[1], sys.argv[2]

try:
    d = pd.read_csv(domain_path)
except Exception as e:
    print(f"Failed to read domain CSV: {e}", file=sys.stderr)
    sys.exit(1)

required = {"gene", "domain", "start", "end"}
missing = required - set(d.columns.str.lower() if hasattr(d.columns, "str") else d.columns)
if missing:
    # normalize case-insensitive check
    lower_map = {c.lower(): c for c in d.columns}
    missing = required - set(lower_map)
    if missing:
        print(f"Domain CSV missing required columns: {', '.join(sorted(missing))}", file=sys.stderr)
        sys.exit(1)

lower_map = {c.lower(): c for c in d.columns}
d = d.rename(columns={lower_map[k]: k for k in lower_map})
d = d.rename(columns={"start": "start", "end": "end"})

for col in ("gene", "domain", "start", "end"):
    if col not in d.columns:
        print(f"Domain CSV missing normalized column: {col}", file=sys.stderr)
        sys.exit(1)

if d.empty:
    print("Domain CSV is empty", file=sys.stderr)
    sys.exit(1)

try:
    d["start"] = pd.to_numeric(d["start"])
    d["end"] = pd.to_numeric(d["end"])
except Exception:
    print("Columns 'start' and 'end' must be numeric", file=sys.stderr)
    sys.exit(1)

if ((d["start"] <= 0) | (d["end"] <= 0) | (d["end"] < d["start"])).any():
    print("Invalid domain intervals detected (require start>0, end>0, end>=start)", file=sys.stderr)
    sys.exit(1)

og = pd.read_csv(og_table_path)
if "gene" not in og.columns:
    print("OG table does not contain 'gene' column", file=sys.stderr)
    sys.exit(1)
valid_genes = set(og["gene"].astype(str))
bad_genes = sorted(set(d["gene"].astype(str)) - valid_genes)
if bad_genes:
    print(
        "Domain CSV contains gene(s) not present in stats/entropy/OG_genes_entropy.csv: "
        + ", ".join(bad_genes),
        file=sys.stderr,
    )
    sys.exit(1)
PY
}

validate_og_names_vs_msa() {
  local msa_dir="$1"
  local og_table_csv="$2"
  local tmp_msa="$TEMP_DIR/msa_ogs.txt"
  local tmp_tbl="$TEMP_DIR/table_ogs.txt"

  : > "$tmp_msa"
  shopt -s nullglob
  for f in "$msa_dir"/OG*.fa; do
    basename "${f%.fa}" >> "$tmp_msa"
  done
  shopt -u nullglob

  if [[ ! -s "$tmp_msa" ]]; then
    log_error "No OG*.fa files found in MSA directory: $msa_dir"
    exit 1
  fi

  awk -F',' 'NR>1 {print $1}' "$og_table_csv" | sort -u > "$tmp_tbl"
  sort -u "$tmp_msa" -o "$tmp_msa"

  local missing_in_table extra_in_table
  missing_in_table="$(comm -23 "$tmp_msa" "$tmp_tbl" || true)"
  extra_in_table="$(comm -13 "$tmp_msa" "$tmp_tbl" || true)"

  if [[ -n "$missing_in_table" ]]; then
    log_error "OG names present in $msa_dir but missing in $og_table_csv:"
    while IFS= read -r og; do [[ -n "$og" ]] && log_error "  - $og"; done <<< "$missing_in_table"
    exit 1
  fi

  if [[ -n "$extra_in_table" ]]; then
    log_warn "Some OGs in $og_table_csv are not present in $msa_dir (continuing)."
  fi
  log_info "Validated OG names: MSA files are covered by $OG_ENTROPY_TABLE"
}

detect_concat_phy() {
  local seq_lc="$1"
  local candidates=()

  shopt -s nullglob
  for f in "$OUT_DIR"/output/concat_*_"${seq_lc}".phy "$OUT_DIR"/concat_*_"${seq_lc}".phy; do
    [[ -f "$f" ]] && candidates+=("$f")
  done
  shopt -u nullglob

  if [[ ${#candidates[@]} -eq 0 ]]; then
    log_error "Could not find concatenated alignment for seq_type=${seq_lc} (*.phy) under $OUT_DIR"
    exit 1
  fi

  if [[ ${#candidates[@]} -gt 1 ]]; then
    log_warn "Multiple concatenated alignments found. Using first one: ${candidates[0]}" >&2
  fi
  echo "${candidates[0]}"
}

if [[ $# -eq 0 ]]; then
  usage
  log_info "Try '$PROGNAME --help' for more information."
  exit 1
fi

log_info "Script invoked with: $PROGNAME $*\n"
log_info "========== Step 3.1: Validating parameters =========="

while [[ $# -gt 0 ]]; do
  case "$1" in
    --o2t_out)
      ROOT_DIR="${2%/}"
      shift 2
      ;;
    -l|--label)
      LABEL="$2"
      shift 2
      ;;
    -m|--metadata)
      METADATA_FILE="$2"
      shift 2
      ;;
    --seq_type)
      SEQ_MODE="${2,,}"
      shift 2
      ;;
    -T|--threads)
      THREADS="$2"
      shift 2
      ;;
    --bootstrap)
      BOOTSTRAP="$2"
      shift 2
      ;;
    --temp_dir)
      TEMP_DIR="${2%/}"
      shift 2
      ;;
    --debug)
      DEBUG=true
      shift
      ;;
    --exclude_pattern)
      EXCLUDE_PATTERN="$2"
      shift 2
      ;;
    --filter_column)
      FILTER_COLUMN="$2"
      shift 2
      ;;
    --filter_value)
      FILTER_VALUE="$2"
      shift 2
      ;;
    --group_by)
      shift
      if [[ $# -eq 0 || "$1" =~ ^- ]]; then
        log_error "--group_by requires at least one column name"
        exit 1
      fi
      while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
        GROUP_BY+=("$1")
        shift
      done
      ;;
    --min_samples)
      MIN_SAMPLES="$2"
      shift 2
      ;;
    --exclude_gaps)
      EXCLUDE_GAPS=true
      shift
      ;;
    --add_domain)
      ADD_DOMAIN_FILE="$2"
      shift 2
      ;;
    -h|--help)
      show_help
      ;;
    *)
      log_error "Unknown parameter passed: $1"
      usage
      log_info "Try '$PROGNAME --help' for more information."
      exit 1
      ;;
  esac
done

if [[ -z "$ROOT_DIR" ]]; then
  log_error "The user must specify --o2t_out <dir> (base directory containing O2T_RESULTS)"
  exit 1
fi
if [[ -z "$METADATA_FILE" ]]; then
  log_error "--metadata (-m) is required for visualization and entropy metadata-aware steps"
  exit 1
fi

case "$SEQ_MODE" in
  aa) SEQ_LC="aa"; SEQ_UC="AA" ;;
  dna) SEQ_LC="dna"; SEQ_UC="DNA" ;;
  *)
    log_error "Invalid --seq_type '$SEQ_MODE'. Use 'aa' or 'dna'"
    exit 1
    ;;
esac

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
  log_error "--threads must be a positive integer"
  exit 1
fi
if ! [[ "$BOOTSTRAP" =~ ^[0-9]+$ ]]; then
  log_error "--bootstrap must be a non-negative integer"
  exit 1
fi
if (( BOOTSTRAP > 0 && BOOTSTRAP < 1000 )); then
  log_error "--bootstrap must be 0 or >= 1000 for IQ-TREE ultrafast bootstrap"
  exit 1
fi
if [[ -n "$MIN_SAMPLES" ]] && ! [[ "$MIN_SAMPLES" =~ ^[1-9][0-9]*$ ]]; then
  log_error "--min_samples must be a positive integer"
  exit 1
fi
if [[ -n "$FILTER_COLUMN" && -z "$FILTER_VALUE" ]] || [[ -z "$FILTER_COLUMN" && -n "$FILTER_VALUE" ]]; then
  log_error "--filter_column and --filter_value must be provided together"
  exit 1
fi

check_dependencies
log_info "Checked system dependencies"

require_file "$METADATA_FILE" "Metadata file not found or empty"
METADATA_FILE="$(realpath "$METADATA_FILE")"

mkdir -p "$ROOT_DIR"
cd "$ROOT_DIR"
ROOT_DIR="$(pwd -P)"
OUT_DIR="$ROOT_DIR/O2T_RESULTS"
log_info "Using root directory: $ROOT_DIR"
log_info "Using read2tree output directory: $OUT_DIR"

require_dir "$OUT_DIR" "Expected step1/step2 read2tree output directory not found"
require_dir "marker_genes" "Required marker_genes directory not found"
require_file "dna_ref.fa" "Required dna_ref.fa not found"
require_file "five_letter_taxon.tsv" "Required five_letter_taxon.tsv not found"
require_file "$OG_ENTROPY_TABLE" "Required entropy OG mapping file not found (step1 should generate it)"

if [[ -z "$TEMP_DIR" ]]; then
  TEMP_DIR="$(mktemp -d)"
  log_info "Created temp directory at '$TEMP_DIR'"
else
  mkdir -p "$TEMP_DIR"
  TEMP_DIR="$(realpath "$TEMP_DIR")"
  log_info "Using temp directory: '$TEMP_DIR'"
fi

if [[ "$DEBUG" == false ]]; then
  trap '[[ -n "${TEMP_DIR:-}" && -d "${TEMP_DIR:-}" ]] && rm -rf "$TEMP_DIR"' EXIT
else
  log_info "Debug mode enabled, keeping temporary directory: '$TEMP_DIR'"
fi

if [[ -n "$ADD_DOMAIN_FILE" ]]; then
  require_file "$ADD_DOMAIN_FILE" "Domain CSV not found or empty"
  ADD_DOMAIN_FILE="$(realpath "$ADD_DOMAIN_FILE")"
  validate_domain_csv "$ADD_DOMAIN_FILE" "$OG_ENTROPY_TABLE" || {
    log_error "Domain CSV validation failed"
    exit 1
  }
  log_info "Validated domain annotations file: $ADD_DOMAIN_FILE"
fi

log_info "========== Step 3.1b: Validating metadata (early) =========="
VALIDATE_META_CMD=(python3 "$UTILS_DIR/validate_metadata.py"
  -m "$METADATA_FILE"
  --five_letter "$FIVE_LETTER_FILE"
  --o2t_results "$OUT_DIR")
run_cmd "${VALIDATE_META_CMD[@]}"

log_info "========== Step 3.2: Running Read2tree (step 3 combine) =========="
READ2TREE_CMD=(read2tree --step 3combine --standalone_path marker_genes --dna_reference dna_ref.fa --output_path "$OUT_DIR" --tree --debug)
run_cmd "${READ2TREE_CMD[@]}"

log_info "========== Step 3.3: Running IQ-TREE =========="
CONCAT_PHY="$(detect_concat_phy "$SEQ_LC")"
require_file "$CONCAT_PHY" "Concatenated alignment not found"
IQTREE_CMD=(iqtree -T "$THREADS" -s "$CONCAT_PHY" -bb "$BOOTSTRAP")
IQTREE_CMD+=( -pre "$CONCAT_PHY" )
run_cmd "${IQTREE_CMD[@]}"

TREEFILE_IQTREE="${CONCAT_PHY}.treefile"
TREEFILE_FALLBACK="$OUT_DIR/tree_merge.nwk"
VIEW_TREE_INPUT="$TREEFILE_IQTREE"
if [[ ! -s "$VIEW_TREE_INPUT" ]]; then
  if [[ -s "$TREEFILE_FALLBACK" ]]; then
    log_warn "IQ-TREE .treefile not found; using fallback tree: $TREEFILE_FALLBACK"
    VIEW_TREE_INPUT="$TREEFILE_FALLBACK"
  else
    log_error "Neither IQ-TREE treefile nor fallback tree_merge.nwk is available"
    exit 1
  fi
fi

log_info "========== Step 3.4: Preparing visualization inputs =========="
mkdir -p stats visualization
VIEW_TREE_OUTPUT="$OUT_DIR/concat_merge_view_${SEQ_LC}.phy.treefile"
VIEW_META_OUTPUT="visualization/metadata_o2t_view.csv"
PREP_META_CMD=(python3 "$UTILS_DIR/prepare_metadata_o2t_view.py"
  -m "$METADATA_FILE"
  --five_letter "$FIVE_LETTER_FILE"
  --in_nwk "$VIEW_TREE_INPUT"
  --out_nwk "$VIEW_TREE_OUTPUT"
  --out_meta "$VIEW_META_OUTPUT")
run_cmd "${PREP_META_CMD[@]}"

LABEL_SAFE="$(sanitize_label_id "$LABEL")"
VIEW_OUT_PREFIX="visualization/omni2treeview_${LABEL_SAFE}"
# NOTE (intentional, no changes applied to omni2treeview.py):
# - The generated metadata uses sample_id=accession and keeps label/source as metadata variables.
# - omni2treeview.py matches tree nodes mainly via sample_id (and may treat 'label' as a regular variable).
# - If the tree is relabeled with long labels, omni2treeview.py can still truncate/transform node names
#   depending on its internal leaf-name parsing. The viewer code itself is NOT modified here.
# - What would need changing (if desired later): omni2treeview.py matching logic and/or leaf-name parsing.
# - What is NOT affected: step3 combine/iqtree/entropy table generation.
log_warn "Viewer invocation note: omni2treeview.py is intentionally unmodified."
log_warn "With sample_id=accession, 'label' may be treated as metadata (not unique ID), and relabeled Newick names may be truncated by viewer parsing."
OMNIVIEW_CMD=(python3 "$VIEW_DIR/omni2treeview.py"
  -n "$VIEW_TREE_OUTPUT"
  -m "$VIEW_META_OUTPUT"
  -t "$VIEW_DIR/template_v5.html"
  -l "$LABEL"
  -o "$VIEW_OUT_PREFIX")
run_cmd "${OMNIVIEW_CMD[@]}"

log_info "========== Step 3.5: Entropy Step 1 (MSA -> position table) =========="
MSA_DIR="$OUT_DIR/06_align_merge_${SEQ_LC}"
require_dir "$MSA_DIR" "MSA directory not found for selected seq_type"
validate_og_names_vs_msa "$MSA_DIR" "$OG_ENTROPY_TABLE"

ENTROPY_DIR="stats/entropy"
mkdir -p "$ENTROPY_DIR"
POSITIONS_CSV="$ENTROPY_DIR/${SEQ_LC}_positions.csv"
ENTROPY_CSV="$ENTROPY_DIR/${SEQ_LC}_entropy.csv"
ENTROPY_PLOTS_DIR="$ENTROPY_DIR"

MSA2POS_CMD=(python3 "$MAIN_DIR/msa_to_position_table.py"
  --msa_dir "$MSA_DIR"
  --og_table "$OG_ENTROPY_TABLE"
  --output "$POSITIONS_CSV"
  --seq_type "$SEQ_UC"
  --metadata "$VIEW_META_OUTPUT"
  --metadata_match_column "label"
  --five_letter "$FIVE_LETTER_FILE")

if [[ -n "$EXCLUDE_PATTERN" ]]; then
  MSA2POS_CMD+=(--exclude_pattern "$EXCLUDE_PATTERN")
fi
if [[ -n "$FILTER_COLUMN" ]]; then
  MSA2POS_CMD+=(--filter_column "$FILTER_COLUMN" --filter_value "$FILTER_VALUE")
fi

run_cmd "${MSA2POS_CMD[@]}"

log_info "========== Step 3.6: Entropy Step 2 (Shannon entropy) =========="
CALC_ENTROPY_CMD=(python3 "$MAIN_DIR/calculate_entropy.py"
  --input "$POSITIONS_CSV"
  --output "$ENTROPY_CSV"
  --metadata "$VIEW_META_OUTPUT")

if [[ ${#GROUP_BY[@]} -gt 0 ]]; then
  CALC_ENTROPY_CMD+=(--group_by "${GROUP_BY[@]}")
fi
if [[ -n "$MIN_SAMPLES" ]]; then
  CALC_ENTROPY_CMD+=(--min_samples "$MIN_SAMPLES")
fi
if [[ "$EXCLUDE_GAPS" == true ]]; then
  CALC_ENTROPY_CMD+=(--exclude_gaps)
fi

run_cmd "${CALC_ENTROPY_CMD[@]}"

if ! awk 'NR>1 {found=1; exit 0} END {exit found ? 0 : 1}' "$ENTROPY_CSV"; then
  log_error "Entropy output contains no data rows after filtering. Adjust --filter_column/--filter_value and/or lower --min_samples."
  exit 1
fi

log_info "========== Step 3.7: Entropy Step 3 (Plot entropy) =========="
# NOTE (intentional, no changes applied to plot_entropy.R):
# - plot_entropy.R remains unchanged and does not consume the optional domain CSV.
# - If --add_domain is provided, the file is validated earlier but NOT passed to Rscript.
# - What would need changing (if desired later): plot_entropy.R argument parsing + domain overlay logic.
# - What is NOT affected: msa_to_position_table.py and calculate_entropy.py outputs.
log_warn "plot_entropy.R invocation note: R script is intentionally unmodified and will not process --add_domain CSV."
PLOT_CMD=(Rscript "$MAIN_DIR/plot_entropy.R" "$ENTROPY_CSV" "$ENTROPY_PLOTS_DIR" "$SEQ_UC")
if [[ -n "$ADD_DOMAIN_FILE" ]]; then
  log_warn "--add_domain was provided ($ADD_DOMAIN_FILE), but it is ignored because plot_entropy.R was intentionally left unchanged."
fi
run_cmd "${PLOT_CMD[@]}"

log_info "Step 3 completed successfully."
log_info "Main outputs:"
log_info "  Tree/alignment: $CONCAT_PHY"
log_info "  Visualization HTML: ${VIEW_OUT_PREFIX}.html"
log_info "  Visualization metadata: $VIEW_META_OUTPUT"
log_info "  Entropy positions: $POSITIONS_CSV"
log_info "  Entropy table: $ENTROPY_CSV"
log_info "  Entropy plots dir: $ENTROPY_PLOTS_DIR"
