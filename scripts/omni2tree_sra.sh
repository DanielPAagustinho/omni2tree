#!/usr/bin/env bash
set -euo pipefail

##############################################################################
# This script downloads SRA data, handling "run" (SRR, ERR, DRR) or
# "experiment" (SRX, ERX, DRX) accessions from a file, one SPECIES per line.
#
# It splits the accessions into chunks to retrieve bulk runinfo once, then:
#   - For RUNs, it finds LAYOUT from runinfo with grep (unless forced).
#   - For EXPERIMENTs, it grabs *all* RUN lines for that experiment,
#     and downloads each RUN individually.
#
##############################################################################

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

PROGNAME="$(basename "$0")"

log_info() {
  echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]${NC} $*"
}

log_warn() {
  echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]${NC} $*"
}

log_error() {
  echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR]${NC} $*" >&2
}

########################################
# Functions
########################################

check_dependencies() {
  #Initialize array to map tools to messages
  local missing=0
  declare -A tools=(
    ["prefetch"]="SRA Toolkit"
    ["efetch"]="Entrez Direct utilities - efetch"
    ["esearch"]="Entrez Direct utilities - esearch"

  )

  #log_info "Checking system dependencies..."
  #Looping through the array to test for commands and detect possible missing tools  
  for cmd in "${!tools[@]}"; do
    if ! command -v "$cmd" &>/dev/null; then
      log_error "Missing requirement: ${tools[$cmd]}"
      ((missing++))
    #else
      #log_info "Found: ${tools[$cmd]}"
    fi
  done
}

# For species dir name: allow only alnum and replace spaces with underscores
sanitize_name() {
	  shopt -s extglob 
	  local input="$1"
	  # Quitar espacios al inicio y final
	  input="${input##+([[:space:]])}"
	  input="${input%%+([[:space:]])}"
     	  input="${input//+([[:space:]])/_}"
	  input="${input//[^a-zA-Z0-9_]/}"
	  echo "$input"
 }


is_run_prefix() {
  [[ "$1" == SRR* || "$1" == ERR* || "$1" == DRR* ]]
}

is_experiment_prefix() {
  [[ "$1" == SRX* || "$1" == ERX* || "$1" == DRX* ]]
}

# Split an array into multiple lines (chunks), each up to $CHUNK_SIZE items
chunk_array() {
  local array_name="$1"
  local -n _arr="$array_name"

  local total="${#_arr[@]}"
  local start=0

  while (( start < total )); do
    local end=$(( start + CHUNK_SIZE ))
    if (( end > total )); then
      end="$total"
    fi
    local chunk=("${_arr[@]:start:end-start}")
    echo "${chunk[*]}"
    start="$end"
  done
}

usage() {
  log_info "Usage: $(basename "$0") -i <input_file> [-o <output_dir>] [options]\n"
  #echo "Try '$0 --help' for more information."
}

show_help() {
  cat <<EOF
$(usage)

Required:
  -i, --input         Input file containing SRA IDs with one taxon per line:
                      <species_name>,SRA_ID1,SRA_ID2,SRA_ID3,...
                      <species_name2>,SRA_ID4,SRA_ID5,SRA_ID6,...
                      Where <species_name> is the taxon name, and SRA_IDs for each taxon must be either all SRA RUNs (SRR, ERR, DRR) 
                      or all SRA EXPERIMENTs (SRX, ERX, DRX) per line.

Optional:
  -o, --outdir        Output directory (default: current dir)
  -c, --chunk-size    Number of SRA IDs in each chunk when fetching metadata using 
                      esearch and efetch (default: 350)
  -w, --sleep-secs    Seconds to sleep between chunks (default: 1)
  -l, --layout        Force layout (SINGLE or PAIRED) for all runs, if SRA IDs correspond to RUNS, it skips metadata fetching
  -d, --debug         Avoid removing of species temporary directory (default: off)
  -h, --help          Show help


Example:
  $PROGNAME -i species_accessions.txt -o results --chunk-size 100 --sleep-secs 2

EOF
}

########################################
# Defaults
########################################
INPUT_FILE=""
OUTPUT_DIR="$(pwd)"
CHUNK_SIZE=350
SLEEP_SECS=1
USER_LAYOUT=""    # If the user determines SINGLE or PAIRED
DEBUG=false
global_ok_runs=()      # todas las runs/experimentos exitosos
global_fail_runs=()    # todas las runs/experimentos fallidos
species_reports=()     # textos con el resumen de cada especie
########################################
# Parse arguments
########################################

if [[ $# -eq 0 ]]; then
  usage
  log_info "Try '$PROGNAME --help' for more information."
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      INPUT_FILE="$2"
      shift 2
      ;;
    -o|--outdir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    -c|--chunk-size)
      CHUNK_SIZE="$2"
      shift 2
      ;;
    -w|--sleep-secs)
      SLEEP_SECS="$2"
      shift 2
      ;;
    -l|--layout)
      USER_LAYOUT="$2"
      shift 2
      ;;
    -d|--debug)
      DEBUG=true
      shift
      ;;
    -h|--help)
      show_help
      exit 0
      ;;
    *)
      log_error "Unknown argument: $1"
      usage
      log_info "Try '$PROGNAME --help' for more information."
      exit 1
      ;;
  esac
done

check_dependencies
log_info "Checked system dependencies"

if [[ -z "$INPUT_FILE" ]]; then
  log_error "-i/--input is required."
  usage
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  log_error "input file '$INPUT_FILE' does not exist."
  exit 1
fi

if [[ -n "$USER_LAYOUT" ]]; then
  LAYOUT_UPPER="$(echo "$USER_LAYOUT" | tr '[:lower:]' '[:upper:]')"
  if [[ "$LAYOUT_UPPER" != "SINGLE" && "$LAYOUT_UPPER" != "PAIRED" ]]; then
    log_error "Forced layout must be SINGLE or PAIRED."
    exit 1
  fi
  USER_LAYOUT="$LAYOUT_UPPER"
fi

if [[ "$DEBUG" == true ]]; then
  log_info "Debug mode enabled, keeping species temporary directory"
fi

OUTPUT_DIR="${OUTPUT_DIR%/}"
mkdir -p "$OUTPUT_DIR"
SUMMARY_FILE="${OUTPUT_DIR}/summary_download.txt"
: > "$SUMMARY_FILE"   # truncar/crear archivo de resúmenes


########################################
# Main
########################################

mapfile -t lines < "$INPUT_FILE"

for line in "${lines[@]}"; do

  success_runs=()  fail_runs=()
  success_exp=()   fail_exp=()
  declare -A fail_by_exp=()
  # Skip empty or commented lines
  [[ -z "$line" || "$line" =~ ^# ]] && continue
  # Parse the line by commas, first part is the taxon names, and from the second on are the accessions
  IFS=',' read -ra fields <<< "$line"
  local_name="${fields[0]}"
  accessions=()
  for field in "${fields[@]:1}"; do
    clean_field="${field//[[:space:]]/}"  # Delete any spaces or tab
    accessions+=("$clean_field")
  done
  if (( ${#accessions[@]} == 0 )); then
    log_warn "No SRA IDs found for species '$local_name' in input. Check the CSV line."
    continue
  fi
  species_dir="$(sanitize_name "$local_name")"
  species_outdir="${OUTPUT_DIR}/${species_dir}"
  mkdir -p "$species_outdir"

  log_info "==========================================================="
  log_info "Processing species: $local_name"
  log_info "Output dir:         $species_outdir"
  log_info "SRA IDs:         ${accessions[*]}"

  # 1) Determine if RUN or EXPERIMENT by first accession
  first_acc="${accessions[0]}"
  if is_run_prefix "$first_acc"; then
    experiment_flag="false"
    #echo "Enabling RUN mode" 
  elif is_experiment_prefix "$first_acc"; then
    experiment_flag="true"
    #echo "Enabling EXPERIMENT mode" 
  else
    log_error "First SRA ID '$first_acc' is neither RUN nor EXPERIMENT."
    exit 1
  fi

  # 2) Validate that all accessions share are or run or experiment
  for a in "${accessions[@]}"; do
    if [[ "$experiment_flag" == "false" ]]; then
      if ! is_run_prefix "$a"; then
        log_error "Expected RUN prefix but found: $a"
        exit 1
      fi
    else
      if ! is_experiment_prefix "$a"; then
        log_error "Expected EXPERIMENT prefix but found: $a"
        exit 1
      fi
    fi
  done

  if [[ -z "$USER_LAYOUT" || "$experiment_flag" == "true" ]] ; then 
    # 3) Build a single runinfo CSV for all accessions of this species IF layout doesn't exist (variable void) or it exists, but it is a experiment
    tmp_runinfo="${species_outdir}/tmp_runinfo_${species_dir}.csv"
    : > "$tmp_runinfo"   
    log_info "Getting metadata for SRA IDs CHUNKS"
    while IFS= read -r chunk_line; do
      # chunk_line might be "ACC1 ACC2"
      IFS=' ' read -ra chunk_accs <<< "$chunk_line"
      # Build "acc1 OR acc2" ... #replacen spaces with parameter expansion. use printf instead of sed
      printf -v query_str '%s OR ' "${chunk_accs[@]}"  # Join elements with " OR "
      query_str="${query_str% OR }"
  
      log_info "  [CHUNK] Acc: ${chunk_accs[*]}"
      esearch -db sra -query "$query_str" < /dev/null \
        | efetch -format runinfo >> "$tmp_runinfo" 2>> "${species_outdir}/error_runinfo.log"

      log_info "  Sleeping $SLEEP_SECS second(s)..."
      sleep "$SLEEP_SECS"

    done < <(chunk_array accessions)
  else
    log_info "The user has forced the layout to be ${USER_LAYOUT} and it is RUN mode"
  fi
  # 4) Download step
  if [[ "$experiment_flag" == "false" ]]; then
    ################################################
    # RUN MODE
    ################################################
    log_info "RUN mode: each accession is already a RUN (SRR/ERR/DRR)."

    # For each run accession
    for acc in "${accessions[@]}"; do

      # Determine layout
      layout="UNKNOWN"
      if [[ -n "$USER_LAYOUT" ]]; then
        layout="$USER_LAYOUT"
      else
        # Try to grep for the acc in $tmp_runinfo
        run_line="$(grep "$acc" "$tmp_runinfo" || true)"
        if [[ -n "$run_line" ]]; then
          # parse col16
          raw_layout="$(echo "$run_line" | cut -d',' -f16 | tr '[:upper:]' '[:lower:]')"
          if [[ "$raw_layout" == "paired" ]]; then
            layout="PAIRED"
          elif [[ "$raw_layout" == "single" ]]; then
            layout="SINGLE"
          fi
        fi
      fi

      log_info "---------------------------------------------------"
      log_info " SRA ID (run): $acc"
      log_info " Layout:          $layout"

      # Download with prefetch
      log_info "   Downloading .sra with prefetch..."
      if ! prefetch --output-directory "${species_outdir}" "$acc" </dev/null; then
        log_warn "   Prefetch failed for $acc"
        fail_runs+=("$acc")
        continue
      fi

      # Look for the .sra
      SRA_PATH="${species_outdir}/${acc}"
      if [[ ! -d "$SRA_PATH" ]]; then
        log_warn "   SRA directory not found for $acc. Skipping..."
        continue
      fi

            # Convert to FASTQ
      log_info "   Converting to FASTQ..."
      cmd=(fasterq-dump)
      [[ "$layout" == "PAIRED" ]] && cmd+=(--split-files)
      if "${cmd[@]}" "$SRA_PATH" -O "$species_outdir" </dev/null; then
        moved=0
        shopt -s nullglob
        for ext in fastq fastq.gz fq fq.gz; do
            for f in "${species_outdir}/${acc}"*.${ext}; do
              if mv -n "$f" "${OUTPUT_DIR}/${species_dir}_$(basename "$f")"; then
                (( moved+=1 ))
              else
                log_warn "$f could not be moved to output directory."
              fi
            done
        done
        shopt -u nullglob
        if [[ $moved -eq 0 ]]; then
          log_warn "No FASTQ files found for ${acc}. Not moving to output directory."
        fi
     
        success_runs+=("$acc")
        log_info "   [OK] Finished $acc"
      else
        fail_runs+=("$acc")
        log_warn "   Conversion failed for $acc"
      fi

    done

  else
    ################################################
    # EXPERIMENT MODE
    ################################################
    log_info "EXPERIMENT mode: each accession can map to multiple RUNs."

    for exp_id in "${accessions[@]}"; do
      log_info "---------------------------------------------------"
      log_info " Experiment: $exp_id"

      # grep all lines with that experiment ID.
      # Often, experiment is in column #11 (Experiment).
      # But the order can vary; we rely on the fact that the experiment ID
      # typically appears *somewhere* on the line, so we do a substring grep.
      matched_lines="$(grep "$exp_id" "$tmp_runinfo" || true)"
      if [[ -z "$matched_lines" ]]; then
        log_warn "  No runs found for experiment '$exp_id' in $tmp_runinfo."
        continue
      fi

      # For each matching line, parse out col1 (Run) and col16 (Layout).
      # We can have multiple lines (multiple runs).
      while IFS= read -r line; do
        run_id="$(echo "$line" | cut -d',' -f1)"

        layout="UNKNOWN"
        if [[ -n "$USER_LAYOUT" ]]; then
          layout="$USER_LAYOUT"
        else
          raw_layout="$(echo "$line" | cut -d',' -f16 | tr '[:upper:]' '[:lower:]')"
          if [[ "$raw_layout" == "paired" ]]; then
            layout="PAIRED"
          elif [[ "$raw_layout" == "single" ]]; then
            layout="SINGLE"
          fi
        fi

        log_info "   => Found RUN: $run_id"
        log_info "      Layout:    $layout"

        # Download with prefetch
        log_info "   Downloading .sra with prefetch..."
        if ! prefetch --output-directory "${species_outdir}" "$run_id" </dev/null; then
          log_warn "      Prefetch failed for $run_id"
          fail_exp+=("$exp_id")
          fail_by_exp["$exp_id"]+="$run_id "
          continue
        fi
        # Locate .sra
        SRA_PATH="${species_outdir}/${run_id}"
        if [[ ! -d "$SRA_PATH" ]]; then
          log_warn "   SRA directory not found for $run_id. Skipping..."
          continue
        fi

        log_info "      Converting .sra to FASTQ..."
        cmd=(fasterq-dump)
        [[ "$layout" == "PAIRED" ]] && cmd+=(--split-files)
        if "${cmd[@]}" "$SRA_PATH" -O "$species_outdir" </dev/null; then
          moved=0
          shopt -s nullglob
          for ext in fastq fastq.gz fq fq.gz; do
              for f in "${species_outdir}/${run_id}"*.${ext}; do
                if mv -f "$f" "${OUTPUT_DIR}/${species_dir}_$(basename "$f")"; then
                  (( moved+=1 ))
                else
                  log_warn "$f could not be moved to output directory."
                fi
              done
          done
          shopt -u nullglob
          if [[ $moved -eq 0 ]]; then
            log_warn "No FASTQ files found for ${run_id}. Not moving to output directory."
          fi
          success_exp+=("$exp_id")
          success_runs+=("$run_id")   # <- para unificar métricas por RUN
          log_info "      Finished run $run_id for experiment $exp_id"
        else
          fail_exp+=("$exp_id")
          fail_runs+=("$run_id")      # <- para unificar métricas por RUN
          fail_by_exp["$exp_id"]+="$run_id "
          log_warn "      FASTQ conversion failed for $run_id"
        fi

      done <<< "$matched_lines"
    done
  fi

  log_info "Done with $local_name"
  if [[ "$DEBUG" == false ]]; then
    log_info "Removing intermediate taxon directory: $species_outdir"
    rm -rf -- "$species_outdir"
  fi
  log_info "### Summary for $local_name ###"
  log_info "successful RUNs:     ${#success_runs[@]}  -> ${success_runs[*]}"
  log_info "failed RUNs: ${#fail_runs[@]}      -> ${fail_runs[*]}"
  if [[ "$experiment_flag" == "true" ]]; then
    for eid in "${!fail_by_exp[@]}"; do
      log_warn "Exp $eid failed in RUNs: ${fail_by_exp[$eid]}"
    done
  fi
  species_report="### Summary for $local_name ###
  RUNs OK:     ${#success_runs[@]} -> ${success_runs[*]:-none}
  RUNs failed: ${#fail_runs[@]} -> ${fail_runs[*]:-none}"
  if [[ "$experiment_flag" == "true" && ${#fail_by_exp[@]} -gt 0 ]]; then
    for eid in "${!fail_by_exp[@]}"; do
      species_report+="
  Exp $eid failed in RUNs: ${fail_by_exp[$eid]}"
    done
  fi
  species_reports+=("$species_report")
  global_ok_runs+=( "${success_runs[@]}" )
  global_fail_runs+=( "${fail_runs[@]}" )
  {
  printf 'SRA IDs: %s\n' "${accessions[*]}"
  printf '%b\n\n' "$species_report"
  } >> "$SUMMARY_FILE"
done

log_info "### Global Summary ###"
log_info "Total RUNs OK:     $(( ${#global_ok_runs[@]} ))"
log_info "Total RUNs failed: $(( ${#global_fail_runs[@]} ))"
log_info "Processing complete. Outputs are in '$OUTPUT_DIR'."
###

log_info "Summaries saved to '$SUMMARY_FILE'."