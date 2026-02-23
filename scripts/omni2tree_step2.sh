#!/usr/bin/env bash
set -euo pipefail
#set -x #bug fixes

############################################
# Default values
############################################
PROGNAME="$(basename "$0")"
READ_TYPE="single"
MINIMAP2_OPTIONS="-ax map-ont"
READS=()         # Array to store multiple input reads
TEMP_DIR=""      # If empty, we will create a new one with mktemp
OUT_DIR=""
DEDUP=false
DEDUP_L=""       # optional argument for czid-dedup
DOWNSAMPLE=false
# Downsampling method can be coverage & genome_size, or number_of_bases, or number_of_reads
COVERAGE=0.0
GENOME_SIZE=""
NUM_BASES=0
NUM_READS=0
THREADS=4        # you can adjust a default or parse it
DEBUG=false
STATS_FILE=""
ROOT_DIR=""
STATS_DIR=""

############################################
# Functions
############################################
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

check_dependencies() {
  #Initialize array to map tools to messages
  if ! command -v read2tree &>/dev/null; then
    log_error "Missing requirement: read2tree"
    exit 1
  fi
  if [[ "$DEDUP" == true ]]; then
    if ! command -v czid-dedup &>/dev/null; then
      log_error "Missing requirement: czid-dedup"
      exit 1
    fi
  fi
  if [[ "$DOWNSAMPLE" == true ]]; then
    if ! command -v rasusa &>/dev/null; then
      log_error "Missing requirement: rasusa (because --downsample was requested)"
      exit 1
    fi
  fi
  if ! command -v bc &>/dev/null; then
  log_error "Missing requirement: bc command"
  exit 1
  fi
  if ! command -v flock &>/dev/null; then
  log_error "Missing requirement: flock"
  exit 1
  fi


}

usage() {
  log_info "Usage: ${PROGNAME} -r <read_file1> [read_file2...] -t <single|paired> [-map_op|--minimap2_options \"<opts>\"] [options]\n"
  #echo "Try '$0 --help' for more information."
}

show_help() {
  cat << EOF
$(usage)

  Required:
    -r, --reads <file1> [file2 file3 ...]            Input fastq/fastq.gz

  Optional:
    -t, --read_type <paired|single>    Generic type of reads [default: single]
    -map_op, --minimap2_options <options>   Options for minimap2 when mapping read set to the reference [default: "-ax map-ont"]
                                            Suggested: -ax map-sr (for short single or paired end short reads), 
                                            -ax map-pb (for PacBio reads), -ax map-hifi (for HiFi reads)
    -T, --threads <int>      Threads to use [default: 4]
    --temp_dir <dir>         Specify temp directory (otherwise mktemp -d is used)
    --root_dir <dir>         Specify root directory containing step 1 result; all the outputs will be saved here [default: current directory]
    --out_dir <dir>          Specify read2tree step 1 output directory (absolute or relative to --root_dir). Must exist from step 1 [default: root_dir/read2tree_output]
    --stats_file <file>      Specify read statistics output file (absolute or relative to --root_dir) [default: root_directory/stats/reads_statistics.tsv]
    --dedup                  Run czid-dedup
    --dedup_l <int>          Provide '-l <int>' to czid-dedup (requires --dedup)
    --downsample             Run rasusa
    --coverage <float>       Coverage to target (requires --genome_size and --downsample)
    --genome_size <size>     Genome size to calculate coverage with respect to. E.g., 4.3kb, 7Tb, 9000, 4.1MB (requires --coverage and --downsample)
    --num_bases <int>        Target number of bases (requires --downsample)
    --num_reads <int>        Target number of reads (requires --downsample)
    --debug                  Keep temporary directory with intermediate read preprocessing files
    -h, --help               Show this help

Examples:
  $PROGNAME -t single -r sample.fastq.gz --dedup --downsample --num_reads 10000 --minimap2_options "-ax map-ont"
  $PROGNAME --read_type paired --reads R1.fastq R2.fastq --coverage 30 --genome_size 5MB --minimap2_options "-ax sr"

EOF
  exit 0
}

append_stats() {
  local stats_block="$1"
  (
    # Acquire an exclusive lock on them very same stats file (200 is the file descriptor)
    flock -x 200
    if [[ ! -s "$STATS_FILE" ]]; then
      printf "%s\t%s\t%s\t%s\n" "FASTQ File" "Num Reads" "Avg Length" "Total Bases" >> "$STATS_FILE"
    fi
    printf "%s\n" "$stats_block" >> "$STATS_FILE"
  ) 200>>"$STATS_FILE"
}

#function to analyze each fastq
analyze_fastq() {
    local fastq_file="$1"
    #local output_file="$2"
    local total_lines total_reads total_bases avg_length

    # Verify if the file is compressed
    if [[ "$fastq_file" == *.gz ]]; then
        total_lines=$(zcat "$fastq_file" | wc -l)
        total_bases=$(zcat "$fastq_file" | awk 'NR % 4 == 2 {total += length($0)} END {print total}')
    else
        total_lines=$(wc -l < "$fastq_file")
        total_bases=$(awk 'NR % 4 == 2 {total += length($0)} END {print total}' "$fastq_file")
    fi

    # Each read correspond to four lines
    total_reads=$((total_lines / 4))

    #Get mean length
      if [[ "$total_reads" -gt 0 ]]; then
          avg_length=$((total_bases / total_reads))
      else
          avg_length=0
      fi

      # Print results in tabular format
      printf "%s\t%'d\t%'d\t%'d" "$(basename "$fastq_file")" "$total_reads" "$avg_length" "$total_bases"
  }

############################################
# Parse arguments
############################################

if [[ $# -eq 0 ]]; then
  usage
  log_info "Try '$PROGNAME --help' for more information."
  exit 1
fi

log_info "Script invoked with: $PROGNAME $*\n"

log_info "========== Step 2.1: Validating parameters =========="

#Not positional arguments different from the options below
while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--read_type)
      READ_TYPE="$2"
      shift 2
      ;;
    -map_op|--minimap2_options)
      MINIMAP2_OPTS_ARR=()
      shift
      while [[ $# -gt 0 ]]; do
        case "$1" in
          -t|--read_type|-r|--reads|-T|--threads|--temp_dir|--root_dir|--out_dir|--stats_file|\
          --dedup|--dedup_l|--downsample|--coverage|--genome_size|--num_bases|--num_reads|\
          --debug|-h|--help|-map_op|--minimap2_options)
            break;;
          *)
            MINIMAP2_OPTS_ARR+=("$1"); shift;;
        esac
      done
      MINIMAP2_OPTIONS="${MINIMAP2_OPTS_ARR[*]}"
      ;;
    -r|--reads)
      shift
      while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
        READS+=("$1")
        shift
      done
      ;;
    -T|--threads)
      THREADS="$2"
      shift 2
      ;;
    --temp_dir)
      TEMP_DIR="${2%/}"
      shift 2
      ;;
    --root_dir)
      ROOT_DIR="${2%/}"
      shift 2
      ;;
    --out_dir)
      OUT_DIR="${2%/}"
      shift 2
      ;;
    --stats_file)
      STATS_FILE="$2"
      shift 2
      ;;
    --dedup)
      DEDUP=true
      shift
      ;;
    --dedup_l)
      DEDUP_L="$2"
      shift 2
      ;;
    --downsample)
      DOWNSAMPLE=true
      shift
      ;;
    --coverage)
      COVERAGE="$2"
      shift 2
      ;;
    --genome_size)
      GENOME_SIZE="$2"
      shift 2
      ;;
    --num_bases)
      NUM_BASES="$2"
      shift 2
      ;;
    --num_reads)
      NUM_READS="$2"
      shift 2
      ;;
    --debug)
      DEBUG=true
      shift
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

############################################
# Validate required parameters
############################################

check_dependencies
log_info "Checked system dependencies"


if [[ ${#READS[@]} -eq 0 ]]; then
  log_error "Please specify at least one input file with --reads"
  usage
  exit 1
else
  for file in "${READS[@]}"; do
    if [[ ! -f "$file" ]] || [[ ! -s "$file" ]]; then
      log_error "File not found or is empty: $file"
      exit 1
    fi
  done
fi

for i in "${!READS[@]}"; do
  READS[$i]="$(realpath "${READS[$i]}")"
done

case "$READ_TYPE" in
  single|paired) ;;
  *)
    log_error "Invalid --read_type: $READ_TYPE. Must be 'single' or 'paired'."
    exit 1
    ;;
esac

if [[ "$READ_TYPE" == "paired" && ${#READS[@]} -ne 2 ]]; then
  log_error "Paired-end reads requires exactly 2 files. You provided ${#READS[@]}"
  usage
  exit 1
fi

if [[ -n "$DEDUP_L" ]]; then
    if [[ "$DEDUP" != true ]]; then
        log_error "--dedup is required when using --dedup_l."
        exit 1
    fi
    if [[ "$DEDUP_L" -le 0 ]]; then
        log_error "Argument --dedup_l must be a positive integer"
        exit 1
    fi
fi

#Verify the use of just one of the 3 rasusa options
if [[ "$DOWNSAMPLE" == true ]]; then
  method_count=0
  #echo "$GENOME_SIZE"
  if [[ $(echo "$COVERAGE > 0" | bc -l) -eq 1 && -n "$GENOME_SIZE" ]]; then
    #echo "here in if"
    (( method_count+=1 )) 
    #echo "here after adding 1"
  fi
  if [[ $NUM_BASES -gt 0 ]]; then
    (( method_count+=1 ))
  fi
  if [[ $NUM_READS -gt 0 ]]; then
    (( method_count+=1 ))
  fi
  if [[ $method_count -eq 0 ]]; then
    log_error "Downsampling requires --coverage + --genome_size OR --num_bases OR --num_reads."
    exit 1
  fi
  if [[ $method_count -gt 1 ]]; then
    log_error "Please specify only ONE downsampling method (coverage+genome_size) OR num_bases OR num_reads."
    exit 1
  fi
else
  if [[ $(echo "$COVERAGE != 0" | bc -l) -eq 1 || -n "$GENOME_SIZE" || $NUM_BASES -gt 0 || $NUM_READS -gt 0 ]]; then
    log_error "Downsampling parameters provided without specifying --downsample. Please include --downsample."
    exit 1
  fi
fi

############################################
# Create or verify root, r2t and temp directories
############################################

if [[ -n "$ROOT_DIR" ]]; then
  if [[ ! -d "$ROOT_DIR" ]]; then
    mkdir -p "$ROOT_DIR"
  fi
  cd "$ROOT_DIR"
  ROOT_DIR="$(pwd -P)"
  log_info "Using root directory: $ROOT_DIR"
else
  ROOT_DIR="$(pwd -P)"
  log_info "No --root_dir provided. Using current directory as root: $ROOT_DIR"
fi

if [[ -z "$OUT_DIR" ]]; then
  OUT_DIR="$ROOT_DIR/read2tree_output"
  log_info "No --out_dir specified; assuming $OUT_DIR"
else
  [[ "$OUT_DIR" = /* ]] || OUT_DIR="$ROOT_DIR/$OUT_DIR"
  if [[ -d "$OUT_DIR" ]]; then
    OUT_DIR="$(realpath "$OUT_DIR")"
  else
    log_error "read2tree output directory not found: $OUT_DIR. Please run step 1 first."
    exit 1
  fi
fi
if [[ -z "$TEMP_DIR" ]]; then
  TEMP_DIR="$(mktemp -d)"
  # TEMP_DIR_CREATED=true
  log_info "Created temp directory at $TEMP_DIR"
else
  # Validate if it is a directory
  if mkdir "$TEMP_DIR" 2>/dev/null; then
    # TEMP_DIR_CREATED=true
    log_info "Using temp directory: $TEMP_DIR"
  else
    log_error "Temp directory already exists: $TEMP_DIR. Refuse to continue."
    exit 1
  fi
fi

# Set default stats file if not provided (using an absolute path)
if [[ -z "$STATS_FILE" ]]; then
  STATS_DIR="$ROOT_DIR/stats"
  mkdir -p "$STATS_DIR"
  STATS_FILE="$STATS_DIR/reads_statistics.tsv"
  log_info "No --stats_file specified; using $STATS_FILE"
else
  STATS_DIR="$(dirname "$STATS_FILE")"
  mkdir -p "$STATS_DIR"
  STATS_FILE="$(realpath "$STATS_FILE")"
fi


if [[ "$DEBUG" == false ]]; then
  trap '[[ -n "$TEMP_DIR" && -d "$TEMP_DIR" ]] && rm -rf "$TEMP_DIR"' EXIT
else
  log_info "Debug mode enabled, keeping temporary directory: $TEMP_DIR"
fi

############################################
# Decompression and Concatenation if necessary
############################################

log_info "========== Step 2.2: File decompression, concatenation and moving to the temporary directory =========="

FINAL_READS=()
if [[ ${#READS[@]} -gt 1 ]]; then
  if [[ ${READ_TYPE} == "paired" ]]; then
    for READ_SAMPLE in "${READS[@]}"; do
      filename=$(basename "$READ_SAMPLE")
      output_file=$TEMP_DIR/${filename%.gz}
      if gzip -t "$READ_SAMPLE" 2>/dev/null; then
        zcat "$READ_SAMPLE" > "$output_file"
      else
        cat "$READ_SAMPLE" > "$output_file"
      fi
      FINAL_READS+=("$output_file")
    done
  else
    log_info "Concatenating multiple files..."
    filename=$(basename "${READS[0]}")
    output_file=$TEMP_DIR/${filename%.gz}
    for READ_SAMPLE in "${READS[@]}"; do
      if gzip -t "$READ_SAMPLE" 2>/dev/null; then
        zcat "$READ_SAMPLE" >> "$output_file"
      else
        cat "$READ_SAMPLE" >> "$output_file"
      fi
    done
    FINAL_READS=("$output_file")
  fi
else
  filename=$(basename "${READS[0]}")
  output_file=$TEMP_DIR/${filename%.gz}
  if gzip -t "${READS[0]}" 2>/dev/null; then
    zcat "${READS[0]}" > "$output_file"
  else
    cat "${READS[0]}" > "$output_file"
  fi
  FINAL_READS=("$output_file")
fi

############################################
# Deduplication with czid-dedup (optional) 
############################################


if [[ "$DEDUP" == true ]]; then
  log_info "========== Step 2.3: Deduplication with czid-dedup (optional) =========="
  log_info "Preparing input for czid-dedup..."
  if [[ "$READ_TYPE" == "paired" ]]; then
    input_file_1="${FINAL_READS[0]}"
    input_file_2="${FINAL_READS[1]}"
    filename_1=$(basename "$input_file_1")
    filename_2=$(basename "$input_file_2")
    output_file_1="$TEMP_DIR/${filename_1%.fastq}_dedup.fastq"
    output_file_2="$TEMP_DIR/${filename_2%.fastq}_dedup.fastq"
    DEDUP_CMD=(czid-dedup -i "$input_file_1" -i "$input_file_2" -o "$output_file_1" -o "$output_file_2")
    if [[ -n "$DEDUP_L" ]]; then
      DEDUP_CMD+=(-l "$DEDUP_L")
    fi
    log_info "Executing czid-dedup command: ${DEDUP_CMD[*]}"
    "${DEDUP_CMD[@]}"
    FINAL_READS=("$output_file_1" "$output_file_2")
  else
    log_info "Running czid-dedup on reads..."
    # single file
    input_file="${FINAL_READS[0]}"
    filename=$(basename "$input_file")
    output_file="$TEMP_DIR/${filename%.fastq}_dedup.fastq"
    DEDUP_CMD=(czid-dedup -i "$input_file" -o "$output_file")
    if [[ -n "$DEDUP_L" ]]; then
      DEDUP_CMD+=(-l "$DEDUP_L")
    fi
    log_info "Executing czid-dedup command: ${DEDUP_CMD[*]}"
    "${DEDUP_CMD[@]}"
    FINAL_READS=("$output_file")
  fi
else
  log_info "No deduplication requested. Skipping step 2.3."
fi

############################################
# (Optional) Downsampling with rasusa
############################################
# The user may specify coverage+genome_size, or num_bases, or num_reads


if [[ "$DOWNSAMPLE" == true ]]; then
  log_info "========== Step 2.4: Downsampling with rasusa (optional) =========="

  log_info "Preparing input for rasusa..."

  # We'll generate new file(s):
  if [[ "$READ_TYPE" == "paired" ]]; then
    # we have two files in FINAL_READS
    input_file_1="${FINAL_READS[0]}"
    input_file_2="${FINAL_READS[1]}"
    filename_1=$(basename "$input_file_1")
    filename_2=$(basename "$input_file_2")
    output_file_1="$TEMP_DIR/${filename_1%.fastq}_ds.fastq"
    output_file_2="$TEMP_DIR/${filename_2%.fastq}_ds.fastq"
    # Decide which rasusa flags to use:
    RASUSA_CMD=(rasusa reads)
    if [[ $(echo "$COVERAGE > 0" | bc -l) -eq 1 && -n "$GENOME_SIZE" ]]; then
      RASUSA_CMD+=(--coverage "$COVERAGE" --genome-size "$GENOME_SIZE")
    elif [[ $NUM_BASES -gt 0 ]]; then
      RASUSA_CMD+=(--bases "$NUM_BASES")
    elif [[ $NUM_READS -gt 0 ]]; then
      RASUSA_CMD+=(--num "$NUM_READS")
    fi
    RASUSA_CMD+=(-o "$output_file_1" -o "$output_file_2" "$input_file_1" "$input_file_2")
    #Here show commands before execute
    log_info "Executing rasusa command: ${RASUSA_CMD[*]}"
    # run it
    "${RASUSA_CMD[@]}"
    FINAL_READS=("$output_file_1" "$output_file_2")
  else
    # single file
    input_file="${FINAL_READS[0]}"
    filename=$(basename "$input_file")
    output_file="$TEMP_DIR/${filename%.fastq}_ds.fastq"
    RASUSA_CMD=(rasusa reads)
    if [[ $(echo "$COVERAGE > 0" | bc -l) -eq 1 && -n "$GENOME_SIZE" ]]; then
      RASUSA_CMD+=(--coverage "$COVERAGE" --genome-size "$GENOME_SIZE")
    elif [[ $NUM_BASES -gt 0 ]]; then
      RASUSA_CMD+=(--bases "$NUM_BASES")
    elif [[ $NUM_READS -gt 0 ]]; then
      RASUSA_CMD+=(--num "$NUM_READS")
    fi
    RASUSA_CMD+=(-o "$output_file" "$input_file")
    log_info "Executing rasusa command: ${RASUSA_CMD[*]}"
    "${RASUSA_CMD[@]}"
    FINAL_READS=("$output_file")
  fi
else
  log_info "No downsampling requested. Skipping Step 2.4."
fi

###########################################
# Create summary table with read_statistics
###########################################
log_info "========== Step 2.5: Create summary table with read_statistics =========="

export -f analyze_fastq
export STATS_FILE
#Analyze fastq in parallel, keep the parallelizing?

input_file="${FINAL_READS[0]}"
filename=$(basename "$input_file")
sample_code=$(echo "${filename}" |sed -E 's/(_dedup|_ds|_dedup_ds)?\.fastq//')

if [[ "$READ_TYPE" == "paired" ]]; then
    input_file2="${FINAL_READS[1]}"
    filename2=$(basename "$input_file2")
    sample_code2=$(echo "${filename2}" | sed -E 's/(_dedup|_ds|_dedup_ds)?\.fastq//')
    log_info "Sample code 1: $sample_code"
    log_info "Sample code 2: $sample_code2"

    #log_info "What find is finding"
    #find "$TEMP_DIR" -type f \( -name "${sample_code}*.fastq" -o -name "${sample_code2}*.fastq" \) -print0 | sort -z

    sample_stats=""
    while IFS= read -r -d '' fastq_file; do
      log_info "Processing: $fastq_file"
      stats_line=$(analyze_fastq "$fastq_file")
      if [ -z "$sample_stats" ]; then
        sample_stats="$stats_line"
      else
        sample_stats+=$'\n'"$stats_line"
      fi
    done < <(find "$TEMP_DIR" -type f \( -name "${sample_code}*.fastq" -o -name "${sample_code2}*.fastq" \) -print0 | sort -z)
    
    # Write all the statistics for the sample at once
    if [[ -n "$sample_stats" ]]; then
      append_stats "$sample_stats"
    fi
else
    log_info "Sample code: $sample_code"
    #log_info "What find is finding"
    #find "$TEMP_DIR" -type f \( -name "${sample_code}*.fastq" \) -print0 | sort -z

    sample_stats=""
    while IFS= read -r -d '' fastq_file; do
      log_info "Processing: $fastq_file"
      stats_line=$(analyze_fastq "$fastq_file")
      if [ -z "$sample_stats" ]; then
        sample_stats="$stats_line"
      else
        sample_stats+=$'\n'"$stats_line"
      fi
    done < <(find "$TEMP_DIR" -type f -name "${sample_code}*.fastq" -print0 | sort -z)
    
    if [[ -n "$sample_stats" ]]; then
      append_stats "$sample_stats"
    fi
fi

log_info "Read statistics have been saved in file: $STATS_FILE"

############################################
# Finally, run read2tree
############################################
log_info "========== Step 2.6: Running Read2tree (step 2 map) =========="

log_info "Running read2tree --step 2map with read_type=$READ_TYPE, minimap2_options='$MINIMAP2_OPTIONS', threads=$THREADS"

# If we have paired => pass two files as arguments
# Otherwise pass single file
READ2TREE_CMD=( read2tree --step 2map
                --standalone_path marker_genes
                --dna_reference dna_ref.fa
                --threads "$THREADS"
                --output_path "$OUT_DIR"
                --debug )

READ2TREE_CMD+=( --read_type "$MINIMAP2_OPTIONS" )
if [[ "$READ_TYPE" == "paired" ]]; then
  READ2TREE_CMD+=( --reads "${FINAL_READS[0]}" "${FINAL_READS[1]}" )
else
  READ2TREE_CMD+=( --reads "${FINAL_READS[0]}" )
fi

log_info "Executing read2tree command: ${READ2TREE_CMD[*]}"
"${READ2TREE_CMD[@]}"

log_info "read2tree step 2 map completed successfully."
# No longer needed to remove the lock file, since we are using the same stats file for all samples.
# rm -f "${STATS_FILE}.lock"


# Optionally do step 3combine
# log_info "Running read2tree --step 3combine..."
# read2tree --step 3combine \
#           --standalone_path marker_genes \
#           --dna_reference dna_ref.fa \
#           --output_path read2tree_output \
#           --tree \
#           --debug

#log_info "Done. If you want to run iqtree or similar, you can do so now."
