#!/usr/bin/env bash
set -euo pipefail
# usage: ./uninstall.sh [/path/to/prefix]
REPO_ROOT="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAVED_PREFIX=$(cat "$REPO_ROOT/.install_prefix" 2>/dev/null || echo "")
PREFIX="${SAVED_PREFIX:-$HOME/.local/bin}"

if [[ $# -gt 0 && "$1" != -* ]]; then
  PREFIX="$1"
  shift
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      echo "Usage: ./uninstall.sh [/path/to/prefix]"
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown option: $1"
      exit 1
      ;;
  esac
  shift
done

SYMLINKS=("o2t-step1" "o2t-step2" "o2t-sra")

echo "[INFO] Removing symlinks from: $PREFIX"
for name in "${SYMLINKS[@]}"; do
  target="$PREFIX/$name"
  if [[ -L "$target" ]]; then
    rm -f "$target"
    echo "[OK] Removed $target"
  elif [[ -e "$target" ]]; then
    echo "[WARN] $target exists but is not a symlink. Skipping."
  else
    echo "[INFO] $target not found."
  fi
done

removed_any=false
comp_files=("omni2tree")
comp_symlinks=("o2t-step1" "o2t-step2" "o2t-sra" "v2t-step1" "v2t-step2" "v2t-sra")

for base_path in "$HOME/.local/share/bash-completion/completions" "/usr/share/bash-completion/completions"; do
  for comp_file in "${comp_files[@]}"; do
    if [[ -e "$base_path/$comp_file" ]]; then
      if rm -f "$base_path/$comp_file" 2>/dev/null; then
        echo "[OK] Removed $base_path/$comp_file"
        removed_any=true
      else
        echo "[WARN] Could not remove $base_path/$comp_file (permission denied?)."
      fi
    fi
  done
  
  for link in "${comp_symlinks[@]}"; do
    if [[ -L "$base_path/$link" ]]; then
      if rm -f "$base_path/$link" 2>/dev/null; then
        echo "[OK] Removed completion symlink $base_path/$link"
        removed_any=true
      else
        echo "[WARN] Could not remove completion symlink $base_path/$link (permission denied?)."
      fi
    fi
  done
done

echo "[DONE] Omni2tree uninstall complete."
