#!/usr/bin/env bash
set -euo pipefail
# usage: ./install.sh [/path/to/prefix]
PREFIX="${1:-/usr/local/bin}"  #prefix is the place where symlinks will be created 
APP_NAME="omni2tree"
BASENAMES=("omni2tree_step1.sh:o2t-step1"
           "omni2tree_step2.sh:o2t-step2"
           "omni2tree_sra.sh:o2t-sra")

# Resolve install.sh absolute directory, following symlinks
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
REPO_ROOT="$(cd -P "$(dirname "$SOURCE")" && pwd)"
SCRIPTS_DIR="$REPO_ROOT/scripts" #here will be the scripts to be installed

# Check the three targets exist and are executable
for pair in "${BASENAMES[@]}"; do
  src="${pair%%:*}"
  target="$SCRIPTS_DIR/$src"

  # Must exist
  [[ -f "$target" ]] || { echo "[ERROR] Missing file: $target"; exit 1; }

  # Ensure executable bit; try to fix if not
  if [[ ! -x "$target" ]]; then
    if [[ -w "$target" ]]; then
      chmod u+x "$target" || true
    fi
  fi

  # Final check
  [[ -x "$target" ]] || {
    echo "[ERROR] Not executable and could not fix perms: $target"
    echo "        Try: chmod +x \"$target\""
    exit 1
  }
done

# Create the prefix if it's a user-owned location
if ! mkdir -p "$PREFIX" 2>/dev/null; then
  echo "[ERROR] Can't create $PREFIX (need permissions?)."
  echo "        Try: sudo ./install.sh   or   ./install.sh \"\$HOME/.local/bin\""
  exit 1
fi
if [[ ! -w "$PREFIX" ]]; then
  echo "[ERROR] $PREFIX is not writable."
  echo "        Try: sudo ./install.sh   or   ./install.sh \"\$HOME/.local/bin\""
  exit 1
fi

# Create absolute symlinks
for pair in "${BASENAMES[@]}"; do
  src="${pair%%:*}"
  dst="${pair##*:}"
  # refuse to clobber a regular file
  if [[ -e "$PREFIX/$dst" && ! -L "$PREFIX/$dst" ]]; then
    echo "[ERROR] $PREFIX/$dst exists and is not a symlink. Remove it first."
    exit 1
  fi
  ln -sfn "$SCRIPTS_DIR/$src" "$PREFIX/$dst"
  echo "[OK] $PREFIX/$dst -> $SCRIPTS_DIR/$src"
done

# Install bash completions
COMP_SRC="$REPO_ROOT/share/bash-completion/completions/$APP_NAME"
COMP_INSTALLED=false
# install bash completions 
if [[ -f "$COMP_SRC" ]]; then
  if [[ "${EUID:-$(id -u)}" -eq 0 ]]; then
    COMP_SHARE="/usr/share/bash-completion/completions"
  else
    COMP_SHARE="$HOME/.local/share/bash-completion/completions"
  fi
  
  mkdir -p "$COMP_SHARE" 2>/dev/null || true
  
  if [[ -d "$COMP_SHARE" && -w "$COMP_SHARE" ]]; then
    if cp -f "$COMP_SRC" "$COMP_SHARE/$APP_NAME" 2>/dev/null; then
      echo "[OK] Bash completions installed to $COMP_SHARE/$APP_NAME"
      COMP_INSTALLED=true
      
      # Symlinks so bash-completion can reuse the same completion file per command
      for cmd in o2t-step1 o2t-step2 o2t-sra; do
        ln -sf "$APP_NAME" "$COMP_SHARE/$cmd" 2>/dev/null || true
      done
      echo "[OK] Created completion symlinks: o2t-step1, o2t-step2, o2t-sra"
    fi
  fi
fi

# Warn if the prefix is not in PATH
case ":$PATH:" in
  *":$PREFIX:"*) echo "[OK] $PREFIX is in PATH.";;
  *) echo "[WARN] $PREFIX is not in PATH. Add this to your shell rc:"
     echo "export PATH=\"$PREFIX:\$PATH\"";;
esac

if ! (printf '%s\n' "$PREFIX" > "$REPO_ROOT/.install_prefix") 2>/dev/null; then
  echo "[WARN] Could not write $REPO_ROOT/.install_prefix (permission denied?)."
  echo "[WARN] Uninstall may not auto-detect the prefix; pass it explicitly to ./uninstall.sh."
fi

echo "[DONE] Installed Omni2tree executables: o2t-step1, o2t-step2, o2t-sra"
if [[ "$COMP_INSTALLED" == true ]]; then
  echo "[INFO] To enable completions now, run:"
  echo "  source ~/.bashrc"
fi
echo "[INFO] To uninstall, run: ./uninstall.sh \"$PREFIX\""
