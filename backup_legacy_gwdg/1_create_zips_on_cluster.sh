#!/usr/bin/env bash
set -euo pipefail

# --- CONFIGURATION ---
SOURCE_DIR="$HOME"                        # Where your data is now
BACKUP_DIR="/scratch/users/$USER/backup_staging" # Fast storage for the zips
EXCLUDE_ITEMS=("scratch" "tmp" "core" "Anaconda3" "miniconda3") # Items to skip

# --- SETUP ---
mkdir -p "$BACKUP_DIR"
# Create a sub-folder for loose files so they don't get mixed with the zips
LOOSE_FILES_DIR="$BACKUP_DIR/loose_files"
mkdir -p "$LOOSE_FILES_DIR"

echo "Starting Smart Backup Job..."
cd "$SOURCE_DIR" || exit 1

# Loop through everything in Home
for item in *; do
    # Exclusions
    if [[ " ${EXCLUDE_ITEMS[*]} " =~ " ${item} " ]]; then
        echo "SKIPPING (Excluded): $item"
        continue
    fi

    # DIRECTORIES → tar.gz + md5
    if [[ -d "$item" ]]; then
        archive_name="${item}.tar.gz"
        dest_path="$BACKUP_DIR/$archive_name"
        md5_path="$dest_path.md5"

        # Skip if archive already exists
        if [[ -f "$dest_path" ]]; then
            echo "SKIPPING (Archive exists): $archive_name"
            # Optionally regenerate checksum if missing
            if [[ ! -f "$md5_path" ]]; then
                echo "Generating missing checksum for: $archive_name"
                (cd "$BACKUP_DIR" && md5sum "$archive_name" > "$md5_path")
            fi
            continue
        fi

        echo "[DIR] Compressing: $item"
        if tar -I 'pigz -9' -cf "$dest_path" "$item" 2>> "$BACKUP_DIR/backup_errors.log"; then
            echo "Success: $archive_name"
            echo "Generating checksum..."
            (cd "$BACKUP_DIR" && md5sum "$archive_name" > "$md5_path")
        else
            echo "ERROR: Failed to compress $item"
        fi

    # FILES → copy as-is
    elif [[ -f "$item" ]]; then
        echo "[FILE] Copying: $item"
        cp -p "$item" "$LOOSE_FILES_DIR/"

    else
        echo "WARNING: $item is not regular file/dir, skipping."
    fi
done

echo "------------------------------------------------"
echo "Backup staging complete."
echo " - Archives + checksums: $BACKUP_DIR"
echo " - Loose files:          $LOOSE_FILES_DIR"