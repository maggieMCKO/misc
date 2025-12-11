#!/usr/bin/env bash
set -euo pipefail

# --- CONFIGURATION ---
SOURCE_DIR="$HOME"
BACKUP_DIR="/scratch/users/$USER/backup_staging"
VERIFY_LOG="$BACKUP_DIR/backup_verify.log"
ERROR_LOG="$BACKUP_DIR/backup_errors.log"

echo "Recreating corrupted archives based on: $VERIFY_LOG"
cd "$SOURCE_DIR" || exit 1

# Ensure logs exist
touch "$ERROR_LOG"

# Grep lines with 'Archive CORRUPT:' and extract the archive name
# Example line: "Archive CORRUPT: AMY.tar.gz"
grep 'Archive CORRUPT:' "$VERIFY_LOG" | while read -r line; do
    archive_name=$(echo "$line" | awk '{print $3}')
    # Safety: skip if empty
    [[ -n "${archive_name:-}" ]] || continue

    # Strip .tar.gz to get directory name
    dir_name="${archive_name%.tar.gz}"

    echo "------------------------------------------------"
    echo "Recreating archive for: $dir_name ($archive_name)"

    # Check that source directory exists
    if [[ ! -d "$dir_name" ]]; then
        echo "WARNING: Source directory not found: $dir_name" | tee -a "$ERROR_LOG"
        continue
    fi

    dest_path="$BACKUP_DIR/$archive_name"
    md5_path="$dest_path.md5"

    # Remove old (corrupt) archive and checksum if present
    rm -f "$dest_path" "$md5_path"

    # Recreate archive
    echo "Compressing $dir_name -> $dest_path"
    if tar -I 'pigz -9' -cf "$dest_path" "$dir_name" 2>>"$ERROR_LOG"; then
        echo "Success: $archive_name"

        # Verify archive structure
        if tar -tzf "$dest_path" > /dev/null 2>>"$ERROR_LOG"; then
            echo "Archive OK after recreate: $archive_name"
            # Generate checksum
            (cd "$BACKUP_DIR" && md5sum "$archive_name" > "$md5_path")
        else
            echo "Archive still CORRUPT after recreate: $archive_name" | tee -a "$ERROR_LOG"
        fi
    else
        echo "ERROR: Failed to recompress $dir_name" | tee -a "$ERROR_LOG"
    fi
done

echo "Done recreating corrupted archives."
