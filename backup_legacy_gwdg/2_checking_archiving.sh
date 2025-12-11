#!/usr/bin/env bash
set -euo pipefail

# --- CONFIGURATION ---
BACKUP_DIR="/scratch/users/$USER/backup_staging"

LOG_FILE="$BACKUP_DIR/backup_verify.log"
ERROR_LOG="$BACKUP_DIR/backup_errors.log"

echo "Starting archive verification in: $BACKUP_DIR"
echo "Log:   $LOG_FILE"
echo "Error: $ERROR_LOG"

# Empty or create logs
: > "$LOG_FILE"
touch "$ERROR_LOG"

cd "$BACKUP_DIR" || exit 1

# Find all .tar.gz archives (non-recursive; add -maxdepth if you prefer)
for archive_path in "$BACKUP_DIR"/*.tar.gz; do
    # Handle case where glob doesn't match anything
    [[ -e "$archive_path" ]] || { echo "No .tar.gz archives found."; break; }

    archive_name="$(basename "$archive_path")"
    echo "Checking: $archive_name" | tee -a "$LOG_FILE"

    if tar -tzf "$archive_path" > /dev/null 2>>"$ERROR_LOG"; then
        echo "Archive OK: $archive_name" | tee -a "$LOG_FILE"
    else
        echo "Archive CORRUPT: $archive_name" | tee -a "$LOG_FILE"
        echo "Archive CORRUPT: $archive_name" >>"$ERROR_LOG"
    fi
done

echo "Verification finished."
