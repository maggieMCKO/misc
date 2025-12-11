# CLEAN UP

Remove unimportant and unnecessary files before continuing.

# STEP 1: create archives (tar.gz)

After login in the terminal, make a new screen, which will keep the command running and does not require the terminal to stay open constantly:

```bash
screen -S NAME # to create one
```

In the created screen, navigate to `$HOME` and run `1_create_zips_on_cluster.sh`
This script will create a tar.gz for each folder in the /scratch

```bash
cd $HOME
bash 1_create_zips_on_cluster.sh
```

Depending on the amount of data, this will take some time. You can detach the screen (while keeps the script running).
To detach a screen, use the following keys

```bash
Ctrl + A + D (careful, not Ctrl + D, which kills the screen, then you need to start over)
```

You can check the progress from time to time

```bash
# To resume a detached screen:
screen -r NAME
```


# STEP2: downloading to hard disk

# STEP3: check whether the archives are properly generated

in screen because this also takes time.

```bash
bash 2_checking_archiving.sh
```

see the log generated `backup_verify.log`
if all archives are `OK`
if there are some `CORRUPT`, rerun for those folders using

```bash
bash 3_repeat_for_corrupt.sh
```

These scripts help you archive and back up data on the GWDG legacy system.

# CLEAN UP

Remove unimportant or unnecessary files before continuing.

# STEP 1: Create archives (tar.gz)

After logging in on the terminal, start a new screen session so the command keeps running without requiring the terminal to stay open:

```bash
screen -S NAME   # create a new screen session
```

In the new screen, go to `$HOME` and run `1_create_zips_on_cluster.sh`.
This script creates a `.tar.gz` archive for each folder under your home directory and writes them to `/scratch/users/$USER/backup_staging` (or the path you configured).

```bash
cd "$HOME"
bash 1_create_zips_on_cluster.sh
```

Depending on the amount of data, this can take a long time. You can detach from the screen session while the script continues to run.
To detach from a screen session:

```bash
Ctrl + A, then D   # careful: Ctrl + D alone will close the shell in the screen
```

To check progress later:

```bash
# Resume a detached screen session:
screen -r NAME
```


# STEP 2: Download to external hard disk

Use the download script (for example `2_download_backups_locally.sh`) from your local machine to copy the archives from the cluster to your external hard disk, typically via `rsync` over SSH.

# STEP 3: Check archives and fix corrupt ones

Run the archive-checking script in a screen session, because this can also take some time:

```bash
bash 2_checking_archiving.sh
```

Inspect the log file `backup_verify.log`.
If all archives are marked `Archive OK`, no further action is needed.
If some archives are marked `Archive CORRUPT`, rerun the archive creation only for those folders using:

```bash
bash 3_repeat_for_corrupt.sh
```

HOPEFULLY DONE.
