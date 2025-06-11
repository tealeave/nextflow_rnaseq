#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

def harmonize_filenames(data_dir: Path, dry_run: bool):
    """
    Finds and renames inconsistently named FASTQ files in a data directory.
    This handles two cases:
    1. Renaming raw sequencing core files (e.g., 4R044-...)
    2. Standardizing the .fq.gz extension to .fastq.gz
    """
    if not data_dir.is_dir():
        print(f"Error: Directory not found at '{data_dir}'", file=sys.stderr)
        sys.exit(1)

    print(">>> Starting filename harmonization...")
    if dry_run:
        print(">>> Mode: DRY RUN (no files will be changed)")
    else:
        print(">>> Mode: COMMIT (files will be renamed)")
    print("-" * 50)

    processed_count = 0
    # Iterate through each item in the data directory
    for sample_dir in sorted(data_dir.iterdir()):
        if not sample_dir.is_dir():
            continue

        # --- CASE 1: Handle raw sequencing core files ---
        read1_files = list(sample_dir.glob("4R044-L6-*.txt.gz"))
        read2_files = list(sample_dir.glob("4R044-L7-*.txt.gz"))

        if len(read1_files) == 1 and len(read2_files) == 1:
            processed_count += 1
            old_r1_path = read1_files[0]
            old_r2_path = read2_files[0]
            sample_name = sample_dir.name
            new_r1_path = sample_dir / f"{sample_name}_1.fastq.gz"
            new_r2_path = sample_dir / f"{sample_name}_2.fastq.gz"

            print(f"Found '4R044' files in: '{sample_name}'")
            if dry_run:
                print(f"  [DRY RUN] Would rename '{old_r1_path.name}'\n    -> to '{new_r1_path.name}'")
                print(f"  [DRY RUN] Would rename '{old_r2_path.name}'\n    -> to '{new_r2_path.name}'")
            else:
                try:
                    old_r1_path.rename(new_r1_path)
                    old_r2_path.rename(new_r2_path)
                    print(f"  SUCCESS: Renamed files for '{sample_name}'")
                except OSError as e:
                    print(f"  ERROR: Could not rename files in '{sample_name}'. Reason: {e}", file=sys.stderr)
            print("-" * 50)
            continue # Move to the next directory

        # --- CASE 2: Handle files with .fq.gz extension ---
        fq_files = list(sample_dir.glob("*.fq.gz"))
        
        if fq_files:
            processed_count += 1
            print(f"Found '.fq.gz' files in: '{sample_dir.name}'")
            for old_path in fq_files:
                # Replace the extension safely
                new_name = old_path.name.replace(".fq.gz", ".fastq.gz")
                new_path = old_path.with_name(new_name)

                if dry_run:
                    print(f"  [DRY RUN] Would rename '{old_path.name}' to '{new_name}'")
                else:
                    try:
                        old_path.rename(new_path)
                        print(f"  SUCCESS: Renamed '{old_path.name}' to '{new_name}'")
                    except OSError as e:
                        print(f"  ERROR: Could not rename '{old_path.name}'. Reason: {e}", file=sys.stderr)
            print("-" * 50)


    if processed_count == 0:
        print("No files matching any renaming patterns were found.")
    else:
        print(f"Scan complete. Found {processed_count} sample directories to process.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Harmonize FASTQ filenames from a sequencing core to a standard format.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "data_dir",
        type=Path,
        help="The path to the top-level data directory containing sample subfolders."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the rename operations that would be performed without actually changing any files."
    )
    args = parser.parse_args()
    harmonize_filenames(data_dir=args.data_dir, dry_run=args.dry_run)