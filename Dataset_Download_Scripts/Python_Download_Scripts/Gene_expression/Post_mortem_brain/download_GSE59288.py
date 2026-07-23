#!/usr/bin/env python3
"""
Download and stage public repository inputs for GSE59288.
Converted from R download script to standalone Python script.
"""

import os
import sys
import time
import shutil
import urllib.request
import hashlib
import tarfile
import argparse
import pandas as pd
from datetime import datetime, timezone

DOWNLOADS = [{'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59288/suppl/GSE59288_exp_mRNA.txt.gz', 'destination': '01_Raw_Data/Gene_expression/Post_mortem_brain/GSE59288/GSE59288_exp_mRNA.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 2120725, 'expected_sha256': '3d7ab94695cc09609a19826efad46eab214370ec7da760ff28e23a545c681a21'}]

EXPECTED_FILES = [{'path': '01_Raw_Data/Gene_expression/Post_mortem_brain/GSE59288/GSE59288_exp_mRNA.txt.gz', 'expected_bytes': 2120725, 'expected_sha256': '3d7ab94695cc09609a19826efad46eab214370ec7da760ff28e23a545c681a21'}]

RESOURCES = {}

def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest().lower()

def download_one(url, destination, overwrite=False, dry_run=False, retries=3, timeout=7200):
    if dry_run:
        print(f"[dry-run] {url} -> {destination}")
        return "dry_run"
    
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    if os.path.exists(destination) and os.path.getsize(destination) > 0 and not overwrite:
        print(f"[present] {destination}")
        return "already_present"
    
    part_file = destination + ".part"
    if os.path.exists(part_file):
        os.remove(part_file)
        
    last_error = None
    for attempt in range(1, retries + 1):
        print(f"[download {attempt}/{retries}] {url}")
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Python-Dataset-Downloader/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(part_file, 'wb') as out_f:
                shutil.copyfileobj(resp, out_f)
            if os.path.exists(part_file) and os.path.getsize(part_file) > 0:
                if os.path.exists(destination):
                    os.remove(destination)
                os.rename(part_file, destination)
                return "downloaded"
        except Exception as e:
            last_error = str(e)
            if os.path.exists(part_file):
                os.remove(part_file)
            if attempt < retries:
                time.sleep(min(30, 2 ** attempt))
                
    raise RuntimeError(f"Download failed for {url}: {last_error}")

def main():
    parser = argparse.ArgumentParser(description="Download repository inputs for GSE59288.")
    parser.add_argument("--output-root", default=os.path.join(os.getcwd(), "downloaded_public_inputs"),
                        help="Root directory for downloaded inputs")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files")
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without downloading")
    parser.add_argument("--verify-sha256", action="store_true", help="Verify SHA256 checksums")
    parser.add_argument("--retries", type=int, default=3, help="Number of download retries")
    parser.add_argument("--timeout", type=int, default=7200, help="Download timeout in seconds")
    
    args = parser.parse_args()
    
    items_to_download = DOWNLOADS or RESOURCES
    if not items_to_download:
        print("No download targets defined.")
        return

    log_rows = []
    for row in items_to_download:
        dest_rel = os.path.normpath(row['destination'])
        dest_full = os.path.join(args.output_root, dest_rel)
        
        status = download_one(
            url=row['remote_url'],
            destination=dest_full,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            retries=args.retries,
            timeout=args.timeout
        )
        
        action = row.get('action', 'none')
        ext_root_rel = row.get('extraction_root')
        if status in ("downloaded", "already_present") and action and action != "none" and ext_root_rel:
            ext_root_full = os.path.join(args.output_root, os.path.normpath(ext_root_rel))
            os.makedirs(ext_root_full, exist_ok=True)
            if not args.dry_run and tarfile.is_tarfile(dest_full):
                print(f"[extract] Extracting {dest_full} to {ext_root_full}...")
                with tarfile.open(dest_full, 'r:*') as tar:
                    tar.extractall(path=ext_root_full)
                    
        log_rows.append({
            'accession': 'GSE59288',
            'remote_url': row['remote_url'],
            'destination': dest_full,
            'status': status,
            'checked_at_utc': datetime.now(timezone.utc).isoformat()
        })
        
    if not args.dry_run and EXPECTED_FILES:
        checks = []
        for row in EXPECTED_FILES:
            file_rel = os.path.normpath(row['path'])
            file_full = os.path.join(args.output_root, file_rel)
            present = os.path.exists(file_full)
            obs_bytes = os.path.getsize(file_full) if present else None
            exp_bytes = row.get('expected_bytes')
            size_match = present and (exp_bytes is None or obs_bytes == exp_bytes)
            
            exp_sha256 = row.get('expected_sha256')
            obs_sha256 = compute_sha256(file_full) if (present and args.verify_sha256) else None
            sha256_match = (obs_sha256 == exp_sha256.lower()) if (present and args.verify_sha256 and exp_sha256) else None
            
            checks.append({
                'file': file_full,
                'present': present,
                'expected_bytes': exp_bytes,
                'observed_bytes': obs_bytes,
                'size_match': size_match,
                'expected_sha256': exp_sha256,
                'observed_sha256': obs_sha256,
                'sha256_match': sha256_match
            })
            
        log_dir = os.path.join(args.output_root, "04_Download_Logs")
        os.makedirs(log_dir, exist_ok=True)
        pd.DataFrame(log_rows).to_csv(os.path.join(log_dir, "GSE59288_download_log.csv"), index=False)
        checks_df = pd.DataFrame(checks)
        checks_df.to_csv(os.path.join(log_dir, "GSE59288_file_check.csv"), index=False)
        
        failures = checks_df[~checks_df['present'] | ~checks_df['size_match']]
        if not failures.empty:
            print(f"WARNING: {len(failures)} expected file(s) failed verification.")
        else:
            print(f"Dataset GSE59288 staged and verified successfully.")

if __name__ == "__main__":
    main()
