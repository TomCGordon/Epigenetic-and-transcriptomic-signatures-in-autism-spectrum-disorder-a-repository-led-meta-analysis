#!/usr/bin/env python3
"""
Download and stage public repository inputs for GSE77103.
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

DOWNLOADS = [{'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77103/suppl/GSE77103_RAW.tar', 'destination': '01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW.tar', 'action': 'extract_archive', 'extraction_root': '01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW', 'expected_bytes': 25323520, 'expected_sha256': 'f6d148b348bcc80cdb8991778106c395ae1563e700d5032cf2576f9d6a8a10e2'}]

EXPECTED_FILES = [{'path': '01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW.tar', 'expected_bytes': 25323520, 'expected_sha256': 'f6d148b348bcc80cdb8991778106c395ae1563e700d5032cf2576f9d6a8a10e2'}, {'path': '01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044363_C2.txt.gz', 'expected_bytes': 3194167, 'expected_sha256': '9dda3d40b82dc9e6ae06cc2b89c1355bc93bb3d8fb79b576ced6cf1cb20f0519'}, {'path': '01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044364_C4.txt.gz', 'expected_bytes': 3190286, 'expected_sha256': '5135754a7945252e87883da7875aa555cb40fcfaf0eead82078f166cdfb4dce7'}, {'path': '01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044365_C5.txt.gz', 'expected_bytes': 3177295, 'expected_sha256': 'd63daa354ed8c6aa95e267b92859760f3ea7d0a42d9efccc0bf9d3ff486aa240'}, {'path': '01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044366_C6.txt.gz', 'expected_bytes': 3167885, 'expected_sha256': '17351fb88706445a672d00d19db1570d5852632153f2d0ba9cd76745a78c3089'}, {'path': '01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044367_A2.txt.gz', 'expected_bytes': 3142461, 'expected_sha256': 'c574e1f72d107a171c458eaa0749e2d1447d85dd21f27f81562590530d98fee3'}, {'path': '01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044368_A3.txt.gz', 'expected_bytes': 3133875, 'expected_sha256': 'e0ca8536ec8e8810f2417d87137c8b36ca1eedb3cfc191be1c39246dacbbfa25'}, {'path': '01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044369_A4.txt.gz', 'expected_bytes': 3154385, 'expected_sha256': '9649a8042fe10b3a5f14c15332760bcd72039a4db66a286df4f45d169506662f'}, {'path': '01_Raw_Data/Gene_expression/Blood_family/GSE77103/GSE77103_RAW/GSM2044370_A5.txt.gz', 'expected_bytes': 3147378, 'expected_sha256': 'daecc71df818ef0aacf7a24ecebeec0125e16990a1ded26e11fc201c3c372b23'}]

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
    parser = argparse.ArgumentParser(description="Download repository inputs for GSE77103.")
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
            'accession': 'GSE77103',
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
        pd.DataFrame(log_rows).to_csv(os.path.join(log_dir, "GSE77103_download_log.csv"), index=False)
        checks_df = pd.DataFrame(checks)
        checks_df.to_csv(os.path.join(log_dir, "GSE77103_file_check.csv"), index=False)
        
        failures = checks_df[~checks_df['present'] | ~checks_df['size_match']]
        if not failures.empty:
            print(f"WARNING: {len(failures)} expected file(s) failed verification.")
        else:
            print(f"Dataset GSE77103 staged and verified successfully.")

if __name__ == "__main__":
    main()
