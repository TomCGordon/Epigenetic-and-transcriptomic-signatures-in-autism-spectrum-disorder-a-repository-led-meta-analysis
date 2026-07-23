#!/usr/bin/env python3
"""
Download and stage public repository inputs for GSE81541.
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

DOWNLOADS = [{'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156966/suppl/GSM2156966_JLKD002.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156966.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 128901381, 'expected_sha256': '1622c8d12db3f92f25b7fb1ee458d455bbf2de77b30b5076852b7165ebd1eff1'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156967/suppl/GSM2156967_JLKD003.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156967.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 130362321, 'expected_sha256': '894abfccb41c910448fcd169a776db4567c76e712512299bff940ba974c5be08'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156968/suppl/GSM2156968_JLKD005.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156968.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 143713191, 'expected_sha256': '9e457a35139575d12b58f8208e603adb42581a4367d442cbcb9c5e8dc49fd9d8'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156969/suppl/GSM2156969_JLKD001.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156969.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 127132779, 'expected_sha256': '1a5877598d07e3900ff27254c858ae5530255d659d1d178919c02874474637c3'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156970/suppl/GSM2156970_JLKD014.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156970.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 167672527, 'expected_sha256': 'a058eed930cdfeda2a0d565759225005f07396c2b2646160bebb6b1fcaea2fc9'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156971/suppl/GSM2156971_JLKD004.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156971.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 136040105, 'expected_sha256': '491c47460e7b668e2917b8616da994e9133533ccae599571e897af934d0f5894'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156972/suppl/GSM2156972_JLKD040.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156972.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 182705486, 'expected_sha256': 'e4b13aa2b3088d870346f519693d368422a83a5fad1266d5785716d0296d251e'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156973/suppl/GSM2156973_JLKD041.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156973.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 182715004, 'expected_sha256': 'a51b31c5d7a34bf546f5d4af66c19b5cd38ed6d7cbb4ee1e9361623c2d73789e'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156974/suppl/GSM2156974_JLKD042.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156974.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 184013401, 'expected_sha256': 'bc62e3fc610ee5dab8f623328000d844f8cfe9ffd60a058a410adae282cb2afd'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156975/suppl/GSM2156975_JLKD026.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156975.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 158416633, 'expected_sha256': 'cbf02086fdfd946e1de3aa36ce36716d8863a4194bccd5c0a727d4ab6485aaf1'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156976/suppl/GSM2156976_JLKD028.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156976.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 144597767, 'expected_sha256': 'b13f62fad2f4d9cb757c52661286a71ffa08b50c7c4581f9aa4e922fdbfd1f10'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156998/suppl/GSM2156998_JLKD009.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156998.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 175755501, 'expected_sha256': '1c6935cdd3caa6681139ea0dc81f0ba5e3057ba76dc92667843eee3c6362bd35'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2156nnn/GSM2156999/suppl/GSM2156999_JLKD013.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156999.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 172212170, 'expected_sha256': '04fc3ca24b3d6c4aa0c95d3a3accacffa0801ee3da802af2b971c2329513ff43'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2157nnn/GSM2157000/suppl/GSM2157000_JLKD016.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157000.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 172110276, 'expected_sha256': 'b7237f7abb9c985337e79a873ff9b85702b26ea1ce3c6dbd738d3d572bfe35cc'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2157nnn/GSM2157001/suppl/GSM2157001_JLKD018.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157001.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 180508054, 'expected_sha256': '41f76a0daf88c9c320725e26d64b0cac4907def21d08bb062b078b6f605e2991'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2157nnn/GSM2157002/suppl/GSM2157002_JLKD019.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157002.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 175676620, 'expected_sha256': '899bd1a29f040fe9a3ea6c3c640ec0c182152bcc352c0cd475f9f162bfa53dce'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2157nnn/GSM2157003/suppl/GSM2157003_JLKD020.bed.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157003.bed.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 177292075, 'expected_sha256': 'ff5c40e09e2c1758e10c301139618ace40cfcc1cd0a537f8c6bf7bfb04920e68'}]

EXPECTED_FILES = [{'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156966.bed.gz', 'expected_bytes': 128901381, 'expected_sha256': '1622c8d12db3f92f25b7fb1ee458d455bbf2de77b30b5076852b7165ebd1eff1'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156967.bed.gz', 'expected_bytes': 130362321, 'expected_sha256': '894abfccb41c910448fcd169a776db4567c76e712512299bff940ba974c5be08'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156968.bed.gz', 'expected_bytes': 143713191, 'expected_sha256': '9e457a35139575d12b58f8208e603adb42581a4367d442cbcb9c5e8dc49fd9d8'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156969.bed.gz', 'expected_bytes': 127132779, 'expected_sha256': '1a5877598d07e3900ff27254c858ae5530255d659d1d178919c02874474637c3'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156970.bed.gz', 'expected_bytes': 167672527, 'expected_sha256': 'a058eed930cdfeda2a0d565759225005f07396c2b2646160bebb6b1fcaea2fc9'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156971.bed.gz', 'expected_bytes': 136040105, 'expected_sha256': '491c47460e7b668e2917b8616da994e9133533ccae599571e897af934d0f5894'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156972.bed.gz', 'expected_bytes': 182705486, 'expected_sha256': 'e4b13aa2b3088d870346f519693d368422a83a5fad1266d5785716d0296d251e'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156973.bed.gz', 'expected_bytes': 182715004, 'expected_sha256': 'a51b31c5d7a34bf546f5d4af66c19b5cd38ed6d7cbb4ee1e9361623c2d73789e'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156974.bed.gz', 'expected_bytes': 184013401, 'expected_sha256': 'bc62e3fc610ee5dab8f623328000d844f8cfe9ffd60a058a410adae282cb2afd'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156975.bed.gz', 'expected_bytes': 158416633, 'expected_sha256': 'cbf02086fdfd946e1de3aa36ce36716d8863a4194bccd5c0a727d4ab6485aaf1'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156976.bed.gz', 'expected_bytes': 144597767, 'expected_sha256': 'b13f62fad2f4d9cb757c52661286a71ffa08b50c7c4581f9aa4e922fdbfd1f10'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156998.bed.gz', 'expected_bytes': 175755501, 'expected_sha256': '1c6935cdd3caa6681139ea0dc81f0ba5e3057ba76dc92667843eee3c6362bd35'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2156999.bed.gz', 'expected_bytes': 172212170, 'expected_sha256': '04fc3ca24b3d6c4aa0c95d3a3accacffa0801ee3da802af2b971c2329513ff43'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157000.bed.gz', 'expected_bytes': 172110276, 'expected_sha256': 'b7237f7abb9c985337e79a873ff9b85702b26ea1ce3c6dbd738d3d572bfe35cc'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157001.bed.gz', 'expected_bytes': 180508054, 'expected_sha256': '41f76a0daf88c9c320725e26d64b0cac4907def21d08bb062b078b6f605e2991'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157002.bed.gz', 'expected_bytes': 175676620, 'expected_sha256': '899bd1a29f040fe9a3ea6c3c640ec0c182152bcc352c0cd475f9f162bfa53dce'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE81541/WGBS/GSE81541/GSM2157003.bed.gz', 'expected_bytes': 177292075, 'expected_sha256': 'ff5c40e09e2c1758e10c301139618ace40cfcc1cd0a537f8c6bf7bfb04920e68'}]

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
    parser = argparse.ArgumentParser(description="Download repository inputs for GSE81541.")
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
            'accession': 'GSE81541',
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
        pd.DataFrame(log_rows).to_csv(os.path.join(log_dir, "GSE81541_download_log.csv"), index=False)
        checks_df = pd.DataFrame(checks)
        checks_df.to_csv(os.path.join(log_dir, "GSE81541_file_check.csv"), index=False)
        
        failures = checks_df[~checks_df['present'] | ~checks_df['size_match']]
        if not failures.empty:
            print(f"WARNING: {len(failures)} expected file(s) failed verification.")
        else:
            print(f"Dataset GSE81541 staged and verified successfully.")

if __name__ == "__main__":
    main()
