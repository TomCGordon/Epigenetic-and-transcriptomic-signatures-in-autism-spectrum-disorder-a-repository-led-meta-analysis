#!/usr/bin/env python3
"""
Download and stage public repository inputs for GSE109875.
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

DOWNLOADS = [{'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971944/suppl/GSM2971944_JLKD051_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971944.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 264284960, 'expected_sha256': '4a3c5c1401c9f3dea2326c053a81da508e1b52d15a8873395eb2e1746f945e07'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971945/suppl/GSM2971945_JLKD052_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971945.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 266852074, 'expected_sha256': '0d864da2d8827420adde8a4c32ec9df5a73fe66418284632639da7caea687cb8'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971946/suppl/GSM2971946_JLKD054_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971946.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 268162085, 'expected_sha256': '2749ced8526facd5993e07ad6b12f95bc9337a5d2e43c2b8df0e42bef39df1c6'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971947/suppl/GSM2971947_JLKD055_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971947.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 266937349, 'expected_sha256': 'eba1f625556de2a50938b6b5a21253ad833b58e1e23e3e8f2583c18f90a09453'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971948/suppl/GSM2971948_JLKD056_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971948.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 265581013, 'expected_sha256': '1f43203efc9353c901cebb4589e7585e2c99e19a5885c9ec3f9ddbbc10aab0ab'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971949/suppl/GSM2971949_JLKD057_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971949.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 262852523, 'expected_sha256': '8935a9a92cc18a6451ecabf275901d6071a3eb30a7218bd7c1db695c77b269fe'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971950/suppl/GSM2971950_JLKD058_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971950.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 262658188, 'expected_sha256': '81a5ffedce6706cee64f0bee99bc1f145042ea1c8fcf8755dbc5f9414e8bd2a2'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971951/suppl/GSM2971951_JLKD059_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971951.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 266345342, 'expected_sha256': '944dd645f98a19c1dbeca95de1740d1eb49b4c883562f91947e4320d246074d3'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971952/suppl/GSM2971952_JLKD060_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971952.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 262298260, 'expected_sha256': '735cb407a7c7be2b403555d9b972bac3048275ca89c8d42f48e9ceb82fae460b'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971953/suppl/GSM2971953_JLKD061_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971953.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 244560070, 'expected_sha256': '8c842fc5b49e7e88f301f7919206d09b7e6c5a3f4cc1f28ccdbb882a0e70b8e6'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971954/suppl/GSM2971954_JLKD062_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971954.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 268523047, 'expected_sha256': 'd5b4e564ea2aedee7488453e17171fb921fcbc0288c4b4ffcbddcb2e828aace7'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971955/suppl/GSM2971955_JLKD063_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971955.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 268286531, 'expected_sha256': '2d3443cea658f2a4c45e23c405bd8532cf58b75d59ec9dfc91c6975be8c366c1'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971956/suppl/GSM2971956_JLKD065_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971956.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 269230091, 'expected_sha256': '29c1ac2142e6c3b4a009fbd09d5293a1119a7f2c1ed097c0ce770d6d06d0c20d'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971957/suppl/GSM2971957_JLKD066_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971957.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 266725679, 'expected_sha256': 'dec0a66bf0b8cd9cb2d4d155cfff8755bbfb76b96cedad83039e2f2985433720'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971958/suppl/GSM2971958_JLKD067_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971958.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 268041654, 'expected_sha256': 'db3a481c398af7339e050c3b9bc2a712cbe4b0e2e007309286e1611b9bcd9967'}, {'remote_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2971nnn/GSM2971959/suppl/GSM2971959_JLKD069_filtered_trimmed_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz', 'destination': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971959.cpg_report.txt.gz', 'action': 'none', 'extraction_root': None, 'expected_bytes': 269663171, 'expected_sha256': 'a435a07f4773fa2b678b423b8e38853c217c03acfd5344bf1ce16a585d2ce4ac'}]

EXPECTED_FILES = [{'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971944.cpg_report.txt.gz', 'expected_bytes': 264284960, 'expected_sha256': '4a3c5c1401c9f3dea2326c053a81da508e1b52d15a8873395eb2e1746f945e07'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971945.cpg_report.txt.gz', 'expected_bytes': 266852074, 'expected_sha256': '0d864da2d8827420adde8a4c32ec9df5a73fe66418284632639da7caea687cb8'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971946.cpg_report.txt.gz', 'expected_bytes': 268162085, 'expected_sha256': '2749ced8526facd5993e07ad6b12f95bc9337a5d2e43c2b8df0e42bef39df1c6'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971947.cpg_report.txt.gz', 'expected_bytes': 266937349, 'expected_sha256': 'eba1f625556de2a50938b6b5a21253ad833b58e1e23e3e8f2583c18f90a09453'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971948.cpg_report.txt.gz', 'expected_bytes': 265581013, 'expected_sha256': '1f43203efc9353c901cebb4589e7585e2c99e19a5885c9ec3f9ddbbc10aab0ab'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971949.cpg_report.txt.gz', 'expected_bytes': 262852523, 'expected_sha256': '8935a9a92cc18a6451ecabf275901d6071a3eb30a7218bd7c1db695c77b269fe'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971950.cpg_report.txt.gz', 'expected_bytes': 262658188, 'expected_sha256': '81a5ffedce6706cee64f0bee99bc1f145042ea1c8fcf8755dbc5f9414e8bd2a2'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971951.cpg_report.txt.gz', 'expected_bytes': 266345342, 'expected_sha256': '944dd645f98a19c1dbeca95de1740d1eb49b4c883562f91947e4320d246074d3'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971952.cpg_report.txt.gz', 'expected_bytes': 262298260, 'expected_sha256': '735cb407a7c7be2b403555d9b972bac3048275ca89c8d42f48e9ceb82fae460b'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971953.cpg_report.txt.gz', 'expected_bytes': 244560070, 'expected_sha256': '8c842fc5b49e7e88f301f7919206d09b7e6c5a3f4cc1f28ccdbb882a0e70b8e6'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971954.cpg_report.txt.gz', 'expected_bytes': 268523047, 'expected_sha256': 'd5b4e564ea2aedee7488453e17171fb921fcbc0288c4b4ffcbddcb2e828aace7'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971955.cpg_report.txt.gz', 'expected_bytes': 268286531, 'expected_sha256': '2d3443cea658f2a4c45e23c405bd8532cf58b75d59ec9dfc91c6975be8c366c1'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971956.cpg_report.txt.gz', 'expected_bytes': 269230091, 'expected_sha256': '29c1ac2142e6c3b4a009fbd09d5293a1119a7f2c1ed097c0ce770d6d06d0c20d'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971957.cpg_report.txt.gz', 'expected_bytes': 266725679, 'expected_sha256': 'dec0a66bf0b8cd9cb2d4d155cfff8755bbfb76b96cedad83039e2f2985433720'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971958.cpg_report.txt.gz', 'expected_bytes': 268041654, 'expected_sha256': 'db3a481c398af7339e050c3b9bc2a712cbe4b0e2e007309286e1611b9bcd9967'}, {'path': '01_Raw_Data/DNA_methylation/Post_mortem_brain/GSE109875/WGBS/GSE109875/GSM2971959.cpg_report.txt.gz', 'expected_bytes': 269663171, 'expected_sha256': 'a435a07f4773fa2b678b423b8e38853c217c03acfd5344bf1ce16a585d2ce4ac'}]

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
    parser = argparse.ArgumentParser(description="Download repository inputs for GSE109875.")
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
            'accession': 'GSE109875',
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
        pd.DataFrame(log_rows).to_csv(os.path.join(log_dir, "GSE109875_download_log.csv"), index=False)
        checks_df = pd.DataFrame(checks)
        checks_df.to_csv(os.path.join(log_dir, "GSE109875_file_check.csv"), index=False)
        
        failures = checks_df[~checks_df['present'] | ~checks_df['size_match']]
        if not failures.empty:
            print(f"WARNING: {len(failures)} expected file(s) failed verification.")
        else:
            print(f"Dataset GSE109875 staged and verified successfully.")

if __name__ == "__main__":
    main()
