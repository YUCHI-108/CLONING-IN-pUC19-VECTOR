import requests
from bs4 import BeautifulSoup
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import time, os, glob
 
BASE_URL = "http://metabolomics.jp"
LISTING_URL = "http://metabolomics.jp/mediawiki/index.php?title=Special:WhatLinksHere/Index:FL&namespace=0&limit=500"
 
checkpoint_dir = "checkpoints"
os.makedirs(checkpoint_dir, exist_ok=True)
 
skipped_entries = []
 
def get_listing_pages():
    """Collect all listing page URLs by following 'next 500' links."""
    urls = [LISTING_URL]
    while True:
        r = requests.get(urls[-1])
        soup = BeautifulSoup(r.text, "html.parser")
        next_link = soup.find("a", string=lambda x: x and "next 500" in x.lower())
        if not next_link:
            break
        urls.append(BASE_URL + next_link["href"])
    return urls
 
def extract_ids(listing_url):
    """Extract valid flavonoid IDs from a listing page."""
    r = requests.get(listing_url)
    soup = BeautifulSoup(r.text, "html.parser")
    ids = []
    for li in soup.select("ul li a"):
        fid = li.text.strip()
        if len(fid) == 12 and fid.startswith("FL"):
            ids.append(fid)
    return ids
 
def fetch_entry(fid):
    """Fetch details for a single flavonoid ID."""
    url = f"{BASE_URL}/wiki/{fid}"
    try:
        r = requests.get(url)
        soup = BeautifulSoup(r.text, "html.parser")
        smiles, systematic_name, avg_mass = "", "", ""
 
        for row in soup.select("table tr"):
            cells = [c.get_text(" ", strip=True) for c in row.find_all(["td", "th"])]
            if len(cells) >= 2:
                key, value = cells[0].lower(), cells[1]
                if "smiles" in key:
                    smiles = value
                elif "systematic name" in key:
                    systematic_name = value
                elif "average mass" in key:
                    try:
                        avg_mass = str(int(float(value.split()[0])))
                    except:
                        pass
 
        if smiles or systematic_name:
            return {"Flavonoid_ID": fid,
                    "Systematic_Name": systematic_name,
                    "SMILES": smiles,
                    "Average_Mass": avg_mass}
        else:
            skipped_entries.append(fid)
            return None
    except Exception:
        skipped_entries.append(fid)
        return None
 
def save_checkpoint(batch_id, data):
    """Save a checkpoint file for safety."""
    if not data:
        return
    df = pd.DataFrame(data)
    path = os.path.join(checkpoint_dir, f"checkpoint_{batch_id}.csv")
    df.to_csv(path, index=False)
    print(f"[Checkpoint] Saved {len(df)} rows → {path}")
 
def load_completed_batches():
    """Find which batches were already completed from checkpoint files."""
    files = sorted(glob.glob(os.path.join(checkpoint_dir, "checkpoint_*.csv")))
    completed = set()
    for f in files:
        try:
            batch_num = int(os.path.basename(f).split("_")[1].split(".")[0])
            completed.add(batch_num)
        except:
            continue
    return completed
 
def main():
    start_time = time.time()
 
    print("Collecting listing pages...")
    listing_pages = get_listing_pages()
    print(f"Found {len(listing_pages)} listing pages.")
 
    all_ids = []
    for page in listing_pages:
        ids = extract_ids(page)
        all_ids.extend(ids)
    print(f"Total flavonoid IDs collected: {len(all_ids)}")
 
    # Resume: find already completed checkpoints
    completed_batches = load_completed_batches()
    print(f"Already completed batches: {sorted(completed_batches)}")
 
    batch_size = 500
    all_data = []
 
    for i in range(0, len(all_ids), batch_size):
        batch_num = i // batch_size + 1
        batch_ids = all_ids[i:i+batch_size]
 
        if batch_num in completed_batches:
            print(f"Skipping batch {batch_num} (already done)")
            # Load existing data for merging later
            df = pd.read_csv(os.path.join(checkpoint_dir, f"checkpoint_{batch_num}.csv"))
            all_data.extend(df.to_dict("records"))
            continue
 
        print(f"\nProcessing batch {batch_num} ({len(batch_ids)} entries)...")
        batch_data = []
        with ThreadPoolExecutor(max_workers=20) as executor:
            futures = {executor.submit(fetch_entry, fid): fid for fid in batch_ids}
            for j, future in enumerate(as_completed(futures), 1):
                result = future.result()
                if result:
                    batch_data.append(result)
 
                if j % 50 == 0 or j == len(batch_ids):
                    print(f"  Progress: {j}/{len(batch_ids)}")
 
        save_checkpoint(batch_num, batch_data)
        all_data.extend(batch_data)
 
    # Merge all checkpoints into final CSV
    final_file = "flavonoids_final.csv"
    df = pd.DataFrame(all_data)
    df.to_csv(final_file, index=False)
 
    if skipped_entries:
        pd.DataFrame(skipped_entries, columns=["Skipped_IDs"]).to_csv("skipped_entries.csv", index=False)
 
    print("\n Done!")
    print(f"Collected entries: {len(all_data)}")
    print(f"Skipped entries: {len(skipped_entries)}")
    print(f"Final merged CSV: {final_file}")
    print(f"Time taken: {time.time() - start_time:.2f} sec")
 
if __name__ == "__main__":
    main()
