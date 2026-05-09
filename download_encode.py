import requests
import os
import csv
import gzip

base_url = "https://encodeproject.org/search/"
encode_url = "https://www.encodeproject.org"

def download_encode_data(params):
    try:
        # get json with file/experiment names from search parameters
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        print(f"Failed to retrieve data: {e}")
        return {}
    

def download_records(records, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    for i, rec in enumerate(records, start=1):
        # file to be written to
        target = os.path.join(out_dir, rec["filename"])
        if os.path.exists(target):
            print(f"[{i}/{len(records)}] Skipping existing file: {os.path.basename(target)}")
            continue

        print(f"[{i}/{len(records)}] Downloading {rec['accession']}: {os.path.basename(target)}")

        try:
            with requests.get(rec["url"], stream=True, timeout=120) as response:
                response.raise_for_status()
                with open(target, "wb") as f:
                    # read in chunks instead of writing all at once
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)

        except requests.RequestException as e:
            print(f"Failed to download {rec['accession']}: {e}")


def get_file_records(data):
    records = []
    for item in data.get("@graph", []):
        href = item.get("href")
        accession = item.get("accession")
        if not href or not accession:
            continue

        # get file name and other useful metadata
        filename = href.split("/@@download/")[-1]
        records.append(
            {
                "accession": accession,
                "url": encode_url + href,
                "filename": filename,
                "assembly": item.get("assembly", ""),
                "file_format": item.get("file_format", ""),
                "file_format_type": item.get("file_format_type", ""),
                "output_type": item.get("output_type", ""),
                "output_category": item.get("output_category", ""),
                "dataset": item.get("dataset", ""),
            }
        )
    return records


def write_records_csv(records, csv_path):
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    # manifest columns to keep
    fieldnames = [
        "accession",
        "filename",
        "url",
        "assembly",
        "file_format",
        "file_format_type",
        "output_type",
        "output_category",
        "dataset",
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            writer.writerow({k: rec.get(k, "") for k in fieldnames})


def merge_wgbs_to_csv(records, download_dir, output_csv, max_files=None, max_rows=None, zip=True):
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    merged_files = 0
    merged_rows = 0

    # handle writing to zipped and unzipped files
    if output_csv.endswith(".gz") and zip:
        out_handle = gzip.open(output_csv, "wt", newline="", encoding="utf-8")
    else:
        out_handle = open(output_csv, "w", newline="", encoding="utf-8")

    with out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(
            [
                "accession",
                "chrom",
                "start",
                "end",
                "strand",
                "methylation_percent"  # likely the last column according to https://www.encodeproject.org/data-standards/wgbs/
            ]
        )

        for rec in records:
            filename = rec["filename"]

            # can use bigBed files with UCSC's bigBedToBed, but we already have a lot
            # of WGBS data so i'll use only bed files for now
            if not (filename.endswith(".bed") or filename.endswith(".bed.gz")):
                continue

            # set download limit
            if max_files is not None and merged_files >= max_files:
                break

            input_path = os.path.join(download_dir, filename)

            if not os.path.exists(input_path):
                continue

            # handle reading from zipped and unzipped files
            if input_path.endswith(".gz"):
                in_handle = gzip.open(input_path, "rt", encoding="utf-8", errors="replace")
            else:
                in_handle = open(input_path, "r", encoding="utf-8", errors="replace")

            with in_handle:
                for line in in_handle:
                    # skip empty lines, comments, and track/browser lines
                    if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                        continue

                    fields = line.rstrip("\n").split("\t")
                    # bed files must have chrom, start, end - skip if malformed
                    if len(fields) < 3:
                        continue

                    if not (fields[1].isdigit() and fields[2].isdigit()):
                        continue

                    # get strand and methylation info if it exists
                    strand = fields[5] if len(fields) > 5 else ""
                    methylation_percent = fields[10] if len(fields) > 10 else ""

                    writer.writerow(
                        [
                            rec["accession"],
                            fields[0],
                            fields[1],
                            fields[2],
                            strand,
                            methylation_percent,
                        ]
                    )
                    merged_rows += 1

                    if max_rows is not None and merged_rows >= max_rows:
                        break

            merged_files += 1

            if max_rows is not None and merged_rows >= max_rows:
                break


def merge_rnaseq_to_csv(records, download_dir, output_csv, max_files=None, max_rows=None, zip=True):
    # very similar to merging WGBS except for parsing columns
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    merged_files = 0
    merged_rows = 0

    if output_csv.endswith(".gz") and zip:
        out_handle = gzip.open(output_csv, "wt", newline="", encoding="utf-8")
    else:
        out_handle = open(output_csv, "w", newline="", encoding="utf-8")

    with out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["accession", "gene_id", "expected_count", "tpm", "fpkm"])

        for rec in records:
            filename = rec["filename"]
            # should have only tsv files
            if not (filename.endswith(".tsv") or filename.endswith(".tsv.gz")):
                continue

            if max_files is not None and merged_files >= max_files:
                break

            input_path = os.path.join(download_dir, filename)
            if not os.path.exists(input_path):
                continue

            if input_path.endswith(".gz"):
                in_handle = gzip.open(input_path, "rt", encoding="utf-8", errors="replace")
            else:
                in_handle = open(input_path, "r", encoding="utf-8", errors="replace")

            with in_handle:
                header = None
                for line in in_handle:
                    if not line.strip() or line.startswith("#"):
                        continue
                    header = line.strip().split()
                    break

                if not header:
                    continue

                index = {name: i for i, name in enumerate(header)}
                gene_idx = index.get("gene_id", 0)
                counts_idx = index.get("expected_count", 1 if len(header) > 1 else 0)
                tpm_idx = index.get("TPM", -1)
                fpkm_idx = index.get("FPKM", -1)

                for line in in_handle:
                    if not line.strip() or line.startswith("#"):
                        continue

                    fields = line.strip().split()
                    if len(fields) <= gene_idx:
                        continue

                    gene_id = fields[gene_idx]
                    expected_count = fields[counts_idx] if len(fields) > counts_idx else ""
                    tpm = fields[tpm_idx] if tpm_idx >= 0 and len(fields) > tpm_idx else ""
                    fpkm = fields[fpkm_idx] if fpkm_idx >= 0 and len(fields) > fpkm_idx else ""

                    writer.writerow([rec["accession"], gene_id, expected_count, tpm, fpkm])
                    merged_rows += 1

                    if max_rows is not None and merged_rows >= max_rows:
                        break

            merged_files += 1

            if max_rows is not None and merged_rows >= max_rows:
                break


rna_seq_params = {
    "type": "File",
    "searchTerm": "RNA-seq",
    "output_category": "quantification",
    "output_type": "gene quantifications",
    "biosample_ontology.cell_slims": "fibroblast",
    "assembly": "GRCh38",
    "status": "released",
    "limit": "all",
    "format": "json"
}

wgbs_params = {
    "type": "File",
    "searchTerm": "WGBS",
    "output_category": "quantification",
    "output_type": "methylation state at CpG",
    "biosample_ontology.cell_slims": "fibroblast",
    "assembly": "GRCh38",
    "status": "released",
    "limit": "all",
    "format": "json"
}

rna_seq_data = download_encode_data(rna_seq_params)
# wgbs_data = download_encode_data(wgbs_params)

rna_seq_file_ids = [
    item["accession"] for item in rna_seq_data.get("@graph", []) if "accession" in item
]
# wgbs_file_ids = [
#     item["accession"] for item in wgbs_data.get("@graph", []) if "accession" in item
# ]

rna_seq_records = get_file_records(rna_seq_data)
# wgbs_records = get_file_records(wgbs_data)

print(f"RNA-seq file IDs found: {len(rna_seq_file_ids)}")
print(rna_seq_file_ids[:10])

# print(f"WGBS file IDs found: {len(wgbs_file_ids)}")
# print(wgbs_file_ids[:10])

write_records_csv(rna_seq_records, "exports/rna_seq_manifest.csv")
# write_records_csv(wgbs_records, "exports/wgbs_manifest.csv")
print("Wrote manifest CSV files to exports/")


# uncomment values to produce the test files
MERGE_MAX_FILES = 5  # None  # 1
MERGE_MAX_ROWS = 100000  # None


download_records(rna_seq_records, "downloads/rna_seq")
# download_records(wgbs_records, "downloads/wgbs")

# merge_wgbs_to_csv(
#     wgbs_records,
#     "downloads/wgbs",
#     "exports/wgbs_merged.csv",
#     max_files=MERGE_MAX_FILES,
#     max_rows=MERGE_MAX_ROWS,
#     zip=True
# )

merge_rnaseq_to_csv(
    rna_seq_records,
    "downloads/rna_seq",
    "exports/rna_seq_merged.csv",
    max_files=None,
    max_rows=MERGE_MAX_ROWS,
    zip=True
)

# The WGBS files generated by this script ended up not being used in the old analysis.
