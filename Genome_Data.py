import requests
import csv
from Bio import Entrez, SeqIO
import time
import os
import pickle
from tqdm import tqdm

# NCBI configuration
Entrez.email = "marcfrancomagana@gmail.com"

def get_phages_data(use_cache=True, cache_file="cache/phages_data.pkl", csv_file="csv/phages_data.csv"):
    """
    Retrieves phage data from the PhagesDB API with enhanced error handling and caching.
    """
    base_url = "https://phagesdb.org/api/phages/"
    all_phages = {}
    page = base_url
    failed_attempts = 0
    max_failures = 5

    if use_cache and os.path.exists(cache_file):
        try:
            with open(cache_file, "rb") as f:
                print("📁 Loading cached PhagesDB data...")
                return pickle.load(f)
        except Exception as e:
            print(f"⚠️ Error loading cache: {e}. Fetching fresh data...")

    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
    os.makedirs(os.path.dirname(csv_file), exist_ok=True)

    print("🔄 Fetching data from PhagesDB API...")

    while page and failed_attempts < max_failures:
        try:
            response = requests.get(page, timeout=30)
            response.raise_for_status()
            data = response.json()

            results = data.get("results", [])
            if not results:
                print(f"⚠️ No results found on page: {page}")
                break

            for phage in results:
                name = phage.get("phage_name", "unknown")
                if name != "unknown":
                    all_phages[name] = phage

            page = data.get("next")
            failed_attempts = 0

            print(f"📄 Processed page with {len(results)} phages. Total: {len(all_phages)}")
            time.sleep(0.2)

        except requests.exceptions.RequestException as e:
            print(f"[WARNING] Failed to access: {page}\nReason: {e}")
            failed_attempts += 1

            if failed_attempts < max_failures:
                print(f"⏳ Retrying... ({failed_attempts}/{max_failures})")
                time.sleep(2 ** failed_attempts)
            else:
                print(f"[ERROR] Maximum retry attempts ({max_failures}) reached.")
                break

    print(f"✅ Successfully fetched {len(all_phages)} phages from PhagesDB")

    if not all_phages:
        print("❌ No phage data obtained. Check API connectivity.")
        return {}

    try:
        with open(csv_file, "w", newline="", encoding='utf-8') as csvfile:
            fieldnames = [
                "phage_name", "genbank_accession", "morphotype", "found_city", "found_country",
                "found_gps_lat", "found_gps_lon", "isolation_temperature", "genome_length", "gc_percent"
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for name, data in all_phages.items():
                row_data = {}
                for field in fieldnames:
                    value = data.get(field.replace("phage_", ""), "NA")
                    if field == "phage_name":
                        value = name
                    row_data[field] = value if value not in [None, "", " "] else "NA"

                writer.writerow(row_data)
        print(f"📄 Data saved to CSV: {csv_file}")
    except Exception as e:
        print(f"⚠️ Error saving CSV: {e}")

    if use_cache:
        try:
            with open(cache_file, "wb") as f:
                pickle.dump(all_phages, f)
            print(f"💾 Data cached to: {cache_file}")
        except Exception as e:
            print(f"⚠️ Error saving cache: {e}")

    return all_phages

def get_ncbi_data(phages_dict, use_cache=True, cache_file="cache/cached_ncbi.pkl", csv_file="output/ncbi_data.csv"):
    """
    Enrich phage data with information from NCBI GenBank.
    """
    # Check if input data is provided
    if not phages_dict:
        print("❌ No phage data provided for NCBI enrichment")
        return {}
    
    # Check cache
    if use_cache and os.path.exists(cache_file):
        try:
            print("📁 Loading cached NCBI data...")
            with open(cache_file, "rb") as f:
                cached_data = pickle.load(f)
                print(f"✅ Loaded {len(cached_data)} cached NCBI records")
                return cached_data
        except Exception as e:
            print(f"⚠️ Error loading NCBI cache: {e}. Fetching fresh data...")

    # Create directories
    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
    os.makedirs(os.path.dirname(csv_file), exist_ok=True)

    ncbi_data = {}
    processed = 0
    errors = 0
    
    # Filter phages with valid accessions
    valid_phages = {name: info for name, info in phages_dict.items() 
                    if info.get("genbank_accession") not in [None, "", "NA"]}
    
    print(f"🔬 Processing {len(valid_phages)} phages with valid GenBank accessions out of {len(phages_dict)} total...")
    
    for name, info in tqdm(valid_phages.items(), desc="Fetching NCBI data"):
        accession = info.get("genbank_accession")
        
        try:
            # Fetch GenBank record
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            # Pause to respect NCBI rate limits
            time.sleep(0.4)

            # Initialize variables
            host = source = lat_lon = country = "NA"
            tail_length = locus_tag = "NA"
            cds_count = 0

            # Process record features
            for feature in record.features:
                if feature.type == "source":
                    host = feature.qualifiers.get("host", feature.qualifiers.get("lab_host", ["NA"]))[0]
                    source = feature.qualifiers.get("isolation_source", ["NA"])[0]
                    lat_lon = feature.qualifiers.get("lat_lon", ["NA"])[0]
                    country = feature.qualifiers.get("country", ["NA"])[0]

                elif feature.type == "CDS":
                    cds_count += 1
                    product = feature.qualifiers.get("product", [""])[0].lower()
                    
                    # Look for tail-related genes
                    tail_keywords = ["tape", "measure", "tail", "neck", "sheath"]
                    if any(term in product for term in tail_keywords):
                        locus_tag = feature.qualifiers.get("locus_tag", ["NA"])[0]
                        translation = feature.qualifiers.get("translation", [""])[0]
                        if translation:
                            tail_length = len(translation) * 3
                        else:
                            # If translation is missing, use sequence length
                            tail_length = len(feature.location)

            # Calculate genome statistics
            genome_length = len(record.seq)
            gc_percent = "NA"
            if genome_length > 0:
                gc_count = record.seq.count("G") + record.seq.count("C")
                gc_percent = round(gc_count / genome_length * 100, 2)

            # Retrieve host taxonomy
            taxonomy = fetch_taxonomy_from_host(host) if host != "NA" else "NA"

            # Store enriched data
            ncbi_data[name] = {
                "phage_name": name,
                "genbank_accession": accession,
                "host": host,
                "taxonomy": taxonomy,
                "source": source,
                "latitude_longitude": lat_lon,
                "country": country,
                "genome_length": genome_length,
                "GC%": gc_percent,  # Duplicate column for compatibility
                "CDS_count": cds_count,
                "host_temperature": "NA",  # Can be enriched later
                "tail_locus_tag": locus_tag,
                "tail_length_nt": tail_length
            }
            
            processed += 1

        except Exception as e:
            print(f"[ERROR] Failed to process {name} ({accession}): {e}")
            errors += 1
            continue

    print(f"✅ Successfully processed {processed} phages from NCBI ({errors} errors)")

    # Save to CSV if data exists
    if ncbi_data:
        try:
            with open(csv_file, "w", newline="", encoding='utf-8') as f:
                fieldnames = [
                    "phage_name", "genbank_accession", "host", "taxonomy", "source",
                    "latitude_longitude", "country", "genome_length", "gc_percent",
                    "GC%", "CDS_count", "host_temperature", "tail_locus_tag", "tail_length_nt"
                ]
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for data in ncbi_data.values():
                    writer.writerow(data)
            print(f"📄 NCBI data saved to CSV: {csv_file}")
        except Exception as e:
            print(f"⚠️ Error saving NCBI CSV: {e}")

        # Save cache
        if use_cache:
            try:
                with open(cache_file, "wb") as f:
                    pickle.dump(ncbi_data, f)
                print(f"💾 NCBI data cached to: {cache_file}")
            except Exception as e:
                print(f"⚠️ Error saving NCBI cache: {e}")

    return ncbi_data

def fetch_taxonomy_from_host(host_name: str) -> str:
    """
    Retrieves host taxonomic information from NCBI.
    """
    if not host_name or host_name in ["NA", "", None]:
        return "NA"

    try:
        clean_host = host_name.strip()
        if not clean_host:
            return "NA"

        search_handle = Entrez.esearch(db="taxonomy", term=clean_host, retmax=1)
        search_record = Entrez.read(search_handle)
        search_handle.close()

        if search_record["IdList"]:
            tax_id = search_record["IdList"][0]

            fetch_handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            fetch_record = Entrez.read(fetch_handle)
            fetch_handle.close()

            if fetch_record and len(fetch_record) > 0:
                lineage = fetch_record[0].get("Lineage", "Not found")
                return lineage if lineage else "Not found"

        return "Not found"

    except Exception as e:
        print(f"⚠️ Taxonomy lookup failed for '{host_name}': {e}")
        return "Error"

def merge_phages_and_ncbi_data(phages_data, ncbi_data, output_csv="output/merged_data.csv"):
    """
    Combines PhagesDB data with enriched NCBI data.
    """
    if not phages_data and not ncbi_data:
        print("❌ No data to merge")
        return {}

    combined_data = {}

    print(f"🔗 Merging {len(phages_data)} PhagesDB entries with {len(ncbi_data)} NCBI entries...")

    for phage_name in phages_data:
        phage_info = phages_data[phage_name]
        combined_info = {
            "phage_name": phage_name,
            "genbank_accession": phage_info.get("genbank_accession", "NA"),
            "morphotype": phage_info.get("morphotype", "NA"),
            "found_city": phage_info.get("found_city", "NA"),
            "found_country": phage_info.get("found_country", "NA"),
            "found_gps_lat": phage_info.get("found_gps_lat", "NA"),
            "found_gps_lon": phage_info.get("found_gps_lon", "NA"),
            "isolation_temperature": phage_info.get("isolation_temperature", "NA"),
            "genome_length": phage_info.get("genome_length", "NA"),
            "gc_percent": phage_info.get("gc_percent", "NA"),
        }

        ncbi_info = ncbi_data.get(phage_name, {})
        combined_info.update({
            "host": ncbi_info.get("host", "NA"),
            "taxonomy": ncbi_info.get("taxonomy", "NA"),
            "source": ncbi_info.get("source", "NA"),
            "latitude_longitude": ncbi_info.get("latitude_longitude", "NA"),
            "country_ncbi": ncbi_info.get("country", "NA"),
            "genome_length_ncbi": ncbi_info.get("genome_length", "NA"),
            "GC%": ncbi_info.get("GC%", "NA"),
            "CDS_count": ncbi_info.get("CDS_count", "NA"),
            "host_temperature": ncbi_info.get("host_temperature", "NA"),
            "tail_locus_tag": ncbi_info.get("tail_locus_tag", "NA"),
            "tail_length_nt": ncbi_info.get("tail_length_nt", "NA")
        })

        combined_data[phage_name] = combined_info

    print(f"✅ Merged dataset with {len(combined_data)} entries.")

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    try:
        with open(output_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(next(iter(combined_data.values())).keys()))
            writer.writeheader()
            for entry in combined_data.values():
                writer.writerow(entry)
        print(f"📄 Merged data saved to: {output_csv}")
    except Exception as e:
        print(f"⚠️ Error saving merged CSV: {e}")

    return combined_data


def filter_phages_data(phages_dict, file_name="filtered_phages.csv"):
    """
    Filters a dictionary of phages by removing those that lack values in:
    'host', 'genome_length', 'tail_length_nt', 'GC%', 'latitude_longitude', 'morphotype'
    and saves the result to a CSV file.

    :param phages_dict: Dictionary containing phage data
    :param file_name: Output CSV file name
    :return: Filtered dictionary with only phages that have all required fields
    """

    required_keys = ['host', 'genome_length', 'tail_length_nt', 'GC%', 'latitude_longitude', 'morphotype']
    filtered_phages = {}

    for phage_id, data in phages_dict.items():
        if all(data.get(key) not in [None, '', []] for key in required_keys):
            filtered_phages[phage_id] = data

    if filtered_phages:
        # Get all possible keys for the CSV header
        all_keys = list(next(iter(filtered_phages.values())).keys())

        with open(file_name, mode='w', newline='', encoding='utf-8') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=['phage_id'] + all_keys)
            writer.writeheader()

            for phage_id, data in filtered_phages.items():
                row = {'phage_id': phage_id}
                row.update(data)
                writer.writerow(row)

        print(f"CSV file saved as: {file_name}")
    else:
        print("No phages met the criteria. No file was generated.")

    return filtered_phages
