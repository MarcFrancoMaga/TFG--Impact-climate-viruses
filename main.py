from NOAA_data import *
from Genome_Data import *
from Data_analyze import *

def main():
    NOAA_API_TOKEN = "PbydWwtEqspgAwvyJyPwZvfbDZWSlIde"

    print("\n📡 Step 1: Downloading PhagesDB data...")
    phages_data = get_phages_data()
    print("\n🧬 Step 2: Fetching Genbank Accession with NCBI data...")
    ncbi_data = get_ncbi_data(phages_data)
    print("\n🔗 Step 3: Cleaning phages datasets...")
    gene_data = merge_phages_and_ncbi_data(phages_data, ncbi_data)
    gene_filtred_data = filter_phages_data(gene_data)
    
    location = generate_phage_world_map(gene_data,100)
    print("Location with most phages:", location)

    latitud = 40.445
    longitud = -79.955
    start_date = "2024-01-01"
    end_date = "2024-12-31"

    
    data_NOAA = download_noaa_climate_data(latitud,longitud,100,start_date,end_date,NOAA_API_TOKEN)
    
    print("Analyzing dataset and extracting plots...")
    print("Analyzing tail length with weather data...")
    analyze_tail_length_vs_weather(gene_filtred_data, data_NOAA)
    print("Analyzing genome length with weather data...")
    analyze_genome_length_vs_weather(gene_filtred_data, data_NOAA)
    print("Analyzing morphotype with weather data...")
    analyze_morphotype_vs_weather(gene_filtred_data, data_NOAA)
    
   
if __name__ == "__main__":
    main()
