import matplotlib.pyplot as plt
import os
import folium
import seaborn as sns
import pandas as pd
import math
from collections import defaultdict
from scipy.stats import pearsonr
from collections import Counter


VALID_MORPHOTYPES = ['SIPHO', 'PODO', 'MYO']

def clean_weather_data(dict_weather):
    """Elimina entradas no válidas del conjunto de datos climáticos."""
    return [r for r in dict_weather if isinstance(r, dict)]


def mapping_stations(station_id, estaciones):
    """
    Genera un mapa interactivo con todas las estaciones devueltas por `find_stations`,
    marcando la estación principal con un punto verde.

    Parámetros:
    - station_id: str. ID de la estación principal.
    - estaciones: lista de dicts, resultado de `find_stations`.
    """
    # Buscar la estación principal y su posición
    principal = next((e for e in estaciones if e["id"] == station_id), None)
    if not principal:
        print("Estación principal no encontrada en la lista.")
        return

    # Crear mapa centrado en la estación principal
    mapa = folium.Map(
        location=[principal["latitude"], principal["longitude"]],
        zoom_start=6,
        tiles="OpenStreetMap"
    )

    # Añadir marcador verde para la estación principal
    folium.Marker(
        location=[principal["latitude"], principal["longitude"]],
        popup=f"Estación principal: {principal['name']}",
        icon=folium.Icon(color="green")
    ).add_to(mapa)

    # Añadir marcadores azules para las demás estaciones
    for est in estaciones:
        if est["id"] != station_id:
            folium.Marker(
                location=[est["latitude"], est["longitude"]],
                popup=est["name"],
                icon=folium.Icon(color="blue")
            ).add_to(mapa)

    # Guardar el mapa
    mapa.save("mapa_estaciones.html")
    print("Mapa guardado como 'mapa_estaciones.html'. Ábrelo en tu navegador.")


def haversine_distance(lat1, lon1, lat2, lon2):
    """Calcula la distancia entre dos puntos geográficos en kilómetros."""
    R = 6371  # Radio de la Tierra en km
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    d_phi = math.radians(lat2 - lat1)
    d_lambda = math.radians(lon2 - lon1)

    a = math.sin(d_phi / 2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(d_lambda / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    return R * c


def generate_phage_world_map(dict_gene, range_km, output_filename="phage_world_map1.html"):
    """Genera un mapa mundial con la localización de los fagos y devuelve la localización con más fagos dentro de un rango."""
    
    color_map = {
        "SIPHO": "red",
        "MYO": "blue",
        "PODO": "green"
    }

    # Crear mapa satelital
    world_map = folium.Map(location=[0, 0], zoom_start=2, tiles="Esri.WorldImagery")
    
    coords_list = []  # Lista para almacenar las coordenadas válidas (lat, lon)

    for phage, data in dict_gene.items():
        latlon = data.get("latitude_longitude")
        morpho = data.get("morphotype", "").upper()

        if not latlon or morpho not in color_map:
            continue

        try:
            parts = latlon.strip().split()
            if len(parts) != 4:
                continue

            lat = float(parts[0]) * (1 if parts[1].upper() == 'N' else -1)
            lon = float(parts[2]) * (1 if parts[3].upper() == 'E' else -1)

            coords_list.append((lat, lon))

            popup_text = f"Phage: {phage}\nMorphotype: {morpho}"
            folium.Marker(
                location=[lat, lon],
                icon=folium.Icon(color=color_map[morpho], icon="map-marker", prefix='fa'),
                popup=popup_text
            ).add_to(world_map)

        except (ValueError, TypeError, IndexError):
            continue

    if not coords_list:
        print("No se encontraron fagos con coordenadas válidas y morphotype reconocido.")
        return None

    # Contar fagos en un rango alrededor de cada punto
    max_count = 0
    location_with_most_phages = None

    for i, (lat1, lon1) in enumerate(coords_list):
        count = 0
        for lat2, lon2 in coords_list:
            distance = haversine_distance(lat1, lon1, lat2, lon2)
            if distance <= range_km:
                count += 1
        if count > max_count:
            max_count = count
            location_with_most_phages = (lat1, lon1)

    # Añadir leyenda
    legend_html = '''
     <div style="position: fixed; 
                 bottom: 50px; left: 50px; width: 160px; height: 130px; 
                 background-color: white; border:2px solid grey; z-index:9999; 
                 font-size:14px; padding: 10px;">
     <b>Morphotype</b><br>
     <i style="background:red; width:10px; height:10px; display:inline-block;"></i> SIPHO<br>
     <i style="background:blue; width:10px; height:10px; display:inline-block;"></i> MYO<br>
     <i style="background:green; width:10px; height:10px; display:inline-block;"></i> PODO<br>
     </div>
     '''
    world_map.get_root().html.add_child(folium.Element(legend_html))

    os.makedirs("maps", exist_ok=True)
    filepath = os.path.join("maps", output_filename)
    world_map.save(filepath)
    print(f"Mapa generado y guardado en {filepath}")

    print(f"La localización con más fagos dentro de un rango de {rango_km} km es: {location_with_most_phages} con {max_count} fagos.")

    return location_with_most_phages

def analyze_tail_length_vs_weather(dict_phage, weather_dict, output_dir="tail_weather_analysis"):
    """
    Analiza la correlación entre tail_length de los fagos y los datos meteorológicos (TMAX, TMIN, PRCP, RHAV),
    utilizando directamente el diccionario weather_dict devuelto por download_climate_data.

    Args:
        dict_phage (dict): Diccionario con los datos de fagos. Debe contener la llave 'tail_length_nt'.
        weather_dict (dict): Diccionario devuelto por download_climate_data (con datos por estación).
        output_dir (str): Carpeta donde se guardarán los gráficos.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Extraer todos los tail_length válidos (convertidos a numérico, ignorando 'NA' o valores no numéricos)
    tail_lengths = []
    for phage_data in dict_phage.values():
        raw_value = phage_data.get("tail_length_nt")
        try:
            tail_length = float(raw_value)
            if tail_length > 0:
                tail_lengths.append(tail_length)
        except (TypeError, ValueError):
            continue  # Ignorar 'NA', None o cualquier valor no numérico

    if len(tail_lengths) < 2:
        print(f"❌ No hay suficientes valores de tail_length para el análisis. Valores válidos encontrados: {len(tail_lengths)}")
        return

    print(f"✅ Se encontraron {len(tail_lengths)} valores válidos de tail_length para el análisis.")

    # Preparar datos meteorológicos combinados (todos los registros de todas las estaciones)
    combined_metrics = {"TMAX": [], "TMIN": [], "PRCP": [], "RHAV": []}

    for station_id, station_info in weather_dict.items():
        organized_data = station_info.get("organized_data", {})
        for date_data in organized_data.values():
            for metric in combined_metrics:
                value = date_data.get(metric)
                try:
                    value = float(value)
                    combined_metrics[metric].append(value)
                except (TypeError, ValueError):
                    continue  # Ignorar valores no numéricos o faltantes

    # Analizar cada métrica climática
    for metric, values in combined_metrics.items():
        valid_values = [v for v in values if isinstance(v, (int, float))]

        if len(valid_values) < 2:
            print(f"⚠️ Datos insuficientes para {metric}: se encontraron {len(valid_values)} valores válidos.")
            continue

        # Emparejar tamaños de muestras
        min_len = min(len(valid_values), len(tail_lengths))
        x = valid_values[:min_len]
        y = tail_lengths[:min_len]

        corr, p_value = pearsonr(x, y)

        # Graficar
        plt.figure(figsize=(8, 6))
        sns.regplot(x=x, y=y, ci=95, line_kws={'color': 'red'})
        plt.xlabel(metric)
        plt.ylabel("Tail Length (nt)")
        plt.title(f"Tail Length vs {metric}\nPearson r = {corr:.2f}, p = {p_value:.3e}")
        plt.tight_layout()

        filename = os.path.join(output_dir, f"tail_vs_{metric}.png")
        plt.savefig(filename)
        plt.close()
        print(f"✅ Gráfico guardado: {filename}")

         
def analyze_genome_length_vs_weather(dict_phage, weather_dict, output_dir="genome_weather_analysis"):
    """
    Analiza la correlación entre genome_length de los fagos y los datos meteorológicos (TMAX, TMIN, PRCP, RHAV),
    utilizando directamente el diccionario weather_dict devuelto por download_climate_data.

    Args:
        dict_phage (dict): Diccionario con los datos de fagos. Debe contener la llave 'genome_length'.
        weather_dict (dict): Diccionario devuelto por download_climate_data (con datos por estación).
        output_dir (str): Carpeta donde se guardarán los gráficos.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Extraer todos los genome_length válidos (convertidos a numérico, ignorando 'NA' o valores no numéricos)
    genome_lengths = []
    for phage_data in dict_phage.values():
        raw_value = phage_data.get("genome_length_ncbi")
        try:
            genome_length = float(raw_value)
            if genome_length > 0:
                genome_lengths.append(genome_length)
        except (TypeError, ValueError):
            continue  # Ignorar 'NA', None o cualquier valor no numérico

    if len(genome_lengths) < 2:
        print(f"❌ No hay suficientes valores de genome_length para el análisis. Valores válidos encontrados: {len(genome_lengths)}")
        return

    print(f"✅ Se encontraron {len(genome_lengths)} valores válidos de genome_length para el análisis.")

    # Preparar datos meteorológicos combinados (todos los registros de todas las estaciones)
    combined_metrics = {"TMAX": [], "TMIN": [], "PRCP": [], "RHAV": []}

    for station_id, station_info in weather_dict.items():
        organized_data = station_info.get("organized_data", {})
        for date_data in organized_data.values():
            for metric in combined_metrics:
                value = date_data.get(metric)
                try:
                    value = float(value)
                    combined_metrics[metric].append(value)
                except (TypeError, ValueError):
                    continue  # Ignorar valores no numéricos o faltantes

    # Analizar cada métrica climática
    for metric, values in combined_metrics.items():
        valid_values = [v for v in values if isinstance(v, (int, float))]

        if len(valid_values) < 2:
            print(f"⚠️ Datos insuficientes para {metric}: se encontraron {len(valid_values)} valores válidos.")
            continue

        # Emparejar tamaños de muestras
        min_len = min(len(valid_values), len(genome_lengths))
        x = valid_values[:min_len]
        y = genome_lengths[:min_len]

        corr, p_value = pearsonr(x, y)

        # Graficar
        plt.figure(figsize=(8, 6))
        sns.regplot(x=x, y=y, ci=95, line_kws={'color': 'red'})
        plt.xlabel(metric)
        plt.ylabel("Genome Length (bp)")
        plt.title(f"Genome Length vs {metric}\nPearson r = {corr:.2f}, p = {p_value:.3e}")
        plt.tight_layout()

        filename = os.path.join(output_dir, f"genome_vs_{metric}.png")
        plt.savefig(filename)
        plt.close()
        print(f"✅ Gráfico guardado: {filename}")
        

def analyze_morphotype_vs_weather(phage_dict, weather_dict, output_dir="morphotype_weather_analysis"):
    """
    Analiza la distribución porcentual de los diferentes morphotype (SIHPO, PODO, MYO) 
    en función de las métricas climáticas (TMAX, TMIN, PRCP, RHAV), generando múltiples tipos de gráficos.

    Args:
        phage_dict (dict): Diccionario con los datos de fagos. Debe contener la llave 'morphotype'.
        weather_dict (dict): Diccionario devuelto por download_climate_data (con datos por estación).
        output_dir (str): Carpeta donde se guardarán los gráficos.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    valid_morphotypes = ["SIPHO", "PODO", "MYO"]
    morphotypes = []

    for phage_data in phage_dict.values():
        morph = phage_data.get("morphotype", "").upper()
        if morph in valid_morphotypes:
            morphotypes.append(morph)

    if len(morphotypes) < 2:
        print(f"❌ No hay suficientes valores de morphotype para el análisis. Valores válidos encontrados: {len(morphotypes)}")
        return

    print(f"✅ Se encontraron {len(morphotypes)} fagos válidos con morphotype para el análisis.")

    combined_metrics = {"TMAX": [], "TMIN": [], "PRCP": [], "RHAV": []}

    for station_id, station_info in weather_dict.items():
        organized_data = station_info.get("organized_data", {})
        for date_data in organized_data.values():
            for metric in combined_metrics:
                value = date_data.get(metric)
                try:
                    value = float(value)
                    combined_metrics[metric].append(value)
                except (TypeError, ValueError):
                    continue

    morphotype_counts = Counter(morphotypes)
    total_phages = sum(morphotype_counts.values())

    morphotype_percentages = {k: (v / total_phages) * 100 for k, v in morphotype_counts.items()}

    print("\nDistribución porcentual de morphotypes:")
    for morph, pct in morphotype_percentages.items():
        print(f"  {morph}: {pct:.2f}%")

    for metric, values in combined_metrics.items():
        valid_values = [v for v in values if isinstance(v, (int, float))]

        if len(valid_values) < 2:
            print(f"⚠️ Datos insuficientes para {metric}: se encontraron {len(valid_values)} valores válidos.")
            continue

        min_len = min(len(valid_values), len(morphotypes))
        x = valid_values[:min_len]
        y_raw = morphotypes[:min_len]

        morph_numeric = [1 if morph == "SIPHO" else 2 if morph == "PODO" else 3 for morph in y_raw]
        morph_str = y_raw

        df = pd.DataFrame({
            "Morphotype": morph_str,
            metric: x
        })

        # Boxplot
        plt.figure(figsize=(8, 6))
        sns.boxplot(x="Morphotype", y=metric, data=df)
        plt.title(f"Boxplot: {metric} según morphotype")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"boxplot_morphotype_vs_{metric}.png"))
        plt.close()

        # Violin plot
        plt.figure(figsize=(8, 6))
        sns.violinplot(x="Morphotype", y=metric, data=df)
        plt.title(f"Violin Plot: {metric} según morphotype")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"violinplot_morphotype_vs_{metric}.png"))
        plt.close()

        # Strip plot
        plt.figure(figsize=(8, 6))
        sns.stripplot(x="Morphotype", y=metric, data=df, jitter=True)
        plt.title(f"Strip Plot: {metric} según morphotype")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"stripplot_morphotype_vs_{metric}.png"))
        plt.close()

        print(f"✅ Gráficos de {metric} guardados en: {output_dir}")

    print("\n=== Análisis completado ===")

        



