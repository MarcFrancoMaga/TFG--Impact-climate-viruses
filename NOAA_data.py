import requests
import pandas as pd
import time
from math import radians, cos, sin, asin, sqrt

class NOAAClimateDownloader:
    def __init__(self, token):
        self.token = token
        self.base_url = "https://www.ncdc.noaa.gov/cdo-web/api/v2"
        self.headers = {"token": token}
        self.datatypes = ["TMAX", "TMIN", "PRCP", "RHAV"]

    def haversine_distance(self, lat1, lon1, lat2, lon2):
        lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
        c = 2 * asin(sqrt(a))
        r = 6371  # Earth's radius in km
        return c * r

    def get_stations_in_radius(self, lat, lon, radius_km, limit=1000):
        print(f"Searching for stations within {radius_km} km of ({lat}, {lon})...")

        lat_range = radius_km / 111.0
        lon_range = radius_km / (111.0 * cos(radians(lat)))
        extent = f"{lat - lat_range},{lon - lon_range},{lat + lat_range},{lon + lon_range}"

        url = f"{self.base_url}/stations"
        params = {
            "extent": extent,
            "datasetid": "GHCND",
            "limit": limit,
        }

        try:
            response = requests.get(url, headers=self.headers, params=params)
            response.raise_for_status()
            data = response.json()

            stations_in_radius = []
            if "results" in data:
                for station in data["results"]:
                    if station.get("latitude") and station.get("longitude"):
                        distance = self.haversine_distance(
                            lat, lon,
                            station["latitude"],
                            station["longitude"]
                        )
                        if distance <= radius_km:
                            station["distance_km"] = distance
                            stations_in_radius.append(station)

                stations_in_radius.sort(key=lambda x: x["distance_km"])

            print(f"Found {len(stations_in_radius)} stations within the specified radius.")
            return stations_in_radius

        except requests.exceptions.RequestException as e:
            print(f"Error searching for stations: {e}")
            return []

    def get_station_data(self, station_id, start_date, end_date, datatypes=None):
        if datatypes is None:
            datatypes = self.datatypes

        url = f"{self.base_url}/data"
        all_data = {}

        for datatype in datatypes:
            print(f"  Downloading {datatype} for station {station_id}...")

            params = {
                "datasetid": "GHCND",
                "stationid": station_id,
                "datatypeid": datatype,
                "startdate": start_date,
                "enddate": end_date,
                "limit": 1000,
                "units": "metric"
            }

            try:
                response = requests.get(url, headers=self.headers, params=params)
                response.raise_for_status()
                data = response.json()

                if "results" in data:
                    all_data[datatype] = data["results"]
                    if not data["results"]:
                        print(f"    No data for {datatype}.")
                else:
                    print(f"    No results for {datatype}.")
                    all_data[datatype] = []

                time.sleep(0.2)

            except requests.exceptions.RequestException as e:
                print(f"    Error downloading {datatype}: {e}")
                all_data[datatype] = []

        return all_data

    def download_climate_data(self, lat, lon, radius_km, start_date, end_date, output_filename="climate_data.csv"):
        print("=== Starting NOAA climate data download ===")
        print(f"Center point: ({lat}, {lon})")
        print(f"Radius: {radius_km} km")
        print(f"Period: {start_date} to {end_date}")
        print(f"Data types: {', '.join(self.datatypes)}\n")

        stations = self.get_stations_in_radius(lat, lon, radius_km)
        if not stations:
            print("No stations found in the specified area.")
            return {}

        all_climate_data = {}
        rows_for_csv = []

        for i, station in enumerate(stations, 1):
            station_id = station["id"]
            station_name = station.get("name", "Unnamed")
            station_lat = station.get("latitude", "N/A")
            station_lon = station.get("longitude", "N/A")
            distance = station.get("distance_km", "N/A")

            print(f"\n[{i}/{len(stations)}] Processing station:")
            print(f"  ID: {station_id}")
            print(f"  Name: {station_name}")
            print(f"  Location: ({station_lat}, {station_lon})")
            print(f"  Distance: {distance:.2f} km")

            station_data = self.get_station_data(station_id, start_date, end_date)
            dates_data = {}
            for datatype, records in station_data.items():
                if not records:
                    print(f"  No available data for {datatype} at this station.")
                for record in records:
                    date = record["date"]
                    value = record["value"]
                    if datatype in ["TMAX", "TMIN"]:
                        value = round(value, 1)
                    if date not in dates_data:
                        dates_data[date] = {
                            "station_id": station_id,
                            "station_name": station_name,
                            "station_lat": station_lat,
                            "station_lon": station_lon,
                            "distance_km": distance,
                            "date": date
                        }
                    dates_data[date][datatype] = value

            for date, data in dates_data.items():
                rows_for_csv.append(data)

            all_climate_data[station_id] = {
                "station_info": {
                    "id": station_id,
                    "name": station_name,
                    "latitude": station_lat,
                    "longitude": station_lon,
                    "distance_km": distance
                },
                "climate_data": station_data,
                "organized_data": dates_data
            }

        if rows_for_csv:
            df = pd.DataFrame(rows_for_csv)
            columns_order = ["station_id", "station_name", "station_lat", "station_lon", "distance_km", "date"] + self.datatypes
            existing_columns = [col for col in columns_order if col in df.columns]
            df = df[existing_columns]
            df = df.sort_values(["station_id", "date"])
            df.to_csv(output_filename, index=False)
            print(f"\n=== Data saved to '{output_filename}' ===")
            print(f"Total records: {len(df)}")
            print(f"Stations processed: {len(stations)}")

            station_summary = df.groupby("station_name").size()
            print("\nRecords per station:")
            for station, count in station_summary.items():
                print(f"  {station}: {count} records")

        print("\n=== Download completed ===")
        return all_climate_data

def download_noaa_climate_data(lat, lon, radius_km, start_date, end_date, token, output_filename="climate_data.csv"):
    downloader = NOAAClimateDownloader(token)
    return downloader.download_climate_data(lat, lon, radius_km, start_date, end_date, output_filename)
