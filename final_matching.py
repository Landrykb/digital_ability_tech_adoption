import time
import pandas as pd
from shapely.geometry import Polygon, Point
from geopy.geocoders import Nominatim
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import folium
from folium import GeoJson
import unidecode
import numpy as np
import geopandas as gpd
from dbscan import DBSCAN

# Function to handle geocoding with retries and increasing delay
def do_geocode_with_retry(address, geolocator, attempt=1, max_attempts=5):
    try:
        return geolocator.reverse(address)
    except Exception as e:
        if attempt <= max_attempts:
            time.sleep(2 ** attempt)  # Increase the delay exponentially
            return do_geocode_with_retry(address, geolocator, attempt=attempt + 1)
        raise
    except Exception as e:
        print(f"An error occurred: {e}")
        return None  # Return None or handle the error appropriately


def create_city_polygon_cluster(tech_coordinates_dict, optimal_n_clusters, buffer_factor=0.05, use_dbscan=False):
    all_coords = np.concatenate(list(tech_coordinates_dict.values()))
    all_coords_array = np.array(all_coords)

    if all_coords and optimal_n_clusters:
        if use_dbscan:
            # Use DBSCAN for clustering
            dbscan = DBSCAN(eps=0.1, min_samples=5).fit(all_coords_array)
            dbscan.fit(all_coords_array)
            labels = dbscan.labels_
        else:
            # Use KMeans for clustering
            kmeans = KMeans(n_clusters=optimal_n_clusters, n_init=10, random_state=0)
            kmeans.fit(all_coords_array)
            labels = kmeans.labels_

        # Extract the cluster center
        center_point = np.mean(all_coords_array, axis=0)

        # Create a shapely Point object from the center
        center = Point(center_point)

        # Create a buffer around the point to form the city polygon with expanded boundaries
        city_polygon = center.buffer(buffer_factor)  # Adjust the buffer size as needed

        # Calculate area
        area = city_polygon.area
        return city_polygon, area
    return None, None


# Function to find optimal clusters for technology using silhouette score
def find_optimal_clusters_for_technology(coordinates_dict, use_dbscan=False):
    silhouette_scores = []
    all_coords = []

    for tech, coords_list in coordinates_dict.items():
        if not coords_list or any(len(coords) < 2 for coords in coords_list):
            print(f"No valid data for clustering in {tech}.")
            continue

        # Concatenate all coordinates for the current technology
        all_coords.extend(coords_list)

    if all_coords:
        all_coords = np.concatenate(all_coords)
        all_coords_array = np.array(all_coords)

        for i in range(2, min(11, len(all_coords))):  # Start from 2 clusters, up to the number of samples
            if use_dbscan:
                dbscan = DBSCAN(eps=0.1, min_samples=5).fit(all_coords_array)
                dbscan.fit(all_coords_array)
                labels = dbscan.labels_
            else:
                kmeans = KMeans(n_clusters=i, init='k-means++', max_iter=300, n_init=10, random_state=0)
                kmeans.fit(all_coords_array)
                labels = kmeans.labels_

            # Check if there are enough unique labels for silhouette score calculation
            unique_labels = len(np.unique(labels))
            if 2 <= unique_labels < len(all_coords_array):  # Require 2 to (n_samples - 1) unique labels
                silhouette_scores.append(silhouette_score(all_coords_array, labels))
            else:
                print(
                    f"Skipping silhouette score calculation for {i} clusters as there are not enough unique labels or too many clusters.")

    if silhouette_scores:
        # Plot silhouette scores
        plt.plot(range(2, min(11, len(all_coords))), silhouette_scores, marker='o')  # Update x-axis range
        plt.title('Silhouette Score')
        plt.xlabel('Number of clusters')
        plt.ylabel('Silhouette Score')
        plt.show()

        optimal_n_clusters = int(input("Enter the optimal number of clusters based on the plots: "))
        return optimal_n_clusters
    else:
        print("No valid data for clustering.")
        return None

# Function to resize excessively large polygons
def resize_large_polygon(polygon, threshold_area):
    if polygon.area > threshold_area:
        target_area = threshold_area / 2
        minx, miny, maxx, maxy = polygon.bounds
        new_polygon_coords = [
            (minx, miny),
            (minx, miny + target_area),
            (minx + target_area, miny + target_area),
            (minx + target_area, miny)
        ]
        return Polygon(new_polygon_coords)
    return polygon


# Initialize geolocator
geolocator = Nominatim(user_agent="myGeoapi")

# Load the LSMS individual data
key_farmers_data_for_polygons = pd.read_stata("new_key_farmers_data_with_matching_status223.dta")

# Prepare city polygons
cluster_polygons = {}
cluster_polygon_areas = {}
#
# # Function to create adoption clusters and polygons for a given technology
# def create_adoption_clusters(data, technology_column, buffer_factor=0.05):
#     adopters = data[data[technology_column] == 1]
#     adopter_coords = list(zip(adopters['GPS_Long'], adopters['GPS_Latitude']))
#
#     if adopter_coords:
#         city_polygon, area = create_city_polygon_cluster(adopter_coords, buffer_factor)
#         if city_polygon:
#             cluster_polygons[technology_column] = city_polygon
#             cluster_polygon_areas[technology_column] = area
#             data[f'{technology_column}_Polygon'] = city_polygon
#             data[f'{technology_column}_Polygon_Area'] = area
#             print(f"Fetched adoption cluster for {technology_column}, Area: {area} sq. units")
#         else:
#             print(f"Failed to create adoption cluster for {technology_column}")
#     else:
#         print(f"No adopters found for {technology_column}")
#
# #
# # Function to create adoption clusters for "Any Technology" and "Non-Adopters"
# def create_combined_clusters(data, technologies, cluster_name, buffer_factor=0.05):
#     adopters_coords = list(zip(
#         *[data[data[tech].astype(bool)]['GPS_Long'] for tech in technologies],
#         *[data[data[tech].astype(bool)]['GPS_Latitude'] for tech in technologies]
#     ))
#
#     if adopters_coords:
#         city_polygon, area = create_city_polygon_cluster(adopters_coords, buffer_factor)
#         if city_polygon:
#             cluster_polygons[cluster_name] = city_polygon
#             cluster_polygon_areas[cluster_name] = area
#             data[f'{cluster_name}_Polygon'] = city_polygon
#             data[f'{cluster_name}_Polygon_Area'] = area
#             print(f"Fetched adoption cluster for {cluster_name}, Area: {area} sq. units")
#         else:
#             print(f"Failed to create adoption cluster for {cluster_name}")
#     else:
#         print(f"No adopters found for {cluster_name}")
#

# Function to calculate dynamic buffer size based on the number of farmers per cluster
def calculate_dynamic_buffer_size(farmers_count):
    # Adjust the multiplier as needed
    return 0.01 * farmers_count  # You can change 0.01 to another value depending on your preferences


# Iterate through individual technologies and create adoption clusters
technologies = ['Any_tech', 'SRI_Adoption', 'SmartValley_Adoption', 'RiceAdvice_Adoption']

tech_coordinates_dict = {}  # Dictionary to store all coordinates for each technology

for tech in technologies:
    # Extract coordinates for the current technology
    tech_coords = key_farmers_data_for_polygons[key_farmers_data_for_polygons[tech] == 1][
        ['GPS_Long', 'GPS_Latitude']]

    # Convert tech_coords DataFrame to list of coordinates
    tech_coordinates_list = tech_coords.values.tolist()  # Convert DataFrame to list

    # Add coordinates to the dictionary
    tech_coordinates_dict[tech] = tech_coordinates_list
    print(tech_coordinates_dict[tech])

# Iterate through each technology
for tech in technologies:
    tech_cluster_column_name = f'{tech}_Cluster'
    tech_buffer_column_name = f'{tech}_Buffer'
    tech_area_column_name = f'{tech}_Area'
    tech_farmers_count_column_name = f'{tech}_Farmers_Count'
    tech_adjusted_polygon_column_name = f'{tech}_Adjusted_Polygon'

    # Initialize list to store coordinates for the current technology

    # Loop to create city polygons
    tech_cluster_polygons = {}
    tech_cluster_polygon_areas = {}
    tech_cluster_farmers_count = {}



    # Print data types and sample row values
    print("Data types:", key_farmers_data_for_polygons[['GPS_Long', 'GPS_Latitude']].dtypes)
    print("Sample row values:", key_farmers_data_for_polygons[['GPS_Long', 'GPS_Latitude']].iloc[0])
    print("Sample row types:", key_farmers_data_for_polygons[['GPS_Long', 'GPS_Latitude']].iloc[0].apply(type))

    # Find the optimal number of clusters using silhouette score
    optimal_n_clusters = find_optimal_clusters_for_technology(tech_coordinates_dict[tech], use_dbscan=True)

    print(f"Optimal number of clusters for {tech}: {optimal_n_clusters}")

    # Iterate through each row in the DataFrame
    for index, row in key_farmers_data_for_polygons.iterrows():
        # Extract relevant information from the current row
        village = unidecode.unidecode(str(row['Area_Village'])).strip().lower()
        country = row['country_match']
        gps_longitude = row['GPS_Long']
        gps_latitude = row['GPS_Latitude']
        coordinates = f"{gps_latitude}, {gps_longitude}"
        city_coordinates = [[float(gps_longitude), float(gps_latitude)]]  # Convert to list of coordinates

        # Create adoption clusters for the current technology with the optimal number of clusters
        if city_coordinates:
            # Apply the clustering
            city_polygon, area = create_city_polygon_cluster(tech_coordinates_dict[tech], optimal_n_clusters)

            if city_polygon:
                # Update dictionaries and DataFrame with cluster information
                tech_cluster_polygons[village] = city_polygon
                tech_cluster_polygon_areas[village] = area
                key_farmers_data_for_polygons.at[index, f'{tech}_Polygon'] = city_polygon
                key_farmers_data_for_polygons.at[index, tech_area_column_name] = area
                print(f"Fetched {tech} city boundary for {village}, Area: {area} sq. units")

                # Add farmer count by cluster
                farmers_count = len(
                    key_farmers_data_for_polygons[
                        key_farmers_data_for_polygons[tech_cluster_column_name] == village])
                tech_cluster_farmers_count[village] = farmers_count
                key_farmers_data_for_polygons.at[index, tech_farmers_count_column_name] = farmers_count

            else:
                print(f"Failed to create {tech} polygon for {village}")

    # Calculate dynamic threshold based on standard deviation
    if tech_cluster_polygon_areas:
        mean_area = sum(tech_cluster_polygon_areas.values()) / len(tech_cluster_polygon_areas)
        std_dev = (sum((area - mean_area) ** 2 for area in tech_cluster_polygon_areas.values()) / len(
            tech_cluster_polygon_areas)) ** 0.5
        dynamic_threshold = mean_area + 3 * std_dev  # Adjust the multiplier as needed

        # Resize excessively large polygons using the dynamic threshold
        for index, row in key_farmers_data_for_polygons.iterrows():
            village = unidecode.unidecode(str(row['Area_Village'])).strip().lower()
            if village in tech_cluster_polygons:
                city_polygon = tech_cluster_polygons[village]
                adjusted_polygon = resize_large_polygon(city_polygon, dynamic_threshold)
                if adjusted_polygon:
                    key_farmers_data_for_polygons.at[index, tech_adjusted_polygon_column_name] = adjusted_polygon
    else:
        print("No valid data for calculating dynamic threshold.")

        if f'{tech}_Polygon' in key_farmers_data_for_polygons:
            key_farmers_data_for_polygons[f'{tech}_Polygon_JSON'] = key_farmers_data_for_polygons[
                f'{tech}_Polygon'].apply(
                lambda x: x.__geo_interface__ if isinstance(x, Polygon) else None
            )
            key_farmers_data_for_polygons[f'{tech}_Polygon_JSON'] = key_farmers_data_for_polygons[
                f'{tech}_Polygon_JSON'].astype(
                'str')

            if f'{tech}_Polygon' in key_farmers_data_for_polygons.columns:
                key_farmers_data_for_polygons_data = key_farmers_data_for_polygons[
                    key_farmers_data_for_polygons[f'{tech}_Polygon'].apply(lambda x: isinstance(x, Polygon))]
                key_farmers_data_for_polygons_data.drop(columns=[f'{tech}_Polygon'], inplace=True)

                # Fix non-string labels in specific columns
                # columns_to_check = ['zae', 'region', 'sousregion', 'zae_string_recode', 'region_string_recode',
                #                     'sousregion_string_recode']
                # for column in columns_to_check:
                #     key_farmers_data_for_polygons_data[column] = key_farmers_data_for_polygons_data[column].astype(str)

                key_farmers_data_for_polygons_data.to_stata('key_farmers_data_for_polygons_data_new2.dta',
                                                            write_index=False,
                                                            version=118)
            else:
                print(f"{tech}_Polygon column does not exist in the dataset.")
        else:
            print(f"{tech}_Polygon column does not exist in the dataset.")

    # Create a Folium Map for each technology
    tech_map = folium.Map(location=[6.5244, 2.684])  # Update with the desired map center coordinates

    # Add GeoJSON layer for the city polygons
    geojson = GeoJson(
        key_farmers_data_for_polygons.set_index('Area_Village')[[f'{tech}_Polygon_JSON']],
        style_function=lambda feature: {
            'fillColor': 'green',
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.5,
        },
    ).add_to(tech_map)

    # Save the map as an HTML file
    tech_map.save(f'{tech}_map.html')
    print(f"{tech} map created and saved as {tech}_map.html")