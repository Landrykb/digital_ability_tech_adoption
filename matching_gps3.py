import time
import pandas as pd
from shapely.geometry import shape
from geopy.geocoders import Nominatim
from shapely.geometry import Point, Polygon
import folium
from selenium import webdriver
from geopy.exc import GeocoderTimedOut
import unidecode
from geopy.distance import geodesic
import geopandas as gpd
from folium.plugins import MeasureControl
import json
from shapely.wkt import loads
import pandas as pd
import plotly.graph_objects as go
from shapely.wkt import loads


# Function to handle geocoding with retries and increasing delay
def do_geocode_with_retry(address, geolocator, attempt=1, max_attempts=5):
    try:
        return geolocator.reverse(address)
    except GeocoderTimedOut as e:
        if attempt <= max_attempts:
            time.sleep(2 ** attempt)  # Increase the delay exponentially
            return do_geocode_with_retry(address, geolocator, attempt=attempt + 1)
        raise
    except Exception as e:
        print(f"An error occurred: {e}")
        # Handle other specific exceptions as needed
        return None  # Return None or handle the error appropriately


# Function to encode city name and handle Unicode characters
def encode_city_name(city_name):
    if city_name:
        city_name = unidecode.unidecode(city_name)
        return city_name.strip().lower()
    return None


# Function to fetch city boundaries from OpenStreetMap using Nominatim
def fetch_city_boundary(city_name, geolocator):
    location = do_geocode_with_retry(city_name, geolocator)
    if location:
        return location.raw.get('boundingbox')
    return None


# Function to check if coordinates fall within a polygon boundary
def check_within_boundary(point, polygon):
    return polygon.contains(point)


# Initialize geolocator
geolocator = Nominatim(user_agent="myGeoapi")

lsms_data = pd.read_stata("ben_civ_co_dataset_farmers_aggDigLit.dta")

filtered_data = lsms_data[(lsms_data['country'] == 'BEN') | (lsms_data['country'] == 'CIV')]

data_ben = filtered_data[filtered_data['country'] == 'BEN']
data_civ = filtered_data[filtered_data['country'] == 'CIV']

unique_data_ben = data_ben.drop_duplicates(subset=['country', 'grappe']).copy()
lenght1 = len(unique_data_ben)
print(lenght1)
unique_data_civ = data_civ.drop_duplicates(subset=['country', 'grappe']).copy()
lenght2 = len(unique_data_civ)
print(lenght2)


columns_to_check = ['coordonnes_gps__Longitude', 'coordonnes_gps__Latitude']
unique_data_ben = unique_data_ben.dropna(subset=columns_to_check)
unique_data_civ = unique_data_civ.dropna(subset=columns_to_check)

unique_lsms_data = pd.concat([unique_data_ben, unique_data_civ], ignore_index=True)

city_polygons = {}
for index, row in unique_lsms_data.iterrows():
    lsms_city = encode_city_name(row['sousregion_string'])
    country = row['country']  # Assuming 'country' column contains country code
    gps_longitude = row['coordonnes_gps__Longitude']
    gps_latitude = row['coordonnes_gps__Latitude']
    coordinates = f"{gps_latitude}, {gps_longitude}"  # Create a coordinate pair

    location = do_geocode_with_retry(coordinates, geolocator)
    if location:
        fetched_city = encode_city_name(location.raw.get('address', {}).get('city'))

        city_boundary = fetch_city_boundary(coordinates, geolocator)
        if city_boundary:
            city_polygon_coords = [
                (float(city_boundary[2]), float(city_boundary[0])),
                (float(city_boundary[2]), float(city_boundary[1])),
                (float(city_boundary[3]), float(city_boundary[1])),
                (float(city_boundary[3]), float(city_boundary[0]))
            ]
            city_polygons[lsms_city] = Polygon(city_polygon_coords)
            city_polygon = Polygon(city_polygon_coords)
            unique_lsms_data.at[index, 'City_Polygon'] = city_polygon
            print(f"Fetched with GPS coordinates; {lsms_city} and village '{fetched_city}' Matched")
        else:
            coordinates_with_country = f"{coordinates}, {country}"
            location_with_country = do_geocode_with_retry(coordinates_with_country, geolocator)
            if location_with_country:
                city_boundary = fetch_city_boundary(coordinates_with_country, geolocator)
                if city_boundary:
                    city_polygons[lsms_city] = Polygon([(float(city_boundary[2]), float(city_boundary[0])),
                                                        (float(city_boundary[2]), float(city_boundary[1])),
                                                        (float(city_boundary[3]), float(city_boundary[1])),
                                                        (float(city_boundary[3]), float(city_boundary[0]))])
                    unique_lsms_data.at[index, 'City_Polygon'] = Polygon(
                        [(float(city_boundary[2]), float(city_boundary[0])),
                         (float(city_boundary[2]), float(city_boundary[1])),
                         (float(city_boundary[3]), float(city_boundary[1])),
                         (float(city_boundary[3]), float(city_boundary[0]))])
                    print(f"Fetched with city/village names; {lsms_city}, {country} and {fetched_city}, {country} matched")
                else:
                    print(
                        f"No city boundary found for GPS coordinates with country; {lsms_city} and {fetched_city} NOT "
                        f"Matched")
            else:
                print(f"No city information fetched for GPS coordinates; {lsms_city} and {fetched_city}")
    else:
        print(f"No city information fetched for GPS coordinates; {gps_latitude}, {gps_longitude}")

if 'City_Polygon' in unique_lsms_data:
    unique_lsms_data['City_Polygon_JSON'] = unique_lsms_data['City_Polygon'].apply(
        lambda x: x.__geo_interface__ if isinstance(x, Polygon) else None
    )
    unique_lsms_data['City_Polygon_JSON'] = unique_lsms_data['City_Polygon_JSON'].astype('str')

    unique_lsms_data = unique_lsms_data[unique_lsms_data['City_Polygon'].apply(lambda x: isinstance(x, Polygon))]
    unique_lsms_data.drop(columns=['City_Polygon'], inplace=True)
else:
    print("City_Polygon column does not exist in the dataset.")

unique_lsms_data.to_stata('lsms_data_with_polygons_full.dta', write_index=False, version=118)

# Load the dataset containing city polygons

# Load the dataset containing city polygons
lsms_data_with_polygons = pd.read_stata("lsms_data_with_polygons_full.dta")

# Convert City_Polygon_JSON strings back to Polygon objects
lsms_data_with_polygons['City_Polygon'] = lsms_data_with_polygons['City_Polygon_JSON'].apply(lambda x: shape(eval(x)))

# Load Key farmers dataset
key_farmers_data = pd.read_stata("3_Dataset_StataFile.dta")

# Set 'Matching_Status' and 'Matched_Polygon' columns
key_farmers_data['Matching_Status'] = 0
key_farmers_data['Matched_Polygon'] = None
key_farmers_data['country_match'] = None
key_farmers_data['grappe_match'] = None

# Initialize a dictionary to store country, grappe, and matching polygons
match_info = {}

# Access the City_Polygon column (now containing Polygon objects)
city_polygons = lsms_data_with_polygons['City_Polygon']

for index, row in key_farmers_data.iterrows():
    farmer_coords = Point(row['GPS_Long'], row['GPS_Latitude'])
    farmer_village = row['Area_Village'].strip().lower()

    potential_matches = [city_polygon for city_polygon in city_polygons if farmer_coords.within(city_polygon)]

    # Count the number of potential matches
    num_potential_matches = len(potential_matches)
    key_farmers_data.at[index, 'Num_Potential_Matches'] = num_potential_matches  # Add count to the dataset

    if potential_matches:
        # Check containment within each potential match
        contained_polygons = [polygon for polygon in potential_matches if farmer_coords.within(polygon)]

        if contained_polygons:
            # Farmer coordinates are within one of the polygons
            distances_from_center = [farmer_coords.distance(polygon.centroid) for polygon in contained_polygons]
            closest_polygon_center = contained_polygons[distances_from_center.index(min(distances_from_center))]

            # Assign the matched polygon to the farmer
            key_farmers_data.at[index, 'Matching_Status'] = 1
            key_farmers_data.at[index, 'Matched_Polygon'] = closest_polygon_center

            matching_row = lsms_data_with_polygons[lsms_data_with_polygons['City_Polygon'] == closest_polygon_center]

            if not matching_row.empty:
                country_value = matching_row['country'].values[0]
                grappe_value = matching_row['grappe'].values[0]

                match_info[index] = {
                    'country_match': country_value,
                    'grappe_match': grappe_value,
                    'matched_polygon': closest_polygon_center
                }

            print(f"Coordinates for '{farmer_village}' within one of the polygons, assigned based on center distance")
        else:
            # Calculate distances from edges and area proportions for potential matches
            distances = [farmer_coords.distance(polygon.exterior) for polygon in potential_matches]
            areas = [polygon.area for polygon in potential_matches]

            # Normalize areas to get proportion
            total_area = sum(areas)
            proportions = [area / total_area for area in areas]

            # Weighted distance considering both distance to edges and area proportion
            weighted_distances = [distance * proportion for distance, proportion in zip(distances, proportions)]

            # Find the nearest boundary considering both distance and area proportion
            nearest_polygon = potential_matches[weighted_distances.index(min(weighted_distances))]

            key_farmers_data.at[index, 'Matching_Status'] = 1
            key_farmers_data.at[index, 'Matched_Polygon'] = nearest_polygon  # Assign the matched polygon
            print(f"Coordinates for '{farmer_village}' assigned considering proximity and area weighting")

            matching_row = lsms_data_with_polygons[lsms_data_with_polygons['City_Polygon'] == nearest_polygon]

            if not matching_row.empty:
                country_value = matching_row['country'].values[0]
                grappe_value = matching_row['grappe'].values[0]

                match_info[index] = {
                    'country_match': country_value,
                    'grappe_match': grappe_value,
                    'matched_polygon': nearest_polygon
                }

    else:
        print(f"Coordinates did NOT match within the city or within close range of any city boundary")



# Update key_farmers_data using match_info dictionary after the loop
for index, info in match_info.items():
    key_farmers_data.at[index, 'country_match'] = info['country_match']
    key_farmers_data.at[index, 'grappe_match'] = info['grappe_match']
    key_farmers_data.at[index, 'Matched_Polygon'] = info['matched_polygon']



# ... Your initial data loading and processing code ...

# Convert 'grappe_match' column to string
key_farmers_data['grappe_match'] = key_farmers_data['grappe_match'].astype(str)

# Filter out rows with null 'Matched_Polygon' values
key_farmers_data = key_farmers_data.dropna(subset=['Matched_Polygon'])

# Convert 'Matched_Polygon' to different formats
formats_to_try = ['wkt', 'json']  # Add more formats as needed

for format_type in formats_to_try:
    try:
        if format_type == 'wkt':
            key_farmers_data['Matched_Polygon_WKT'] = key_farmers_data['Matched_Polygon'].apply(lambda x: x.wkt if x else None)
        elif format_type == 'json':
            key_farmers_data['Matched_Polygon_JSON'] = key_farmers_data['Matched_Polygon'].apply(lambda x: x.__geo_interface__ if x else None)

        # Drop the original 'Matched_Polygon' column after conversion
        key_farmers_data.drop(columns=['Matched_Polygon'], inplace=True)

        # Save the updated dataset as CSV
        key_farmers_data.to_csv('key_farmers_data_with_matching_status22.csv', index=False)

        # Save the updated dataset as STATA
        key_farmers_data.to_stata('key_farmers_data_with_matching_status22.dta', write_index=False, version=118)

        print(f"Data exported successfully in {format_type.upper()} format!")
        # Drop the remaining problematic columns after successful export
        # key_farmers_data.drop(columns=['Matched_Polygon_WKT', 'Matched_Polygon_JSON'], inplace=True)
        break  # Exit the loop if successful
    except Exception as e:
        print(f"Error occurred during export attempt in {format_type.upper()} format: {e}")
        continue  # Try the next format
else:
    raise RuntimeError("All export attempts failed! Check the data and try again.")

key_farmers_data.to_stata('key_farmers_data_with_matching_status22.dta', write_index=False, version=118)




lsms_data_with_polygons = pd.read_stata("lsms_data_with_polygons_full22.dta")
key_farmers_data_with_matching_status = pd.read_stata("key_farmers_data_with_matching_status22.dta")


def check_wkt_format(wkt_string):
    try:
        geom = loads(wkt_string)
        # If the load succeeds, it's a valid WKT string
        return True
    except Exception as e:
        # If the load fails, it's likely due to a format issue
        print(f"Error: {e}")
        return False


# Apply the function to your dataset column containing WKT strings
for wkt_string in key_farmers_data_with_matching_status['Matched_Polygon_WKT']:
    is_valid = check_wkt_format(wkt_string)
    if not is_valid:
        print(f"Invalid WKT: {wkt_string}")
key_farmers_data_with_matching_status['Matched_Polygon'] = key_farmers_data_with_matching_status[
    'Matched_Polygon_WKT'].apply(lambda x: shape(eval(x)) if x else None)
# key_farmers_data_with_matching_status.drop(columns=['Matched_Polygon_JSON'], inplace=True)
key_farmers_data_with_matching_status.to_stata('key_farmers_data_with_matching_status22.dta', write_index=False,
                                               version=118)

key_farmers_data_matched = pd.read_stata("key_farmers_data_with_matching_status22.dta")
m = folium.Map(location=[6.5244, 2.6844], zoom_start=8)

crs = 'EPSG:4326'  # Replace with your desired CRS
gdf_polygons = gpd.GeoDataFrame(key_farmers_data_with_matching_status, geometry='Matched_Polygon', crs=crs)

# Convert GeoDataFrame to GeoJSON and add to the map
gdf_polygons_json = gdf_polygons.to_json()


folium.GeoJson(gdf_polygons_json, name='Matched City Polygons').add_to(m)

# # Add markers for matched farmers
matched_farmers = key_farmers_data_matched[key_farmers_data_matched['Matching_Status'] == 1]
# for idx, row in matched_farmers.iterrows():
#     folium.Marker([row['GPS_Latitude'], row['GPS_Long']], popup=row['Area_Village'],
#                   icon=folium.Icon(color='green')).add_to(m)
# Add circle markers for matched farmers
for idx, row in matched_farmers.iterrows():
    folium.Circle(
        location=[row['GPS_Latitude'], row['GPS_Long']],
        radius=50,  # Adjust the radius as needed
        color='#191970',  # Color of the circle
        fill=True,
        fill_color='#191970'  # Fill color of the circle
    ).add_to(m)


# Add a legend with scale
m.add_child(folium.LatLngPopup())
m.add_child(folium.LayerControl())
m.add_child(folium.plugins.MeasureControl(position='topleft', primary_length_unit='meters',
                                          secondary_length_unit='kilometers'))

measure_control = MeasureControl(position='topleft', primary_length_unit='meters', secondary_length_unit='kilometers')
m.add_child(measure_control)

# Save the map as HTML
m.save('map_with_matched_farmers22.html')

# Open the HTML map in a browser using Selenium
browser = webdriver.Chrome()  # Adjust this based on your browser and its driver
browser.get('file:///path_to/map_with_matched_farmers22.html')
# Assuming key_farmers_data_with_matching_status has 'Matched_Polygon_WKT' and necessary location columns



# Convert WKT to Shapely geometries using Shapely's loads method
key_farmers_data_with_matching_status['Matched_Polygon'] = key_farmers_data_with_matching_status[
    'Matched_Polygon_WKT'].apply(lambda x: loads(x) if x else None)

# Create GeoDataFrame
crs = 'EPSG:4326'
gdf_polygons = gpd.GeoDataFrame(key_farmers_data_with_matching_status, geometry='Matched_Polygon', crs=crs)

# Create a Folium Map
m = folium.Map(location=[6.5244, 2.6844], zoom_start=8)

# Add the polygons to the map with only boundaries (without filling)
for _, row in gdf_polygons.iterrows():
    if row['Matched_Polygon']:
        geo_json = folium.GeoJson(row['Matched_Polygon'], style_function=lambda x: {'fillColor': 'transparent', 'color': 'darksalmom', 'weight': 2})
        m.add_child(geo_json)

# Add markers for matched farmers
matched_farmers = key_farmers_data_with_matching_status[key_farmers_data_with_matching_status['Matching_Status'] == 1]
for idx, row in matched_farmers.iterrows():
    folium.Circle(
        location=[row['GPS_Latitude'], row['GPS_Long']],
        radius=50,
        color='#191970',
        fill=True,
        fill_color='#191970'
    ).add_to(m)

# Add a custom scale bar using HTML
m.get_root().html.add_child(folium.Element('<div style="position:fixed; bottom:20px; left:20px; background:white; z-index:1000; border:2px solid grey; padding:5px">Scale: 0 to 100 km</div>'))

# Save the map as HTML
m.save('map_with_matched_farmers_and_polygons.html')

# Open the HTML map in a browser using Selenium
browser = webdriver.Chrome()  # Use the relevant driver based on your browser
browser.get('file:///path_to/map_with_matched_farmers_and_polygons.html')


import plotly.graph_objects as go
from shapely.wkt import loads

# Columns in your dataset (replace these with your columns)
wkt_column = 'Matched_Polygon_WKT'
latitude_column = 'GPS_Latitude'
longitude_column = 'GPS_Long'

# Convert WKT to Shapely geometries using Shapely's loads method
key_farmers_data_with_matching_status['Matched_Polygon'] = key_farmers_data_with_matching_status[wkt_column].apply(lambda x: loads(x) if x else None)

# # Plot polygons using Plotly
# fig = go.Figure()
#
# for polygon in key_farmers_data_with_matching_status['Matched_Polygon']:
#     if polygon:
#         x, y = polygon.exterior.xy
#         x = list(x)  # Convert x coordinates to a list
#         y = list(y)  # Convert y coordinates to a list
#         fig.add_trace(go.Scatter(x=x, y=y, mode='lines', line=dict(color='red', width=2)))
#
# # Plot points on the map
# fig.add_trace(go.Scatter(x=key_farmers_data_with_matching_status[longitude_column], y=key_farmers_data_with_matching_status[latitude_column], mode='markers', marker=dict(size=10, color='blue')))
#
# # Update layout
# fig.update_layout(
#     mapbox=dict(
#         style="carto-positron",
#         zoom=8,
#         center=dict(lat=key_farmers_data_with_matching_status[latitude_column].mean(), lon=key_farmers_data_with_matching_status[longitude_column].mean())  # Centers the map based on the mean of latitudes and longitudes
#     )
# )
#
# fig.show()

import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Assuming 'key_farmers_data_with_matching_status' contains geographical data

# Create a GeoDataFrame with the matched polygons
# key_farmers_data_with_matching_status['Matched_Polygon'] = ...

# Create a figure and axes with Cartopy projection
fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))

# Plot country boundaries
ax.add_feature(cfeature.BORDERS, linestyle='-', edgecolor='black')

# Assuming 'key_farmers_data_with_matching_status' contains city names and coordinates
# Plot cities/villages as points
for idx, row in key_farmers_data_with_matching_status.iterrows():
    ax.plot(row['GPS_Long'], row['GPS_Latitude'], marker='o', markersize=5, color='red', transform=ccrs.PlateCarree())

# Add legend and labels for cities
plt.legend()
for idx, row in key_farmers_data_with_matching_status.iterrows():
    ax.text(row['GPS_Long'], row['GPS_Latitude'], None, transform=ccrs.PlateCarree(), ha='right')

# Set title and show plot
plt.title('Map with Countries, and Matched Polygons')
plt.show()