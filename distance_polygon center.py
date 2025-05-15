import pandas as pd
from shapely.geometry import shape, Point
from shapely.wkt import loads

# Load Key farmers dataset
lsms_data = pd.read_stata("ben_civ_lsms_data_with_matching_status22newfinal.dta")  # Replace with your file name
lsms_data.dropna(subset=['coordonnes_gps__Longitude', 'coordonnes_gps__Latitude'], inplace=True)
# Function to calculate distance
def calculate_distance(row):
    # Access GPS coordinates
    farmer_coords = Point(row['coordonnes_gps__Longitude'], row['coordonnes_gps__Latitude'])

    # Access the polygon from WKT
    polygon_wkt = row['Matched_Polygon_WKT2']
    if pd.notnull(polygon_wkt):  # Check if WKT data exists
        try:
            polygon = loads(polygon_wkt)
            polygon_center = polygon.centroid

            # Calculate and return the distance
            return farmer_coords.distance(polygon_center)
        except Exception as e:
            print(f"Error parsing WKT for row {row.name}: {e}")
            return None
    else:
        return None  # Return None for missing or invalid WKT data

# Apply the function to calculate distances and create a new column
lsms_data['Distance_to_Adopters'] = lsms_data.apply(calculate_distance, axis=1)

print(lsms_data[['grappe', 'Distance_to_Adopters']])

# Save the updated dataset
lsms_data.to_stata('ben_civ_lsms_data_with_distance_Adopter.dta', write_index=False, version=118)
# import pandas as pd
# from shapely.geometry import shape, Point
# from shapely.wkt import loads
#
# # Load Key farmers dataset
# key_farmers_data = pd.read_stata("key_farmers_data_with_matching_status22.dta")  # Replace with your file name
#
#
# # Function to calculate distance
# def calculate_distance(row):
#     # Access GPS coordinates
#     farmer_coords = Point(row['GPS_Long'], row['GPS_Latitude'])
#
#     # Access the polygon from WKT
#     polygon_wkt = row['Matched_Polygon_WKT']
#     if pd.notnull(polygon_wkt):  # Check if WKT data exists
#         polygon = loads(polygon_wkt)
#         polygon_center = polygon.centroid
#
#         # Calculate and return the distance
#         return farmer_coords.distance(polygon_center)
#     else:
#         return None  # Return None for missing or invalid WKT data
#
# # Apply the function to calculate distances and create a new column
# key_farmers_data['Distance_to_Polygon_Center'] = key_farmers_data.apply(calculate_distance, axis=1)
#
# print(key_farmers_data[['Codes', 'Distance_to_Polygon_Center']])
#
# # Save the updated dataset
# key_farmers_data.to_stata('key_farmers_data_with_matching_status223.dta', write_index=False, version=118)

# # ... (Previous code remains unchanged) ...
#
# # Iterate through the matched farmers and calculate the distance to the matched polygon's center
# for index, info in match_info.items():
#     matched_polygon_center = info['matched_polygon'].centroid
#     farmer_coords = Point(key_farmers_data.loc[index, 'GPS_Long'], key_farmers_data.loc[index, 'GPS_Latitude'])
#
#     # Calculate distance and assign it to a new column
#     key_farmers_data.at[index, 'Distance_to_Matched_Polygon_Center'] = farmer_coords.distance(matched_polygon_center)
#
# # Save the updated dataset as CSV
# key_farmers_data.to_csv('key_farmers_data_with_distances.csv', index=False)
