import pandas as pd
from shapely.geometry import Point
from shapely.wkt import loads



farmers_data = pd.read_stata("ben_civ_rice_lsms_data_with_polygons_clusters_.dta")  # Replace with your file name



# Load Key farmers dataset
key_farmers_data = pd.read_stata("key_farmers_data_after_5km_10km_matching_use.dta")  # Replace with your file name

def calculate_distance(row):
    column_names = ['WKT2_MP_Any_tech_Adoption_5km', 'WKT2_MP_RiceAdvice_Adoption_5km', 'WKT2_MP_SRI_Adoption_5km', 'WKT2_MP_SmartValley_Adoption_5km', 'WKT2_MP_Any_tech_Adoption_10km', 'WKT2_MP_RiceAdvice_Adoption_10km', 'WKT2_MP_SRI_Adoption_10km', 'WKT2_MP_SmartValley_Adoption_10k']
    # Access GPS coordinates
    farmer_coords = Point(row['GPS_Long'], row['GPS_Latitude'])

    distances = []
    for column_name in column_names:
        # Access the polygon from WKT
        polygon_wkt = row[column_name]
        if pd.notnull(polygon_wkt):  # Check if WKT data exists
            try:
                polygon = loads(polygon_wkt)
                polygon_center = polygon.centroid

                # Calculate and append the distance
                distance_degree = farmer_coords.distance(polygon_center)
                distance_km = distance_degree * 111.32  # Convert Haversine distance to kilometers
                distances.append(distance_km)
            except Exception as e:
                print(f"Error processing WKT in column '{column_name}' for row {row.name}: {e}")
                distances.append(None)  # Append None for problematic WKT
        else:
            distances.append(None)  # Append None for missing or invalid WKT data

    return distances

# Apply the function to calculate distances in kilometers and divide by 111.32
key_farmers_data[['Distance_Any_tech_5km', 'Distance_RiceAdvice_5km', 'Distance_SRI_5km', 'Distance_SmartValley_5km',
                   'Distance_Any_tech_10km', 'Distance_RiceAdvice_10km', 'Distance_SRI_10km', 'Distance_SmartValley_10km']] = \
    pd.DataFrame(key_farmers_data.apply(calculate_distance, axis=1).to_list(), index=key_farmers_data.index)
# Apply the function to calculate distances and create new columns
key_farmers_data[['Distance_Any_tech_5km', 'Distance_RiceAdvice_5km', 'Distance_SRI_5km', 'Distance_SmartValley_5km',
                   'Distance_Any_tech_10km', 'Distance_RiceAdvice_10km', 'Distance_SRI_10km', 'Distance_SmartValley_10km']] = \
    pd.DataFrame(key_farmers_data.apply(calculate_distance, axis=1).to_list(), index=key_farmers_data.index)

print(key_farmers_data[['Codes', 'Distance_Any_tech_5km', 'Distance_RiceAdvice_5km', 'Distance_SRI_5km', 'Distance_SmartValley_5km',
                        'Distance_Any_tech_10km', 'Distance_RiceAdvice_10km', 'Distance_SRI_10km', 'Distance_SmartValley_10km']])

# Save the updated dataset
key_farmers_data.to_stata('key_farmers_data_after_5km_10km_matching_use_distance.dta', write_index=False, version=118)







#
#
# # Load Key farmers dataset
# key_farmers_data = pd.read_stata("key_farmers_data_after_5km_10km_matching_use.dta")  # Replace with your file name
#
# def calculate_distance(row):
#     column_names = ['WKT2_MP_Any_tech_Adoption_5km', 'WKT2_MP_RiceAdvice_Adoption_5km', 'WKT2_MP_SRI_Adoption_5km', 'WKT2_MP_SmartValley_Adoption_5km', 'WKT2_MP_Any_tech_Adoption_10km', 'WKT2_MP_RiceAdvice_Adoption_10km', 'WKT2_MP_SRI_Adoption_10km', 'WKT2_MP_SmartValley_Adoption_10k']
#     # Access GPS coordinates
#     farmer_coords = Point(row['GPS_Long'], row['GPS_Latitude'])
#
#     distances = []
#     for column_name in column_names:
#         # Access the polygon from WKT
#         polygon_wkt = row[column_name]
#         if pd.notnull(polygon_wkt):  # Check if WKT data exists
#             try:
#                 polygon = loads(polygon_wkt)
#                 polygon_center = polygon.centroid
#
#                 # Calculate and append the distance
#                 distance_degree = farmer_coords.distance(polygon_center)
#                 distance_km = distance_degree * 111.32  # Convert Haversine distance to kilometers
#                 distances.append(distance_km)
#             except Exception as e:
#                 print(f"Error processing WKT in column '{column_name}' for row {row.name}: {e}")
#                 distances.append(None)  # Append None for problematic WKT
#         else:
#             distances.append(None)  # Append None for missing or invalid WKT data
#
#     return distances
#
# # Apply the function to calculate distances in kilometers and divide by 111.32
# key_farmers_data[['Distance_Any_tech_5km', 'Distance_RiceAdvice_5km', 'Distance_SRI_5km', 'Distance_SmartValley_5km',
#                    'Distance_Any_tech_10km', 'Distance_RiceAdvice_10km', 'Distance_SRI_10km', 'Distance_SmartValley_10km']] = \
#     pd.DataFrame(key_farmers_data.apply(calculate_distance, axis=1).to_list(), index=key_farmers_data.index)
# # Apply the function to calculate distances and create new columns
# key_farmers_data[['Distance_Any_tech_5km', 'Distance_RiceAdvice_5km', 'Distance_SRI_5km', 'Distance_SmartValley_5km',
#                    'Distance_Any_tech_10km', 'Distance_RiceAdvice_10km', 'Distance_SRI_10km', 'Distance_SmartValley_10km']] = \
#     pd.DataFrame(key_farmers_data.apply(calculate_distance, axis=1).to_list(), index=key_farmers_data.index)
#
# print(key_farmers_data[['Codes', 'Distance_Any_tech_5km', 'Distance_RiceAdvice_5km', 'Distance_SRI_5km', 'Distance_SmartValley_5km',
#                         'Distance_Any_tech_10km', 'Distance_RiceAdvice_10km', 'Distance_SRI_10km', 'Distance_SmartValley_10km']])
#
# # Save the updated dataset
# key_farmers_data.to_stata('key_farmers_data_after_5km_10km_matching_use_distance.dta', write_index=False, version=118)
#












# import pandas as pd
# from shapely.geometry import Point
# from shapely.wkt import loads
#
# # Load Key farmers dataset
# key_farmers_data = pd.read_stata("key_farmers_data_after_5km_10km_matching_use.dta")  # Replace with your file name
#
# # Columns to consider
# column_names = ['WKT2_MP_Any_tech_Adoption_5km', 'WKT2_MP_RiceAdvice_Adoption_5km', 'WKT2_MP_SRI_Adoption_5km', 'WKT2_MP_SmartValley_Adoption_5km', 'WKT2_MP_Any_tech_Adoption_10km', 'WKT2_MP_RiceAdvice_Adoption_10km', 'WKT2_MP_SRI_Adoption_10km', 'WKT2_MP_SmartValley_Adoption_10k']
# ms_columns = ['MS_Any_tech_Adoption_5km', 'MS_RiceAdvice_Adoption_5km', 'MS_SRI_Adoption_5km', 'MS_SmartValley_Adoption_5km', 'MS_Any_tech_Adoption_10km', 'MS_RiceAdvice_Adoption_10km', 'MS_SRI_Adoption_10km', 'MS_SmartValley_Adoption_10km']
# # Function to calculate distance
# def calculate_distance(row):
#     distances = []
#     for column_name in column_names:
#         # Check if the corresponding 'MS' column is 1
#         ms_column_name = column_name.replace('WKT2_MP_', 'MS_')
#         if row[ms_column_name] == 1:
#             # Access GPS coordinates
#             farmer_coords = Point(row['GPS_Long'], row['GPS_Latitude'])
#
#             # Access the polygon from WKT
#             polygon_wkt = row[column_name]
#             if pd.notnull(polygon_wkt):  # Check if WKT data exists
#                 polygon = loads(polygon_wkt)
#                 polygon_center = polygon.centroid
#
#                 # Calculate and append the distance
#                 distances.append(farmer_coords.distance(polygon_center))
#             else:
#                 distances.append(None)  # Append None for missing or invalid WKT data
#         else:
#             distances.append(None)  # Append None for rows where 'MS' column is not 1
#
#     return distances
#
# # Apply the function to calculate distances and create new columns
# key_farmers_data[['Distance_Any_tech_5km', 'Distance_RiceAdvice_5km', 'Distance_SRI_5km', 'Distance_SmartValley_5km',
#                    'Distance_Any_tech_10km', 'Distance_RiceAdvice_10km', 'Distance_SRI_10km', 'Distance_SmartValley_10km']] = \
#     pd.DataFrame(key_farmers_data.apply(calculate_distance, axis=1).to_list(), index=key_farmers_data.index)
#
# print(key_farmers_data[['Codes', 'Distance_Any_tech_5km', 'Distance_RiceAdvice_5km', 'Distance_SRI_5km', 'Distance_SmartValley_5km',
#                         'Distance_Any_tech_10km', 'Distance_RiceAdvice_10km', 'Distance_SRI_10km', 'Distance_SmartValley_10km']])
#
# # Save the updated dataset
# key_farmers_data.to_stata('key_farmers_data_after_5km_10km_matching_use_distance.dta', write_index=False, version=118)
#
#
