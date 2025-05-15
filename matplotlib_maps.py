

import geopandas as gpd
import matplotlib.pyplot as plt
import json
import numpy as np
import pandas as pd
from shapely.wkt import loads
from matplotlib.colors import Normalize
import pandas as pd

# key_farmers_dict = pd.read_json('key_farmers_dict.json')
# lsms_farmers_dict = pd.read_json('lsms_farmers_dict.json')



# Load shapefiles for Benin and Cote D'Ivoire
benin_shapefile_0 = gpd.read_file("/Users/apple/Downloads/gadm40_BEN_shp/gadm40_BEN_0.shp")
benin_shapefile_1 = gpd.read_file("/Users/apple/Downloads/gadm40_BEN_shp/gadm40_BEN_1.shp")
benin_shapefile_2 = gpd.read_file("/Users/apple/Downloads/gadm40_BEN_shp/gadm40_BEN_2.shp")

# cotedivoire_shapefile_0 = gpd.read_file("/Users/apple/Downloads/gadm40_CIV_shp/gadm40_CIV_0.shp")
# cotedivoire_shapefile_1 = gpd.read_file("/Users/apple/Downloads/gadm40_CIV_shp/gadm40_CIV_1.shp")
# cotedivoire_shapefile_2 = gpd.read_file("/Users/apple/Downloads/gadm40_CIV_shp/gadm40_CIV_2.shp")

# Load key farmers data with polygons
key_farmers_data = pd.read_stata("key_farmers_data_after_10km_matching.dta")
lsms_data = pd.read_stata("new_lsms_data_for_tech_adopt_readyforregressions_newcelltower.dta")

filtered_data = key_farmers_data[(key_farmers_data['Country'] == 'Bénin') | (key_farmers_data['Country'] == 'Côte d’Ivoire')]

key_farmers_data = filtered_data[filtered_data['Country'] == 'Bénin']

lsms_data['Matched_Polygon'] = lsms_data['Matched_Polygon_WKT'].apply(lambda x: loads(x) if x else None)
crs = 'EPSG:4326'
gdf_polygons = gpd.GeoDataFrame(lsms_data, geometry='Matched_Polygon', crs=crs)


tech_categories = ['Any_tech_Adoption']

# Create consistent colors for each technology
tech_colors = {
    'Any_tech_Adoption': '#4a9e28'
}

# Loop through each tech category
for tech in tech_categories:
    fig, ax = plt.subplots(figsize=(10, 8))
    benin_shapefile_0.plot(ax=ax, color='white', edgecolor='black')
    benin_shapefile_1.boundary.plot(ax=ax, linewidth=1.5)
    benin_shapefile_2.boundary.plot(ax=ax, linewidth=1)

    # Plot WKT polygons and add them to the legend
    gdf_polygons.plot(ax=ax, color='none', edgecolor='#B70409', linewidth=1, label='LSMS 10km Polygon Cluster')

    # Plot adopters from the current tech category
    Adopters = key_farmers_data[(key_farmers_data[tech] == 1)]
    sc = ax.scatter(
        x=Adopters['GPS_Long'],
        y=Adopters['GPS_Latitude'],
        label=f'{tech} (Census)',
        c=tech_colors[tech],
        s=20,
        marker='.'  # Use dot marker
    )

    # Create a custom legend
    handles, labels = ax.get_legend_handles_labels()
    polygon_patch = plt.Line2D([0], [0], color='#B70409', linewidth=1)
    handles.append(polygon_patch)
    labels.append('LSMS 10km Polygon Cluster')
    ax.legend(handles=handles, labels=labels)

    plt.title(f'Census {tech} Adopters - Polygons (clusters)')
    plt.show()






# Load key farmers data with polygons
lsms_data = pd.read_stata("new_lsms_data_for_tech_adopt_readyforregressions_newcelltower.dta")
lsms_data['Matched_Polygon'] = lsms_data['WKT2_MP_Any_tech_Adoption_10km'].apply(lambda x: loads(x) if x else None)
crs = 'EPSG:4326'
gdf_polygons = gpd.GeoDataFrame(lsms_data, geometry='Matched_Polygon', crs=crs)


tech_categories_itt = ['Vill_Anytech_adopt']
# Create consistent colors for each technology
tech_colors = {
    'Vill_Anytech_adopt': '#4a9e28'
}

# Loop through each tech category
for tech in tech_categories_itt:
    fig, ax = plt.subplots(figsize=(10, 8))
    benin_shapefile_0.plot(ax=ax, color='white', edgecolor='black')
    benin_shapefile_1.boundary.plot(ax=ax, linewidth=1.5)
    benin_shapefile_2.boundary.plot(ax=ax, linewidth=1)
    # cotedivoire_shapefile_0.boundary.plot(ax=ax, color='lightgray', edgecolor='black')
    # cotedivoire_shapefile_1.boundary.plot(ax=ax, linewidth=1)
    # cotedivoire_shapefile_2.boundary.plot(ax=ax, linewidth=1)


    # Plot adopters from the current tech category
    treated_farmers = lsms_data[(lsms_data[tech] == 1)]
    control_farmers = lsms_data[(lsms_data[tech] == 0)]

    # Plot treated farmers with triangular shapes
    ax.scatter(
        x=treated_farmers['coordonnes_gps__Longitude'],
        y=treated_farmers['coordonnes_gps__Latitude'],
        c=tech_colors[tech],
        label=f'Treated - {tech}',
        s=20,
        marker='^'  # Use triangular marker
    )

    # Plot control farmers with triangular shapes
    ax.scatter(
        x=control_farmers['coordonnes_gps__Longitude'],
        y=control_farmers['coordonnes_gps__Latitude'],
        c='gray',
        label=f'Control - {tech}',
        s=20,
        marker='^'  # Use triangular marker
    )

    plt.title(f'{tech} Adopters - Polygons (clusters)')
    plt.legend()
    plt.show()


# import geopandas as gpd
# import matplotlib.pyplot as plt
# import pandas as pd
# from shapely.wkt import loads
#
# # Load shapefiles for Benin
# benin_shapefile_0 = gpd.read_file("/Users/apple/Downloads/gadm40_BEN_shp/gadm40_BEN_0.shp")
# benin_shapefile_1 = gpd.read_file("/Users/apple/Downloads/gadm40_BEN_shp/gadm40_BEN_1.shp")
# benin_shapefile_2 = gpd.read_file("/Users/apple/Downloads/gadm40_BEN_shp/gadm40_BEN_2.shp")
#
# # Load key farmers data with polygons
# key_farmers_data = pd.read_stata("key_farmers_data_after_10km_matching.dta")
# lsms_data = pd.read_stata("new_lsms_data_for_tech_adopt_readyforregressions_newcelltower.dta")
#
# filtered_data = key_farmers_data[(key_farmers_data['Country'] == 'Bénin') | (key_farmers_data['Country'] == 'Côte d’Ivoire')]
# key_farmers_data = filtered_data[filtered_data['Country'] == 'Bénin']
#
# lsms_data['Matched_Polygon'] = lsms_data['Matched_Polygon_WKT'].apply(lambda x: loads(x) if x else None)
# crs = 'EPSG:4326'
# gdf_polygons = gpd.GeoDataFrame(lsms_data, geometry='Matched_Polygon', crs=crs)
#
# tech_categories = ['Any_tech_Adoption']
#
# # Create consistent colors for each technology
# tech_colors = {
#     'Any_tech_Adoption': '#4a9e28'
# }
#
# # Loop through each tech category
# for tech in tech_categories:
#     fig, ax = plt.subplots(figsize=(10, 8))
#     benin_shapefile_0.plot(ax=ax, color='white', edgecolor='black')
#     benin_shapefile_1.boundary.plot(ax=ax, linewidth=1.5)
#     benin_shapefile_2.boundary.plot(ax=ax, linewidth=1)
#
#     # Plot WKT polygons and add them to the legend
#     gdf_polygons.plot(ax=ax, color='none', edgecolor='#13103d', linewidth=1, label='LSMS 10km Polygon Cluster')
#
#     # Plot adopters from the current tech category
#     Adopters = key_farmers_data[(key_farmers_data[tech] == 1)]
#     sc = ax.scatter(
#         x=Adopters['GPS_Long'],
#         y=Adopters['GPS_Latitude'],
#         label=f'{tech} (Census)',
#         c=tech_colors[tech],
#         s=20,
#         marker='.'  # Use dot marker
#     )
#
#     # Create a custom legend
#     handles, labels = ax.get_legend_handles_labels()
#     polygon_patch = plt.Line2D([0], [0], color='#13103d', linewidth=1)
#     handles.append(polygon_patch)
#     labels.append('LSMS 10km Polygon Cluster')
#     ax.legend(handles=handles, labels=labels)
#
#     plt.title(f'Census {tech} Adopters - Polygons (clusters)')
#     plt.show()


