import os
import json
import ast
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from shapely.geometry import shape


output_dir = (
    "/Users/apple/MyDocs/Mydocuments/"
    "GSAResEcon(RDFlab)/PhDProgram2022_2025/"
    "Datasets_4/New regressions Digital Ability"
)
os.makedirs(output_dir, exist_ok=True)


# LSMS data (Stata)
lsms_data = pd.read_stata("DA_Tech_Adopt_regressions.dta")

# Benin shapefiles
ben0 = gpd.read_file("ben_shapefile/ben_admbnda_adm0_1m_salb_20190816.shp")
ben1 = gpd.read_file("ben_shapefile/ben_admbnda_adm1_1m_salb_20190816.shp")
ben2 = gpd.read_file("ben_shapefile/ben_admbnda_adm2_1m_salb_20190816.shp")


agg = (
    lsms_data
        .groupby("region_string")["Dig_Ability_norm"]
        .mean()
        .reset_index()
)
# normalize names to lowercase and fix known mismatches
agg["region_string"] = agg["region_string"].str.lower().replace({"atacora": "atakora"})
ben2["adm1_name"] = ben2["adm1_name"].str.lower()

# merge into adm2 shapefile (which also has adm1_name column)
ben2 = ben2.merge(
    agg,
    left_on="adm1_name",
    right_on="region_string",
    how="left"
)
# fill missing with sentinel
ben2["Dig_Ability_norm"] = ben2["Dig_Ability_norm"].fillna(-1)


colors = ["#dadade", "#a8a8a8", "#99a3ae", "#acbcc4", "#43596b"]
# get maximum of real data (excluding -1)
max_val = ben2.loc[ben2["Dig_Ability_norm"] != -1, "Dig_Ability_norm"].max()
# create bin boundaries: first is 0, then equally spaced up to max_val
boundaries = [0] + list(np.linspace(0, max_val, len(colors) - 1))
# extended: include -1 for "no data"
extended_cmap = mcolors.ListedColormap(["white"] + colors)
extended_boundaries = [-1] + boundaries
extended_norm = mcolors.BoundaryNorm(extended_boundaries, extended_cmap.N)


tech_categories = {
    "EA_Anytech_adopt": "Technology (Anytech)",
    "EA_Vill_RiceAdvice_adopt": "RiceAdvice",
    "EA_Vill_SRI_adopt": "SRI",
    "EA_Vill_SmartValley_adopt": "SmartValley",
}
tech_colors = {
    "EA_Anytech_adopt": "#34AA1F",
    "EA_Vill_RiceAdvice_adopt": "#34AA1F",
    "EA_Vill_SRI_adopt": "#34AA1F",
    "EA_Vill_SmartValley_adopt": "#34AA1F",
}



def parse_geom(obj):
    """
    Turn a dict or a stringified dict into a Shapely geometry.
    """
    if isinstance(obj, dict):
        geom_dict = obj
    elif isinstance(obj, str):
        # try JSON first
        try:
            geom_dict = json.loads(obj)
        except ValueError:
            # fallback to Python literal (single quotes, tuples, etc.)
            geom_dict = ast.literal_eval(obj)
    else:
        return None

    return shape(geom_dict)


# apply parser to build geometry column
city_geoms = lsms_data["City_Polygon_JSON"].apply(parse_geom)
city_gdf = gpd.GeoDataFrame(lsms_data, geometry=city_geoms, crs="EPSG:4326")
# reproject to match ben2 CRS
city_gdf = city_gdf.to_crs(ben2.crs)


for tech_col, tech_label in tech_categories.items():
    fig, ax = plt.subplots(figsize=(12, 10))

    # draw country and admin boundaries
    ben0.boundary.plot(ax=ax, linewidth=1, color="black")
    ben1.boundary.plot(ax=ax, linewidth=1, color="black")
    # choropleth of digital ability
    ben2.plot(
        column="Dig_Ability_norm",
        ax=ax,
        cmap=extended_cmap,
        norm=extended_norm,
        linewidth=0.8,
        edgecolor="0.8",
        legend=False
    )

    # mask treatment vs control
    is_treat = city_gdf[tech_col] == 1
    is_control = city_gdf[tech_col] == 0

    # plot treatment footprints (filled)
    city_gdf[is_treat].plot(
        ax=ax,
        facecolor=tech_colors[tech_col],
        edgecolor="black",
        linewidth=0.5,
        alpha=0.6,
        label=f"Treatment – {tech_label}"
    )
    # plot control footprints (hatches)
    city_gdf[is_control].plot(
        ax=ax,
        facecolor="none",
        edgecolor="#C52700",
        linewidth=1.0,
        hatch="///",
        label=f"Control – {tech_label}"
    )

    # title
    ax.set_title(f"Adopters’ Enumeration Area – {tech_label}", fontsize=14)

    # colorbar
    sm = plt.cm.ScalarMappable(cmap=extended_cmap, norm=extended_norm)
    sm.set_array([])
    cbar = plt.colorbar(
        sm,
        ax=ax,
        orientation="horizontal",
        fraction=0.03,
        pad=0.10
    )
    cbar.set_label("Aggregated Digital Ability (LSMS Rice farmers)", fontsize=10)
    # ticks & labels
    tick_locs = extended_boundaries
    tick_lbls = ["No data"] + [f"{v:.2f}" for v in np.linspace(0, max_val, len(colors))]
    cbar.set_ticks(tick_locs)
    cbar.ax.set_xticklabels(tick_lbls, fontsize=8)

    # legend outside
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        fontsize="small"
    )

    plt.tight_layout()
    plt.show()

    # —— Save figure ——
    fname = f"{tech_label}.png"
    save_path = os.path.join(output_dir, fname)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)  # free memory

    print(f"Saved: {save_path}")
