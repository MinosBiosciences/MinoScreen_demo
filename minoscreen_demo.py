import os
import sys
import glob
from collections import OrderedDict
import numpy as np
import pandas as pd
import dash
from dash import html
from dash import dcc
from dash.dash_table.Format import Format, Scheme
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from dash.dependencies import Input, Output

# Setting paths and specific module imports
# from config import resources_repo_path, title
# sys.path.insert(0, f"{resources_repo_path}Scripts")
#from global_variables import samples, bigbi_path, human_color, mouse_color, grey_color, cell_type_colors, ggyor_colors, Minos_dark_color, Minos_light_color

# Global_variables variables
Minos_dark_color = "#1e8ba5"
mouse_color = "#ff00ff"
human_color = "#b1d700"
grey_color = "#abb2b9"
cell_type_colors = {"Human": human_color,
                    "H_Ramos": human_color,
                    "human cell": human_color,  # TO BE REMOVED
                    "Mouse": mouse_color,
                    "M_B3Z": mouse_color,
                    "mouse cell": mouse_color,  # TO BE REMOVED
                    "Undetermined": grey_color,
                    "undetermined": "#a9a9a9",  # TO BE REMOVED
                    "Undetermined_cells": "#a9a9a9"}  # TO BE REMOVED
bigbi_path = "assets/Minos_space/"


def get_single_cells(chip_ID: str):
    try:
        # Load ref cells file
        cells_file = os.path.join(f"{bigbi_path}Chips/Chip{chip_ID}/{chip_ID}_data/{chip_ID}.20X.T0.cells.tsv")
        df_cells = pd.read_table(cells_file)

        # Load ref cages file
        cages_file = os.path.join(f"{bigbi_path}Chips/Chip{chip_ID}/{chip_ID}_data/{chip_ID}.20X.T0.cages.tsv")
        df_cages = pd.read_table(cages_file)

        # Filter cells file to only include rows where Cell_Index_Global equals Unique_Cell (ref cell version for redundant ones)
        df_cells_filtered = df_cells[df_cells['Cell_Index_Global'] == df_cells['Unique_Cell']]

        # Filter cells file further to only include rows where the Cage_ID in the cage file has Cell_Count_Cage with a value of 1 (single-cell cages)
        single_cell_cages = df_cages[df_cages['Cell_Count_Cage'] == 1]['Cage_ID'].tolist()
        df_cells_filtered = df_cells_filtered[df_cells_filtered['Cage_ID'].isin(single_cell_cages)]

    except Exception as e:
        print(f"Error loading or processing files: {e}")
        return {}
    return df_cells_filtered

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

server = app.server

# Hardcoded samples
samples_list = ["Chip204_M049-PoC79", "Chip206_M050-PoC80"]
selected_sample = samples_list[0]  # Default to first sample

# Creating dropdown options
samples_dropdown = []
for sample in samples_list:
    samples_dropdown.append({"label": sample, "value": sample})

#############################################################
#############################################################
#                                                           #
#                         LAYOUT                            #
#                                                           #
#############################################################
#############################################################

app.layout = html.Div([
    dbc.Row(dbc.Col(html.H1(id="H1", children=f"MinoScreen: Minos coupled datasets explorer", className="app-header--title"))),

    # Sample selector (dropdown)
    dbc.Row(dbc.Col(dcc.Dropdown(id="dropdown", style={"width": "330px"}, options=samples_dropdown, value=selected_sample), width=4), justify="center", className="mb-4"),

    # First row with table, pie chart, and images
    dbc.Row([
        # Single table on the left (template)
        dbc.Col(html.Div([
            html.H3("Experiment summary"),
            html.Div(id="sample_data_table", className="info_table")
        ]), width=4),

        # Two cell crops grids side by side on the right
        dbc.Col(html.Div([
            html.H3("Example cell crops"),
            html.Div(id="cells_per_cage_spatial")]), width=8)
    ]),

    # Second row
    dbc.Row([
        # Cell tye pie chart on the left
        dbc.Col(html.Div([
            html.H3("Cell Type Distribution"), dcc.Graph(id="cell_types_piechart")
        ]), width=3),

        # Fluo pairplot on the right
        dbc.Col(html.Div([
            html.H3("Fluo intensity pairplot (log values)"), dcc.Graph(id="cell_fluo_intensity_pairplot")
        ]), width=9),
    ])
])


#############################################################
#############################################################
#                                                           #
#                       CALLBACKS                           #
#                                                           #
#############################################################
#############################################################


# Sample data table callback
@app.callback(Output("sample_data_table", "children"),
              [Input("dropdown", "value")])
def update_sample_table(selected_sample):
    # Create a template table with the specified structure
    # This is just a placeholder - you'll fill this with actual data
    df = pd.DataFrame({
        selected_sample: ['Cage count', 'Global UMI count', 'Median cage UMI count', 'Global unique gene count', 'Median cage unique gene count'],
        'All single-cell cages': ['', '', '', '', ''],
        'Human single-cell cages': ['', '', '', '', ''],
        'Murine single-cell cages': ['', '', '', '', '']
    })

    return html.Div(
        [dash.dash_table.DataTable(
            df.to_dict('records'),
            [{"name": i, "id": i} for i in df.columns],
            style_cell={
                'padding': "5px", 
                "font-family": "Avenir",
                'maxWidth': '100px',  # Set maximum width for all columns
                'overflow': 'hidden',
                'textOverflow': 'ellipsis',
                'whiteSpace': 'normal'  # Enable text wrapping
            },
            style_header={"fontWeight": "bold"},
            style_data_conditional=[
                {'if': {'row_index': 'odd'},
                 'backgroundColor': 'rgb(245, 245, 245)'}
            ],
            # Style the first column as a header
            style_cell_conditional=[
                {'if': {'column_id': selected_sample},
                 'fontWeight': 'bold',
                 'textAlign': 'left',
                 'width': '120px'  # Set specific width for the first column
                },
                {'if': {'column_id': 'All single-cell cages'},
                 'width': '100px'  # Set specific width for this column
                },
                {'if': {'column_id': 'Human single-cell cages'},
                 'width': '100px'  # Set specific width for this column
                },
                {'if': {'column_id': 'Murine single-cell cages'},
                 'width': '100px'  # Set specific width for this column
                }
            ]
        )]
    )

# Cell crops grids callback
@app.callback(Output("cells_per_cage_spatial", "children"),
              [Input("dropdown", "value")])
def update_cells_per_cage(selected_sample):
    parts = selected_sample.replace("Chip", "").split("_")
    chip_ID = parts[0]
    grid_path_Cy5 = glob.glob(f"assets/Minos_space/Chips/Chip{chip_ID}/{chip_ID}_images/20x/Grids/T0/{chip_ID}.10x10_grid.cells.species.Grid*FITC-_Cy5+.RGB.best.tif")[0]
    grid_path_FITC = glob.glob(f"assets/Minos_space/Chips/Chip{chip_ID}/{chip_ID}_images/20x/Grids/T0/{chip_ID}.10x10_grid.cells.species.Grid*FITC+_Cy5-.RGB.best.tif")[0]
    return [html.Div([html.Img(id="image1", src=app.get_asset_url(grid_path_Cy5.lstrip("assets/")), alt="mesc", style={"width": "100%", "height": "auto", "border": "1px solid #ddd"}),
                      html.P("Human cells (Cy5 marker)")], style={"width": "48%", "display": "inline-block", "marginRight": "4%"}),
            html.Div([html.Img(id="image2", src=app.get_asset_url(grid_path_FITC.lstrip("assets/")), style={"width": "100%", "height": "auto", "border": "1px solid #ddd"}),
                      html.P("Murine cells (FITC marker)")], style={"width": "48%", "display": "inline-block"})
            ]

# Cell type pie chart callback
@app.callback(Output("cell_types_piechart", "figure"),
              [Input("dropdown", "value")])
def update_cell_pie_chart(selected_sample):
    # Extract chip_ID and seqID from the format "Chip204_M049-PoC79"
    parts = selected_sample.replace("Chip", "").split("_")
    chip_ID = parts[0]
    seqID = parts[1].split("-")[0]  # Extract M049 from M049-PoC79

    # Keep onyl single cells in cages
    df_sc = get_single_cells(chip_ID)

    # Computing the cell type count for each type
    type_list = ["H_Ramos", "M_B3Z", "Undetermined"]
    cell_type_count = [len(df_sc.query("Cell_Type==@kind")) for kind in type_list]
    fig = go.Figure(data=[go.Pie(values=cell_type_count, labels=type_list, textinfo='label+percent', hole=0.5)])
    fig.update_traces(marker=dict(colors=type_list, line=dict(color="black", width=2)))

    fig.update_layout(
        # title_text=f"Chip {chip_ID}: Cell type repartition",
        legend_title='Cell type',
        title_font=dict(family='Avenir', size=24, color=Minos_dark_color),
        title_x=0.5,
        font_family="Avenir",
        font_color="black",
        legend_title_font_color="black",
        height=500,
        margin=dict(
            l=50,
            r=50,
            b=100,
            t=100,
            pad=40
        ),
        legend=dict(
            orientation="v",
            yanchor="bottom",
            y=-0.9,
            xanchor="right",
            x=0.6,
        ),
        font=dict(
            size=15)
    )

    return fig


@app.callback(Output("cell_fluo_intensity_pairplot", "figure"),
              [Input("dropdown", component_property="value")])
def update_cell_fluo_pairplot(dropdown_value):
    # Extract chip_ID from the selected sample value
    parts = selected_sample.replace("Chip", "").split("_")
    chip_ID = parts[0]

    # Keep onyl single cells in cages
    df_sc = get_single_cells(chip_ID)

    # Checking whether cell types were determined
    if "Cell_Type" not in df_sc.columns:
        return {}

    # Column indices for fluo intensity
    intensity_cols = [col_name for col_name in df_sc.columns if col_name.startswith("Median_intensity_") and not df_sc[col_name].isna().all()]
    # log transformation
    df_sc.loc[:, intensity_cols] = np.log10(df_sc.loc[:, intensity_cols])
    # Get all cols
    all_cols = intensity_cols + ["Cage_ID", "Cell_Type", "X_Cell_Center_Global", "Y_Cell_Center_Global", "X_Cell_Center_Local", "Y_Cell_Center_Local"]
    fig = px.scatter_matrix(df_sc[all_cols], opacity=0.45, dimensions=intensity_cols, color="Cell_Type", color_discrete_map=cell_type_colors)

    fig.update_traces(showupperhalf=False)
    fig.update_traces(selector=dict(mode='markers', marker=dict(line=dict(color="black", width=2))))

    fig.update_layout(
        legend_title='Cell type',
        title_font=dict(family='Avenir', size=24, color=Minos_dark_color),
        title_x=0.5,
        font_family="Avenir",
        font_color="black",
        legend_title_font_color="black",
        height=900,
        legend=dict(
            yanchor="top",
            y=0.5,
            xanchor="right",
            x=1.3,
        ),
        font=dict(
            size=15)
    )
    return fig


if __name__ == '__main__':
    import argparse  # Fixed a bug with gunicorn

    # Managing command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", dest="port", type=int, default=8050, help="Port number to display on.")
    args = parser.parse_args()
    # args, unknown = parser.parse_known_args()  # Fixed a bug with gunicorn

    app.run(debug=True, port=args.port)
