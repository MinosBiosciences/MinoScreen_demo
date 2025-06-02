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
from dash.exceptions import PreventUpdate

# Setting paths and specific module imports
# from config import resources_repo_path, title
# sys.path.insert(0, f"{resources_repo_path}Scripts")
#from global_variables import samples, bigbi_path, human_color, mouse_color, grey_color, cell_type_colors, ggyor_colors, Minos_dark_color, Minos_light_color

# Global_variables variables
Minos_dark_color = "#1e8ba5"
# mouse_color = "#ff00ff"
mouse_color = "#47d348"  # Temp modif
# human_color = "#b1d700"
human_color = "#dc4142"  # Temp modif
grey_color = "#abb2b9"
cell_type_colors = {"Human": human_color,
                    "H_Ramos": human_color,
                    "Ramos_Human": human_color,
                    "human cell": human_color,  # TO BE REMOVED
                    "Mouse": mouse_color,
                    "M_B3Z": mouse_color,
                    "B3Z_murine": mouse_color,
                    "mouse cell": mouse_color,  # TO BE REMOVED
                    "Undetermined": grey_color,
                    "undetermined": "#a9a9a9",  # TO BE REMOVED
                    "Undetermined_cells": "#a9a9a9"}  # TO BE REMOVED
bigbi_path = "assets/Minos_space/"


def get_single_cells(chip_ID: str):
    try:
        # Load ref cells file
        cells_file = f"{bigbi_path}Chips/Chip{chip_ID}/{chip_ID}_data/{chip_ID}.20X.T0.cells.tsv"
        df_cells = pd.read_table(cells_file)

        # Load ref cages file
        cages_file = os.path.join(f"{bigbi_path}Chips/Chip{chip_ID}/{chip_ID}_data/{chip_ID}.20X.T0.cages.tsv")
        df_cages = pd.read_table(cages_file)

        # Filter cells file further to only include rows where the Cage_ID in the cage file has Cell_Count_Cage with a value of 1 (single-cell cages)
        single_cell_cages = df_cages[df_cages['Cell_Count_Cage'] == 1]['Cage_ID'].tolist()
        df_cells_filtered = df_cells[df_cells['Cage_ID'].isin(single_cell_cages)]

    except Exception as e:
        print(f"Error loading or processing files: {e}")
        return {}
    return df_cells_filtered

def get_cell_crop_path(chip_ID, cage_ID, mag, cells, channel):
    try:
        # print(chip_ID)
        cells_info = cells.loc[cage_ID]
        tile_ID = cells_info.loc['Tile_ID']
        if chip_ID == "159":
            cell_folder = f"Cell_crops2/{channel}"
        else:
            cell_folder = f"Cells/Crops/{channel}-T0"
            tile_ID = str(tile_ID).zfill(4)
        cell_path = f"Minos_space/Chips/Chip{chip_ID}/{chip_ID}_images/{mag.lower()}/{cell_folder}/chip{chip_ID}_tile{tile_ID}_cell{cells_info.loc['Cell_Index_Global']}_local{cells_info.loc['Cell_Index_Local']}.png"
        # print(cell_path)
    except:
        cell_path = ""
    return cell_path

def get_top_correlated_pairs(corr_matrix, mode, n, thresh, feat_type):
    # Get only selected feature type
    if feat_type == "genes_genes":
        corr_matrix = corr_matrix[(corr_matrix.Variable1_type == "Genes") & (corr_matrix.Variable2_type == "Genes")]
    elif feat_type == "img_img":
        corr_matrix = corr_matrix[(corr_matrix.Variable1_type.isin(["CP features"]) & corr_matrix.Variable2_type.isin(["DL features"])) | (corr_matrix.Variable1_type.isin(["DL features"]) & corr_matrix.Variable2_type.isin(["CP features"]))]
    elif feat_type == "genes_img":
        corr_matrix = corr_matrix[(corr_matrix.Variable1_type.isin(["CP features", "DL features"]) & corr_matrix.Variable2_type.isin(["Genes"])) | (corr_matrix.Variable1_type.isin(["Genes"]) & corr_matrix.Variable2_type.isin(["CP features", "DL features"]))]
    # Get the top correlated pairs
    if mode == "top_n":
        positive_pairs = corr_matrix[corr_matrix['Correlation'] > 0].nlargest(n, 'Correlation')
        negative_pairs = corr_matrix[corr_matrix['Correlation'] < 0].nsmallest(n, 'Correlation')
    elif mode == "threshold":
        positive_pairs = corr_matrix[(corr_matrix['Correlation'] > thresh)]
        negative_pairs = corr_matrix[(corr_matrix['Correlation'] < -thresh)]
    positive_pairs["CorrelationType"] = "Positive"
    negative_pairs["CorrelationType"] = "Negative"
    top_corr_pairs = pd.concat([positive_pairs, negative_pairs])
    top_corr_pairs = top_corr_pairs.reset_index()

    return top_corr_pairs.sort_values(by='Correlation', ascending=False)


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

server = app.server

# Hardcoded samples
samples_list = ["Chip204_M049-PoC78", "Chip206_M051-PoC80"]
selected_sample = samples_list[1]  # Default to first sample

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
    dbc.Row([
        dbc.Col(html.Label("Sample to analyse:", style={"margin-right": "10px", "display": "inline-block", "vertical-align": "middle", "fontWeight": "bold"}), width="auto"),
        dbc.Col(dcc.Dropdown(id="dropdown", style={"width": "330px", "vertical-align": "middle"}, options=samples_dropdown, value=selected_sample), width="auto")
    ], justify="center", className="mb-4", align="center"),

    # First row with table, pie chart, and images
    dbc.Row([
        # Single table on the left (template)
        dbc.Col(html.Div([
            html.H3("Experiment summary"),
            html.Div(id="sample_data_table", className="info_table")
        ]), width=4),

        # Two cell crops grids side by side on the right
        dbc.Col(html.Div([
            html.H3("Example cell crops from chip"),
            html.Div(id="cell_crops_grids")]), width=8)
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
    ]),

    # Third row
    dbc.Row([
        # Omics, bimodal and imaging clustering
        dbc.Col(html.Div([
            html.H3("Clustering"), dcc.Graph(id="clustering_sc", clear_on_unhover=True)
        ]))
    ]),

    # Fourth row
    dbc.Row([
        # Scatter plot UMI on the left
        dbc.Col(html.Div([
            html.H3("Cells omics vs imaging"), dcc.Graph(id="scatterplot_UMI_sc", clear_on_unhover=True),
        ]), width=3),

        # Image display on the center
        dbc.Col(html.Div([
            html.H3("Cells display", className="text-center"),
            html.Div(id="cell_display")
        ]), width=5),# justify="center", style={"height": "120px"},

        # Features pairplot on the right
        dbc.Col(html.Div([
            html.H3("Feat-feat plot"),
            html.Div(dcc.RadioItems(id='features-option-radio',
                                            options=[
                                                {'label': ' Genes vs Genes', 'value': 'genes_genes'},
                                                {'label': ' Genes vs imaging features', 'value': 'genes_img'},
                                                {'label': ' Imaging features: classic vs DL', 'value': 'img_img'}
                                            ],
                                            value='genes_genes',
                                            labelStyle={'display': 'inline-block', 'margin-right': '10px'}
                                        ),
                                 className='d-flex justify-content-center'
                         ),
                          html.Div(
                              dcc.RadioItems(id='filter-option-radio',
                                            options=[
                                                {'label': ' Top N', 'value': 'top_n'},
                                                {'label': ' Threshold', 'value': 'threshold'}
                                            ],
                                            value='top_n',
                                            labelStyle={'display': 'inline-block', 'margin-right': '10px'}
                                        ),
                                  className='d-flex justify-content-center'
                          ),
                          html.Div([
                              html.Div(
                                  dcc.Input(
                                      id='top-n-input',
                                      type='number',
                                      placeholder='Entrez la valeur de N',
                                      value=10,
                                      style={'margin-right': '10px'}
                                  ),
                                className='d-flex justify-content-center'
                              ),
                              html.Div([
                                  dcc.Slider(
                                      id='threshold-input',
                                      min=0.3,
                                      max=1,
                                      step=0.01,
                                      value=0.8,
                                      marks={i / 10: f'{i / 10:.1f}' for i in range(0, 11)},
                                      tooltip={"placement": "bottom", "always_visible": True}
                                  )
                              ], id='threshold-input-container', style={'display': 'none'})
                          ], id='filter-input-container'),
                         html.Br(),
        dcc.Dropdown(id='correlation-pairs-dropdown', placeholder="Select one of the most correlated or anti-correlated pair of variables", style={'font-size': '12px', 'width': '95%', 'margin': '0 auto'}),
        html.Div(dcc.Graph(id='scatter-plot', style={'width': '580px', 'height': '580px'}, clear_on_unhover=True), className='d-flex justify-content-center'),
        ]), width=4),
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
    # Extract chip_ID and seq_ID from the format the selected sample
    parts = selected_sample.replace("Chip", "").split("_")
    chip_ID = parts[0]
    seq_ID = parts[1]
    # Read the TSV file
    df = pd.read_csv(f"{bigbi_path}Chips/Chip{chip_ID}/{seq_ID}/{seq_ID}_data/{chip_ID}_{seq_ID}.MinoScreen_info.tsv", sep='\t', skiprows=[3, 5])  # For now we only show some info
    # Rename the first column (which has an empty or 'Unnamed: 0' header) to selected sample
    first_col_name = df.columns[0]  # Get the name of the first column
    df = df.rename(columns={first_col_name: selected_sample})

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
                 'backgroundColor': 'rgb(245, 245, 245)'},
                # Make first column bold
                {'if': {'column_id': selected_sample},
                 'fontWeight': 'bold'}
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
@app.callback(Output("cell_crops_grids", "children"),
              [Input("dropdown", "value")])
def update_cells_per_cage(selected_sample):
    parts = selected_sample.replace("Chip", "").split("_")
    chip_ID = parts[0]
    grid_path_Cy5 = glob.glob(f"{bigbi_path}Chips/Chip{chip_ID}/{chip_ID}_images/20x/Grids/T0/{chip_ID}.10x10_grid.cells.species.Grid*FITC-_Cy5+.RGB.best.png")[0]
    grid_path_FITC = glob.glob(f"{bigbi_path}Chips/Chip{chip_ID}/{chip_ID}_images/20x/Grids/T0/{chip_ID}.10x10_grid.cells.species.Grid*FITC+_Cy5-.RGB.best.png")[0]
    return [html.Div([html.Img(id="image1", src=app.get_asset_url(grid_path_Cy5.lstrip("assets/")), alt="mesc", style={"width": "100%", "height": "auto", "border": "1px solid #ddd"}),
                      html.P("Ramos cells (Human, Cy5 marker)")], style={"width": "48%", "display": "inline-block", "marginRight": "4%"}),
            html.Div([html.Img(id="image2", src=app.get_asset_url(grid_path_FITC.lstrip("assets/")), style={"width": "100%", "height": "auto", "border": "1px solid #ddd"}),
                      html.P("B3Z cells (murine, FITC marker)")], style={"width": "48%", "display": "inline-block"})
            ]

# Cell type pie chart callback
@app.callback(Output("cell_types_piechart", "figure"),
              [Input("dropdown", "value")])
def update_cell_pie_chart(selected_sample):
    # Extract chip_ID from the selected sample value
    parts = selected_sample.replace("Chip", "").split("_")
    chip_ID = parts[0]

    # Keep onyl single cells in cages
    df_sc = get_single_cells(chip_ID)
    df_sc.loc[df_sc.Cell_Type == "H_Ramos", "Cell_Type"] = "Ramos_Human"  # Temp modif
    df_sc.loc[df_sc.Cell_Type == "M_B3Z", "Cell_Type"] = "B3Z_murine"  # Temp modif

    # Computing the cell type count for each type
    # type_list = ["H_Ramos", "M_B3Z", "Undetermined"]
    type_list = ["Ramos_Human", "B3Z_murine", "Undetermined"]
    colors = [cell_type_colors[cell_type] for cell_type in type_list]
    cell_type_count = [len(df_sc.query("Cell_Type==@kind")) for kind in type_list]
    fig = go.Figure(data=[go.Pie(values=cell_type_count, labels=type_list, textinfo='label+percent', hole=0.5)])
    fig.update_traces(marker=dict(colors=colors, line=dict(color="black", width=2)))

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
def update_cell_fluo_pairplot(selected_sample):
    # Extract chip_ID from the selected sample value
    parts = selected_sample.replace("Chip", "").split("_")
    chip_ID = parts[0]

    # Keep onyl single cells in cages
    df_sc = get_single_cells(chip_ID)
    df_sc.loc[df_sc.Cell_Type == "H_Ramos", "Cell_Type"] = "Ramos_Human"  # Temp modif
    df_sc.loc[df_sc.Cell_Type == "M_B3Z", "Cell_Type"] = "B3Z_murine"  # Temp modif

    # Checking whether cell types were determined
    if "Cell_Type" not in df_sc.columns:
        return {}

    # Column indices for fluo intensity
    intensity_cols = [col_name for col_name in df_sc.columns if col_name.startswith("Median_intensity_") and not df_sc[col_name].isna().all()]
    # log transformation
    df_sc.loc[:, intensity_cols] = np.log10(df_sc.loc[:, intensity_cols])
    # Get all cols
    all_cols = intensity_cols + ["Cage_ID", "Cell_Type"]
    try:
        fig = px.scatter_matrix(df_sc[all_cols], opacity=0.45, dimensions=intensity_cols, color="Cell_Type", color_discrete_map=cell_type_colors)
    except ValueError:  # Hack to avoid a bug (https://community.plotly.com/t/valueerror-invalid-value-in-basedatatypes-py/55993/6)
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

# Omics, bimodal and imaging clustering callback
@app.callback(
    Output("clustering_sc", "figure"),
    Input("dropdown", "value")
)
def update_clustering(selected_sample):
    if selected_sample is None:
        # If no dropdown value or hover data is selected, clear the image.
        return None

    chip_ID, seq_ID = selected_sample.replace("Chip", "").split("_")

    df_clustering = pd.read_table(f"{bigbi_path}Chips/Chip{chip_ID}/{seq_ID}/{seq_ID}_data/{chip_ID}_{seq_ID}.clustering.tsv", index_col=0)
    df_clustering.loc[df_clustering.Cell_type_img == "Human", "Cell_type_img"] = "Ramos_Human"  # Temp modif
    df_clustering.loc[df_clustering.Cell_type_img == "Mouse", "Cell_type_img"] = "B3Z_murine"  # Temp modif
    df_clustering.loc[df_clustering.Cell_type_img == "Undetermined_cells", "Cell_type_img"] = "Undetermined"  # Temp modif
    df_clustering.loc[df_clustering.Cell_type_img == "Undetermined_cells", "Cell_type_img"] = "Undetermined"
    df_clustering["ID"] = df_clustering.index

    try:
        fig_omics = px.scatter(df_clustering, x="UMAP1_sc", y="UMAP2_sc", color="Cell_type_img", opacity=0.7, custom_data=["ID"],
                      color_discrete_map=cell_type_colors, category_orders={"Cell_type_img": ["Undetermined", "Mouse", "Human"]})
    except ValueError:
        fig_omics = px.scatter(df_clustering, x="UMAP1_sc", y="UMAP2_sc", color="Cell_type_img", opacity=0.7, custom_data=["ID"],
                               color_discrete_map=cell_type_colors, category_orders={"Cell_type_img": ["Undetermined", "Mouse", "Human"]})

    fig_omics.update_traces(hovertemplate="UMAP1: %{x}<br>UMAP2: %{y}<br>Cage_ID: %{customdata[0]}<extra></extra>")

    fig_bimodal = px.scatter(df_clustering, x="UMAP1_bimodal_sc", y="UMAP2_bimodal_sc", color="Cell_type_img", opacity=0.7, custom_data=["ID"],
                      color_discrete_map=cell_type_colors, category_orders={"Cell_type_img": ["Undetermined", "Mouse", "Human"]})
    fig_bimodal.update_traces(hovertemplate="UMAP1: %{x}<br>UMAP2: %{y}<br>Cage_ID: %{customdata[0]}<extra></extra>")

    fig_img = px.scatter(df_clustering, x="UMAP1_CP_sc", y="UMAP2_CP_sc", color="Cell_type_img", opacity=0.7, custom_data=["ID"],
                      color_discrete_map=cell_type_colors, category_orders={"Cell_type_img": ["Undetermined", "Mouse", "Human"]})
    fig_img.update_traces(hovertemplate="UMAP1: %{x}<br>UMAP2: %{y}<br>Cage_ID: %{customdata[0]}<extra></extra>")

    subplots_fig = make_subplots(rows=1, cols=3, subplot_titles=("Omics-based clustering", "Bimodal clustering", "Imaging-based clustering"),
                                 horizontal_spacing=0.08)

    for i in range(len(fig_omics["data"])):
        if i == 0:
            fig_omics['data'][i]["legendgrouptitle_text"] = "Cell type from imaging"
        subplots_fig.add_trace(fig_omics["data"][i], col=1, row=1)
    for i in range(len(fig_bimodal["data"])):
        if i == 0:
            fig_bimodal['data'][i]["legendgrouptitle_text"] = "Cell type from imaging"
        subplots_fig.add_trace(fig_bimodal["data"][i], col=2, row=1)
    for i in range(len(fig_img["data"])):
        if i == 0:
            fig_img['data'][i]["legendgrouptitle_text"] = "Cell type from imaging"
        subplots_fig.add_trace(fig_img["data"][i], col=3, row=1)

    for trace in subplots_fig.data:
        trace.update(marker_size=6, marker_line_color="white")

    subplots_fig.update_layout(#plot_bgcolor="white",
                               xaxis1_range=[df_clustering.UMAP1_sc.min() - 1, df_clustering.UMAP1_sc.max() + 1],
                               xaxis2_range=[df_clustering.UMAP1_bimodal_sc.min() - 1, df_clustering.UMAP1_bimodal_sc.max() + 1],
                               xaxis3_range=[df_clustering.UMAP1_CP_sc.min() - 1, df_clustering.UMAP1_CP_sc.max() + 1],
                               yaxis1_range=[df_clustering.UMAP2_sc.min() - 1, df_clustering.UMAP2_sc.max() + 1],
                               yaxis2_range=[df_clustering.UMAP2_bimodal_sc.min() - 1, df_clustering.UMAP2_bimodal_sc.max() + 1],
                               yaxis3_range=[df_clustering.UMAP2_CP_sc.min() - 1, df_clustering.UMAP2_CP_sc.max() + 1],
                               font_size=15,
                               margin=dict(
                                   b=30
                               ))

    subplots_fig.update_xaxes(showgrid=True, showticklabels=False)
    subplots_fig.update_yaxes(showgrid=True, showticklabels=False)

    names = set()
    subplots_fig.for_each_trace(
        lambda trace:
        trace.update(showlegend=False)
        if (trace.name in names) else names.add(trace.name))

    return subplots_fig

# Hover data highlight in the 3 clustering callback
"""
@app.callback(
    Output("clustering_sc", "figure", allow_duplicate=True),
    Input("clustering_sc", "figure"),
    Input("clustering_sc", "hoverData"),
    prevent_initial_call=True
)
def update_hover_point_pca_clustering(fig, hoverData):
    if not fig or "data" not in fig:
        return fig

    hovered_id = None
    if hoverData and "points" in hoverData and hoverData["points"]:
        try:
            hovered_id = hoverData['points'][0]['customdata'][0]
        except (KeyError, IndexError, TypeError):
            hovered_id = None

    for scat in fig["data"]:
        if "customdata" in scat and scat["customdata"]:
            ids = [x[0] for x in scat["customdata"]]

            # Index du point survolé (si hover actif)
            idx = -1
            if hovered_id in ids:
                idx = ids.index(hovered_id)

            # Réinitialiser tous les points
            size_array = [12 if j == idx else 6 for j in range(len(scat["x"]))]
            line_color_array = ['black' if j == idx else 'white' for j in range(len(scat["x"]))]

            # Mise à jour des propriétés
            scat.setdefault("marker", {})
            scat["marker"]["size"] = size_array
            scat["marker"].setdefault("line", {})
            scat["marker"]["line"]["color"] = line_color_array

    return fig
"""
# Cell crop display on hover clustering callback
@app.callback(
    Output("cell_display", "children", allow_duplicate=True),
    [Input('dropdown', 'value'),
     Input("clustering_sc", "hoverData")],
    prevent_initial_call=True
)
def cell_crop_on_hover_clustering(selected_sample, hoverData):
    """
    Display cell crop of sc cages.
    """
    if hoverData is None or selected_sample is None:
        # If no dropdown value or hover data is selected, clear the image.
        return None

    chip_ID, seq_ID = selected_sample.replace("Chip", "").split("_")
    mag="20X"

    channel_available = np.unique([c.split("-T")[0] for c in sorted(os.listdir(f"{bigbi_path}/Chips/Chip{chip_ID}/{chip_ID}_images/{mag.lower()}/Cells/Crops/"))])
    # Fix RGB image in first
    if "RGB" in channel_available:
        channel_available = np.delete(channel_available, np.where(channel_available == "RGB")[0])
        channel_available = np.insert(channel_available, 0, "RGB")

    cage_ID = hoverData["points"][0]["customdata"][0]
    df_cells_sc = get_single_cells(chip_ID)
    df_cells_sc.set_index("Cage_ID", inplace=True)

    cell_items = []
    for channel in channel_available:
        cell_path = get_cell_crop_path(chip_ID, cage_ID, mag, df_cells_sc, channel)
        # print(cell_path)
        if os.path.exists(cell_path.replace("Minos_space/", bigbi_path)):
            cell_items.append(html.Div([
                html.Div(channel, style={"text-align": "center", "font-size": "12px", "color": "#555"}),
                html.Img(
                    src=dash.get_asset_url(cell_path),
                    alt=f"Image not available (chip {chip_ID}, {seq_ID})",
                    style={"width": "70px", "margin-bottom": "5px"}
                )
            ], style={"margin-right": "10px", "margin-top": "0px", "display": "inline-block", "text-align": "center"}))

    try:
        if not cell_items:
            raise Exception
        return html.Div(cell_items, style={"display": "flex", "justify-content": "center", "align-items": "flex-start"})
    except Exception as error:
        print(error)
        raise PreventUpdate

# UMI scatterplot callback
@app.callback(Output("scatterplot_UMI_sc", "figure"),
              [Input("dropdown", "value")])
def update_scatterplot_UMI_sc(selected_sample):
    if selected_sample is None:
        # If no dropdown value or hover data is selected, clear the image.
        return None

    chip_ID, seq_ID = selected_sample.replace("Chip", "").split("_")

    df_scatter = pd.read_table(f"{bigbi_path}Chips/Chip{chip_ID}/{seq_ID}/{seq_ID}_data/{chip_ID}_{seq_ID}.scatter_UMI.tsv", index_col=0)
    maxi = max(max(df_scatter.H_UMIs_cage), max(df_scatter.M_UMIs_cage))
    df_scatter.loc[df_scatter.Cell_type_img == "Human", "Cell_type_img"] = "Ramos_Human"  # Temp modif
    df_scatter.loc[df_scatter.Cell_type_img == "Mouse", "Cell_type_img"] = "B3Z_murine"  # Temp modif
    df_scatter.loc[df_scatter.Cell_type_img == "Undetermined_cells", "Cell_type_img"] = "Undetermined"  # Temp modif

    fig = go.Figure()
    hovertemplate = "<b>CageID</b>: %{text}<br /><i>Ramos (Human) UMI count</i>: %{x}<br /><i>B3Z (murine) UMI count: %{y}"
    for spe in df_scatter.Cell_type_img.unique():
        df = df_scatter.query("Cell_type_img==@spe")
        fig.add_trace(go.Scatter(x=df.H_UMIs_cage, y=df.M_UMIs_cage, name=f"{spe} cell", mode="markers", marker_opacity=0.55, marker_size=4, marker_color=cell_type_colors[spe], text=df.index, hovertemplate=hovertemplate))
    fig.update_xaxes(range=[-30, maxi + 10], title="Ramos (Human) UMI counts")
    fig.update_yaxes(range=[-30, maxi + 10], title="B3Z (murine) UMI counts")
    fig.update_layout(height=400, width=500, margin=dict(t=10),
                      legend=dict(
                        title="Cell type from imaging",
                        font_size=10,
                        yanchor="top",
                        y=0.99,
                        xanchor="right",
                        x=0.99,
                        bgcolor='rgba(0,0,0,0)',
                     ), xaxis=dict(
                            tick0=0,
                            dtick=200,
                        ),
                        yaxis=dict(
                            tick0=0,
                            dtick=200,
                        )
                     )
    return fig

# Cell crop display on hover scatterplot UMI callback
@app.callback(
    Output("cell_display", "children", allow_duplicate=True),
    [Input('dropdown', 'value'),
     Input("scatterplot_UMI_sc", "hoverData")],
    prevent_initial_call=True
)
def cell_crop_on_hover_scatterUMI(selected_sample, hoverData):
    """
    Display cell crop of sc cages.
    """
    if hoverData is None or selected_sample is None:
        # If no dropdown value or hover data is selected, clear the image.
        return None

    chip_ID, seq_ID = selected_sample.replace("Chip", "").split("_")
    mag="20X"

    channel_available = np.unique([c.split("-T")[0] for c in sorted(os.listdir(f"{bigbi_path}/Chips/Chip{chip_ID}/{chip_ID}_images/{mag.lower()}/Cells/Crops/"))])
    # Fix RGB image in first
    if "RGB" in channel_available:
        channel_available = np.delete(channel_available, np.where(channel_available == "RGB")[0])
        channel_available = np.insert(channel_available, 0, "RGB")

    cage_ID = hoverData["points"][0]["text"]
    df_cells_sc = get_single_cells(chip_ID)
    df_cells_sc.set_index("Cage_ID", inplace=True)

    cell_items = []
    for channel in channel_available:
        cell_path = get_cell_crop_path(chip_ID, cage_ID, mag, df_cells_sc, channel)
        if os.path.exists(cell_path.replace("Minos_space/", bigbi_path)):
            cell_items.append(html.Div([
                html.Div(channel, style={"text-align": "center", "font-size": "12px", "color": "#555"}),
                html.Img(
                    src=dash.get_asset_url(cell_path),
                    alt=f"Image not available (chip {chip_ID}, {seq_ID})",
                    style={"width": "80px", "margin-bottom": "5px"}
                )
            ], style={"margin-right": "10px", "margin-top": "0px", "display": "inline-block", "text-align": "center"}))

    try:
        if not cell_items:
            raise Exception
        return html.Div(cell_items, style={"display": "flex", "justify-content": "center", "align-items": "flex-start"})
    except Exception as error:
        print(error)
        raise PreventUpdate


@app.callback(
    Output('top-n-input', 'style'),
    Output('threshold-input-container', 'style'),
    [Input('filter-option-radio', 'value')]
)
# Allows to display int input or threshold according to the selected filter option
def update_filter_input(filter_option):
    if filter_option == 'top_n':
        return {'margin-right': '10px'}, {'display': 'none'}
    elif filter_option == 'threshold':
        return {'display': 'none'}, {'margin-right': '10px'}
    return {}, {}

@app.callback(
    Output('correlation-pairs-dropdown', 'options'),
    [Input('dropdown', 'value'),
    Input('features-option-radio', 'value'),
    Input('filter-option-radio', 'value'),
    Input('top-n-input', 'value'),
    Input('threshold-input', 'value')]
)
# Get the pair of most correlated variables
def update_correlation_dropdown(selected_sample, feat_type, filter_option, top_n, threshold):
    # If any seq is selected, or if cell_type is selected for data type selection but this cell type isn't yet choose, return empty
    if selected_sample is None:
        return []

    chip_ID, seq_ID = selected_sample.replace("Chip", "").split("_")

    corr_matrix = pd.read_table(f"{bigbi_path}Chips/Chip{chip_ID}/{seq_ID}/{seq_ID}_data/{chip_ID}_{seq_ID}.correlation_greater_0.3.tsv", index_col=0)

    top_corr_pairs = get_top_correlated_pairs(corr_matrix, filter_option, top_n, threshold, feat_type)
    top_corr_pairs["Correlation"].astype(float)

    # Add the top correlated pairs into a dropdown list
    dropdown_options = [
        {
            'label': f"{i + 1}. {row['Variable1']} ({t1}) vs {row['Variable2']} ({t2})"
                     + f" (R: {row['Correlation']:.2f}, "
                     + f"p-value: {np.format_float_scientific(row['P-value'], precision=4) if row['P-value'] >= np.nextafter(0, 1) else '<1e-324'})",
            'value': f"{row['Variable1']},{row['Variable2']}"
        }
        for (i, row), t1, t2 in
        zip(top_corr_pairs.iterrows(), top_corr_pairs.Variable1_type, top_corr_pairs.Variable2_type)
    ]
    return dropdown_options

@app.callback(
    Output('scatter-plot', 'figure'),
    [Input('correlation-pairs-dropdown', 'value'),
     Input('dropdown', 'value')]
)
# Create the scatter plot of the 2 variables selected from the top correlated pairs
def update_scatter_plot(correlation_value, selected_sample):
    # If any seq is selected, or if the top pair isn't yet chosen, return empty
    if correlation_value is None or selected_sample is None:
        return {}

    chip_ID, seq_ID = selected_sample.replace("Chip", "").split("_")

    df_sc = get_single_cells(chip_ID)
    df_sc.loc[df_sc.Cell_Type == "H_Ramos", "Cell_Type"] = "Ramos_Human"  # Temp modif
    df_sc.loc[df_sc.Cell_Type == "M_B3Z", "Cell_Type"] = "B3Z_murine"  # Temp modif

    # Get df_count
    df_count = pd.read_table(f"{bigbi_path}Chips/Chip{chip_ID}/{seq_ID}/{seq_ID}_data/{chip_ID}_{seq_ID}.counts.tsv", index_col=0, dtype=str)
    df_count = df_count.loc[df_count.index != 'Type']
    df_count = df_count.astype(float)

    df_count["Cell_type_img"] = df_sc.set_index("Cage_ID").loc[df_count.index, "Cell_Type"]

    # Get the 2 variables of the selected top correlated pair from the dropdown list
    var1, var2 = correlation_value.split(',')
    # Add the coloration of the cages if needed
    fig = px.scatter(df_count, y=var1, x=var2, opacity=0.5, color="Cell_type_img", color_discrete_map=cell_type_colors, title=f"Correlation scatter plot of {var1} vs {var2}", custom_data=[df_count.index])
    #fig.update_yaxes(range=[min(df_count[var1].min(), df_count[var2].min()) - 1, max(df_count[var1].max(), df_count[var2].max()) + 1])
    #fig.update_xaxes(range=[min(df_count[var2].min(), df_count[var1].min()) - 1, max(df_count[var2].max(), df_count[var1].max()) + 1])
    fig.update_traces(
        hovertemplate="<b>Cage_ID: %{customdata[0]}</b><br>" +
                      "%{yaxis.title.text}: %{y}<br>" +
                      "%{xaxis.title.text}: %{x}<extra></extra>")
    fig.update_layout(legend=dict(
        title="Cell type from imaging",
        font_size=8,
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=0.01,
        bgcolor='rgba(0,0,0,0)',
    ))

    return fig


@app.callback(
    Output("cell_display", "children", allow_duplicate=True),
    [Input('dropdown', 'value'),
     Input('correlation-pairs-dropdown', 'value'),
     Input("scatter-plot", "hoverData")],
    prevent_initial_call=True
)
def cell_crop_on_hover_scatter_plot(selected_sample, correlation_dropdown, hoverData):
    """
    Display cell crop of sc cages.
    """
    if correlation_dropdown is None or hoverData is None or selected_sample is None:
        # If no dropdown value or hover data is selected, clear the image.
        return None

    chip_ID, seq_ID = selected_sample.replace("Chip", "").split("_")
    mag="20X"

    channel_available = np.unique([c.split("-T")[0] for c in sorted(os.listdir(f"{bigbi_path}/Chips/Chip{chip_ID}/{chip_ID}_images/{mag.lower()}/Cells/Crops/"))])
    # Fix RGB image in first
    if "RGB" in channel_available:
        channel_available = np.delete(channel_available, np.where(channel_available == "RGB")[0])
        channel_available = np.insert(channel_available, 0, "RGB")

    cage_ID = hoverData["points"][0]["customdata"][0]
    df_cells_sc = get_single_cells(chip_ID)
    df_cells_sc.set_index("Cage_ID", inplace=True)

    cell_crops_section = [html.Div("Cell crops:", style={"text-align": "center", "font-size": "12px", "color": "#555"})]
    for channel in channel_available:
        cell_path = get_cell_crop_path(chip_ID, cage_ID, mag, df_cells_sc, channel)
        if os.path.exists(cell_path.replace("Minos_space/", bigbi_path)):
            cell_crops_section.append(html.Div([
                html.Div(channel, style={"text-align": "center", "font-size": "12px", "color": "#555"}),
                html.Img(
                    src=dash.get_asset_url(cell_path),
                    alt=f"Image not available (chip {chip_ID}, {seq_ID})",
                    style={"width": "80px", "margin-bottom": "5px"}
                )
            ], style={"margin-right": "10px", "margin-top": "0px", "display": "inline-block", "text-align": "center"}))

    try:
        if not cell_crops_section:
            raise Exception
        return html.Div(cell_crops_section, style={"display": "flex", "flex-wrap": "wrap", "justify-content": "center"})
    except Exception as error:
        print(error)
        raise PreventUpdate


if __name__ == '__main__':
    import argparse  # Fixed a bug with gunicorn

    # Managing command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", dest="port", type=int, default=8050, help="Port number to display on.")
    args = parser.parse_args()

    app.run(debug=True, port=args.port)
