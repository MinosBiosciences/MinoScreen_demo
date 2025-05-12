import argparse
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import pandas as pd
import plotly.express as px

# Sample data
df = pd.DataFrame({
    'Sample': ['Sample1', 'Sample1', 'Sample2', 'Sample2', 'Sample3', 'Sample3'],
    'X': [1, 2, 3, 4, 5, 6],
    'Y': [10, 11, 12, 13, 14, 15]
})

# Initialize the Dash app
app = dash.Dash(__name__)

# Define app layout
app.layout = html.Div([
    html.Label("Select Samples:"),
    dcc.Dropdown(
        id='sample-dropdown',
        options=[{'label': sample, 'value': sample} for sample in df['Sample'].unique()],
        multi=True,
        value=['Sample1', 'Sample2', 'Sample3']
    ),
    html.Div([
        dcc.Graph(id='scatter-plot-sample1'),
        dcc.Graph(id='scatter-plot-sample2')
    ], style={'display': 'flex'})
])

# Define callback to update scatter plots based on selected samples
@app.callback(
    [Output('scatter-plot-sample1', 'figure'),
     Output('scatter-plot-sample2', 'figure')],
    [Input('sample-dropdown', 'value')]
)
def update_scatter_plots(selected_samples):
    figures = []

    for sample in selected_samples:
        filtered_df = df[df['Sample'] == sample]
        scatter_plot = px.scatter(filtered_df, x='X', y='Y', title=f'Scatter Plot - {sample}')
        figures.append(scatter_plot)

    return figures

# Run the app
# if __name__ == '__main__':
#     app.run_server(debug=True)

# Managing command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--port", dest="port", type=int, default=8050, help="Port number to display on.")
args = parser.parse_args()


if __name__ == '__main__':
    app.run_server(debug=True, port=args.port)