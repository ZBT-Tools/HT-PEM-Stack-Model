import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from pemfc.src import interpolation as ip
from pemfc.src import stack_simulation
from pemfc.gui import data_transfer
from flask_caching import Cache
import pickle

#external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]
external_stylesheets = [dbc.themes.BOOTSTRAP]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

CACHE_CONFIG = {
    "DEBUG": True,          # some Flask specific configs
    "CACHE_TYPE": "simple",  # Flask-Caching related configs
    "CACHE_DEFAULT_TIMEOUT": 300
}
cache = Cache()
cache.init_app(app.server, config=CACHE_CONFIG)

cell_input = dbc.FormGroup(
    [
        dbc.Label("Cell Number", html_for="cell_number", width=2),
        dbc.Col(
            dbc.Input(
                id="cell_number", type="number", min=0, max=100, step=1,
                value=1,
                placeholder="Enter number of cells in stack"
            ),
            width=2,
        ),
    ],
    row=True,
)


# output = dbc.Row()
app.layout = html.Div([
    cell_input, html.Button('Run Simulation', id='run_button'),
    html.P(id="output_stack_power"),
    dcc.Dropdown(id='results-dropdown'),
    html.Div(dcc.Graph(id='graph')),
    # hidden signal value
    html.Div(id='signal', style={'display': 'none'})])
# app.layout = html.Div(
#     [
#         # html.I("Inputs"),
#         # html.Br(),
#         # dcc.Input(id="input1", type="text", placeholder=""),
#         # dcc.Input(id="input2", type="text", placeholder="", debounce=True),
#         # html.Div(id="output"),
#         cell_input
#     ]
# )


@cache.memoize()
def simulation_store(cell_number):
    values = {'cell number': {'sim_name': ['stack', 'cell_number'],
                              'gui_name': 'Cell Number:', 'value': cell_number}}
    data_transfer.gui_to_sim_transfer(values, data_transfer.sim_dict)
    # print(data_transfer.sim_dict)
    global_data, local_data, sim = stack_simulation.main()
    return [global_data, local_data]


@app.callback(
    Output('signal', 'children'),
    [Input("cell_number", "value"), Input("run_button", "n_clicks")]
)
def compute_simulation(cell_number, n_clicks):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if 'run_button' in changed_id and n_clicks is not None:
        simulation_store(cell_number)
    else:
        raise PreventUpdate
    return cell_number


@app.callback(
    Output('results-dropdown', 'options'), [Input('signal', 'children')]
)
def get_dropdown_options(value):
    results = simulation_store(value)
    local_data = results[1]
    return [{'label': key, 'value': key} for key in local_data]


@app.callback(
    Output("graph", "figure"),
    [Input('signal', 'children'), Input('results-dropdown', 'value')])
def update_graph(value, dropdown_key):
    results = simulation_store(value)
    global_data = results[0]
    local_data = results[1]
    result = 'Stack power [W/m²]: {:05.2f}'.format(
        global_data['Stack Power']['value'])

    #z_key = 'Current Density'
    z_key = dropdown_key
    x_key = 'Channel Location'
    y_key = 'Cells'
    xvalues = ip.interpolate_1d(local_data[x_key]['value'])
    yvalues = local_data[y_key]['value']
    # fig = px.imshow(local_data[z_key]['value'],
    #                 labels=dict(x=x_key, y=y_key, color=z_key),
    #                 x=xvalues, y=yvalues, width=1000, height=600,
    #                 aspect='auto')
    fig = go.Figure(go.Heatmap(z=local_data[z_key]['value'],
                    x=xvalues, y=yvalues, xgap=1, ygap=1))
    fig.update_xaxes(showgrid=True, tickmode='array',
                     tickvals=local_data[x_key]['value'])
    fig.update_yaxes(showgrid=True, tickmode='array', tickvals=yvalues)
    # fig.update_layout(ygap=2)
    return fig


if __name__ == "__main__":
    app.run_server(debug=True)
