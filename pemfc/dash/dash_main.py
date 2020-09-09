import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
from pemfc.src import stack_simulation
from pemfc.gui import data_transfer

#external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]
external_stylesheets = [dbc.themes.BOOTSTRAP]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

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

output = dbc.Row()
app.layout = dbc.Container(
    [cell_input, html.Button('Run Simulation', id='run_simulation'),
     html.P(id="output_stack_power")])
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


@app.callback(
    Output("output_stack_power", "children"),
    [Input("cell_number", "value"), Input("run_simulation", "n_clicks")],
)
def update_output(cell_number, n_clicks):
    if not isinstance(cell_number, int):
        result = 'Integer number must be provided'
    else:
        if n_clicks is not None:
            values = {'cell number': {'sim_name': ['stack', 'cell_number'],
                      'gui_name': 'Cell Number:', 'value': cell_number}}
            data_transfer.transfer(values, data_transfer.sim_dict)
            print(data_transfer.sim_dict)
            avg_icd = stack_simulation.main()
            result = 'Stack power [W/m²]: {}'.format(avg_icd)
        else:
            result = None
    return result


if __name__ == "__main__":
    app.run_server(debug=True)
