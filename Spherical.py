### Import all extra functions ###

import plotly.graph_objects as go
import numpy as numpy
import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import time
import Utils ### Custom function scripts, Thanks to Tom Irons of the Teale group (UoN School of Chemistry) for molecule generation code. 19/7/22 ###


### Variables ###
work_dir = "" ### Folder with all relevant files. Needs to end in / ###
server = "" ### Set to login001.augusta.nottingham.ac.uk for HPC access, if "LOCAL", no server will be used. ###
qm_package = "QChem" ### Set to either QChem or GAUSSIAN ###
host_work_dir = "" ### Folder on server for calculations. Needs to end in / ###
start_coordinate = "" ### Currently only supports .xyz files, although new file types can be added easily ###
calculate_data = True ### Whether to run the data calculations (if you're just doing analysis, set to False) ###
qm_functional = ""
qm_basis_set = ""
qm_dispersion = ""
system_spin = "" ### Spin/multiplicity in the format the QM package requires, no calculation is performed on this and is used "as is" ###
system_charge = ""
spherical_datapoints = 200 ### Number of single point calculations required ###
charge_value = 3 ### Value of point charges around the sphere ###
radius = 30 ### Radius of sphere ###
threads = 4 ### Number of threads to give each single point calculation ###
mem = 2 ### Ammount of memory to give each single point calculation ###
random = True ### Decides whether spherical points are random or uniform ###
gen_newcoord = True ### Decides whether to generate new coordinate file in format for later visualisation ###
recenter = True ### Decides whether to recenter the coordinates if not already in the center. ###
center_type = "" ### either simple, mass weighted, atom, or translate ###
center_info = [] ### if atom, put atom number in array, if translate, put coordinates in array ([x,y,z]) ###
coord_filename = "coords_formatted.txt" ### Coord file, with format X, Y, Z, Atomic number (hint: H = 1, C = 6, N = 7, O = 8), (only change if adding own file) ###
energies_filename = "Energy.xyzc" ### Spherical energy file, with format X, Y, Z, C (only change if you're adding your own file) ###
resolution= ### Resolution of sphere, increase to improve contours but decrease for speed, Suggest 100 ###


### Calculation setup ###
if calculate_data == True:
    (type_atom, x_atom, y_atom, z_atom) = Utils.CoordGet(work_dir,start_coordinate) ### loads in initial coordinates ###

    (x_atom, y_atom, z_atom) = Utils.CoordCheck(type_atom, x_atom, y_atom, z_atom, recenter, center_type, center_info) ### Checks the coordinates to make sure they exist and locates the center of the system ###

    if gen_newcoord == True:    Utils.FormatCoord(type_atom, x_atom, y_atom, z_atom,work_dir) ### Generates formatted coordinate file if requested ###

    if random == True: Utils.RandomCoords(work_dir,spherical_datapoints,radius) ### Generates random spherical points ###
    else: ### Generates non-random points about the sphere ###
        print("Warning, this method is currently unsupported and probably doesnt work. DONT USE UNLESS YOU HAVE FIXED IT")
        (px,py,pz) = Utils.SphereGen(0,0,0,radius,spherical_datapoints)
        with open(work_dir + "ChargeP.xyz",'w') as f:
            for i in range(len(px)):
                print(str(px[i])+"\t"+str(py[i])+"\t"+str(pz[i]),file=f)
        with open(work_dir + "ChargeN.xyz",'w') as f:
            for i in range(len(px)):
                print(str(-px[i])+"\t"+str(-py[i])+"\t"+str(-pz[i]),file=f)

    print("Generating input files")

    Utils.GenSP(type_atom, x_atom, y_atom, z_atom, work_dir, charge_value, mem, threads, qm_functional, qm_basis_set, system_charge, system_spin, qm_dispersion, qm_package)

    Utils.SPRun(work_dir,mem,threads,spherical_datapoints,host_work_dir,server,qm_package)
else:
    print("Data should already be calculated.")
    
failed = Utils.OutputExtract(work_dir,spherical_datapoints,threads,qm_package)


start_time = time.perf_counter()
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

### Extract energy data ###
(x_data, y_data, z_data, en_data) =  Utils.EnergyExtract(work_dir,energies_filename)
MolData = numpy.genfromtxt(work_dir + coord_filename)

data_loaded = time.perf_counter()

print(f"Time taken to import data is {data_loaded - start_time} seconds")

### Generate spheres ###
(x_sphere, y_sphere, z_sphere) = Utils.SphereGen(0,0,0,1,resolution)

colour_start = time.perf_counter()

print(f"Time taken to generate sphere is {colour_start - data_loaded} seconds")

### Generate colors and contours ###
(col_sphere, x_contour, y_contour, z_contour, val_contour) = Utils.GenColors(x_sphere, y_sphere, z_sphere, x_data, y_data, z_data, en_data, radius)

color_end = time.perf_counter()

print(f"Time taken to interpolate and generate color data is {color_end - colour_start} seconds")

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.layout = html.Div([
#     html.Div(children=[html.H1(children="Spherical visualisation tool",style={'text-align' : 'center'})]),
#     html.Div(children=[
#         html.Div(children=[
#             html.H3(children="Sphere opacity",style={'display': 'inline-block'}),
#             dcc.Slider(id='alphval', min=0, max=1, step=0.05, value=0.8,style={'display': 'inline-block'})
#         ],),
#     ],),
#     html.Div(children=[
#         html.Div(children=[
#             html.H3(children="Interpolation data",style={'display': 'inline-block'}),
#             dcc.RadioItems(id='interp',options=["None","Sphere","Contours"],value="Sphere",inline=True,style={'display': 'inline-block'}),
#         ]),
#     ]),
    html.H1(children="Spherical visualisation tool",style={'text-align' : 'center'}),
    html.H3("Sphere opacity"),
    dcc.Slider(id='alphval', min=0, max=1, step=0.05, value=0.8),
    html.H3("Interpolation data"),
    dcc.RadioItems(id='interp',options=["None","Sphere","Contours"],value="Sphere",inline=True),
    html.H3("Molecule texture style"),
    dcc.RadioItems(id='texture',options=["matte","shiny","orbs"],value="orbs",inline=True),
    html.H3("Molecule draw style"),
    dcc.RadioItems(id='draw_type',options=["ball_and_stick","tubes","wireframe","spacefilling"],value="tubes",inline=True),
    dcc.Graph(id="graph"),
])


@app.callback(Output('graph','figure'),Input('alphval','value'),Input('interp','value'),Input('texture','value'),Input('draw_type','value'))
def update_graph(alphval,interp,texture,draw_type):
    Fig = go.Figure(layout = go.Layout(title="Visualisation",uirevision='camera'))
    Fig.update_layout(autosize=True, width=1000, height=1000)
    Molecule = Utils.DrawMolecule(MolData,texture,draw_type)
    for Bond in Molecule['bond_list']: Fig.add_trace(Bond)
    for Atom in Molecule['atom_list']: Fig.add_trace(Atom)
    Fig.update_layout(Utils.GetLayout(None))
    # MinRange = numpy.min(Molecule['geometry'])
    # MaxRange = numpy.max(Molecule['geometry'])
    # Fig.update_layout(Utils.GetRange(MinRange,MaxRange))
    Sphere = go.Surface(x=x_sphere*radius, y=y_sphere*radius, z=z_sphere*radius, customdata=col_sphere, opacity=alphval,surfacecolor=col_sphere,contours={"x":{"show":True},"y":{"show":True},"z":{"show":True}},colorscale='Turbo')
    Energies = go.Scatter3d(x=x_data,y=y_data,z=z_data,mode='markers',marker=dict(size=5,color=en_data,colorscale='Turbo',opacity=1))
    Contours = go.Scatter3d(x=x_contour,y=y_contour,z=z_contour,mode='markers',marker=dict(size=3,color='black'),hovertext=val_contour)
    Fig.add_trace(Energies)
    if interp == "Sphere" or "Contours": Fig.add_trace(Sphere)
    if interp == "Contours":  Fig.add_trace(Contours)
    return Fig

app.run_server(debug=True, use_reloader=False, host='0.0.0.0')
