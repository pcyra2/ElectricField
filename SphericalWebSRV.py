### Import all extra functions ###

from  plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as numpy
import dash
from dash import dcc
from dash import html, ctx
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import time
import Utils ### Custom function scripts, Thanks to Tom Irons of the Teale group (UoN School of Chemistry) for molecule generation code. 19/7/22 ###
from collections import OrderedDict
import multiprocessing
import pandas as pd
from IPython.display import display, HTML
import json
import webbrowser
from threading import Timer
display(HTML(""))

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.LUX])
app.layout = html.Div([
    html.H1(children="Spherical visualisation tool",
            style={'text-align': 'center'},
            className="bg-primary text-white m-1 p-3 text-decoration-underline"),
    dbc.Row([
        html.H2("Global Variables", style={'text-align': 'center'}, className="bg-secondary text-black m-1 p-3 text-decoration-underline"),
        dbc.Row([
            html.H3("Working Directory"),
            dcc.Input(id='work_dir', type='text',value="/home/pcyra2/TEST/1/", required=True, placeholder='~/Data/')
        ]),
        dbc.Row([
            html.H3("Coordinate File"),
            dcc.Input(id='file', type='text',value="H2O.xyz", required=True, placeholder='coords.xyz')
        ]),
    ], className="border-dark border-bottom border-3 m-1 p-3 "),
    dbc.Row([
        dbc.Col([
            html.Div([
                html.H2("Computational Variables", style={'text-align': 'center'}, className="bg-secondary text-black m-1 p-3 text-decoration-underline"),
                dbc.Row([
                    dbc.Col([
                        html.H3("QM Variables:"),
                        html.H5("QM Package"),
                        dcc.RadioItems(id='package', options=['QChem', 'Gaussian'], value='QChem', inline=True),
                        html.H5("DFT Functional"),
                        dcc.Input(id='functional', type='text',value="B3LYP", required=True, placeholder='B3LYP')
                        ]),
                    dbc.Col([
                        html.H5("DFT Basis Set"),
                        dcc.Input(id='basis', type='text',value="6-31+G*", required=True, placeholder='6-31+G*')
                        ]),
                    dbc.Col([
                        html.H5("DFT Dispersion "),
                        dcc.Input(id='dispersion', type='text',value="EMPIRICAL_GRIMME", required=True, placeholder='6-31+G*')
                        ]),
                    ], className="border-dark border border-1 p-1 m-2"),
                dbc.Row([
                    html.H3("System:"),
                    html.H5("recenter"),
                    dcc.RadioItems(id='recenter', options=['None', 'simple', 'mass weighted', 'atom', 'translate'], value='simple', inline=True),
                    dbc.Col([
                        html.H5("Charge"),
                        dcc.Input(id='charge', type='number', required=True, value=0)
                    ]),
                    dbc.Col([
                        html.H5("Spin"),
                        dcc.Input(id='spin', type='number', required=True, value=1)
                    ]),
                ], className="border-dark border border-1 p-1 m-2"),
                dbc.Row([
                    html.H3("Sphere:"),
                    # html.H6("leave blank if using previous data positions"),
                    dbc.Row([
                            html.H5("Use existing charge positions"),
                            dcc.RadioItems(id='new_chg', options=['True', 'False'],
                                           value='False', inline=True)
                        ]),
                    dbc.Col([
                        html.H5("Charge"),
                        dcc.Input(id='sphere-charge', type='number', required=True, value=3)
                    ], className="w-30"),
                    dbc.Col([
                        html.H5("Radius"),
                        dcc.Input(id='sphere-radius', type='number', required=True, value=30)
                    ], className="w-30"),
                    dbc.Col([
                        html.H5("Points"),
                        dcc.Input(id='num-points', type='number', required=True, value=300)
                    ], className="w-30"),
                ], className="border-dark border border-1 p-1 m-2"),
                dbc.Row([
                    html.H3("Calculation:"),
                    dbc.Col([
                        html.H5("Threads"),
                        dcc.Input(id='threads', type='number', required=True, value=4)
                    ]),
                    dbc.Col([
                        html.H5("Memory"),
                        dcc.Input(id='mem', type='number', required=True, value=1)
                    ]),
                ], className="border-dark border border-1 p-1 m-2"),
                dbc.Row([
                    html.H3("Server"),
                    dcc.RadioItems(id='SRV',options=['True', 'False'], value='True', inline=True),
                    dbc.Col([
                        dcc.Input(id="server", type="text", required=True, value="login001.augusta.nottingham.ac.uk", className="w-50"),
                        dcc.Input(id="srvDir", type='text', required=True, value='~/tmp/', className="w-50")
                        ])
                    ], className="border-dark border border-1 p-1 m-2"),
                dbc.Row([
                    html.Button("Run Calculation", id="calc_run", n_clicks=0)
                ]),
                dcc.Loading(
                    id="loading-1",
                    type="default",
                    children=html.H1(id="loading-output-1"),
                    fullscreen=True,
                ),
            ], className=""),
        ], className="m-1 p-3 align-middle "),
        dbc.Col([
            html.Div([
                html.H2("Interpolation Variables", style={'text-align': 'center'},
                        className="bg-secondary text-black m-1 p-3 text-decoration-underline"),
                html.Div([
                    dbc.Row([
                        html.H3("Calculate data"),
                        dcc.RadioItems(id='calc-data', options=['True', 'False'],
                                           value='False', inline=True)
                        ]),
                        dbc.Row([
                            html.H3("Save Interpolation"),
                            dcc.RadioItems(id='interp save', options=['True', 'False'],
                                           value='True', inline=True)
                        ]),
                        dbc.Col([
                            html.H5("Visualisation Resolution"),
                            dcc.Input(id='res', type='number', required=True, value=400)
                        ]),
                ],className="border-dark border border-1 p-1 m-2"),
                dbc.Row([
                    html.Button("Run Interpolation", id="interp_run", n_clicks=0)
                ]),
                dcc.Loading(
                    id="loading-2",
                    type="default",
                    children=html.H1(id="loading-output-2"),
                    fullscreen=True,
                ),
            ], className="m-1 p-3 align-middle "),
            html.Div([
                html.H2("Visualisation Variables", style={'text-align': 'center'}, className="bg-secondary text-black m-1 p-3 text-decoration-underline"),
                html.Div([
                    dbc.Row([
                        html.H3("Sphere opacity"),
                        dcc.Slider(id='alphval', min=0, max=1, step=0.05, value=0.8),
                    ]),
                    dbc.Row([
                        html.H3("Interpolation data"),
                        dcc.RadioItems(id='interp',options=["None","Sphere","Contours"],value="Sphere",inline=True),
                    ]),
                    dbc.Row([
                        html.H3("Molecule texture style"),
                        dcc.RadioItems(id='texture',options=["matte","shiny","orbs"],value="orbs",inline=True),
                    ]),
                    dbc.Row([
                        html.H3("Molecule draw style"),
                        dcc.RadioItems(id='draw_type',options=["ball_and_stick","tubes","wireframe","spacefilling"],value="tubes",inline=True),
                    ]),
                ],className="border-dark border border-1 p-1 m-2"),
                dbc.Row([
                   html.Button("Run Visualisation", id="vis_run", n_clicks=0)
                ]),
                dcc.Loading(
                    id="loading-3",
                    type="default",
                    children=html.Div(id="loading-output-3"),
                    fullscreen=True,
                ),
            ],className="m-1 p-3 align-middle"),
        ], className=""),
    ], className="border-dark border-bottom border-3 m-1 p-3 "),
    html.Div(dcc.Graph(id="graph"), className="align-middle border-dark border-bottom border-3 m-1 p-3"),
    dbc.Row([
        html.H2("Data Comparison", style={'text-align': 'center'},
                className="bg-secondary text-black m-1 p-3 text-decoration-underline"),
        html.Div([
            dbc.Row([
                html.H3("Dataset 1"),
                dcc.Input(id='work_dir1', type='text', value="/home/pcyra2/TEST/1/",
                          required=True, placeholder='~/Data/')
            ]),
            dbc.Row([
                html.H3("Dataset 2"),
                dcc.Input(id='work_dir2', type='text', value="/home/pcyra2/TEST/2/", required=True,)
            ]),
            dbc.Row([
                html.H3("Save Location"),
                dcc.Input(id='save_location', type='text', value="/home/pcyra2/TEST/", required=True,)
            ]),
        ],className="border-dark border border-1 p-1 m-2"),
        dbc.Row([
            html.Button("Run Comparison", id="Comp_run", n_clicks=0)]),
            dcc.Loading(
                id="loading-4",
                type="default",
                children=html.Div(id="loading-output-4"),
                fullscreen=True,),
    ], className="border-dark border-bottom border-3 m-1 p-3 "),
    html.Button("Run Visualisation", id="Comp_vis", n_clicks=0),
    html.Div(dcc.Graph(id="graph2"),className="align-middle border-dark border-bottom border-3 m-1 p-3"),
])


@app.callback(Output('work_dir','value'),
              Output('file','value'),
              Output('package','value'),
              Output('functional','value'),
              Output('basis','value'),
              Output('dispersion','value'),
              Output('recenter','value'),
              Output('charge','value'),
              Output('spin','value'),
              Output('new_chg','value'),
              Output('sphere-charge','value'),
              Output('sphere-radius','value'),
              Output('num-points','value'),
              Output('threads','value'),
              Output('mem','value'),
              Output('SRV','value'),
              Output('server','value'),
              Output('srvDir','value'),
              Output('calc-data','value'),
              Output('interp save','value'),
              Output('res','value'),
              Output('alphval','value'),
              Output('interp','value'),
              Output('texture','value'),
              Output('draw_type','value'),
              Output('work_dir1','value'),
              Output('work_dir2','value'),
              Output('save_location','value'),
              Input('calc_run','n_clicks'),
              Input('interp_run','n_clicks'),
              Input('vis_run','n_clicks'),
              Input('Comp_run','n_clicks'),
              State('work_dir','value'),
              State('file','value'),
              State('package','value'),
              State('functional','value'),
              State('basis','value'),
              State('dispersion','value'),
              State('recenter','value'),
              State('charge','value'),
              State('spin','value'),
              State('new_chg','value'),
              State('sphere-charge','value'),
              State('sphere-radius','value'),
              State('num-points','value'),
              State('threads','value'),
              State('mem','value'),
              State('SRV','value'),
              State('server','value'),
              State('srvDir','value'),
              State('calc-data','value'),
              State('interp save','value'),
              State('res','value'),
              State('alphval','value'),
              State('interp','value'),
              State('texture','value'),
              State('draw_type','value'),
              State('work_dir1', 'value'),
              State('work_dir2', 'value'),
              State('save_location', 'value'),
              )
def retain_state(nclicks1,nclicks2,nclicks3,nclicks4,
                work_dir,file,package,functional,basis,disp,recenter,chg,spin,new_chg,sphchg,sphr,num_points,threads,mem,SRV,server,SRVDIR,calc_data,interpsv,res,alphval,interp,text,draw, wd1, wd2, save_loc):
    if str(work_dir).endswith("/") == False:
        work_dir =  str(work_dir)+"/"
    if str(wd1).endswith("/") == False:
        wd1 =  str(wd1)+"/"
    if str(wd2).endswith("/") == False:
        wd2 =  str(wd2)+"/"
    if str(save_loc).endswith("/") == False:
        save_loc =  str(save_loc)+"/"
    if str(SRVDIR).endswith("/") == False:
            SRVDIR = str(SRVDIR) + "/"
    if nclicks1 != 0 or nclicks2 != 0 or nclicks3 != 0 or nclicks4 != 0:
        vars = {'work_dir' : work_dir,
            'file' : file,
            'package' : package,
            'functional' : functional,
            'basis' : basis,
            'dispersion' : disp,
            'recenter' : recenter,
            'charge' : chg,
            'spin' : spin,
            'new_chg' : new_chg,
            'sphere-charge' : sphchg,
            'sphere-radius' : sphr,
            'num_points' : num_points,
            'threads' : threads,
            'mem' : mem,
            'SRV' : SRV,
            'server' : server,
            'srvDir' : SRVDIR,
            'calc-data' : calc_data,
            'interp save' : interpsv,
            'res' : res,
            'alphval' : alphval,
            'interp' : interp,
            'texture' : text,
            'draw_type' : draw,
            'work_dir1' : wd1,
            'work_dir2' : wd2,
            'save_location' : save_loc,
        }
        json.dump(vars, open("cache.dict",'w'))
        return (work_dir,file,package,functional,basis,disp,recenter,chg,spin,new_chg,sphchg,sphr,num_points,threads,mem,SRV,server,SRVDIR,calc_data,interpsv,res,alphval,interp,text,draw, wd1, wd2, save_loc)
    else:
        try:
            cache = json.load(open("cache.dict", 'r'))
        except FileNotFoundError:
            print("No Cache file")
            return (
            work_dir, file, package, functional, basis, disp, recenter, chg,
            spin, new_chg, sphchg, sphr, num_points, threads, mem, SRV,
            server, SRVDIR, calc_data, interpsv, res, alphval, interp, text,
            draw)
        work_dir = cache['work_dir']
        file = cache['file']
        package = cache['package']
        functional = cache['functional']
        basis = cache['basis']
        disp = cache['dispersion']
        recenter = cache['recenter']
        chg = cache['charge']
        spin = cache['spin']
        new_chg = cache['new_chg']
        sphchg = cache['sphere-charge']
        sphr = cache['sphere-radius']
        num_points = cache['num_points']
        threads = cache['threads']
        mem = cache['mem']
        SRV = cache['SRV']
        server = cache['server']
        SRVDIR = cache['srvDir']
        calc_data = cache['calc-data']
        interpsv = cache['interp save']
        res = cache['res']
        alphval = cache['alphval']
        interp = cache['interp']
        text = cache['texture']
        draw = cache['draw_type']
        wd1 = cache['work_dir1']
        wd2 = cache['work_dir2']
        save_loc = cache['save_location']
        return (work_dir, file, package, functional, basis, disp, recenter, chg, spin, new_chg, sphchg, sphr, num_points, threads, mem, SRV, server, SRVDIR, calc_data, interpsv, res, alphval, interp, text, draw, wd1, wd2, save_loc)


@app.callback(Output('loading-1','children'),
                Input('calc_run','n_clicks'),
                State('work_dir', 'value'),
                State('file', 'value'),
                State('functional', 'value'),
                State('basis', 'value'),
                State('charge', 'value'),
                State('spin', 'value'),
                State('sphere-charge', 'value'),
                State('sphere-radius', 'value'),
                State('num-points', 'value'),
                State('threads', 'value'),
                State('mem', 'value'),
                State('server','value'),
                State('srvDir','value'),
                State('SRV','value'),
                State('package', 'value'),
                State('dispersion','value'),
                State('recenter', 'value'),
                )
def Calculate_Data(nclicks, work_dir, start_coordinate, qm_functional, qm_basis_set, system_charge, system_spin, charge_value, radius, spherical_datapoints, threads, mem,server, host_work_dir, locRun, qm_package, qm_dispersion, recenter):
    if nclicks > 0:
        if str(work_dir).endswith("/") == False:
            work_dir = str(work_dir) + "/"
        if str(host_work_dir).endswith("/") == False:
            host_work_dir = str(host_work_dir) + "/"
        if str(locRun) =="True":
            HPC="Y"
        else:
            HPC="N"
        gen_newcoord = True
        if str(recenter) == "None":
            recenter = False
            center_type = ""
        else:
            center_type = str(recenter)
            recenter = True
            print("Recentering data")
        calculate_data = True
        run_clean = True
        center_info=[]
        random=True ### For future implementation of systematic datapoints.
        if calculate_data == True:
            if run_clean == True:
                Utils.RunClean(work_dir, server, host_work_dir)
            (type_atom, x_atom, y_atom, z_atom) = Utils.CoordGet(str(work_dir), str(start_coordinate))  ### loads in initial coordinates ###
            (x_atom, y_atom, z_atom) = Utils.CoordCheck(type_atom, x_atom, y_atom, z_atom, recenter, center_type, center_info)  ### Checks the coordinates to make sure they exist and locates the center of the system ###
            if gen_newcoord == True:
                Utils.FormatCoord(type_atom, x_atom, y_atom, z_atom, work_dir)  ### Generates formatted coordinate file if requested ###
            if random == True:
                Utils.RandomCoords(str(work_dir), spherical_datapoints, radius)  ### Generates random spherical points ###
            else:  ### Generates non-random points about the sphere ###
                print(
                    "Warning, this method is currently unsupported and probably doesnt work. DONT USE UNLESS YOU HAVE FIXED IT")
                (px, py, pz) = Utils.SphereGen(0, 0, 0, radius, spherical_datapoints)
                with open(str(work_dir) + "ChargeP.xyz", 'w') as f:
                    for i in range(len(px)):
                        print(str(px[i]) + "\t" + str(py[i]) + "\t" + str(pz[i]), file=f)
                with open(str(work_dir) + "ChargeN.xyz", 'w') as f:
                    for i in range(len(px)):
                        print( str(-px[i]) + "\t" + str(-py[i]) + "\t" + str(-pz[i]), file=f)
            Utils.GenSP(type_atom, x_atom, y_atom, z_atom, str(work_dir), charge_value,
                        mem, threads, str(qm_functional), str(qm_basis_set), system_charge,
                        system_spin, str(qm_dispersion), str(qm_package))
            Utils.SPRun(str(work_dir), mem, threads, spherical_datapoints,
                        str(host_work_dir), str(server), str(qm_package))
        else:
            print("Data should already be calculated.")
        return "Calculation Complete"


@app.callback(Output('loading-2','children'),
              Input('interp_run','n_clicks'),
              State('work_dir','value'),
              State('file','value'),
              State('res','value'),
              State('sphere-radius','value'),
              State('interp save','value')
              )
def Interpolate_Data(nclicks, work_dir, file, resolution, radius, save):
    if nclicks > 0:
        if str(work_dir).endswith("/") == False:
            work_dir = str(work_dir) + "/"
        start_time = time.perf_counter()
        energies_filename = "Energy.xyzc"
        coord_filename = "coords_formatted.txt"
        (x_data, y_data, z_data, en_data) =  Utils.EnergyExtract(str(work_dir),energies_filename)
        MolData = numpy.genfromtxt(str(work_dir) + coord_filename)
        data_loaded = time.perf_counter()
        print(f"Time taken to import data is {data_loaded - start_time} seconds")
        ### Generate spheres ###
        (x_sphere, y_sphere, z_sphere) = Utils.SphereGen(0,0,0,1,resolution)
        colour_start = time.perf_counter()
        print(f"Time taken to generate sphere is {colour_start - data_loaded} seconds")
        ### Generate colors and contours ###
        (col_sphere, x_contour, y_contour, z_contour, val_contour) = Utils.GenColorsFaster(x_sphere, y_sphere, z_sphere, x_data, y_data, z_data, en_data, radius)
        color_end = time.perf_counter()
        print(f"Time taken to interpolate and generate color data is {color_end - colour_start} seconds")
        if save == "True":
            data = {"x_sphere" : x_sphere.tolist(),
                    "y_sphere" : y_sphere.tolist(),
                    "z_sphere" : z_sphere.tolist(),
                    "col_sphere" : col_sphere.tolist(),
                    "x_contour" : x_contour,
                    "y_contour" : y_contour,
                    "z_contour" : z_contour,
                    "val_contour" : val_contour,
                    }
            # with open(str(work_dir)+"interpolation.dict", 'w') as file:
            #     print(data, file=file)
            json.dump(data, open(str(work_dir)+"interpolation.dict", 'w'))
        return "Interpolation completed"


@app.callback(Output('loading-4','children'),
              Input('Comp_run','n_clicks'),
              State('work_dir1','value'),
              State('work_dir2','value'),
              State('res','value'),
              State('sphere-radius','value'),
              State('interp save','value'),
              State('save_location', 'value'),
    )
def Compare_Data(nclicks, work_dir1, work_dir2, resolution, radius, save, save_loc):
    if nclicks > 0:
        if str(work_dir1).endswith("/") == False:
            work_dir1 = str(work_dir1) + "/"
        if str(work_dir2).endswith("/") == False:
            work_dir2 = str(work_dir2) + "/"
        if str(save_loc).endswith("/") == False:
            save_loc = str(save_loc) + "/"
        start_time = time.perf_counter()
        energies_filename = "Energy.xyzc"
        coord_filename = "coords_formatted.txt"
        (x_data, y_data, z_data, en_data1) =  Utils.EnergyExtract(str(work_dir1),energies_filename)
        (x_data, y_data, z_data, en_data2) =  Utils.EnergyExtract(str(work_dir2),energies_filename)
        MolData = numpy.genfromtxt(str(work_dir1) + coord_filename)
        en_data_ar=numpy.subtract(numpy.array(en_data1),numpy.array(en_data2))
        en_data = list(en_data_ar)
        with open(str(save_loc)+"Energy.xyzc",'w') as f:
            for i in range(len(en_data)):
                print(str(x_data[i])+"\t"+str(y_data[i])+"\t"+str(z_data[i])+"\t"+str(en_data[i]), file=f)
        data_loaded = time.perf_counter()
        print(f"Time taken to import data is {data_loaded - start_time} seconds")
        ### Generate sphere ###
        (x_sphere, y_sphere, z_sphere) = Utils.SphereGen(0,0,0,1,resolution)
        colour_start = time.perf_counter()
        print(f"Time taken to generate sphere is {colour_start - data_loaded} seconds")
        ### Generate colors and contours ###
        (col_sphere, x_contour, y_contour, z_contour, val_contour) = Utils.GenColorsFaster(x_sphere, y_sphere, z_sphere, x_data, y_data, z_data, en_data, radius)
        color_end = time.perf_counter()
        print(f"Time taken to interpolate and generate color data is {color_end - colour_start} seconds")
        if save == "True":
            data = {"x_sphere" : x_sphere.tolist(),
                    "y_sphere" : y_sphere.tolist(),
                    "z_sphere" : z_sphere.tolist(),
                    "col_sphere" : col_sphere.tolist(),
                    "x_contour" : x_contour,
                    "y_contour" : y_contour,
                    "z_contour" : z_contour,
                    "val_contour" : val_contour,
                    }
            # with open(str(work_dir)+"interpolation.dict", 'w') as file:
            #     print(data, file=file)
            json.dump(data, open(str(save_loc)+"interpolation_comparison.dict", 'w'))
        return "Interpolation completed"

@app.callback(Output('graph', 'figure'),
              Input('vis_run', 'n_clicks'),
              State('alphval', 'value'),
              State('interp', 'value'),
              State('texture', 'value'),
              State('draw_type', 'value'),
              State('work_dir','value'),
              State('sphere-radius','value'),
)
def Visualise_Data(n_clicks,alphval,interp,texture,draw_type,work_dir,radius):
    if n_clicks > 0:
        if str(work_dir).endswith("/") == False:
            work_dir = str(work_dir) + "/"
        energies_filename = "Energy.xyzc"
        MolData = numpy.genfromtxt(str(work_dir) + "coords_formatted.txt")
        sphdat=json.load(open(str(work_dir)+"interpolation.dict"))
        x_sphere = numpy.array(sphdat["x_sphere"])
        y_sphere = numpy.array(sphdat["y_sphere"])
        z_sphere = numpy.array(sphdat["z_sphere"])
        col_sphere = numpy.array(sphdat["col_sphere"])
        x_contour = sphdat["x_contour"]
        y_contour = sphdat["y_contour"]
        z_contour = sphdat["z_contour"]
        val_contour = sphdat["val_contour"]
        (x_data, y_data, z_data, en_data) = Utils.EnergyExtract(str(work_dir),
                                                                energies_filename)
        Fig = go.Figure(layout=go.Layout(title="Visualisation", uirevision='camera'))
        # Fig.update_layout(autosize=True, width=1000, height=1000)
        Molecule = Utils.DrawMolecule(MolData, texture, draw_type)
        for Bond in Molecule['bond_list']: Fig.add_trace(Bond)
        for Atom in Molecule['atom_list']: Fig.add_trace(Atom)
        #
        Fig.update_layout(Utils.GetLayout(None))
        # MinRange = numpy.min(Molecule['geometry'])
        # MaxRange = numpy.max(Molecule['geometry'])
        # Fig.update_layout(Utils.GetRange(MinRange,MaxRange))
        Sphere = go.Surface(x=x_sphere * radius, y=y_sphere * radius,
                            z=z_sphere * radius, customdata=col_sphere,
                            opacity=alphval, surfacecolor=col_sphere,
                            contours={"x": {"show": True}, "y": {"show": True},
                                      "z": {"show": True}}, colorscale='Turbo')
        Energies = go.Scatter3d(x=x_data, y=y_data, z=z_data, mode='markers',
                                marker=dict(size=5, color=en_data,
                                            colorscale='Turbo', opacity=1))
        Contours = go.Scatter3d(x=x_contour, y=y_contour, z=z_contour,
                                mode='markers',
                                marker=dict(size=3, color='black'),
                                hovertext=val_contour)
        Fig.add_trace(Energies)
        if interp == "Sphere" or "Contours": Fig.add_trace(Sphere)
        if interp == "Contours":  Fig.add_trace(Contours)
        # print(Fig)
        return Fig
    else:
        return go.Figure(layout=go.Layout(title="Visualisation", uirevision='camera'))


@app.callback(Output('graph2', 'figure'),
              Input('Comp_vis', 'n_clicks'),
              State('alphval', 'value'),
              State('interp', 'value'),
              State('texture', 'value'),
              State('draw_type', 'value'),
              State('save_location','value'),
              State('sphere-radius','value'),
              State('work_dir1','value'),
)
def Visualise_Interp(n_clicks,alphval,interp,texture,draw_type,work_dir,radius,crd_dir):
    if n_clicks > 0:
        if str(work_dir).endswith("/") == False:
            work_dir = str(work_dir) + "/"
        if str(crd_dir).endswith("/") == False:
            crd_dir = str(crd_dir) + "/"
        energies_filename = "Energy.xyzc"
        MolData = numpy.genfromtxt(str(crd_dir) + "coords_formatted.txt")
        sphdat=json.load(open(str(work_dir)+"interpolation_comparison.dict"))
        x_sphere = numpy.array(sphdat["x_sphere"])
        y_sphere = numpy.array(sphdat["y_sphere"])
        z_sphere = numpy.array(sphdat["z_sphere"])
        col_sphere = numpy.array(sphdat["col_sphere"])
        x_contour = sphdat["x_contour"]
        y_contour = sphdat["y_contour"]
        z_contour = sphdat["z_contour"]
        val_contour = sphdat["val_contour"]
        (x_data, y_data, z_data, en_data) = Utils.EnergyExtract(str(work_dir),
                                                                energies_filename)
        Fig = go.Figure(layout=go.Layout(title="Visualisation", uirevision='camera'))
        # Fig.update_layout(autosize=True, width=1000, height=1000)
        Molecule = Utils.DrawMolecule(MolData, texture, draw_type)
        for Bond in Molecule['bond_list']: Fig.add_trace(Bond)
        for Atom in Molecule['atom_list']: Fig.add_trace(Atom)
        #
        Fig.update_layout(Utils.GetLayout(None))
        # MinRange = numpy.min(Molecule['geometry'])
        # MaxRange = numpy.max(Molecule['geometry'])
        # Fig.update_layout(Utils.GetRange(MinRange,MaxRange))
        Sphere = go.Surface(x=x_sphere * radius, y=y_sphere * radius,
                            z=z_sphere * radius, customdata=col_sphere,
                            opacity=alphval, surfacecolor=col_sphere,
                            contours={"x": {"show": True}, "y": {"show": True},
                                      "z": {"show": True}}, colorscale='Turbo')
        Energies = go.Scatter3d(x=x_data, y=y_data, z=z_data, mode='markers',
                                marker=dict(size=5, color=en_data,
                                            colorscale='Turbo', opacity=1))
        # Contours = go.Scatter3d(x=x_contour, y=y_contour, z=z_contour,
        #                         mode='markers',
        #                         marker=dict(size=3, color='black'),
        #                         hovertext=val_contour)
        Fig.add_trace(Energies)
        if interp == "Sphere" or "Contours": Fig.add_trace(Sphere)
        # if interp == "Contours":  Fig.add_trace(Contours)
        # print(Fig)
        return Fig
    else:
        return go.Figure(layout=go.Layout(title="Visualisation", uirevision='camera'))

def Open_Browser():
    webbrowser.open("0.0.0.0:8050")


if __name__ == '__main__':
    Timer(1, Open_Browser).start();
    app.run_server(debug=True, use_reloader=False, host='0.0.0.0', )


