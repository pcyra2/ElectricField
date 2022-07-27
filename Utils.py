import numpy as numpy
import plotly.graph_objects as go
import os
import subprocess
# from subprocess import DEVNULL, STDOUT, PIPE, Popen
import time
from tqdm import tnrange, tqdm_notebook

Ang2Bohr = 1.889725989
Bohr2Ang = 0.529177249
Amu2SiKg = 1.660539040e-27

Gaussian_module = "gaussian-uon/avx2/g16"
QChem_HPC_module = "qchem-uon/foss-2020a/5.3"
QChem_CHEM_module = "qchem/qchem-5.4.2-openmp"


CovRadii = {'H' : 0.32, 'He': 0.46, 'Li': 1.33, 'Be': 1.02, 'B' : 0.85, 'C' : 0.75,
            'N' : 0.71, 'O' : 0.63, 'F' : 0.64, 'Ne': 0.67, 'Na': 1.55, 'Mg': 1.39,
            'Al': 1.26, 'Si': 1.16, 'P' : 1.11, 'S' : 1.03, 'Cl': 0.99, 'Ar': 0.96,
            'K' : 1.96, 'Ca': 1.71, 'Sc': 1.48, 'Ti': 1.36, 'V' : 1.34, 'Cr': 1.22,
            'Mn': 1.19, 'Fe': 1.16, 'Co': 1.11, 'Ni': 1.10, 'Cu': 1.12, 'Zn': 1.18,
            'Ga': 1.24, 'Ge': 1.21, 'As': 1.21, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.17,
            'Rb': 2.10, 'Sr': 1.85, 'Y' : 1.63, 'Zr': 1.54, 'Nb': 1.47, 'Mo': 1.38,
            'Tc': 1.28, 'Ru': 1.25, 'Rh': 1.25, 'Pd': 1.20, 'Ag': 1.28, 'Cd': 1.36,
            'In': 1.42, 'Sn': 1.40, 'Sb': 1.40, 'Te': 1.36, 'I' : 1.33, 'Xe': 1.31,
            'Cs': 2.32, 'Ba': 1.96, 'La': 1.80, 'Ce': 1.63, 'Pr': 1.76, 'Nd': 1.74,
            'Pm': 1.73, 'Sm': 1.72, 'Eu': 1.68, 'Gd': 1.69, 'Tb': 1.68, 'Dy': 1.67,
            'Ho': 1.66, 'Er': 1.65, 'Tm': 1.64, 'Yb': 1.70, 'Lu': 1.62, 'Hf': 1.52,
            'Ta': 1.46, 'W' : 1.37, 'Re': 1.31, 'Os': 1.29, 'Ir': 1.22, 'Pt': 1.23,
            'Au': 1.24, 'Hg': 1.33, 'Tl': 1.44, 'Pb': 1.44, 'Bi': 1.51, 'Po': 1.45,
            'At': 1.47, 'Rn': 1.42, 'Fr': 2.23, 'Ra': 2.01, 'Ac': 1.86, 'Th': 1.75,
            'Pa': 1.69, 'U' : 1.70, 'Np': 1.71, 'Pu': 1.72, 'Am': 1.66, 'Cm': 1.66,
            'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65, 'Fm': 1.67, 'Md': 1.73, 'No': 1.76,
            'Lr': 1.61, 'Rf': 1.57, 'Db': 1.49, 'Sg': 1.43, 'Bh': 1.41, 'Hs': 1.34,
            'Mt': 1.29, 'Ds': 1.28, 'Rg': 1.21, 'Gh': 0.00}

Mass = [0.0, 1.00782503207, 4.00260325415, 7.016004548, 9.012182201, 11.009305406,
        12, 14.00307400478, 15.99491461956, 18.998403224, 19.99244017542,
        22.98976928087, 23.985041699, 26.981538627, 27.97692653246, 30.973761629,
        31.972070999, 34.968852682, 39.96238312251, 38.963706679, 39.962590983,
        44.955911909, 47.947946281, 50.943959507, 51.940507472, 54.938045141,
        55.934937475, 58.933195048, 57.935342907, 62.929597474, 63.929142222,
        68.925573587, 73.921177767, 74.921596478, 79.916521271, 78.918337087,
        85.910610729, 84.911789737, 87.905612124, 88.905848295, 89.904704416,
        92.906378058, 97.905408169, 98.906254747, 101.904349312, 102.905504292,
        105.903485715, 106.90509682, 113.90335854, 114.903878484, 119.902194676,
        120.903815686, 129.906224399, 126.904472681, 131.904153457, 132.905451932,
        137.905247237, 138.906353267, 139.905438706, 140.907652769, 141.907723297,
        144.912749023, 151.919732425, 152.921230339, 157.924103912, 158.925346757,
        163.929174751, 164.93032207, 165.930293061, 168.93421325, 173.938862089,
        174.940771819, 179.946549953, 180.947995763, 183.950931188, 186.955753109,
        191.96148069, 192.96292643, 194.964791134, 196.966568662, 201.970643011,
        204.974427541, 207.976652071, 208.980398734, 208.982430435, 210.987496271,
        222.017577738, 222.01755173, 228.031070292, 227.027752127, 232.038055325,
        231.03588399, 238.050788247, 237.048173444, 242.058742611, 243.06138108,
        247.07035354, 247.07030708, 251.079586788, 252.082978512, 257.095104724,
        258.098431319, 255.093241131, 260.105504, 263.112547, 255.107398,
        259.114500, 262.122892, 263.128558, 265.136151, 281.162061, 272.153615]

Symbol = ["X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al",
          "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co",
          "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb",
          "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs",
          "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
          "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
          "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk",
          "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg"]

Symb2Mass = dict(zip(Symbol,Mass))
Atno2Mass = dict(zip(list(range(11)),Mass))
Atno2Symb = dict(zip(list(range(11)),Symbol))
Symb2Atno = dict(zip(Symbol,list(range(11))))


AtomColours = {"H" : ["rgba(255, 255, 255, 1.0)", 1],
               "He": ["rgba(255, 191, 203, 1.0)", 2],
               "Li": ["rgba(178, 034, 033, 1.0)", 3],
               "Be": ["rgba(255, 019, 147, 1.0)", 4],
               "B" : ["rgba(000, 251, 002, 1.0)", 5],
               "C" : ["rgba(198, 198, 198, 1.0)", 6],
               "N" : ["rgba(140, 142, 251, 1.0)", 7],
               "O" : ["rgba(236,   0,   0, 1.0)", 8],
               "F" : ["rgba(217, 164,  32, 1.0)", 9],
               "Ne": ["rgba(253,  20, 146, 1.0)", 10],
               "Na": ["rgba(  0,   0, 253, 1.0)", 11],
               "Mg": ["rgba( 30, 124,  30, 1.0)", 12],
               "Al": ["rgba(128, 128, 144, 1.0)", 13],
               "Si": ["rgba(211, 159,  31, 1.0)", 14],
               "P" : ["rgba(241, 156,   0, 1.0)", 15],
               "S" : ["rgba(250, 196,  49, 1.0)", 16],
               "Cl": ["rgba(  0, 241,   0, 1.0)", 17],
               "Ar": ["rgba(251, 20,  145, 1.0)", 18],
               "K" : ["rgba(255, 255, 255, 1.0)", 19],
               "Ca": ["rgba(255, 191, 203, 1.0)", 20],
               "Sc": ["rgba(178,  34,  33, 1.0)", 21],
               "Ti": ["rgba(255,  19, 147, 1.0)", 22],
               "V" : ["rgba(  0, 251,   2, 1.0)", 23],
               "Cr": ["rgba(198, 198, 198, 1.0)", 24],
               "Mn": ["rgba(140, 142, 251, 1.0)", 25],
               "Fe": ["rgba(236,   0,   0, 1.0)", 26],
               "Co": ["rgba(217, 164,  32, 1.0)", 27],
               "Ni": ["rgba(253,  20, 146, 1.0)", 28],
               "Cu": ["rgba(  0,   0, 253, 1.0)", 29],
               "Zn": ["rgba( 30, 124,  30, 1.0)", 30],
               "Ga": ["rgba(128, 128, 144, 1.0)", 31],
               "Ge": ["rgba(211, 159,  31, 1.0)", 32],
               "As": ["rgba(241, 156,   0, 1.0)", 33],
               "Se": ["rgba(250, 196,  49, 1.0)", 34],
               "Br": ["rgba(  0, 241,   0, 1.0)", 35],
               "Kr": ["rgba(251, 20,  145, 1.0)", 36],
               "Rb": ["rgba(255, 255, 255, 1.0)", 37],
               "Sr": ["rgba(255, 191, 203, 1.0)", 38],
               "Y" : ["rgba(178,  34,  33, 1.0)", 39],
               "Zr": ["rgba(255,  19, 147, 1.0)", 40],
               "Nb": ["rgba(  0, 251,   2, 1.0)", 41],
               "Mo": ["rgba(198, 198, 198, 1.0)", 42],
               "Tc": ["rgba(140, 142, 251, 1.0)", 43],
               "Ru": ["rgba(236,   0,   0, 1.0)", 44],
               "Rh": ["rgba(217, 164,  32, 1.0)", 45],
               "Pd": ["rgba(253,  20, 146, 1.0)", 46],
               "Ag": ["rgba(  0,   0, 253, 1.0)", 47],
               "Cd": ["rgba( 30, 124, 30,  1.0)", 48],
               "In": ["rgba(128, 128, 144, 1.0)", 49],
               "Sn": ["rgba(211, 159,  31, 1.0)", 50],
               "Sb": ["rgba(241, 156,   0, 1.0)", 51],
               "Te": ["rgba(250, 196,  49, 1.0)", 52],
               "I" : ["rgba(  0, 241,   0, 1.0)", 53],
               "Xe": ["rgba(251,  20, 145, 1.0)", 54]}

SurfaceStyles = {'matte': {'ambient'             : 0.60,
                           'diffuse'             : 0.35,
                           'fresnel'             : 0.05,
                           'specular'            : 0.03,
                           'roughness'           : 0.05,
                           'facenormalsepsilon'  : 1e-15,
                           'vertexnormalsepsilon': 1e-15},
                 'shiny': {'ambient'             : 0.30,
                           'diffuse'             : 0.85,
                           'fresnel'             : 0.10,
                           'specular'            : 0.70,
                           'roughness'           : 0.05,
                           'facenormalsepsilon'  : 1e-15,
                           'vertexnormalsepsilon': 1e-15},
                 'orbs' : {'ambient'             : 0.30,
                           'diffuse'             : 0.60,
                           'fresnel'             : 0.10,
                           'specular'            : 0.70,
                           'roughness'           : 0.90,
                           'facenormalsepsilon'  : 1e-15,
                           'vertexnormalsepsilon': 1e-15},
                 'glass': {'ambient'             : 0.30,
                           'diffuse'             : 0.10,
                           'fresnel'             : 0.45,
                           'specular'            : 0.70,
                           'roughness'           : 0.10,
                           'facenormalsepsilon'  : 1e-15,
                           'vertexnormalsepsilon': 1e-15}}

def CoordGet(work_dir, file):
    with open(work_dir+file,'r') as f:
        data = f.readlines()
    x,y,z,at,=[],[],[],[]
    for line in data:
        words = line.split()
        if len(words) == 4:
            at.append(str(words[0]))
            x.append(float(words[1]))
            y.append(float(words[2]))
            z.append(float(words[3]))
    return (at, x, y, z)

def CoordCheck(at,x,y,z,recenter,center_type,data):
    print(x)
    x_center = average(x)
    y_center = average(y)
    z_center = average(z)
    print(f"Center of coordinates is {x_center} {y_center} {z_center}")
    if round(x_center) != 0 or round(y_center) != 0 or round(z_center) != 0:
        print(f"Warning, your system is not centered about the origin using simple metrics")
        if recenter == True:
            if center_type == "simple":
                (x,y,z) = Recenter(x,y,z,x_center,y_center,z_center)
            elif center_type == "mass":
                print("Mass weighted centering is not currently tested, Check centering correctly!")
                newat = []
                for i in range(len(at)):
                    for j in range(len(Symbol)):
                        if at[i] == Symbol[j]: newat.append(float(j))
                x_center = average(x*newat)
                y_center = average(y*newat)
                z_center = average(z*newat)
                (x,y,z) = Recenter(x,y,z,x_center,y_center,z_center)
            elif center_type == "atom":
                cx = x[data[0]-1]
                cy = y[data[0]-1]
                cz = z[data[0]-1]
                (x,y,z) = Recenter(x,y,z,cx,cy,cz)
            elif center_type == "translate":
               (x,y,z) = Recenter(x,y,z,data[0],data[1],data[2]) 
    return(x, y, z)

def Recenter(x,y,z,cx,cy,cz):
    for i in range(len(x)):
        x[i] = x[i] - cx
        y[i] = y[i] - cy
        z[i] = z[i] - cz
    return(x,y,z)

def SPRun(work_dir, sp_memory, sp_threads, NumPoints, hostdir, server="loginchem01.nottingham.ac.uk",HPC="Y",QMPACKAGE="QChem"):
    if HPC == "Y":
        print("Forming HPC Array JOB")
        with open(work_dir+"array_job.sh",'w') as f:
            print("#!/bin/bash",file=f)
            print("RUNLINE=$(cat $ARRAY_TASKFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)",file=f)
            print("eval $RUNLINE",file=f)
        if QMPACKAGE == "GAUSSIAN":
            if server == "loginchem01.nottingham.ac.uk":
                print("This cluster does not currently have gaussian installed")
                return ValueError
            elif server == "login001.augusta.nottingham.ac.uk":
                module_file = Gaussian_module
                scratch_loc = "~/SCRATCH"
                partition = "defq"
            with open(work_dir+"sub.sh",'w') as f:
                print("#!/bin/bash",file=f)
                print("export ARRAY_JOBFILE=array_job.sh",file=f)
                print("export ARRAY_TASKFILE=jobs.tmp",file=f)
                print("export ARRAY_NTASKS=$(cat $ARRAY_TASKFILE | wc -l)",file=f)
                print("module load "+str(module_file),file=f)
                print("export "+str(scratch_loc),file=f)
                print("CAL000=$(sbatch -J \"SinglePoint\" --mem="+str(sp_memory)+"G --nodes=1 --cpus-per-task="+str(sp_threads) +" -p \"defq\" -t 20 --get-user-env --parsable --array=1-$ARRAY_NTASKS $ARRAY_JOBFILE)",file=f)
                print("echo $CAL000",file=f)
            with open(work_dir+"jobs.tmp",'w') as f:
                for i in range(0,NumPoints):
                    text = "cd ./"+str(i)+" ; g16 SP.com"
                    print(str(text),file=f)
        if QMPACKAGE == "QChem":
            if server == "loginchem01.nottingham.ac.uk":
                module_file = QChem_CHEM_module
                scratch_loc = "~/PERSONAL_SCRATCH"
                partition = "dumbledore"
            elif server == "login001.augusta.nottingham.ac.uk":
                module_file = QChem_HPC_module
                scratch_loc = "~/SCRATCH"
                partition = "defq"
            with open(work_dir + "sub.sh", 'w') as f:
                print("#!/bin/bash", file=f)
                print("export ARRAY_JOBFILE=array_job.sh", file=f)
                print("export ARRAY_TASKFILE=jobs.tmp", file=f)
                print("export ARRAY_NTASKS=$(cat $ARRAY_TASKFILE | wc -l)", file=f)
                print("module load " + str(module_file), file=f)
                print("export QCSCRATCH=" + str(scratch_loc), file=f)
                print("CAL000=$(sbatch -J \"SinglePoint\" --mem=" + str(
                    sp_memory) + "G --nodes=1 --cpus-per-task=" + str(
                    sp_threads) + " -p \""+partition+"\" -t 20 --get-user-env --parsable --array=1-$ARRAY_NTASKS $ARRAY_JOBFILE)",
                      file=f)
                print("echo $CAL000", file=f)
            with open(work_dir + "jobs.tmp", 'w') as f:
                for i in range(0,NumPoints):
                    text = "cd ./" + str(i) + " ; qchem -nt " + str(sp_threads) + " SP.inp SP.out"
                    print(str(text), file=f)
        print("Submitting HPC ARRAY JOB")
        Submit(server,hostdir,work_dir)
        print("HPC ArrayJob submitted")
        FinishCheck(server,hostdir,work_dir,NumPoints,sp_threads)
    else:
        print("Starting Calculations locally")
        line="seq 0 "+str(NumPoints)+" | parallel -j "+str(sp_threads)+ " g16 ./{}/SP.com"
        os.system(line)
        print("Single Point Calculations Complete")

def average(lst):
    if len(lst)==0:
        print("coordinates are not loaded, breaking")
        raise IndexError
    return sum(lst)/len(lst)

def FormatCoord(at, x, y, z,work_dir,file="coords_formatted.txt"):
    newat = []
    for i in range(len(at)):
        for j in range(len(Symbol)):
            if at[i] == Symbol[j]: newat.append(str(j))
    with open(work_dir+file,'w') as f:
        for i in range(len(at)):
            print(str(x[i])+"\t"+str(y[i])+"\t"+str(z[i])+"\t"+str(newat[i]),file=f)

def sample_spherical(npoints):          ###Dont call (Script uses this not you)
    vec = numpy.random.randn(3, npoints)
    vec /= numpy.linalg.norm(vec, axis=0)
    return vec

def RandomCoords(work_dir,num_points,sphere_radius,file="ChargeP.xyz"):
    phi = numpy.linspace(0, numpy.pi, 20)
    theta = numpy.linspace(0, 2 * numpy.pi, 40)
    x = numpy.outer(numpy.sin(theta), numpy.cos(phi)) #+ x0
    y = numpy.outer(numpy.sin(theta), numpy.sin(phi)) #+ y0
    z = numpy.outer(numpy.cos(theta), numpy.ones_like(phi)) #+ z0
    xi, yi, zi = sample_spherical(num_points)
    pxi =sphere_radius * xi #+ x0
    pyi =sphere_radius * yi #+ y0
    pzi =sphere_radius * zi #+ z0
    nxi =0 - sphere_radius*xi
    nyi =0 - sphere_radius*yi
    nzi =0 - sphere_radius*zi
    with open(str(work_dir)+file,'w') as f:
        for i in range(len(pxi)):
            print(str(pxi[i])+"\t"+str(pyi[i])+"\t"+str(pzi[i]),file=f)
    with open(str(work_dir)+file.replace("P","N"),'w') as f:
        for i in range(len(nxi)):
            print(str(nxi[i])+"\t"+str(nyi[i])+"\t"+str(nzi[i]),file=f)

def GenSP(opt_sys_at, opt_sys_x, opt_sys_y, opt_sys_z, work_dir, point_charge_value, sp_memory,
          sp_threads, sp_method, sp_basis, sys_charge, sys_spin, Dispersion = "D3",QMPACKAGE  = "QChem", charge_file = "ChargeP.xyz"):       ##generates Single point input files
    sys_center_x = average(opt_sys_x)
    sys_center_y = average(opt_sys_y)
    sys_center_z = average(opt_sys_z)
    if len(opt_sys_at)==0:
        print("coordinates are not loaded, breaking")
        raise IndexError
    cp=[]
    cn=[]
    with open(str(work_dir)+charge_file) as chargep:
        rline = chargep.readlines()
    for line in rline:
        cp.append(str(line.replace("\n","\t")))
    with open(str(work_dir)+charge_file.replace("P","N")) as chargen:
        rline = chargen.readlines()
    for line in rline:
        cn.append(str(line.replace("\n","\t")))
    for i in range(0,len(cp)):
        directory = str(work_dir) + str(i)
        try:
            os.mkdir(directory)
        except FileExistsError:
            print("directory "+ str(i) + " already exists")
        pc = str(cp[i]) + " " + str(point_charge_value)
        nc = str(cn[i]) + " -" + str(point_charge_value)
        if QMPACKAGE == "Gaussian":
            sp_file = str(directory) + "/SP.com"
            with open(sp_file,'w') as f:
                print("%mem=" + str(sp_memory) + "GB", file=f)
                print("%nprocshared=" + str(sp_threads), file=f)
                print(" #p OPT " + str(sp_method + "/" + str(sp_basis)), file=f)
                print("", file=f)
                print(str(inp) + " opt", file=f)
                print("", file=f)
                print(str(sys_charge) + " " + str(sys_spin), file=f)
                for j in range(len(opt_sys_at)):
                    line = str(opt_sys_at[j]) + "\t" + str(opt_sys_x[j]-sys_center_x) + "\t" + str(opt_sys_y[j]-sys_center_y) + "\t" + str(opt_sys_z[j]-sys_center_z)
                    print(str(line), file=f)
                print("",file=f)
                print(str(pc),file=f)
                print(str(nc),file=f)
                print("", file=f)
                print("", file=f)
        elif QMPACKAGE == "QChem":
            sp_file = str(directory) + "/SP.inp"
            with open(sp_file,'w') as f:
                print("$molecule",file=f)
                print(str(sys_charge) + " " + str(sys_spin), file=f)
                for j in range(len(opt_sys_at)):
                    line = str(opt_sys_at[j]) + "\t" + str(opt_sys_x[j]-sys_center_x) + "\t" + str(opt_sys_y[j]-sys_center_y) + "\t" + str(opt_sys_z[j]-sys_center_z)
                    print(str(line), file=f)
                print("$end", file=f)
                print("",file=f)
                print("$rem", file=f)
                print("JOBTYPE = SP", file=f)
                print("METHOD = " + str(sp_method), file=f)
                print("BASIS = " + str(sp_basis), file=f)
                print("DFT_D = "+ str(Dispersion),file=f)
                print("SYM_IGNORE = TRUE",file=f)
                print("SYMMETRY = FALSE",file=f)
                print("$end", file=f)
                print("", file=f)
                print("$external_charges", file=f)
                print(str(pc),file=f)
                print(str(nc),file=f)
                print("$end", file=f)

def Submit(hostname,hostdir,workdir):
    print("syncing data")
    print("Hostname = "+ hostname + ", Workdir = "+ workdir+", Serverdir = "+hostdir)
    # rsync = subprocess.Popen(['rsync','-azP', workdir+"*", hostname+":"+hostdir+"."], stderr=subprocess.STDOUT, stdout=subprocess.PIPE,shell=True).communicate()[0]
    print("sync initiated, waiting for sync to complete")
    sync = subprocess.run(['rsync', '-azP',  workdir, hostname+':'+hostdir],capture_output=True, text=True)
    # sync = Popen(['rsync', '-azP', workdir, hostname+':'+hostdir+'.',],stdout=PIPE,stdin=PIPE)#.communicate()[0]
#     print("sync initiated, waiting for sync to complete")
    # sync.wait()
    # status  = sync#.decode().split("\n")
    # print(status)
    print("Data synced, running jobs")
    #run = Popen(['ssh' , hostname, 'cd '+hostdir+' ;sbatch sub.sh'], stdout=PIPE, stderr=PIPE, stdin=PIPE)#.communicate()[0]
    run = subprocess.run(['ssh', hostname, 'cd '+ hostdir + ' ; sbatch sub.sh'], capture_output = True, text = True)
    # run.wait()
    # run = run.communicate()[0]
    # print(run)
    status = run.stdout.split()
    # print(status)
    pid = status[3]
    print("job id = "+pid)

def FinishCheck(hostname,hostdir,workdir,njobs,threads):
    print("Waiting for calculations to be completed")
    done = False
    counter=0
    while done == False:
        counter = counter + 1
        print("wait cycle: "+str(counter))
        # time.sleep(10)
        # rsync = subprocess.Popen(['rsync', '-azP', hostname+':'+hostdir+'*',  workdir+'.'], stdout=PIPE, stderr=PIPE, stdin=PIPE)
        print("Waiting for rsync")
        rysync = subprocess.run(['rsync', '-azP', hostname+':'+hostdir,  workdir], capture_output = True, text = True)
        #rsync.wait()
        print("Rsync done")
        # slurmcheck = subprocess.Popen(['ls '+workdir+'slurm*'],stdout=PIPE, stderr=PIPE, stdin=PIPE,shell=True)#.communicate()[0]
        # slurmcheck.wait()
        slurmcheck = subprocess.run(['ls '+ workdir+'slurm*'], shell=True, text=True, capture_output = True)
        print("Counting slurm files")
        # slurmcheck = slurmcheck.communicate()[0]
        # print(slurmcheck.stdout.split("\n"))
        slurms = slurmcheck.stdout.split("\n")
        print("Number of jobs started = "+str(len(slurms)-2)+" out of "+str(njobs))
        if len(slurms) < njobs+1:
            print("Waiting for jobs to run")
            time.sleep(120)
        else:
            try:
                failed = OutputExtract(workdir,njobs,threads)
            except IndexError:
                time.sleep(120)
                print("Running a final check for calculations")
            if failed > 0 :
                rysync = subprocess.run(['rsync', '-azP', hostname+':'+hostdir,  workdir], capture_output = True, text = True)
                print("Final rsync performed")
                failed = OutputExtract(workdir,njobs,threads)
            if failed == 0: done = True
            else:
                print("There are failed calculations, please run them manually or check the outputs.")
                return ValueError
    print("Jobs completed")

def OutputExtract(work_dir,num_points,sp_threads,QMPACKAGE  = "QChem"):
    with open(str(work_dir) + "ChargeP.xyz") as chargep:
        rline = chargep.readlines()
    xyzc=[]
    failed=[]
    fail_tot=0
    if QMPACKAGE == "GAUSSIAN":
        for i in range(0,num_points):
            file_dir=str(work_dir)+str(i)
            #with open(file_dir+inp+".log") as f:
            with open(file_dir+"/SP.log") as f:
                lines=f.readlines()
            for line in lines:
                if "SCF Done" in line:
                    words=line.split()
                    energy=words[4]
                    break
            xyzc.append(str(rline[i].replace("\n","\t"))+str(energy))
    elif QMPACKAGE == "QChem":
        # print("qchem energy")
        for i in range(0,num_points):
            file_dir=str(work_dir)+str(i)
            #with open(file_dir+inp+".log") as f:
            with open(file_dir+"/SP.out") as f:
                lines=f.readlines()
            if (len(lines) < 100):
                print("probably a problem with file number "+ str(i))
            energy="0"
            # if (i == 125):
            #     print(file_dir)
            for line in lines:
                # if (i == 125):
                #     print(line)
                if "Total energy in the " in line:
                    words=line.split()
                    energy=words[8]
                    # print(str(i) + " " + str(energy))
                    # print(line)
                    break
            if energy == "0":
                # print(i)
                print("Calculation "+ str(i)+ " has failed.")
                failed.append(int(i))
                fail_tot=fail_tot+1
            xyzc.append(str(rline[i].replace("\n","\t"))+str(energy))
    print("Total number of failed jobs is "+str(fail_tot))
    with open(str(work_dir)+"Energy.xyzc",'w') as f:
        for i in range(len(xyzc)):
            j=i-1
            if j not in failed:
                #print(xyzc[i])
                print(xyzc[i],file=f)
    with open(str(work_dir)+"ReRunJobs.txt",'w') as f:
        for i in failed:
            text = "cd ./" + str(i) + " ; qchem -nt " + str(sp_threads) + " SP.com"
            print(str(text),file=f)
    with open(str(work_dir)+"Failed_Energy.csv",'w') as f:
        print(xyzc[0],file=f)
        for i in failed:
            print(i)
            j=i+1
            print(str(xyzc[j]),file=f)
    return fail_tot

def EnergyExtract(work_dir, energies_filename):
    with open(work_dir + energies_filename,'r') as f:
        data = f.readlines()
    X,Y,Z,C=[],[],[],[]
    for line in data:
        words = line.split()
        X.append(float(words[0]))
        Y.append(float(words[1]))
        Z.append(float(words[2]))
        C.append(float(words[3]))
    return (X,Y,Z,C)

def SphereGen(x, y, z, radius, resolution=10):
    u, v = numpy.mgrid[0:2*numpy.pi:resolution*2j, 0:numpy.pi:resolution*1j]
    X = radius * numpy.cos(u)*numpy.sin(v) + x
    Y = radius * numpy.sin(u)*numpy.sin(v) + y
    Z = radius * numpy.cos(v) + z
    return (X, Y, Z)

def GenColors(x_sphere, y_sphere, z_sphere,x_data, y_data, z_data, en_data,radius):
    ### Generate array for colors
    col = numpy.zeros([len(x_sphere[:,1]),len(x_sphere[1])])
    ### Populate array of colors
    cX,cY,cZ,contVal=[],[],[],[]
    iso_values = set(numpy.around(numpy.arange(round(numpy.max(en_data),3),round(numpy.min(en_data),3),-0.001),4))
    # print(round(numpy.min(en_data),3),round(numpy.max(en_data),3))
    # print(iso_values)
    # for i in range(len(x_sphere[:,1])):
    for i in tqdm_notebook(range(len(x_sphere[:,1])), desc="Color progress"):
        # print("Color: "+str(i/len(x_sphere[:,1])*100)+"%")
        for j in tqdm_notebook(range(len(y_sphere[1])),desc="Step "+ str(i) + " progress", leave = False):
            top = 0
            bottom = 0
            for k in range(len(x_data)):
                dx = (x_sphere[i,j] * radius) * x_data[k]
                dy = (y_sphere[i,j] * radius) * y_data[k]
                dz = (z_sphere[i,j] * radius) * z_data[k]
                dP = numpy.arccos((dx+dy+dz)/(numpy.power(radius,2)))*radius
                top = top + en_data[k] * numpy.exp(-numpy.power(dP,2))
                bottom = bottom + numpy.exp(-numpy.power(dP,2))
            # print(top, bottom)
            cVal = round(top/bottom,4)
            # print(cVal)
            # for k in numpy.arange(round(numpy.min(en_data),3),round(numpy.max(en_data),3)-0.001):
            if cVal in iso_values:
                cX.append(x_sphere[i,j] * radius)
                cY.append(y_sphere[i,j] * radius)
                cZ.append(z_sphere[i,j] * radius)
                contVal.append((round(cVal,4)+456.080)*627.5)
                # print((round(cVal,4)+456.080)*627.5)
            col[i,j] = top/bottom
    print("Coloring done")
    return (col,cX,cY,cZ,contVal)

def GetRange(min_range, max_range, overage=2.5):
    axis = {'autorange': False,
            'range': (min_range * overage, max_range * overage)}

    layout = {'scene': {'xaxis': axis, 'yaxis': axis, 'zaxis': axis},
              'scene_camera': {'up': {'x': 0, 'y': 0, 'z': 1},
                               'center': {'x': 0, 'y': 0, 'z': 0},
                               'eye': {'x': 2.2 / max_range,
                                       'y': 2.2 / max_range,
                                       'z': 2.2 / max_range}}}
    return layout

def GetLayout(figsize=None):
    axis = {'showgrid': False,
            'zeroline': False,
            'showline': False,
            'title': '',
            'ticks': '',
            'showticklabels': False,
            'showbackground': False,
            'showspikes': False}

    layout = {'scene_aspectmode': 'manual',
              'scene_aspectratio': {'x': 1, 'y': 1, 'z': 1},
              'scene_xaxis_showticklabels': False,
              'scene_yaxis_showticklabels': False,
              'scene_zaxis_showticklabels': False,
              'dragmode': 'orbit',
              'template': 'plotly_white',
              'showlegend': False,
              'hovermode': False,
              'margin': {'t': 0, 'l': 0, 'b': 0, 'r': 0},
              'scene': {'xaxis': axis, 'yaxis': axis, 'zaxis': axis}}

    if figsize is not None:
        layout['height'] = figsize[0]
        layout['width'] = figsize[1]

    return layout

def Connectivity(molecule):
    T = []
    n = molecule.shape[0]
    for i in range(n):
        for j in range(i):
            Rij = numpy.linalg.norm(molecule[i, :3] - molecule[j, :3])
            ri = CovRadii[Atno2Symb[int(molecule[i, 3])]]
            rj = CovRadii[Atno2Symb[int(molecule[j, 3])]]
            ithr = 0.7 * Ang2Bohr * (ri + rj)
            if Rij != 0.0 and Rij < ithr:
                T.append((i, j))
    return T

def GetAtoms(geometry, atomic_numbers, symbols, style, surface,
           quality='high'):
    numpyts = {'low': 5, 'medium': 10, 'high': 20}
    trace_list = []
    sphere = Sphere(1, numpyts[quality])
    for atom, xyz in enumerate(geometry):
        if style == 'ball_and_stick':
            reshaped_sphere = sphere * (atomic_numbers[atom] / 30.0 + 0.6) * 0.6
        elif style == 'tubes':
            reshaped_sphere = sphere * 0.3
        elif style == 'wireframe':
            reshaped_sphere = sphere * 0.06
        elif style == 'spacefilling':
            reshaped_sphere = sphere * (atomic_numbers[atom] / 20.0 + 1.5)
        else:
            sys.exit('Molecule style {} not available'.format(style))
        trace_list.append(
            GetAtomMesh(reshaped_sphere, symbols[atom], xyz, surface))
    return trace_list

def GetBonds(geometry, symbols, bonds, style, surface, quality='high'):
    numpyts = {'low': 20, 'medium': 50, 'high': 100}
    trace_list = []
    for i, j in bonds:
        vi = geometry[i]
        vj = geometry[j]
        dij = numpy.linalg.norm(vj - vi)
        R = RotationMatrix(numpy.array([0.0, 0.0, 1.0]), vj - vi)
        if (style == 'ball_and_stick') or (style == 'tubes'):
            r = 0.15, #0.3
        elif style == 'wireframe':
            r = 0.06
        elif style == 'spacefilling':
            return []
        else:
            sys.exit('Molecule style {} not available'.format(style))
        if symbols[i] == symbols[j]:
            cyl = Cylinder(r, numpyts[quality])
            cyl[:, 2] *= dij
            cyl = numpy.dot(R, cyl.transpose()).transpose()
            cyl += vi
            mesh = GetBondMesh(cyl, i, symbols, surface)
            trace_list.append(mesh)
        else:
            cyl = Cylinder(r, numpyts[quality])
            cyl[:, 2] *= (0.5 * dij)
            cyl = numpy.dot(R, cyl.transpose()).transpose()
            cyl_i = cyl + vi
            cyl_j = cyl + 0.5 * (vi + vj)
            mesh = GetBondMesh(cyl_i, i, symbols, surface)
            trace_list.append(mesh)
            mesh = GetBondMesh(cyl_j, j, symbols, surface)
            trace_list.append(mesh)
    return trace_list

def GetAtomMesh(sphere, symbol, xyz, surface):
    mesh = go.Mesh3d({'x': sphere[0] + xyz[0],
                      'y': sphere[1] + xyz[1],
                      'z': sphere[2] + xyz[2],
                      'color': AtomColours[symbol][0],
                      'alphahull': 0,
                      'flatshading': False,
                      'cmin': -7,
                      'lighting': SurfaceStyles[surface],
                      'lightposition': {'x': 100, 'y': 200, 'z': 0}})
    return mesh

def GetBondMesh(cylinder, bond, symbols, surface):
    mesh = go.Mesh3d({'x': cylinder[:, 0],
                      'y': cylinder[:, 1],
                      'z': cylinder[:, 2],
                      'color': AtomColours[symbols[bond][0]][0],
                      'alphahull': 0,
                      'flatshading': True,
                      'cmin': -7,
                      'lighting': SurfaceStyles[surface],
                      'lightposition': {'x': 100, 'y': 200, 'z': 0}})
    return mesh

def Cylinder(r=1.0, n=100):
    phi = numpy.linspace(0, 2.0 * numpy.pi, n)
    x = r * numpy.cos(phi)
    y = r * numpy.sin(phi)
    pts = numpy.zeros((2 * len(phi), 3))
    pts[:len(phi), 0] = x
    pts[len(phi):, 0] = x
    pts[:len(phi), 1] = y
    pts[len(phi):, 1] = y
    pts[:len(phi), 2] = 0.0
    pts[len(phi):, 2] = 1.0
    return pts

def Sphere(r=1.0, n=20):
    phi = numpy.linspace(0, 2.0 * numpy.pi, 2 * n)
    theta = numpy.linspace(-0.5 * numpy.pi, 0.5 * numpy.pi, n)
    phi, theta = numpy.meshgrid(phi[1:], theta)
    pts = numpy.zeros((3, phi.size))
    pts[0] = r * (numpy.cos(theta) * numpy.sin(phi)).flatten()
    pts[1] = r * (numpy.cos(theta) * numpy.cos(phi)).flatten()
    pts[2] = r * (numpy.sin(theta)).flatten()
    return pts

def RotationMatrix(a, b):
    """
        Find the rotation matrix that aligns a to b
    """
    u = a / numpy.linalg.norm(a)
    v = b / numpy.linalg.norm(b)
    c = numpy.cross(u, v)
    d = numpy.dot(u, v)
    s = numpy.linalg.norm(c)
    K = numpy.asarray(
        [[0.0, -c[2], c[1]], [c[2], 0.0, -c[0]], [-c[1], c[0], 0.0]])
    if s <= 1.0e-10:
        R = (1.0 - numpy.linalg.norm(u - v)) * numpy.identity(3)
    else:
        R = numpy.identity(3) + K + numpy.dot(K, K) * (1.0 - d) / (s * s)
    return R

def DrawMolecule(MolData,txture,style):
    Molecule = {}
    Molecule['geometry'] = MolData[:,:3].copy()
    Molecule['symbols'] = [Atno2Symb[int(MolData[i,3])] for i in range(MolData.shape[0])]
    Molecule['atomic_numbers'] = [int(MolData[i,3]) for i in range(MolData.shape[0])]
    Molecule['bonds'] = Connectivity(MolData)
    Molecule['bond_list'] = GetBonds(Molecule['geometry'],Molecule['symbols'],Molecule['bonds'],style,txture)
    Molecule['atom_list'] = GetAtoms(Molecule['geometry'],Molecule['atomic_numbers'],Molecule['symbols'],style,txture)
    return Molecule
