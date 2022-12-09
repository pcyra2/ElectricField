# Welcome to the Spherical charge electric field script

## Install & Configure:
1. Download source files from GIT
2. Configure conda environment using: \
  conda env create -f CondaEnv.yml
3. Activate conda environment using: \
  conda activate ElectricField
4. Edit the server.config file if not using UoN clusters.

## Run single calculation
1. Run script using: \
  python SphericalWebSRV.py \
   (If on mac, also open a browser window and type 127.0.0.1:8050)
2. Edit user variables
3. Run Calculation
4. Run Interpolation
5. Run Visualisation

## Run comparison
1. Ensure single calculations are run on both systems to be compared
2. Ensure the same charge points are used (ChargeP and ChargeN)
3. Set DATASET1 as the folder with coordinates to be visualised and initial charge data
4. Set DATASET2 as the folder with the comparison charge data
5. Set SAVE LOCATION as the folder to save interpolation cache file. (Usually the parent directory)
6. Run Visualisation (uses the same visualisation parameters from initial visualisation)



# *** Jupyter version now depreciated ***

## Install and general information:

Install Jupyter Notebook if you havent already
Code is dependant on conda packages: numpy, plotly, dash, os, subprocess, time, tqdm.

For cluster/server calculations, passwordless SSH is required so that python can communicate with the slurm submit node.

If on the University of Nottingham network, a version of the visualization may be running on duip76109.nottingham.ac.uk:8050 however this is not always true.

### Usage:

1. Create a working directory on your local machine and your server
2. Place your coordinate file (.xyz) in your working directory and decide on how it should be centered.
3. Open the Jupyter Notebook called "Spherical.ipynb" and set all variables/parameters
4. Run the full script using the run all command
5. Load up the visualization using the links printed


### Running Jupyter notebook via ssh:

Create a ssh tunnel to the port 8080 on the machine 

{ssh -N -L 8080:localhost:8080 XXX@YYYY.nottingham.ac.uk}


In a new terminal, run jupyter notebook from the folder containing the notebook

{jupyter notebook --no-browser --port=8080}


On your browser on your local machine, connect to the notebook through the link:

{localhost:8080/notebooks}


To run the visualisation, you also need to tunnel to port 8050 on the machine:

{ssh -N -L 8080:localhost:8050 XXX@YYYY.nottingham.ac.uk}
