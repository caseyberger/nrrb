#! usr/bin/env python
import os, string, h5py
import pandas as pd
import numpy as np


def get_io_files(work_dir):
    all_files = os.listdir(work_dir)
    input_filename = None
    output_filename = None
    for filename in all_files:        
        if filename.startswith("output"):
            output_filename = filename
        elif filename.startswith("input"):
            input_filename = filename
        else:
            pass
    return input_filename, output_filename

def extract_parameters(work_dir,input_filename):
    #takes input file generated by AMReX and reads from it the various parameters. 
    #returns a dictionary of all the values
    input_file = work_dir+'/'+input_filename
    fo = open(input_file, "r")
    parameters = dict()
    for line in fo.readlines():
        #print(line)
        if line.startswith("nsteps"):
            temp = line.split()
            parameters["nL"] = int(temp[-1])
        elif line.startswith("n_cell"):
            temp = line.split()
            parameters["Nx"] = int(temp[2])
            parameters["Nt"] = int(temp[-1])
        elif line.startswith("nrrb.l"):
            temp = line.split()
            parameters["lambda"] = float(temp[-1])
        elif line.startswith("nrrb.mu"):
            temp = line.split()
            parameters["mu"] = float(temp[-1])
        elif line.startswith("nrrb.m"):
            temp = line.split()
            parameters["m"] = float(temp[-1])
        elif line.startswith("nrrb.w"):
            temp = line.split()
            if temp[0] == "nrrb.w_t":
                parameters["wtr"] = float(temp[-1])
            else:
                parameters["wz"] = float(temp[-1])
        elif line.startswith("nrrb.dtau"):
            temp = line.split()
            parameters["dt"] = float(temp[-1])
        elif line.startswith("nrrb.eps"):
            temp = line.split()
            parameters["eps"] = float(temp[-1])
    parameters["beta"] = parameters["dt"]*parameters["Nt"]
    parameters["dim"] = float(input_filename.split("_")[1][-1])
    fo.close()
    return parameters


def average_new_data(csv_filename, data_dir, obs_list):
    df_old = pd.read_csv(csv_filename)
    df = df_old.drop('tL', axis = 1)
    params_list = ['dim','Nx','Nt','dt','beta','nL','eps','m','mu','wtr','wz','lambda']
    df_params = df[params_list] #takes only the parameters from the original data file
    df_params=df_params.sort_index(axis=1)
    file_header = "nrrb"
    count = 0
    subdirs = os.listdir(data_dir)
    for subdir in subdirs:
        if subdir.startswith(file_header):
            input_filename, output_filename = get_io_files(data_dir+subdir)
            p = extract_parameters(data_dir+subdir,input_filename)
            matches = df_params[(df_params==p).all(axis=1)]
            if len(matches.index.values)==0:
                f = h5py.File(data_dir+subdir+"/"+output_filename)
                df_temp = pd.DataFrame(p,index=[0])
                for obs in obs_list:
                    Re_obs = f["Observables"][obs]["Re"]
                    Im_obs = f["Observables"][obs]["Im"]
                    therm_step = int(0.2*len(Re_obs))
                    df_temp[obs+" (Re)"] = np.mean(Re_obs[therm_step:])
                    df_temp[obs+" Err (Re)"] = np.std(Re_obs[therm_step:])
                    df_temp[obs+" (Im)"] = np.mean(Im_obs[therm_step:])
                    df_temp[obs+" Err (Im)"] = np.std(Im_obs[therm_step:])
                df.append(df_temp,sort=False,ignore_index = True)
                count = count + 1
    tL = df["eps"]*df["nL"]
    df["tL"] = tL
    print(str(count)+" new data points added to data frame")
    return df


def average_data(data_dir, obs_list):
    subdirs = os.listdir(data_dir)
    params_list = ['dim','Nx','Nt','dt','beta','nL','eps','m','mu','wtr','wz','lambda']
    file_header = "nrrb"
    count = 0
    for subdir in subdirs:
        if subdir.startswith(file_header):
            input_filename, output_filename = get_io_files(data_dir+subdir)
            p = extract_parameters(data_dir+subdir,input_filename)
            df_temp = pd.DataFrame(p,index=[0])
            f = h5py.File(data_dir+subdir+"/"+output_filename)
            for obs in obs_list:
                Re_obs = f["Observables"][obs]["Re"]
                Im_obs = f["Observables"][obs]["Im"]
                therm_step = int(0.2*len(Re_obs))
                df_temp[obs+" (Re)"] = np.mean(Re_obs[therm_step:])
                df_temp[obs+" Err (Re)"] = np.std(Re_obs[therm_step:])
                df_temp[obs+" (Im)"] = np.mean(Im_obs[therm_step:])
                df_temp[obs+" Err (Im)"] = np.std(Im_obs[therm_step:])
            if count == 0:
                df = df_temp
            else:
                df = df.append(df_temp,sort=False,ignore_index = True)
            count = count + 1
    tL = df["eps"]*df["nL"]
    df["tL"] = tL
    return df

