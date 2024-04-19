#!/usr/bin/env python3

import sys
import matplotlib
import h5py
matplotlib.use('TkAgg')
import numpy as np
import random
import torch
import joblib as jl
from cebra import CEBRA
           
def main(modelParams_filename):
    #modelParams
    with h5py.File(modelParams_filename, 'r') as hdf:
        model_params = {
            'model_architecture': hdf.attrs['model_architecture'],
            'batch_size': int(hdf.attrs['batch_size']),
            'learning_rate': float(hdf.attrs['learning_rate']),
            'temperature': int(hdf.attrs['temperature']),
            'output_dimension': int(hdf.attrs['output_dimension']),
            'max_iterations': int(hdf.attrs['max_iterations']),
            'distance': hdf.attrs['distance'],
            'conditional': hdf.attrs['conditional'],
            'time_offsets': int(hdf.attrs['time_offsets']),
            'hybrid': hdf.attrs['hybrid'],
            'verbose': hdf.attrs['verbose']
        }
        group_field         = hdf.attrs['group_field']
        behavior_field      = hdf.attrs['behavior_field']
        neural_field        = hdf.attrs['neural_field']
        model_filename      = hdf.attrs['model_filename']
        neural_filename     = hdf.attrs['neural_filename']
        behavior_filename   = hdf.attrs['behavior_filename']
        seed                = int(hdf.attrs['seed'])
        maxI                = int(hdf.attrs['maxI'])

    # seed 
    if seed!=0:
        random.seed(seed)                                   # random
        np.random.seed(random.randint(1,maxI))              # numpy
        # seed PyTorch
        torch.manual_seed(random.randint(1,maxI))
        torch.cuda.manual_seed_all(random.randint(1,maxI))  # multi-GPU
        if torch.backends.cudnn.enabled:
            torch.backends.cudnn.deterministic  = True      # may reduce performance
            torch.backends.cudnn.benchmark      = False

    # load behavior
    with h5py.File(behavior_filename, 'r') as hdf:
        behavior_data  = hdf[group_field + '/' +  behavior_field][:]
    # load neural 
    with h5py.File(neural_filename, 'r') as hdf:
        neural_data    = hdf[group_field + '/' +  neural_field][:]

    # fit model
    cebra_model = CEBRA(**model_params)
    cebra_model.fit(neural_data,behavior_data)
 
    # save model
    torch.save(cebra_model, model_filename)
    print(f"Model saved at {model_filename}")

    return
    
if __name__=="__main__":
    modelParams_filename = sys.argv[1]
    main(modelParams_filename)