#!/usr/bin/env python3

#import fcca
from FCCA.fcca import LQGComponentsAnalysis as LQGCA
#import dca
import numpy as np
#import base
# Define dimensions
# nTimes      = 100  # Number of rows
# nChannels   = 5  # Number of columns

# Create a random matrix with values between 0 and 1
# X           = np.random.rand(nTimes, nChannels)
# lqgca_model = LQGCA(T=4, d=2, n_init=10).fit(X)
# Projection matrix corresponding to the optimized projection
# print(lqgca_model)
# print(lqgca_model.coef_)

import sys
import h5py
import numpy as np
import random
import torch
bool_vals = {'true': True, 'false': False}           
def main(modelParams_filename):
    #modelParams
    with h5py.File(modelParams_filename, 'r') as hdf:
        # model_params = {
        #     "d": int(hdf.attrs['d']),
        #     "T": int(hdf.attrs['T']),
        #     "chunk_cov_estimate":None,
        #     "dtype": torch.float64,
           
        model_params = {
            "d": int(hdf.attrs['d']),
            "T": int(hdf.attrs['T']),
            "init": hdf.attrs['init'],
            "n_init": int(hdf.attrs['n_init']),
            "stride": int(hdf.attrs['stride']),
            "tol": float(hdf.attrs['tol']),
            "ortho_lambda": float(hdf.attrs['ortho_lambda']),
            "verbose": bool_vals[hdf.attrs['verbose'].lower()],
            "device": hdf.attrs['device'],
            "rng_or_seed": int(hdf.attrs['rng_or_seed']),
        }
        group_field         = hdf.attrs['group_field']
        neural_field        = hdf.attrs['neural_field']
        model_field         = hdf.attrs['model_field']
        neural_filename     = hdf.attrs['neural_filename']
        model_filename      = hdf.attrs['model_filename']
        #seed                = int(hdf.attrs['seed'])
        #maxI                = int(hdf.attrs['maxI'])

    # seed 
    #if seed!=0:
    #    random.seed(seed)                                   # random
    #    np.random.seed(random.randint(1,maxI))              # numpy
        # seed PyTorch
    #    torch.manual_seed(random.randint(1,maxI))
    #    torch.cuda.manual_seed_all(random.randint(1,maxI))  # multi-GPU
    #    if torch.backends.cudnn.enabled:
    #        torch.backends.cudnn.deterministic  = True      # may reduce performance
    #        torch.backends.cudnn.benchmark      = False

    # load neural 
    #  https://stackoverflow.com/questions/21624653/python-created-hdf5-dataset-transposed-in-matlab
    # note: same hdf5 file result transposed 
    # in matlab: Channels x nTimes
    # in python:   nTimes x nChannels 

    with h5py.File(neural_filename, 'r') as hdf:
        neural_data    = hdf[group_field + '/' +  neural_field][:]
    # fit model
    #lqgca_model = LQGCA(d=2, T=4,  
    #             init="random_ortho", n_init=100, stride=1,
    #             chunk_cov_estimate=None, tol=1e-6, ortho_lambda=10., verbose=False,
    #             device="cpu", dtype=torch.float64, rng_or_seed=None).fit(X)
    # lqgca_model = LQGCA(** model_params).fit(np.transpose(neural_data))
    print(neural_data.shape)
    # input('neural data size')
    lqgca_model = LQGCA(** model_params).fit(neural_data)
    Wfcaa=lqgca_model.coef_;
    print(Wfcaa.shape)
    # input('Wfcaa size')

    # save data in manifold_filename 
    with h5py.File(model_filename, 'w') as out_hdf:
        group = out_hdf.create_group(group_field)
        group.create_dataset(model_field, data=Wfcaa, compression='gzip', compression_opts=9)
        print(f"{model_field} saved in {model_filename}")   

    # # save model
    # torch.save(Wfcaa, model_filename)
    # Wfcaa=torch.load(model_filename)
    # print(f"Model saved at {model_filename}")
    # print(Wfcaa)
    return
    
if __name__=="__main__":
    modelParams_filename = sys.argv[1]
    main(modelParams_filename)