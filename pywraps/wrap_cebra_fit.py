# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from cebra import CEBRA
from scipy.io import loadmat
#from cebra.datasets.hippocampus import *
import scipy.io
import joblib

from scipy.io import loadmat
import torch

# Check # input arguments
if len(sys.argv) != 3:
    print("Usage: python script.py <input_directory> <output_directory>")
    sys.exit(1)

input_directory     = sys.argv[1]
output_directory    = sys.argv[2]
################## load neural and behavioral data #################################
data_mat            = loadmat(os.path.join(input_directory, 'rat_n.mat'))
data_               = data_mat['rat_n']
label_mat           = loadmat(os.path.join(input_directory, 'rat_b.mat'))
label_              = label_mat['rat_b']

##################### load model parameters #############################
m_params            = loadmat(os.path.join(input_directory,'params.mat'))
model_params        = m_params['params']
    
try:
    mod_type        = model_params['model_type'][0].item() 
    mod_type        = mod_type[0]
    
    mod_arch        = model_params['mod_arch'][0].item() 
    mod_arch        = mod_arch[0]
    
    dist            = model_params['distance'][0].item() 
    dist            = dist[0]
    
    cond            = model_params['conditional'][0].item() 
    cond            = cond[0]
    
    temp            = int(model_params['temperature'][0][0]);

    time_off        = int(model_params['time_offsets'][0][0])

    max_iter        = int(model_params['max_iter'][0][0])

    max_adapt_iter  = int(model_params['max_adapt_iter'][0][0])

    b_size          = int(model_params['batch_size'][0][0])

    learn_rate      = int(model_params['learning_rate'][0][0])
    
    out_dim         = int(model_params['output_dimension'][0][0])
    
    verb            = model_params['verbose'][0].item() 
    verb            = verb[0]
    
    n_h_u           =int(model_params['num_hidden_units'][0][0])
    
    pad_before_transform_= model_params['pad_before_transform'][0].item() 
    # Convert string of 'pad_before_transform' in logical value
    if pad_before_transform_ == "True" or pad_before_transform_ == "true":
        p_b_t = True
    else:
        p_b_t = False
        
    hybrid_= model_params['hybrid'][0].item() 
    if hybrid_ == "True" or hybrid_== "true":
        hyb = True
    else:
        hyb = False

except FileNotFoundError:
   print("model_params file is missing!!!")

#############################################################
# Verify if GPU is available
if torch.cuda.is_available():
    # Print info on GPU
    print(f"GPU available: {torch.cuda.get_device_name(0)}")
else:
    print("no GPU available switching to CPU.")


################################## CNN MODEL FIT ###############################
torch.cuda.empty_cache()
cebra_target_model = CEBRA(model_architecture=mod_arch,
                           distance=dist,
                           conditional= cond,
                           temperature=temp,
                           time_offsets=time_off,
                           max_iterations=max_iter,
                           max_adapt_iterations=max_adapt_iter,
                           batch_size=b_size,
                           learning_rate=learn_rate,
                           output_dimension=out_dim,
                           verbose=verb,
                           num_hidden_units=n_h_u,
                           pad_before_transform=p_b_t,
                           hybrid=hyb,
                           device='cuda_if_available')


def run_model(model, data, labels, model_type):
        ## CEBRA BEHAVIOUR
        if model_type == "hypothesis":
            # If the model is in supervised mode, use both data and labels
            model.fit(data, labels)
        ## CEBRA TIME      
        elif model_type == "discovery":
             model.fit(data)
        elif model_type == "shuffle":
            shuffled_labels = np.random.permutation(labels)
            model.fit(data, shuffled_labels)
         # save model in the output directory 
        joblib.dump(model, os.path.join(output_directory,'fitted_model.pkl'))
        with torch.no_grad():
            return model.transform(data), model.model_.state_dict()
    
cebra_output, ceb_model = run_model(cebra_target_model, data_, label_, mod_type)

### generate and save output (to be deleted)
cebra_mat={'cebra_output':cebra_output}

numpy_model_    = {key.replace('.', '_'): value for key, value in ceb_model.items()}
numpy_model     = {key: value.detach().cpu().numpy() for key, value in numpy_model_.items()}

# save fitted model and output
scipy.io.savemat(os.path.join(output_directory, 'cebra_output.mat'), cebra_mat)
scipy.io.savemat(os.path.join(output_directory, 'model_struct.mat'), numpy_model)

