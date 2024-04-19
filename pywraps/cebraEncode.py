import sys
import h5py
import torch
from cebra import CEBRA

def main(modelParams_filename):
    #modelParams
    with h5py.File(modelParams_filename, 'r') as hdf:
        model_filename      = hdf.attrs['model_filename']
        group_field         = hdf.attrs['group_field']
        neural_field        = hdf.attrs['neural_field']
        neural_filename     = hdf.attrs['neural_filename']
        manifold_field      = hdf.attrs['manifold_field']
        manifold_filename   = hdf.attrs['manifold_filename']

    # load neural 
    with h5py.File(neural_filename, 'r') as hdf:
        neural_data    = hdf[group_field + '/' +  neural_field][:]

    # load model
    device = "cuda" if torch.cuda.is_available() else "cpu"
    cebra_model  = torch.load(model_filename,map_location=torch.device(device))
    cebra_model.to(device)
    
    # encode new data with model fitted
    manifold     = cebra_model.transform(neural_data)
    
    # save data in 
    with h5py.File(manifold_filename, 'w') as out_hdf:
        group = out_hdf.create_group(group_field)
        group.create_dataset(manifold_field, data=manifold, compression='gzip', compression_opts=9)
        print(f"{manifold_field} saved in {manifold_filename}")   
    return
    
if __name__=="__main__":
    modelParams_filename = sys.argv[1]
    main(modelParams_filename)