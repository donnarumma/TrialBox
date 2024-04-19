import h5py
import cebra.datasets
from cebra import CEBRA

model_params = {
    "model_architecture": 'offset10-model',
    "batch_size": "512",
    "learning_rate": "3e-4",
    "temperature": "1",
    "output_dimension": "3",
    "max_iterations": "10000",
    "distance": "cosine",
    "conditional": "time_delta",
    "hybrid": "False",
    "time_offsets": "10",
    "verbose": "true",
    "seed": "42",
    "maxI": "9999",
    "model_filename": "cebra_model.pkl",
    "group_field": "data",       
    "neural_filename": "neural.hd5",
    "neural_field": "neural",
    "behavior_filename": "behavior.hd5",
    "behavior_field": "behavior"    
    }

encode_params={
    "model_filename": "cebra_model.pkl",
    "group_field": "data",
    "neural_filename": "neural.hd5",
    "neural_field": "neural",
    "manifold_filename": "manifold.hd5",
    "manifold_field": "manifold"
}

plot_params={
    "group_field": "data",
    "behavior_filename": "behavior.hd5",
    "behavior_field": "behavior",
    "manifold_filename": "manifold.hd5",
    "manifold_field": "manifold"
}

modelParams_filename    = "cebraModelParams.hd5"
encodeParams_filename   = "cebraEncodeParams.hd5"
plotCebraParams_filename= "plotCebraParams.hd5"

def main():
    # get dataset
    hippocampus_pos     = cebra.datasets.init('rat-hippocampus-single-achilles')
    neural_data         = hippocampus_pos.neural
    behavior_data       = hippocampus_pos.continuous_index.numpy()
    
    group_field         = model_params.get("group_field");
    behavior_filename   = model_params.get("behavior_filename");
    behavior_field      = model_params.get("behavior_field");
    neural_filename     = model_params.get("neural_filename");
    neural_field        = model_params.get("neural_field");

    # save neural and behavior in hd5
    with h5py.File(behavior_filename, 'w') as out_hdf:
        group = out_hdf.create_group(group_field)
        group.create_dataset(behavior_field, data=behavior_data,compression='gzip', compression_opts=9)
        print(f"{behavior_field} saved in {behavior_filename}")
    with h5py.File(neural_filename, 'w') as out_hdf:
        group = out_hdf.create_group(group_field)
        group.create_dataset(neural_field, data=neural_data,compression='gzip', compression_opts=9)
        print(f"{neural_field} saved in {neural_filename}")
    # save model params
    with h5py.File(modelParams_filename, 'w') as hdf:
        for param,   value in model_params.items():    
            hdf.attrs[param] = value
        print(f"modelParams saved in {modelParams_filename}")
    # save encode params
    with h5py.File(encodeParams_filename, 'w') as hdf:
        for param,   value in encode_params.items():    
            hdf.attrs[param] = value
        print(f"encodeParams saved in {encodeParams_filename}")
    # save plot params
    with h5py.File(plotCebraParams_filename, 'w') as hdf:
        for param,   value in plot_params.items():    
            hdf.attrs[param] = value
        print(f"plotParams saved in {plotCebraParams_filename}")
    return
    
if __name__=="__main__":
    main()
