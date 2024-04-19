import sys
import matplotlib
import h5py
matplotlib.use('TkAgg')
from fig_cebra_1 import plot_cebra

def plot_results (plotParams_filename):
    with h5py.File(plotParams_filename, 'r') as hdf:
        group_field         = hdf.attrs['group_field']
        behavior_field      = hdf.attrs['behavior_field']
        behavior_filename   = hdf.attrs['behavior_filename']
        manifold_field      = hdf.attrs['manifold_field']
        manifold_filename   = hdf.attrs['manifold_filename']

    with h5py.File(behavior_filename, 'r') as hdf:
        behavior    = hdf[group_field + '/' + behavior_field][:]
    with h5py.File(manifold_filename, 'r') as hdf:
        manifold    = hdf[group_field + '/' + manifold_field][:]
    
    plot_cebra(manifold, behavior)

if __name__ == "__main__":
    plotParams_filename = sys.argv[1] 
    plot_results(plotParams_filename)