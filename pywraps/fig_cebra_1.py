# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 17:47:10 2024

@author: zlollo2
"""
#sys.path.insert(0,'/path/to/your/directory')
import os
os.getcwd()
import sys
#### Mettere la directory di interesse la stessa di matlab
from pathlib import Path
#base_dir=r'/home/zlollo/CNR/Cebra_for_all'
#os.chdir(base_dir)
import time

#!pip install ripser
#import ripser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib as jl
import cebra.datasets
from cebra import CEBRA
from scipy.io import loadmat
from scipy.io import savemat
#from dataset import SingleRatDataset  # Assumendo che il codice sia in 'dataset.py'
from matplotlib.collections import LineCollection
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
import sklearn.metrics
import inspect
import torch
from cebra.datasets.hippocampus import *
#import tensorflow as tf

def plot_cebra(emb, label):
    fig = plt.figure(figsize=(10, 5))
        
        
    # Variabili per gestire la colorbar unica
    cmap1 = plt.get_cmap('cool')
    cmap2 = plt.get_cmap('summer')
    #label = method_viz['hypothesis']["label"]
    norm = plt.Normalize(vmin=label[:, 0].min(), vmax=label[:, 0].max())
    
     #    norm = plt.Normalize(vmin=label[:, 0].min(), vmax=label[:, 0].max())
        
         # Ciclo per ogni modello
    # for i, model in enumerate(['Hypothesis: position', 
    #    'Discovery: time only','Hybrid: time + behavior']):
    #    # Preparazione del grafico
       # fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection="3d")
     
    # emb = method_viz[model]["embedding"]
    #    label = method_viz[model]["label"]
    idx1, idx2, idx3 = (0, 1, 2)
    #    if i == 3:
    #        idx1, idx2, idx3 = (1, 2, 0)
     
    r = label[:, 1] == 1
    l = label[:, 2] == 1
     
       # Plot per sinistra e destra
    ax.scatter(emb[l, idx1], emb[l, idx2], emb[l, idx3], c=label[l, 0], 
                 cmap=cmap1, norm=norm, s=0.1, label='Left')
    ax.scatter(emb[r, idx1], emb[r, idx2], emb[r, idx3], c=label[r, 0], 
                  cmap=cmap2, norm=norm, s=0.1, label='Right')
     
       # Axes Removal
    #ax.axis("off")
     
       # Titolo
    ax.set_title("{Cebra}", fontsize=20)
     
       # Aggiunta delle annotazioni di testo per indicare le direzioni
       #ax.text2D(0.05, 0.95, "Left", transform=ax.transAxes)
       #ax.text2D(0.95, 0.95, "Right", transform=ax.transAxes)
     
    # Aggiungi una colorbar per ciascun grafico
    sm = plt.cm.ScalarMappable(cmap=cmap1, norm=norm)
    sm.set_array([])
       # # Colorbar 1
    cbar = plt.colorbar(sm, ax=fig.axes, orientation='vertical', fraction=0.01, pad=0.15)
    cbar.set_label('left')
     
    sm2 = plt.cm.ScalarMappable(cmap=cmap2, norm=norm)
    sm2.set_array([])
    cbar2 = plt.colorbar(sm2, ax=fig.axes, orientation='vertical', fraction=0.01, pad=0.15)  # Imposta pad a un valore diverso
    cbar2.set_label('right')
    
      # save_fig(fig,f"plot_{model}", tight_layout=False, fig_extension="eps")        
       # Mostra il grafico
    plt.show()