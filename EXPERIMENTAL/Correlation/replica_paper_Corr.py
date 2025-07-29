# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 21:04:04 2025

@author: loren
"""

# encoding_simulation.py


# Hasson_simulation.py
"""
Simplified and documented replication of the encoding model and
inter-subject correlation analysis (as in Hasson Zada et al.2024),
adapted for synthetic data.
"""

import numpy as np
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from scipy.spatial.distance import cdist

# Simulate data
np.random.seed(0)
n_words, n_feat, n_elecs, n_lags = 50, 100, 32, 10
 # Shared embeddings
E = np.random.randn(n_words, n_feat) 
# Speaker neural data
S_speaker = np.random.randn(n_words, n_elecs, n_lags)  
 # Listener neural data
S_listener = np.random.randn(n_words, n_elecs, n_lags) 

def encoding_sim(S, E, alphas=np.logspace(-3, 3, 5), n_folds=5):
    n_words, n_elecs, n_lags = S.shape
    preds_all = np.zeros((n_words, n_elecs, n_lags))
    true_all = np.zeros((n_words, n_elecs, n_lags))
    coefs = np.zeros((E.shape[1], n_elecs, n_lags))

    kf = KFold(n_splits=n_folds)
    for train_idx, test_idx in kf.split(E):
        E_train, E_test = E[train_idx], E[test_idx]
        S_train, S_test = S[train_idx], S[test_idx]

        scaler_E = StandardScaler()
        scaler_S = StandardScaler()
        E_train = scaler_E.fit_transform(E_train)
        E_test = scaler_E.transform(E_test)
        S_train = scaler_S.fit_transform(S_train.reshape(len(train_idx), -1)).reshape(S_train.shape)
        S_test = scaler_S.transform(S_test.reshape(len(test_idx), -1)).reshape(S_test.shape)

        for e in range(n_elecs):
            for l in range(n_lags):
                y_train = S_train[:, e, l]
                y_test = S_test[:, e, l]

                model = RidgeCV(alphas=alphas)
                model.fit(E_train, y_train)
                y_pred = model.predict(E_test)

                preds_all[test_idx, e, l] = y_pred
                true_all[test_idx, e, l] = y_test
                coefs[:, e, l] += model.coef_

    return {
        "preds": preds_all,
        "true": true_all,
        "coefs": coefs / n_folds,
    }

def get_lag_lag_matrix(a, b):
    # average over electrodes
    a = a.mean(1)  
    b = b.mean(1)
    # shape: (n_lags, n_lags)
    return 1 - cdist(a.T, b.T, metric='correlation')  

# Step 1: train on speaker, predict on speaker embeddings
res_speaker = encoding_sim(S_speaker, E)
res_listener = encoding_sim(S_listener, E)

# Step 2: compute S2L (speaker preds vs listener true)
s2l_matrix = get_lag_lag_matrix(res_speaker["preds"], res_listener["true"])

# Step 3: compute L2S (listener preds vs speaker true)
l2s_matrix = get_lag_lag_matrix(res_listener["preds"], res_speaker["true"])

# Step 4: average
ise_matrix = (s2l_matrix + l2s_matrix.T) / 2

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(5, 5))
im = ax.imshow(ise_matrix, origin="lower", cmap="coolwarm", vmin=-1, vmax=1)
ax.set_title("Avg Corr (S2L + L2S average)")
ax.set_xlabel("Listener lag")
ax.set_ylabel("Speaker lag")
fig.colorbar(im, ax=ax)
plt.tight_layout()
plt.show()
