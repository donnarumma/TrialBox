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
from scipy.stats import pearsonr
import plotly.express as px
import matplotlib.pyplot as plt

def encoding_sim(S, E, alphas=np.logspace(-3, 3, 10), n_folds=5):
    """
    Simulates the encoding process using ridge regression with cross-validation.

    Args:
        S: np.ndarray, shape (n_words, n_elecs, n_lags)
            Neural activity (per word, per electrode, per lag.
        E: np.ndarray, shape (n_words, n_feat)
            Word embeddings (e.g., GPT-2 vectors).
        alphas: list or np.ndarray
            Regularization parameters to test for RidgeCV.
        n_folds: int
            Number of folds for cross-validation.

    Returns:
        A dictionary with:
            - "preds": predicted neural responses, shape (n_words, n_elecs, n_lags)
            - "true": actual neural responses, same shape
            - "coefs": model weights (averaged across folds), shape (n_feat, n_elecs, n_lags)
    """
    n_words, n_elecs, n_lags = S.shape
    preds_all = np.zeros((n_words, n_elecs, n_lags))
    true_all = np.zeros((n_words, n_elecs, n_lags))
    coefs = np.zeros((E.shape[1], n_elecs, n_lags))

    kf = KFold(n_splits=n_folds)
    for train_idx, test_idx in kf.split(E):
        E_train, E_test = E[train_idx], E[test_idx]
        S_train, S_test = S[train_idx], S[test_idx]

        # Z-score standardization
        scaler_E = StandardScaler()
        scaler_S = StandardScaler()
        E_train = scaler_E.fit_transform(E_train)
        E_test = scaler_E.transform(E_test)

        # Flatten electrodes × lags into single feature dimension for standardization
        S_train = scaler_S.fit_transform(S_train.reshape(len(train_idx), -1)).reshape(S_train.shape)
        S_test = scaler_S.transform(S_test.reshape(len(test_idx), -1)).reshape(S_test.shape)

        # Train a separate model for each electrode and lag
        for e in range(n_elecs):
            for l in range(n_lags):
                y_train = S_train[:, e, l]
                y_test = S_test[:, e, l]

                model = RidgeCV(alphas=alphas, store_cv_values=False)
                model.fit(E_train, y_train)

                y_pred = model.predict(E_test)
                preds_all[test_idx, e, l] = y_pred
                true_all[test_idx, e, l] = y_test
                coefs[:, e, l] += model.coef_  # accumulate across folds

    return {
        "preds": preds_all,
        "true": true_all,
        "coefs": coefs / n_folds,  # average across folds
    }


def get_lag_correlation(preds, true):
    """
    Computes a lag-by-lag correlation matrix between predicted and true neural responses,
    similar to the intersubject encoding lag analysis in Hasson et al.

    Args:
        preds: np.ndarray, shape (n_words, n_elecs, n_lags)
            Predicted neural activity.
        true: np.ndarray, shape (n_words, n_elecs, n_lags)
            Actual observed neural activity.

    Returns:
        corr_matrix: np.ndarray, shape (n_lags, n_lags)
            Matrix of correlation values between each prediction lag and each observation lag.
    """
    n_words, n_elecs, n_lags = preds.shape
    corr_mat = np.zeros((n_lags, n_lags))

    # Average over electrodes (like Hasson et al. do)
    # shape: (n_words, n_lags)
    preds_mean = preds.mean(axis=1)
    # shape: (n_words, n_lags)
    true_mean = true.mean(axis=1)    

    # Compute correlation for every lag pair (i, j)
    for i in range(n_lags):
        for j in range(n_lags):
            corr, _ = pearsonr(preds_mean[:, i], true_mean[:, j])
            corr_mat[i, j] = corr

    return corr_mat


def average_bidirectional(matrix1, matrix2):
    """
    Computes symmetric average of two lag-lag correlation matrices.

    Used to combine:
        - speaker-to-listener encoding (matrix1)
        - listener-to-speaker encoding (matrix2, transposed to match)

    Args:
        matrix1: np.ndarray, shape (n_lags, n_lags)
        matrix2: np.ndarray, shape (n_lags, n_lags)

    Returns:
        np.ndarray: symmetric average, shape (n_lags, n_lags)
    """
    return (matrix1 + matrix2.T) / 2


# Example
Lags = 10
Electrodes = 100
Words = 50
GPT_dim = 768
# synthetic neural data
S = np.random.randn(Words, Electrodes, Lags)     
# synthetic embeddings
E = np.random.randn(Words, GPT_dim)              

results = encoding_sim(S, E)
corr_mat = get_lag_correlation(results["preds"], results["true"])

# Definiamo lags (esempio: da -4s a +4s con 10 step)
n_lags = corr_mat.shape[0]
lags = np.linspace(-4, 4, n_lags)

# Plot
fig, ax = plt.subplots()
im = ax.imshow(corr_mat, origin="lower", cmap="coolwarm", vmin=-1, vmax=1)
plt.show()

