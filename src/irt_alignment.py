import warnings
# from .util import timestamped_echo

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

import click
import os
import sys
import pathlib
import posixpath, ntpath
import re
import operator
import numpy as np
import pandas as pd

# alignment
from sklearn import preprocessing
import sklearn.isotonic
import sklearn.linear_model
import statsmodels.api as sm
from scipy.interpolate import interp1d

try:
    
    from pyprophet.stats import pemp, qvalue, pi0est
    from pyprophet.ipf import compute_model_fdr
except ModuleNotFoundError:
    pass

from scipy.stats import gaussian_kde
from numpy import linspace, concatenate
from seaborn import lmplot

from datetime import datetime
import click

def timestamped_echo(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    click.echo(f"{timestamp} - {message}")

    
def lowess_iso(x, y, lowess_frac):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='invalid value encountered in ', category=RuntimeWarning)
        lwf = sm.nonparametric.lowess(y, x.ravel(), frac=lowess_frac)
    while pd.isna(lwf[:, 1]).any():
        lowess_frac *= 2
        lwf = sm.nonparametric.lowess(y, x.ravel(), frac=lowess_frac)
    lwf_x = lwf[:, 0]
    ir = sklearn.isotonic.IsotonicRegression()  # make the regression strictly increasing
    lwf_y = ir.fit_transform(lwf_x, lwf[:, 1])
    mask = np.concatenate([[True], np.diff(lwf_y) != 0])  # remove non increasing points
    try:
        return interp1d(lwf_x[mask], lwf_y[mask], bounds_error=False, fill_value="extrapolate")
    except ValueError as e:
        timestamped_echo(e)
    return interp1d(lwf_x, lwf_y, bounds_error=False, fill_value="extrapolate")

class LowessIsoEstimator:
    def __init__(self, lowess_frac):
        self.lowess_frac = lowess_frac

    def fit(self, x, y):
        self.lwi = lowess_iso(x, y, self.lowess_frac)
        return self

    def get_params(self, deep=False):
        return {'lowess_frac': self.lowess_frac}

    def set_params(self, lowess_frac):
        self.lowess_frac = lowess_frac
        return self

    def score(self, x, y):
        resid = self.lwi(x.ravel()) - y
        return 1 / resid.dot(resid)

    def predict(self, x):
        return self.lwi(x.ravel())

    def __repr__(self):
        return str(self.get_params())

def lowess_iso_predictor(filename, x, y, xpred):
    gsc = sklearn.model_selection.GridSearchCV(LowessIsoEstimator(None), {'lowess_frac': [0.01, 0.02, 0.04, 0.08]},
                                             cv=sklearn.model_selection.KFold(4, shuffle=True, random_state=0),
                                             n_jobs=min(os.cpu_count(), 61))

    gsc.fit(x.reshape(-1, 1), y)
    timestamped_echo(f'Info: {filename}; Lowess fraction used: {gsc.best_params_["lowess_frac"]}.')
    return gsc.best_estimator_.predict(xpred)
      
def lowess(run, reference_run, xcol, ycol, lowess_frac, psm_fdr_threshold, min_peptides, filename, main_path):
  # Filter alignment data
    run_alignment = run[run['q_value'] < psm_fdr_threshold] if 'q_value' in run else run
    
    if 'q_value' in reference_run:
        reference_run_alignment = reference_run[reference_run['q_value'] < psm_fdr_threshold]
    else:
        reference_run_alignment = reference_run

    dfm = pd.merge(run_alignment, reference_run_alignment[['modified_peptide','precursor_charge',ycol]], on=['modified_peptide','precursor_charge'])
    timestamped_echo(f'Info: {filename}; Peptide overlap between run and reference: {dfm.shape[0]}.')
    if dfm.shape[0] <= min_peptides:
        timestamped_echo(f'Info: {filename}; Skipping run because not enough peptides could be found for alignment.')
        return pd.DataFrame()

    if dfm.shape[0] < 50:  # use linear regression for small reference size
        linreg = sklearn.linear_model.LinearRegression().fit(dfm[xcol].to_numpy().reshape(-1, 1), dfm[ycol])
        run[ycol] = linreg.predict(run[xcol].to_numpy().reshape(-1, 1))
    else:
    # Fit and apply the lowess model
        run[ycol] = lowess_iso_predictor(filename, dfm[xcol].to_numpy(), dfm[ycol].to_numpy(), run[xcol].to_numpy()) \
        if lowess_frac == 0 else \
        lowess_iso(dfm[xcol].to_numpy(), dfm[ycol].to_numpy(), lowess_frac)(run[xcol].to_numpy())

    return run

def submod(text):
    t1 = re.sub('M\[147\]', 'M(UniMod:35)', text)
    t2 = re.sub('S\[167\]', 'S(UniMod:21)', t1)
    t3 = re.sub('T\[181\]', 'T(UniMod:21)', t2)
    t4 = re.sub('Y\[243\]', 'Y(UniMod:21)', t3)
    t5 = re.sub('N\[115\]', 'N(UniMod:7)', t4)
    t6 = re.sub('Q\[129\]', 'Q(UniMod:7)', t5)
    t7 = re.sub('S\[129\]', 'Q(UniMod:1)', t6)
    t8 = re.sub('n\[43\]', '(UniMod:1)', t7)
    return t8 

def lowess2(run, reference_run, xcol, ycol, lowess_frac, psm_fdr_threshold, min_peptides):
  # Filter alignment data
    run_alignment = run
    
    reference_run_alignment = reference_run

    dfm = pd.merge(run_alignment, reference_run_alignment, on=['modified_peptide','precursor_charge'])
    timestamped_echo(f'Info: Peptide overlap between SysteMHC and reference: {dfm.shape[0]}.')

    if dfm.shape[0] < 50:  # use linear regression for small reference size
        linreg = sklearn.linear_model.LinearRegression().fit(dfm[xcol].to_numpy().reshape(-1, 1), dfm[ycol])
        run[ycol] = linreg.predict(run[xcol].to_numpy().reshape(-1, 1))
    else:
    # Fit and apply the lowess model
        run[ycol] = lowess_iso_predictor(filename, dfm[xcol].to_numpy(), dfm[ycol].to_numpy(), run[xcol].to_numpy()) \
        if lowess_frac == 0 else \
        lowess_iso(dfm[xcol].to_numpy(), dfm[ycol].to_numpy(), lowess_frac)(run[xcol].to_numpy())
    return run