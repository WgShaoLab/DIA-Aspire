import warnings
warnings.filterwarnings("ignore")

from pathlib import Path
import pandas as pd
import numpy as np
import os
from .rescore_ML import Percolator_ML
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def PerformRescore(psm_df: pd.DataFrame, output_dir: str, model: str):
    
    psm_all = psm_df.copy()
    psm_all['ml_score'] = 1-psm_all['fdr']

    col1 = ['raw_name','Precursor.Id','Stripped.Sequence','precursor_intensity','Protein.Group']
    col2 = ['decoy','ml_score','fdr']

    ms2_rt_features = [
    'ms1_profile_corr', 'evidence', 'pep',
    'cos', 'sa', 'spc','pcc', 
    'merr_weighted_frag_score', 'pred_weighted_frag_score',
    'matched_frag_num', 'matched_frag_ratio', 
    'both_matched_pred_frag_to_matched', 
    'matched_not_pred_frag_ratio',
    'matched_frag_rel_to_pred', # good
    'pred_frag_rel_to_matched', # good
    'merr_weighted_bion_score',
    'merr_weighted_yion_score', 
    'pred_yion_rel_to_matched', 
    'abs_rt_delta', 'rt_ratio', 
        ]
    xic_features = [
        "cos_corr_top3_median",
        "cos_corr_top3_mean",
        "pearson_corr_top3_median",
        "pearson_corr_top3_mean",
        "cos_corr_filtered_top3_median",
        "cos_corr_filtered_top3_mean",
        "pearson_corr_filtered_top3_median",
        "pearson_corr_filtered_top3_mean",
        "fragment_filter_ratio",
        "peak_symmetry",
        "peak_jaggedness",
        "peak_modality",
        "apex_rt_ratio",
        "peak_tailing_factor",
        "co_elution_score",
        "main_peak_consistency",
        ]
    
    features = ms2_rt_features + xic_features
    # features = ms2_generator.feature_names + rt_generator.feature_names + ['ms1_profile_corr', 'evidence', 'pep']
    col3 = col1 + features + col2
    psm_all1 = psm_all[col3]
    data = psm_all1.copy()
 
    def convert_matrix(data):
        col = ['Precursor.Id','raw_name','precursor_intensity']
        data1 = data[col]
        res = data1.pivot_table(index=['Precursor.Id'],columns=['raw_name'],values=['precursor_intensity']).reset_index()
        res1 = res.set_index(['Precursor.Id'])
        res1.columns = res1.columns.droplevel(0)
        res1 = res1.reset_index().rename_axis([None], axis=1)
        return res1

    def rescoring(data,backend,model):
        perc_config= {
                    "percolator_backend": backend,
                    "percolator_model": model,
                    "min_train_sample": 1000,
                    "max_train_sample": 30000,
                    "cv_fold": 3,
                    "iter_num":8,
                    "fdr": 0.01,
                    "fdr_level": "precursor",
                    "per_raw_fdr": False,
                    "feature_list": features,
                    "outpath":output_dir,
                }
        perc = Percolator_ML(**perc_config)
        dr = perc.re_score(data)
        dr1 = dr[(dr['fdr']<=0.01) & (dr['decoy']==0)]
        dr2 = dr1[dr1['precursor_intensity']>0]

        pre_target_num = dr1['Precursor.Id'].nunique()

        dr1d = dr[(dr['fdr']<=0.01) & (dr['decoy']==1)]
        pre_decoy_num = dr1d['Precursor.Id'].nunique()

        outputfile = os.path.join(output_dir, f"{model}_target_precursor.csv")
        dr2.to_csv(outputfile,index=False)

        rescore_df = convert_matrix(dr1)
        rescore_df2 = dr2[['Stripped.Sequence','Protein.Group','Precursor.Id']].drop_duplicates()
        dr_matrix = rescore_df.merge(rescore_df2,on='Precursor.Id',how='inner')

        outputfile2 = os.path.join(output_dir, f"rescore_pr_matrix.csv")
        dr_matrix.to_csv(outputfile2,index=False)

        logger.info(f"model: {model}")
        logger.info(f"target_precursor: {pre_target_num}")
        logger.info(f"decoy_precursor: {pre_decoy_num}")
        logger.info(f"result_precursor: {pre_target_num - pre_decoy_num}")
        logger.info(f"final result saved to: {outputfile}")
        return True


    # for backend in ['sklearn','pytorch']:
    #     # if backend == 'sklearn':
    #     #     for model in ['SVM']:  # ,'RF',,'LR'
    #     #         rescoring(data,backend,model)

    #     if backend == 'pytorch':
    #         for model in ['DNN']:   # 'NNlinear','CNN'
    #             rescoring(data,backend,model)

    modelx = model
    if modelx == 'DNN':
        backend = 'pytorch'
    elif modelx == 'SVM':
        backend = 'sklearn'
    rescoring(data,backend,modelx)
