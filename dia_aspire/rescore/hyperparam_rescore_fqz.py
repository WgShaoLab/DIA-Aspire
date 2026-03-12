import numpy as np
import pandas as pd
import torch
import os
import sys
import multiprocessing as mp
from tqdm import tqdm
import contextlib
import shap
import matplotlib.pyplot as plt

from .fdr import fdr_from_ref, fdr_to_q_values
from peptdeep.utils import logging
from sklearn.model_selection import GridSearchCV, KFold

import warnings
warnings.filterwarnings("ignore", message="The figure layout has changed to tight")
warnings.filterwarnings("ignore", message="'mode' parameter is deprecated")

NUM_THREADS = 16  # 可以根据服务器情况调整，如4、8、16等

# 限制所有科学计算库的线程数
os.environ['OMP_NUM_THREADS'] = str(NUM_THREADS)
os.environ['MKL_NUM_THREADS'] = str(NUM_THREADS)
os.environ['OPENBLAS_NUM_THREADS'] = str(NUM_THREADS)
os.environ['VECLIB_MAXIMUM_THREADS'] = str(NUM_THREADS)
os.environ['NUMEXPR_NUM_THREADS'] = str(NUM_THREADS)

# 限制PyTorch线程数
torch.set_num_threads(NUM_THREADS)
torch.set_num_interop_threads(NUM_THREADS)

# CUDA内存优化
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True,max_split_size_mb:128'


# ==================== Hyperparameter Grid Definitions ====================
def get_hyperparam_grid(model_type: str) -> dict:
    """Define hyperparameter search space for different models"""
    if model_type == "svm":
        return {
            'C': [0.1, 1, 10],
            'epsilon': [0.01, 0.1, 0.5]
        }
    
    elif model_type == "lr":
        return {
            'C': [0.1, 1, 10],
            'penalty': ['l2'],
            'solver': ['lbfgs', 'liblinear']
        }
    
    elif model_type == "rf":
        return {
            'n_estimators': [100, 200],
            'max_depth': [None, 10, 20],
            'min_samples_split': [10, 20],
            'min_samples_leaf': [5, 10],
            'max_features': ['sqrt', 'log2']
        }
    
    elif model_type == "lda":
        return {
            'solver': ['lsqr', 'eigen'],
            'shrinkage': ['auto', 0.5]
        }
    
    elif model_type == "bayes":
        return {
            'var_smoothing': [1e-9, 1e-8, 1e-7, 1e-6]
        }
    
    elif model_type == "dnn":
        return {
            'hidden_dim1': [32, 64, 128],
            'hidden_dim2': [8, 16, 32],
            'dropout': [0.1, 0.2, 0.3]
        }
    
    elif model_type == "cnn":
        return {
            'hidden_channels': [64, 128, 256],
            'kernel_size': [3, 5],
            'dropout': [0.1, 0.2, 0.3]
        }
    
    else:
        return {}


# ==================== CNN/DNN Models ====================
class CNNModel(torch.nn.Module):
    """CNN model"""
    def __init__(self, input_dim, hidden_channels=64, kernel_size=3, dropout=0.2):
        super().__init__()
        torch.manual_seed(1290)
        
        self.conv1 = torch.nn.Conv1d(1, hidden_channels, kernel_size, padding=kernel_size//2)
        self.bn1 = torch.nn.BatchNorm1d(hidden_channels)
        
        self.conv2 = torch.nn.Conv1d(hidden_channels, hidden_channels*2, kernel_size, padding=kernel_size//2)
        self.bn2 = torch.nn.BatchNorm1d(hidden_channels*2)
    
        self.conv3 = torch.nn.Conv1d(hidden_channels*2, hidden_channels*4, kernel_size, padding=kernel_size//2)
        self.bn3 = torch.nn.BatchNorm1d(hidden_channels*4)
        
        self.global_avg_pool = torch.nn.AdaptiveAvgPool1d(1)
        self.fc = torch.nn.Linear(hidden_channels*4, 1)
        self.dropout = torch.nn.Dropout(dropout)
        self.relu = torch.nn.ReLU()
    
    def forward(self, x):
        if x.dim() == 2:
            x = x.unsqueeze(1)
        elif x.dim() == 3:
            if x.shape[1] != 1:
                if x.shape[2] == 1:
                    x = x.squeeze(-1).unsqueeze(1)
                else:
                    raise ValueError(f"Expected channel dim=1, got {x.shape}")
        elif x.dim() == 4:
            x = x.squeeze()
            if x.dim() == 2:
                x = x.unsqueeze(1)
            elif x.dim() == 3:
                if x.shape[1] == 1:
                    pass
                elif x.shape[2] == 1:
                    x = x.squeeze(-1).unsqueeze(1)
                else:
                    raise ValueError(f"Cannot handle 3D input with shape {x.shape}")
            else:
                raise ValueError(f"Unexpected shape after squeeze: {x.shape}")
        else:
            raise ValueError(f"Unsupported input dimension: {x.dim()}")

        assert x.shape[1] == 1, f"Channel dim must be 1, got {x.shape}"

        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu(x)
        x = self.dropout(x)
        
        x = self.conv2(x)
        x = self.bn2(x)
        x = self.relu(x)
        x = self.dropout(x)
        
        x = self.conv3(x)
        x = self.bn3(x)
        x = self.relu(x)
        x = self.dropout(x)
        
        x = self.global_avg_pool(x)
        x = x.squeeze(-1)
        x = self.fc(x)
        return x


class DNNModel(torch.nn.Module):
    """DNN model"""
    def __init__(self, input_dim, hidden_dim1=None, hidden_dim2=None, dropout=0.2):
        super().__init__()
        torch.manual_seed(1290)
        
        hidden_dim1 = hidden_dim1 or max(32, input_dim // 2)
        hidden_dim2 = hidden_dim2 or max(16, input_dim // 4)
        
        self.model = torch.nn.Sequential(
            torch.nn.Linear(input_dim, hidden_dim1),
            torch.nn.BatchNorm1d(hidden_dim1),
            torch.nn.ReLU(),
            torch.nn.Dropout(dropout),
            
            torch.nn.Linear(hidden_dim1, hidden_dim2),
            torch.nn.BatchNorm1d(hidden_dim2),
            torch.nn.ReLU(),
            torch.nn.Dropout(dropout),
            
            torch.nn.Linear(hidden_dim2, 1)
        )
    
    def forward(self, x):
        return self.model(x)


class NNLinear(torch.nn.Module):
    """Torch-based NN_Linear rescore model"""
    def __init__(self, input_dim, **kwargs):
        super().__init__()
        torch.manual_seed(1290)
        self.linear = torch.nn.Linear(input_dim, 1)

    def forward(self, x):
        return self.linear(x)


# ==================== Model Provider ====================
class RescoreModelProvider:
    def __init__(self):
        self.model_dict = {}
        self.model_dict["nnlinear"] = NNLinear
        self.model_dict["dnn"] = DNNModel
        self.model_dict["cnn"] = CNNModel

    def register(self, model_name, model_class):
        self.model_dict[model_name.lower()] = model_class

    def get_model(self, model_name, input_dim, **kwargs):
        if model_name.lower() not in self.model_dict:
            print(
                "[Rescore] "
                f"PyTorch rescoring model '{model_name}' is not "
                "implemented, switch to 'linear' model."
            )
            return self.model_dict["nnlinear"](input_dim, **kwargs)
        else:
            return self.model_dict[model_name.lower()](input_dim, **kwargs)
            
rescore_model_provider = RescoreModelProvider()


# ==================== NNLinear Rescorer ====================
class NNLinearRescore:
    def __init__(self, num_features, nn_model_type="nnlinear", **model_kwargs):
        self.nn_model = rescore_model_provider.get_model(nn_model_type, num_features, **model_kwargs)
        self.train_batch_size = 10000
        self.predict_batch_size = 100000
        self.nn_model_type = nn_model_type.lower()
        self.num_features = num_features

        self.optimizer = torch.optim.Adam(
            self.nn_model.parameters(), lr=0.001, weight_decay=1e-4
        )
        self.loss_func = torch.nn.BCEWithLogitsLoss()

        if torch.cuda.is_available():
            self.device = torch.device("cuda")
            self.nn_model.to(self.device)
        else:
            self.device = torch.device("cpu")
        self.epoch = 20

    def fit(self, features, labels):
        labels = torch.tensor(labels, dtype=torch.float, device=self.device)
        sample_idxes = np.random.RandomState(1290).permutation(len(features))

        for _ in range(self.epoch):
            for i in range(0, len(features), self.train_batch_size):
                self.optimizer.zero_grad()

                outputs = self.nn_model(
                    torch.tensor(
                        features[sample_idxes[i : i + self.train_batch_size]],
                        dtype=torch.float,
                        device=self.device,
                    )
                )
                loss = self.loss_func(
                    outputs, labels[sample_idxes[i : i + self.train_batch_size]].unsqueeze(-1)
                )
                loss.backward()
                self.optimizer.step()

    def decision_function(self, features):
        outputs = np.empty(len(features))
        for i in range(0, len(features), self.predict_batch_size):
            outputs[i : i + self.predict_batch_size] = (
                torch.sigmoid(self.nn_model(
                    torch.tensor(
                        features[i : i + self.predict_batch_size],
                        dtype=torch.float,
                        device=self.device,))).detach().squeeze(1).cpu().numpy()
            )
        return outputs


# ==================== NN Rescorer with Hyperparameter Search ====================
class NNRescore:
    def __init__(self, num_features, nn_model_type="dnn", **model_kwargs):
        self.num_features = num_features
        self.nn_model_type = nn_model_type.lower()
        self.model_kwargs = model_kwargs
        self.train_batch_size = 10000
        self.predict_batch_size = 10000
        
        # Initialize model (may be overridden by hyperparameter search)
        self.nn_model = rescore_model_provider.get_model(
            nn_model_type, num_features, **model_kwargs
        )
        
        # Device configuration
        if torch.cuda.is_available():
            self.device = torch.device("cuda")
            self.nn_model.to(self.device)
        else:
            self.device = torch.device("cpu")

        self.optimizer = torch.optim.Adam(
            self.nn_model.parameters(), 
            lr=0.001,
            weight_decay=1e-4
        )
        self.loss_func = torch.nn.BCEWithLogitsLoss()
        self.epoch = 20
        
        if self.nn_model_type in ["dnn", "cnn"]:
            self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                self.optimizer, 
                mode='min',
                factor=0.5,
                patience=2,
            )
            self.patience = 3
        else:
            self.scheduler = None

    def grid_search_hyperparams(self, features, labels):
        """Grid search for PyTorch models"""
        param_grid = get_hyperparam_grid(self.nn_model_type)
        
        if not param_grid:
            logging.info(f"[HyperOpt] No hyperparameter grid for {self.nn_model_type}, using defaults")
            return self.model_kwargs
        
        logging.info(f"[HyperOpt] Starting grid search for {self.nn_model_type}...")
        logging.info(f"[HyperOpt] Search space: {param_grid}")
        
        # Generate all parameter combinations
        from itertools import product
        keys = list(param_grid.keys())
        values = list(param_grid.values())
        combinations = [dict(zip(keys, v)) for v in product(*values)]
        
        logging.info(f"[HyperOpt] Total combinations: {len(combinations)}")
        
        # 3-fold CV
        kf = KFold(n_splits=3, shuffle=True, random_state=42)
        
        best_params = None
        best_score = float('inf')
        
        for i, params in enumerate(combinations):
            logging.debug(f"[HyperOpt] Testing combination {i+1}/{len(combinations)}: {params}")
            
            fold_scores = []
            for fold_idx, (train_idx, val_idx) in enumerate(kf.split(features)):
                X_train_fold = features[train_idx]
                y_train_fold = labels[train_idx]
                X_val_fold = features[val_idx]
                y_val_fold = labels[val_idx]
                
                # Create temporary model
                temp_model = rescore_model_provider.get_model(
                    self.nn_model_type, self.num_features, **params
                ).to(self.device)
                
                temp_optimizer = torch.optim.Adam(
                    temp_model.parameters(), lr=0.001, weight_decay=1e-4
                )
                
                # Quick training (few epochs for evaluation)
                temp_model.train()
                X_train_tensor = torch.tensor(X_train_fold, dtype=torch.float32).to(self.device)
                y_train_tensor = torch.tensor(y_train_fold, dtype=torch.float32).to(self.device)
                
                for epoch in range(5):  # Quick evaluation with 5 epochs
                    temp_optimizer.zero_grad()
                    outputs = temp_model(X_train_tensor)
                    loss = self.loss_func(outputs, y_train_tensor.unsqueeze(-1))
                    loss.backward()
                    temp_optimizer.step()
                
                # Validation
                temp_model.eval()
                with torch.no_grad():
                    X_val_tensor = torch.tensor(X_val_fold, dtype=torch.float32).to(self.device)
                    y_val_tensor = torch.tensor(y_val_fold, dtype=torch.float32).to(self.device)
                    val_outputs = temp_model(X_val_tensor)
                    val_loss = self.loss_func(val_outputs, y_val_tensor.unsqueeze(-1)).item()
                    fold_scores.append(val_loss)
                
                del temp_model, temp_optimizer
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            
            # Average validation loss
            avg_score = np.mean(fold_scores)
            logging.debug(f"[HyperOpt]   Avg validation loss: {avg_score:.4f}")
            
            if avg_score < best_score:
                best_score = avg_score
                best_params = params
                logging.info(f"[HyperOpt] New best params: {best_params} (score={best_score:.4f})")
        
        logging.info(f"[HyperOpt] Best hyperparameters: {best_params}")
        logging.info(f"[HyperOpt] Best validation loss: {best_score:.4f}")
        
        return best_params

    def fit(self, features, labels, is_for_shap=False, enable_hyperparam_search=True):
        """Train model with optional hyperparameter search"""
        # Hyperparameter search
        if enable_hyperparam_search and not is_for_shap:
            best_params = self.grid_search_hyperparams(features, labels)
            
            # Reinitialize model with best hyperparameters
            self.nn_model = rescore_model_provider.get_model(
                self.nn_model_type, self.num_features, **best_params
            ).to(self.device)
            
            self.optimizer = torch.optim.Adam(
                self.nn_model.parameters(), lr=0.001, weight_decay=1e-4
            )
            
            if self.nn_model_type in ["dnn", "cnn"]:
                self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                    self.optimizer, mode='min', factor=0.5, patience=2
                )
        
        # Convert data to tensors
        features_tensor = torch.tensor(features, dtype=torch.float)
        labels_tensor = torch.tensor(labels, dtype=torch.float)
        
        # Create train/val split for DNN/CNN
        if self.nn_model_type in ["dnn", "cnn"]:
            dataset_size = len(features_tensor)
            indices = torch.randperm(dataset_size).tolist()
            train_size = int(0.9 * dataset_size)
            train_indices = indices[:train_size]
            val_indices = indices[train_size:]
            
            train_features = features_tensor[train_indices].to(self.device)
            train_labels = labels_tensor[train_indices].to(self.device)
            val_features = features_tensor[val_indices].to(self.device)
            val_labels = labels_tensor[val_indices].to(self.device)
            
            best_val_loss = float('inf')
            patience_counter = 0
            best_model_state = None
       
        # Training loop
        for epoch in range(self.epoch):
            self.nn_model.train()
            train_loss = 0.0
            num_batches = 0
            
            for i in range(0, len(train_features), self.train_batch_size):
                self.optimizer.zero_grad()
                
                batch_features = train_features[i:i+self.train_batch_size]
                batch_labels = train_labels[i:i+self.train_batch_size].unsqueeze(-1)
                
                outputs = self.nn_model(batch_features)
                loss = self.loss_func(outputs, batch_labels)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.nn_model.parameters(), 1.0)
                self.optimizer.step()
                
                train_loss += loss.item()
                num_batches += 1
            
            avg_train_loss = train_loss / num_batches if num_batches > 0 else 0
            
            # Validation
            self.nn_model.eval()
            val_loss = 0.0
            num_val_batches = 0
            
            with torch.no_grad():
                for i in range(0, len(val_features), self.train_batch_size):
                    batch_features = val_features[i:i+self.train_batch_size]
                    batch_labels = val_labels[i:i+self.train_batch_size].unsqueeze(-1)
                    
                    outputs = self.nn_model(batch_features)
                    loss = self.loss_func(outputs, batch_labels)
                    val_loss += loss.item()
                    num_val_batches += 1
            
            avg_val_loss = val_loss / num_val_batches if num_val_batches > 0 else 0
            
            self.scheduler.step(avg_val_loss)
            
            logging.debug(f"[{self.nn_model_type.upper()}] Epoch {epoch+1}/{self.epoch} - "
                         f"Train Loss: {avg_train_loss:.4f}, Val Loss: {avg_val_loss:.4f}")
            
            # Early stopping
            if avg_val_loss < best_val_loss:
                best_val_loss = avg_val_loss
                patience_counter = 0
                best_model_state = self.nn_model.state_dict()
            else:
                patience_counter += 1
                if patience_counter >= self.patience:
                    if is_for_shap:
                        logging.info(f"[SHAP] [{self.nn_model_type.upper()}] Early stopping at epoch {epoch+1}")
                    else:
                        logging.info(f"[Rescore] [{self.nn_model_type.upper()}] Early stopping at epoch {epoch+1}")
                    break
        
        if best_model_state is not None:
            self.nn_model.load_state_dict(best_model_state)

    def decision_function(self, features):
        self.nn_model.eval()
        outputs = np.empty(len(features))
        
        with torch.no_grad():
            for i in range(0, len(features), self.predict_batch_size):
                batch = torch.tensor(
                    features[i:i+self.predict_batch_size], 
                    dtype=torch.float, 
                    device=self.device
                )
                pred = self.nn_model(batch)
                
                if self.nn_model_type in ["dnn", "cnn"]:
                    pred = torch.sigmoid(pred)
                
                outputs[i:i+self.predict_batch_size] = pred.squeeze(1).cpu().numpy()
        return outputs

    def predict_proba(self, X):
        """For SHAP"""
        self.nn_model.eval()
        with torch.no_grad():
            if not isinstance(X, torch.Tensor):
                X = torch.tensor(X, dtype=torch.float32)
            X = X.to(self.device)
            preds = self.nn_model(X)
            preds = torch.sigmoid(preds)
            return preds.cpu().numpy()


# ==================== Main Percolator_ML Class ====================
class Percolator_ML:
    
    def __init__(
        self,
        *,
        percolator_model: str,
        percolator_backend: str,
        min_train_sample: int,
        max_train_sample: int,
        cv_fold: int,
        iter_num: int,
        fdr: float,
        fdr_level: str,  
        per_raw_fdr: bool,
        feature_list: list,
        shap_explain: bool = True,
        outpath: str,
        enable_hyperparam_search: bool = False,
    ):
        
        self.fdr_level = fdr_level
        self.fdr = fdr
        self.cv_fold = cv_fold
        self.iter_num = iter_num
        self.feature_list = feature_list
        self.max_train_sample = max_train_sample
        self.min_train_sample = min_train_sample
        self.per_raw_fdr = per_raw_fdr
        self.shap_explain = shap_explain
        self.outpath = outpath
        self.enable_hyperparam_search = enable_hyperparam_search

        self.init_percolator_model(percolator_model, percolator_backend)

    def init_percolator_model(
        self, percolator_model="svm", percolator_backend="sklearn"
    ):
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.linear_model import LogisticRegression
        from sklearn.naive_bayes import GaussianNB
        from sklearn.svm import LinearSVR
        from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
        from sklearn.preprocessing import StandardScaler

        self.percolator_model = percolator_model.lower()
        self.percolator_backend = percolator_backend.lower()

        # PyTorch models
        if percolator_backend.lower() == "pytorch":
            model_kwargs = {}
            if self.percolator_model == "dnn":
                model_kwargs = {
                    'hidden_dim1': 64,
                    'hidden_dim2': 16,
                    'dropout': 0.2
                }
                self.model = NNRescore(
                    len(self.feature_list), nn_model_type=percolator_model, **model_kwargs
                )

            elif self.percolator_model == "cnn":
                model_kwargs = {
                    'hidden_channels': 128,
                    'kernel_size': 3,
                    'dropout': 0.2
                }
                self.model = NNRescore(
                    len(self.feature_list), nn_model_type=percolator_model, **model_kwargs
                )
                self.predict_batch_size = 1000

            elif self.percolator_model == "nnlinear":
                self.model = NNLinearRescore(
                    len(self.feature_list), nn_model_type=percolator_model, **model_kwargs
                )

        # sklearn models with hyperparameter search
        elif self.percolator_model == "lr":
            self.scaler = StandardScaler()
            base_model = LogisticRegression(max_iter=1000, random_state=42)
            
            if self.enable_hyperparam_search:
                param_grid = get_hyperparam_grid("lr")
                self.model = GridSearchCV(
                    base_model,
                    param_grid=param_grid,
                    cv=KFold(3, shuffle=True, random_state=42),
                    scoring='neg_log_loss',
                    n_jobs=NUM_THREADS,
                    verbose=0
                )
                logging.info(f"[HyperOpt] LR with GridSearchCV: {param_grid}")
            else:
                self.model = base_model

        elif self.percolator_model == "svm":
            # Use dual=True for compatibility
            base_model = LinearSVR(
                dual='auto',
                loss='squared_epsilon_insensitive',
                max_iter=5000,
                random_state=42
            )
            
            if self.enable_hyperparam_search:
                param_grid = get_hyperparam_grid("svm")
                self.model = GridSearchCV(
                    base_model,
                    param_grid=param_grid,
                    cv=KFold(3, shuffle=True, random_state=42),
                    scoring='neg_mean_squared_error',
                    n_jobs=NUM_THREADS,
                    verbose=0
                )
                logging.info(f"[HyperOpt] SVM with GridSearchCV: {param_grid}")
            else:
                self.model = base_model

        elif self.percolator_model == "rf":
            base_model = RandomForestClassifier(random_state=42, n_jobs=NUM_THREADS)
            
            if self.enable_hyperparam_search:
                param_grid = get_hyperparam_grid("rf")
                self.model = GridSearchCV(
                    base_model,
                    param_grid=param_grid,
                    cv=KFold(3, shuffle=True, random_state=42),
                    scoring='neg_log_loss',
                    n_jobs=1,
                    verbose=1,
                    pre_dispatch='1*n_jobs'
                )
                logging.info(f"[HyperOpt] RF with GridSearchCV: {param_grid}")
            else:
                self.model = base_model

        elif self.percolator_model == "lda":
            self.scaler = StandardScaler()
            base_model = LinearDiscriminantAnalysis()
            
            if self.enable_hyperparam_search:
                param_grid = get_hyperparam_grid("lda")
                self.model = GridSearchCV(
                    base_model,
                    param_grid=param_grid,
                    cv=KFold(3, shuffle=True, random_state=42),
                    scoring='neg_log_loss',
                    n_jobs=NUM_THREADS,
                    verbose=0
                )
                logging.info(f"[HyperOpt] LDA with GridSearchCV: {param_grid}")
            else:
                self.model = base_model

        elif self.percolator_model == "bayes":
            self.scaler = StandardScaler()
            base_model = GaussianNB()
            
            if self.enable_hyperparam_search:
                param_grid = get_hyperparam_grid("bayes")
                self.model = GridSearchCV(
                    base_model,
                    param_grid=param_grid,
                    cv=KFold(3, shuffle=True, random_state=42),
                    scoring='neg_log_loss',
                    n_jobs=NUM_THREADS,
                    verbose=0
                )
                logging.info(f"[HyperOpt] Bayes with GridSearchCV: {param_grid}")
            else:
                self.model = base_model

        else:
            logging.info(
                "[Rescore] "
                f"Rescoring model '{percolator_model}' is not "
                "implemented, switch to sklearn 'svm' model."
            )
            self.model = LinearSVR(
                dual='auto',
                loss='squared_epsilon_insensitive',
                max_iter=5000,
                random_state=42
            )
            self.percolator_model = "svm"
            self.percolator_backend = "sklearn"

    def _estimate_fdr(
        self,
        df: pd.DataFrame,
        fdr_level: str = None,
        per_raw_fdr: bool = None,
    ) -> pd.DataFrame:

        if df.empty:
            return df

        df = df.sort_values(by=["ml_score", "decoy"], ascending=[False, True])
        df = df.reset_index(drop=True)
        
        if fdr_level is None:
            fdr_level = self.fdr_level
        if per_raw_fdr is None:
            per_raw_fdr = self.per_raw_fdr
        if per_raw_fdr:
            df_list = []
            for raw_name, df_raw in df.groupby("raw_name"):
                if not df_raw.empty:
                    df_list.append(
                        self._estimate_fdr(df_raw, fdr_level=fdr_level, per_raw_fdr=False)
                    )
            return pd.concat(df_list) if df_list else df
       
        if fdr_level == "precursor":
            group_key = "precursor_id"
        elif fdr_level == "peptide":
            group_key = "sequence"
        else:
            group_key = "precursor_id"

        _df = df.groupby([group_key, "decoy"], as_index=False)["ml_score"].max()
        _df = _df.sort_values(["ml_score", "decoy"], ascending=[False, True])

        target_values = 1 - _df["decoy"].values
        decoy_cumsum = np.cumsum(_df["decoy"].values)
        target_cumsum = np.cumsum(target_values)

        with np.errstate(divide='ignore', invalid='ignore'):
            fdr_values = np.where(
                target_cumsum > 0, 
                decoy_cumsum / target_cumsum, 
                1 
            )
        _df["fdr"] = fdr_to_q_values(fdr_values)
        
        df["fdr"] = fdr_from_ref(
            df["ml_score"].values, _df["ml_score"].values, _df["fdr"].values
        )
        return df

    def _train(self, train_t_df, train_d_df):
        if self.percolator_model in ["cnn", "nnlinear"]:
        # 训练前清理缓存
            torch.cuda.empty_cache()

        if train_t_df.empty or train_d_df.empty:
            logging.warning("[Rescore] Empty training set. Skipping training.")
            return
        
        if len(train_t_df) > self.max_train_sample:
            train_t_df = train_t_df.sample(n=self.max_train_sample, random_state=1337)
        if len(train_d_df) > self.max_train_sample:
            train_d_df = train_d_df.sample(n=self.max_train_sample, random_state=1337)

        train_df = pd.concat((train_t_df, train_d_df))
        train_label = np.ones(len(train_df), dtype=np.int32)
        train_label[len(train_t_df):] = 0

        X_train = train_df[self.feature_list].values

        # Standardization (only for models that need it)
        if self.percolator_model in ["lr", "lda", "bayes"]:
            X_train = self.scaler.fit_transform(X_train)

        # PyTorch models
        if self.percolator_backend == "pytorch":
            self.model.fit(
                X_train, 
                train_label, 
                enable_hyperparam_search=self.enable_hyperparam_search
            )
        
        # sklearn models
        else:
            self.model.fit(X_train, train_label)
            
            # Log best parameters if GridSearchCV was used
            if isinstance(self.model, GridSearchCV):
                logging.info(f"[HyperOpt] Best params: {self.model.best_params_}")
                logging.info(f"[HyperOpt] Best score: {self.model.best_score_:.4f}")

        # Log OOB score for RF
        if self.percolator_model == "rf":
            final_model = self.model.best_estimator_ if isinstance(self.model, GridSearchCV) else self.model
            if hasattr(final_model, 'oob_score_') and hasattr(final_model, 'oob_score') and final_model.oob_score:
                logging.info(f"[RF] OOB Score: {final_model.oob_score_:.4f}")
        
        if self.percolator_model in ["cnn", "nnlinear"]:
            torch.cuda.empty_cache()

    def _predict(self, test_df):
        X_test = test_df[self.feature_list].values

        if self.percolator_backend == "pytorch":
            test_df["ml_score"] = self.model.decision_function(X_test)

        elif self.percolator_model == "lr":
            if hasattr(self, 'scaler'):
                X_test = self.scaler.transform(X_test)
            test_df["ml_score"] = self.model.decision_function(X_test)

        elif self.percolator_model == "bayes":
            if hasattr(self, 'scaler'):
                X_test = self.scaler.transform(X_test)
            final_model = self.model.best_estimator_ if isinstance(self.model, GridSearchCV) else self.model
            test_df["ml_score"] = final_model.predict_proba(X_test)[:, 1]
        
        elif self.percolator_model == "lda":
            if hasattr(self, 'scaler'):
                X_test = self.scaler.transform(X_test)
            final_model = self.model.best_estimator_ if isinstance(self.model, GridSearchCV) else self.model
            test_df["ml_score"] = final_model.predict_proba(X_test)[:, 1]

        elif self.percolator_model == "svm":
            final_model = self.model.best_estimator_ if isinstance(self.model, GridSearchCV) else self.model
            test_df["ml_score"] = final_model.predict(X_test)

        elif self.percolator_model == "rf":
            final_model = self.model.best_estimator_ if isinstance(self.model, GridSearchCV) else self.model
            test_df["ml_score"] = final_model.predict_proba(X_test)[:, 1]
        
        return test_df

    def _cv_score(self, df: pd.DataFrame) -> pd.DataFrame:
        """Original FDR-based training set selection strategy"""
        if df.empty:
            return df
            
        df = df.sample(frac=1, random_state=1337).reset_index(drop=True)
        df_target = df[df.decoy == 0]
        df_decoy = df[df.decoy != 0]
        if (
            np.sum(df_target.fdr < 0.01) < self.min_train_sample * self.cv_fold
            or len(df_decoy) < self.min_train_sample * self.cv_fold
        ):
            logging.info(
                "[Rescore] "
                f"#target={np.sum(df_target.fdr<0.01)} or #decoy={len(df_decoy)} "
                f"< minimal training sample={self.min_train_sample} "
                f"for cv-fold={self.cv_fold}. Skip rescoring!!!"
            )
            return df

        if self.cv_fold > 1:
            test_df_list = []
            for i in range(self.cv_fold):
                t_mask = np.ones(len(df_target), dtype=bool)
                _slice = slice(i, len(df_target), self.cv_fold)
                t_mask[_slice] = False
                cv_df_target = df_target[t_mask]
                train_t_df = cv_df_target[cv_df_target.fdr <= self.fdr]
                test_t_df = df_target[_slice]

                d_mask = np.ones(len(df_decoy), dtype=bool)
                _slice = slice(i, len(df_decoy), self.cv_fold)
                d_mask[_slice] = False
                train_d_df = df_decoy[d_mask]
                test_d_df = df_decoy[_slice]

                self._train(train_t_df, train_d_df)

                test_df = pd.concat((test_t_df, test_d_df))
                test_df_list.append(self._predict(test_df))

            return pd.concat(test_df_list)
        else:
            train_t_df = df_target[df_target.fdr <= self.fdr]

            self._train(train_t_df, df_decoy)
            test_df = pd.concat((df_target, df_decoy))

            return self._predict(test_df)

    def _cv_score_topk(self, df: pd.DataFrame, top_ratio: float = 0.3) -> pd.DataFrame:
        """Top-K training set selection strategy"""
        if df.empty:
            return df
            
        df = df.sample(frac=1, random_state=1337).reset_index(drop=True)
        
        df_target = df[df.decoy == 0].copy()
        df_decoy = df[df.decoy != 0].copy()
        
        if (
            len(df_target) < self.min_train_sample * self.cv_fold
            or len(df_decoy) < self.min_train_sample * self.cv_fold
        ):
            logging.warning(
                "[Rescore] "
                f"#target={len(df_target)} or #decoy={len(df_decoy)} "
                f"< minimal training sample={self.min_train_sample} "
                f"for cv-fold={self.cv_fold}. Skip rescoring!!!"
            )
            return df

        # Split high-quality and low-quality sets
        df_target_sorted = df_target.sort_values('ml_score', ascending=False).reset_index(drop=True)
        n_target_hq = int(len(df_target_sorted) * top_ratio)
        df_target_hq = df_target_sorted.iloc[:n_target_hq].copy()
        df_target_lq = df_target_sorted.iloc[n_target_hq:].copy()
        
        df_decoy_sorted = df_decoy.sort_values('ml_score', ascending=False).reset_index(drop=True)
        n_decoy_hq = int(len(df_decoy_sorted) * top_ratio)
        df_decoy_hq = df_decoy_sorted.iloc[:n_decoy_hq].copy()
        df_decoy_lq = df_decoy_sorted.iloc[n_decoy_hq:].copy()
        
        logging.info(
            f"[Rescore] Data split (top_ratio={top_ratio}):\n"
            f"  Target: HQ={len(df_target_hq)}, LQ={len(df_target_lq)}\n"
            f"  Decoy:  HQ={len(df_decoy_hq)}, LQ={len(df_decoy_lq)}"
        )
        
        test_df_list = []
        
        for i in range(self.cv_fold):
            logging.info(f"[Rescore] Fold {i+1}/{self.cv_fold}")
            
            # Create CV masks for high-quality sets
            hq_t_mask = np.ones(len(df_target_hq), dtype=bool)
            hq_t_slice = slice(i, len(df_target_hq), self.cv_fold)
            hq_t_mask[hq_t_slice] = False
            
            train_t_hq = df_target_hq[hq_t_mask]
            test_t_hq = df_target_hq[~hq_t_mask]
            
            hq_d_mask = np.ones(len(df_decoy_hq), dtype=bool)
            hq_d_slice = slice(i, len(df_decoy_hq), self.cv_fold)
            hq_d_mask[hq_d_slice] = False
            
            train_d_hq = df_decoy_hq[hq_d_mask]
            test_d_hq = df_decoy_hq[~hq_d_mask]
            
            # Create CV masks for low-quality sets (test only)
            lq_t_mask = np.zeros(len(df_target_lq), dtype=bool)
            lq_t_slice = slice(i, len(df_target_lq), self.cv_fold)
            lq_t_mask[lq_t_slice] = True
            
            test_t_lq = df_target_lq[lq_t_mask]
            
            lq_d_mask = np.zeros(len(df_decoy_lq), dtype=bool)
            lq_d_slice = slice(i, len(df_decoy_lq), self.cv_fold)
            lq_d_mask[lq_d_slice] = True
            
            test_d_lq = df_decoy_lq[lq_d_mask]
            
            # Balance training set
            min_train_samples = min(len(train_t_hq), len(train_d_hq))
            min_train_samples = min(min_train_samples, self.max_train_sample)
            
            if min_train_samples < self.min_train_sample:
                logging.warning(
                    f"[Rescore] Fold {i}: Only {min_train_samples} samples available "
                    f"(< {self.min_train_sample}). Using all available HQ samples."
                )
                min_train_samples = min(len(train_t_hq), len(train_d_hq))
            
            train_t_df = train_t_hq.iloc[:min_train_samples]
            train_d_df = train_d_hq.iloc[:min_train_samples]
            
            # Combine test sets
            test_t_df = pd.concat([test_t_hq, test_t_lq])
            test_d_df = pd.concat([test_d_hq, test_d_lq])
            
            logging.info(
                f"  Fold {i+1} split:\n"
                f"    Train: Target={len(train_t_df)}, Decoy={len(train_d_df)}\n"
                f"    Test:  Target={len(test_t_df)} (HQ={len(test_t_hq)}, LQ={len(test_t_lq)}), "
                f"Decoy={len(test_d_df)} (HQ={len(test_d_hq)}, LQ={len(test_d_lq)})"
            )
            
            self._train(train_t_df, train_d_df)
            
            test_df = pd.concat([test_t_df, test_d_df])
            test_df_list.append(self._predict(test_df))
        
        result_df = pd.concat(test_df_list)
        
        assert len(result_df) == len(df), \
            f"Data loss detected! Original: {len(df)}, After CV: {len(result_df)}"
        
        logging.info(f"[Rescore] CV completed. All {len(result_df)} samples tested.")
        
        return result_df

    def explain_model_with_shap(self, df: pd.DataFrame):
        """SHAP model explanation"""
        if not self.shap_explain:
            logging.info("[SHAP] is not used")
            return
        
        try:    
            train_t_df = df[(df.fdr <= self.fdr) & (df.decoy == 0)]
            train_d_df = df[df.decoy != 0]
            
            if len(train_t_df) < self.min_train_sample or len(train_d_df) < self.min_train_sample:
                logging.warning(f"[SHAP] sample is not sufficient (target={len(train_t_df)}, decoy={len(train_d_df)}), skipping SHAP")
                return
                
            n_sample = min(len(train_t_df), len(train_d_df), self.max_train_sample)
            train_t_df = train_t_df.sample(n=n_sample, random_state=42)
            train_d_df = train_d_df.sample(n=n_sample, random_state=42)
            
            train_df = pd.concat([train_t_df, train_d_df])
            train_labels = np.concatenate([
                np.ones(len(train_t_df), dtype=int),
                np.zeros(len(train_d_df), dtype=int)
            ])
            
            logging.info("[SHAP] training final model for explanation...")
            if self.percolator_model in ["cnn", 'dnn']:
                self.model.fit(train_df[self.feature_list].values, train_labels, is_for_shap=True)
            else:
                X_train = train_df[self.feature_list].values
                if self.percolator_model in ["lr", "lda", "bayes"]:
                    X_train = self.scaler.fit_transform(X_train)
                self.model.fit(X_train, train_labels)
            
            explain_df = pd.concat([
                train_t_df.sample(min(500, len(train_t_df))), 
                train_d_df.sample(min(500, len(train_d_df)))
            ])
            shap_values = None

            def silent_call(func):
                """Execute func with stdout/stderr redirected"""
                with open(os.devnull, 'w') as devnull_out, open(os.devnull, 'w') as devnull_err:
                    old_out, old_err = sys.stdout, sys.stderr
                    try:
                        sys.stdout, sys.stderr = devnull_out, devnull_err
                        shap_logger = logging.getLogger('shap')
                        original_shap_level = shap_logger.level
                        shap_logger.setLevel(logging.CRITICAL)
                        
                        np_logger = logging.getLogger('numpy')
                        original_numpy_level = np_logger.level
                        np_logger.setLevel(logging.CRITICAL)

                        os.environ['PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION'] = 'python'
                        os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'
                        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

                        return func()
                    finally:
                        sys.stdout, sys.stderr = old_out, old_err
                        shap_logger.setLevel(original_shap_level)
                        np_logger.setLevel(original_numpy_level)
                
            def _compute_shap():
                if self.percolator_backend == "pytorch":
                    logging.info("[SHAP] using DeepExplainer to explain PyTorch-based model")

                    pytorch_model = self.model.nn_model
                    background_data = train_df[self.feature_list].values[:100]
                    explain_data = explain_df[self.feature_list].values
                    
                    background_tensor = torch.tensor(background_data, dtype=torch.float32).to(self.model.device)
                    explain_tensor = torch.tensor(explain_data, dtype=torch.float32).to(self.model.device)

                    if self.model.nn_model_type == "cnn":
                        background_tensor = background_tensor.unsqueeze(1)
                        explain_tensor = explain_tensor.unsqueeze(1)
    
                    explainer = shap.DeepExplainer(pytorch_model, background_tensor)
                    shap_values = explainer.shap_values(explain_tensor)
                    
                    if isinstance(shap_values, list):
                        shap_values = shap_values[0]
                    
                    if self.model.nn_model_type == "cnn":
                        shap_values = shap_values.squeeze(1)
                    shap_values = shap_values.squeeze(-1)
                    
                elif self.percolator_model == "rf":
                    logging.info("[SHAP] using TreeExplainer to explain tree-based model")
                    final_model = self.model.best_estimator_ if isinstance(self.model, GridSearchCV) else self.model
                    explainer = shap.TreeExplainer(final_model)
                    shap_values = explainer.shap_values(explain_df[self.feature_list].values)
                    
                elif self.percolator_model == "lr":
                    logging.info("[SHAP] using LinearExplainer to explain linear model")
                    final_model = self.model.best_estimator_ if isinstance(self.model, GridSearchCV) else self.model
                    background = self.scaler.transform(train_df[self.feature_list].values)
                    explain_data = self.scaler.transform(explain_df[self.feature_list].values)
                    explainer = shap.LinearExplainer(final_model, background)
                    shap_values = explainer.shap_values(explain_data)

                elif self.percolator_model == "lda":
                    logging.info("[SHAP] using LinearExplainer to explain LDA model")
                    final_model = self.model.best_estimator_ if isinstance(self.model, GridSearchCV) else self.model
                    background = self.scaler.transform(train_df[self.feature_list].values)
                    explain_data = self.scaler.transform(explain_df[self.feature_list].values)
                    explainer = shap.LinearExplainer(
                        final_model, background,
                        feature_perturbation="interventional"
                    )
                    shap_values = explainer.shap_values(explain_data)
                    
                else:
                    logging.info(f"[SHAP] using KernelExplainer to explain {self.percolator_model} model")
                    final_model = self.model.best_estimator_ if isinstance(self.model, GridSearchCV) else self.model
                    explainer = shap.KernelExplainer(
                        final_model.predict_proba if hasattr(final_model, 'predict_proba') else final_model.predict,
                        shap.sample(train_df[self.feature_list].values, 100)
                    )
                    shap_values = explainer.shap_values(explain_df[self.feature_list].values)

                return shap_values
            
            shap_values = silent_call(_compute_shap)
            
            if shap_values is not None:
                self._visualize_shap(shap_values, explain_df)

        except Exception as e:
            logging.error(f"[SHAP] Failed: {str(e)}")
            import traceback
            logging.debug(traceback.format_exc())

    
    def _visualize_shap(self, shap_values, explain_df):
        """Visualizing SHAP with multiple plot types"""
        try:
            # ============================================================
            # 数据预处理 - 处理不同模型的SHAP输出格式
            # ============================================================
            logging.info(f"[SHAP] Raw SHAP values type: {type(shap_values)}")
            
            # 处理PyTorch tensor（先转换，因为后续需要检查shape）
            if torch.is_tensor(shap_values):
                shap_values = shap_values.detach().cpu().numpy()
                logging.info("[SHAP] Converted PyTorch tensor to NumPy array")
            
            # 处理列表格式（某些版本的RF/分类器）
            if isinstance(shap_values, list):
                logging.info(f"[SHAP] SHAP values is a list with {len(shap_values)} elements")
                
                if len(shap_values) == 1:
                    shap_values_for_plot = shap_values[0]
                    logging.info("[SHAP] Using single class SHAP values")
                
                elif len(shap_values) == 2:
                    # 二分类：class 0 = decoy, class 1 = target
                    shap_values_for_plot = shap_values[1]
                    logging.info("[SHAP] Binary classification (list format), using class 1 (target)")
                
                else:
                    shap_values_for_plot = shap_values[1] if len(shap_values) > 1 else shap_values[0]
                    logging.warning(
                        f"[SHAP] Multi-class output ({len(shap_values)} classes), using class 1"
                    )
            
            # 处理NumPy数组格式
            elif isinstance(shap_values, np.ndarray):
                logging.info(f"[SHAP] SHAP values shape: {shap_values.shape}")
                
                # 情况1: 2D数组 (samples, features) - 理想格式
                if shap_values.ndim == 2:
                    shap_values_for_plot = shap_values
                    logging.info("[SHAP] 2D array format, using directly")
                
                # 情况2: 3D数组 (samples, features, classes) - RF常见格式
                elif shap_values.ndim == 3:
                    n_classes = shap_values.shape[2]
                    logging.info(f"[SHAP] 3D array format with {n_classes} classes")
                    
                    if n_classes == 2:
                        # 二分类：取正类（class 1 = target）
                        shap_values_for_plot = shap_values[:, :, 1]
                        logging.info("[SHAP] Binary classification (3D format), extracting class 1 (target)")
                    
                    elif n_classes == 1:
                        # 单类别（降维）
                        shap_values_for_plot = shap_values[:, :, 0]
                        logging.info("[SHAP] Single class in 3D format, squeezing")
                    
                    else:
                        # 多分类（取第一个或最后一个类）
                        shap_values_for_plot = shap_values[:, :, -1]  # 或者用 [:, :, 0]
                        logging.warning(
                            f"[SHAP] Multi-class (3D format, {n_classes} classes), "
                            f"using class {n_classes-1} for visualization"
                        )
                
                # 情况3: 4D数组（CNN可能产生，需要squeeze）
                elif shap_values.ndim == 4:
                    logging.warning(f"[SHAP] 4D array detected: {shap_values.shape}, attempting to squeeze")
                    shap_values_for_plot = shap_values.squeeze()
                    
                    if shap_values_for_plot.ndim == 3:
                        # Squeeze后变成3D，递归处理
                        logging.info(f"[SHAP] After squeeze: {shap_values_for_plot.shape}, extracting class dimension")
                        shap_values_for_plot = shap_values_for_plot[:, :, 1] if shap_values_for_plot.shape[2] == 2 else shap_values_for_plot[:, :, 0]
                    
                    elif shap_values_for_plot.ndim != 2:
                        raise ValueError(f"Cannot handle shape {shap_values_for_plot.shape} after squeeze")
                
                else:
                    raise ValueError(f"Unsupported SHAP values dimension: {shap_values.ndim}")
            
            else:
                raise TypeError(f"Unsupported SHAP values type: {type(shap_values)}")

            # 最终验证
            if shap_values_for_plot.ndim != 2:
                raise ValueError(
                    f"SHAP values must be 2D after processing, "
                    f"but got shape {shap_values_for_plot.shape}"
                )
            
            # 验证特征数量
            if shap_values_for_plot.shape[1] != len(self.feature_list):
                raise ValueError(
                    f"SHAP values feature dimension ({shap_values_for_plot.shape[1]}) "
                    f"does not match feature list length ({len(self.feature_list)})"
                )

            X_display = explain_df[self.feature_list].values
            
            logging.info(f"[SHAP] ✓ Processed SHAP values shape: {shap_values_for_plot.shape}")
            logging.info(f"[SHAP] ✓ X_display shape: {X_display.shape}")
            logging.info(f"[SHAP] ✓ Number of features: {len(self.feature_list)}")
            
            # ============================================================
            # 1. Feature Importance Bar Plot (原有)
            # ============================================================
            logging.info("[SHAP] Generating feature importance plot...")
            plt.figure(figsize=(12, 6))
            shap.summary_plot(
                shap_values_for_plot, 
                X_display,
                feature_names=self.feature_list,
                plot_type="bar", 
                show=False
            )
            plt.title(f"{self.percolator_model.upper()} - Feature Importance", fontsize=14, fontweight='bold')
            plt.tight_layout()
            
            plt.savefig(self.outpath / f"{self.percolator_model}_shap_feature_importance.png", dpi=300, bbox_inches='tight')
            plt.savefig(self.outpath / f"{self.percolator_model}_shap_feature_importance.pdf", dpi=300, bbox_inches='tight')
            plt.close()
            logging.info(f"[SHAP] ✓ Feature importance saved")
            
            # ============================================================
            # 2. Summary Plot (原有)
            # ============================================================
            logging.info("[SHAP] Generating summary plot...")
            plt.figure(figsize=(12, 8))
            shap.summary_plot(
                shap_values_for_plot, 
                X_display,
                feature_names=self.feature_list,
                show=False
            )
            plt.title(f"{self.percolator_model.upper()} - SHAP Values Distribution", fontsize=14, fontweight='bold')
            plt.tight_layout()
            
            plt.savefig(self.outpath / f"{self.percolator_model}_shap_summary_plot.png", dpi=300, bbox_inches='tight')
            plt.savefig(self.outpath / f"{self.percolator_model}_shap_summary_plot.pdf", dpi=300, bbox_inches='tight')
            plt.close()
            logging.info(f"[SHAP] ✓ Summary plot saved")
            
            # ============================================================
            # 3. Dependence Plots (依赖图) - 展示特征与SHAP值的关系
            # ============================================================
            logging.info("[SHAP] Generating dependence plots for top features...")
            
            # 计算特征重要性（按SHAP值绝对值的平均值排序）
            feature_importance = np.abs(shap_values_for_plot).mean(axis=0)
            top_features_idx = np.argsort(feature_importance)[::-1][:5]  # 前5个最重要特征
            
            for idx in top_features_idx:
                feature_name = self.feature_list[idx]
                
                plt.figure(figsize=(10, 6))
                shap.dependence_plot(
                    idx,
                    shap_values_for_plot,
                    X_display,
                    feature_names=self.feature_list,
                    show=False,
                    interaction_index='auto'  # 自动选择最强交互特征
                )
                plt.title(f"{self.percolator_model.upper()} - Dependence Plot: {feature_name}", 
                        fontsize=12, fontweight='bold')
                plt.tight_layout()
                
                safe_feature_name = feature_name.replace('/', '_').replace(' ', '_')
                plt.savefig(
                    self.outpath / f"{self.percolator_model}_shap_dependence_{safe_feature_name}.png", 
                    dpi=300, bbox_inches='tight'
                )
                plt.close()
            
            logging.info(f"[SHAP] ✓ Dependence plots saved for top {len(top_features_idx)} features")
            
            # ============================================================
            # 4. Force Plots (力图) - 展示单个预测的解释
            # ============================================================

            logging.info("[SHAP] Generating force plots...")

            # 计算base_value（所有样本SHAP值的期望）
            if self.percolator_backend == "pytorch":
                base_value = shap_values_for_plot.mean(axis=0).sum()
            else:
                base_value = shap_values_for_plot.mean(axis=0).sum()

            # 选择几个代表性样本
            # 1) 预测为正例且置信度最高的样本
            target_samples = explain_df[explain_df['decoy'] == 0]
            if len(target_samples) > 0:
                high_conf_target_idx = target_samples.nlargest(1, 'ml_score').index[0]
                sample_idx_in_explain = explain_df.index.get_loc(high_conf_target_idx)
                
                plt.figure(figsize=(20, 4))  # 增加高度
                shap.force_plot(
                    base_value,
                    shap_values_for_plot[sample_idx_in_explain],
                    X_display[sample_idx_in_explain],
                    feature_names=self.feature_list,
                    matplotlib=True,
                    show=False,
                    text_rotation=10,  # 文字旋转角度
                    contribution_threshold=0.01  # 只显示贡献度>1%的特征值
                )
                plt.title(f"{self.percolator_model.upper()} - Force Plot: High Confidence Target", 
                        fontsize=12, fontweight='bold', pad=10)
                # 调整字体大小
                ax = plt.gca()
                for text in ax.texts:
                    text.set_fontsize(7)  # 减小字体
                plt.tight_layout()
                plt.savefig(
                    self.outpath / f"{self.percolator_model}_shap_force_high_target.png", 
                    dpi=300, bbox_inches='tight'
                )
                plt.close()

            # 2) 预测为负例（decoy）的样本
            decoy_samples = explain_df[explain_df['decoy'] != 0]
            if len(decoy_samples) > 0:
                high_conf_decoy_idx = decoy_samples.nlargest(1, 'ml_score').index[0]
                sample_idx_in_explain = explain_df.index.get_loc(high_conf_decoy_idx)
                
                plt.figure(figsize=(20, 4))
                shap.force_plot(
                    base_value,
                    shap_values_for_plot[sample_idx_in_explain],
                    X_display[sample_idx_in_explain],
                    feature_names=self.feature_list,
                    matplotlib=True,
                    show=False,
                    text_rotation=10,
                    contribution_threshold=0.01
                )
                plt.title(f"{self.percolator_model.upper()} - Force Plot: High Confidence Decoy", 
                        fontsize=12, fontweight='bold', pad=10)
                ax = plt.gca()
                for text in ax.texts:
                    text.set_fontsize(7)
                plt.tight_layout()
                plt.savefig(
                    self.outpath / f"{self.percolator_model}_shap_force_high_decoy.png", 
                    dpi=300, bbox_inches='tight'
                )
                plt.close()

            # 3) 边界样本（预测分数接近阈值）
            median_score = explain_df['ml_score'].median()
            boundary_samples = explain_df.iloc[(explain_df['ml_score'] - median_score).abs().argsort()[:1]]
            if len(boundary_samples) > 0:
                boundary_idx = boundary_samples.index[0]
                sample_idx_in_explain = explain_df.index.get_loc(boundary_idx)
                
                plt.figure(figsize=(20, 4))
                shap.force_plot(
                    base_value,
                    shap_values_for_plot[sample_idx_in_explain],
                    X_display[sample_idx_in_explain],
                    feature_names=self.feature_list,
                    matplotlib=True,
                    show=False,
                    text_rotation=10,
                    contribution_threshold=0.01
                )
                plt.title(f"{self.percolator_model.upper()} - Force Plot: Boundary Sample", 
                        fontsize=12, fontweight='bold', pad=10)
                ax = plt.gca()
                for text in ax.texts:
                    text.set_fontsize(7)
                plt.tight_layout()
                plt.savefig(
                    self.outpath / f"{self.percolator_model}_shap_force_boundary.png", 
                    dpi=300, bbox_inches='tight'
                )
                plt.close()

            logging.info(f"[SHAP] ✓ Force plots saved")
            # ============================================================
            # 4. Interaction Summary Plot (交互作用摘要图) - 聚焦Top特征
            # ============================================================

            logging.info("[SHAP] Generating interaction summary plot for top features...")

            try:
                # 首先确定Top 5重要特征
                feature_importance = np.abs(shap_values_for_plot).mean(axis=0)
                top_5_indices = np.argsort(feature_importance)[::-1][:5]
                top_5_features = [self.feature_list[i] for i in top_5_indices]
                
                logging.info(f"[SHAP] Top 5 features for interaction analysis:")
                for i, (idx, feat) in enumerate(zip(top_5_indices, top_5_features)):
                    logging.info(f"[SHAP]   {i+1}. {feat} (importance: {feature_importance[idx]:.4f})")
                
                # 计算交互作用SHAP值
                logging.info("[SHAP] Computing interaction values (this may take a while)...")
                
                # 限制样本数量以加快计算（最多1000个样本）
                sample_size = min(1000, len(X_display))
                if sample_size < len(X_display):
                    np.random.seed(42)
                    sample_indices = np.random.choice(len(X_display), sample_size, replace=False)
                    X_sample = X_display[sample_indices]
                    logging.info(f"[SHAP] Using {sample_size} samples for interaction calculation")
                else:
                    X_sample = X_display
                    sample_indices = np.arange(len(X_display))
                
                shap_interaction_values = None
                
                if self.percolator_backend == "pytorch":
                    # PyTorch模型使用DeepExplainer
                    logging.info(f"[SHAP] Using DeepExplainer for PyTorch model")
                    
                    # 获取实际的PyTorch模型
                    pytorch_model = self.model.nn_model
                    
                    # 准备background数据
                    background_data = X_display[:100]
                    background_tensor = torch.tensor(background_data, dtype=torch.float32).to(self.model.device)
                    
                    # 准备解释数据
                    X_sample_tensor = torch.tensor(X_sample, dtype=torch.float32).to(self.model.device)
                    
                    # CNN需要额外处理维度
                    if self.model.nn_model_type == "cnn":
                        background_tensor = background_tensor.unsqueeze(1)
                        X_sample_tensor = X_sample_tensor.unsqueeze(1)
                    
                    # 创建explainer
                    explainer_interaction = shap.DeepExplainer(pytorch_model, background_tensor)
                    
                    # 计算交互值
                    shap_interaction_values = explainer_interaction.shap_interaction_values(X_sample_tensor)
                    
                    # 处理输出格式
                    if isinstance(shap_interaction_values, list):
                        shap_interaction_values = shap_interaction_values[0] if len(shap_interaction_values) == 1 else shap_interaction_values[1]
                    
                    if torch.is_tensor(shap_interaction_values):
                        shap_interaction_values = shap_interaction_values.detach().cpu().numpy()
                    
                    # CNN需要去掉channel维度
                    if self.model.nn_model_type == "cnn" and shap_interaction_values.ndim == 4:
                        shap_interaction_values = shap_interaction_values.squeeze(1)
                        
                else:
                    # sklearn模型 - 获取实际模型（可能包装在GridSearchCV中）
                    if isinstance(self.model, GridSearchCV):
                        final_model = self.model.best_estimator_
                    else:
                        final_model = self.model
                    
                    model_type = type(final_model).__name__
                    
                    if model_type in ['RandomForestClassifier', 'GradientBoostingClassifier', 
                                    'XGBClassifier', 'LGBMClassifier', 'CatBoostClassifier',
                                    'ExtraTreesClassifier', 'HistGradientBoostingClassifier']:
                        # 树模型支持交互计算
                        logging.info(f"[SHAP] Using TreeExplainer for {model_type}")
                        explainer_interaction = shap.TreeExplainer(final_model)
                        
                        # 准备数据（可能需要标准化）
                        if hasattr(self, 'scaler') and self.percolator_model in ["lr", "lda", "bayes"]:
                            X_sample_scaled = self.scaler.transform(X_sample)
                        else:
                            X_sample_scaled = X_sample
                        
                        shap_interaction_values = explainer_interaction.shap_interaction_values(X_sample_scaled)
                        
                        # 如果是多分类，取正类的交互值
                        if isinstance(shap_interaction_values, list):
                            shap_interaction_values = shap_interaction_values[1]
                    
                    elif model_type in ['LogisticRegression', 'LinearDiscriminantAnalysis']:
                        # 线性模型 - 使用近似方法
                        logging.warning(f"[SHAP] {model_type} doesn't support direct interaction calculation")
                        logging.info("[SHAP] Using correlation-based approximation")
                        shap_interaction_values = None
                        
                    else:
                        # 其他模型尝试使用KernelExplainer（很慢）
                        logging.warning(f"[SHAP] Model type {model_type} may not support interaction calculation")
                        logging.info("[SHAP] Attempting KernelExplainer (this will be VERY slow)...")
                        
                        try:
                            # 进一步减少样本
                            kernel_sample_size = min(100, sample_size)
                            kernel_indices = np.random.choice(len(X_sample), kernel_sample_size, replace=False)
                            X_kernel_sample = X_sample[kernel_indices]
                            
                            # 准备数据
                            if hasattr(self, 'scaler') and self.percolator_model in ["lr", "lda", "bayes"]:
                                X_kernel_sample = self.scaler.transform(X_kernel_sample)
                                background_kernel = shap.kmeans(self.scaler.transform(X_display), 10)
                            else:
                                background_kernel = shap.kmeans(X_display, 10)
                            
                            # 定义预测函数
                            if hasattr(final_model, 'predict_proba'):
                                predict_fn = lambda x: final_model.predict_proba(x)[:, 1]
                            elif hasattr(final_model, 'decision_function'):
                                predict_fn = final_model.decision_function
                            else:
                                predict_fn = final_model.predict
                            
                            explainer_interaction = shap.KernelExplainer(predict_fn, background_kernel)
                            
                            # KernelExplainer不支持interaction_values，设为None
                            logging.warning("[SHAP] KernelExplainer doesn't support shap_interaction_values()")
                            shap_interaction_values = None
                            
                        except Exception as e:
                            logging.error(f"[SHAP] KernelExplainer failed: {str(e)}")
                            shap_interaction_values = None
                
                # ========================================================
                # 方法1: 完整的交互作用摘要图（如果成功计算了交互值）
                # ========================================================
                if shap_interaction_values is not None:
                    logging.info("[SHAP] Drawing full interaction summary plot...")
                    
                    plt.figure(figsize=(14, 12))
                    shap.summary_plot(
                        shap_interaction_values,
                        X_sample,
                        feature_names=self.feature_list,
                        show=False,
                        max_display=20
                    )
                    plt.title(f"{self.percolator_model.upper()} - SHAP Interaction Values (All Features)", 
                            fontsize=14, fontweight='bold')
                    plt.tight_layout()
                    plt.savefig(
                        self.outpath / f"{self.percolator_model}_shap_interaction_all.png", 
                        dpi=300, bbox_inches='tight'
                    )
                    plt.savefig(
                        self.outpath / f"{self.percolator_model}_shap_interaction_all.pdf", 
                        dpi=300, bbox_inches='tight'
                    )
                    plt.close()
                    logging.info(f"[SHAP] ✓ Full interaction summary saved")
                    
                    # ========================================================
                    # 方法2: Top 5特征的交互矩阵热图
                    # ========================================================
                    logging.info("[SHAP] Creating interaction matrix for top 5 features...")
                    
                    # 提取Top 5特征的交互矩阵
                    # shap_interaction_values shape: (n_samples, n_features, n_features)
                    interaction_matrix_top5 = np.zeros((5, 5))
                    
                    for i, idx_i in enumerate(top_5_indices):
                        for j, idx_j in enumerate(top_5_indices):
                            # 取所有样本的平均交互强度
                            interaction_matrix_top5[i, j] = np.abs(
                                shap_interaction_values[:, idx_i, idx_j]
                            ).mean()
                    
                    # 绘制热图
                    plt.figure(figsize=(10, 8))
                    import seaborn as sns
                    
                    mask = np.triu(np.ones_like(interaction_matrix_top5, dtype=bool), k=1)  # 只显示下三角
                    
                    sns.heatmap(
                        interaction_matrix_top5,
                        xticklabels=top_5_features,
                        yticklabels=top_5_features,
                        annot=True,
                        fmt='.4f',
                        cmap='YlOrRd',
                        square=True,
                        linewidths=1,
                        cbar_kws={"label": "Mean |SHAP Interaction|"},
                        vmin=0,
                        mask=mask
                    )
                    
                    plt.title(f"{self.percolator_model.upper()} - Feature Interaction Matrix (Top 5 Features)", 
                            fontsize=14, fontweight='bold')
                    plt.xlabel("Feature", fontsize=11, fontweight='bold')
                    plt.ylabel("Feature", fontsize=11, fontweight='bold')
                    plt.xticks(rotation=45, ha='right', fontsize=10)
                    plt.yticks(rotation=0, fontsize=10)
                    
                    plt.text(0.5, -0.15, "Diagonal = Main effects | Off-diagonal = Interaction effects", 
                            ha='center', va='top', transform=plt.gca().transAxes, 
                            fontsize=9, style='italic', color='gray')
                    
                    plt.tight_layout()
                    plt.savefig(
                        self.outpath / f"{self.percolator_model}_shap_interaction_top5_matrix.png", 
                        dpi=300, bbox_inches='tight'
                    )
                    plt.savefig(
                        self.outpath / f"{self.percolator_model}_shap_interaction_top5_matrix.pdf", 
                        dpi=300, bbox_inches='tight'
                    )
                    plt.close()
                    logging.info(f"[SHAP] ✓ Top 5 interaction matrix saved")
                    
                    # ========================================================
                    # 方法3: Top 5特征的交互摘要图（分面显示）
                    # ========================================================
                    logging.info("[SHAP] Drawing interaction summary for top 5 features...")
                    
                    # 只显示Top 5特征的交互
                    shap_interaction_top5 = shap_interaction_values[:, top_5_indices, :][:, :, top_5_indices]
                    X_sample_top5 = X_sample[:, top_5_indices]
                    
                    plt.figure(figsize=(12, 10))
                    shap.summary_plot(
                        shap_interaction_top5,
                        X_sample_top5,
                        feature_names=top_5_features,
                        show=False
                    )
                    plt.title(f"{self.percolator_model.upper()} - SHAP Interaction Summary (Top 5 Features)", 
                            fontsize=14, fontweight='bold')
                    plt.tight_layout()
                    plt.savefig(
                        self.outpath / f"{self.percolator_model}_shap_interaction_top5_summary.png", 
                        dpi=300, bbox_inches='tight'
                    )
                    plt.savefig(
                        self.outpath / f"{self.percolator_model}_shap_interaction_top5_summary.pdf", 
                        dpi=300, bbox_inches='tight'
                    )
                    plt.close()
                    logging.info(f"[SHAP] ✓ Top 5 interaction summary saved")
                    
                    # ========================================================
                    # 方法4: 绘制最强的交互对（从Top 5中选择）
                    # ========================================================
                    logging.info("[SHAP] Identifying strongest interactions among top 5 features...")
                    
                    # 从Top 5中找到最强的交互（排除对角线）
                    interaction_matrix_no_diag = interaction_matrix_top5.copy()
                    np.fill_diagonal(interaction_matrix_no_diag, 0)
                    
                    # 找到最强的3对交互
                    flat_indices = np.argsort(interaction_matrix_no_diag.ravel())[::-1]
                    top_interaction_pairs = []
                    seen_pairs = set()
                    
                    for flat_idx in flat_indices:
                        i, j = np.unravel_index(flat_idx, interaction_matrix_no_diag.shape)
                        pair = tuple(sorted([i, j]))
                        if pair not in seen_pairs and i != j and interaction_matrix_no_diag[i, j] > 0:
                            top_interaction_pairs.append((i, j, interaction_matrix_no_diag[i, j]))
                            seen_pairs.add(pair)
                            if len(top_interaction_pairs) >= 3:
                                break
                    
                    logging.info(f"[SHAP] Top 3 interaction pairs among top 5 features:")
                    for rank, (i, j, strength) in enumerate(top_interaction_pairs):
                        logging.info(f"[SHAP]   {rank+1}. {top_5_features[i]} × {top_5_features[j]} (strength: {strength:.4f})")
                    
                    # 为每对交互绘制详细的依赖图
                    for rank, (i, j, strength) in enumerate(top_interaction_pairs):
                        feat_i_global = top_5_indices[i]
                        feat_j_global = top_5_indices[j]
                        
                        plt.figure(figsize=(10, 7))
                        shap.dependence_plot(
                            (feat_i_global, feat_j_global),
                            shap_interaction_values,
                            X_sample,
                            feature_names=self.feature_list,
                            show=False
                        )
                        plt.title(
                            f"{self.percolator_model.upper()} - Interaction #{rank+1}: "
                            f"{top_5_features[i]} × {top_5_features[j]}\n"
                            f"(Interaction Strength: {strength:.4f})", 
                            fontsize=12, fontweight='bold'
                        )
                        plt.tight_layout()
                        
                        safe_name_i = top_5_features[i].replace('/', '_').replace(' ', '_')
                        safe_name_j = top_5_features[j].replace('/', '_').replace(' ', '_')
                        plt.savefig(
                            self.outpath / f"{self.percolator_model}_shap_interaction_pair_{rank+1}_{safe_name_i}_x_{safe_name_j}.png", 
                            dpi=300, bbox_inches='tight'
                        )
                        plt.close()
                    
                    logging.info(f"[SHAP] ✓ Top 3 interaction pair plots saved")
                    
                    # ========================================================
                    # 方法5: 保存交互强度数值表
                    # ========================================================
                    interaction_df = pd.DataFrame(
                        interaction_matrix_top5,
                        index=top_5_features,
                        columns=top_5_features
                    )
                    interaction_df.to_csv(
                        self.outpath / f"{self.percolator_model}_shap_interaction_top5_matrix.csv"
                    )
                    logging.info(f"[SHAP] ✓ Interaction matrix saved to CSV")
                    
                else:
                    # ========================================================
                    # 备选方案: 使用SHAP值相关性近似交互（当无法计算交互值时）
                    # ========================================================
                    logging.warning("[SHAP] Using SHAP correlation approximation for interaction analysis")
                    
                    interaction_matrix_approx = np.zeros((5, 5))
                    
                    for i, idx_i in enumerate(top_5_indices):
                        for j, idx_j in enumerate(top_5_indices):
                            if i == j:
                                # 对角线：主效应（平均绝对SHAP值）
                                interaction_matrix_approx[i, j] = np.abs(shap_values_for_plot[:, idx_i]).mean()
                            else:
                                # 非对角线：交互（SHAP值的相关性）
                                correlation = np.corrcoef(
                                    shap_values_for_plot[:, idx_i],
                                    shap_values_for_plot[:, idx_j]
                                )[0, 1]
                                interaction_matrix_approx[i, j] = abs(correlation) * np.sqrt(
                                    np.abs(shap_values_for_plot[:, idx_i]).mean() * 
                                    np.abs(shap_values_for_plot[:, idx_j]).mean()
                                )
                    
                    plt.figure(figsize=(10, 8))
                    import seaborn as sns
                    
                    sns.heatmap(
                        interaction_matrix_approx,
                        xticklabels=top_5_features,
                        yticklabels=top_5_features,
                        annot=True,
                        fmt='.4f',
                        cmap='YlOrRd',
                        square=True,
                        linewidths=1,
                        cbar_kws={"label": "Interaction Strength (Approximation)"}
                    )
                    
                    plt.title(f"{self.percolator_model.upper()} - Feature Interaction (Correlation-based Approximation)", 
                            fontsize=13, fontweight='bold')
                    plt.xlabel("Feature", fontsize=11, fontweight='bold')
                    plt.ylabel("Feature", fontsize=11, fontweight='bold')
                    plt.xticks(rotation=45, ha='right', fontsize=10)
                    plt.yticks(rotation=0, fontsize=10)
                    
                    plt.text(0.5, -0.12, "Note: Approximation based on SHAP value correlations", 
                            ha='center', va='top', transform=plt.gca().transAxes, 
                            fontsize=9, style='italic', color='red')
                    
                    plt.tight_layout()
                    plt.savefig(
                        self.outpath / f"{self.percolator_model}_shap_interaction_top5_approximation.png", 
                        dpi=300, bbox_inches='tight'
                    )
                    plt.close()
                    
                    # 保存近似矩阵
                    interaction_df_approx = pd.DataFrame(
                        interaction_matrix_approx,
                        index=top_5_features,
                        columns=top_5_features
                    )
                    interaction_df_approx.to_csv(
                        self.outpath / f"{self.percolator_model}_shap_interaction_top5_approximation.csv"
                    )
                    
                    logging.info(f"[SHAP] ✓ Approximated interaction matrix saved")
                
                logging.info(f"[SHAP] ========== Interaction Analysis Complete ==========")
                
            except Exception as e:
                logging.error(f"[SHAP] Failed in interaction analysis: {str(e)}")
                import traceback
                logging.debug(traceback.format_exc())
            # ============================================================
            # 5. SHAP Heatmap (热图) - 展示多个样本的SHAP值模式
            # ============================================================
            logging.info("[SHAP] Generating SHAP heatmap...")
            
            # 选择最多50个样本进行可视化
            n_samples_heatmap = min(50, len(explain_df))
            
            # 按ml_score排序，选择高分和低分样本
            sorted_indices = explain_df['ml_score'].argsort()
            top_indices = sorted_indices[-n_samples_heatmap//2:]  # 高分样本
            bottom_indices = sorted_indices[:n_samples_heatmap//2]  # 低分样本
            selected_indices = np.concatenate([top_indices, bottom_indices])
            
            plt.figure(figsize=(14, 10))
            shap.plots.heatmap(
                shap.Explanation(
                    values=shap_values_for_plot[selected_indices],
                    base_values=base_value,
                    data=X_display[selected_indices],
                    feature_names=self.feature_list
                ),
                show=False,
                max_display=15  # 最多显示15个特征
            )
            plt.title(f"{self.percolator_model.upper()} - SHAP Heatmap (Top/Bottom Samples)", 
                    fontsize=14, fontweight='bold')
            plt.tight_layout()
            plt.savefig(
                self.outpath / f"{self.percolator_model}_shap_heatmap.png", 
                dpi=300, bbox_inches='tight'
            )
            plt.savefig(
                self.outpath / f"{self.percolator_model}_shap_heatmap.pdf", 
                dpi=300, bbox_inches='tight'
            )
            plt.close()
            logging.info(f"[SHAP] ✓ Heatmap saved")
            
            # ============================================================
            # 6. Waterfall Plot (瀑布图) - 展示单个预测的详细分解
            # ============================================================
            logging.info("[SHAP] Generating waterfall plots...")
            
            # 为代表性样本生成瀑布图
            representative_samples = [
                ("high_target", target_samples.nlargest(1, 'ml_score').index[0] if len(target_samples) > 0 else None),
                ("high_decoy", decoy_samples.nlargest(1, 'ml_score').index[0] if len(decoy_samples) > 0 else None),
            ]
            
            for label, idx in representative_samples:
                if idx is None:
                    continue
                    
                sample_idx_in_explain = explain_df.index.get_loc(idx)
                
                plt.figure(figsize=(10, 8))
                shap.plots.waterfall(
                    shap.Explanation(
                        values=shap_values_for_plot[sample_idx_in_explain],
                        base_values=base_value,
                        data=X_display[sample_idx_in_explain],
                        feature_names=self.feature_list
                    ),
                    show=False,
                    max_display=15
                )
                plt.title(f"{self.percolator_model.upper()} - Waterfall Plot: {label.replace('_', ' ').title()}", 
                        fontsize=12, fontweight='bold')
                plt.tight_layout()
                plt.savefig(
                    self.outpath / f"{self.percolator_model}_shap_waterfall_{label}.png", 
                    dpi=300, bbox_inches='tight'
                )
                plt.close()
            
            logging.info(f"[SHAP] ✓ Waterfall plots saved")
            
            # ============================================================
            # 7. Decision Plot (决策图) - 展示多个样本的决策路径
            # ============================================================
            logging.info("[SHAP] Generating decision plot...")
            
            # 选择一部分样本（target和decoy各一半）
            n_samples_decision = min(20, len(explain_df))
            target_decision = target_samples.sample(min(n_samples_decision//2, len(target_samples)), random_state=42).index if len(target_samples) > 0 else []
            decoy_decision = decoy_samples.sample(min(n_samples_decision//2, len(decoy_samples)), random_state=42).index if len(decoy_samples) > 0 else []
            
            decision_indices = []
            for idx in list(target_decision) + list(decoy_decision):
                decision_indices.append(explain_df.index.get_loc(idx))
            
            if len(decision_indices) > 0:
                plt.figure(figsize=(10, 8))
                shap.decision_plot(
                    base_value,
                    shap_values_for_plot[decision_indices],
                    X_display[decision_indices],
                    feature_names=self.feature_list,
                    show=False,
                    highlight=[i for i in range(len(target_decision))]  # 高亮target样本
                )
                plt.title(f"{self.percolator_model.upper()} - Decision Plot (Blue=Target, Gray=Decoy)", 
                        fontsize=12, fontweight='bold')
                plt.tight_layout()
                plt.savefig(
                    self.outpath / f"{self.percolator_model}_shap_decision.png", 
                    dpi=300, bbox_inches='tight'
                )
                plt.close()
                logging.info(f"[SHAP] ✓ Decision plot saved")
            
            # ============================================================
            # 8. Scatter Plot (散点图) - 特征值vs SHAP值
            # ============================================================
            logging.info("[SHAP] Generating scatter plots for top features...")
            
            for idx in top_features_idx[:3]:  # 前3个最重要特征
                feature_name = self.feature_list[idx]
                
                plt.figure(figsize=(10, 6))
                shap.plots.scatter(
                    shap.Explanation(
                        values=shap_values_for_plot[:, idx],
                        base_values=base_value,
                        data=X_display[:, idx],
                        feature_names=[feature_name]
                    ),
                    show=False
                )
                plt.title(f"{self.percolator_model.upper()} - Scatter Plot: {feature_name}", 
                        fontsize=12, fontweight='bold')
                plt.tight_layout()
                
                safe_feature_name = feature_name.replace('/', '_').replace(' ', '_')
                plt.savefig(
                    self.outpath / f"{self.percolator_model}_shap_scatter_{safe_feature_name}.png", 
                    dpi=300, bbox_inches='tight'
                )
                plt.close()
            
            logging.info(f"[SHAP] ✓ Scatter plots saved")
            
            # ============================================================
            # 9. Beeswarm Plot (蜂群图) - 更现代的summary plot
            # ============================================================
            logging.info("[SHAP] Generating beeswarm plot...")
            
            plt.figure(figsize=(12, 8))
            shap.plots.beeswarm(
                shap.Explanation(
                    values=shap_values_for_plot,
                    base_values=base_value,
                    data=X_display,
                    feature_names=self.feature_list
                ),
                show=False,
                max_display=15
            )
            plt.title(f"{self.percolator_model.upper()} - Beeswarm Plot", 
                    fontsize=14, fontweight='bold')
            plt.tight_layout()
            plt.savefig(
                self.outpath / f"{self.percolator_model}_shap_beeswarm.png", 
                dpi=300, bbox_inches='tight'
            )
            plt.savefig(
                self.outpath / f"{self.percolator_model}_shap_beeswarm.pdf", 
                dpi=300, bbox_inches='tight'
            )
            plt.close()
            logging.info(f"[SHAP] ✓ Beeswarm plot saved")
            
            # ============================================================
            # 10. 生成SHAP值统计摘要
            # ============================================================
            logging.info("[SHAP] Generating summary statistics...")
            
            summary_stats = {
                'Feature': self.feature_list,
                'Mean_Abs_SHAP': np.abs(shap_values_for_plot).mean(axis=0),
                'Mean_SHAP': shap_values_for_plot.mean(axis=0),
                'Std_SHAP': shap_values_for_plot.std(axis=0),
                'Max_SHAP': shap_values_for_plot.max(axis=0),
                'Min_SHAP': shap_values_for_plot.min(axis=0)
            }
            
            summary_df = pd.DataFrame(summary_stats)
            summary_df = summary_df.sort_values('Mean_Abs_SHAP', ascending=False)
            summary_df.to_csv(
                self.outpath / f"{self.percolator_model}_shap_summary_stats.csv", 
                index=False
            )
            logging.info(f"[SHAP] ✓ Summary statistics saved to CSV")
            
            # ============================================================
            # 打印总结信息
            # ============================================================
            logging.info(f"\n[SHAP] ========== Visualization Summary ==========")
            logging.info(f"[SHAP] Total plots generated: 10+ types")
            logging.info(f"[SHAP] 1. Feature Importance (Bar)")
            logging.info(f"[SHAP] 2. Summary Plot (Dot)")
            logging.info(f"[SHAP] 3. Dependence Plots (Top 5 features)")
            logging.info(f"[SHAP] 4. Force Plots (3 representative samples)")
            logging.info(f"[SHAP] 5. Heatmap (Top/Bottom samples)")
            logging.info(f"[SHAP] 6. Waterfall Plots (2 representative samples)")
            logging.info(f"[SHAP] 7. Decision Plot (Target vs Decoy)")
            logging.info(f"[SHAP] 8. Scatter Plots (Top 3 features)")
            logging.info(f"[SHAP] 9. Beeswarm Plot (Modern summary)")
            logging.info(f"[SHAP] 10. Summary Statistics (CSV)")
            logging.info(f"[SHAP] ===============================================\n")
                    
        except Exception as e:
            logging.error(f"[SHAP] Failed in visualization: {str(e)}")
            import traceback
            logging.debug(traceback.format_exc())

    def re_score(self, df: pd.DataFrame) -> pd.DataFrame:
        """Rescoring using FDR strategy"""
        logging.info(
            "[Rescore] "
            f"{np.sum((df.fdr<=self.fdr) & (df.decoy==0))} "
            f"target PSMs at {self.fdr} precursor-level FDR"
        )
        for i in range(self.iter_num):
            logging.info(f"[Rescore] Iteration {i+1} of Rescore ...")
            df = self._cv_score(df)
            df = self._estimate_fdr(df, self.fdr_level, False)

            logging.info(
                f"[Rescore] {len(df[(df.fdr<=self.fdr) & (df.decoy==0)])} "
                f"target PSMs at {self.fdr} precursor-level FDR"
            )
        df = self._estimate_fdr(df)
        dr = df[(df.fdr<=self.fdr) & (df.decoy==0)]
        dr1 = dr.drop_duplicates(subset=['precursor_id'])

        logging.info(
            f"[Rescore] Finished with {dr1.shape[0]} "
            f"target Precursors at {self.fdr} {self.fdr_level}-level FDR"
        )

        if self.shap_explain:
            try:
                self.explain_model_with_shap(df)
            except Exception as e:
                logging.error(f"[SHAP] failed: {str(e)}")
        return df

    def re_score_topk(self, df: pd.DataFrame, top_ratio: float = 0.3) -> pd.DataFrame:
        """Rescoring using Top-K strategy"""
        logging.info(
            "[Rescore] "
            f"{np.sum((df.fdr<=self.fdr) & (df.decoy==0))} "
            f"target PSMs at {self.fdr} precursor-level FDR"
        )
        
        for i in range(self.iter_num):
            logging.info(f"[Rescore] Iteration {i+1} of Rescore (Top-{top_ratio*100}%) ...")
            df = self._cv_score_topk(df, top_ratio)
            df = self._estimate_fdr(df, self.fdr_level, False)

            logging.info(
                f"[Rescore] {len(df[(df.fdr<=self.fdr) & (df.decoy==0)])} "
                f"target PSMs at {self.fdr} precursor-level FDR"
            )
        
        df = self._estimate_fdr(df)
        dr = df[(df.fdr<=self.fdr) & (df.decoy==0)]
        dr1 = dr.drop_duplicates(subset=['precursor_id'])

        logging.info(
            f"[Rescore] Finished with {dr1.shape[0]} "
            f"target Precursors at {self.fdr} {self.fdr_level}-level FDR"
        )

        if self.shap_explain:
            try:
                self.explain_model_with_shap(df)
            except Exception as e:
                logging.error(f"[SHAP] failed: {str(e)}")
        
        return df