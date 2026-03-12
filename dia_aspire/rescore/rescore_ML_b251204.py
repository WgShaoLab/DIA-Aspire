import numpy as np
import pandas as pd
import torch
import os
import sys
import multiprocessing as mp
from tqdm import tqdm
import shap
import matplotlib.pyplot as plt

from rescore.fdr import fdr_from_ref, fdr_to_q_values
from peptdeep.utils import logging

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
        
        # 全局平均池化
        self.global_avg_pool = torch.nn.AdaptiveAvgPool1d(1)
        
        # 全连接层
        self.fc = torch.nn.Linear(hidden_channels*4, 1)
        
        # Dropout
        self.dropout = torch.nn.Dropout(dropout)
        
        # 激活函数
        self.relu = torch.nn.ReLU()
    
    def forward(self, x):

        if x.dim() == 2:
        # [B, D] -> [B, 1, D]
            x = x.unsqueeze(1)

        elif x.dim() == 3:
            # [B, C, D] — 如果 C != 1，可能是错误；但我们假设 C=1
            if x.shape[1] != 1:
                # 如果是 [B, D, 1]，转置并 squeeze
                if x.shape[2] == 1:
                    x = x.squeeze(-1).unsqueeze(1)  # [B, D, 1] -> [B, D] -> [B, 1, D]
                else:
                    raise ValueError(f"Expected channel dim=1, got {x.shape}")
                
        elif x.dim() == 4:
            # 处理 SHAP 可能传入的 [B, 1, 1, D] 或 [B, 1, D, 1]
            # 先 squeeze 所有长度为1的维度（除了 batch 和 feature）
            x = x.squeeze()
            if x.dim() == 2:
                x = x.unsqueeze(1)  # [B, D] -> [B, 1, D]
            elif x.dim() == 3:
                # 可能是 [B, 1, D] 或 [B, D, 1]
                if x.shape[1] == 1:
                    pass  # already [B, 1, D]
                elif x.shape[2] == 1:
                    x = x.squeeze(-1).unsqueeze(1)  # [B, D, 1] -> [B, 1, D]
                else:
                    raise ValueError(f"Cannot handle 3D input with shape {x.shape}")
            else:
                raise ValueError(f"Unexpected shape after squeeze: {x.shape}")
        else:
            raise ValueError(f"Unsupported input dimension: {x.dim()}")

        # 现在 x 一定是 [B, 1, D]
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
        
        # 全局平均池化
        x = self.global_avg_pool(x)
        x = x.squeeze(-1)  # 移除长度维度

        # 全连接层
        x = self.fc(x)
        return x    #.squeeze(1)

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
        return self.model(x)   #.squeeze(1)

class NNLinear(torch.nn.Module):
    """Torch-based NN_Linear rescore model"""

    def __init__(self, input_dim, **kwargs):
        super().__init__()
        torch.manual_seed(1290)
        self.linear = torch.nn.Linear(input_dim, 1)

    def forward(self, x):
        return self.linear(x)     #.squeeze(1)

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

class NNLinearRescore:
    def __init__(self, num_features, nn_model_type="nnlinear", **model_kwargs):
        self.nn_model = rescore_model_provider.get_model(nn_model_type, num_features, **model_kwargs)
        self.train_batch_size = 10000
        self.predict_batch_size = 100000
        self.nn_model_type = nn_model_type.lower() 

        self.optimizer = torch.optim.Adam(
            self.nn_model.parameters(), lr=0.001,weight_decay=1e-4
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

class NNRescore:
    def __init__(self, num_features, nn_model_type="dnn", **model_kwargs):
        
        self.nn_model = rescore_model_provider.get_model(nn_model_type, num_features,**model_kwargs)
        self.train_batch_size = 10000
        self.predict_batch_size = 100000
        
        self.nn_model_type = nn_model_type.lower()  # 保存模型类型
        
        # 设备配置
        if torch.cuda.is_available():
            self.device = torch.device("cuda")
            self.nn_model.to(self.device)
        else:
            self.device = torch.device("cpu")

        # 优化器配置
        self.optimizer = torch.optim.Adam(
            self.nn_model.parameters(), 
            lr=0.001,
            weight_decay=1e-4
        )
        self.loss_func = torch.nn.BCEWithLogitsLoss()
            
        # 训练参数
        self.epoch = 20
        
        # 为DNN模型添加学习率调度器和早停机制
        if self.nn_model_type in ["dnn","cnn"]:
            self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                self.optimizer, 
                mode='min',
                factor=0.5,
                patience=2,
                # verbose=True ,
            )
            self.patience = 3  # 早停耐心值
        else:
            self.scheduler = None

    def fit(self, features, labels, is_for_shap = False):
        # 将数据转换为PyTorch张量
        features_tensor = torch.tensor(features, dtype=torch.float)
        labels_tensor = torch.tensor(labels, dtype=torch.float)
        
        # 为DNN模型创建训练集和验证集
        if self.nn_model_type in ["dnn","cnn"]:
            # 分割数据集（90%训练，10%验证）
            dataset_size = len(features_tensor)
            indices = torch.randperm(dataset_size).tolist()
            train_size = int(0.9 * dataset_size)
            train_indices = indices[:train_size]
            val_indices = indices[train_size:]
            
            train_features = features_tensor[train_indices]
            train_labels = labels_tensor[train_indices]
            val_features = features_tensor[val_indices]
            val_labels = labels_tensor[val_indices]
            
            # 将数据移动到设备
            train_features = train_features.to(self.device)
            train_labels = train_labels.to(self.device)
            val_features = val_features.to(self.device)
            val_labels = val_labels.to(self.device)
            
            best_val_loss = float('inf')
            patience_counter = 0
            best_model_state = None
       
        # 训练循环
        for epoch in range(self.epoch):
            
            # DNN模型：完整训练流程
            self.nn_model.train()
            train_loss = 0.0
            num_batches = 0
            
            # 训练阶段
            for i in range(0, len(train_features), self.train_batch_size):
                self.optimizer.zero_grad()
                
                batch_features = train_features[i:i+self.train_batch_size]
                batch_labels = train_labels[i:i+self.train_batch_size].unsqueeze(-1)
                
                outputs = self.nn_model(batch_features)
                loss = self.loss_func(outputs, batch_labels)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.nn_model.parameters(), 1.0)  # 梯度裁剪
                self.optimizer.step()
                
                train_loss += loss.item()
                num_batches += 1
            
            avg_train_loss = train_loss / num_batches if num_batches > 0 else 0
            
            # 验证阶段
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
            
            # 更新学习率
            self.scheduler.step(avg_val_loss)
            
            # 打印训练状态
            logging.debug(f"[{self.nn_model_type.upper()}] Epoch {epoch+1}/{self.epoch} - "
                         f"Train Loss: {avg_train_loss:.4f}, Val Loss: {avg_val_loss:.4f}")
            
            # 早停机制
            if avg_val_loss < best_val_loss:
                best_val_loss = avg_val_loss
                patience_counter = 0
                # 保存最佳模型
                best_model_state = self.nn_model.state_dict()
            else:
                patience_counter += 1
                if patience_counter >= self.patience:
                    if is_for_shap:
                        logging.info(f"[SHAP] [{self.nn_model_type.upper()}] Early stopping at epoch {epoch+1}")
                    else:
                        logging.info(f"[Rescore] [{self.nn_model_type.upper()}] Early stopping at epoch {epoch+1}")
                    break
        
        # 恢复最佳模型
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
                
                if self.nn_model_type in ["dnn","cnn"]:
                    pred = torch.sigmoid(pred)
                
                outputs[i:i+self.predict_batch_size] = pred.squeeze(1).cpu().numpy()
        return outputs

    # for SHAP
    def predict_proba(self, X):
        """ """
        self.nn_model.eval()  # 设置为评估模式
        with torch.no_grad():
            if not isinstance(X, torch.Tensor):
                X = torch.tensor(X, dtype=torch.float32)
            X = X.to(self.device)
            preds = self.nn_model(X)
            preds = torch.sigmoid(preds)
            return preds.cpu().numpy()

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

        self.init_percolator_model(percolator_model, percolator_backend)

    def init_percolator_model(
        self, percolator_model="svm", percolator_backend="sklearn"
    ):
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.linear_model import LogisticRegression
        from sklearn.naive_bayes import GaussianNB
        from sklearn.svm import LinearSVR, LinearSVC
        from sklearn.calibration import CalibratedClassifierCV
        from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
        from sklearn.preprocessing import StandardScaler 

        self.percolator_model = percolator_model.lower()
        self.percolator_backend = percolator_backend.lower()

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
            if self.percolator_model == "cnn":
                model_kwargs = {
                    'hidden_channels': 128,
                    'kernel_size': 3,
                    'dropout': 0.2
                }
                self.model = NNRescore(
                len(self.feature_list), nn_model_type=percolator_model, **model_kwargs
                )
            
            elif self.percolator_model == "nnlinear":
                self.model = NNLinearRescore(
                len(self.feature_list), nn_model_type=percolator_model,**model_kwargs
            )
                
        elif self.percolator_model == "lr":
            # self.model = LogisticRegression(solver="liblinear",
            #                                 max_iter=1000,random_state=42)
            self.model = LogisticRegression(
                    penalty='elasticnet',   # L1+L2 混合
                    l1_ratio=0.3,           # 偏向 Ridge（更稳定）
                    C=0.1,                  # 较强正则（C 越小正则越强）
                    solver='saga',
                    max_iter=2000,
                    random_state=42
                )
        elif self.percolator_model == "bayes": 
            
            self.scaler = StandardScaler()
            self.model = GaussianNB(var_smoothing=1e-6)  
            # self.model = GaussianNB()
        
        elif self.percolator_model == "lda":
            # LDA 对数值稳定性敏感，建议标准化（虽然理论上不需要，但实践中更稳）
            self.scaler = StandardScaler()
            self.model = LinearDiscriminantAnalysis(
                solver='lsqr', 
                shrinkage= 'auto',
                )
            
        elif self.percolator_model == "rf":
            self.model = RandomForestClassifier(n_estimators=100,
                                                max_depth=4,     
                                                min_samples_split=2000,       
                                                min_samples_leaf=1000,        
                                                max_features="sqrt", 
                                                bootstrap=True,
                                                oob_score=True,       
                                                random_state=42,
                                                n_jobs=-1)
        elif self.percolator_model == "svm":
            self.model = LinearSVR(epsilon=0.1,C=1,dual=False,
                                             max_iter=5000,
                                             random_state=42,
                                             loss='squared_epsilon_insensitive')
            
            # base_svm = LinearSVC(C=0.1, random_state=42, max_iter=2000)
            # self.model = CalibratedClassifierCV(base_svm, method='sigmoid', cv=3)

        else:
            logging.info(
                "[Rescore] "
                f"Rescoring model '{percolator_model}' is not "
                "implemented, switch to sklearn 'svm' model."
            )
            self.model = LinearSVR(epsilon=0.1,C=1,dual=False,
                                             max_iter=5000,
                                             random_state=42,
                                             loss='squared_epsilon_insensitive')
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

        df = df.sort_values(by=["ml_score", "decoy"], ascending=[False,True])
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
            group_key = "Precursor.Id"
        elif fdr_level == "peptide":
            group_key = "Stripped.Sequence"
        else:
            group_key = "Precursor.Id"

        _df = df.groupby([group_key, "decoy"],as_index=False)["ml_score"].max()
        _df = _df.sort_values(["ml_score", "decoy"], ascending=[False,True])

        target_values = 1 - _df["decoy"].values
        decoy_cumsum = np.cumsum(_df["decoy"].values)
        target_cumsum = np.cumsum(target_values)

        with np.errstate(divide='ignore', invalid='ignore'):
            fdr_values = np.where(
                target_cumsum > 0, 
                decoy_cumsum / target_cumsum, 
                1 
            )
        #fdr_values = decoy_cumsum / target_cumsum
        _df["fdr"] = fdr_to_q_values(fdr_values)
        
        df["fdr"] = fdr_from_ref(
                df["ml_score"].values, _df["ml_score"].values, _df["fdr"].values
            )
        return df

    def _train(self, train_t_df, train_d_df):

        if train_t_df.empty or train_d_df.empty:
            logging.warning("[Rescore] Empty training set. Skipping training.")
            return
        if len(train_t_df) > self.max_train_sample:
            train_t_df = train_t_df.sample(n=self.max_train_sample, random_state=1337)
        if len(train_d_df) > self.max_train_sample:
            train_d_df = train_d_df.sample(n=self.max_train_sample, random_state=1337)

        train_df = pd.concat((train_t_df, train_d_df))
        train_label = np.ones(len(train_df), dtype=np.int32)
        train_label[len(train_t_df) :] = 0

        # self.model.fit(train_df[self.feature_list].values, train_label)

        X_train = train_df[self.feature_list].values

        # 标准化（仅对 LDA 和 Bayes）
        if self.percolator_model in ["lda", "bayes"]:
            X_train = self.scaler.fit_transform(X_train)
        self.model.fit(X_train, train_label)

        if self.percolator_model == "rf":
            if hasattr(self.model, 'oob_score_') and self.model.oob_score_:
                logging.info(f"[RF] OOB Score: {self.model.oob_score_:.4f}")

    def _predict(self, test_df):

        X_test = test_df[self.feature_list].values

        if self.percolator_backend == "pytorch":
            test_df["ml_score"] = self.model.decision_function(X_test)

        elif self.percolator_model == "lr":
            test_df["ml_score"] = self.model.decision_function(X_test)

        elif self.percolator_model == "bayes":
            if hasattr(self, 'scaler'):
                X_test = self.scaler.transform(X_test)
            test_df["ml_score"] = self.model.predict_proba(X_test)[:,1]
        
        elif self.percolator_model == "lda":
            if hasattr(self, 'scaler'):
                X_test = self.scaler.transform(X_test)
            test_df["ml_score"] = self.model.predict_proba(X_test)[:,1]

        elif self.percolator_model == "svm":
            test_df["ml_score"] = self.model.predict(X_test)

        elif self.percolator_model == "rf":
            test_df["ml_score"] = self.model.predict_proba(X_test)[:, 1]
        
        return test_df

    def _cv_score(self, df: pd.DataFrame) -> pd.DataFrame:

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
    
    def explain_model_with_shap(self, df: pd.DataFrame):
        """
        Parameters
        ----------
        df : pd.DataFrame
        """
        if not self.shap_explain:
            logging.info("[SHAP] is not used")
            return
        
        try:    
        # 准备训练数据
            train_t_df = df[(df.fdr <= self.fdr) & (df.decoy == 0)]
            train_d_df = df[df.decoy != 0]
            
            if len(train_t_df) < self.min_train_sample or len(train_d_df) < self.min_train_sample:
                logging.warning(f"[SHAP] sample is not sufficient (target={len(train_t_df)}, decoy={len(train_d_df)}), skipping SHAP")
                return
                
            # 平衡采样
            n_sample = min(len(train_t_df), len(train_d_df), self.max_train_sample)
            train_t_df = train_t_df.sample(n=n_sample, random_state=42)
            train_d_df = train_d_df.sample(n=n_sample, random_state=42)
            
            train_df = pd.concat([train_t_df, train_d_df])
            train_labels = np.concatenate([
                np.ones(len(train_t_df), dtype=int),
                np.zeros(len(train_d_df), dtype=int)
            ])
            
            # 训练一个用于解释的最终模型
            logging.info("[SHAP] training final model for explanation...")
            if self.percolator_model in ["cnn",'dnn']:
                self.model.fit(train_df[self.feature_list].values, train_labels, is_for_shap=True)
            else:
                self.model.fit(train_df[self.feature_list].values, train_labels)
            # 准备解释样本 - 使用目标PSM和decoy的混合
            explain_df = pd.concat([
                train_t_df.sample(min(500, len(train_t_df))), 
                train_d_df.sample(min(500, len(train_d_df)))
            ])
            shap_values = None

            try:
                logging.info("[SHAP] calculating values to explain the model ...")
            
                # 保存当前输出设置
                import contextlib
                original_stdout = sys.stdout
                original_stderr = sys.stderr
                sys.stdout = open(os.devnull, 'w')
                sys.stderr = open(os.devnull, 'w')
                
                # 临时禁用SHAP库的日志
                shap_logger = logging.getLogger('shap')
                original_shap_level = shap_logger.level
                shap_logger.setLevel(logging.CRITICAL)
                
                # 临时禁用NumPy的警告
                np_logger = logging.getLogger('numpy')
                original_numpy_level = np_logger.level
                np_logger.setLevel(logging.CRITICAL)

                # 设置环境变量避免进一步的问题
                os.environ['PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION'] = 'python'
                os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'  # 禁用oneDNN警告
                os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'   # 禁用TensorFlow日志

                # 使用上下文管理器确保设置恢复
                with contextlib.redirect_stdout(open(os.devnull, 'w')), \
                    contextlib.redirect_stderr(open(os.devnull, 'w')):

                    if self.percolator_backend == "pytorch":

                        logging.info("[SHAP] using DeepExplainer to explain PyTorch-based model")

                        pytorch_model = self.model.nn_model
                        background_data = train_df[self.feature_list].values[:100]
                        explain_data = explain_df[self.feature_list].values
                        
                        # 转换为tensor并移动到设备
                        background_tensor = torch.tensor(background_data, dtype=torch.float32).to(self.model.device)
                        explain_tensor = torch.tensor(explain_data, dtype=torch.float32).to(self.model.device)

                        if self.model.nn_model_type == "cnn":
                            background_tensor = background_tensor.unsqueeze(1)
                            explain_tensor = explain_tensor.unsqueeze(1)
        
                        explainer = shap.DeepExplainer(pytorch_model, background_tensor)
                        shap_values = explainer.shap_values(explain_tensor)
                        
                        if isinstance(shap_values, list):
                            shap_values = shap_values[0]  # 取类别0的SHAP值（或类别1，根据模型输出）
                        
                        # 移除通道维度（CNN）和最后一维（输出维度）
                        if self.model.nn_model_type == "cnn":
                            shap_values = shap_values.squeeze(1)  # [batch, input_dim, 1] → [batch, input_dim]
                        shap_values = shap_values.squeeze(-1)  # 
                        
                                               
                    elif self.percolator_model == "rf":
                        logging.info("[SHAP] using TreeExplainer to explain tree-based model")
                        explainer = shap.TreeExplainer(self.model)
                        shap_values = explainer.shap_values(explain_df[self.feature_list].values)
                        
                    elif self.percolator_model == "lr":
                        logging.info("[SHAP] using LinearExplainer to explain linear model")
                        explainer = shap.LinearExplainer(self.model, train_df[self.feature_list].values)
                        shap_values = explainer.shap_values(explain_df[self.feature_list].values)

                    elif self.percolator_model == "lda":
                        logging.info("[SHAP] using LinearExplainer to explain LDA model")
                        background = self.scaler.transform(train_df[self.feature_list].values)
                        explain_data = self.scaler.transform(explain_df[self.feature_list].values)

                        explainer = shap.LinearExplainer(self.model, background,
                            feature_perturbation="interventional"  # 或 "correlation_dependent"
                        )
                        shap_values = explainer.shap_values(explain_data)
                        
                    else:
                        logging.info(f"[SHAP] using KernelExplainer to explain {self.percolator_model} model")
                        explainer = shap.KernelExplainer(
                            self.model.predict_proba if hasattr(self.model, 'predict_proba') else self.model.predict,
                            shap.sample(train_df[self.feature_list].values, 100)
                        )
                        shap_values = explainer.shap_values(explain_df[self.feature_list].values)
                
            except Exception as e:
                logging.error(f"[SHAP] Failed: {str(e)}")
                import traceback
                logging.debug(traceback.format_exc())
                return
            
            finally:
                sys.stdout = original_stdout
                sys.stderr = original_stderr
                shap_logger.setLevel(original_shap_level)
                np_logger.setLevel(original_numpy_level)

            self._visualize_shap(shap_values, explain_df)

        except Exception as e:
            logging.error(f"[SHAP] error during processing: {str(e)}")

    import warnings
    warnings.filterwarnings("ignore", message="The figure layout has changed to tight")

    def _visualize_shap(self, shap_values, explain_df):
        """visualizing SHAP"""

        try:
            if isinstance(shap_values, list):
                if len(shap_values) == 1:
                    # PyTorch DeepExplainer 单输出情况: [array]
                    shap_values_for_plot = shap_values[0]
                elif len(shap_values) == 2:
                    # 二分类，取正类（target class）
                    shap_values_for_plot = shap_values[1]
                else:
                    # 多分类，取第一个类别（或根据需求调整）
                    shap_values_for_plot = shap_values[0]
            else:
                # sklearn 等返回的直接是 array
                shap_values_for_plot = shap_values

            # 确保是 NumPy array
            if torch.is_tensor(shap_values_for_plot):
                shap_values_for_plot = shap_values_for_plot.detach().cpu().numpy()

            # 现在 shap_values_for_plot 应该是 [n_samples, n_features]
            assert shap_values_for_plot.shape[1] == len(self.feature_list), \
                f"SHAP values feature dim {shap_values_for_plot.shape[1]} != expected {len(self.feature_list)}"

            # feature importance
            plt.figure(figsize=(12, 6))
            shap.summary_plot(shap_values_for_plot,  explain_df[self.feature_list], plot_type="bar", show=False)
            plt.title("feature importance")
            plt.tight_layout()
            
            plt.savefig(self.percolator_model+"_shap_feature_importance.png",dpi=300)
            plt.savefig(self.percolator_model+"_shap_feature_importance.pdf",dpi=300)
            plt.close()
            logging.info(f"[SHAP] feature-importance saved to: {self.percolator_model}_shap_feature_importance.png")
            
            # SHAP values
            plt.figure(figsize=(12, 8))
            shap.summary_plot(shap_values_for_plot,  explain_df[self.feature_list], show=False)
            plt.title(" SHAP values distribution")
            plt.tight_layout()
            
            plt.savefig(self.percolator_model+"_shap_summary_plot.png",dpi=300)
            plt.savefig(self.percolator_model+"_shap_summary_plot.pdf",dpi=300)
            plt.close()
            logging.info(f"[SHAP] SHAP-value saved to: {self.percolator_model}_shap_summary_plot.png")
                
            # for feature in self.feature_list:  # 只展示前5个重要特征
            #     plt.figure(figsize=(10, 6))
            #     shap.dependence_plot(
            #         feature, 
            #         shap_values_for_plot,  
            #         explain_df[self.feature_list],
            #         interaction_index=None,
            #         show=False
            #     )
            #     plt.title(f"{feature} dependence")
            #     plt.tight_layout()
            #     plt.savefig(f"shap_dependence_{feature}.png",dpi=300)
            #     # plt.savefig(f"shap_dependence_{feature}.pdf")
            #     plt.close()
            # logging.info(f"[SHAP] feature-dependance saved to: {self.percolator_model}_shap_dependence_*.png")
                
        except Exception as e:
            logging.error(f"[SHAP] Failed in visualization: {str(e)}")

    def re_score(self, df: pd.DataFrame) -> pd.DataFrame:
        logging.info(
            "[Rescore] "
            f"{np.sum((df.fdr<=self.fdr) & (df.decoy==0))} "
            f"target PSMs at {self.fdr} precursor-level FDR"
        )
        for i in range(self.iter_num):
            logging.info(f"[Rescore] Iteration {i+1} of Rescore ...")
            df = self._cv_score(df)
            # print(df)
            df = self._estimate_fdr(df, self.fdr_level, False)
            # print(df)

            logging.info(
                f"[Rescore] {len(df[(df.fdr<=self.fdr) & (df.decoy==0)])} "
                f"target PSMs at {self.fdr} precursor-level FDR"
            )
        df = self._estimate_fdr(df)
        dr = df[(df.fdr<=self.fdr) & (df.decoy==0)]
        dr1 = dr.drop_duplicates(subset=['Precursor.Id'])

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
