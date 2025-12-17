import sqlite3
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import joblib
import os


def load_bioactivity_data(db_path):
    query = """
    SELECT DISTINCT
        md.molregno,
        cs.canonical_smiles,
        act.standard_value,
        act.standard_type,
        td.pref_name as target_name,
        td.chembl_id as target_chembl_id
    FROM molecule_dictionary md
    INNER JOIN compound_structures cs ON md.molregno = cs.molregno
    INNER JOIN activities act ON md.molregno = act.molregno
    INNER JOIN assays a ON act.assay_id = a.assay_id
    INNER JOIN target_dictionary td ON a.tid = td.tid
    WHERE act.standard_type = 'IC50'
        AND act.standard_units = 'nM'
        AND act.standard_value IS NOT NULL
        AND act.standard_value > 0
        AND a.assay_type = 'B'
        AND cs.canonical_smiles IS NOT NULL
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    print(f"Loaded {len(df)} bioactivity records from database")
    return df


def clean_bioactivity_data(df):
    initial_count = len(df)
    
    df = df.dropna(subset=['canonical_smiles', 'standard_value'])
    
    valid_smiles_mask = df['canonical_smiles'].apply(lambda x: Chem.MolFromSmiles(x) is not None)
    df = df[valid_smiles_mask].copy()
    print(f"After removing invalid SMILES: {len(df)} records")
    
    df = df.drop_duplicates(subset=['canonical_smiles'], keep='first')
    print(f"After removing duplicates: {len(df)} records")
    
    df['pIC50'] = -np.log10(df['standard_value'] * 1e-9)
    
    df = df[(df['pIC50'] >= 0) & (df['pIC50'] <= 15)].copy()
    print(f"After filtering pIC50 range [0, 15]: {len(df)} records")
    
    print(f"Cleaned dataset: {len(df)} records (removed {initial_count - len(df)} records)")
    return df


def compute_morgan_fingerprints(smiles_list, radius=2, n_bits=2048):
    fingerprints = []
    valid_indices = []
    
    for idx, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            fingerprints.append(np.array(fp))
            valid_indices.append(idx)
    
    X = np.array(fingerprints)
    print(f"Generated {len(fingerprints)} Morgan fingerprints (radius={radius}, nBits={n_bits})")
    return X, valid_indices


def train_bioactivity_model(X_train, y_train, model_type='random_forest'):
    print(f"Training {model_type} model...")
    print(f"Training set size: {len(X_train)} samples")
    
    if model_type == 'random_forest':
        model = RandomForestRegressor(
            n_estimators=100,
            random_state=42,
            n_jobs=-1,
            max_depth=20,
            min_samples_split=5,
            min_samples_leaf=2
        )
    elif model_type == 'gradient_boosting':
        model = GradientBoostingRegressor(
            n_estimators=100,
            random_state=42,
            max_depth=10,
            min_samples_split=5,
            learning_rate=0.1
        )
    
    model.fit(X_train, y_train)
    print("Model training completed")
    return model


def evaluate_model(model, X_test, y_test):
    print(f"Evaluating model on {len(X_test)} test samples...")
    y_pred = model.predict(X_test)
    
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    mae = mean_absolute_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    
    print("=" * 50)
    print("Bioactivity Model Evaluation Results")
    print("=" * 50)
    print(f"RMSE: {rmse:.4f}")
    print(f"MAE: {mae:.4f}")
    print(f"R² Score: {r2:.4f}")
    print("=" * 50)
    
    return {'rmse': rmse, 'mae': mae, 'r2': r2, 'y_pred': y_pred, 'y_true': y_test}


def save_bioactivity_model(model, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    model_path = os.path.join(output_dir, 'bioactivity_model.joblib')
    joblib.dump(model, model_path)
    print(f"Model saved to {model_path}")
    
    config = {
        'model_type': 'pIC50_regressor',
        'fingerprint_type': 'Morgan',
        'fingerprint_radius': 2,
        'fingerprint_bits': 2048,
        'target_column': 'pIC50'
    }
    config_path = os.path.join(output_dir, 'bioactivity_config.joblib')
    joblib.dump(config, config_path)
    print(f"Config saved to {config_path}")
    
    return model_path, config_path


def main():
    print("=" * 60)
    print("Bioactivity Model Training Pipeline")
    print("=" * 60)
    
    db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"
    output_dir = "trained_models"
    
    print("\n[Step 1/5] Loading bioactivity data...")
    df = load_bioactivity_data(db_path)
    
    print("\n[Step 2/5] Cleaning dataset...")
    df_clean = clean_bioactivity_data(df)
    
    print("\n[Step 3/5] Computing molecular fingerprints...")
    smiles_list = df_clean['canonical_smiles'].tolist()
    X, valid_indices = compute_morgan_fingerprints(smiles_list)
    y = df_clean['pIC50'].iloc[valid_indices].values
    
    print(f"Final dataset: {len(X)} molecules")
    
    print("\n[Step 4/5] Splitting and training...")
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    print(f"Training: {len(X_train)}, Test: {len(X_test)}")
    
    model = train_bioactivity_model(X_train, y_train)
    
    print("\n[Step 5/5] Evaluating and saving...")
    results = evaluate_model(model, X_test, y_test)
    save_bioactivity_model(model, output_dir)
    
    print("\nBioactivity model training completed!")
    return model, results


if __name__ == "__main__":
    model, results = main()
