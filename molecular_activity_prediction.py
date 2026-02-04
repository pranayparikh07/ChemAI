import sqlite3
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
import joblib


def load_data(db_path):
    query = """
    SELECT DISTINCT
        md.molregno,
        cs.canonical_smiles,
        act.standard_value
    FROM molecule_dictionary md
    INNER JOIN compound_structures cs ON md.molregno = cs.molregno
    INNER JOIN activities act ON md.molregno = act.molregno
    INNER JOIN assays a ON act.assay_id = a.assay_id
    WHERE act.standard_type = 'IC50'
        AND act.standard_units = 'nM'
        AND act.standard_value IS NOT NULL
        AND a.assay_type = 'B'
        AND cs.canonical_smiles IS NOT NULL
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    print(f"Loaded {len(df)} records from database")
    return df


def clean_data(df):
    initial_count = len(df)
    df = df.dropna(subset=['canonical_smiles', 'standard_value'])
    print(f"After removing nulls: {len(df)} records")
    
    valid_smiles_mask = df['canonical_smiles'].apply(lambda x: Chem.MolFromSmiles(x) is not None)
    df = df[valid_smiles_mask].copy()
    print(f"After removing invalid SMILES: {len(df)} records")
    
    df = df.drop_duplicates(subset=['canonical_smiles'], keep='first')
    print(f"After removing duplicates: {len(df)} records")
    
    df = df[df['standard_value'] > 0].copy()
    df['pIC50'] = -np.log10(df['standard_value'] * 1e-9)
    print(f"After log transformation: {len(df)} records")
    
    df = df[(df['pIC50'] >= 0) & (df['pIC50'] <= 15)].copy()
    print(f"After filtering pIC50 range [0, 15]: {len(df)} records")
    
    print(f"Cleaned dataset: {len(df)} records (removed {initial_count - len(df)} records)")
    return df


def featurize(smiles_list, radius=2, n_bits=2048):
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


def train_model(X_train, y_train, n_estimators=100, random_state=42):
    print(f"Training RandomForestRegressor with {n_estimators} estimators...")
    print(f"Training set size: {len(X_train)} samples")
    
    model = RandomForestRegressor(
        n_estimators=n_estimators,
        random_state=random_state,
        n_jobs=-1,
        max_depth=20,
        min_samples_split=5,
        min_samples_leaf=2
    )
    model.fit(X_train, y_train)
    print("Model training completed")
    return model


def evaluate(model, X_test, y_test):
    print(f"Evaluating model on {len(X_test)} test samples...")
    y_pred = model.predict(X_test)
    
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    r2 = r2_score(y_test, y_pred)
    
    print("=" * 50)
    print("Model Evaluation Results")
    print("=" * 50)
    print(f"RMSE: {rmse:.4f}")
    print(f"R² Score: {r2:.4f}")
    print("=" * 50)
    
    return {'rmse': rmse, 'r2': r2, 'y_pred': y_pred}


def save_model(model, model_path, feature_params_path):
    joblib.dump(model, model_path)
    print(f"Model saved to {model_path}")
    
    feature_params = {
        'radius': 2,
        'n_bits': 2048,
        'fingerprint_type': 'Morgan'
    }
    joblib.dump(feature_params, feature_params_path)
    print(f"Feature parameters saved to {feature_params_path}")


def main():
    print("=" * 60)
    print("Molecular Activity Prediction Pipeline")
    print("Agentic AI-Based Molecular Design and Optimization System")
    print("=" * 60)
    
    db_path = "chembl_36/chembl_36_sqlite/chembl_36.db"
    model_path = "pIC50_random_forest_model.joblib"
    feature_params_path = "feature_generator_params.joblib"
    
    print("\n[Step 1/6] Loading data from ChEMBL database...")
    df = load_data(db_path)
    
    print("\n[Step 2/6] Cleaning dataset...")
    df_clean = clean_data(df)
    
    print("\n[Step 3/6] Generating molecular fingerprints...")
    smiles_list = df_clean['canonical_smiles'].tolist()
    X, valid_indices = featurize(smiles_list)
    
    y = df_clean['pIC50'].iloc[valid_indices].values
    print(f"Final dataset size: {len(X)} molecules with features")
    
    print("\n[Step 4/6] Splitting data into train/test sets (80/20)...")
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    print(f"Training set: {len(X_train)} samples")
    print(f"Test set: {len(X_test)} samples")
    
    print("\n[Step 5/6] Training RandomForest model...")
    model = train_model(X_train, y_train)
    
    print("\n[Step 6/6] Evaluating model performance...")
    results = evaluate(model, X_test, y_test)
    
    print("\nSaving model and feature generator...")
    save_model(model, model_path, feature_params_path)
    
    print("\n" + "=" * 60)
    print("Pipeline completed successfully!")
    print("=" * 60)
    
    return model, results


if __name__ == "__main__":
    model, results = main()
