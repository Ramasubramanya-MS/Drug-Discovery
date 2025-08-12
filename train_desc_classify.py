#!/usr/bin/env python3
# train_desc_classify.py
# Usage (PowerShell/cmd):
#   python train_desc_classify.py "EGFR_with_descriptors.csv" --topN 50 --thresh 0.8 --outdir .
# Requires: pandas numpy scikit-learn

import argparse, os, numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score

DESC_CANDIDATES = [
    "MolWt","LogP","TPSA","HBD","HBA","RB","FractionCSP3",
    "HeavyAtomCount","NumHeteroatoms","NumRings",
    "NumAromaticRings","NumAliphaticRings","NumSaturatedRings","QED"
]

def main():
    ap = argparse.ArgumentParser(description="Descriptors-only LR classification with OOF ranking")
    ap.add_argument("input", help="CSV with columns: smiles, ic50_nM, active + descriptor columns")
    ap.add_argument("--penalty", default="l2", choices=["l2","elasticnet"], help="LR penalty")
    ap.add_argument("--l1_ratio", type=float, default=0.1, help="ElasticNet l1_ratio if penalty=elasticnet")
    ap.add_argument("--folds", type=int, default=5, help="Stratified K folds")
    ap.add_argument("--thresh", type=float, default=0.80, help="P(active) threshold for hits")
    ap.add_argument("--topN", type=int, default=50, help="Top-N rows to export")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--prefix", default="EGFR_DESC_LR", help="Output file prefix")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.input)
    for c in ["smiles","active"]:
        if c not in df.columns: raise ValueError(f"Missing required column: {c}")

    # pick descriptors that exist
    feat_cols = [c for c in DESC_CANDIDATES if c in df.columns]
    if not feat_cols:
        raise ValueError("No descriptor columns found. Expected some of: " + ", ".join(DESC_CANDIDATES))
    X = df[feat_cols].astype(float).values
    y = df["active"].astype(int).values

    # OOF probabilities
    skf = StratifiedKFold(n_splits=args.folds, shuffle=True, random_state=42)
    oof = np.zeros(len(df), dtype=float)

    # LR model (scaled)
    if args.penalty == "elasticnet":
        lr = LogisticRegression(
            penalty="elasticnet", l1_ratio=args.l1_ratio, solver="saga",
            C=1.0, max_iter=5000, class_weight="balanced", random_state=42
        )
    else:
        lr = LogisticRegression(
            penalty="l2", solver="lbfgs", C=1.0,
            max_iter=5000, class_weight="balanced", random_state=42
        )
    pipe = make_pipeline(StandardScaler(), lr)

    for i,(tr,te) in enumerate(skf.split(X,y), 1):
        pipe.fit(X[tr], y[tr])
        oof[te] = pipe.predict_proba(X[te])[:,1]
        print(f"fold {i}/{args.folds} done")

    roc = roc_auc_score(y, oof)
    pr  = average_precision_score(y, oof)
    print(f"OOF ROC-AUC={roc:.3f} | PR-AUC={pr:.3f}")

    # Attach OOF and simple ADMET gates (need those descriptor cols)
    df["p_active_desc_oof"] = oof
    df["Lipinski_OK"] = (df["MolWt"]<=500) & (df["LogP"]<=5) & (df["HBD"]<=5) & (df["HBA"]<=10)
    df["Veber_OK"]    = (df["TPSA"]<=140) & (df["RB"]<=10)
    df["QED_OK"]      = (df["QED"]>=0.40) if "QED" in df.columns else True
    df["ADMET_OK"]    = df["Lipinski_OK"] & df["Veber_OK"] & df["QED_OK"]

    # Save per-molecule predictions
    preds_path = os.path.join(args.outdir, f"{args.prefix}_OOF_preds.csv")
    df.to_csv(preds_path, index=False)
    print("Saved:", preds_path)

    # Make a Top-N list (rank by prob, gate by ADMET)
    hits = df[df["ADMET_OK"] & (df["p_active_desc_oof"] >= args.thresh)].copy()
    cols = ["smiles","p_active_desc_oof","ic50_nM","MolWt","LogP","TPSA","QED","Lipinski_OK","Veber_OK"]
    cols = [c for c in cols if c in hits.columns]
    topN = hits.sort_values("p_active_desc_oof", ascending=False)[cols].head(args.topN)
    top_path = os.path.join(args.outdir, f"{args.prefix}_Top{args.topN}.csv")
    topN.to_csv(top_path, index=False)
    print(f"High-confidence, drug-like hits: {len(hits)} (showing Top-{args.topN})")
    print("Saved:", top_path)

    # Train once on ALL data for interpretability (coef on scaled features)
    pipe.fit(X, y)
    # Pull coefficients back to feature names
    # (coef_ corresponds to positive class; StandardScaler is inside the pipeline)
    coefs = pipe.named_steps["logisticregression"].coef_.ravel()
    coef_df = pd.DataFrame({"feature": feat_cols, "coef": coefs}).sort_values("coef", ascending=False)
    coef_path = os.path.join(args.outdir, f"{args.prefix}_coefficients.csv")
    coef_df.to_csv(coef_path, index=False)
    print("Saved:", coef_path)
    print("Top +features:\n", coef_df.head(5).to_string(index=False))
    print("Top -features:\n", coef_df.tail(5).to_string(index=False))

if __name__ == "__main__":
    main()
