#!/usr/bin/env python3
# chembl_egfr_clean.py
# -------------------------------------------------------------
# CLI tool to clean ChEMBL EGFR (CHEMBL203) bioactivity dumps
# into a QSAR-ready table: smiles, ic50_nM, active
#
# Usage:
#   python chembl_egfr_clean.py INPUT.[csv|tsv|gz] \
#       --outdir . --prefix EGFR \
#       --active-threshold-nm 1000 \
#       [--nsclc-only]
#
# Requires: pandas, numpy, rdkit
#   pip install pandas numpy rdkit-pypi
# -------------------------------------------------------------

import argparse, os, sys, time, gzip, io, re, textwrap
import numpy as np
import pandas as pd

# RDKit for SMILES standardization (strip salts/counter-ions)
try:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
except Exception as e:
    sys.stderr.write("ERROR: RDKit is required. Install with: pip install rdkit-pypi\n")
    raise

def now(): return time.strftime("%Y-%m-%d %H:%M:%S")

def sniff_sep_and_text(path):
    """Detects gzip, reads a small chunk, returns (is_gz, sample_text, sep)"""
    is_gz = False
    with open(path, "rb") as f:
        sig = f.read(2)
        is_gz = (sig == b"\x1f\x8b")
    if is_gz:
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as g:
            sample = g.read(4096)
    else:
        with open(path, "rb") as f:
            raw = f.read(4096)
        try:
            sample = raw.decode("utf-8")
        except UnicodeDecodeError:
            sample = raw.decode("latin-1", errors="ignore")
    # crude delimiter guess: choose the most frequent among ',', '\t', ';'
    counts = {",": sample.count(","), "\t": sample.count("\t"), ";": sample.count(";")}
    sep = max(counts, key=counts.get) if any(counts.values()) else ","
    return is_gz, sample, sep

def read_table_safely(path):
    """Robust reader for CSV/TSV/semicolon CSV; gz ok; returns pandas DataFrame (all str)."""
    _, sample, sep = sniff_sep_and_text(path)
    # If it looks like HTML, stop early
    if "<html" in sample.lower():
        raise ValueError("Input looks like an HTML page, not a data file. Please re-download as CSV/TSV from ChEMBL.")
    df = pd.read_csv(path, sep=sep, engine="python", dtype=str, on_bad_lines="skip")
    return df

def case_insensitive_picker(columns):
    """Returns a helper pick(*aliases) -> actual_column_name or None"""
    lookup = {c.lower(): c for c in columns}
    def pick(*names):
        for n in names:
            if n and n.lower() in lookup:
                return lookup[n.lower()]
        return None
    return pick

def nm_from_value_and_unit(val, unit):
    """Convert any value/unit to nM. Returns float or np.nan."""
    try:
        x = float(str(val).replace(",", "").strip())
    except Exception:
        return np.nan
    u = (str(unit).lower().strip() if unit is not None and str(unit).strip() != "" else "nm")
    if u in ("nm", "nanomolar"):               return x
    if u in ("µm", "um", "μm", "micromolar"):  return x * 1000.0
    if u in ("pm", "picomolar"):               return x / 1000.0
    if u in ("mm", "millimolar", "mmol"):      return x * 1_000_000.0
    # other weird units (% inhibition, rates, etc.) -> unknown
    return np.nan

def standardize_parent_smiles(smi):
    """Strip salts/solvents and uncharge. Return standardized SMILES or None if parse fails."""
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        parent = rdMolStandardize.FragmentParent(m)      # drops ".Cl", ".Na", etc.
        unchg  = rdMolStandardize.Uncharger().uncharge(parent)
        Chem.SanitizeMol(unchg)
        return Chem.MolToSmiles(unchg, isomericSmiles=True)
    except Exception:
        return None

def main():
    ap = argparse.ArgumentParser(
        description="Clean ChEMBL EGFR (CHEMBL203) activity into QSAR-ready CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ap.add_argument("input", help="Input CSV/TSV (semicolon CSV also ok; .gz supported).")
    ap.add_argument("--target", default="CHEMBL203", help="Target ChEMBL ID to keep.")
    ap.add_argument("--std-type", default="IC50", help="Standard Type to keep (e.g., IC50).")
    ap.add_argument("--assay-type", default="B", help="Assay Type to keep (B = Binding).")
    ap.add_argument("--active-threshold-nm", type=float, default=1000.0, help="Active cutoff in nM.")
    ap.add_argument("--nsclc-only", action="store_true", help="Keep rows whose assay_description mentions NSCLC mutants/cell lines.")
    ap.add_argument("--prefix", default="EGFR", help="Output file prefix.")
    ap.add_argument("--outdir", default=".", help="Output directory.")
    args = ap.parse_args()

    t0 = time.time()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"[{now()}] Reading: {args.input}")
    df_raw = read_table_safely(args.input)
    n0 = len(df_raw)
    print(f"  loaded rows={n0}, cols={len(df_raw.columns)}")

    # --- Column mapping to unified names ---
    pick = case_insensitive_picker(df_raw.columns)
    colmap = {
        "canonical_smiles":     pick("canonical_smiles","Smiles","smiles"),
        "standard_value":       pick("standard_value","Standard Value"),
        "standard_units":       pick("standard_units","Standard Units"),
        "standard_relation":    pick("standard_relation","Standard Relation"),
        "pchembl_value":        pick("pchembl_value","pChEMBL Value"),
        "standard_type":        pick("standard_type","Standard Type"),
        "assay_type":           pick("assay_type","Assay Type"),
        "assay_description":    pick("assay_description","Assay Description"),
        "assay_variant_mutation": pick("assay_variant_mutation","Assay Variant Mutation"),
        "target_chembl_id":     pick("target_chembl_id","Target ChEMBL ID"),
        "molecule_chembl_id":   pick("molecule_chembl_id","Molecule ChEMBL ID"),
        "molecule_pref_name":   pick("molecule_pref_name","Molecule Name"),
    }
    keep_src = [v for v in colmap.values() if v]
    missing_core = [k for k in ("canonical_smiles","standard_value","standard_units","standard_type","assay_type","target_chembl_id") if not colmap.get(k)]
    if missing_core:
        sys.stderr.write(f"ERROR: Required columns missing in input: {missing_core}\n")
        sys.exit(2)

    df = df_raw[keep_src].copy()
    # Rename to unified keys (lower snake_case)
    df = df.rename(columns={v:k for k,v in colmap.items() if v})
    # Uppercase a few fields for easy compares
    for c in ("standard_type","assay_type","target_chembl_id"):
        if c in df.columns:
            df[c] = df[c].astype(str).str.upper()

    # --- Optional NSCLC focus (by keywords in assay_description) ---
    if args.nsclc_only and "assay_description" in df.columns:
        patt = r"L858R|T790M|C797S|Del19|Exon 19|H1975|HCC827|PC9|H1650|HCC4006"
        before = len(df)
        df = df[df["assay_description"].fillna("").str.contains(patt, case=False, regex=True)].copy()
        print(f"NSCLC filter: {before} -> {len(df)}")

    # --- Core filters: target + IC50 + Binding ---
    before = len(df)
    df = df[
        (df["target_chembl_id"] == args.target.upper()) &
        (df["standard_type"] == args.std_type.upper()) &
        (df["assay_type"] == args.assay_type.upper())
    ].copy()
    print(f"Core filters (target={args.target}, type={args.std_type}, assay={args.assay_type}): {before} -> {len(df)}")

    # --- Require SMILES and numeric value present ---
    before = len(df)
    df = df[df["canonical_smiles"].notna() & df["standard_value"].notna()].copy()
    print(f"Drop missing SMILES/value: {before} -> {len(df)}")

    # --- Normalize relation and keep '=', '<', '<=' (drop '>', '~', etc.) ---
    if "standard_relation" in df.columns:
        df["standard_relation"] = (
            df["standard_relation"].astype(str)
              .str.replace("'", "", regex=False)
              .str.replace('"', "", regex=False)
              .str.strip()
        )
        before = len(df)
        df = df[df["standard_relation"].isin(["=", "<", "<="])].copy()
        print(f"Keep relations '=', '<', '<=': {before} -> {len(df)}")

    # --- Convert to nM, backfill from pChEMBL if needed ---
    before = len(df)
    df["ic50_nM"] = [nm_from_value_and_unit(v, u) for v, u in zip(df["standard_value"], df["standard_units"])]
    if "pchembl_value" in df.columns:
        miss = df["ic50_nM"].isna()
        if miss.any():
            pc = pd.to_numeric(df.loc[miss, "pchembl_value"], errors="coerce")
            df.loc[miss, "ic50_nM"] = 10 ** (9 - pc)  # pChEMBL -> nM
    df = df.dropna(subset=["ic50_nM"]).copy()
    print(f"Unit conversion to nM (+pChEMBL backfill): {before} -> {len(df)}")

    # --- Standardize SMILES to parent (strip salts, uncharge) ---
    before = len(df)
    df["smiles_std"] = df["canonical_smiles"].apply(standardize_parent_smiles)
    df = df.dropna(subset=["smiles_std"]).copy()
    print(f"SMILES standardization (strip salts/uncharge): {before} -> {len(df)}")

    # --- Dedupe by standardized SMILES; keep best (lowest IC50) ---
    before = len(df)
    df = df.sort_values("ic50_nM", ascending=True).drop_duplicates(subset=["smiles_std"], keep="first")
    print(f"Deduplicate by structure (keep lowest IC50): {before} -> {len(df)}")

    # --- Build minimal QSAR table ---
    out_qsar = df[["smiles_std","ic50_nM"]].rename(columns={"smiles_std":"smiles"}).copy()
    out_qsar["ic50_nM"] = out_qsar["ic50_nM"].astype(float)
    out_qsar["active"] = (out_qsar["ic50_nM"] <= float(args.active_threshold_nm)).astype(int)

    # Basic sanity stats (pIC50 range)
    out_qsar["pIC50_from_nM"] = 9 - np.log10(out_qsar["ic50_nM"])
    pmin, pmed, pmax = out_qsar["pIC50_from_nM"].min(), out_qsar["pIC50_from_nM"].median(), out_qsar["pIC50_from_nM"].max()

    # --- Save files ---
    prefix = args.prefix
    qsar_path = os.path.join(args.outdir, f"{prefix}_EGFR_QSAR_dataset.csv")
    out_qsar.to_csv(qsar_path, index=False)

    # Also save a slim metadata join to help later reporting
    # (if molecule IDs/names exist)
    meta_cols = [c for c in ("molecule_chembl_id","molecule_pref_name","assay_description") if c in df.columns]
    if meta_cols:
        meta = df[["smiles_std"] + meta_cols].drop_duplicates("smiles_std").rename(columns={"smiles_std":"smiles"})
        meta_path = os.path.join(args.outdir, f"{prefix}_EGFR_QSAR_metadata.csv")
        meta.to_csv(meta_path, index=False)
    else:
        meta_path = None

    # --- Print summary ---
    n_final = len(out_qsar)
    n_act = int(out_qsar["active"].sum())
    n_inact = int(n_final - n_act)
    took = time.time() - t0

    print("\n=== SUMMARY ===")
    print(f"Input file         : {args.input}")
    print(f"Output QSAR file   : {qsar_path}")
    if meta_path:
        print(f"Output meta file   : {meta_path}")
    print(f"Target / Type / Assay : {args.target} / {args.std_type} / {args.assay_type}")
    if args.nsclc_only:
        print("Assay context      : NSCLC-focused only (mutants/cell lines keywords)")
    print(f"Final rows         : {n_final} (unique molecules)")
    print(f"Actives/Inactives : {n_act} / {n_inact} (threshold <= {args.active_threshold_nm:.1f} nM)")
    print(f"pIC50 (min/med/max): {pmin:.3f} / {pmed:.3f} / {pmax:.3f}")
    print(f"Time taken         : {took:.1f} s")
    print("================\n")

if __name__ == "__main__":
    main()
