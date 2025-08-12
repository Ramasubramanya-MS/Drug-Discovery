#!/usr/bin/env python3
"""
make_descriptors.py
-------------------
CLI tool: read a QSAR table (needs a SMILES column), compute RDKit descriptors,
and write a new CSV with descriptors appended.

Typical input schema (others are carried through if present):
  - smiles (required)
  - ic50_nM (optional)
  - active  (optional)

Usage:
  python make_descriptors.py INPUT.csv --out OUTPUT.csv
  python make_descriptors.py "C:\path with spaces\EGFR_QSAR_dataset.csv" --standardize

Args:
  INPUT                 Path to CSV/TSV/semicolon-CSV (gz ok). Must contain a SMILES column.

Options:
  --out OUT             Output CSV path (default: <input_basename>_with_descriptors.csv)
  --smiles-col NAME     Column name for SMILES (case-insensitive; default: smiles)
  --no-standardize      Do NOT standardize SMILES (keep as-is). Default is to standardize.
  --quiet               Less console output.

Requires:
  - pandas, numpy, rdkit (pip install rdkit-pypi pandas numpy)
"""

import argparse, os, sys, time, gzip, io
import numpy as np
import pandas as pd

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, QED
    from rdkit.Chem.MolStandardize import rdMolStandardize
except Exception as e:
    sys.stderr.write("ERROR: RDKit is required. Install with: pip install rdkit-pypi\n")
    raise

def sniff_sep_and_read(path: str) -> pd.DataFrame:
    """Read CSV/TSV/semicolon CSV; gz supported. Returns all columns as-is."""
    def _is_gz(p): 
        with open(p, "rb") as f: 
            return f.read(2) == b"\x1f\x8b"
    is_gz = _is_gz(path)
    # peek some text to guess delimiter
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
    if "<html" in sample.lower():
        raise ValueError("Input looks like HTML, not a data file. Re-download as CSV/TSV.")
    counts = {",": sample.count(","), "\t": sample.count("\t"), ";": sample.count(";")}
    sep = max(counts, key=counts.get) if any(counts.values()) else ","
    return pd.read_csv(path, sep=sep, engine="python", dtype=str, on_bad_lines="skip")

def pick_column(df: pd.DataFrame, name: str) -> str:
    """Case-insensitive column picker; returns actual name or raises."""
    lookup = {c.lower(): c for c in df.columns}
    key = name.lower()
    if key in lookup:
        return lookup[key]
    # try common variants for smiles
    if key == "smiles":
        for alt in ("canonical_smiles","smile","SMILES"):
            if alt.lower() in lookup:
                return lookup[alt.lower()]
    raise KeyError(f"Column '{name}' not found (available: {list(df.columns)[:8]}...)")

def to_parent_smiles(smi: str):
    """Strip salts/solvents and uncharge. Return standardized SMILES or None."""
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        p = rdMolStandardize.FragmentParent(m)
        p = rdMolStandardize.Uncharger().uncharge(p)
        Chem.SanitizeMol(p)
        return Chem.MolToSmiles(p, isomericSmiles=True)
    except Exception:
        return None

def compute_desc_row(smi: str):
    """Compute a compact, interpretable descriptor set for one SMILES."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    try:
        return {
            # size / composition
            "MolWt":            Descriptors.MolWt(m),
            "HeavyAtomCount":   Descriptors.HeavyAtomCount(m),
            "NumHeteroatoms":   rdMolDescriptors.CalcNumHeteroatoms(m),
            # lipophilicity / polarity
            "LogP":             Descriptors.MolLogP(m),         # Crippen cLogP
            "TPSA":             rdMolDescriptors.CalcTPSA(m),   # Topological polar surface area
            # hydrogen bonding
            "HBD":              rdMolDescriptors.CalcNumHBD(m),
            "HBA":              rdMolDescriptors.CalcNumHBA(m),
            # flexibility / rings / shape
            "RB":               rdMolDescriptors.CalcNumRotatableBonds(m),
            "FractionCSP3":     rdMolDescriptors.CalcFractionCSP3(m),
            "NumRings":         rdMolDescriptors.CalcNumRings(m),
            "NumAromaticRings": rdMolDescriptors.CalcNumAromaticRings(m),
            "NumAliphaticRings":rdMolDescriptors.CalcNumAliphaticRings(m),
            "NumSaturatedRings":rdMolDescriptors.CalcNumSaturatedRings(m),
            # heuristic drug-likeness
            "QED":              float(QED.qed(m)),
        }
    except Exception:
        return None

def main():
    ap = argparse.ArgumentParser(
        description="Compute RDKit descriptors from a SMILES column and write a new CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ap.add_argument("input", help="Input CSV/TSV (gz ok). Must contain a SMILES column.")
    ap.add_argument("--out", default=None, help="Output CSV path (default: <input_basename>_with_descriptors.csv)")
    ap.add_argument("--smiles-col", default="smiles", help="Name of the SMILES column (case-insensitive).")
    ap.add_argument("--no-standardize", action="store_true", help="Do NOT standardize SMILES (keep as-is).")
    ap.add_argument("--quiet", action="store_true", help="Reduce console output.")
    args = ap.parse_args()

    t0 = time.time()
    if not os.path.exists(args.input):
        sys.stderr.write(f"ERROR: Input not found: {args.input}\n")
        sys.exit(2)

    if not args.quiet:
        print(f"[INFO] Reading: {args.input}")

    df = sniff_sep_and_read(args.input)
    n_in = len(df)
    if not args.quiet:
        print(f"[INFO] Loaded rows={n_in}, cols={len(df.columns)}")

    # locate SMILES column
    try:
        smiles_col = pick_column(df, args.smiles_col)
    except KeyError as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(2)

    # pick optional passthrough columns if present
    passthrough = []
    for c in ("ic50_nM","active"):
        c_found = next((col for col in df.columns if col.lower()==c.lower()), None)
        if c_found: passthrough.append(c_found)

    # choose SMILES to featurize
    if args.no_standardize:
        smi_series = df[smiles_col].astype(str)
        standardizing = False
    else:
        if not args.quiet: print("[INFO] Standardizing SMILES (strip salts/uncharge)...")
        smi_std = df[smiles_col].astype(str).apply(to_parent_smiles)
        ok_mask = smi_std.notna()
        dropped_std = int((~ok_mask).sum())
        if dropped_std and not args.quiet:
            print(f"[WARN] {dropped_std} rows failed standardization and will be dropped.")
        df = df.loc[ok_mask].reset_index(drop=True).copy()
        df["smiles"] = smi_std[ok_mask].values
        smi_series = df["smiles"]
        standardizing = True

    # compute descriptors
    if not args.quiet: print("[INFO] Computing descriptors ...")
    rows = [compute_desc_row(s) for s in smi_series]
    ok_desc = np.array([r is not None for r in rows], dtype=bool)
    n_bad_desc = int((~ok_desc).sum())
    if n_bad_desc and not args.quiet:
        print(f"[WARN] {n_bad_desc} rows failed RDKit parsing and will be dropped.")

    X = pd.DataFrame([r for r in rows if r is not None]).reset_index(drop=True)
    df_ok = df.loc[np.where(ok_desc)[0]].reset_index(drop=True).copy()

    # build output table
    out_cols = ["smiles"] if standardizing else [smiles_col]
    out = df_ok[out_cols + passthrough].copy()
    out = out.rename(columns={smiles_col:"smiles"})  # unify name
    out = pd.concat([out.reset_index(drop=True), X], axis=1)

    # default output path
    if args.out:
        out_path = args.out
    else:
        base = os.path.basename(args.input)
        stem = os.path.splitext(os.path.splitext(base)[0])[0]  # handle .csv.gz
        out_path = os.path.join(os.path.dirname(args.input), f"{stem}_with_descriptors.csv")

    out.to_csv(out_path, index=False)

    # summary
    took = time.time() - t0
    kept = len(out)
    dropped = n_in - kept
    if not args.quiet:
        print("\n=== SUMMARY ===")
        print(f"Input rows        : {n_in}")
        print(f"Kept rows         : {kept}")
        print(f"Dropped (parse/std): {dropped}")
        print(f"Standardized      : {'Yes' if standardizing else 'No'}")
        # quick property medians
        mids = {k: float(out[k].median()) for k in ("MolWt","LogP","TPSA") if k in out.columns}
        print(f"Medians (MolWt/LogP/TPSA): {mids}")
        print(f"Output            : {out_path}")
        print(f"Took              : {took:.1f} s")
        print("===============")

if __name__ == "__main__":
    main()
