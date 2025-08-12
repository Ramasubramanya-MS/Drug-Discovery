from __future__ import annotations

"""
FastAPI backend to run the EGFR pipeline steps from your CLI tools:
  1) chembl_egfr_clean.py      -> /pipeline/clean
  2) make_descriptors.py        -> /pipeline/descriptors
  3) train_desc_classify.py     -> /pipeline/train

Returns JSON with captured logs and downloadable file URLs for your React app.
Place this file (main.py) in the same folder as the three scripts.
Run:  uvicorn main:app --reload --port 8000
"""

import os
import sys
import uuid
import shutil
import traceback
from pathlib import Path
from typing import Optional, Dict, Any

from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles

import subprocess
import time
import pandas as pd

# ----------------------------- Paths & Setup -----------------------------
BASE_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = BASE_DIR  # expects scripts in the same folder as main.py
WORKSPACE_DIR = BASE_DIR / "workspace"
WORKSPACE_DIR.mkdir(parents=True, exist_ok=True)

# Helpful: serve any generated files under /files/* (React can download)
app = FastAPI(title="EGFR Pipeline Backend", version="1.0.0")
app.mount("/files", StaticFiles(directory=str(WORKSPACE_DIR)), name="files")

# CORS: allow local dev
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:5173",
        "http://127.0.0.1:5173",
        "http://localhost:3000",
        "http://127.0.0.1:3000",
        "*",  # relax during dev; tighten for prod
    ],
    allow_methods=["*"],
    allow_headers=["*"],
)

# --------------------------- Utility Functions ---------------------------

def mk_run_dir(prefix: str = "job") -> Path:
    stamp = time.strftime("%Y%m%d-%H%M%S")
    rid = f"{prefix}-{stamp}-{uuid.uuid4().hex[:8]}"
    p = WORKSPACE_DIR / rid
    p.mkdir(parents=True, exist_ok=True)
    return p


def save_upload(dst_dir: Path, up: UploadFile, name: Optional[str] = None) -> Path:
    dst = dst_dir / (name or up.filename)
    with dst.open("wb") as f:
        shutil.copyfileobj(up.file, f)
    return dst


def run_script(script_name: str, args: list[str], cwd: Optional[Path] = None, timeout: int = 3600) -> dict:
    """Run a Python script using the current interpreter, capture stdout/stderr."""
    script_path = (SCRIPTS_DIR / script_name).resolve()
    if not script_path.exists():
        raise HTTPException(status_code=500, detail=f"Script not found: {script_path}")

    cmd = [sys.executable, str(script_path), *[str(a) for a in args]]
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        raise HTTPException(status_code=504, detail=f"Timed out running: {' '.join(cmd)}")

    return {
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "cmd": cmd,
    }


def url_for(path: Path) -> str:
    """Build a browser URL under /files for a workspace path."""
    try:
        rel = path.relative_to(WORKSPACE_DIR)
    except ValueError:
        # if not in workspace, expose nothing
        raise HTTPException(status_code=400, detail="Path is outside workspace")
    # use forward slashes for URLs
    return f"/files/{str(rel).replace(os.sep, '/')}"


def try_read_csv_preview(path: Path, n: int = 10) -> list[dict[str, Any]]:
    try:
        df = pd.read_csv(path)
        return df.head(n).to_dict(orient="records")
    except Exception:
        return []

# --------------------------------- Routes --------------------------------

@app.get("/health")
def health():
    return {"ok": True, "ts": time.time()}


@app.post("/pipeline/clean")
def clean(
    file: UploadFile = File(..., description="Raw ChEMBL dump (csv/tsv/csv.gz)"),
    prefix: str = Form("EGFR"),
    target: str = Form("CHEMBL203"),
    std_type: str = Form("IC50"),
    assay_type: str = Form("B"),
    active_threshold_nm: float = Form(1000.0),
    nsclc_only: bool = Form(False),
):
    run_dir = mk_run_dir("clean")
    raw_path = save_upload(run_dir, file, name=file.filename)

    args = [
        str(raw_path),
        "--outdir", str(run_dir),
        "--prefix", prefix,
        "--target", target,
        "--std-type", std_type,
        "--assay-type", assay_type,
        "--active-threshold-nm", str(active_threshold_nm),
    ]
    if nsclc_only:
        args.append("--nsclc-only")

    res = run_script("chembl_egfr_clean.py", args, cwd=run_dir)

    # Expected outputs
    qsar = run_dir / f"{prefix}_EGFR_QSAR_dataset.csv"
    meta = run_dir / f"{prefix}_EGFR_QSAR_metadata.csv"

    payload = {
        "ok": res["returncode"] == 0,
        "workdir": url_for(run_dir),
        "logs": res["stdout"] + ("\n" + res["stderr"] if res["stderr"] else ""),
        "qsar_url": url_for(qsar) if qsar.exists() else None,
        "meta_url": url_for(meta) if meta.exists() else None,
    }
    return JSONResponse(payload, status_code=200 if payload["ok"] else 500)


@app.post("/pipeline/descriptors")
def descriptors(
    # Either upload a QSAR CSV, or pass an existing workspace file URL
    file: Optional[UploadFile] = File(None, description="QSAR CSV from /pipeline/clean"),
    input_path: Optional[str] = Form(None, description="Existing /files/... URL or relative workspace path"),
    smiles_col: str = Form("smiles"),
    no_standardize: bool = Form(False),
):
    run_dir = mk_run_dir("desc")

    # Resolve input
    if file is not None:
        src_path = save_upload(run_dir, file, name=file.filename)
    elif input_path:
        # allow both "/files/<...>" and raw relative paths under workspace
        input_path = input_path.replace("\\", "/")
        if input_path.startswith("/files/"):
            input_path = input_path[len("/files/"):]
        src_path = WORKSPACE_DIR / input_path
        if not src_path.exists():
            raise HTTPException(status_code=404, detail=f"Input not found: {src_path}")
    else:
        raise HTTPException(status_code=400, detail="Provide either 'file' upload or 'input_path'.")

    # Compute output filename
    stem = src_path.stem
    if stem.endswith(".csv"):
        stem = stem[:-4]
    out_path = run_dir / f"{stem}_with_descriptors.csv"

    args = [str(src_path), "--out", str(out_path), "--smiles-col", smiles_col]
    if no_standardize:
        args.append("--no-standardize")

    res = run_script("make_descriptors.py", args, cwd=run_dir)

    payload = {
        "ok": res["returncode"] == 0,
        "workdir": url_for(run_dir),
        "logs": res["stdout"] + ("\n" + res["stderr"] if res["stderr"] else ""),
        "descriptors_url": url_for(out_path) if out_path.exists() else None,
        "preview": try_read_csv_preview(out_path, n=10),
    }
    return JSONResponse(payload, status_code=200 if payload["ok"] else 500)


@app.post("/pipeline/train")
def train(
    input_path: str = Form(..., description="/files/... URL or relative workspace path to descriptors CSV"),
    topN: int = Form(50),
    thresh: float = Form(0.80),
    folds: int = Form(5),
    penalty: str = Form("l2"),
    l1_ratio: float = Form(0.1),
    prefix: str = Form("EGFR_DESC_LR"),
):
    run_dir = mk_run_dir("train")

    input_path = input_path.replace("\\", "/")
    if input_path.startswith("/files/"):
        input_path = input_path[len("/files/"):]
    src_path = WORKSPACE_DIR / input_path
    if not src_path.exists():
        raise HTTPException(status_code=404, detail=f"Input not found: {src_path}")

    args = [
        str(src_path),
        "--outdir", str(run_dir),
        "--topN", str(topN),
        "--thresh", str(thresh),
        "--folds", str(folds),
        "--penalty", penalty,
        "--l1_ratio", str(l1_ratio),
        "--prefix", prefix,
    ]

    res = run_script("train_desc_classify.py", args, cwd=run_dir)

    preds = run_dir / f"{prefix}_OOF_preds.csv"
    top = run_dir / f"{prefix}_Top{topN}.csv"
    coefs = run_dir / f"{prefix}_coefficients.csv"

    payload = {
        "ok": res["returncode"] == 0,
        "workdir": url_for(run_dir),
        "logs": res["stdout"] + ("\n" + res["stderr"] if res["stderr"] else ""),
        "preds_url": url_for(preds) if preds.exists() else None,
        "top_url": url_for(top) if top.exists() else None,
        "coefs_url": url_for(coefs) if coefs.exists() else None,
        "top_preview": try_read_csv_preview(top, n=10),
    }
    return JSONResponse(payload, status_code=200 if payload["ok"] else 500)


# Optional: simple pipeline convenience combining all three steps in one call
# You can wire this to a single "Run All" button if desired.
@app.post("/pipeline/run-all")
def run_all(
    file: UploadFile = File(..., description="Raw ChEMBL dump (csv/tsv/csv.gz)"),
    prefix: str = Form("EGFR"),
    target: str = Form("CHEMBL203"),
    std_type: str = Form("IC50"),
    assay_type: str = Form("B"),
    active_threshold_nm: float = Form(1000.0),
    nsclc_only: bool = Form(False),
    topN: int = Form(50),
    thresh: float = Form(0.80),
    folds: int = Form(5),
    penalty: str = Form("l2"),
    l1_ratio: float = Form(0.1),
    model_prefix: str = Form("EGFR_DESC_LR"),
):
    """Run Clean -> Descriptors -> Train in one request and return all artifacts."""
    # 1) Clean
    clean_dir = mk_run_dir("clean")
    raw_path = save_upload(clean_dir, file, name=file.filename)
    clean_args = [
        str(raw_path), "--outdir", str(clean_dir), "--prefix", prefix,
        "--target", target, "--std-type", std_type, "--assay-type", assay_type,
        "--active-threshold-nm", str(active_threshold_nm),
    ]
    if nsclc_only:
        clean_args.append("--nsclc-only")
    r1 = run_script("chembl_egfr_clean.py", clean_args, cwd=clean_dir)
    qsar = clean_dir / f"{prefix}_EGFR_QSAR_dataset.csv"
    if not qsar.exists():
        raise HTTPException(status_code=500, detail="QSAR dataset not produced during clean stage.")

    # 2) Descriptors
    desc_dir = mk_run_dir("desc")
    out_desc = desc_dir / f"{prefix}_EGFR_QSAR_dataset_with_descriptors.csv"
    r2 = run_script("make_descriptors.py", [str(qsar), "--out", str(out_desc)], cwd=desc_dir)
    if not out_desc.exists():
        raise HTTPException(status_code=500, detail="Descriptor CSV not produced.")

    # 3) Train
    train_dir = mk_run_dir("train")
    r3 = run_script(
        "train_desc_classify.py",
        [
            str(out_desc), "--outdir", str(train_dir), "--topN", str(topN),
            "--thresh", str(thresh), "--folds", str(folds),
            "--penalty", penalty, "--l1_ratio", str(l1_ratio),
            "--prefix", model_prefix,
        ],
        cwd=train_dir,
    )

    preds = train_dir / f"{model_prefix}_OOF_preds.csv"
    top = train_dir / f"{model_prefix}_Top{topN}.csv"
    coefs = train_dir / f"{model_prefix}_coefficients.csv"

    return {
        "ok": (r1["returncode"] == 0 and r2["returncode"] == 0 and r3["returncode"] == 0),
        "clean": {
            "workdir": url_for(clean_dir),
            "logs": r1["stdout"] + ("\n" + r1["stderr"] if r1["stderr"] else ""),
            "qsar_url": url_for(qsar),
        },
        "descriptors": {
            "workdir": url_for(desc_dir),
            "logs": r2["stdout"] + ("\n" + r2["stderr"] if r2["stderr"] else ""),
            "descriptors_url": url_for(out_desc),
            "preview": try_read_csv_preview(out_desc, 10),
        },
        "train": {
            "workdir": url_for(train_dir),
            "logs": r3["stdout"] + ("\n" + r3["stderr"] if r3["stderr"] else ""),
            "preds_url": url_for(preds) if preds.exists() else None,
            "top_url": url_for(top) if top.exists() else None,
            "coefs_url": url_for(coefs) if coefs.exists() else None,
            "top_preview": try_read_csv_preview(top, 10),
        },
    }


# ------------------------------ Error Handler ----------------------------
@app.exception_handler(Exception)
async def unhandled_exception_handler(request, exc):
    # Keep responses JSON for the frontend
    tb = ''.join(traceback.format_exception(type(exc), exc, exc.__traceback__))
    return JSONResponse(
        status_code=500,
        content={
            "ok": False,
            "error": str(exc),
            "traceback": tb,
        },
    )
