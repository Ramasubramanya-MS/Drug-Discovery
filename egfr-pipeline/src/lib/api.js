// src/lib/api.js
import { API_BASE } from "./config";

const throwIfNotOk = async (res) => {
  const body = await res.json().catch(() => null);
  if (!res.ok) {
    const msg = body?.error || body?.detail || body?.logs || res.statusText;
    const e = new Error(typeof msg === "string" ? msg : "Request failed");
    e.body = body;
    throw e;
  }
  return body;
};

export async function runUpload(payload) {
  if (!payload?.upload?.file) throw new Error("No file selected.");
  return { message: `File "${payload.upload.name}" staged in UI.` };
}

export async function runClean(payload) {
  const f = payload?.upload?.file;
  if (!f) throw new Error("No file selected.");

  const fd = new FormData();
  fd.append("file", f);
  fd.append("prefix", "EGFR");
  fd.append("target", "CHEMBL203");
  fd.append("std_type", "IC50");
  fd.append("assay_type", "B");
  fd.append("active_threshold_nm", "1000");

  const res = await fetch(`${API_BASE}/pipeline/clean`, { method: "POST", body: fd });
  const data = await throwIfNotOk(res);
  return { message: "Cleaned dataset created.", qsar_url: data.qsar_url, logs: data.logs };
}

export async function runLearn(payload, ctx) {
  const qsar = ctx?.qsar_url;
  if (!qsar) throw new Error("Missing QSAR URL from Clean step.");
  const topN = (payload?.gate?.topN ?? 50);

  // Descriptors
  const fd1 = new FormData();
  fd1.append("input_path", qsar);
  fd1.append("smiles_col", "smiles");
  const r1 = await fetch(`${API_BASE}/pipeline/descriptors`, { method: "POST", body: fd1 });
  const desc = await throwIfNotOk(r1);

  // Train + Top-N
  const fd2 = new FormData();
  fd2.append("input_path", desc.descriptors_url);
  fd2.append("topN", String(topN));
  fd2.append("thresh", "0.80");
  fd2.append("folds", "5");
  fd2.append("penalty", "l2");
  fd2.append("l1_ratio", "0.1");
  fd2.append("prefix", "EGFR_DESC_LR");

  const r2 = await fetch(`${API_BASE}/pipeline/train`, { method: "POST", body: fd2 });
  const trn = await throwIfNotOk(r2);

  return {
    message: `Descriptors + Train complete. Selected Top-${topN}.`,
    descriptors_url: desc.descriptors_url,
    preview: desc.preview,
    logs: `${desc.logs}\n${trn.logs}`,
    top_url: trn.top_url,
    preds_url: trn.preds_url,
    coefs_url: trn.coefs_url,
    top_preview: trn.top_preview,
  };
}

export async function runOutput(payload, ctx) {
  if (!ctx?.top_url) throw new Error("Run Learn first to generate Top-N.");
  return { message: "Artifacts ready.", ...ctx };
}
