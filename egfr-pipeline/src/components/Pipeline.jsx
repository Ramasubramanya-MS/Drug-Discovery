// src/components/Pipeline.jsx
import * as React from "react";
import Box from "@mui/material/Box";
import Tooltip from "@mui/material/Tooltip";
import IconButton from "@mui/material/IconButton";
import Button from "@mui/material/Button";
import TextField from "@mui/material/TextField";
import Snackbar from "@mui/material/Snackbar";
import Alert from "@mui/material/Alert";
import Chip from "@mui/material/Chip";
import CircularProgress from "@mui/material/CircularProgress";
import CloudUploadIcon from "@mui/icons-material/CloudUpload";
import CleaningServicesIcon from "@mui/icons-material/CleaningServices";
import ScienceIcon from "@mui/icons-material/Science";
import OutboxIcon from "@mui/icons-material/Outbox";
import PlayArrowIcon from "@mui/icons-material/PlayArrow";
import CheckCircleIcon from "@mui/icons-material/CheckCircle";
import ErrorIcon from "@mui/icons-material/Error";
import { runUpload, runClean, runLearn, runOutput } from "../lib/api.js";

import UploadCard from "./UploadCard.jsx";
import ResultsTable from "./ResultsTable.jsx";
import styles from "./pipeline.module.css";
import { API_BASE } from "../lib/config";

const baseSteps = [
  {
    key: "upload",
    label: "Upload",
    icon: <CloudUploadIcon fontSize="small" />,
    tooltip:
      "Upload a CSV/SDF for CHEMBL203 (EGFR). Required fields for CSV: smiles, activity/label.",
    action: runUpload,
  },
  {
    key: "clean",
    label: "Clean",
    icon: <CleaningServicesIcon fontSize="small" />,
    tooltip:
      "Standardize SMILES, deduplicate by structure, convert IC50 to nM, keep valid assay rows.",
    action: runClean,
  },
  {
    key: "learn",
    label: "Learn",
    icon: <ScienceIcon fontSize="small" />,
    tooltip:
      "Compute descriptors then train the classifier; returns Top-N predictions and coefficients.",
    action: runLearn, // this performs descriptors -> train
  },
  {
    key: "output",
    label: "Output",
    icon: <OutboxIcon fontSize="small" />,
    tooltip:
      "Download Top-N shortlist (CSV) and other artifacts (predictions, coefficients).",
    action: runOutput,
  },
];

export default function Pipeline() {
  // per-step status
  const [statuses, setStatuses] = React.useState({
    upload: "idle",
    clean: "idle",
    learn: "idle",
    output: "idle",
  });

  // messages & toast
  const [messages, setMessages] = React.useState([]);
  const [toast, setToast] = React.useState({
    open: false,
    message: "",
    severity: "success",
  });

  // inputs
  const [topN, setTopN] = React.useState(50);
  const [uploadPayload, setUploadPayload] = React.useState(null);

  // context shared across steps (urls, previews)
  const [ctx, setCtx] = React.useState({});

  const stepIndex = (k) => baseSteps.findIndex((s) => s.key === k);
  const isEnabled = (k) => {
    const idx = stepIndex(k);
    for (let i = 0; i < idx; i++) {
      if (statuses[baseSteps[i].key] !== "done") return false;
    }
    return true;
  };

  const updateStatus = (key, s) =>
    setStatuses((prev) => ({ ...prev, [key]: s }));

  const runStep = async (stepKey) => {
    if (!isEnabled(stepKey)) return;
    const step = baseSteps.find((s) => s.key === stepKey);
    try {
      updateStatus(stepKey, "running");

      const payload = {
        upload: uploadPayload,
        gate: { topN }, // reused name to avoid refactor in api.js
      };

      const res = await step.action(payload, ctx);

      // Merge returned context fields (only those present)
      const {
        qsar_url,
        descriptors_url,
        top_url,
        preds_url,
        coefs_url,
        preview,
        top_preview,
        logs,
      } = res || {};

      setCtx((prev) => ({
        ...prev,
        ...(qsar_url ? { qsar_url } : {}),
        ...(descriptors_url ? { descriptors_url } : {}),
        ...(top_url ? { top_url } : {}),
        ...(preds_url ? { preds_url } : {}),
        ...(coefs_url ? { coefs_url } : {}),
        ...(preview ? { preview } : {}),
        ...(top_preview ? { top_preview } : {}),
      }));

      if (logs) setMessages((m) => [...m, logs]);

      updateStatus(stepKey, "done");
      const msg = `${step.label}: ${res?.message ?? "Completed."}`;
      setMessages((m) => [...m, msg]);
      setToast({ open: true, message: msg, severity: "success" });
    } catch (e) {
      console.error(e);
      updateStatus(stepKey, "error");
      const extra = e?.body?.logs || e?.body?.traceback;
      if (extra) setMessages((m) => [...m, `${step.label} logs:\n${extra}`]);
      setToast({
        open: true,
        message: `${step.label}: Failed. ${e?.message ?? ""}`,
        severity: "error",
      });
    }
  };

  const StatusChip = ({ s }) => {
    if (s === "running")
      return (
        <Chip
          size="small"
          icon={<CircularProgress size={14} />}
          label="Running…"
        />
      );
    if (s === "done")
      return (
        <Chip
          size="small"
          color="success"
          icon={<CheckCircleIcon />}
          label="Done"
        />
      );
    if (s === "error")
      return <Chip size="small" color="error" icon={<ErrorIcon />} label="Error" />;
    return <Chip size="small" variant="outlined" label="Idle" />;
  };

  return (
    <>
      <Box className={styles.pipeline}>
        {/* Upload */}
        <Tooltip title={baseSteps[0].tooltip} placement="top">
          <Box className={styles.card} data-enabled={isEnabled("upload")} data-state={statuses.upload}>
            <Box className={styles.header}>
              <Box className={styles.title}>
                {baseSteps[0].icon}
                <span>{baseSteps[0].label}</span>
              </Box>
              <StatusChip s={statuses.upload} />
            </Box>

            <UploadCard
              disabled={statuses.upload === "running"}
              onChange={(payload) => setUploadPayload(payload)}
            />

            <Box className={styles.actions}>
              <IconButton
                aria-label="Run Upload"
                onClick={() => runStep("upload")}
                disabled={
                  !isEnabled("upload") ||
                  statuses.upload === "running" ||
                  !uploadPayload
                }
              >
                <PlayArrowIcon />
              </IconButton>
            </Box>
          </Box>
        </Tooltip>

        <Box className={styles.connector} />

        {/* Clean */}
        <Tooltip title={baseSteps[1].tooltip} placement="top">
          <Box className={styles.card} data-enabled={isEnabled("clean")}>
            <Box className={styles.header}>
              <Box className={styles.title}>
                {baseSteps[1].icon}
                <span>{baseSteps[1].label}</span>
              </Box>
              <StatusChip s={statuses.clean} />
            </Box>

            <Box className={styles.body}>
              Standardize, deduplicate, convert units; produce QSAR-ready dataset.
            </Box>

            <Box className={styles.actions}>
              <IconButton
                aria-label="Run Clean"
                onClick={() => runStep("clean")}
                disabled={!isEnabled("clean") || statuses.clean === "running"}
              >
                <PlayArrowIcon />
              </IconButton>
            </Box>
          </Box>
        </Tooltip>

        <Box className={styles.connector} />

        {/* Learn (descriptors -> train + Top-N) */}
        <Tooltip title={baseSteps[2].tooltip} placement="top">
          <Box className={styles.card} data-enabled={isEnabled("learn")}>
            <Box className={styles.header}>
              <Box className={styles.title}>
                {baseSteps[2].icon}
                <span>{baseSteps[2].label}</span>
              </Box>
              <StatusChip s={statuses.learn} />
            </Box>

            <Box className={styles.body} sx={{ display: "grid", gap: 12 }}>
              Compute descriptors and train classifier; show Top-N shortlist.
              <TextField
                label="Top-N"
                type="number"
                size="small"
                value={topN}
                onChange={(e) => setTopN(Number(e.target.value))}
                inputProps={{ min: 1 }}
              />
            </Box>

            <Box className={styles.actions}>
              <IconButton
                aria-label="Run Learn"
                onClick={() => runStep("learn")}
                disabled={!isEnabled("learn") || statuses.learn === "running"}
              >
                <PlayArrowIcon />
              </IconButton>
            </Box>
          </Box>
        </Tooltip>

        <Box className={styles.connector} />

        {/* Output */}
        <Tooltip title={baseSteps[3].tooltip} placement="top">
          <Box className={styles.card} data-enabled={isEnabled("output")}>
            <Box className={styles.header}>
              <Box className={styles.title}>
                {baseSteps[3].icon}
                <span>{baseSteps[3].label}</span>
              </Box>
              <StatusChip s={statuses.output} />
            </Box>

            <Box className={styles.body}>
              Download Top-N CSV (and optional predictions/coefficients).
            </Box>

            <Box className={styles.actions} sx={{ gap: 8 }}>
              <IconButton
                aria-label="Run Output"
                onClick={() => runStep("output")}
                disabled={!isEnabled("output") || statuses.output === "running"}
              >
                <PlayArrowIcon />
              </IconButton>

              <Button
                size="small"
                variant="outlined"
                startIcon={<OutboxIcon />}
                onClick={() => {
                  if (!ctx?.top_url) {
                    setToast({
                      open: true,
                      message: "Run Learn first to generate Top-N.",
                      severity: "info",
                    });
                    return;
                  }
                  window.open(`${API_BASE}${ctx.top_url}`, "_blank");
                }}
              >
                Download Top-N
              </Button>
            </Box>
          </Box>
        </Tooltip>
      </Box>

      {/* Top-N table appears right after Learn */}
      {ctx?.top_url && (
        <ResultsTable csvUrl={ctx.top_url} title={`Top-${topN} Shortlist`} />
      )}

      {/* Run Log */}
      <Box sx={{ mt: 3, p: 2, bgcolor: "rgba(255,255,255,0.04)", borderRadius: 2 }}>
        <Box sx={{ fontWeight: 600, mb: 1 }}>Run Log</Box>
        {messages.length === 0 ? (
          <Box sx={{ opacity: 0.7 }}>No runs yet.</Box>
        ) : (
          <Box component="ul" sx={{ m: 0, pl: 3, whiteSpace: "pre-wrap" }}>
            {messages.map((m, i) => (
              <li key={i}>{m}</li>
            ))}
          </Box>
        )}
      </Box>

      <Snackbar
        open={toast.open}
        autoHideDuration={3200}
        onClose={() => setToast((t) => ({ ...t, open: false }))}
      >
        <Alert
          severity={toast.severity}
          onClose={() => setToast((t) => ({ ...t, open: false }))}
          variant="filled"
          sx={{ width: "100%" }}
        >
          {toast.message}
        </Alert>
      </Snackbar>
    </>
  );
}
