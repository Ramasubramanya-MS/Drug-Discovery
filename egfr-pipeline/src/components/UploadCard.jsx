import * as React from "react";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Typography from "@mui/material/Typography";
import CloudUploadIcon from "@mui/icons-material/CloudUpload";

export default function UploadCard({ disabled, onChange }) {
  const [fileName, setFileName] = React.useState("");

  const handlePick = (e) => {
    const f = e.target.files?.[0];
    if (!f) return;
    setFileName(f.name);
    onChange?.({ file: f, name: f.name });
  };

  return (
    <Box sx={{ display: "grid", gap: 12 }}>
      <Typography variant="body2" sx={{ opacity: 0.8 }}>
        Accepted: .csv (e.g., smiles,activity), .sdf
      </Typography>
      <Button
        component="label"
        variant="contained"
        startIcon={<CloudUploadIcon />}
        disabled={disabled}
      >
        Choose File
        <input type="file" hidden accept=".csv,.sdf" onChange={handlePick} />
      </Button>
      {fileName && (
        <Typography variant="caption" sx={{ opacity: 0.9 }}>
          Selected: {fileName}
        </Typography>
      )}
    </Box>
  );
}
