import * as React from "react";
import Container from "@mui/material/Container";
import Typography from "@mui/material/Typography";
import Box from "@mui/material/Box";
import Pipeline from "./components/Pipeline.jsx";

export default function App() {
  return (
    <Box sx={{ bgcolor: "#0f1115", minHeight: "100vh", color: "#e8eaf6", py: 6 }}>
      <Container maxWidth="xl">
        <Typography variant="h4" sx={{ fontWeight: 700, mb: 1 }}>
          EGFR Inhibitor Screening – AI Pipeline
        </Typography>
        <Typography variant="body1" sx={{ opacity: 0.8, mb: 4 }}>
          Upload ChEMBL data for CHEMBL203 (human EGFR), clean it, learn structure–activity,
          gate by drug-likeness, and emit a Top-N shortlist for literature cross-check
          (NSCLC mutants L858R/T790M).
        </Typography>
        <Pipeline />
      </Container>
    </Box>
  );
}
