import * as React from "react";
import Box from "@mui/material/Box";
import Typography from "@mui/material/Typography";
import { DataGrid } from "@mui/x-data-grid";
import Papa from "papaparse";
import { API_BASE } from "../lib/config";

export default function ResultsTable({ csvUrl, title = "Top-50 Predictions" }) {
  const [rows, setRows] = React.useState([]);
  const [cols, setCols] = React.useState([]);

  React.useEffect(() => {
    if (!csvUrl) return;
    fetch(`${API_BASE}${csvUrl}`)
      .then(r => r.text())
      .then(txt => {
        const parsed = Papa.parse(txt, { header: true, dynamicTyping: true, skipEmptyLines: true });
        const fields = parsed.meta.fields || [];
        const columns = fields.map(f => ({ field: f, headerName: f, flex: 1, minWidth: 120 }));
        const data = (parsed.data || []).map((r, i) => ({ id: i, ...r }));
        setCols(columns);
        setRows(data);
      });
  }, [csvUrl]);

  if (!csvUrl) return null;

  return (
    <Box sx={{ mt: 3 }}>
      <Typography variant="h6" sx={{ mb: 1 }}>{title}</Typography>
      <div style={{ height: 560, width: "100%" }}>
        <DataGrid
          rows={rows}
          columns={cols}
          pageSizeOptions={[10, 25, 50]}
          initialState={{ pagination: { paginationModel: { pageSize: 50 } } }}
        />
      </div>
    </Box>
  );
}
