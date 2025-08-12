import { createTheme } from "@mui/material/styles";

const theme = createTheme({
  palette: {
    mode: "dark",
    background: { default: "#0A0F1C", paper: "rgba(255,255,255,0.04)" },
    primary:   { main: "#4F8BFF" },
    secondary: { main: "#9D7DFF" },
    success:   { main: "#27D39F" },
    text:      { primary: "rgba(255,255,255,0.92)", secondary: "rgba(255,255,255,0.72)" },
    divider:   "rgba(255,255,255,0.08)",
  },
  typography: {
    fontFamily: "InterVariable, Inter, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial",
    h4: { fontWeight: 800, letterSpacing: -0.2 },
    h6: { fontWeight: 700 },
    button: { textTransform: "none", fontWeight: 700, letterSpacing: 0.2 },
  },
  shape: { borderRadius: 16 },
  components: {
    MuiCssBaseline: {
      styleOverrides: {
        body: {
          backgroundAttachment: "fixed",
          backgroundImage:
            "radial-gradient(1200px 600px at 10% -10%, rgba(79,139,255,0.18), transparent 60%)," +
            "radial-gradient(900px 500px at 90% 0%, rgba(157,125,255,0.14), transparent 55%)",
        },
      },
    },
    MuiPaper: {
      styleOverrides: {
        root: {
          backdropFilter: "blur(8px)",
          border: "1px solid rgba(255,255,255,0.08)",
        },
      },
    },
    MuiTooltip: { styleOverrides: { tooltip: { fontSize: "0.9rem" } } },
    MuiButton:  { styleOverrides: { root: { borderWidth: 2 } } },
    MuiChip:    { styleOverrides: { root: { fontWeight: 700 } } },
  },
});

export default theme;
