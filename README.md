# Drug-Discovery

## In egfr-pipeline
npm i @mui/material @mui/icons-material @emotion/react @emotion/styled \
      @mui/x-data-grid papaparse @fontsource-variable/inter

## Frontend
cd frontend
npm i
npm run dev

## Backend (Windows/venv)
cd backend
python -m venv .venv
.\.venv\Scripts\Activate
pip install -r requirements.txt
python -m uvicorn main:app --reload --port 8000
