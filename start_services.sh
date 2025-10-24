#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BACKEND_DIR="$PROJECT_ROOT/backend"
FRONTEND_DIR="$PROJECT_ROOT/frontend"
PYTHON_BIN="python3"

log() {
  echo "[start] $1"
}

# --- Backend env setup -----------------------------------------------------
if [[ ! -d "$BACKEND_DIR/.venv" ]]; then
  log "Creating backend virtual environment (.venv)"
  "$PYTHON_BIN" -m venv "$BACKEND_DIR/.venv"
fi

# shellcheck disable=SC1091
source "$BACKEND_DIR/.venv/bin/activate"

log "Installing backend dependencies"
pip install --upgrade pip >/dev/null
pip install -r "$BACKEND_DIR/requirements.txt" >/dev/null

# --- Frontend dependencies --------------------------------------------------
if [[ ! -d "$FRONTEND_DIR/node_modules" ]]; then
  log "Installing frontend dependencies"
  (cd "$FRONTEND_DIR" && npm install)
fi

# --- Service startup --------------------------------------------------------
log "Building frontend assets"
(cd "$FRONTEND_DIR" && npm run build >/dev/null)

cd "$BACKEND_DIR"
log "Starting FastAPI service on http://localhost:8080"
uvicorn app.main:app --host 0.0.0.0 --port 8080 --reload &
BACKEND_PID=$!
cd "$PROJECT_ROOT"

trap 'log "Stopping service"; kill $BACKEND_PID >/dev/null 2>&1 || true' EXIT
log "Service running. Press Ctrl+C to stop."
wait $BACKEND_PID
