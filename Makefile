.PHONY: demo-on demo-propka backend frontend

backend:
	npm --prefix frontend install
	npm --prefix frontend run build
	cd backend && USE_REAL_PROPka=$${USE_REAL_PROPka:-false} USE_REAL_DOCKING=$${USE_REAL_DOCKING:-false} uvicorn app.main:app --reload --port 8080

frontend:
	cd frontend && npm run dev

demo-on:
	USE_REAL_PROPka=false USE_REAL_DOCKING=false make backend

demo-propka:
	USE_REAL_PROPka=true USE_REAL_DOCKING=false make backend
