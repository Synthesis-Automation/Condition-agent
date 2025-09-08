install:
	python -m pip install -r requirements.txt
run:
	uvicorn app.main:app --reload --port 8000
test:
	pytest -q

# Resolve items from the registry
# Examples:
#   make registry Q="Toluene"
#   make registry Q="108-88-3 PhMe" JSONL=1
#   make registry FILE=queries.txt PRETTY=1
#   make registry DEMO=1 FIELD=uid
registry:
ifeq (,$(strip $(REG_ARGS)))
REG_ARGS :=
endif
ifdef Q
REG_ARGS += $(Q)
endif
ifdef FILE
REG_ARGS += -f $(FILE)
endif
ifdef DEMO
REG_ARGS += --demo
endif
ifdef FIELD
REG_ARGS += --field $(FIELD)
endif
ifdef JSONL
REG_ARGS += --jsonl
endif
ifdef PRETTY
REG_ARGS += --pretty
endif
	python -m chemtools.cli.registry $(REG_ARGS)

# Build DRFP NPZ index for fast warm start
drfp-index:
	python scripts/precompute_drfp.py --out artifacts/drfp_index.npz

# Build higher-fidelity 4096-bit index across data/reaction_dataset
drfp-index-4096:
	python scripts/precompute_drfp.py --source data/reaction_dataset --out artifacts/all_drfp_4096.npz --n-bits 4096 --radius 3

# Run API with DRFP index preloaded (prefers 4096-bit if present)
run-drfp:
	@if [ -f artifacts/all_drfp_4096.npz ]; then \
	  CHEMTOOLS_DRFPPATH=artifacts/all_drfp_4096.npz uvicorn app.main:app --reload --port 8000; \
	elif [ -f artifacts/all_drfp.npz ]; then \
	  CHEMTOOLS_DRFPPATH=artifacts/all_drfp.npz uvicorn app.main:app --reload --port 8000; \
	else \
	  echo "No DRFP NPZ found. Build with 'make drfp-index-4096' or 'make drfp-index' first."; exit 1; \
	fi

# Launch Gradio UI (auto-loads dataset if present)
ui:
	python scripts/ui_gradio.py
