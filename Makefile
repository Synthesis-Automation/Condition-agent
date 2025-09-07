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
