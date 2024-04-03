# Define Python interpreter (adjust if needed)
PYTHON := python3

# Target to install all dependencies
all: install

# Target to install dependencies using pip
install:
	@$(PYTHON) -m pip install -r requirements.txt

# Phony target for cleaning (no actual file to clean)
.PHONY: clean

# Target to clean installed dependencies
clean:
	@$(PYTHON) -m pip uninstall -r requirements.txt