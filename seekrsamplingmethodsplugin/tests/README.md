# Running Tests

## Basic Commands

```bash
# Run all tests (excludes slow tests by default)
pytest

# Run all tests including slow scenarios
pytest --full

# Run a specific test file
pytest test_prepare.py

# Run a specific test function
pytest test_prepare.py::test_function_name
```

## Test times

Show the slowest test/setup timings at the end of the run:

```bash
pytest seekr/tests -vv --full --durations=0 --durations-min=0
```

## Available Marks

```bash
# Run only tests that need CUDA
pytest -m needs_cuda

# Run only tests that need OpenMM
pytest -m needs_openmm

# Run slow tests explicitly
pytest -m slow

# Exclude specific marks
pytest -m "not needs_cuda"

# Combine marks with boolean expressions
pytest -m "needs_openmm and not slow"
```

## Coverage

Run from the `seekrsamplingmethodsplugin/tests/` directory:

```bash
# Generate HTML report, then open it
pytest --full --cov=seekrsamplingmethodsplugin --cov-report=html:htmlcov
xdg-open htmlcov/index.html # Linux

# Terminal report
pytest --full --cov=seekrsamplingmethodsplugin --cov-report=term-missing
```

## Marks Defined

- `needs_cuda` - Tests requiring CUDA
- `needs_openmm` - Tests requiring OpenMM
- `slow` - Slow tests (skipped by default, run with `--full`)
