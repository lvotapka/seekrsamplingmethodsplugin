"""
conftest.py

Shared pytest fixtures.

Parametrized pipeline fixtures
-------------------------------
  seekr_prepared  -- runs prepare.prepare() for each scenario
  seekr_run       -- runs all non-BD stages for each scenario


"""

import os
import copy

import pytest

import seekr.modules.structures as structures
import seekr.prepare as prepare
import seekr.run as run
import seekr.status as status

import seekrsamplingmethodsplugin.tests.create_seekr_input as create_seekr_input

# ---------------------------------------------------------------------------
# Scenario registry
# ---------------------------------------------------------------------------

SCENARIO_CONFIGS: dict = {
    # toy scenarios
    "toy_scenario":         dict(toy=True),
    # host-guest scenarios
    "hg_scenario":          dict(),
    
}

SCENARIO_PARAMS = [
    pytest.param("toy_scenario"),
    pytest.param("hg_scenario"),
]

# ---------------------------------------------------------------------------
# Custom CLI option
# ---------------------------------------------------------------------------

def pytest_addoption(parser):
    parser.addoption(
        "--full", action="store_true", default=False,
        help="Run the full test suite including slow scenarios "
             "(default: only fast scenarios)."
    )


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--full"):
        skip_slow = pytest.mark.skip(reason="slow scenario; use --full to run")
        for item in items:
            if "slow" in [m.name for m in item.iter_markers()]:
                item.add_marker(skip_slow)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _create_seekr_input(scenario_name: str, root_dir: str):
    cfg = SCENARIO_CONFIGS[scenario_name].copy()
    if cfg.pop("toy", False):
        return create_seekr_input.create_toy_input(root_dir, **cfg)
    return create_seekr_input.create_host_guest_input(root_dir, **cfg)


def _prepare_model(seekr_input):
    old_dir = os.getcwd()
    try:
        model, model_xml_path = prepare.prepare(
            seekr_input, force_overwrite=False)
        model.directory = os.path.abspath(os.path.dirname(model_xml_path))
    finally:
        os.chdir(old_dir)
    return model, model_xml_path


# ---------------------------------------------------------------------------
# Parametrized pipeline fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session", params=SCENARIO_PARAMS)
def seekr_prepared(request, tmp_path_factory):
    """Prepare a SEEKR model for a given scenario. Returns (model, xml_path, name)."""
    scenario_name = request.param
    root_dir = str(tmp_path_factory.mktemp(scenario_name))
    seekr_input = _create_seekr_input(scenario_name, root_dir)
    model, model_xml_path = _prepare_model(seekr_input)
    return model, model_xml_path, scenario_name


@pytest.fixture(scope="session")
def seekr_run(seekr_prepared):
    """
    Run all non-brownian-dynamics stages sequentially.
    Calls status() after each stage as a smoke test.
    Returns the same (model, xml_path, name) tuple.
    """
    model, model_xml_path, scenario_name = seekr_prepared
    structures.load_plugins(model.plugins)
    old_dir = os.getcwd()
    try:
        for stage in model.stages:
            if stage.scale_type == "brownian_dynamics":
                continue
            run.run(model, stage.name)
            status.status(model, "progress", print_json=False)
    finally:
        os.chdir(old_dir)
    return model, model_xml_path, scenario_name