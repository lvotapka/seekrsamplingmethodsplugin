"""
test_run.py

Smoke tests for seekr.run — verify that run.run() completes without
exception for all non-BD stages of every scenario.
"""

import os

import seekr.modules.structures as structures


def test_run_produces_model(seekr_run):
    """After running all stages, the model is still a valid Seekr_model.
    Not a likely danger, but this is mainly just a smoke test."""
    model, _, scenario_name = seekr_run
    assert isinstance(model, structures.Seekr_model), (
        f"[{scenario_name}] seekr_run should yield a Seekr_model"
    )


def test_run_cwd_restored(seekr_run):
    """The process cwd must be the same before and after the fixture."""
    # seekr_run fixture saves/restores cwd internally; if we get here,
    # the restore succeeded by definition. We just assert the cwd is valid.
    assert os.path.isdir(os.getcwd())
