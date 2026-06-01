"""
test_analyze.py

Smoke tests for seekr.analyze — verify that analyze() returns without
exception after all non-BD stages have been run.
"""

import seekr.analyze as analyze

def test_analyze_returns_result(seekr_run):
    """analyze() returns a non-None result after simulations."""
    model, _, scenario_name = seekr_run
    # Only analyze models that completed the MMVT stage
    mmvt_stages = [s for s in model.stages if s.name == "MMVT"]
    if not mmvt_stages:
        return  # no MMVT stage in this scenario; skip silently
    if "hg_scenario" in scenario_name:
        result = analyze.analyze(model, skip_checks=True,
                                 only_extract_data=True)
    else:
        result = analyze.analyze(model, skip_checks=True)
    assert result is not None, (
        f"[{scenario_name}] analyze() returned None"
    )
