from __future__ import annotations
"""
modules/engines/runner.py

Functions usable by plugin engine(s) for running.
"""

import seekr.modules.structures as structures

def run_simulation(
        model: structures.Seekr_model, 
        stage_index: int, 
        anchor_index: int | None = None, 
        device_index: str | None = None,
        force_overwrite: bool = False
        ) -> None:
    """
    Determine which md engine to run and run it
    """
    #stage = model.stages[stage_index]
    #engine = model.get_engine_by_scale_type(stage.scale_type)
    #if engine.type == "custom_engine":
    #    import seekrsamplingmethodsplugin.modules.engines.custom_engine.runner \
    #        as runner_custom_engine
    #    runner_custom_engine.run_simulation(
    #        model, stage_index, anchor_index, 
    #        device_index, force_overwrite)
    #else:
    #    pass

    return