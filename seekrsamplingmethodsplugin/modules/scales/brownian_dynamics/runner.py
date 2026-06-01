"""
modules/scales/brownian_dynamics/runner.py

Evaluating whether bd_scale stages are complete.
"""

import seekr.modules.structures as structures
import seekr.modules.scales.structures as scales_structures

def is_finished_bd_with_message(
        model: structures.Seekr_model,
        stage_bd: scales_structures.Stage,
        use_json: bool = True
    ) -> tuple[bool, str | dict]:
    """

    """
    # Currently no completion criteria for Brownian Dynamics stages.
    return False, None