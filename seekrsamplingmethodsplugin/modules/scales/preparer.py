from __future__ import annotations
"""
modules/scales/preparer.py

Objects for the prepare stage.
"""

import typing

import seekr.modules.structures as structures

import seekrsamplingmethodsplugin.modules.scales.logistic\
    .swarm_structure_from_reporter_logistic_stage \
        as swarm_structure_from_reporter_logistic_stage
import seekrsamplingmethodsplugin.modules.scales.logistic\
    .egress_pathway_logistic_stage \
        as egress_pathway_logistic_stage
# This dictionary is needed to structure/unstructure prepare stage objects.

def generate_files(
        model: structures.Seekr_model, 
        directory: str,
        anchor_set_name: str,
        scale_type: str
        ) -> bool:
    """
    Determine the scale engines and prepare the necessary files for them.
    Returns a bool based on whether this plugin contained the desired scale.
    """
    #if scale_type == "custom_scale":
    #    # Do something to generate starting files for this scale?
    #    return True
    #else:
    #    return False
    return False

#stage_inputs: typing.Dict[str, typing.Any] \
#    = {"custom_stage": custom_scale_preparer.Input_stage_custom, }
stage_inputs: typing.Dict[str, typing.Any] \
    = {"swarm_structures_from_reporter_stage": \
        swarm_structure_from_reporter_logistic_stage\
            .Input_stage_swarm_structures_from_reporter, 
        "egress_pathway_stage": \
            egress_pathway_logistic_stage.Input_stage_egress_pathway, }

def register_structure_hooks(
        converter: typing.Any
        ) -> None:
    """
    Register the structure hooks for all scales and their sampling methods.
    """
    #custom_scale_preparer.register_structure_hooks(converter)
    return
    