from __future__ import annotations
"""
modules/scales/structures.py

Structures for scales within seekrsamplingmethodsplugin - handle 
structuring of the scale structures for this plugin.
"""

import typing

import seekrsamplingmethodsplugin.modules.scales.logistic\
    .swarm_structure_from_reporter_logistic_stage \
        as swarm_structure_from_reporter_logistic_stage
import seekrsamplingmethodsplugin.modules.scales.logistic\
    .structures as logistic_structures

# This dictionary is needed by seekr to structure/unstructure data for the model file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.
#scales: typing.Dict[str, typing.Any] \
#    = {"custom_scale": custom_scale_structures.Scale_custom, }

scales: typing.Dict[str, typing.Any] \
    = {}
stages: typing.Dict[str, typing.Any] \
    = {"swarm_structures_from_reporter_stage": \
        swarm_structure_from_reporter_logistic_stage\
            .Stage_swarm_structures_from_reporter,}

def register_structure_hooks(
        converter: typing.Any
        ) -> None:
    """
    Register the structure hooks for all scales and their sampling methods.
    """
    logistic_structures.register_structure_hooks(converter)
    return