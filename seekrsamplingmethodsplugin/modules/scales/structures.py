from __future__ import annotations
"""
modules/scales/structures.py

Structures for scales within seekrsamplingmethodsplugin - handle structuring
of the scale structures for this plugin.
"""

import typing

#import seekrsamplingmethodsplugin.modules.scales.custom_scale.structures \
#    as custom_scale_structures

# This dictionary is needed by seekr to structure/unstructure data for the model file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.
#scales: typing.Dict[str, typing.Any] \
#    = {"custom_scale": custom_scale_structures.Scale_custom, }
#stages: typing.Dict[str, typing.Any] \
#    = {"custom_stage": custom_scale_structures.Stage_custom, }

scales: typing.Dict[str, typing.Any] \
    = {}
stages: typing.Dict[str, typing.Any] \
    = {}