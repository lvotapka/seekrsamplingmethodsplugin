from __future__ import annotations
"""
modules/engines/preparer.py

Import input classes and functions from the different engines
for easy access here.
"""
import os
import typing
from attrs import define, field, validators, Factory

import seekr.modules.base as base

#import seekrsamplingmethodsplugin.modules.engines.custom_engine.preparer as custom_engine_preparer
            
#engine_inputs: typing.Dict[str, typing.Dict[str, typing.Any]] \
#    = {"custom": custom_engine_preparer.Custom_engine_settings_input, }

engine_inputs: typing.Dict[str, typing.Dict[str, typing.Any]] \
    = {}