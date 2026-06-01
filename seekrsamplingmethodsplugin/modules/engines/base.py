from __future__ import annotations
"""
modules/engines/base.py

Provide base structures common to engines in the plugin.
"""

from attrs import define, field, validators, Factory

import seekr.modules.base as base

#@define
#class CUSTOM_unpartitioned_settings:
#    # To be used for a common set of parameters for all anchors
#    # within CD simulations
#    directory: str = field(default=base.UNPARTITIONED_DIRECTORY, 
#                           validator=validators.instance_of(str))
#    stage_parameters: dict[str, dict | None] = field(default=Factory(dict),)

#@define
#class CUSTOM_unpartitioned_settings_input:
#    # To be used for a common set of parameters for all anchors
#    # within CD simulations
#    directory: str = field(default=base.UNPARTITIONED_DIRECTORY, 
#                           validator=validators.instance_of(str))