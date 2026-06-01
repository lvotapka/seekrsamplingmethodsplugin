"""
modules/anchor_schemes/preparer.py

In seekr, particular regions of configuration space are called 'anchors'
(a holdover from milestoning theory). Code in this module is used to
define the partitioning of these anchors throughout configuration 
space. Placement of the anchors allow for the definition of 'milestone',
which, in general, are surfaces in the configuration space. If a 
sampling method provided in this plugin requires dividing up configuration
space somehow, anchors and milestones should be determined in this module
for inclusion within the model.
"""

from attrs import define

# other imports here...

# @define
# class Anchor_scheme_class(anchor_scheme_preparer.Base_anchor_scheme):
#     """
#     This class defines a custom 'anchor_scheme' - a way to divide up
#     configuration space, as needed for an enhanced sampling method.
#     A 'make_anchors()' method must be devised here which is called by
#     seekr's modules/preparer.py
#     """
#     def make_anchors(
#         self,
#         ...):

# A new anchor scheme must be named and defined in the dictionary below
#anchor_schemes = {
#    "YOUR_ANCHOR_SCHEME_NAME": Anchor_scheme_class
#}

# Leaving empty for now...
anchor_schemes = {}
