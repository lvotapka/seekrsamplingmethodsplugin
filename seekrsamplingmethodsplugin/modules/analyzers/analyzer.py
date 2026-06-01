"""
modules/analyzers/analyzer.py

The general hub where the correct analyzer is selected based on the model type,
the analysis is performed, and results are returned in the form of an analysis
object for this plugin.

Typically, a new sampling method could have an analyzer supplied for it to allow
some sort of analysis of the sampling method's output. For instance, a sampling
method may allow one to perform free energy analyses or residence time 
estimations. Such analysis would be implemented in scripts in this directory
that could be imported and called in this script.

Alternatively, this program could perform additional analyses on top of the
analyses already bundled in the main seekr program. For example, if SEEKR
is used to run the string method, a program could be written here that would
allow visualization and quantification of those results.
"""

# This list must contain any strings named with the 'type' attribute of a 
# sampling method to be analyzed.
AVAILABLE_ANALYZERS = []

import typing

import seekr.modules.structures as structures
#  Place other imports here
# import seekrsamplingmethodsplugin.modules.analyzers.SOMETHING_analyzer as plugin_analyzers_SOMETHING_analyzer

def analyze(
        model: structures.Seekr_model,
        stage_index: int | None = None,
        skip_checks: bool = False,
        min_time_ps: float = 0.0,
        max_time_ps: float | None = None,
        error_analysis_method: str = "bootstrap",
        only_extract_data: bool = False
    ) -> tuple[bool, typing.Any]: #| plugin_analyzers_SOMETHING_analyzer.SOMETHING_analyzer
    """
    A call to analyze.py in seekr will automatically go to here if the plugin is
    included in the model. If the used sampling method is found in
    AVAILABLE_ANALYZERS, then the analyses put here will be run.
    """
    stage = model.stages[stage_index-1]
    is_found = False
    analysis = None
    #if stage.sampling_method.type == "SOMETHING":
    #    is_found = True
    #    analysis = plugin_analyzers_SOMETHING_analyzer.analyze_SOMETHING(...)

    return is_found, analysis
