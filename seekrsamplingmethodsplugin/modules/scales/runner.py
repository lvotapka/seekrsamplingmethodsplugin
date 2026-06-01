from __future__ import annotations
"""
modules/scales/runner.py

Differentiate between different scales in this plugin and call 
the appropriate runner functions for each scale.
"""

import seekr.modules.structures as structures
import seekr.modules.scales.base as scales_base
import seekr.modules.engines.structures as engines_structures
import seekr.modules.scales.structures as scales_structures

#import seekrsamplingmethodsplugin.modules.engines.base as plugin_engines_base
#import seekrsamplingmethodsplugin.modules.scales.custom_scale.runner as custom_scale_runner


def is_finished_with_message(
        model: structures.Seekr_model,
        stage: scales_base.Base_stage,
        anchor_index: str | None = None,
        use_json: bool = True
    ) -> tuple[bool, bool, list[str]]:
    """
    Check if the given stage is finished, and return a tuple of (is_finished,
    message).
    """
    if use_json:
        final_message = {}
    else:
        final_message = ""
    found = False
    finished = False
    #if stage.scale_type == "custom_scale":
    #    found = True
    #    custom_engine = model.get_engine_by_scale_type("custom_scale")
    #    anchor_or_unpartitioned_index_list = custom_scale_runner\
    #        .construct_anchor_or_unpartitioned_index_list(model, stage, anchor_index)
    #    anchor_or_unpartitioned: plugin_engines_base.CD_unpartitioned_settings \
    #        | structures.Anchor
    #    if len(anchor_or_unpartitioned_index_list) == 0:
    #        finished = False
    #    else:
    #        finished = True
    #    for anchor_or_unpartitioned_index in anchor_or_unpartitioned_index_list:
    #        if anchor_or_unpartitioned_index is None:
    #            assert custom_engine.unpartitioned_settings is not None
    #            anchor_or_unpartitioned = custom_engine.unpartitioned_settings
    #        else:
    #            anchor_or_unpartitioned \
    #                = model.anchors[anchor_or_unpartitioned_index]
    #        tmp_finished, message = custom_scale_runner.is_finished_custom_with_message(
    #            model, stage, anchor_or_unpartitioned,
    #            use_json=use_json)
    #        if not tmp_finished:
    #            finished = False
    #        if use_json:
    #            final_message[anchor_or_unpartitioned_index] = message
    #        else:
    #            if anchor_or_unpartitioned_index is None:
    #                final_message += f"    Unpartitioned: {message}\n"
    #            else:
    #                final_message += f"    Anchor {anchor_or_unpartitioned_index}: {message}\n"

    return found, finished, final_message

def construct_anchor_or_unpartitioned_index_list(
        model: structures.Seekr_model,
        stage: scales_base.Base_stage,
        anchor_index: str
        ) -> list[int] | None:
    """
    Construct a list of anchor or unpartitioned indices based on the stage and
    anchor index.
    """
    found = False
    result = None
    #if stage.scale_type == "custom_scale":
    #    result = cd_runner.construct_anchor_or_unpartitioned_index_list(
    #        model, stage, anchor_index)
    #    found = True
        
    return found, result

def choose_next_simulation(
        model: structures.Seekr_model,
        stage: scales_base.Base_stage,
        choose_params: dict,
        force_overwrite: bool,
        swarm_id: int | None = None,
        use_json: bool = True,
        silent: bool = False,
        benchmark: bool = False
        ) -> list:
    """
    Choose the next simulation to run based on the stage scale, and return
    a list of instructions to run in this round.
    """
    found = False
    run_instructions = []
    #if stage.scale_type == "custom_scale":
    #    anchor_index_instruction = choose_params["anchor_index_instruction"]
    #    device_index = choose_params["device_index"]
    #    run_instructions = custom_scale_runner.choose_next_simulation_cd(
    #        model, stage, anchor_index_instruction, force_overwrite, 
    #        device_index, use_json, silent=silent, benchmark=benchmark)
    #    found = True

    return found, run_instructions

def get_calculation_directory(
        engine: engines_structures.Base_engine,
        root_directory: str):
    """
    Return the calculation directory for a given custom scale stage index.
    """
    found = False
    full_directory = ""
    #if engine.scale_type == "custom_scale":
    #    full_directory = custom_scale_runner.get_calculation_directory_cd(
    #        engine, root_directory)
    #    found = True

    return found, full_directory

def do_run_instruction(
        model: structures.Seekr_model,
        force_overwrite: bool,
        stage_index: int,
        instruction: list,
        device_index: str | None,
        benchmark: bool = False
    ) -> bool:
    """
    Execute the run instructions for a given scale type.
    """
    found = False
    stage = model.stages[stage_index]
    #if stage.scale_type == "custom_scale":
    #    custom_scale_runner.do_run_instruction_cd(
    #        model, force_overwrite, stage_index=stage_index, 
    #        instruction=instruction, device_index=device_index,
    #        benchmark=benchmark)
    #    found = True
    return found
