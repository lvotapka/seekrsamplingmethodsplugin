from __future__ import annotations
"""
modules/scales/logistic/runner.py

Perform runner functions specific to logistic stages.
"""

import cattrs
import logging
import importlib
import importlib.util

import seekr.modules.base as base
import seekr.modules.structures as structures
import seekr.modules.scales.logistic.base as logistic_base

import seekrsamplingmethodsplugin.modules.scales.logistic\
    .swarm_structure_from_reporter_logistic_stage \
        as swarm_structure_from_reporter_logistic_stage
import seekrsamplingmethodsplugin.modules.scales.logistic\
    .egress_pathway_logistic_stage \
        as egress_pathway_logistic_stage

def is_finished_logistic_with_message(
        model: structures.Seekr_model,
        stage: logistic_base.Stage_logistic,
        use_json: bool = True
    ) -> tuple[bool, str | dict]:
    """

    """
    found = False
    completed = False
    if use_json:
        message = {}
    else:
        message = ""

    if stage.type == "swarm_structures_from_reporter_stage":
        completed, message \
            = swarm_structure_from_reporter_logistic_stage\
                .is_finished_swarm_structures_from_reporter_logistic_with_message(
                    model, stage, use_json=use_json)
        found = True
    
    elif stage.type == "egress_pathway_stage":
        completed, message \
            = egress_pathway_logistic_stage\
                .is_finished_egress_pathway_logistic_with_message(
                    model, stage, use_json=use_json)
        found = True
    
    return found, completed, message

def choose_next_simulation_logistic(
        model: structures.Seekr_model, 
        stage: logistic_base.Stage_logistic,
        force_overwrite: bool,
        use_json: bool = True,
        silent: bool = False
        ) -> list:
    found = False
    instructions = None
    if stage.type == "swarm_structures_from_reporter_stage":
        instructions = swarm_structure_from_reporter_logistic_stage\
            .choose_next_simulation_swarm_structures_from_reporter_logistic(
                model, stage, force_overwrite, use_json, silent=silent)
        found = True
    elif stage.type == "egress_pathway_stage":
        instructions = egress_pathway_logistic_stage\
            .choose_next_simulation_egress_pathway_logistic(
                model, stage, force_overwrite, use_json, silent=silent)
        found = True
    
    assert instructions is not None, "No instructions found for stage type: " + stage.type
    if not silent:
        stage.write_to_log(model, instructions)
    return found, instructions
    
def do_run_instruction_logistic(
        model: structures.Seekr_model,
        force_overwrite: bool,
        stage_index: int,
        instruction: list
        ) -> None:
    found = False
    stage = model.stages[stage_index]
    if stage.type == "swarm_structures_from_reporter_stage":
        swarm_structure_from_reporter_logistic_stage\
            .do_run_instruction_swarm_structures_from_reporter_logistic(
                model, force_overwrite, stage_index=stage_index,
                instruction=instruction)
        found = True
    
    elif stage.type == "egress_pathway_stage":
        egress_pathway_logistic_stage\
            .do_run_instruction_egress_pathway_logistic(
                model, force_overwrite, stage_index=stage_index,
                instruction=instruction)
        found = True
    
    return found