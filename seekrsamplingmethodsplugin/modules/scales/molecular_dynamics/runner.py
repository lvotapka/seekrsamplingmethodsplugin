"""
modules/scales/molecular_dynamics/runner.py

Evaluating whether md_scale stages are complete.
"""

import os
import typing

import seekr.modules.structures as structures
import seekr.modules.engines.structures as engines_structures
import seekr.modules.engines.runner as engines_runner
import seekr.modules.completion_criteria.structures \
    as completion_criteria_structures
import seekr.modules.scales.molecular_dynamics.structures as md_structures

def is_finished_md(
        model: structures.Seekr_model,
        stage: md_structures.Stage_md,
        anchor_or_unpartitioned: engines_structures.MD_unpartitioned_settings | structures.Anchor,
        current_step: int,
        swarm_id: int = 0,
        cv_values: list[float] | None = None
    ) -> bool:
    """
    Check if the simulation is finished.
    """
    found, ending, _ = is_finished_md_with_message(
        model, stage, anchor_or_unpartitioned, current_step, swarm_id, cv_values)
    return found, ending

def is_finished_md_with_message(
        model: structures.Seekr_model,
        stage: md_structures.Stage_md,
        anchor_or_unpartitioned: engines_structures.MD_unpartitioned_settings \
            | structures.Anchor,
        current_step: int,
        swarm_id: int = 0,
        cv_values: list[float] | None = None,
        use_json: bool = False
    ) -> tuple[bool, bool, typing.Union[None, str, dict]]:
    """
    Check if the simulation is finished and return with message.
    """
    found = False
    ending = False
    message: typing.Union[None, str, dict] = None
    completion_criteria: completion_criteria_structures.Base_md_completion_criteria = \
        stage.completion_criteria
    if completion_criteria.type == "cv_value_attained":
        cv = model.get_cv_by_name(completion_criteria.cv_name)
        if cv_values is None:
            checkpoint_information_basename = engines_runner\
                .get_checkpoint_information_basename(swarm_id)
            stage_parameter = anchor_or_unpartitioned.stage_parameters[stage.name]
            stage_directory = stage_parameter["directory"]
            full_directory = os.path.join(
                model.directory, anchor_or_unpartitioned.directory, stage_directory)
            checkpoint_information_filename = os.path.join(
                full_directory, checkpoint_information_basename)
            if not os.path.exists(checkpoint_information_filename):
                # Simulation hasn't run yet
                return True, False, "No simulation yet run - criteria not satisfied."
            checkpoint_information = engines_runner.load_checkpoint_information(
                checkpoint_information_filename)
            cv_values = checkpoint_information.cv_values
            assert not cv_values is None, \
                "There is a problem - previous restart had no CV values provided"
            
        cv_index = cv.index
        cv_value = cv_values[cv_index]
        if use_json:
            incomplete_message = {
                "completion_criteria": completion_criteria.type,
                "attained": False,
                "cv_value_reached": cv_value,
                "cv_value_needed_for_completion": completion_criteria.cv_value
            }
            complete_message = {
                "completion_criteria": completion_criteria.type,
                "attained": True,
                "cv_value_reached": cv_value,
                "cv_value_needed_for_completion": completion_criteria.cv_value
            }
        else:
            incomplete_message = f"System reached a CV value of {cv_value}. "\
                f"Need to attain a value of {completion_criteria.cv_value} "\
                f"for completion."
            complete_message = f"CV value of {cv_value} attained. Required "\
            f"value of {completion_criteria.cv_value}."
        if completion_criteria.cv_inequality == "greater_than":
            if cv_value > completion_criteria.cv_value:
                return True, True, complete_message
            else:
                return True, False, incomplete_message
        elif completion_criteria.cv_inequality == "less_than":
            if cv_value < completion_criteria.cv_value:
                return True, True, complete_message
            else:
                return True, False, incomplete_message
    return found, ending, message