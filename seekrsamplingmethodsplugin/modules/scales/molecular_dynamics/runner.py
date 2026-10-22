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
import seekr.modules.scales.molecular_dynamics.runner as scales_md_runner

def is_finished_md(
        model: structures.Seekr_model,
        stage: md_structures.Stage_md,
        anchor_or_unpartitioned: engines_structures.MD_unpartitioned_settings | structures.Anchor,
        swarm_id: int = 0,
        completion_criteria: completion_criteria_structures\
            .Base_completion_criteria | None = None,
        cv_values: list[float] | None = None
    ) -> bool:
    """
    Check if the simulation is finished.
    """
    found, ending, _ = is_finished_md_with_message(
        model, stage, anchor_or_unpartitioned, swarm_id, 
        completion_criteria=completion_criteria, cv_values=cv_values)
    return found, ending

def is_finished_md_with_message(
        model: structures.Seekr_model,
        stage: md_structures.Stage_md,
        anchor_or_unpartitioned: engines_structures.MD_unpartitioned_settings \
            | structures.Anchor,
        swarm_id: int = 0,
        completion_criteria: completion_criteria_structures\
            .Base_completion_criteria | None = None,
        cv_values: list[float] | None = None,
        use_json: bool = False
    ) -> tuple[bool, bool, typing.Union[None, str, dict]]:
    """
    Check if the simulation is finished and return with message.
    """
    found = False
    ending = False
    message: typing.Union[None, str, dict] = None
    completion_criteria: completion_criteria_structures.Base_md_completion_criteria
    if completion_criteria is None:
        completion_criteria = stage.completion_criteria
    
    if completion_criteria.type == "cv_value_attained":

        cv = model.get_cv_by_name(completion_criteria.cv_name)
        found = True
        if cv_values is None:
            checkpoint_information_basename = engines_runner\
                .get_checkpoint_information_basename(swarm_id)
            stage_parameter = anchor_or_unpartitioned.stage_parameters[stage.name]
            stage_directory = stage_parameter["directory"]
            if stage.scope == "partitioned":
                assert anchor_or_unpartitioned is not None
                full_directory = os.path.join(
                    model.directory, "anchor_sets", stage.set_name, 
                    anchor_or_unpartitioned.directory, stage_directory)
            else:
                full_directory = os.path.join(
                    model.directory, anchor_or_unpartitioned.directory, stage_directory)
            checkpoint_information_filename = os.path.join(
                full_directory, checkpoint_information_basename)
            if not os.path.exists(checkpoint_information_filename):
                # Simulation hasn't run yet
                if use_json:
                    return True, False, {
                        "completion_criteria": completion_criteria.type,
                        "finished": False,
                        "progress": 0.0,
                        "cv_value_reached": None,
                        "cv_value_needed_for_completion": completion_criteria.cv_value
                    }
                else:
                    return True, False, "No simulation yet run - criteria not satisfied."
            checkpoint_information = engines_runner.load_checkpoint_information(
                checkpoint_information_filename)
            cv_values = checkpoint_information.cv_values
            
        cv_index = cv.index
        if cv_values is None:
            cv_value = 0.0
        else:
            cv_value = cv_values[cv_index]

        done_value = completion_criteria.cv_value
        #progress = max(0.0, min(cv_value / done_value, 1.0))
        if completion_criteria.cv_inequality == "greater_than":
            finished = cv_value > done_value
            if finished:
                progress = 1.0
            elif done_value > 0:
                progress = max(0.0, min(cv_value / done_value, 1.0))
            elif done_value == 0:
                progress = 0.0
            else:
                # Target is negative and cv_value is still at or below it.
                # cv/target is >= 1 on that side, so invert the ratio.
                progress = max(0.0, min(done_value / cv_value, 1.0))
        elif completion_criteria.cv_inequality == "less_than":
            finished = cv_value < done_value
            if finished:
                progress = 1.0
            elif done_value > 0:
                progress = max(0.0, min(done_value / cv_value, 1.0))
            elif done_value == 0:
                progress = 0.0
            else:
                # Target is negative and cv_value has not gone below it.
                if cv_value >= 0:
                    progress = 0.0
                else:
                    progress = max(0.0, min(cv_value / done_value, 1.0))

        else:
            raise ValueError(
                "cv_inequality must be 'greater_than' or 'less_than', "
                f"got {completion_criteria.cv_inequality!r}.")
                
        if use_json:
            if finished:
                message = {
                    "completion_criteria": completion_criteria.type,
                    "finished": True,
                    "progress": progress,
                    "cv_value_reached": cv_value,
                    "cv_value_needed_for_completion": completion_criteria.cv_value
                }
            else:
                message = {
                    "completion_criteria": completion_criteria.type,
                    "finished": False,
                    "progress": progress,
                    "cv_value_reached": cv_value,
                    "cv_value_needed_for_completion": completion_criteria.cv_value
                }
        else:
            if finished:
                message = f"CV value of {cv_value} attained. Required "\
                    f"value of {completion_criteria.cv_value}."
            else:
                message = f"System reached a CV value of {cv_value}. "\
                    f"Need to attain a value of {completion_criteria.cv_value} "\
                    f"for completion."
            
        return found, finished, message
    
    elif completion_criteria.type == "combined":
        found = True
        messages = []
        endings = []
        operator = completion_criteria.operator
        for completion_criterion in completion_criteria.completion_criteria:
            this_ending, this_message = scales_md_runner.is_finished_md_with_message(
                model, stage, anchor_or_unpartitioned, swarm_id, use_json,
                completion_criteria=completion_criterion, cv_values=cv_values)
            messages.append(this_message)
            endings.append(this_ending)
        
        if operator == "and":
            if all(endings):
                ending = True
            else:
                ending = False
        elif operator == "or":
            if any(endings):
                ending = True
            else:
                ending = False
        if use_json:
            total_progress = sum([message["progress"] for message in messages]) \
                / len(messages) if len(messages) > 0 else 0.0
            message = {
                "completion_criteria": completion_criteria.type,
                "finished": ending,
                "progress": total_progress,
                "messages": messages
            }
        else:
            message = "\n".join(messages)

    return found, ending, message