"""
modules/completion_criteria/openmm_interface.py

OpenMM-specific functions for seekrsamplingmethodsplugin completion criteria.
"""

import typing
import logging

from attrs import define
import seekr.modules.structures as structures
import seekr.modules.completion_criteria.structures \
    as completion_criteria_structures
import seekr.modules.engines.structures as engines_structures
import seekr.modules.completion_criteria.openmm_interface \
    as compl_crit_openmm_interface
import seekr.modules.scales.molecular_dynamics.structures \
    as md_structures

import seekrsamplingmethodsplugin.modules.completion_criteria\
    .structures as plugin_compl_crit_struct
import seekrsamplingmethodsplugin.modules.scales.molecular_dynamics\
    .runner as plugin_md_runner

@define
class CV_value_attained_completion_criteria_openmm_interface(
        compl_crit_openmm_interface\
            .Base_completion_criteria_openmm_interface):
    """
    Interface for the CV-value-attained completion criteria using
    the OpenMM engine.
    """
    _completion_criteria: plugin_compl_crit_struct\
        .CV_value_attained_completion_criteria
    _model: structures.Seekr_model
    _anchor_or_unpartitioned: engines_structures\
        .MD_unpartitioned_settings | structures.Anchor
    _stage: md_structures.Stage_md

    def criteria_satisfied(
            self,
            **kwargs: typing.Any
            ) -> bool:
        # simulation, self._cv_force
        """
        Evaluate the current simulation state to see if the CV is 
        attained.
        """
        print_update: bool = kwargs["print_update"]
        current_step: int = kwargs["current_step"]
        cv_values: list[float] = kwargs["cv_values"]
        swarm_id: int = kwargs["swarm_id"]
        return_message: bool
        if print_update:
            found, ending, message = plugin_md_runner\
                .is_finished_md_with_message(
                    self._model, self._stage, 
                    self._anchor_or_unpartitioned, 
                    current_step, cv_values=cv_values, 
                    swarm_id=swarm_id)
            logging.info(message)
        else:
            found, ending = plugin_md_runner.is_finished_md(
                self._model, self._stage, self._anchor_or_unpartitioned, 
                current_step, cv_values=cv_values, swarm_id=swarm_id)

        return ending
        
    def needs_cv_values(self) -> bool:
        """
        This completion criteria does need the CV values in order 
        to be evaluated.
        """
        return True

def assign_completion_interface(
        completion_criteria: completion_criteria_structures\
            .Base_completion_criteria, 
        model: structures.Seekr_model,
        stage: md_structures.Stage_md,
        anchor_or_unpartitioned: engines_structures\
            .MD_unpartitioned_settings | structures.Anchor
        )-> tuple[bool, typing.Optional[
            compl_crit_openmm_interface\
                .Base_completion_criteria_openmm_interface]]:
    """
    Given this plugin's completion criteria object, assign the 
    appropriate OpenMM interface.
    """
    is_found = False
    openmm_interface = None
    # Detect if it's a completion criteria defined in this plugin
    if completion_criteria.type == "cv_value_attained":
        is_found = True
        openmm_interface \
            = CV_value_attained_completion_criteria_openmm_interface(
                completion_criteria, model, anchor_or_unpartitioned, 
                stage,)
    
    return is_found, openmm_interface