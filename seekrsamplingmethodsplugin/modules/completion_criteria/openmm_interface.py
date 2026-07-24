"""
modules/completion_criteria/openmm_interface.py

OpenMM-specific functions for seekrsamplingmethodsplugin completion criteria.
"""

import typing
import logging

from attrs import define, field
import seekr.modules.structures as structures
import seekr.modules.completion_criteria.structures \
    as completion_criteria_structures
import seekr.modules.engines.structures as engines_structures
import seekr.modules.completion_criteria.openmm_interface \
    as compl_crit_openmm_interface
import seekr.modules.scales.molecular_dynamics.structures \
    as md_structures
import seekr.modules.reporters.openmm_interface as reporters_openmm_interface

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
        cv_values: list[float] = kwargs["cv_values"]
        swarm_id: int = kwargs["swarm_id"]
        return_message: bool
        if print_update:
            found, ending, message = plugin_md_runner\
                .is_finished_md_with_message(
                    self._model, self._stage, 
                    self._anchor_or_unpartitioned, 
                    swarm_id=swarm_id,
                    completion_criteria=self._completion_criteria,
                    cv_values=cv_values)
            logging.info(message)
        else:
            found, ending = plugin_md_runner.is_finished_md(
                self._model, self._stage, self._anchor_or_unpartitioned, 
                completion_criteria=self._completion_criteria,
                cv_values=cv_values, swarm_id=swarm_id)

        return ending
        
    def needs_cv_values(self) -> bool:
        """
        This completion criteria does need the CV values in order 
        to be evaluated.
        """
        return True

@define
class Combined_completion_criteria_openmm_interface(
        compl_crit_openmm_interface\
            .Base_completion_criteria_openmm_interface):
    """
    Interface for the combined completion criteria using
    the OpenMM engine.
    """
    _completion_criteria: plugin_compl_crit_struct\
        .Combined_completion_criteria
    _model: structures.Seekr_model
    _anchor_or_unpartitioned: engines_structures\
        .MD_unpartitioned_settings | structures.Anchor
    _stage: md_structures.Stage_md
    _reporter_interfaces: list[reporters_openmm_interface.Base_reporter_openmm_interface] | None = None
    _children_interfaces: list[compl_crit_openmm_interface\
        .Base_completion_criteria_openmm_interface] = field(factory=list)

    def _ensure_children(self):
        if self._children_interfaces:
            return
        for cc in self._completion_criteria.completion_criteria:
            self._children_interfaces.append(
                compl_crit_openmm_interface.assign_completion_interface(
                    cc, self._model, self._stage, self._anchor_or_unpartitioned, 
                    self._reporter_interfaces))
    
    def __attrs_post_init__(self):
        self._ensure_children()

    def criteria_satisfied(
            self,
            **kwargs: typing.Any
            ) -> bool:
        """
        Evaluate the current simulation state to see if the combined 
        completion criteria is satisfied.
        """
        assert len(self._completion_criteria.completion_criteria) > 0, \
            "Combined completion criteria must have at least one completion criterion."
        operator = self._completion_criteria.operator
        for i, completion_criterion in enumerate(
                self._completion_criteria.completion_criteria):
            completion_criterion_openmm_interface = self._children_interfaces[i]
            if operator == "and":
                if not completion_criterion_openmm_interface.criteria_satisfied(**kwargs):
                    return False
            elif operator == "or":
                if completion_criterion_openmm_interface.criteria_satisfied(**kwargs):
                    return True
        if operator == "and":
            return True
        elif operator == "or":
            return False
        
    def needs_cv_values(self) -> bool:
        """
        This completion criteria does need the CV values in order 
        to be evaluated.
        """
        for child_interface in self._children_interfaces:
            if child_interface.needs_cv_values():
                return True
        return False

def assign_completion_interface(
        completion_criteria: completion_criteria_structures\
            .Base_completion_criteria, 
        model: structures.Seekr_model,
        stage: md_structures.Stage_md,
        anchor_or_unpartitioned: engines_structures\
            .MD_unpartitioned_settings | structures.Anchor,
        reporter_interfaces: list | None = None
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
    elif completion_criteria.type == "combined":
        is_found = True
        openmm_interface \
            = Combined_completion_criteria_openmm_interface(
                completion_criteria, model, anchor_or_unpartitioned, 
                stage, reporter_interfaces)
    
    return is_found, openmm_interface