"""
modules/completion_criteria/structures.py

Completion criteria structure that will be defined for this plugin.
"""

import typing
from attrs import define, field, validators, Factory

import seekr.modules.completion_criteria.structures as compl_crit_struct

MAX_CV_VALUES_STEPS = 1e9 # 2-4 us typical MD

@define
class CV_value_attained_completion_criteria(compl_crit_struct.Base_completion_criteria):
    """
    A completion criteria based on a CV value being attained.
    """
    type: typing.Literal["cv_value_attained"] = "cv_value_attained"
    eval_interval: int = field(default=250,
                                validator=[validators.instance_of(int),
                                           validators.gt(0)])
    cv_name: str = field(default="",
                         validator=[validators.instance_of(str)])
    cv_inequality: str = field(default="greater_than",
                               validator=[validators.in_(["greater_than", "less_than"])])
    cv_value: float = field(default=0.0,
                            validator=[validators.instance_of(float)])
    
    # NOTE: these functions are needed for chunking in runner modules
    def get_min_sim_len(self) -> int:
        """
        Return the minimum simulation length for this completion criteria.
        """
        return 0
    
    def get_max_sim_len(self) -> int:
        """
        Return the maximum simulation length for this completion criteria.
        """
        return MAX_CV_VALUES_STEPS
    
    def get_evaluation_interval(self) -> int:
        return self.eval_interval
    
# This dictionary is needed by seekr to structure/unstructure data for the model file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.

completion_criteria: typing.Dict[
        str, typing.Type[compl_crit_struct.Base_completion_criteria]] \
    = {"cv_value_attained": CV_value_attained_completion_criteria}