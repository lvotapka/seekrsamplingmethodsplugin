"""
modules/reporters/conditions.py

Condition definitions for Conditional_position_reporter.
"""

from __future__ import annotations

import typing

from attrs import define, field, validators


Inequalities = typing.Literal[
    "greater_than",
    "less_than",
    "greater_than_or_equal_to",
    "less_than_or_equal_to",
]


@define
class CV_threshold_condition:
    """
    True when a named CV satisfies an inequality against a threshold.
    """
    cv_name: str = field(validator=validators.instance_of(str))
    cv_value: float = field(validator=validators.instance_of(float))
    type: typing.Literal["cv_threshold"] = "cv_threshold"
    cv_inequality: Inequalities = field(
        default="greater_than",
        validator=validators.in_([
            "greater_than",
            "less_than",
            "greater_than_or_equal_to",
            "less_than_or_equal_to",
        ]))

    def is_satisfied(self, cv_value: float) -> bool:
        if self.cv_inequality == "greater_than":
            return cv_value > self.cv_value
        if self.cv_inequality == "less_than":
            return cv_value < self.cv_value
        if self.cv_inequality == "greater_than_or_equal_to":
            return cv_value >= self.cv_value
        if self.cv_inequality == "less_than_or_equal_to":
            return cv_value <= self.cv_value
        raise ValueError(f"Unknown inequality: {self.cv_inequality}")


conditions: typing.Dict[str, typing.Type[CV_threshold_condition]] = {
    "cv_threshold": CV_threshold_condition,
}


def structure_condition(obj: typing.Any) -> CV_threshold_condition:
    """
    Convert a dict or condition object into a CV_threshold_condition.
    """
    if isinstance(obj, CV_threshold_condition):
        return obj
    if isinstance(obj, dict):
        typ = obj.get("type", "cv_threshold")
        assert typ in conditions, f"Unknown condition type: {typ}"
        kwargs = {k: v for k, v in obj.items() if k != "type"}
        return conditions[typ](**kwargs)
    raise TypeError(f"Cannot structure condition from {type(obj)}: {obj!r}")


def structure_condition_list(
        objs: typing.Sequence[typing.Any] | None,
        ) -> list[CV_threshold_condition]:
    """
    Convert a list of condition dicts/objects into condition instances.
    """
    if objs is None:
        return []
    return [structure_condition(obj) for obj in objs]
