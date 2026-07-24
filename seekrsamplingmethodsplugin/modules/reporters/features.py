"""
modules/reporters/features.py

Geometric feature definitions for the Feature_reporter.
"""

from __future__ import annotations

import typing

import numpy as np
from attrs import define, field, validators


class Base_feature(typing.Protocol):
    """
    Base protocol for all Feature_reporter features.
    """
    type: typing.Literal["base"] = "base"
    name: str | None

    def column_name(self) -> str:
        """
        Return the CSV column header for this feature.
        """
        raise NotImplementedError("Subclasses must implement this method.")

    def compute(
            self,
            positions: np.ndarray,
            box_vectors: np.ndarray | None = None,
            masses: np.ndarray | None = None,
            ) -> float:
        """
        Compute the feature value from atomic positions (nm).
        """
        raise NotImplementedError("Subclasses must implement this method.")


def _as_index_list(indices: int | list[int]) -> list[int]:
    if isinstance(indices, int):
        return [indices]
    assert len(indices) > 0, "Atom index list must be non-empty."
    return list(indices)


def _format_index_group(indices: int | list[int]) -> str:
    if isinstance(indices, int):
        return str(indices)
    return "-".join(str(i) for i in indices)


def _com(
        positions: np.ndarray,
        indices: list[int],
        masses: np.ndarray | None,
        ) -> np.ndarray:
    coords = positions[indices]
    if len(indices) == 1:
        return coords[0].copy()
    assert masses is not None, \
        "Masses are required to compute a center of mass for multi-atom groups."
    group_masses = masses[indices].astype(float)
    total_mass = float(np.sum(group_masses))
    assert total_mass > 0.0, "Total mass of atom group must be positive."
    return np.sum(coords * group_masses[:, None], axis=0) / total_mass


def _delta_min_image(
        delta: np.ndarray,
        box_vectors: np.ndarray | None,
        ) -> np.ndarray:
    if box_vectors is None:
        return delta
    # Fractional coordinates in the box basis, then wrap to [-0.5, 0.5).
    box = np.asarray(box_vectors, dtype=float)
    inv_box = np.linalg.inv(box)
    frac = delta @ inv_box
    frac -= np.floor(frac + 0.5)
    return frac @ box


def _dihedral(
        p0: np.ndarray,
        p1: np.ndarray,
        p2: np.ndarray,
        p3: np.ndarray,
        box_vectors: np.ndarray | None = None,
        ) -> float:
    b0 = _delta_min_image(p0 - p1, box_vectors)
    b1 = _delta_min_image(p2 - p1, box_vectors)
    b2 = _delta_min_image(p3 - p2, box_vectors)
    b1_norm = np.linalg.norm(b1)
    assert b1_norm > 0.0, "Dihedral middle bond has zero length."
    b1_unit = b1 / b1_norm
    v = b0 - np.dot(b0, b1_unit) * b1_unit
    w = b2 - np.dot(b2, b1_unit) * b1_unit
    x = np.dot(v, w)
    y = np.dot(np.cross(b1_unit, v), w)
    return float(np.arctan2(y, x))


def _validate_int_or_list_of_int(instance, attribute, value) -> None:
    if isinstance(value, int):
        return
    if isinstance(value, list) and len(value) > 0 \
            and all(isinstance(v, int) for v in value):
        return
    raise TypeError(
        f"{attribute.name} must be an int or a non-empty list of ints, "
        f"got {value!r}")


@define
class Pair_distance_feature:
    """
    Distance between two atoms or mass-weighted centers of mass (nm).
    """
    atom_i: int | list[int] = field(validator=_validate_int_or_list_of_int)
    atom_j: int | list[int] = field(validator=_validate_int_or_list_of_int)
    type: typing.Literal["pair_distance"] = "pair_distance"
    name: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))

    def column_name(self) -> str:
        if self.name is not None:
            return self.name
        return (
            f"pair_dist_{_format_index_group(self.atom_i)}_"
            f"{_format_index_group(self.atom_j)}"
        )

    def compute(
            self,
            positions: np.ndarray,
            box_vectors: np.ndarray | None = None,
            masses: np.ndarray | None = None,
            ) -> float:
        i_indices = _as_index_list(self.atom_i)
        j_indices = _as_index_list(self.atom_j)
        ri = _com(positions, i_indices, masses)
        rj = _com(positions, j_indices, masses)
        delta = _delta_min_image(rj - ri, box_vectors)
        return float(np.linalg.norm(delta))


@define
class Dihedral_feature:
    """
    Torsion angle defined by four atoms (radians).
    """
    atom_i: int = field(validator=validators.instance_of(int))
    atom_j: int = field(validator=validators.instance_of(int))
    atom_k: int = field(validator=validators.instance_of(int))
    atom_l: int = field(validator=validators.instance_of(int))
    type: typing.Literal["dihedral"] = "dihedral"
    name: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))

    def column_name(self) -> str:
        if self.name is not None:
            return self.name
        return (
            f"dihedral_{self.atom_i}_{self.atom_j}_"
            f"{self.atom_k}_{self.atom_l}"
        )

    def compute(
            self,
            positions: np.ndarray,
            box_vectors: np.ndarray | None = None,
            masses: np.ndarray | None = None,
            ) -> float:
        return _dihedral(
            positions[self.atom_i],
            positions[self.atom_j],
            positions[self.atom_k],
            positions[self.atom_l],
            box_vectors=box_vectors,
        )


features: typing.Dict[str, typing.Type[Base_feature]] = {
    "pair_distance": Pair_distance_feature,
    "dihedral": Dihedral_feature,
}


def structure_feature(obj: typing.Any) -> Base_feature:
    """
    Convert a dict or feature object into a concrete feature instance.
    """
    if isinstance(obj, (Pair_distance_feature, Dihedral_feature)):
        return obj
    if isinstance(obj, dict):
        typ = obj.get("type")
        assert typ in features, f"Unknown feature type: {typ}"
        return features[typ](**{k: v for k, v in obj.items() if k != "type"})
    raise TypeError(f"Cannot structure feature from {type(obj)}: {obj!r}")


def structure_feature_list(
        objs: typing.Sequence[typing.Any] | None,
        ) -> list[Base_feature]:
    """
    Convert a list of feature dicts/objects into feature instances.
    """
    if objs is None:
        return []
    return [structure_feature(obj) for obj in objs]
