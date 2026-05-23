"""Unit tests for src/particle/external_field.py.

Covers:
- ExternalField.uniform  factory
- ExternalField.non_uniform factory
- ExternalField.__init__ direct construction
- ExternalField.get_energy (both modes)
- ExternalField.set_energy
- ExternalField.nodes
- ExternalField.__repr__
- Input validation (TypeError / ValueError guards)
"""

import math

import pytest

from src.particle.external_field import ExternalField, FieldMode

# ---------------------------------------------------------------------------
# Constants shared across tests
# ---------------------------------------------------------------------------

_ORIGIN_2D: tuple[int, ...] = (0, 0)
_POINT_3D: tuple[int, ...] = (1, 2, 3)
_FAR_AWAY: tuple[int, ...] = (99, 99, 99)


# ===========================================================================
# ExternalField.uniform
# ===========================================================================


class TestUniformFactory:
    def test_mode_is_uniform(self):
        field = ExternalField.uniform(strength=-1.0)
        assert field.mode is FieldMode.UNIFORM

    def test_default_energy_equals_strength(self):
        field = ExternalField.uniform(strength=2.5)
        assert field.default_energy == 2.5

    def test_get_energy_returns_strength_for_any_coord(self):
        strength = -0.75
        field = ExternalField.uniform(strength=strength)
        assert field.get_energy(_ORIGIN_2D) == strength
        assert field.get_energy(_POINT_3D) == strength
        assert field.get_energy(_FAR_AWAY) == strength

    def test_zero_strength(self):
        field = ExternalField.uniform(strength=0.0)
        assert field.get_energy((5, -3)) == 0.0

    def test_large_positive_strength(self):
        field = ExternalField.uniform(strength=1e12)
        assert field.get_energy((0,)) == pytest.approx(1e12)

    def test_negative_strength(self):
        field = ExternalField.uniform(strength=-100.0)
        assert field.get_energy((0, 0)) == -100.0

    def test_nodes_returns_empty_dict_by_default(self):
        field = ExternalField.uniform(strength=1.0)
        assert field.nodes() == {}

    def test_nan_strength_raises(self):
        with pytest.raises(ValueError, match="finite"):
            ExternalField.uniform(strength=math.nan)

    def test_inf_strength_raises(self):
        with pytest.raises(ValueError, match="finite"):
            ExternalField.uniform(strength=math.inf)


# ===========================================================================
# ExternalField.non_uniform
# ===========================================================================


class TestNonUniformFactory:
    def test_mode_is_non_uniform(self):
        field = ExternalField.non_uniform({_ORIGIN_2D: -2.0})
        assert field.mode is FieldMode.NON_UNIFORM

    def test_explicit_node_energy(self):
        field = ExternalField.non_uniform({_ORIGIN_2D: -3.0, _POINT_3D: 1.5})
        assert field.get_energy(_ORIGIN_2D) == -3.0
        assert field.get_energy(_POINT_3D) == 1.5

    def test_missing_node_returns_default(self):
        field = ExternalField.non_uniform({_ORIGIN_2D: -2.0}, default_energy=0.0)
        assert field.get_energy(_FAR_AWAY) == 0.0

    def test_custom_default_energy(self):
        field = ExternalField.non_uniform({}, default_energy=-0.5)
        assert field.default_energy == -0.5
        assert field.get_energy((1, 1)) == -0.5

    def test_empty_map_uses_default_everywhere(self):
        field = ExternalField.non_uniform({}, default_energy=7.0)
        assert field.get_energy((0, 0)) == 7.0

    def test_overwrite_existing_node(self):
        field = ExternalField.non_uniform({_ORIGIN_2D: 1.0})
        field.set_energy(_ORIGIN_2D, 99.0)
        assert field.get_energy(_ORIGIN_2D) == 99.0

    def test_nodes_returns_copy(self):
        energy_map = {_ORIGIN_2D: -2.0, _POINT_3D: 1.0}
        field = ExternalField.non_uniform(energy_map)
        returned = field.nodes()
        assert returned == energy_map
        # Mutating the returned dict must not affect the internal state.
        returned[_FAR_AWAY] = 999.0
        assert _FAR_AWAY not in field.nodes()

    def test_non_dict_map_raises_type_error(self):
        with pytest.raises(TypeError, match="dict"):
            ExternalField.non_uniform([(0, 0), -1.0])  # type: ignore[arg-type]

    def test_nan_in_map_raises_value_error(self):
        with pytest.raises(ValueError, match="non-finite"):
            ExternalField.non_uniform({_ORIGIN_2D: math.nan})

    def test_inf_in_map_raises_value_error(self):
        with pytest.raises(ValueError, match="non-finite"):
            ExternalField.non_uniform({_ORIGIN_2D: math.inf})

    def test_negative_inf_in_map_raises_value_error(self):
        with pytest.raises(ValueError, match="non-finite"):
            ExternalField.non_uniform({_ORIGIN_2D: -math.inf})

    def test_nan_default_energy_raises_value_error(self):
        with pytest.raises(ValueError, match="finite"):
            ExternalField.non_uniform({}, default_energy=math.nan)


# ===========================================================================
# ExternalField.get_energy – input validation
# ===========================================================================


class TestGetEnergyValidation:
    def test_non_tuple_raises_type_error_uniform(self):
        field = ExternalField.uniform(strength=1.0)
        with pytest.raises(TypeError, match="tuple"):
            field.get_energy([0, 0])  # type: ignore[arg-type]

    def test_non_tuple_raises_type_error_non_uniform(self):
        field = ExternalField.non_uniform({_ORIGIN_2D: -1.0})
        with pytest.raises(TypeError, match="tuple"):
            field.get_energy("0,0")  # type: ignore[arg-type]

    def test_1d_coord_accepted(self):
        field = ExternalField.uniform(strength=-1.0)
        assert field.get_energy((5,)) == -1.0

    def test_high_dimensional_coord_accepted(self):
        field = ExternalField.uniform(strength=3.0)
        coord = tuple(range(10))  # 10-dimensional
        assert field.get_energy(coord) == 3.0

    def test_negative_coord_components_accepted(self):
        field = ExternalField.non_uniform({(-1, -1): 5.0})
        assert field.get_energy((-1, -1)) == 5.0


# ===========================================================================
# ExternalField.set_energy
# ===========================================================================


class TestSetEnergy:
    def test_set_energy_adds_new_node(self):
        field = ExternalField.non_uniform({})
        field.set_energy((3, 4), -7.0)
        assert field.get_energy((3, 4)) == -7.0

    def test_set_energy_overwrites_existing_node(self):
        field = ExternalField.non_uniform({_ORIGIN_2D: 1.0})
        field.set_energy(_ORIGIN_2D, -5.0)
        assert field.get_energy(_ORIGIN_2D) == -5.0

    def test_set_energy_on_uniform_overrides_specific_node(self):
        field = ExternalField.uniform(strength=1.0)
        field.set_energy(_ORIGIN_2D, -99.0)
        # The overridden node now returns the new value even in UNIFORM mode.
        assert field.get_energy(_ORIGIN_2D) == -99.0
        # Other nodes are unaffected.
        assert field.get_energy(_POINT_3D) == 1.0

    def test_set_energy_nan_raises(self):
        field = ExternalField.non_uniform({})
        with pytest.raises(ValueError, match="finite"):
            field.set_energy(_ORIGIN_2D, math.nan)

    def test_set_energy_inf_raises(self):
        field = ExternalField.non_uniform({})
        with pytest.raises(ValueError, match="finite"):
            field.set_energy(_ORIGIN_2D, math.inf)

    def test_set_energy_non_tuple_coord_raises(self):
        field = ExternalField.non_uniform({})
        with pytest.raises(TypeError, match="tuple"):
            field.set_energy([0, 0], 1.0)  # type: ignore[arg-type]


# ===========================================================================
# ExternalField.nodes
# ===========================================================================


class TestNodes:
    def test_nodes_returns_dict(self):
        field = ExternalField.non_uniform({_ORIGIN_2D: 1.0})
        assert isinstance(field.nodes(), dict)

    def test_nodes_after_set_energy(self):
        field = ExternalField.non_uniform({})
        field.set_energy((1, 2), 5.0)
        assert field.nodes() == {(1, 2): 5.0}

    def test_nodes_copy_independence(self):
        field = ExternalField.non_uniform({_ORIGIN_2D: 1.0})
        copy = field.nodes()
        copy[_FAR_AWAY] = 999.0
        assert _FAR_AWAY not in field.nodes()


# ===========================================================================
# ExternalField.__repr__
# ===========================================================================


class TestRepr:
    def test_repr_contains_mode_name(self):
        field = ExternalField.uniform(strength=1.0)
        assert "UNIFORM" in repr(field)

    def test_repr_contains_default_energy(self):
        field = ExternalField.uniform(strength=-3.5)
        assert "-3.5" in repr(field)

    def test_repr_non_uniform_shows_node_count(self):
        field = ExternalField.non_uniform({_ORIGIN_2D: 1.0, _POINT_3D: 2.0})
        assert "2" in repr(field)


# ===========================================================================
# ExternalField.__init__ direct construction
# ===========================================================================


class TestDirectConstruction:
    def test_direct_init_with_valid_args(self):
        field = ExternalField(
            mode=FieldMode.NON_UNIFORM,
            energy_map={_ORIGIN_2D: -1.0},
            default_energy=0.0,
        )
        assert field.get_energy(_ORIGIN_2D) == -1.0
        assert field.get_energy(_FAR_AWAY) == 0.0

    def test_direct_init_non_dict_raises(self):
        with pytest.raises(TypeError, match="dict"):
            ExternalField(
                mode=FieldMode.UNIFORM,
                energy_map="bad",  # type: ignore[arg-type]
                default_energy=0.0,
            )

    def test_direct_init_non_finite_default_raises(self):
        with pytest.raises(ValueError, match="finite"):
            ExternalField(
                mode=FieldMode.UNIFORM,
                energy_map={},
                default_energy=math.inf,
            )
