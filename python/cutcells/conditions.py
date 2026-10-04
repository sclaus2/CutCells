# Copyright (c) 2026 ONERA
# Authors: Susanne Claus
# This file is part of CutCells
# SPDX-License-Identifier: MIT
"""Conditions: selections built from the level sets themselves.

A level set compared with 0 is a condition, and ``&`` (and) and ``|`` (or)
combine conditions, as boolean masks do in pandas. Strings stay as the
shorthand, as ``df.query`` does in pandas::

    result = cutcells.cut(mesh, [body, fluid], names=["body", "fluid"])
    solid = result[body < 0]
    around = result[(body > 0) & (fluid < 0)]
    surface = result[body == 0]
    result["body > 0 and fluid < 0"]   # the same as `around`
    result["body == 0"]                # strings take == as well as =

- ``<``, ``>`` and ``==`` against 0 give conditions; against any other number
  they raise.
- ``&`` binds more tightly than ``<``, so each comparison needs parentheses.
  There is no ``~``: ``body > 0`` is the opposite of ``body < 0``.
- A condition has no truth value, so ``(body < 0) and (fluid > 0)`` raises.
- ``==`` gives a condition only against a number: level sets still compare by
  identity (``body == fluid`` is False, ``body in [fluid, body]`` works) and
  hash by identity, so they work in sets and as dict keys. ``!=`` against a
  number raises.

The result writes a condition as its selection string, naming each level set
as the cut named it, matched by identity, so the C++ side keeps taking strings.
"""

import numbers

from ._cutcellscpp import AnalyticLevelSet, LevelSetFunction_float32, LevelSetFunction_float64


class Condition:
    """A selection built from level sets: or-terms of and-clauses, each clause
    a level set and its relation to 0, "<", ">" or "=". ``&`` distributes over
    the terms, so the text needs no parentheses."""

    __slots__ = ("terms",)

    def __init__(self, terms):
        self.terms = tuple(tuple(term) for term in terms)

    def __and__(self, other):
        if not isinstance(other, Condition):
            return NotImplemented
        return Condition(a + b for a in self.terms for b in other.terms)

    def __or__(self, other):
        if not isinstance(other, Condition):
            return NotImplemented
        return Condition(self.terms + other.terms)

    def __bool__(self):
        raise ValueError(
            "a condition has no truth value: combine conditions with & and |, not with and/or, "
            "with each comparison in parentheses, e.g. (body < 0) & (fluid > 0)")

    def __eq__(self, other):
        if not isinstance(other, Condition):
            return NotImplemented
        return self.terms == other.terms

    def __hash__(self):
        return hash(self.terms)

    def text(self, names=None):
        """The selection string, e.g. "a < 0 and b = 0 or c > 0". names[level_set]
        gives each level set's name; by default its own name, or a placeholder
        for a level set without one (an AnalyticLevelSet)."""
        name = _own_name if names is None else names.__getitem__
        return " or ".join(" and ".join(f"{name(ls)} {relation} 0" for ls, relation in term) for term in self.terms)

    def __str__(self):
        return self.text()

    def __repr__(self):
        return f"Condition({self.text()!r})"


def _own_name(level_set):
    """A level set's own name, or a placeholder for one without (an AnalyticLevelSet)."""
    name = getattr(level_set, "name", None)
    return name if name else f"<{type(level_set).__name__} at {id(level_set):#x}>"


def _clause(level_set, relation, value):
    """The condition level_set <relation> 0; NotImplemented unless value is a number."""
    if not isinstance(value, numbers.Number):
        return NotImplemented
    if value != 0:
        raise ValueError(f"a condition compares a level set with 0, not {value!r}")
    return Condition([[(level_set, relation)]])


class _Comparisons:
    """What the level-set classes get: compared with 0 they give conditions;
    among themselves they compare and hash by identity, as before."""

    def __lt__(self, other):
        return _clause(self, "<", other)

    def __gt__(self, other):
        return _clause(self, ">", other)

    def __eq__(self, other):
        return _clause(self, "=", other)

    def __ne__(self, other):
        if isinstance(other, numbers.Number):
            raise ValueError("level_set != 0 is not a selection: (level_set < 0) | (level_set > 0) "
                             "is everything off its zero set")
        return NotImplemented

    __hash__ = object.__hash__


for _cls in (LevelSetFunction_float32, LevelSetFunction_float64, AnalyticLevelSet):
    for _name in ("__lt__", "__gt__", "__eq__", "__ne__", "__hash__"):
        setattr(_cls, _name, getattr(_Comparisons, _name))


def selection_text(selection, level_sets, names):
    """A selection as the string that a cut of level_sets, named names, parses:
    a string as it is, a condition with each of its level sets named as in the
    cut, found by identity among level_sets."""
    if isinstance(selection, str):
        return selection
    if not isinstance(selection, Condition):
        raise TypeError(f"a selection is a condition such as phi < 0, or a string such as \"phi < 0\", "
                        f"not {type(selection).__name__}")
    found = {}
    for level_set, name in zip(level_sets, names):
        found.setdefault(level_set, []).append(name)
    used = {}
    for term in selection.terms:
        for level_set, _ in term:
            given = found.get(level_set, [])
            if not given:
                raise ValueError(f"the condition's level set {_own_name(level_set)} is not one of the cut's "
                                 f"level sets: {', '.join(names)}")
            if len(given) > 1:
                raise ValueError(f"the condition's level set {_own_name(level_set)} was cut twice, as "
                                 f"{' and '.join(given)}: select it with a string")
            used[level_set] = given[0]
    return selection.text(used)
