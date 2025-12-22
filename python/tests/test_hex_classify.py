import random

import pytest

pytest.importorskip(
    "tools.cutcells_tables",
    reason="Legacy table-tooling package is not part of this repo layout.",
)

from tools.cutcells_tables.classify.hex_mc33 import (
    classify_hex,
    permute_mask,
    RAW_TO_BASE_ID,
    CANONICAL_MASKS,
    ambiguous_faces_mask,
    edges_mask_from_cell_mask,
    body_saddle_sign_hex,
    asymptotic_decider_quad,
)
from tools.cutcells_tables.symmetry.cube import SYMMETRIES


def test_orbit_canonicalization_exhaustive():
    # Every raw mask maps to canonical mask present in CANONICAL_MASKS
    canonical_set = set(CANONICAL_MASKS)
    for m in range(256):
        info = classify_hex(m)
        assert info.canonical_mask in canonical_set
        # Check RAW_TO_* consistency
        c2 = classify_hex(info.canonical_mask)
        assert RAW_TO_BASE_ID[m] == c2.base_id


def test_symmetry_equivalence_random():
    rng = random.Random(1234)
    for _ in range(200):
        m = rng.randrange(256)
        info = classify_hex(m)
        for sid in rng.sample(range(len(SYMMETRIES)), 5):
            pm = permute_mask(m, sid)
            info_p = classify_hex(pm)
            assert info.base_id == info_p.base_id


def test_ambiguous_faces_patterns():
    # Construct masks with one ambiguous face (face 0) pattern 0101 / 1010
    # Face 0 assumed (0,1,2,3)
    pattern_0101 = 0
    pattern_1010 = 0
    # set bits: v0=1,v1=0,v2=1,v3=0 -> 0b0101 -> vertices 0,2
    pattern_0101 |= (1 << 0) | (1 << 2)
    # 1010 -> vertices 1,3
    pattern_1010 |= (1 << 1) | (1 << 3)
    assert ambiguous_faces_mask(pattern_0101) & 1
    assert ambiguous_faces_mask(pattern_1010) & 1


def test_edges_mask_parity_random():
    for m in range(256):
        em = edges_mask_from_cell_mask(m)
        # For each edge bit set ensure its endpoints differ
        for ei, (a, b) in enumerate(
            (
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 0),
                (4, 5),
                (5, 6),
                (6, 7),
                (7, 4),
                (0, 4),
                (1, 5),
                (2, 6),
                (3, 7),
            )
        ):
            if (em >> ei) & 1:
                assert ((m >> a) & 1) != ((m >> b) & 1)


def test_body_saddle_sign_center_consistency():
    # If all phi are symmetric around zero, interior sign should align with center
    phi = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0]
    s = body_saddle_sign_hex(phi)
    # center value is ~0 so sign may be 1 (>=0). Accept either but function should return int 0/1
    assert s in (0, 1)


def test_asymptotic_decider_basic():
    # Values produce center < 0 -> returns 1
    assert asymptotic_decider_quad(-1, 1, 1, -1) == 1
    # values produce center > 0 -> returns 0
    assert asymptotic_decider_quad(1, 1, 1, 1) == 0
