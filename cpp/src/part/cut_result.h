// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <concepts>
#include <cstdint>
#include <map>
#include <span>
#include <string>
#include <vector>

#include "../cell_flags.h"
#include "../cell_topology.h"
#include "../cell_types.h"
#include "../level_set.h"
#include "../mesh_view.h"
#include "entities.h"

/// The front end: every cell classified by every level set, mesh parts
/// selected by expressions, and their quadrature and visualisation by backend.
namespace cutcells::part
{

/// The facets of a cell: its faces in 3D, its edges in 2D, its vertices in 1D.
inline int num_facets(cell::type type)
{
    return cell::num_sub_entities(type, cell::get_tdim(type) - 1);
}

/// The vertices of facet @p f (Basix numbering of the cell and of the facet's
/// vertices).
inline std::span<const int> facet_vertices(cell::type type, int f)
{
    return cell::sub_entity_vertices(type, cell::get_tdim(type) - 1, f);
}

inline cell::type facet_type(cell::type type, int f)
{
    return cell::sub_entity_type(type, cell::get_tdim(type) - 1, f);
}

/// Options of cut().
struct ClassifyOptions
{
    /// Bisections of a cell's box before its sign is given up as unknown; the
    /// cell then counts as cut. Bernstein bounds rarely need more than a few,
    /// Taylor models of distance functions more.
    int max_depth = 12;
};

/// Every cell classified by every level set as inside (phi < 0), outside
/// (phi > 0) or cut, and the facets, edges and vertices on which a level set
/// vanishes.
///
/// A level set made of level sets is replaced by its own: level_sets and
/// level_set_names hold them, named after it ("body.hole"), and insides keeps
/// its inside for selections.
template <std::floating_point T, std::integral I = int>
struct CutResult
{
    const MeshView<T, I>* mesh = nullptr;                  ///< not owned
    std::vector<const LevelSetFunction<T, I>*> level_sets; ///< not owned
    std::vector<std::string> level_set_names;
    /// The level sets made of level sets, by name ("body", or "body.part"
    /// when nested), with their insides as given: selections of their own
    /// level sets, by the names those have within them.
    std::map<std::string, std::string> insides;
    int num_cells = 0;
    /// domains[l * num_cells + c]: inside, outside or intersected
    std::vector<cell::domain> domains;
    /// Cells that some level set cuts, ascending.
    std::vector<I> cut_cells;
    /// The vertices, edges and facets on which a level set vanishes,
    /// entities[d] those of dimension d < tdim, with all the cells that share
    /// each, their local numbers and their sides of each level set vanishing
    /// there.
    std::array<EntityTable<I>, 3> entities;
    /// Facets (faces in 3D, edges in 2D, vertices of intervals) in the zero
    /// set of a level set, from entities[tdim - 1]: each facet once per level
    /// set in its zero mask, with an owning cell and the facet in that cell
    /// (Basix numbering). The owner is the cell on the negative side, or the
    /// lowest sharing cell if none or several are. For codes that take a
    /// facet's rule from one cell; rules on facets come from the entities
    /// (entity_rules, part::EntityPolicy).
    std::vector<int> zero_face_level_sets;
    std::vector<I> zero_face_cells;
    std::vector<std::int8_t> zero_face_local;

    int n_level_sets() const { return static_cast<int>(level_sets.size()); }
    int tdim() const { return num_cells > 0 ? cell::get_tdim(mesh->cell_type(I(0))) : 0; }
    int n_zero_faces() const { return static_cast<int>(zero_face_cells.size()); }
    cell::domain domain(int level_set, I cell_id) const
    {
        return domains[static_cast<std::size_t>(level_set) * static_cast<std::size_t>(num_cells)
                       + static_cast<std::size_t>(cell_id)];
    }
};

/// @brief Bit l set if level set l cuts @p cell_id.
template <std::floating_point T, std::integral I>
std::uint64_t cut_mask(const CutResult<T, I>& result, I cell_id)
{
    std::uint64_t mask = 0;
    for (int l = 0; l < result.n_level_sets(); ++l)
        if (result.domain(l, cell_id) == cell::domain::intersected)
            mask |= std::uint64_t(1) << l;
    return mask;
}

/// @brief Classify every cell of @p mesh by every level set, and find the
/// facets, edges and vertices on which they vanish.
///
/// The cells quadrays takes (intervals, triangles, quadrilaterals, tetrahedra,
/// hexahedra, prisms and pyramids, in a space of their dimension or a higher
/// one, as facets are) are classified by bounds of the level set itself:
/// Bernstein coefficients of a Pk level set, or Taylor models of an analytic
/// one, over the cell and its halves (cell_sign in classify.h). A cell whose
/// sign the bounds cannot prove counts as cut, which the backends handle.
/// Other cells use the signs of their Bernstein coefficients.
///
/// A level set vanishes on a facet, edge or vertex if the bound of its
/// restriction there (the restriction of its Bernstein form, Taylor bounds
/// over the entity itself, a value at a vertex) is within one threshold:
/// touch_tolerance (64 eps) times the level set's largest magnitude over the
/// cells at its zero set (those it cuts or is about 0 at a vertex of). The
/// same threshold decides touching zero there: a cell sharing such an entity
/// and found cut is classified again with values within it counted as
/// touching, so that a level set does not both vanish on an entity and cut a
/// cell next to it. Each sharing cell's side of a vanishing level set comes
/// from its domain, or for a cut cell from the level set just inside it, next
/// to the entity. A vanishing level set counts on the entity (EntityTable::zero)
/// where it changes sign across it, or on the mesh boundary.
///
/// A level set made of level sets (LevelSetFunction::level_sets and inside)
/// is replaced by its own level sets, depth first, named after it: "body.a",
/// "body.b", or "body.c.d" if c is made of level sets in turn.
///
/// @param level_sets  Pk level sets (dof values) or analytic ones, or level
///                    sets made of them; they and @p mesh must outlive the
///                    result
/// @throws std::invalid_argument for level sets with neither, more than 64
///         level sets (counting those level sets made of level sets are made
///         of), analytic level sets on cells quadrays does not take, two
///         level sets of one name, a name that reads as one of the level
///         sets of a level set made of level sets ("body.x" next to body),
///         and level sets made of level sets that check_level_sets rejects
template <std::floating_point T, std::integral I>
CutResult<T, I> cut(const MeshView<T, I>& mesh, std::span<const LevelSetFunction<T, I>> level_sets,
                    const ClassifyOptions& options = {});

} // namespace cutcells::part
