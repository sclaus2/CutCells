// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "cut_mesh.h"




#include "utils.h"
#include "cell_topology.h"
#include "reference_cell.h"
#include <map>
#include <algorithm>
#include <cstdint>

namespace
{
  struct VertexKey
  {
    // 0: original vertex (global vertex id)
    // 1: edge intersection (global vertex id pair)
    // 2: special point (parent cell id + sid)
    uint8_t kind = 0;
    int32_t a = 0;
    int32_t b = 0;

    bool operator==(const VertexKey& other) const noexcept
    {
      return kind == other.kind && a == other.a && b == other.b;
    }
  };

  struct VertexKeyHash
  {
    std::size_t operator()(const VertexKey& k) const noexcept
    {
      // Basic mix; good enough for small integer keys.
      const std::size_t h0 = std::hash<int32_t>{}(k.a);
      const std::size_t h1 = std::hash<int32_t>{}(k.b);
      const std::size_t hk = std::hash<uint8_t>{}(k.kind);
      return h0 ^ (h1 + 0x9e3779b97f4a7c15ULL + (h0 << 6) + (h0 >> 2)) ^ (hk << 1);
    }
  };

  template <std::floating_point T>
  inline void append_vertex_coords(std::vector<T>& out, const std::vector<T>& in,
                   int local_vertex_id, int gdim)
  {
    const int base = local_vertex_id * gdim;
    out.insert(out.end(), in.data() + base, in.data() + base + gdim);
  }
}

namespace cutcells::mesh
{
    template <std::floating_point T>
    void str(const CutCells<T> &cut_mesh)
    {
        std::cout << "CutCells: " << std::endl;
        int cnt = 0;

        for(auto &cell: cut_mesh._cut_cells)
        {
            std::cout << "Cut Cell ";
            std::cout << cnt << ": ";
            cnt++;
            std::cout << "vertex coordinates=[";
            for(auto &i: cell._vertex_coords)
            {
                std::cout << i << ", ";
            }
            std::cout << "]" << std::endl;
            std::cout << "connectivity=[";
            const int ncells = cutcells::cell::num_cells(cell);
            for(int i=0;i<ncells;i++)
            {
              const auto verts = cutcells::cell::cell_vertices(cell, i);
                std::cout << i << ": ";
              for(int j=0;j<verts.size();j++)
                {
                std::cout << verts[j] << ", ";
                }
            }
            std::cout << "]" << std::endl;
        }

        std::cout << "Cell Types= ";
        for(auto &type : cut_mesh._types)
        {
            std::cout << cell_type_to_str(type) << ", ";
        }
        std::cout << std::endl;

        std::cout << "Parent map=[";
        for(int i=0;i<cut_mesh._parent_map.size();i++)
        {
            std::cout << i << ": " << cut_mesh._parent_map[i] << std::endl;
        }
        std::cout << "]" << std::endl;
    }

    std::unordered_map<int, std::vector<int>> create_parent_cut_cells_map(const std::span<int> parent_map)
    {
        std::unordered_map<int, std::vector<int>> parent_cut_cell_map;

        for(std::size_t i=0; i< parent_map.size();i++)
        {
            parent_cut_cell_map[parent_map[i]].push_back(i);
        }

        return parent_cut_cell_map;
    }

    template <std::floating_point T>
    int get_num_cells(const cutcells::mesh::CutCells<T>& cut_mesh)
    {
      int num_cells = 0;

      for(auto & cut_cell : cut_mesh._cut_cells)
      {
        num_cells += cutcells::cell::num_cells(cut_cell);
      }
      return num_cells;
    }

    //create cutmesh from cutcells by merging cut cells
    template <std::floating_point T>
    cutcells::mesh::CutMesh<T> create_cut_mesh(CutCells<T>& cut_cells)
    {
      CutMesh<T> cut_mesh;

      if (cut_cells._cut_cells.empty())
      {
        cut_mesh._gdim = 0;
        cut_mesh._tdim = 0;
        cut_mesh._num_cells = 0;
        cut_mesh._num_vertices = 0;
        return cut_mesh;
      }

      std::size_t gdim = cut_cells._cut_cells[0]._gdim;
      std::size_t tdim = cut_cells._cut_cells[0]._tdim;

      cut_mesh._gdim = gdim;
      cut_mesh._tdim = tdim;

      // Count cells and connectivity in a single pass.
      int num_cells = 0;
      int num_connectivity = 0;
      for (const auto& cut_cell : cut_cells._cut_cells)
      {
        const int nc = cutcells::cell::num_cells(cut_cell);
        num_cells += nc;
        if (nc > 0)
          num_connectivity += cut_cell._offset.back();
      }

      cut_mesh._offset.resize(num_cells + 1);
      cut_mesh._connectivity.resize(num_connectivity);
      // Reserve vertex coord storage up-front to avoid repeated reallocation
      // during the merge loop.  num_connectivity is an over-estimate of the
      // number of unique vertices (dedup reduces it), but it is O(correct).
      cut_mesh._vertex_coords.reserve(num_connectivity * gdim);

      // either two or one; allow missing parent map
      int num_parents = 0;
      if (!cut_cells._parent_map.empty())
        num_parents = static_cast<int>(cut_cells._parent_map.size() / cut_cells._cut_cells.size());
      if (num_parents <= 0)
        num_parents = 1;

      cut_mesh._num_cells = num_cells;
      cut_mesh._types.resize(num_cells);
      cut_mesh._parent_map.resize(num_cells * num_parents, -1);

      int sub_cell_offset = 0;
      int element_offset = 0;
      int cnt = 0;

      // Adaptive dedup map: flat linear-scan vector for small meshes (typ. a
      // handful of cut cells — fast due to cache locality and no hashing),
      // falling back to an unordered_map for larger meshes.
      // Both are thread_local so their allocations are amortised across calls.
      constexpr int kFlatMapThreshold = 128;
      const bool use_flat = (num_connectivity <= kFlatMapThreshold);

      thread_local std::vector<std::pair<VertexKey, int>> flat_vertex_ids;
      thread_local std::unordered_map<VertexKey, int, VertexKeyHash> hash_vertex_ids;

      flat_vertex_ids.clear();
      if (!use_flat)
      {
        hash_vertex_ids.clear();
        hash_vertex_ids.reserve(static_cast<std::size_t>(num_connectivity));
      }

      // Helper lambdas that abstract over the two map implementations.
      // find_vertex: returns the merged id, or -1 if not yet seen.
      // insert_vertex: adds a new entry; only called when find returns -1.
      auto find_vertex = [&](const VertexKey& key) -> int
      {
        if (use_flat)
        {
          for (const auto& [k, v] : flat_vertex_ids)
            if (k == key) return v;
          return -1;
        }
        const auto it = hash_vertex_ids.find(key);
        return it == hash_vertex_ids.end() ? -1 : it->second;
      };
      auto insert_vertex = [&](const VertexKey& key, int id)
      {
        if (use_flat)
          flat_vertex_ids.emplace_back(key, id);
        else
          hash_vertex_ids.emplace(key, id);
      };

      //all cutcells in vector above should have the same gdim and tdim
      for(auto & cut_cell : cut_cells._cut_cells)
      {
        if(cut_cell._vertex_coords.size()==0)
        {
          continue;
        }
        //check that current cut_cell has same dimensions and parent
        if((cut_cell._gdim!=cut_mesh._gdim)||(cut_cell._tdim!=cut_mesh._tdim))
        {
          std::cout << "gdim: (" << gdim << ", " << cut_cell._gdim << ")" << std::endl;
          std::cout << "tdim: (" << tdim << ", " << cut_cell._tdim << ")" << std::endl;
          throw std::runtime_error ("Error in merging cutcells into mesh as they have differing dimensions");
        }

        int num_cut_cell_vertices = cut_cell._vertex_coords.size()/gdim;

        int local_num_cells = cutcells::cell::num_cells(cut_cell);

        // Map from vertex id in current cutcell to merged cutmesh
        // thread_local: reuse the underlying storage each call to avoid per-cell allocation.
        thread_local std::vector<int> local_merged_vertex_ids;
        local_merged_vertex_ids.assign(num_cut_cell_vertices, -1);

        // Determine whether we can use fast token-based dedup for this cut cell.
        const bool has_tokens = (static_cast<int>(cut_cell._vertex_parent_entity.size()) == num_cut_cell_vertices);
        const int parent_vertices = cell::get_num_vertices(cut_cell._parent_cell_type);
        const bool has_parent_ids = (static_cast<int>(cut_cell._parent_vertex_ids.size()) == parent_vertices);
        const int parent_cell_id = cut_cells._parent_map.empty()
                                       ? -1
                                       : cut_cells._parent_map[cnt * num_parents + 0];
        const bool can_fast_dedup = has_tokens && has_parent_ids && (parent_cell_id >= 0);

        if (can_fast_dedup)
        {
          const auto parent_edges = cell::edges(cut_cell._parent_cell_type);

          for (int local_id = 0; local_id < num_cut_cell_vertices; ++local_id)
          {
            const int32_t token = cut_cell._vertex_parent_entity[local_id];
            VertexKey key;

            if (token >= 100 && token < 200)
            {
              const int ref_vid = static_cast<int>(token - 100);
              if (ref_vid < 0 || ref_vid >= parent_vertices)
              {
                // Fallback for unexpected token.
                key.kind = 2;
                key.a = parent_cell_id;
                key.b = token;
              }
              else
              {
                key.kind = 0;
                key.a = static_cast<int32_t>(cut_cell._parent_vertex_ids[ref_vid]);
                key.b = 0;
              }
            }
            else if (token >= 200)
            {
              key.kind = 2;
              key.a = parent_cell_id;
              key.b = static_cast<int32_t>(token - 200);
            }
            else
            {
              const int edge_id = static_cast<int>(token);
              if (edge_id < 0 || edge_id >= static_cast<int>(parent_edges.size()))
              {
                key.kind = 2;
                key.a = parent_cell_id;
                key.b = token;
              }
              else
              {
                // token is a VTK edge index; cell_topology.h uses Basix ordering.
                const int basix_eid = cell::vtk_to_basix_edge(cut_cell._parent_cell_type, edge_id);
                const int bv0 = parent_edges[basix_eid][0];
                const int bv1 = parent_edges[basix_eid][1];
                // _parent_vertex_ids is indexed by VTK vertex; convert.
                const int vtk_v0 = cell::basix_to_vtk_vertex(cut_cell._parent_cell_type, bv0);
                const int vtk_v1 = cell::basix_to_vtk_vertex(cut_cell._parent_cell_type, bv1);
                const int32_t gv0 = static_cast<int32_t>(cut_cell._parent_vertex_ids[vtk_v0]);
                const int32_t gv1 = static_cast<int32_t>(cut_cell._parent_vertex_ids[vtk_v1]);
                key.kind = 1;
                key.a = std::min(gv0, gv1);
                key.b = std::max(gv0, gv1);
              }
            }

            const int existing = find_vertex(key);
            if (existing == -1)
            {
              const int merged_vertex_id = static_cast<int>(cut_mesh._vertex_coords.size() / gdim);
              append_vertex_coords<T>(cut_mesh._vertex_coords, cut_cell._vertex_coords, local_id, static_cast<int>(gdim));
              insert_vertex(key, merged_vertex_id);
              local_merged_vertex_ids[local_id] = merged_vertex_id;
            }
            else
            {
              local_merged_vertex_ids[local_id] = existing;
            }
          }
        }
        else
        {
          // Slow fallback: coordinate-based dedup (kept for compatibility).
          for (int local_id = 0; local_id < num_cut_cell_vertices; ++local_id)
          {
            const int id = cutcells::utils::vertex_exists<T>(cut_mesh._vertex_coords, cut_cell._vertex_coords,
                                                            local_id, static_cast<int>(gdim));
            if (id == -1)
            {
              const int merged_vertex_id = static_cast<int>(cut_mesh._vertex_coords.size() / gdim);
              append_vertex_coords<T>(cut_mesh._vertex_coords, cut_cell._vertex_coords, local_id, static_cast<int>(gdim));
              local_merged_vertex_ids[local_id] = merged_vertex_id;
            }
            else
            {
              local_merged_vertex_ids[local_id] = id;
            }
          }
        }

        for(int i=0;i<local_num_cells;i++)
        {
          auto vertices = cutcells::cell::cell_vertices(cut_cell, i);
          int num_vertices = vertices.size();
          for(int j=0;j<num_vertices;j++)
          {
            int64_t index = vertices[j];
              cut_mesh._connectivity[element_offset+j] = local_merged_vertex_ids[static_cast<int>(index)];
          }

          cut_mesh._offset[sub_cell_offset+i] = element_offset;
          element_offset += num_vertices;
          // one type per cell
          cut_mesh._types[sub_cell_offset+i]=cut_cell._types[i];

          for(int j=0;j<num_parents;j++)
          {
            const int out_parent_index = (sub_cell_offset + i) * num_parents + j;
            const int in_parent_index = cnt * num_parents + j;
            cut_mesh._parent_map[out_parent_index] = cut_cells._parent_map[in_parent_index];
          }
        }

        sub_cell_offset+=local_num_cells;
        cnt++;
      }

      cut_mesh._offset[num_cells]=element_offset; //last offset

      cut_mesh._num_vertices = cut_mesh._vertex_coords.size()/cut_mesh._gdim;

      return cut_mesh;
    }

    //extract the cell ids that are in a given domain
    // either phi<0, phi = 0, or phi>0
    template <std::floating_point T>
    std::vector<int> locate_cells(std::span<const T> ls_vals, std::span<const T> points,
                                    std::span<const int> connectivity, std::span<const int> offset,
                                    std::span<const int> vtk_type,
                                    cell::cut_type ctype)
    {
      std::vector<int> cells;
      std::vector<T> level_set_values;

      int num_cells = vtk_type.size();
      for(std::size_t i=0;i<num_cells;i++)
      {
        int cell_offset = offset[i];
        cell::vtk_types vtype = static_cast<cell::vtk_types>(vtk_type[i]);
        cell::type cell_type = cell::map_vtk_type_to_cell_type(vtype);

        int num_vertices = cell::get_num_vertices(cell_type);

        if (static_cast<int>(level_set_values.size()) != num_vertices)
          level_set_values.resize(num_vertices);

        for(std::size_t j=0;j<num_vertices;j++)
        {
          int vertex_id = connectivity[cell_offset+j];
          level_set_values[j] = ls_vals[vertex_id];
        }

        cell::domain domain_type = cell::classify_cell_domain<T>(std::span(level_set_values.data(),level_set_values.size()));

        if((domain_type == cell::domain::inside) && (ctype == cell::cut_type::philt0))
        {
          cells.push_back(i);
        }
        if((domain_type == cell::domain::outside) && (ctype == cell::cut_type::phigt0))
        {
          cells.push_back(i);
        }
        if((domain_type == cell::domain::intersected) && (ctype == cell::cut_type::phieq0))
        {
          cells.push_back(i);
        }
      }

      return cells;
    }

    // cut a vtk mesh
    // input: level set values at points
    // connectivity of cells with offsets and type
    template <std::floating_point T>
    cutcells::mesh::CutMesh<T> cut_vtk_mesh(std::span<const T> ls_vals, std::span<const T> points,
                                            std::span<const int> connectivity, std::span<const int> offset,
                                            std::span<const int> vtk_type,
                                            const std::string& cut_type_str)
    {
      return cut_vtk_mesh<T>(ls_vals, points,
                             connectivity, offset,
                             vtk_type,
                             cut_type_str,
                             cell::TriangulationStrategy::classical);
    }

    template <std::floating_point T>
    cutcells::mesh::CutMesh<T> cut_vtk_mesh(std::span<const T> ls_vals, std::span<const T> points,
                                            std::span<const int> connectivity, std::span<const int> offset,
                                            std::span<const int> vtk_type,
                                            const std::string& cut_type_str,
                                            cell::TriangulationStrategy strategy)
    {
      // -----------------------------------------------------------------------
      // Fused single-pass implementation.
      //
      // All scratch buffers are thread_local so their heap allocations are
      // amortised across calls (vectors keep their capacity; the hash map keeps
      // its bucket array).  No intermediate CutCells accumulator is created,
      // and no CutCell is ever copy-assigned — the cut writes directly into a
      // reusable scratch object that is immediately merged into the output.
      // -----------------------------------------------------------------------

      // --- thread-local scratch (reused across invocations) ------------------
      thread_local std::vector<int>   tl_intersected;
      thread_local std::vector<T>     tl_vtx_buf;
      thread_local std::vector<T>     tl_ls_buf;
      thread_local cell::CutCell<T>   tl_scratch;
      thread_local std::vector<int>   tl_local_merged;
      thread_local std::unordered_map<VertexKey, int, VertexKeyHash> tl_dedup;

      // --- Pass 1: locate intersected cells ----------------------------------
      tl_intersected.clear();
      {
        const int ncells = static_cast<int>(vtk_type.size());
        for (int i = 0; i < ncells; ++i)
        {
          const int coff = offset[i];
          const cell::type ctype = cell::map_vtk_type_to_cell_type(
              static_cast<cell::vtk_types>(vtk_type[i]));
          const int nv = cell::get_num_vertices(ctype);
          tl_ls_buf.resize(nv);
          for (int j = 0; j < nv; ++j)
            tl_ls_buf[j] = ls_vals[connectivity[coff + j]];
          if (cell::classify_cell_domain<T>(std::span<const T>(tl_ls_buf.data(), nv))
                == cell::domain::intersected)
            tl_intersected.push_back(i);
        }
      }

      if (tl_intersected.empty())
        return CutMesh<T>{};

      const int n_cut = static_cast<int>(tl_intersected.size());

      // --- Prepare output mesh -----------------------------------------------
      // Heuristic reserves: each intersected cell produces ~8 sub-cells,
      // ~4 vertices each (tet / triangle).
      CutMesh<T> cm;
      cm._gdim         = 3;   // VTK meshes always use 3-D coordinates
      cm._tdim         = 0;   // filled from first non-empty cut result
      cm._num_cells    = 0;
      cm._num_vertices = 0;
      cm._types.reserve(n_cut * 8);
      cm._offset.reserve(n_cut * 8 + 1);
      cm._connectivity.reserve(n_cut * 8 * 4);
      cm._parent_map.reserve(n_cut * 8);
      cm._vertex_coords.reserve(n_cut * 12);
      cm._offset.push_back(0);

      // --- Clear dedup map (bucket array is retained across calls) -----------
      tl_dedup.clear();
      tl_dedup.reserve(static_cast<std::size_t>(n_cut * 10));

      int total_conn  = 0;
      int total_cells = 0;

      // --- Pass 2: cut each cell and immediately merge into cm ---------------
      for (int i = 0; i < n_cut; ++i)
      {
        const int ci   = tl_intersected[i];
        const int coff = offset[ci];
        const cell::type ctype = cell::map_vtk_type_to_cell_type(
            static_cast<cell::vtk_types>(vtk_type[ci]));
        const int nv = cell::get_num_vertices(ctype);

        // Gather vertex coords and level-set values
        tl_vtx_buf.resize(nv * 3);
        tl_ls_buf.resize(nv);
        for (int j = 0; j < nv; ++j)
        {
          const int vid = connectivity[coff + j];
          tl_vtx_buf[j * 3 + 0] = points[vid * 3 + 0];
          tl_vtx_buf[j * 3 + 1] = points[vid * 3 + 1];
          tl_vtx_buf[j * 3 + 2] = points[vid * 3 + 2];
          tl_ls_buf[j] = ls_vals[vid];
        }

        // Cut into scratch — all CutCell fields are reset by the cutter
        cell::cut<T>(ctype, tl_vtx_buf, 3, tl_ls_buf, cut_type_str, tl_scratch, strategy);

        // Override parent vertex IDs with context-global mesh indices
        // (the individual cutters default them to local 0..nv-1)
        tl_scratch._parent_cell_type = ctype;
        tl_scratch._parent_vertex_ids.resize(nv);
        for (int j = 0; j < nv; ++j)
          tl_scratch._parent_vertex_ids[j] = connectivity[coff + j];

        const int n_local_verts = static_cast<int>(tl_scratch._vertex_coords.size())
                                  / tl_scratch._gdim;
        const int n_local_cells = cell::num_cells(tl_scratch);
        if (n_local_cells == 0)
          continue;

        // Latch output mesh dimensions from the first non-empty result
        if (total_cells == 0)
          cm._tdim = tl_scratch._tdim;

        // Build local → global vertex id map
        tl_local_merged.assign(n_local_verts, -1);

        const bool has_tokens = (static_cast<int>(tl_scratch._vertex_parent_entity.size())
                                 == n_local_verts);

        if (has_tokens)
        {
          // Fast token-based dedup: shared vertices between neighbouring cut
          // cells (parent vertices and edge intersections) are identified via
          // their VertexKey and mapped to a single global vertex.
          const auto parent_edges = cell::edges(ctype);

          for (int lv = 0; lv < n_local_verts; ++lv)
          {
            const int32_t token = tl_scratch._vertex_parent_entity[lv];
            VertexKey key;

            if (token >= 100 && token < 200)
            {
              const int ref = token - 100;
              if (ref >= 0 && ref < nv)
              {
                key.kind = 0;
                key.a    = static_cast<int32_t>(tl_scratch._parent_vertex_ids[ref]);
                key.b    = 0;
              }
              else { key.kind = 2; key.a = ci; key.b = token; }
            }
            else if (token >= 200)
            {
              key.kind = 2;
              key.a    = ci;
              key.b    = static_cast<int32_t>(token - 200);
            }
            else
            {
              const int eid = static_cast<int>(token);
              if (eid >= 0 && eid < static_cast<int>(parent_edges.size()))
              {
                // token is a VTK edge index; cell_topology.h uses Basix ordering.
                const int basix_eid = cell::vtk_to_basix_edge(ctype, eid);
                const int bv0 = parent_edges[basix_eid][0];
                const int bv1 = parent_edges[basix_eid][1];
                const int vtk_v0 = cell::basix_to_vtk_vertex(ctype, bv0);
                const int vtk_v1 = cell::basix_to_vtk_vertex(ctype, bv1);
                const int32_t gv0 = static_cast<int32_t>(
                    tl_scratch._parent_vertex_ids[vtk_v0]);
                const int32_t gv1 = static_cast<int32_t>(
                    tl_scratch._parent_vertex_ids[vtk_v1]);
                key.kind = 1;
                key.a    = std::min(gv0, gv1);
                key.b    = std::max(gv0, gv1);
              }
              else { key.kind = 2; key.a = ci; key.b = token; }
            }

            const int next_id = static_cast<int>(cm._vertex_coords.size()) / cm._gdim;
            auto [it, inserted] = tl_dedup.emplace(key, next_id);
            tl_local_merged[lv] = it->second;
            if (inserted)
            {
              const int base = lv * tl_scratch._gdim;
              cm._vertex_coords.insert(cm._vertex_coords.end(),
                                       tl_scratch._vertex_coords.data() + base,
                                       tl_scratch._vertex_coords.data() + base + tl_scratch._gdim);
            }
          }
        }
        else
        {
          // Fallback: tokens not available — emit all local vertices without
          // inter-cell dedup (unusual path; occurs only without token support)
          for (int lv = 0; lv < n_local_verts; ++lv)
          {
            tl_local_merged[lv] = static_cast<int>(cm._vertex_coords.size()) / cm._gdim;
            const int base = lv * tl_scratch._gdim;
            cm._vertex_coords.insert(cm._vertex_coords.end(),
                                     tl_scratch._vertex_coords.data() + base,
                                     tl_scratch._vertex_coords.data() + base + tl_scratch._gdim);
          }
        }

        // Emit sub-cells
        for (int lc = 0; lc < n_local_cells; ++lc)
        {
          const auto verts = cell::cell_vertices(tl_scratch, lc);
          const int  vcnt  = static_cast<int>(verts.size());
          for (int v = 0; v < vcnt; ++v)
            cm._connectivity.push_back(tl_local_merged[verts[v]]);
          total_conn += vcnt;
          cm._offset.push_back(total_conn);
          cm._types.push_back(tl_scratch._types[lc]);
          cm._parent_map.push_back(ci);
          ++total_cells;
        }
      }

      cm._num_cells    = total_cells;
      cm._num_vertices = static_cast<int>(cm._vertex_coords.size()) / cm._gdim;

      return cm;
    }

    template <std::floating_point T>
    cutcells::mesh::CutMesh<T> cut_vtk_mesh(std::span<const T> ls_vals, std::span<const T> points,
                                            std::span<const int> connectivity, std::span<const int> offset,
                                            std::span<const int> vtk_type,
                                            const std::string& cut_type_str,
                                            bool triangulate)
    {
      return cut_vtk_mesh<T>(
          ls_vals, points, connectivity, offset, vtk_type, cut_type_str,
          cell::triangulation_strategy_from_bool(triangulate));
    }

//-----------------------------------------------------------------------------
    template void str(const CutCells<double> &cut_mesh);
    template void str(const CutCells<float> &cut_mesh);

    template cutcells::mesh::CutMesh<double> create_cut_mesh(CutCells<double>& cut_cells);
    template cutcells::mesh::CutMesh<float> create_cut_mesh(CutCells<float>& cut_cells);

    template
    std::vector<int> locate_cells(std::span<const double> ls_vals, std::span<const double> points,
                                    std::span<const int> connectivity, std::span<const int> offset,
                                    std::span<const int> vtk_type,
                                    cell::cut_type ctype);

    template
    std::vector<int> locate_cells(std::span<const float> ls_vals, std::span<const float> points,
                                    std::span<const int> connectivity, std::span<const int> offset,
                                    std::span<const int> vtk_type,
                                    cell::cut_type ctype);

    template
    cutcells::mesh::CutMesh<double> cut_vtk_mesh(std::span<const double> ls_vals, std::span<const double> points,
                                            std::span<const int> connectivity, std::span<const int> offset,
                                            std::span<const int> vtk_type,
                                            const std::string& cut_type_str);
    template
    cutcells::mesh::CutMesh<float> cut_vtk_mesh(std::span<const float> ls_vals, std::span<const float> points,
                                            std::span<const int> connectivity, std::span<const int> offset,
                                            std::span<const int> vtk_type,
                                            const std::string& cut_type_str);

    template
    cutcells::mesh::CutMesh<double> cut_vtk_mesh(std::span<const double> ls_vals, std::span<const double> points,
                        std::span<const int> connectivity, std::span<const int> offset,
                        std::span<const int> vtk_type,
                        const std::string& cut_type_str,
                        cell::TriangulationStrategy strategy);
    template
    cutcells::mesh::CutMesh<float> cut_vtk_mesh(std::span<const float> ls_vals, std::span<const float> points,
                        std::span<const int> connectivity, std::span<const int> offset,
                        std::span<const int> vtk_type,
                        const std::string& cut_type_str,
                        cell::TriangulationStrategy strategy);

    template
    cutcells::mesh::CutMesh<double> cut_vtk_mesh(std::span<const double> ls_vals, std::span<const double> points,
                        std::span<const int> connectivity, std::span<const int> offset,
                        std::span<const int> vtk_type,
                        const std::string& cut_type_str,
                        bool triangulate);
    template
    cutcells::mesh::CutMesh<float> cut_vtk_mesh(std::span<const float> ls_vals, std::span<const float> points,
                        std::span<const int> connectivity, std::span<const int> offset,
                        std::span<const int> vtk_type,
                        const std::string& cut_type_str,
                        bool triangulate);

//-----------------------------------------------------------------------------
}
