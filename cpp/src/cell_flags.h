// Copyright (c) 2022 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <span>
#include <cmath>
#include <limits>
#include <concepts>
#include <string>

namespace cutcells
{
    namespace cell
    {
        enum class domain
        {
            unset = -1,
            inside = 0,
            intersected = 1,
            outside = 2
        };

        enum class cut_type
        {
            unset = -1,
            philt0 = 0,
            phieq0 = 1,
            phigt0 = 2
        };

        inline cut_type string_to_cut_type(const std::string& type_str)
        {
            if(type_str=="phi<0")
                return cut_type::philt0;
            else if(type_str=="phi>0")
                return cut_type::phigt0;
            else if(type_str=="phi=0")
                return cut_type::phieq0;
            else
                return cut_type::unset;
        }

    inline std::string domain_type_to_string(const domain& cell_domain)
    {
        if(cell_domain==domain::unset)
            return "unset";
        else if(cell_domain==domain::intersected)
            return "intersected";
        else if(cell_domain==domain::inside)
            return "inside";
        else if(cell_domain==domain::outside)
            return "outside";
        return "unknown";
    }

        template <std::floating_point T>
        inline int get_entity_flag(const std::span<const T> ls_values, bool outside)
        {
            int index = 0;
            int multiplier = 1;

            int ulp = 2;
            //approximate tolerance
            T tol = std::numeric_limits<T>::epsilon() * std::fabs(ls_values[0]) * ulp;

            for(int i=0;i<ls_values.size();i++)
            {
                if(outside)
                {
                   //corner case categorized as phi<0 -> 0
                    if (ls_values[i] > tol)
                    {
                        index |= multiplier;
                    }
                }
                else
                {
                  //corner case categorized as phi>0 -> 0
                    if (ls_values[i] <= -tol) // -tol
                    {
                        index |= multiplier;
                    }
                }

                multiplier *=2;
            }

            return index;
        }

        /// @brief
        /// @param ls_values
        /// @return
        template <std::floating_point T>
        inline bool is_intersected(const std::span<const T> ls_values, const T tol = 1e-14)
        {
            // Check if there is a sign change 
            for(std::size_t i=1;i<ls_values.size();i++)
            {
                //check if sign in current node is different from first node value
                //sign change is only registered if smaller than -tol
                if(ls_values[0]*ls_values[i]<= -tol)
                {
                    return true;
                }
            }

            return false;
        }

        //catch corner cases where absolute level set values are smaller than tolerance abs(ls)<tol
        template <std::floating_point T>
        inline bool is_corner_case(const std::span<const T> ls_values, const T tol = 1e-14)
        {
            // Check if there is a sign change 
            for(std::size_t i=1;i<ls_values.size();i++)
            {
                //check if sign in current node is different from first node value
                if(std::fabs(ls_values[i])<tol)
                {
                    return true;
                }
            }

            return false;
        }

        // Flag all nodes with abs(ls_val) < tol
        template <std::floating_point T>
        inline int get_entity_corner_case_flag(const std::span<const T> ls_values, [[maybe_unused]] bool outside)
        {
            int index = 0;
            int multiplier = 1;

            int ulp = 2;
            //approximate tolerance
            T tol = std::numeric_limits<T>::epsilon() * std::fabs(ls_values[0]) * ulp;

            //todo: consider to provide insideout flag to get outside cut
            for(std::size_t i=0;i<ls_values.size();i++)
            {
                if(std::fabs(ls_values[i])<tol)
                {
                    index |= multiplier;
                }

                multiplier *=2;
            }

            return index;
        }

        // determine if cell is inside, outside or intersected
        // whether a cell is inside (phi<0), outside (phi>0) or intersected (phi changes sign) depends on if the level set values in the nodes change sign 
        // or not

        // classify domain cells using strict criteria that exclude corner cases
        // as corner cases have intersection in nodes/edges/faces
        template <std::floating_point T>
        inline domain classify_cell_domain(std::span<const T> ls_values)
        {
            if (ls_values.empty())
                return domain::unset;

            T scale = T(0);
            for (const T v : ls_values)
                scale = std::max(scale, std::abs(v));
            const T tol = std::numeric_limits<T>::epsilon()
                        * std::max(T(1), scale) * T(8);

            int n_neg = 0;
            int n_pos = 0;
            int n_zero = 0;
            for (const T v : ls_values)
            {
                if (std::abs(v) <= tol)
                    ++n_zero;
                else if (v < T(0))
                    ++n_neg;
                else
                    ++n_pos;
            }

            if (n_neg > 0 && n_pos > 0)
                return domain::intersected;
            if (n_neg > 0)
                return domain::inside;
            if (n_pos > 0)
                return domain::outside;

            // All-zero degenerate touch case: keep deterministic non-cut classification.
            return domain::inside;
        }

    }
}
