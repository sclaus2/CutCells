// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    LGPL-3.0-or-later
#pragma once

#include <span>

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

        inline cut_type string_to_cut_type(std::string type_str)
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
        }

        inline int get_entity_flag(const std::span<const double> ls_values, bool outside)
        {
            int index = 0; 
            int multiplier = 1; 

            //todo: consider to provide insideout flag to get outside cut
            for(int i=0;i<ls_values.size();i++)
            {
                if(outside)
                {
                    if (ls_values[i] >= 0.0)
                    {
                        index |= multiplier; 
                    } 
                }
                else
                {
                    if (ls_values[i] < 0.0)
                    {
                        index |= multiplier; 
                    } 
                } 
                multiplier *=2;
            }

            return index;
        }

        inline bool is_intersected(const std::span<const double> ls_values)
        {
            // Check if there is a sign change 
            for(int i=1;i<ls_values.size();i++)
            {
                //check if sign in current node is different from first node value
                if(ls_values[0]*ls_values[i]<=0)
                {
                    return true;
                }
            }

            return false;
        }

        // determine if cell is inside, outside or intersected 
        // whether a cell is inside (phi<0), outside (phi>0) or intersected (phi changes sign) depends on if the level set values in the nodes change sign 
        // or not 
        inline domain classify_cell_domain(const std::span<const double> ls_values)
        {
            domain marker = domain::unset;

            if(is_intersected(ls_values))
            {
                marker = domain::intersected;
                return marker;
            }

            if(marker == domain::unset)
            {
                if(ls_values[0]>0)
                {
                    marker = domain::outside;
                }
                else
                {
                    marker = domain::inside;
                }
            }

            return marker;
        }

    }
}