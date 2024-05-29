// Copyright (c) 2022-2024 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <span>
#include <string>
#include <sstream>
#include <cassert>

namespace cutcells
{
    namespace math
    {
      inline double dot(std::span<const double> a, std::span<const double> b)
      {
        assert(a.size()==b.size());
        double result = 0;

        for(std::size_t i=0;i<a.size();i++)
          result +=a[i]*b[i];

        return result;
      }

      inline std::vector<double> cross(std::span<const double> a, std::span<const double> b)
      {
        //cross product only implemented in 3D
        assert(a.size()==3);
        assert(a.size()==b.size());
        std::vector<double> result(3);

        result[0] = a[1]*b[2]-a[2]*b[1];
        result[1] = a[2]*b[0]-a[0]*b[2];
        result[2] = a[0]*b[1]-a[1]*b[0];

        return result;
      }

      inline std::vector<double> subtract(std::span<const double> a, std::span<const double> b)
      {
        assert(a.size()==b.size());
        std::vector<double> result(a.size());
        for(std::size_t i=0;i<a.size();i++)
          result[i] =a[i]-b[i];
        return result;
      }

      inline std::vector<double> add(std::span<const double> a, std::span<const double> b)
      {
        assert(a.size()==b.size());
        std::vector<double> result(a.size());
        for(std::size_t i=0;i<a.size();i++)
          result[i] =a[i]+b[i];
        return result;
      }

      inline double distance(std::span<const double> a, std::span<const double> b)
      {
        assert(a.size()==b.size());
        double result = 0;
        for(std::size_t i=0;i<a.size();i++)
          result +=(b[i]-a[i])*(b[i]-a[i]);

        result = sqrt(result);

        //account for numerical round-off error
        if(result < 1e-14)
          result = 0.0;

        return result;
      }

      static std::string print(std::span<const double> a)
      {
        std::ostringstream str_out;

        for(std::size_t i=0; i< a.size();i++)
        {
          str_out << a[i] << ", ";
        }
        str_out << std::endl;
        return str_out.str();
      }
    }
}