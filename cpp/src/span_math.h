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
#include <concepts>

namespace cutcells
{
    namespace math
    {
      template <std::floating_point T>
      inline T dot(std::span<const T> a, std::span<const T> b)
      {
        assert(a.size()==b.size());
        T result = 0;

        for(std::size_t i=0;i<a.size();i++)
          result +=a[i]*b[i];

        return result;
      }

      template <std::floating_point T>
      inline std::vector<T> cross(std::span<const T> a, std::span<const T> b)
      {
        //cross product only implemented in 3D
        assert(a.size()==3);
        assert(a.size()==b.size());
        std::vector<T> result(3);

        result[0] = a[1]*b[2]-a[2]*b[1];
        result[1] = a[2]*b[0]-a[0]*b[2];
        result[2] = a[0]*b[1]-a[1]*b[0];

        return result;
      }

      template <std::floating_point T>
      inline std::vector<T> subtract(std::span<const T> a, std::span<const T> b)
      {
        assert(a.size()==b.size());
        std::vector<T> result(a.size());
        for(std::size_t i=0;i<a.size();i++)
          result[i] =a[i]-b[i];
        return result;
      }

      template <std::floating_point T>
      inline std::vector<T> add(std::span<const T> a, std::span<const T> b)
      {
        assert(a.size()==b.size());
        std::vector<T> result(a.size());
        for(std::size_t i=0;i<a.size();i++)
          result[i] =a[i]+b[i];
        return result;
      }

      template <std::floating_point T>
      inline T distance(std::span<const T> a, std::span<const T> b)
      {
        assert(a.size()==b.size());
        T result = 0;
        for(std::size_t i=0;i<a.size();i++)
          result +=(b[i]-a[i])*(b[i]-a[i]);

        result = sqrt(result);

        //account for numerical round-off error
        if(result < 1e-14)
          result = 0.0;

        return result;
      }

      template <std::floating_point T>
      static std::string print(std::span<const T> a)
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