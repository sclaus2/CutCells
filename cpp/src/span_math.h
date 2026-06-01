// Copyright (c) 2022-2024 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <array>
#include <span>
#include <string>
#include <sstream>
#include <cassert>
#include <concepts>

namespace cutcells
{
    namespace math
    {
      template<typename T, std::size_t N>
      using Vec = std::array<T, N>;

      template<typename T>
      using Vec2 = std::array<T, 2>;

      template<typename T>
      using Vec3 = std::array<T, 3>;

      template <std::size_t N, std::floating_point T>
      inline Vec<T, N> to_vec(std::span<const T> a)
      {
        assert(a.size() == N);
        Vec<T, N> result{};
        for (std::size_t i = 0; i < N; ++i)
          result[i] = a[i];
        return result;
      }

      template <std::floating_point T>
      inline Vec2<T> to_vec2(std::span<const T> a)
      {
        return to_vec<2, T>(a);
      }

      template <std::floating_point T>
      inline Vec3<T> to_vec3(std::span<const T> a)
      {
        return to_vec<3, T>(a);
      }

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
      inline Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b)
      {
        Vec3<T> result{};

        result[0] = a[1]*b[2]-a[2]*b[1];
        result[1] = a[2]*b[0]-a[0]*b[2];
        result[2] = a[0]*b[1]-a[1]*b[0];

        return result;
      }

      template <std::floating_point T>
      inline T cross(const Vec2<T>& a, const Vec2<T>& b)
      {
        return a[0] * b[1] - a[1] * b[0];
      }

      template <std::size_t N, std::floating_point T>
      inline Vec<T, N> subtract(const Vec<T, N>& a, const Vec<T, N>& b)
      {
        Vec<T, N> result{};
        for (std::size_t i = 0; i < N; ++i)
          result[i] = a[i] - b[i];
        return result;
      }

      template <std::floating_point T>
      inline Vec3<T> subtract(const Vec3<T>& a, const Vec3<T>& b)
      {
        return subtract<3, T>(a, b);
      }

      template <std::size_t N, std::floating_point T>
      inline Vec<T, N> add(const Vec<T, N>& a, const Vec<T, N>& b)
      {
        Vec<T, N> result{};
        for (std::size_t i = 0; i < N; ++i)
          result[i] = a[i] + b[i];
        return result;
      }

      template <std::floating_point T>
      inline Vec3<T> add(const Vec3<T>& a, const Vec3<T>& b)
      {
        return add<3, T>(a, b);
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