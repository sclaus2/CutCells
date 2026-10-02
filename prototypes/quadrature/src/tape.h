// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cutcells::proto
{

/// One instruction of a ShapeForest tape (shapeforest/ir/tape.py); opcodes as in
/// shapeforest/_cpp/opcodes.hpp.
struct TapeInstr
{
    int op = 0;
    int a = -1, b = -1, c = -1, out = -1;
    double imm = 0;
};

/// A level set given by a ShapeForest tape: registers, one output, inputs x, y, z.
struct Tape
{
    std::vector<TapeInstr> code;
    int n_regs = 0;
    int out = -1;
    int x = -1, y = -1, z = -1;
};

/// Read the text format written by tools/export_tape.py.
inline Tape read_tape(const std::string& path)
{
    std::ifstream in(path);
    std::string magic;
    int version = 0;
    if (!(in >> magic >> version) || magic != "shapeforest-tape" || version != 1)
        throw std::runtime_error("read_tape: not a shapeforest-tape file: " + path);
    Tape t;
    int n = 0;
    in >> t.n_regs >> t.out >> t.x >> t.y >> t.z >> n;
    t.code.resize(n);
    for (TapeInstr& i : t.code)
        in >> i.op >> i.a >> i.b >> i.c >> i.out >> i.imm;
    if (!in)
        throw std::runtime_error("read_tape: truncated file: " + path);
    return t;
}

// ============================================================================
// Scalar types
// ============================================================================
//
// evaluate<T> runs the tape on any T with arithmetic, sqrt, sin, cos, exp, log
// and the branch helpers below: double for values, Dual<T, N> for derivatives,
// algoim::Interval<N> (tape_algoim.h) for bounds over a box, Dual of an interval
// for bounds of the gradient. Like algoim's templated level-set functors.

/// Sign that holds everywhere T can take values: +1, -1, or 0 if not certain.
inline int certain_sign(double v) { return (v > 0) - (v < 0); }
/// True if a < b everywhere.
inline bool certainly_less(double a, double b) { return a < b; }
/// A T containing both a and b (for double only reached on ties, a == b).
inline double hull(double a, double) { return a; }
inline double tabs(double a) { return std::abs(a); }
inline double tsqrt(double a) { return std::sqrt(a); }
inline double tmin(double a, double b) { return std::min(a, b); }
inline double tmax(double a, double b) { return std::max(a, b); }

/// Value and first derivatives with respect to N variables (forward mode).
template <typename T, int N>
struct Dual
{
    T v;
    std::array<T, N> d;

    Dual() : v(0.0) { d.fill(T(0.0)); }
    explicit Dual(double c) : v(c) { d.fill(T(0.0)); }
};

template <typename T, int N>
Dual<T, N> operator-(const Dual<T, N>& a)
{
    Dual<T, N> r;
    r.v = -a.v;
    for (int k = 0; k < N; ++k)
        r.d[k] = -a.d[k];
    return r;
}

template <typename T, int N>
Dual<T, N> operator+(const Dual<T, N>& a, const Dual<T, N>& b)
{
    Dual<T, N> r;
    r.v = a.v + b.v;
    for (int k = 0; k < N; ++k)
        r.d[k] = a.d[k] + b.d[k];
    return r;
}

template <typename T, int N>
Dual<T, N> operator-(const Dual<T, N>& a, const Dual<T, N>& b)
{
    Dual<T, N> r;
    r.v = a.v - b.v;
    for (int k = 0; k < N; ++k)
        r.d[k] = a.d[k] - b.d[k];
    return r;
}

template <typename T, int N>
Dual<T, N> operator*(const Dual<T, N>& a, const Dual<T, N>& b)
{
    Dual<T, N> r;
    r.v = a.v * b.v;
    for (int k = 0; k < N; ++k)
        r.d[k] = a.v * b.d[k] + b.v * a.d[k];
    return r;
}

template <typename T, int N>
Dual<T, N> operator/(const Dual<T, N>& a, const Dual<T, N>& b)
{
    Dual<T, N> r;
    r.v = a.v / b.v;
    for (int k = 0; k < N; ++k)
        r.d[k] = (a.d[k] - r.v * b.d[k]) / b.v;
    return r;
}

template <typename T, int N>
Dual<T, N> tsqrt(const Dual<T, N>& a)
{
    Dual<T, N> r;
    r.v = tsqrt(a.v);
    const T twice = T(2.0) * r.v;
    for (int k = 0; k < N; ++k)
        r.d[k] = a.d[k] / twice;
    return r;
}

template <typename T, int N>
Dual<T, N> sin(const Dual<T, N>& a)
{
    using std::cos;
    using std::sin;
    Dual<T, N> r;
    r.v = sin(a.v);
    const T c = cos(a.v);
    for (int k = 0; k < N; ++k)
        r.d[k] = c * a.d[k];
    return r;
}

template <typename T, int N>
Dual<T, N> cos(const Dual<T, N>& a)
{
    using std::cos;
    using std::sin;
    Dual<T, N> r;
    r.v = cos(a.v);
    const T s = -sin(a.v);
    for (int k = 0; k < N; ++k)
        r.d[k] = s * a.d[k];
    return r;
}

template <typename T, int N>
Dual<T, N> exp(const Dual<T, N>& a)
{
    using std::exp;
    Dual<T, N> r;
    r.v = exp(a.v);
    for (int k = 0; k < N; ++k)
        r.d[k] = r.v * a.d[k];
    return r;
}

template <typename T, int N>
Dual<T, N> log(const Dual<T, N>& a)
{
    using std::log;
    Dual<T, N> r;
    r.v = log(a.v);
    for (int k = 0; k < N; ++k)
        r.d[k] = a.d[k] / a.v;
    return r;
}

template <typename T, int N>
int certain_sign(const Dual<T, N>& a)
{
    return certain_sign(a.v);
}

template <typename T, int N>
bool certainly_less(const Dual<T, N>& a, const Dual<T, N>& b)
{
    return certainly_less(a.v, b.v);
}

template <typename T, int N>
Dual<T, N> hull(const Dual<T, N>& a, const Dual<T, N>& b)
{
    Dual<T, N> r;
    r.v = hull(a.v, b.v);
    for (int k = 0; k < N; ++k)
        r.d[k] = hull(a.d[k], b.d[k]);
    return r;
}

/// |a|; where the sign is not certain, the derivative is enclosed by both branches.
template <typename T, int N>
Dual<T, N> tabs(const Dual<T, N>& a)
{
    const int s = certain_sign(a.v);
    if (s > 0)
        return a;
    if (s < 0)
        return -a;
    Dual<T, N> r = hull(a, -a);
    r.v = tabs(a.v);
    return r;
}

template <typename T, int N>
Dual<T, N> tmin(const Dual<T, N>& a, const Dual<T, N>& b)
{
    if (certainly_less(a.v, b.v))
        return a;
    if (certainly_less(b.v, a.v))
        return b;
    Dual<T, N> r = hull(a, b);
    r.v = tmin(a.v, b.v);
    return r;
}

template <typename T, int N>
Dual<T, N> tmax(const Dual<T, N>& a, const Dual<T, N>& b)
{
    if (certainly_less(a.v, b.v))
        return b;
    if (certainly_less(b.v, a.v))
        return a;
    Dual<T, N> r = hull(a, b);
    r.v = tmax(a.v, b.v);
    return r;
}

/// a^p for an integer p.
template <typename T>
T integer_power(const T& a, long p)
{
    T result(1.0), base = a;
    for (long e = p < 0 ? -p : p; e > 0; e >>= 1)
    {
        if (e & 1)
            result = result * base;
        base = base * base;
    }
    return p < 0 ? T(1.0) / result : result;
}

// ============================================================================
// Evaluation
// ============================================================================

/// Run the tape at x on registers r (resized as needed); returns the output.
/// Throws std::runtime_error on opcodes it does not support, and whatever T
/// throws (algoim's intervals throw std::domain_error where a bound is lost).
/// Square roots go through tsqrt: algoim's own sqrt for its intervals leaves out
/// part of the remainder (tape_algoim.h).
template <typename T>
T evaluate(const Tape& t, const std::array<T, 3>& x, std::vector<T>& r)
{
    using std::cos;
    using std::exp;
    using std::log;
    using std::sin;
    r.assign(static_cast<std::size_t>(t.n_regs), T(0.0));
    // constant values per register, for exponents given in a register
    thread_local std::vector<double> constant;
    constant.assign(static_cast<std::size_t>(t.n_regs), std::numeric_limits<double>::quiet_NaN());
    for (const TapeInstr& in : t.code)
    {
        T& o = r[in.out];
        switch (in.op)
        {
        case 0: // CONST
            o = T(in.imm);
            constant[in.out] = in.imm;
            break;
        case 1: // INPUT_X
            o = x[0];
            break;
        case 2: // INPUT_Y
            o = x[1];
            break;
        case 3: // INPUT_Z
            o = x[2];
            break;
        case 10: // NEG
            o = -r[in.a];
            break;
        case 11: // ABS
            o = tabs(r[in.a]);
            break;
        case 12: // SQRT
            o = tsqrt(r[in.a]);
            break;
        case 13: // RSQRT
            o = T(1.0) / tsqrt(r[in.a]);
            break;
        case 14: // SIN
            o = sin(r[in.a]);
            break;
        case 15: // COS
            o = cos(r[in.a]);
            break;
        case 18: // EXP
            o = exp(r[in.a]);
            break;
        case 19: // LOG
            o = log(r[in.a]);
            break;
        case 30: // ADD
            o = r[in.a] + r[in.b];
            break;
        case 31: // SUB
            o = r[in.a] - r[in.b];
            break;
        case 32: // MUL
            o = r[in.a] * r[in.b];
            break;
        case 33: // DIV
            o = r[in.a] / r[in.b];
            break;
        case 34: // MIN
            o = tmin(r[in.a], r[in.b]);
            break;
        case 35: // MAX
            o = tmax(r[in.a], r[in.b]);
            break;
        case 36: // POW
        {
            const double p = in.b >= 0 ? constant[in.b] : in.imm;
            if (!(std::floor(p) == p))
                throw std::runtime_error("evaluate: only integer constant exponents are supported");
            o = integer_power(r[in.a], static_cast<long>(p));
            break;
        }
        case 39: // LENGTH3
            o = tsqrt(r[in.a] * r[in.a] + r[in.b] * r[in.b] + r[in.c] * r[in.c]);
            break;
        case 40: // HYPOT2
            o = tsqrt(r[in.a] * r[in.a] + r[in.b] * r[in.b]);
            break;
        case 41: // FMA
            o = r[in.a] * r[in.b] + r[in.c];
            break;
        case 42: // MIN3
            o = tmin(r[in.a], tmin(r[in.b], r[in.c]));
            break;
        case 43: // MAX3
            o = tmax(r[in.a], tmax(r[in.b], r[in.c]));
            break;
        default:
            throw std::runtime_error("evaluate: unsupported tape opcode " + std::to_string(in.op));
        }
    }
    return r[t.out];
}

/// Value and gradient at a point.
inline double value_gradient(const Tape& t, const std::array<double, 3>& x, std::array<double, 3>& g)
{
    thread_local std::vector<Dual<double, 3>> regs;
    std::array<Dual<double, 3>, 3> in;
    for (int i = 0; i < 3; ++i)
    {
        in[i] = Dual<double, 3>(x[i]);
        in[i].d[i] = 1.0;
    }
    const Dual<double, 3> r = evaluate(t, in, regs);
    g = r.d;
    return r.v;
}

inline double tape_value(const Tape& t, const std::array<double, 3>& x)
{
    thread_local std::vector<double> regs;
    return evaluate(t, x, regs);
}

} // namespace cutcells::proto
