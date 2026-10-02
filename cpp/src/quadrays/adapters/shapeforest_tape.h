// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "../analytic.h"
#include "../taylor.h"

/// ShapeForest tapes as analytic level sets. ShapeForest lowers a shape to
/// register code (shapeforest/ir/tape.py); its arrays reach this adapter from
/// Python (cutcells.shapeforest) or from a text file, so CutCells does not
/// depend on ShapeForest. One interpreter, evaluate, runs the tape on every
/// scalar type of taylor.h, so the tape is an algoim-style functor.
namespace cutcells::quadrays::shapeforest
{

/// Opcodes of ShapeForest's tapes (shapeforest/_cpp/opcodes.hpp, table
/// version 5) that the adapter evaluates. TAN, ATAN, ATAN2, FLOOR, MOD and the
/// geometry opcodes are refused.
enum class Op : std::uint8_t
{
    constant = 0,
    input_x = 1,
    input_y = 2,
    input_z = 3,
    neg = 10,
    abs = 11,
    sqrt = 12,
    rsqrt = 13,
    sin = 14,
    cos = 15,
    exp = 18,
    log = 19,
    add = 30,
    sub = 31,
    mul = 32,
    div = 33,
    min = 34,
    max = 35,
    pow = 36,
    length3 = 39,
    hypot2 = 40,
    fma = 41,
    min3 = 42,
    max3 = 43,
    clamp = 44,
    clamp_smooth = 45,
    select_positive = 46,
};

/// A tape: instruction i writes register out[i] from registers a[i], b[i],
/// c[i] (-1 if unused) and the immediate imm[i]. Registers inputs[0..2] start
/// as x, y and z, registers extra_registers as extra_values (ShapeForest's
/// parameters); the result is register output.
struct Tape
{
    std::vector<std::uint8_t> op;
    std::vector<std::int32_t> a, b, c, out;
    std::vector<double> imm;
    int n_registers = 0;
    int output = -1;
    std::array<int, 3> inputs = {-1, -1, -1};
    std::vector<std::int32_t> extra_registers;
    std::vector<double> extra_values;
    /// Integer exponent of each POW instruction, filled by prepare_tape.
    std::vector<long> exponents;

    int n_instructions() const { return static_cast<int>(op.size()); }
};

/// @brief Check @p tape and resolve the exponents of its POW instructions,
/// which must be integer constants.
/// @throws std::invalid_argument on unsupported opcodes or registers out of range
inline void prepare_tape(Tape& tape)
{
    const std::size_t n = tape.op.size();
    if (tape.a.size() != n || tape.b.size() != n || tape.c.size() != n || tape.out.size() != n
        || tape.imm.size() != n)
        throw std::invalid_argument("shapeforest tape: instruction arrays differ in length");
    if (tape.extra_registers.size() != tape.extra_values.size())
        throw std::invalid_argument("shapeforest tape: extra registers and values differ in length");
    const int nr = tape.n_registers;
    auto check = [nr](int r, bool used)
    {
        if (used && (r < 0 || r >= nr))
            throw std::invalid_argument("shapeforest tape: register out of range");
    };
    check(tape.output, true);
    for (const int r : tape.inputs)
        check(r, r >= 0);
    for (const std::int32_t r : tape.extra_registers)
        check(r, true);

    // constant value per register as the tape runs, for exponents
    const double unknown = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> constant(static_cast<std::size_t>(std::max(nr, 0)), unknown);
    for (std::size_t e = 0; e < tape.extra_registers.size(); ++e)
        constant[tape.extra_registers[e]] = tape.extra_values[e];
    tape.exponents.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i)
    {
        int arity = 0;
        switch (static_cast<Op>(tape.op[i]))
        {
        case Op::constant:
        case Op::input_x:
        case Op::input_y:
        case Op::input_z:
            break;
        case Op::neg:
        case Op::abs:
        case Op::sqrt:
        case Op::rsqrt:
        case Op::sin:
        case Op::cos:
        case Op::exp:
        case Op::log:
        case Op::clamp_smooth:
            arity = 1;
            break;
        case Op::add:
        case Op::sub:
        case Op::mul:
        case Op::div:
        case Op::min:
        case Op::max:
        case Op::hypot2:
            arity = 2;
            break;
        case Op::pow:
        {
            arity = 1;
            double p = tape.imm[i];
            if (tape.b[i] >= 0)
            {
                check(tape.b[i], true);
                p = constant[tape.b[i]];
            }
            if (!(std::floor(p) == p) || std::abs(p) > 1024)
                throw std::invalid_argument("shapeforest tape: POW needs a constant integer exponent");
            tape.exponents[i] = static_cast<long>(p);
            break;
        }
        case Op::length3:
        case Op::fma:
        case Op::min3:
        case Op::max3:
        case Op::clamp:
        case Op::select_positive:
            arity = 3;
            break;
        default:
            throw std::invalid_argument("shapeforest tape: unsupported opcode "
                                        + std::to_string(static_cast<int>(tape.op[i])));
        }
        check(tape.a[i], arity >= 1);
        check(tape.b[i], arity >= 2);
        check(tape.c[i], arity >= 3);
        check(tape.out[i], true);
        constant[tape.out[i]] = static_cast<Op>(tape.op[i]) == Op::constant ? tape.imm[i] : unknown;
    }
}

/// @brief Read a tape from the text format of tools/export_tape.py of the
/// quadrature prototype: "shapeforest-tape 1", then "n_regs out x y z n", then
/// one line "op a b c out imm" per instruction.
inline Tape read_tape(const std::string& path)
{
    std::ifstream in(path);
    std::string magic;
    int version = 0;
    if (!(in >> magic >> version) || magic != "shapeforest-tape" || version != 1)
        throw std::runtime_error("read_tape: not a shapeforest-tape file: " + path);
    Tape t;
    int n = 0;
    in >> t.n_registers >> t.output >> t.inputs[0] >> t.inputs[1] >> t.inputs[2] >> n;
    if (!in || n < 0)
        throw std::runtime_error("read_tape: bad header: " + path);
    t.op.resize(n);
    t.a.resize(n);
    t.b.resize(n);
    t.c.resize(n);
    t.out.resize(n);
    t.imm.resize(n);
    for (int i = 0; i < n; ++i)
    {
        int op = 0;
        in >> op >> t.a[i] >> t.b[i] >> t.c[i] >> t.out[i] >> t.imm[i];
        t.op[i] = static_cast<std::uint8_t>(op);
    }
    if (!in)
        throw std::runtime_error("read_tape: truncated file: " + path);
    prepare_tape(t);
    return t;
}

/// @brief Run a prepared tape at @p x on the registers @p r (resized as
/// needed). V is double, Dual<double, 3> or Dual<Taylor<double, N>, N>;
/// Taylor models throw std::domain_error where a bound is lost.
template <typename V>
V evaluate(const Tape& t, const std::array<V, 3>& x, std::vector<V>& r)
{
    using std::abs;
    using std::cos;
    using std::exp;
    using std::log;
    using std::max;
    using std::min;
    using std::sin;
    using std::sqrt;
    // a square root of a value rounded below 0 is 0, as in ShapeForest
    auto root = [](const V& v)
    {
        if constexpr (std::is_floating_point_v<V>)
            return sqrt(max(v, V(0.0)));
        else
            return sqrt(v);
    };
    r.assign(static_cast<std::size_t>(t.n_registers), V(0.0));
    for (int i = 0; i < 3; ++i)
        if (t.inputs[i] >= 0)
            r[t.inputs[i]] = x[i];
    for (std::size_t e = 0; e < t.extra_registers.size(); ++e)
        r[t.extra_registers[e]] = V(t.extra_values[e]);
    const std::size_t n = t.op.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::int32_t A = t.a[i], B = t.b[i], C = t.c[i];
        V v;
        switch (static_cast<Op>(t.op[i]))
        {
        case Op::constant:
            v = V(t.imm[i]);
            break;
        case Op::input_x:
            v = x[0];
            break;
        case Op::input_y:
            v = x[1];
            break;
        case Op::input_z:
            v = x[2];
            break;
        case Op::neg:
            v = -r[A];
            break;
        case Op::abs:
            v = abs(r[A]);
            break;
        case Op::sqrt:
            v = root(r[A]);
            break;
        case Op::rsqrt:
            v = V(1.0) / sqrt(r[A]);
            break;
        case Op::sin:
            v = sin(r[A]);
            break;
        case Op::cos:
            v = cos(r[A]);
            break;
        case Op::exp:
            v = exp(r[A]);
            break;
        case Op::log:
            v = log(r[A]);
            break;
        case Op::add:
            v = r[A] + r[B];
            break;
        case Op::sub:
            v = r[A] - r[B];
            break;
        case Op::mul:
            v = r[A] * r[B];
            break;
        case Op::div:
            v = r[A] / r[B];
            break;
        case Op::min:
            v = min(r[A], r[B]);
            break;
        case Op::max:
            v = max(r[A], r[B]);
            break;
        case Op::pow:
            v = integer_power(r[A], t.exponents[i]);
            break;
        case Op::length3:
            v = root(r[A] * r[A] + r[B] * r[B] + r[C] * r[C]);
            break;
        case Op::hypot2:
            v = root(r[A] * r[A] + r[B] * r[B]);
            break;
        case Op::fma:
            v = r[A] * r[B] + r[C];
            break;
        case Op::min3:
            v = min(r[A], min(r[B], r[C]));
            break;
        case Op::max3:
            v = max(r[A], max(r[B], r[C]));
            break;
        case Op::clamp:
            v = min(max(r[A], r[B]), r[C]);
            break;
        case Op::clamp_smooth:
            v = min(max(r[A], V(0.0)), V(1.0));
            break;
        case Op::select_positive:
        {
            const int s = certain_sign(r[A]);
            v = s > 0 ? r[B] : (s < 0 ? r[C] : hull(r[B], r[C]));
            break;
        }
        default:
            throw std::runtime_error("shapeforest tape: unsupported opcode; call prepare_tape first");
        }
        r[t.out[i]] = v;
    }
    return r[t.output];
}

/// A prepared tape as an algoim-style functor, for analytic_level_set.
struct TapeLevelSet
{
    const Tape* tape = nullptr;

    template <typename V>
    V operator()(const std::array<V, 3>& x) const
    {
        thread_local std::vector<V> registers;
        return evaluate(*tape, x, registers);
    }
};

} // namespace cutcells::quadrays::shapeforest
