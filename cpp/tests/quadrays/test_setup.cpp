// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Set-up check: the quadrays headers are found as <cutcells/quadrays/...>,
// and a test links against the cutcells library.

#include <cutcells/quadrays/rules.h>
#include <cutcells/selection_expr.h>

#include <cstdio>
#include <string>
#include <vector>

int main()
{
    cutcells::SelectionExpr expr = cutcells::parse_selection_expr("phi < 0");
    cutcells::compile_selection_expr(expr, std::vector<std::string>{"phi"});

    if (expr.terms.size() != 1 || expr.terms[0].negative_required != 1u
        || expr.terms[0].zero_required != 0u
        || expr.terms[0].positive_required != 0u)
    {
        std::fprintf(stderr, "quadrays_setup: \"phi < 0\" compiled wrongly\n");
        return 1;
    }
    std::printf("quadrays_setup: ok\n");
    return 0;
}
