# Local patches

This vendored Algoim snapshot is based on the commit recorded in
`UPSTREAM_COMMIT`.

CutCells applies the following local compatibility patch:

- `algoim/quadrature_general.hpp`: replace removed C++20 `std::result_of`
  usage with `std::invoke_result_t`, and include `<type_traits>`.
