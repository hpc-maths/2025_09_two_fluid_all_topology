// Copyright 2021 SAMURAI TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
//
// Author: Giuseppe Orlando, 2025
//
#pragma once

#include <samurai/bc.hpp>

#include "flux_base.hpp"

// Default boundary condition
//
template<class Field>
struct Default: public samurai::Bc<Field> {
  INIT_BC(Default, samurai::Flux<Field>::stencil_size)

  inline stencil_t get_stencil(constant_stencil_size_t) const override {
    return samurai::line_stencil_from<Field::dim, 0, samurai::Flux<Field>::stencil_size>(0);
  }

  inline apply_function_t get_apply_function(constant_stencil_size_t, const direction_t&) const override {
    return [](Field& U, const stencil_cells_t& cells, const value_t& value) {
      U[cells[1]] = value;
    };
  }
};

// Reflective boundary conditions
//
template<class Field>
auto Reflecting(const Field& Q) {
  // Specify the use of this struct where we just store the indices
  using Indices = samurai::Flux<Field>::Indices;

  // Reflective bc
  return[&Q]
  (const auto& normal, const auto& cell_in, const auto& /*coord*/)
  {
    /*--- Impose reflective bc ---*/
    xt::xtensor_fixed<typename Field::value_type, xt::xshape<Field::n_comp>> Q_ghost;
    auto m1u1_dot_n = static_cast<typename Field::value_type>(0.0);
    auto m2u2_dot_n = static_cast<typename Field::value_type>(0.0);
    for(std::size_t d = 0; d < Field::dim; ++d) {
      m1u1_dot_n += Q[cell_in](Indices::ALPHA1_RHO1_U1_INDEX + d)*normal[d];
      m2u2_dot_n += Q[cell_in](Indices::ALPHA2_RHO2_U2_INDEX + d)*normal[d];
    }

    for(std::size_t d = 0; d < Field::dim; ++d) {
      Q_ghost[Indices::ALPHA1_RHO1_U1_INDEX + d] = Q[cell_in](Indices::ALPHA1_RHO1_U1_INDEX + d)
                                                 - static_cast<typename Field::value_type>(2.0)*m1u1_dot_n*normal[d];
      Q_ghost[Indices::ALPHA2_RHO2_U2_INDEX + d] = Q[cell_in](Indices::ALPHA2_RHO2_U2_INDEX + d)
                                                 - static_cast<typename Field::value_type>(2.0)*m2u2_dot_n*normal[d];
    }

    /*--- Complete the ghost state with the inner state ---*/
    Q_ghost[Indices::ALPHA1_INDEX]         = Q[cell_in][Indices::ALPHA1_INDEX];
    Q_ghost[Indices::ALPHA1_RHO1_INDEX]    = Q[cell_in][Indices::ALPHA1_RHO1_INDEX];
    Q_ghost[Indices::ALPHA2_RHO2_INDEX]    = Q[cell_in][Indices::ALPHA2_RHO2_INDEX];
    Q_ghost[Indices::ALPHA1_RHO1_E1_INDEX] = Q[cell_in][Indices::ALPHA1_RHO1_E1_INDEX];
    Q_ghost[Indices::ALPHA2_RHO2_E2_INDEX] = Q[cell_in][Indices::ALPHA2_RHO2_E2_INDEX];

    return Q_ghost;
  };
}
