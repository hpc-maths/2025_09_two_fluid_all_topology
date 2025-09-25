// Copyright 2021 SAMURAI TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
//
// Author: Giuseppe Orlando, 2025
//
#pragma once

#include <samurai/schemes/fv.hpp>

#include "auxiliary_structs.hpp"
#include "eos.hpp"

namespace samurai {
  /**
    * Generic class to compute the flux between a left and right state
    */
  template<class Field>
  class Flux {
  public:
    /*--- Definitions and sanity checks ---*/
    static constexpr std::size_t field_size = Field::n_comp;
    using Indices = EquationData<Field::dim>;
    static_assert(field_size == Indices::NVARS, "The number of elements in the state does not correpsond to the number of equations");
    static constexpr std::size_t output_field_size = field_size;
    static constexpr std::size_t stencil_size = 2;

    using cfg = FluxConfig<SchemeType::NonLinear, output_field_size, stencil_size, Field>;

    using Number = typename Field::value_type; /*--- Shortcut for the arithmetic type ---*/

    Flux(const EOS<Number>& EOS_phase1_,
         const EOS<Number>& EOS_phase2_); /*--- Constructor which accepts in input
                                                the equations of state of the two phases ---*/

  protected:
    const EOS<Number>& EOS_phase1; /*--- Pass it by reference because pure virtual (not so nice, maybe moving to pointers) ---*/
    const EOS<Number>& EOS_phase2; /*--- Pass it by reference because pure virtual (not so nice, maybe moving to pointers) ---*/

    FluxValue<cfg> evaluate_continuous_flux(const FluxValue<cfg>& q,
                                            const std::size_t curr_d); /*--- Evaluate the 'continuous' flux for the state q
                                                                             along direction curr_d ---*/
  };

  // Class constructor in order to be able to work with the equation of state
  //
  template<class Field>
  Flux<Field>::Flux(const EOS<Number>& EOS_phase1_,
                    const EOS<Number>& EOS_phase2_):
    EOS_phase1(EOS_phase1_), EOS_phase2(EOS_phase2_) {}

  // Evaluate the 'continuous flux' along direction 'curr_d'
  //
  template<class Field>
  FluxValue<typename Flux<Field>::cfg>
  Flux<Field>::evaluate_continuous_flux(const FluxValue<cfg>& q,
                                        const std::size_t curr_d) {
    /*--- Sanity check in terms of dimensions ---*/
    assert(curr_d < Field::dim);

    /*--- Initialize with the state ---*/
    FluxValue<cfg> res = q;

    /*--- Pre-fetch variables that will be used several times so as to exploit possible vectorization ---*/
    const auto alpha1 = q(Indices::ALPHA1_INDEX);
    const auto m1     = q(Indices::ALPHA1_RHO1_INDEX);
    const auto m1E1   = q(Indices::ALPHA1_RHO1_E1_INDEX);
    const auto m2     = q(Indices::ALPHA2_RHO2_INDEX);
    const auto m2E2   = q(Indices::ALPHA2_RHO2_E2_INDEX);

    /*--- Compute density, velocity (along the dimension) and internal energy of phase 1 ---*/
    const auto rho1   = m1/alpha1; /*--- TODO: Add treatment for vanishing volume fraction ---*/
    const auto inv_m1 = static_cast<Number>(1.0)/m1; /*--- TODO: Add treatment for vanishing volume fraction ---*/
    auto e1           = m1E1*inv_m1; /*--- TODO: Add treatment for vanishing volume fraction ---*/
    for(std::size_t d = 0; d < Field::dim; ++d) {
      e1 -= static_cast<Number>(0.5)*
            ((q(Indices::ALPHA1_RHO1_U1_INDEX + d)*inv_m1)*
             (q(Indices::ALPHA1_RHO1_U1_INDEX + d)*inv_m1)); /*--- TODO: Add treatment for vanishing volume fraction ---*/
    }
    const auto pres1  = EOS_phase1.pres_value_Rhoe(rho1, e1);
    const auto vel1_d = q(Indices::ALPHA1_RHO1_U1_INDEX + curr_d)*inv_m1; /*--- TODO: Add treatment for vanishing volume fraction ---*/

    /*--- Compute the flux for the equations "associated" to phase 1 ---*/
    res(Indices::ALPHA1_INDEX) = static_cast<Number>(0.0);
    res(Indices::ALPHA1_RHO1_INDEX) *= vel1_d;
    for(std::size_t d = 0; d < Field::dim; ++d) {
      res(Indices::ALPHA1_RHO1_U1_INDEX + d) *= vel1_d;
    }
    res(Indices::ALPHA1_RHO1_U1_INDEX + curr_d) += alpha1*pres1;
    res(Indices::ALPHA1_RHO1_E1_INDEX) *= vel1_d;
    res(Indices::ALPHA1_RHO1_E1_INDEX) += alpha1*pres1*vel1_d;

    /*--- Compute density, velocity (along the dimension) and internal energy of phase 2 ---*/
    const auto alpha2 = static_cast<Number>(1.0) - alpha1;
    const auto rho2   = m2/alpha2; /*--- TODO: Add treatment for vanishing volume fraction ---*/
    const auto inv_m2 = static_cast<Number>(1.0)/m2; /*--- TODO: Add treatment for vanishing volume fraction ---*/
    auto e2           = m2E2*inv_m2; /*--- TODO: Add treatment for vanishing volume fraction ---*/
    for(std::size_t d = 0; d < Field::dim; ++d) {
      e2 -= static_cast<Number>(0.5)*
            ((q(Indices::ALPHA2_RHO2_U2_INDEX + d)*inv_m2)*
             (q(Indices::ALPHA2_RHO2_U2_INDEX + d)*inv_m2)); /*--- TODO: Add treatment for vanishing volume fraction ---*/
    }
    const auto pres2  = EOS_phase2.pres_value_Rhoe(rho2, e2);
    const auto vel2_d = q(Indices::ALPHA2_RHO2_U2_INDEX + curr_d)*inv_m2; /*--- TODO: Add treatment for vanishing volume fraction ---*/

    /*--- Compute the flux for the equations "associated" to phase 2 ---*/
    res(Indices::ALPHA2_RHO2_INDEX) *= vel2_d;
    for(std::size_t d = 0; d < Field::dim; ++d) {
      res(Indices::ALPHA2_RHO2_U2_INDEX + d) *= vel2_d;
    }
    res(Indices::ALPHA2_RHO2_U2_INDEX + curr_d) += alpha2*pres2;
    res(Indices::ALPHA2_RHO2_E2_INDEX) *= vel2_d;
    res(Indices::ALPHA2_RHO2_E2_INDEX) += alpha2*pres2*vel2_d;

    return res;
  }

} // end namespace samurai
