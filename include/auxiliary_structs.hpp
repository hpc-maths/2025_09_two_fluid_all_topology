// Copyright 2021 SAMURAI TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
//
// Author: Giuseppe Orlando, 2025
//
#pragma once

// Auxiliary function to convert unsigned to string
//
template<typename T>
std::string unsigned_to_string(const T value, const unsigned digits = 5) {
  std::string lc_string = std::to_string(value);

  if(lc_string.size() < digits) {
    // We have to add the padding zeroes in front of the number
    const unsigned int padding_position = (lc_string[0] == '-') ? 1 : 0;

    const std::string padding(digits - lc_string.size(), '0');
    lc_string.insert(padding_position, padding);
  }

  return lc_string;
}

// Auxiliary struct to memorize indices
//
template<std::size_t dim>
struct EquationData {
  /*--- Declare suitable static variables for the sake of generalities in the indices ---*/
  static constexpr std::size_t ALPHA1_INDEX         = 0;
  static constexpr std::size_t ALPHA1_RHO1_INDEX    = 1;
  static constexpr std::size_t ALPHA1_RHO1_U1_INDEX = 2;
  static constexpr std::size_t ALPHA1_RHO1_E1_INDEX = ALPHA1_RHO1_U1_INDEX + dim;
  static constexpr std::size_t ALPHA2_RHO2_INDEX    = ALPHA1_RHO1_E1_INDEX + 1;
  static constexpr std::size_t ALPHA2_RHO2_U2_INDEX = ALPHA2_RHO2_INDEX + 1;
  static constexpr std::size_t ALPHA2_RHO2_E2_INDEX = ALPHA2_RHO2_U2_INDEX + dim;

  static constexpr std::size_t NVARS = ALPHA2_RHO2_E2_INDEX + 1;
};
