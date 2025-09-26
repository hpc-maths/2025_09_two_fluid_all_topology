// Copyright 2021 SAMURAI TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
//
// Author: Giuseppe Orlando, 2025
//
#pragma once

// Declare a struct with the simulation parameters
// (domain, levels, final time, and Courant number)
//
template<typename T = double>
struct Simulation_Parameters {
  /*--- Physical parameters ---*/
  double xL;
  double xR;
  double yL;
  double yR;

  T t0;
  T Tf;

  /*--- Numerical parameters ---*/
  T Courant;
  T dt;

  /*--- Mesh parameters ---*/
  std::size_t min_level;
  std::size_t max_level;

  /*--- Output parameters ---*/
  std::size_t nfiles;
};

// Declare a struct with EOS parameters
//
template<typename T = double>
struct EOS_Parameters {
  /*--- SG-EOS parameters phase 1 ---*/
  T gamma_1;
  T pi_infty_1;
  T q_infty_1;
  T c_v_1;

  /*--- SG-EOS parameters phase 2 ---*/
  T gamma_2;
  T pi_infty_2;
  T q_infty_2;
  T c_v_2;
};

// Declare a struct with Riemann problem parameters
//
template<typename T = double>
struct Riemann_Parameters {
  /*--- Initial discontinuity location ---*/
  T xd;
  T yd;

  /*--- North-West state ---*/
  T alpha1NW;
  T rho1NW;
  T p1NW;
  T u1NW;
  T v1NW;
  T rho2NW;
  T p2NW;
  T u2NW;
  T v2NW;

  /*--- South-West state ---*/
  T alpha1SW;
  T rho1SW;
  T p1SW;
  T u1SW;
  T v1SW;
  T rho2SW;
  T p2SW;
  T u2SW;
  T v2SW;

  /*--- North-East ---*/
  T alpha1NE;
  T rho1NE;
  T p1NE;
  T u1NE;
  T v1NE;
  T rho2NE;
  T p2NE;
  T u2NE;
  T v2NE;

  /*--- South-East ---*/
  T alpha1SE;
  T rho1SE;
  T p1SE;
  T u1SE;
  T v1SE;
  T rho2SE;
  T p2SE;
  T u2SE;
  T v2SE;
};
