// Copyright 2021 SAMURAI TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
//
// Author: Giuseppe Orlando, 2025
//
#include <CLI/CLI.hpp>

#include <nlohmann/json.hpp>

#include "test_configuration/All_topology_solver.hpp"

// Main function to run the program
//
int main(int argc, char* argv[]) {
  using json = nlohmann::json;

  auto& app = samurai::initialize("Rusanov flux with non-conservative terms for the all-topology model", argc, argv);

  std::ifstream ifs("input.json"); // Read a JSON file
  json input = json::parse(ifs);

  /*--- Set and declare simulation parameters ---*/
  const std::size_t dim = 1; // Spatial dimension
  using Number = All_Topology_Solver<dim>::Number;
  Simulation_Parameters<Number> sim_param;

  // Physical parameters
  sim_param.xL = input.value("xL", 0.0);
  sim_param.xR = input.value("xR", 1.0);

  sim_param.t0 = input.value("t0", 0.0);
  sim_param.Tf = input.value("Tf", 0.007);

  // Numerical parameters
  sim_param.Courant = input.value("cfl", 0.2);
  sim_param.dt      = input.value("dt", 1e-8);

  // Mesh parameters
  sim_param.min_level = input.value("min-level", static_cast<std::size_t>(10));
  sim_param.max_level = input.value("max-level", static_cast<std::size_t>(10));

  // Output parameters
  sim_param.nfiles = input.value("nfiles", static_cast<std::size_t>(10));

  /*--- Allow for parsing from command line ---*/
  // Physical parameters
  app.add_option("--xL", sim_param.xL, "x Left-end of the domain")->capture_default_str()->group("Physical parameters");
  app.add_option("--xR", sim_param.xR, "x Right-end of the domain")->capture_default_str()->group("Physical parameters");

  app.add_option("--t0", sim_param.t0, "Initial time")->capture_default_str()->group("Physical parameters");
  app.add_option("--Tf", sim_param.Tf, "Final time")->capture_default_str()->group("Physical parameters");

  // Numerical parameters
  app.add_option("--cfl", sim_param.Courant, "The Courant number")->capture_default_str()->group("Simulation parameters");
  app.add_option("--dt", sim_param.dt, "The time step")->capture_default_str()->group("Simulation parameters");

  // Mesh parameters
  app.add_option("--min-level", sim_param.min_level, "Minimum level of the AMR")->capture_default_str()->group("Mesh parameter");
  app.add_option("--max-level", sim_param.max_level, "Maximum level of the AMR")->capture_default_str()->group("Mesh parameter");

  // Output parameters
  app.add_option("--nfiles", sim_param.nfiles, "Number of output files")->capture_default_str()->group("Ouput parameters");

  /*--- Set and declare simulation parameters related to EOS ---*/
  EOS_Parameters<Number> eos_param;

  eos_param.gamma_1    = input.value("gamma_1", 3.0);
  eos_param.pi_infty_1 = input.value("pi_infty_1", 1e2);
  eos_param.q_infty_1  = input.value("q_infty_1", 0.0);
  eos_param.c_v_1      = input.value("c_v_1", 1.040e3);

  eos_param.gamma_2    = input.value("gamma_2", 1.4);
  eos_param.pi_infty_2 = input.value("pi_infty_2", 0.0);
  eos_param.q_infty_2  = input.value("q_infty_2", 0.0);
  eos_param.c_v_2      = input.value("c_v_2", 1.040e3);

  /*--- Allow for parsing from command line ---*/
  app.add_option("--gammma_1", eos_param.gamma_1, "gamma_1")->capture_default_str()->group("EOS parameters");
  app.add_option("--pi_infty_1", eos_param.pi_infty_1, "pi_infty_1")->capture_default_str()->group("EOS parameters");
  app.add_option("--q_infty_1", eos_param.q_infty_1, "q_infty_1")->capture_default_str()->group("EOS parameters");
  app.add_option("--c_v_1", eos_param.c_v_1, "c_v_1")->capture_default_str()->group("EOS parameters");

  app.add_option("--gammma_2", eos_param.gamma_2, "gamma_2")->capture_default_str()->group("EOS parameters");
  app.add_option("--pi_infty_2", eos_param.pi_infty_2, "pi_infty_2")->capture_default_str()->group("EOS parameters");
  app.add_option("--q_infty_2", eos_param.q_infty_2, "q_infty_2")->capture_default_str()->group("EOS parameters");
  app.add_option("--c_v_2", eos_param.c_v_2, "c_v_2")->capture_default_str()->group("EOS parameters");

  /*--- Set and declare simulation parameters related to initial condition ---*/
  Riemann_Parameters<Number> Riemann_param;

  Riemann_param.xd      = input.value("xd", 0.8);

  Riemann_param.alpha1L = input.value("alpha1L", 0.8);
  Riemann_param.rho1L   = input.value("rho1L", 1.0);
  Riemann_param.p1L     = input.value("p1L", 1e3);
  Riemann_param.T1L     = input.value("T1L", 363.0);
  Riemann_param.u1L     = input.value("u1L", -19.59716);

  Riemann_param.rho2L   = input.value("rho2L", 1.0);
  Riemann_param.p2L     = input.value("p2L", 1e3);
  Riemann_param.T2L     = input.value("T2L", 1000.0);
  Riemann_param.u2L     = input.value("u2L", -19.59741);

  Riemann_param.alpha1R = input.value("alpha1R", 0.3);
  Riemann_param.rho1R   = input.value("rho1R", 1.0);
  Riemann_param.p1R     = input.value("p1R", 1e-1);
  Riemann_param.T1R     = input.value("T1R", 363.0);
  Riemann_param.u1R     = input.value("u1R", -19.59741);

  Riemann_param.rho2R   = input.value("rho2R", 1.0);
  Riemann_param.p2R     = input.value("p2R", 1e-1);
  Riemann_param.T2R     = input.value("T2R", 1000.0);
  Riemann_param.u2R     = input.value("u2R", -19.59741);

  app.add_option("--xd", Riemann_param.xd, "Initial discontinuity location")->capture_default_str()->group("Initial conditions");

  app.add_option("--alpha1L", Riemann_param.alpha1L, "Initial volume fraction at left")->capture_default_str()->group("Initial conditions");
  app.add_option("--rho1L", Riemann_param.rho1L, "Initial density phase 1 at left")->capture_default_str()->group("Initial conditions");
  app.add_option("--p1L", Riemann_param.p1L, "Initial pressure phase 1 at left")->capture_default_str()->group("Initial conditions");
  app.add_option("--T1L", Riemann_param.T1L, "Initial temperature phase 1 at left")->capture_default_str()->group("Initial conditions");
  app.add_option("--u1L", Riemann_param.u1L, "Initial horizontal velocity phase 1 at left")->capture_default_str()->group("Initial conditions");
  app.add_option("--v1L", Riemann_param.v1L, "Initial vertical velocity phase 1 at left")->capture_default_str()->group("Initial conditions");

  app.add_option("--rho2L", Riemann_param.rho2L, "Initial density phase 2 at left")->capture_default_str()->group("Initial conditions");
  app.add_option("--p2L", Riemann_param.p2L, "Initial pressure phase 2 at left")->capture_default_str()->group("Initial conditions");
  app.add_option("--T2L", Riemann_param.T2L, "Initial temperature phase 2 at left")->capture_default_str()->group("Initial conditions");
  app.add_option("--u2L", Riemann_param.u2L, "Initial horizontal velocity phase 2 at left")->capture_default_str()->group("Initial conditions");
  app.add_option("--v2L", Riemann_param.v2L, "Initial vertical velocity phase 2 at left")->capture_default_str()->group("Initial conditions");

  app.add_option("--alpha1R", Riemann_param.alpha1R, "Initial volume fraction at right")->capture_default_str()->group("Initial conditions");
  app.add_option("--rho1R", Riemann_param.rho1R, "Initial density phase 1 at right")->capture_default_str()->group("Initial conditions");
  app.add_option("--p1R", Riemann_param.p1R, "Initial pressure phase 1 at right")->capture_default_str()->group("Initial conditions");
  app.add_option("--T1R", Riemann_param.T1R, "Initial temperature phase 1 at right")->capture_default_str()->group("Initial conditions");
  app.add_option("--u1R", Riemann_param.u1R, "Initial horizotnal velocity phase 1 at right")->capture_default_str()->group("Initial conditions");
  app.add_option("--v1R", Riemann_param.v1R, "Initial vertical velocity phase 1 at right")->capture_default_str()->group("Initial conditions");

  app.add_option("--rho2R", Riemann_param.rho2R, "Initial density phase 2 at right")->capture_default_str()->group("Initial conditions");
  app.add_option("--p2R", Riemann_param.p2R, "Initial pressure phase 2 at right")->capture_default_str()->group("Initial conditions");
  app.add_option("--T2R", Riemann_param.T2R, "Initial temperature phase 2 at right")->capture_default_str()->group("Initial conditions");
  app.add_option("--u2R", Riemann_param.u2R, "Initial horizontal velocity phase 2 at right")->capture_default_str()->group("Initial conditions");
  app.add_option("--v2R", Riemann_param.v2R, "Initial vertical velocity phase 2 at right")->capture_default_str()->group("Initial conditions");

  /*--- Create the instance of the class to perform the simulation ---*/
  CLI11_PARSE(app, argc, argv);
  xt::xtensor_fixed<double, xt::xshape<dim>> min_corner = {sim_param.xL};
  xt::xtensor_fixed<double, xt::xshape<dim>> max_corner = {sim_param.xR};
  auto All_Topology_Solver_Sim = All_Topology_Solver(min_corner, max_corner,
                                                     sim_param, eos_param, Riemann_param);

  All_Topology_Solver_Sim.run(sim_param.nfiles);

  samurai::finalize();

  return 0;
}
