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

  auto& app = samurai::initialize("Non-conservative Rusanov scheme for the 2D all-topology model", argc, argv);

  std::ifstream ifs("input.json"); // Read a JSON file
  json input = json::parse(ifs);

  /*--- Set and declare simulation parameters ---*/
  const std::size_t dim = 2; // Space dimension
  using Number = All_Topology_Solver<dim>::Number;
  Simulation_Parameters<Number> sim_param;

  // Physical parameters
  sim_param.xL = input.value("xL", 0.0);
  sim_param.xR = input.value("xR", 1.0);
  sim_param.yL = input.value("yL", 0.0);
  sim_param.yR = input.value("yR", 1.0);

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
  app.add_option("--yL", sim_param.yL, "y Left-end of the domain")->capture_default_str()->group("Physical parameters");
  app.add_option("--yR", sim_param.yR, "y Right-end of the domain")->capture_default_str()->group("Physical parameters");

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

  Riemann_param.xd = input.value("xd", 0.0);
  Riemann_param.yd = input.value("yd", 0.0);

  Riemann_param.alpha1NW = input.value("alpha1NW", 0.4);
  Riemann_param.rho1NW   = input.value("rho1NW", 1.0);
  Riemann_param.p1NW     = input.value("p1NW", 1.0);
  Riemann_param.u1NW     = input.value("u1NW", 0.0);
  Riemann_param.v1NW     = input.value("v1NW", 0.0);
  Riemann_param.rho2NW   = input.value("rho2NW", 0.5);
  Riemann_param.p2NW     = input.value("p2NW", 1.0);
  Riemann_param.u2NW     = input.value("u2NW", 0.0);
  Riemann_param.v2NW     = input.value("v2NW", 0.0);

  Riemann_param.alpha1NE = input.value("alpha1NE", 0.8);
  Riemann_param.rho1NE   = input.value("rho1NE", 2.0);
  Riemann_param.p1NE     = input.value("p1NE", 2.0);
  Riemann_param.u1NE     = input.value("u1NE", 0.0);
  Riemann_param.v1NE     = input.value("v1NE", 0.0);
  Riemann_param.rho2NE   = input.value("rho2NE", 1.5);
  Riemann_param.p2NE     = input.value("p2NE", 2.0);
  Riemann_param.u2NE     = input.value("u2NE", 0.0);
  Riemann_param.v2NE     = input.value("v2NE", 0.0);

  Riemann_param.alpha1SW = input.value("alpha1SW", 0.8);
  Riemann_param.rho1SW   = input.value("rho1SW", 2.0);
  Riemann_param.p1SW     = input.value("p1SW", 2.0);
  Riemann_param.u1SW     = input.value("u1SW", 0.0);
  Riemann_param.v1SW     = input.value("v1SW", 0.0);
  Riemann_param.rho2SW   = input.value("rho2SW", 1.5);
  Riemann_param.p2SW     = input.value("p2SW", 2.0);
  Riemann_param.u2SW     = input.value("u2SW", 0.0);
  Riemann_param.v2SW     = input.value("v2SW", 0.0);

  Riemann_param.alpha1SE = input.value("alpha1SE", 0.4);
  Riemann_param.rho1SE   = input.value("rho1SE", 1.0);
  Riemann_param.p1SE     = input.value("p1SE", 1.0);
  Riemann_param.u1SE     = input.value("u1SE", 0.0);
  Riemann_param.v1SE     = input.value("v1SE", 0.0);
  Riemann_param.rho2SE   = input.value("rho2SE", 0.5);
  Riemann_param.p2SE     = input.value("p2SE", 1.0);
  Riemann_param.u2SE     = input.value("u2SE", 0.0);
  Riemann_param.v2SE     = input.value("v2SE", 0.0);

  app.add_option("--xd", Riemann_param.xd, "Initial discontinuity location along x")->capture_default_str()->group("Initial conditions");
  app.add_option("--yd", Riemann_param.yd, "Initial discontinuity location along y")->capture_default_str()->group("Initial conditions");

  app.add_option("--alpha1NW", Riemann_param.alpha1NW, "Initial volume fraction at top-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--rho1NW", Riemann_param.rho1NW, "Initial density phase 1 at top-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--p1NW", Riemann_param.p1NW, "Initial pressure phase 1 at top-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--u1NW", Riemann_param.u1NW, "Initial horizontal velocity phase 1 at top-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--v1NW", Riemann_param.v1NW, "Initial vertical velocity phase 1 at top-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--rho2NW", Riemann_param.rho2NW, "Initial density phase 2 at top-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--p2NW", Riemann_param.p2NW, "Initial pressure phase 2 at top-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--u2NW", Riemann_param.u2NW, "Initial horizontal velocity phase 2 at top-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--v2NW", Riemann_param.v2NW, "Initial vertical velocity phase 2 at top-left")->capture_default_str()->group("Initial conditions");

  app.add_option("--alpha1NE", Riemann_param.alpha1NE, "Initial volume fraction at top-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--rho1NE", Riemann_param.rho1NE, "Initial density phase 1 at top-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--p1NE", Riemann_param.p1NE, "Initial pressure phase 1 at top-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--u1NE", Riemann_param.u1NE, "Initial horizontal velocity phase 1 at top-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--v1NE", Riemann_param.v1NE, "Initial vertical velocity phase 1 at top-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--rho2NE", Riemann_param.rho2NE, "Initial density phase 2 at top-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--p2NE", Riemann_param.p2NE, "Initial pressure phase 2 at top-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--u2NE", Riemann_param.u2NE, "Initial horizontal velocity phase 2 at top-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--v2NE", Riemann_param.v2NE, "Initial vertical velocity phase 2 at top-right")->capture_default_str()->group("Initial conditions");

  app.add_option("--alpha1SW", Riemann_param.alpha1SW, "Initial volume fraction at bottom-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--rho1SW", Riemann_param.rho1SW, "Initial density phase 1 at bottom-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--p1SW", Riemann_param.p1SW, "Initial pressure phase 1 at bottom-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--u1SW", Riemann_param.u1SW, "Initial horizontal velocity phase 1 at bottom-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--v1SW", Riemann_param.v1SW, "Initial vertical velocity phase 1 at bottom-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--rho2SW", Riemann_param.rho2SW, "Initial density phase 2 at bottom-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--p2SW", Riemann_param.p2SW, "Initial pressure phase 2 at bottom-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--u2SW", Riemann_param.u2SW, "Initial horizontal velocity phase 2 at bottom-left")->capture_default_str()->group("Initial conditions");
  app.add_option("--v2SW", Riemann_param.v2SW, "Initial vertical velocity phase 2 at bottom-left")->capture_default_str()->group("Initial conditions");

  app.add_option("--alpha1SE", Riemann_param.alpha1SE, "Initial volume fraction at bottom-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--rho1SE", Riemann_param.rho1SE, "Initial density phase 1 at bottom-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--p1SE", Riemann_param.p1SE, "Initial pressure phase 1 at bottom-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--u1SE", Riemann_param.u1SE, "Initial horizontal velocity phase 1 at bottom-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--v1SE", Riemann_param.v1SE, "Initial vertical velocity phase 1 at bottom-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--rho2SE", Riemann_param.rho2SE, "Initial density phase 2 at bottom-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--p2SE", Riemann_param.p2SE, "Initial pressure phase 2 at bottom-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--u2SE", Riemann_param.u2SE, "Initial horizontal velocity phase 2 at bottom-right")->capture_default_str()->group("Initial conditions");
  app.add_option("--v2SE", Riemann_param.v2SE, "Initial vertical velocity phase 2 at bottom-right")->capture_default_str()->group("Initial conditions");

  /*--- Create the instance of the class to perform the simulation ---*/
  CLI11_PARSE(app, argc, argv);
  xt::xtensor_fixed<double, xt::xshape<dim>> min_corner = {sim_param.xL, sim_param.yL};
  xt::xtensor_fixed<double, xt::xshape<dim>> max_corner = {sim_param.xR, sim_param.yR};
  auto All_Topology_Solver_Sim = All_Topology_Solver(min_corner, max_corner, sim_param, eos_param, Riemann_param);

  All_Topology_Solver_Sim.run(sim_param.nfiles);

  samurai::finalize();

  return 0;
}
