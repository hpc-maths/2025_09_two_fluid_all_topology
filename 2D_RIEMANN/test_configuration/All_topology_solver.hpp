// Copyright 2021 SAMURAI TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
//
// Author: Giuseppe Orlando, 2025
//
#include <samurai/algorithm/update.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/bc.hpp>
#include <samurai/box.hpp>
#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>

#include <filesystem>
namespace fs = std::filesystem;

/*--- Add header with auxiliary structs ---*/
#include "containers.hpp"

/*--- Add header with the boundary conditions ---*/
#include "user_bc.hpp"

/*--- Include the headers with the numerical fluxes ---*/
#include "Rusanov_flux.hpp"
#include "non_conservative_flux.hpp"

// This is the class for the simulation of the all-topology model
//
template<std::size_t dim>
class All_Topology_Solver {
public:
  using Config  = samurai::MRConfig<dim, 2>;
  using Field   = samurai::VectorField<samurai::MRMesh<Config>,
                                       double,
                                       Utilities::EquationData<dim>::NVARS,
                                       false>;
  using Number  = samurai::Flux<Field>::Number;  /*--- Define the shortcut for the arithmetic type ---*/
  using Indices = samurai::Flux<Field>::Indices; /*--- Shortcut for the indices storage ---*/

  All_Topology_Solver() = default; /*--- Default constructor. This will do nothing
                                         and basically will never be used ---*/

  All_Topology_Solver(const xt::xtensor_fixed<double, xt::xshape<dim>>& min_corner,
                      const xt::xtensor_fixed<double, xt::xshape<dim>>& max_corner,
                      const Simulation_Parameters<Number>& sim_param,
                      const EOS_Parameters<Number>& eos_param,
                      const Riemann_Parameters<Number>& Riemann_param); /*--- Class constrcutor with the arguments related
                                                                              to the grid and to the physics. ---*/

  void run(const std::size_t nfiles = 10); /*--- Function which actually executes the temporal loop ---*/

  template<class... Variables>
  void save(const std::string& suffix,
            const Variables&... fields); /*--- Routine to save the results ---*/

private:
  /*--- Now we declare some relevant variables ---*/
  const samurai::Box<double, dim> box;

  samurai::MRMesh<Config> mesh; /*--- Variable to store the mesh ---*/

  using Field_Scalar = samurai::ScalarField<decltype(mesh), Number>;
  using Field_Vect   = samurai::VectorField<decltype(mesh), Number, dim, false>;

  const Number t0; /*--- Initial time of the simulation ---*/
  const Number Tf; /*--- Final time of the simulation ---*/

  Number cfl; /*--- Courant number of the simulation so as to compute the time step ---*/
  Number dt;  /*--- Time-step (in general modified according to CFL) ---*/

  const SG_EOS<Number> EOS_phase1; /*--- Equation of state of phase 1 ---*/
  const SG_EOS<Number> EOS_phase2; /*--- Equation of state of phase 2 ---*/

  samurai::RusanovFlux<Field> numerical_flux_cons; /*--- function to compute the numerical flux for the conservative part
                                                         (this is necessary to call 'make_flux') ---*/

  samurai::NonConservativeFlux<Field> numerical_flux_non_cons; /*--- function to compute the numerical flux for the non-conservative part
                                                                     (this is necessary to call 'make_flux') ---*/

  fs::path    path;     /*--- Auxiliary variable to store the output directory ---*/
  std::string filename; /*--- Auxiliary variable to store the name of output ---*/

  Field conserved_variables; /*--- The variable which stores the conserved variables,
                                   namely the variables for which we solve a PDE system ---*/

  /*-- Now we declare a bunch of fields which depend from the state, but it is useful
       to have it for the output ---*/
  Field_Scalar rho,
               p,
               rho1,
               p1,
               c1,
               rho2,
               p2,
               c2,
               alpha2,
               Y1,
               Y2,
               T1,
               T2,
               s1,
               s2,
               delta_pres,
               delta_temp;

  Field_Vect vel1,
             vel2,
             vel,
             delta_vel;

  /*--- Now, it's time to declare some member functions that we will employ ---*/
  void create_fields(); /*--- Auxiliary routine to initialize the fileds to the mesh ---*/

  void init_variables(const Riemann_Parameters<Number>& Riemann_param); /*--- Routine to initialize the variables
                                                                              (both conserved and auxiliary, this is problem dependent) ---*/

  void apply_bcs(); /*--- Auxiliary routine for the boundary conditions ---*/

  Number get_max_lambda() const; /*--- Compute the estimate of the maximum eigenvalue ---*/

  void update_auxiliary_fields(); /*--- Routine to update auxiliary fields for output and time step update ---*/
};

//////////////////////////////////////////////////////////////
/*---- START WITH THE IMPLEMENTATION OF THE CONSTRUCTOR ---*/
//////////////////////////////////////////////////////////////

// Implement class constructor
//
template<std::size_t dim>
All_Topology_Solver<dim>::All_Topology_Solver(const xt::xtensor_fixed<double, xt::xshape<dim>>& min_corner,
                                              const xt::xtensor_fixed<double, xt::xshape<dim>>& max_corner,
                                              const Simulation_Parameters<Number>& sim_param,
                                              const EOS_Parameters<Number>& eos_param,
                                              const Riemann_Parameters<Number>& Riemann_param):
  box(min_corner, max_corner),
  t0(sim_param.t0), Tf(sim_param.Tf),
  cfl(sim_param.Courant),
  EOS_phase1(eos_param.gamma_1, eos_param.pi_infty_1, eos_param.q_infty_1, eos_param.c_v_1),
  EOS_phase2(eos_param.gamma_2, eos_param.pi_infty_2, eos_param.q_infty_2, eos_param.c_v_2),
  numerical_flux_cons(EOS_phase1, EOS_phase2),
  numerical_flux_non_cons(EOS_phase1, EOS_phase2)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) {
      std::cout << "Initializing variables " << std::endl;
      std::cout << std::endl;
    }

    /*--- Attach the fields to the mesh ---*/
    create_fields();

    /*--- Initialize the fields ---*/
    mesh = {box, sim_param.min_level, sim_param.max_level, {{false, false}}};
    init_variables(Riemann_param);

    /*--- Apply boundary conditions ---*/
    apply_bcs();
  }

// Auxiliary routine to create the fields
//
template<std::size_t dim>
void All_Topology_Solver<dim>::create_fields() {
  /*--- Create conserved and auxiliary fields ---*/
  conserved_variables = samurai::make_vector_field<Number, Field::n_comp>("conserved", mesh);

  rho  = samurai::make_scalar_field<Number>("rho", mesh);
  p    = samurai::make_scalar_field<Number>("p", mesh);

  rho1 = samurai::make_scalar_field<Number>("rho1", mesh);
  p1   = samurai::make_scalar_field<Number>("p1", mesh);
  c1   = samurai::make_scalar_field<Number>("c1", mesh);

  rho2 = samurai::make_scalar_field<Number>("rho2", mesh);
  p2   = samurai::make_scalar_field<Number>("p2", mesh);
  c2   = samurai::make_scalar_field<Number>("c2", mesh);

  vel1 = samurai::make_vector_field<Number, dim>("vel1", mesh);
  vel2 = samurai::make_vector_field<Number, dim>("vel2", mesh);
  vel  = samurai::make_vector_field<Number, dim>("vel", mesh);

  alpha2 = samurai::make_scalar_field<Number>("alpha2", mesh);
  Y1     = samurai::make_scalar_field<Number>("Y1", mesh);
  Y2     = samurai::make_scalar_field<Number>("Y2", mesh);

  T1 = samurai::make_scalar_field<Number>("T1", mesh);
  T2 = samurai::make_scalar_field<Number>("T2", mesh);

  s1 = samurai::make_scalar_field<Number>("s1", mesh);
  s2 = samurai::make_scalar_field<Number>("s2", mesh);

  delta_pres = samurai::make_scalar_field<Number>("delta_pres", mesh);
  delta_temp = samurai::make_scalar_field<Number>("delta_temp", mesh);
  delta_vel  = samurai::make_vector_field<Number, dim>("delta_vel", mesh);
}

// Initialization of conserved and auxiliary variables
//
template<std::size_t dim>
void All_Topology_Solver<dim>::init_variables(const Riemann_Parameters<Number>& Riemann_param) {
  /*--- Reisze fields since now mesh has been created ---*/
  conserved_variables.resize();
  rho.resize();
  p.resize();
  rho1.resize();
  p1.resize();
  c1.resize();
  rho2.resize();
  p2.resize();
  c2.resize();
  vel1.resize();
  vel2.resize();
  vel.resize();
  alpha2.resize();
  Y1.resize();
  Y2.resize();
  T1.resize();
  T2.resize();
  s1.resize();
  s2.resize();
  delta_pres.resize();
  delta_temp.resize();
  delta_vel.resize();

  /*--- Initialize the fields with a loop over all cells ---*/
  samurai::for_each_cell(mesh,
                         [&](const auto& cell)
                            {
                              const auto center = cell.center();
                              const auto x      = static_cast<Number>(center[0]);
                              const auto y      = static_cast<Number>(center[1]);

                              // South-West state (primitive variables)
                              if(x <= Riemann_param.xd && y <= Riemann_param.yd) {
                                conserved_variables[cell][Indices::ALPHA1_INDEX] = Riemann_param.alpha1SW;

                                p1[cell]      = Riemann_param.p1SW;
                                vel1[cell][0] = Riemann_param.u1SW;
                                vel1[cell][1] = Riemann_param.v1SW;
                                rho1[cell]    = Riemann_param.rho1SW;

                                p2[cell]      = Riemann_param.p2SW;
                                vel2[cell][0] = Riemann_param.u2SW;
                                vel2[cell][1] = Riemann_param.v2SW;
                                rho2[cell]    = Riemann_param.rho2SW;
                              }
                              // North-West state (primitive variables)
                              else if(x <= Riemann_param.xd && y > Riemann_param.yd) {
                                conserved_variables[cell][Indices::ALPHA1_INDEX] = Riemann_param.alpha1NW;

                                p1[cell]      = Riemann_param.p1NW;
                                vel1[cell][0] = Riemann_param.u1NW;
                                vel1[cell][1] = Riemann_param.v1NW;
                                rho1[cell]    = Riemann_param.rho1NW;

                                p2[cell]      = Riemann_param.p2NW;
                                vel2[cell][0] = Riemann_param.u2NW;
                                vel2[cell][1] = Riemann_param.v2NW;
                                rho2[cell]    = Riemann_param.rho2NW;
                              }
                              // South-East state (primitive varaibles)
                              else if(x > Riemann_param.xd && y <= Riemann_param.yd) {
                                conserved_variables[cell][Indices::ALPHA1_INDEX] = Riemann_param.alpha1SE;

                                p1[cell]      = Riemann_param.p1SE;
                                vel1[cell][0] = Riemann_param.u1SE;
                                vel1[cell][1] = Riemann_param.v1SE;
                                rho1[cell]    = Riemann_param.rho1SE;

                                p2[cell]      = Riemann_param.p2SE;
                                vel2[cell][0] = Riemann_param.u2SE;
                                vel2[cell][1] = Riemann_param.v2SE;
                                rho2[cell]    = Riemann_param.rho2SE;
                              }
                              // North-East state (primitive variables)
                              else {
                                conserved_variables[cell][Indices::ALPHA1_INDEX] = Riemann_param.alpha1NE;

                                p1[cell]      = Riemann_param.p1NE;
                                vel1[cell][0] = Riemann_param.u1NE;
                                vel1[cell][1] = Riemann_param.v1NE;
                                rho1[cell]    = Riemann_param.rho1NE;

                                p2[cell]      = Riemann_param.p2NE;
                                vel2[cell][0] = Riemann_param.u2NE;
                                vel2[cell][1] = Riemann_param.v2NE;
                                rho2[cell]    = Riemann_param.rho2NE;
                              }

                              T1[cell] = EOS_phase1.T_value_RhoP(rho1[cell], p1[cell]);

                              conserved_variables[cell][Indices::ALPHA1_RHO1_INDEX] = conserved_variables[cell][Indices::ALPHA1_INDEX]*
                                                                                      rho1[cell];

                              const auto e1   = EOS_phase1.e_value_RhoP(rho1[cell], p1[cell]);
                              auto norm2_vel1 = static_cast<Number>(0.0);
                              for(std::size_t d = 0; d < dim; ++d) {
                                conserved_variables[cell][Indices::ALPHA1_RHO1_U1_INDEX + d] =
                                conserved_variables[cell][Indices::ALPHA1_RHO1_INDEX]*vel1[cell][d];
                                norm2_vel1 += vel1[cell][d]*vel1[cell][d];
                              }

                              conserved_variables[cell][Indices::ALPHA1_RHO1_E1_INDEX] =
                              conserved_variables[cell][Indices::ALPHA1_RHO1_INDEX]*
                              (e1 + static_cast<Number>(0.5)*norm2_vel1);

                              T2[cell] = EOS_phase2.T_value_RhoP(rho2[cell], p2[cell]);

                              alpha2[cell] = static_cast<Number>(1.0)
                                           - conserved_variables[cell][Indices::ALPHA1_INDEX];
                              conserved_variables[cell][Indices::ALPHA2_RHO2_INDEX] = alpha2[cell]*rho2[cell];

                              const auto e2   = EOS_phase2.e_value_RhoP(rho2[cell], p2[cell]);
                              auto norm2_vel2 = static_cast<Number>(0.0);
                              for(std::size_t d = 0; d < dim; ++d) {
                                conserved_variables[cell][Indices::ALPHA2_RHO2_U2_INDEX + d] =
                                conserved_variables[cell][Indices::ALPHA2_RHO2_INDEX]*vel2[cell][d];
                                norm2_vel2 += vel2[cell][d]*vel2[cell][d];
                              }
                              conserved_variables[cell][Indices::ALPHA2_RHO2_E2_INDEX] =
                              conserved_variables[cell][Indices::ALPHA2_RHO2_INDEX]*
                              (e2 + static_cast<Number>(0.5)*norm2_vel2);

                              c1[cell] = EOS_phase1.c_value_RhoP(rho1[cell], p1[cell]);

                              c2[cell] = EOS_phase2.c_value_RhoP(rho2[cell], p2[cell]);

                              rho[cell] = conserved_variables[cell][Indices::ALPHA1_RHO1_INDEX]
                                        + conserved_variables[cell][Indices::ALPHA2_RHO2_INDEX];

                              p[cell] = conserved_variables[cell][Indices::ALPHA1_INDEX]*p1[cell]
                                      + alpha2[cell]*p2[cell];

                              Y1[cell] = conserved_variables[cell][Indices::ALPHA1_RHO1_INDEX]/rho[cell];
                              Y2[cell] = conserved_variables[cell][Indices::ALPHA2_RHO2_INDEX]/rho[cell];

                              for(std::size_t d = 0; d < dim; ++d) {
                                vel[cell][d] = Y1[cell]*vel1[cell][d]
                                             + Y2[cell]*vel2[cell][d];
                              }

                              // Save deltas
                              delta_pres[cell] = p1[cell] - p2[cell];
                              delta_temp[cell] = T1[cell] - T2[cell];
                              for(std::size_t d = 0; d < dim; ++d) {
                                delta_vel[cell][d] = vel1[cell][d] - vel2[cell][d];
                              }

                              // Compute entropies
                              s1[cell] = EOS_phase1.s_value_Rhoe(rho1[cell], e1);
                              s2[cell] = EOS_phase2.s_value_Rhoe(rho2[cell], e2);
                            }
                        );
}

// Auxiliary routine to impose the boundary conditions
//
template<std::size_t dim>
void All_Topology_Solver<dim>::apply_bcs() {
  samurai::make_bc<Default>(conserved_variables, Reflecting(conserved_variables));
}

//////////////////////////////////////////////////////////////
/*---- FOCUS NOW ON THE AUXILIARY FUNCTIONS ---*/
/////////////////////////////////////////////////////////////

// Compute the estimate of the maximum eigenvalue for CFL condition
//
template<std::size_t dim>
typename All_Topology_Solver<dim>::Number All_Topology_Solver<dim>::get_max_lambda() const {
  auto local_res = static_cast<Number>(0.0);

  samurai::for_each_cell(mesh,
                         [&](const auto& cell)
                            {
                              for(std::size_t d = 0; d < dim; ++d) {
                                local_res = std::max(std::max(std::abs(vel1[cell][d]) + c1[cell],
                                                              std::abs(vel2[cell][d]) + c2[cell]),
                                                     local_res);
                              }
                            }
                        );

  double global_res;
  MPI_Allreduce(&local_res, &global_res, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  return global_res;
}

// Update auxiliary fields after solution of the system
//
template<std::size_t dim>
void All_Topology_Solver<dim>::update_auxiliary_fields() {
  /*--- Loop over all cells ---*/
  samurai::for_each_cell(mesh,
                         [&](const auto& cell)
                            {
                              // Pre-fetch varaibles used multiple times in order to exploit (possible) vectorization
                              // as well as to enhance readability
                              const auto alpha1_loc = conserved_variables[cell][Indices::ALPHA1_INDEX];
                              const auto m1_loc     = conserved_variables[cell][Indices::ALPHA1_RHO1_INDEX];
                              const auto m1E1_loc   = conserved_variables[cell][Indices::ALPHA1_RHO1_E1_INDEX];
                              const auto m2_loc     = conserved_variables[cell][Indices::ALPHA2_RHO2_INDEX];
                              const auto m2E2_loc   = conserved_variables[cell][Indices::ALPHA2_RHO2_E2_INDEX];

                              // Compute the fields
                              const auto rho_loc     = m1_loc + m2_loc;
                              const auto inv_rho_loc = static_cast<Number>(1.0)/rho_loc;
                              rho[cell]              = rho_loc;

                              const auto rho1_loc   = m1_loc/alpha1_loc; /*--- TODO: Add treatment for vanishing volume fraction ---*/
                              rho1[cell]            = rho1_loc;
                              const auto inv_m1_loc = static_cast<Number>(1.0)/m1_loc;
                                                      /*--- TODO: Add treatment for vanishing volume fraction ---*/
                              auto e1_loc           = m1E1_loc*inv_m1_loc; /*--- TODO: Add treatment for vanishing volume fraction ---*/
                              for(std::size_t d = 0; d < dim; ++d) {
                                vel1[cell][d] = conserved_variables[cell][Indices::ALPHA1_RHO1_U1_INDEX + d]*inv_m1_loc;
                                                /*--- TODO: Add treatment for vanishing volume fraction ---*/
                                e1_loc -= static_cast<Number>(0.5)*(vel1[cell][d]*vel1[cell][d]);
                              }
                              const auto p1_loc = EOS_phase1.pres_value_Rhoe(rho1_loc, e1_loc);
                              p1[cell]          = p1_loc;
                              c1[cell]          = EOS_phase1.c_value_RhoP(rho1_loc, p1_loc);
                              const auto T1_loc = EOS_phase1.T_value_RhoP(rho1_loc, p1_loc);
                              T1[cell]          = T1_loc;

                              const auto alpha2_loc = static_cast<Number>(1.0) - alpha1_loc;
                              const auto rho2_loc   = m2_loc/alpha2_loc; /*--- TODO: Add treatment for vanishing volume fraction ---*/
                              rho2[cell]            = rho2_loc;
                              const auto inv_m2_loc = static_cast<Number>(1.0)/m2_loc;
                                                      /*--- TODO: Add treatment for vanishing volume fraction ---*/
                              auto e2_loc           = m2E2_loc*inv_m2_loc; /*--- TODO: Add treatment for vanishing volume fraction ---*/
                              for(std::size_t d = 0; d < dim; ++d) {
                                vel2[cell][d] = conserved_variables[cell][Indices::ALPHA2_RHO2_U2_INDEX + d]*inv_m2_loc;
                                                /*--- TODO: Add treatment for vanishing volume fraction ---*/
                                e2_loc -= static_cast<Number>(0.5)*(vel2[cell][d]*vel2[cell][d]);
                              }
                              const auto p2_loc = EOS_phase2.pres_value_Rhoe(rho2_loc, e2_loc);
                              p2[cell]          = p2_loc;
                              c2[cell]          = EOS_phase2.c_value_RhoP(rho2_loc, p2_loc);
                              const auto T2_loc = EOS_phase2.T_value_RhoP(rho2_loc, p2_loc);
                              T2[cell]          = T2_loc;

                              // Compute remaining useful fields
                              alpha2[cell] = alpha2_loc;

                              const auto Y1_loc = m1_loc*inv_rho_loc;
                              Y1[cell]          = Y1_loc;
                              const auto Y2_loc = static_cast<Number>(1.0) - Y1_loc;
                              Y2[cell]          = Y2_loc;

                              for(std::size_t d = 0; d < dim; ++d) {
                                vel[cell][d] = Y1_loc*vel1[cell][d]
                                             + Y2_loc*vel2[cell][d];
                              }

                              p[cell] = alpha1_loc*p1_loc
                                      + alpha2_loc*p2_loc;

                              // Save deltas
                              delta_pres[cell] = p1_loc - p2_loc;
                              delta_temp[cell] = T1_loc - T2_loc;
                              for(std::size_t d = 0; d < dim; ++d) {
                                delta_vel[cell][d] = vel1[cell][d] - vel2[cell][d];
                              }

                              // Compute entropies
                              s1[cell] = EOS_phase1.s_value_Rhoe(rho1_loc, e1_loc);
                              s2[cell] = EOS_phase2.s_value_Rhoe(rho2_loc, e2_loc);
                            }
                        );
}

// Save desired fields and info
//
template<std::size_t dim>
template<class... Variables>
void All_Topology_Solver<dim>::save(const std::string& suffix,
                                    const Variables&... fields) {
  auto level_ = samurai::make_scalar_field<std::size_t>("level", mesh);

  if(!fs::exists(path)) {
    fs::create_directory(path);
  }

  samurai::for_each_cell(mesh,
                         [&](const auto& cell)
                            {
                              level_[cell] = cell.level;
                            }
                        );

  samurai::save(path, fmt::format("{}{}", filename, suffix), mesh, fields..., level_);
}

//////////////////////////////////////////////////////////////
/*---- IMPLEMENT THE FUNCTION THAT EFFECTIVELY SOLVES THE PROBLEM ---*/
/////////////////////////////////////////////////////////////

// Implement the function that effectively performs the temporal loop
//
template<std::size_t dim>
void All_Topology_Solver<dim>::run(const std::size_t nfiles) {
  /*--- Default output arguemnts ---*/
  path = fs::current_path();
  filename = "Rusanov_Flux_order1";

  const auto dt_save = Tf/static_cast<Number>(nfiles);

  /*--- Auxiliary variables to save updated fields ---*/
  auto conserved_variables_np1 = samurai::make_vector_field<Number, Field::n_comp>("conserved_np1", mesh);

  /*--- Create the flux variables ---*/
  auto Rusanov_flux         = numerical_flux_cons.make_flux();
  auto NonConservative_flux = numerical_flux_non_cons.make_flux();

  /*--- Save the initial condition ---*/
  const std::string suffix_init = (nfiles != 1) ? "_ite_" + Utilities::unsigned_to_string(0) : "";
  save(suffix_init, conserved_variables,
                    rho, p, vel,
                    vel1, rho1, p1, c1, T1, s1, Y1,
                    vel2, rho2, p2, c2, T2, s2, alpha2, Y2,
                    delta_pres, delta_temp, delta_vel);

  /*--- Set mesh size ---*/
  const auto dx = static_cast<Number>(mesh.cell_length(mesh.max_level()));
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  using mesh_id_t = typename decltype(mesh)::mesh_id_t;
  const auto n_elements_per_subdomain = mesh[mesh_id_t::cells].nb_cells();
  unsigned n_elements;
  MPI_Allreduce(&n_elements_per_subdomain, &n_elements, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  if(rank == 0) {
    std::cout << "Number of initial elements = " <<  n_elements << std::endl;
    std::cout << std::endl;
  }

  /*--- Start the loop ---*/
  std::size_t nsave = 0;
  std::size_t nt    = 0;
  auto t            = static_cast<Number>(t0);
  while(t != Tf) {
    // Compute time step
    samurai::update_ghost_mr(conserved_variables);
    auto Cons_Flux    = Rusanov_flux(conserved_variables);
    auto NonCons_Flux = NonConservative_flux(conserved_variables);
    dt = std::min(Tf - t, cfl*dx/get_max_lambda());
    t += dt;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) {
      std::cout << fmt::format("Iteration {}: t = {}, dt = {}", ++nt, t, dt) << std::endl;
    }

    // Apply the numerical scheme
    conserved_variables_np1.resize();
    conserved_variables_np1 = conserved_variables - dt*Cons_Flux - dt*NonCons_Flux;
    std::swap(conserved_variables.array(), conserved_variables_np1.array());

    // Save the results
    update_auxiliary_fields();
    if(t >= static_cast<Number>(nsave + 1)*dt_save || t == Tf) {
      const std::string suffix = (nfiles != 1) ? "_ite_" + Utilities::unsigned_to_string(++nsave) : "";
      save(suffix, conserved_variables,
                   rho, p, vel,
                   vel1, rho1, p1, c1, T1, s1, Y1,
                   vel2, rho2, p2, c2, T2, s2, alpha2, Y2,
                   delta_pres, delta_temp, delta_vel);
    }
  }
}
