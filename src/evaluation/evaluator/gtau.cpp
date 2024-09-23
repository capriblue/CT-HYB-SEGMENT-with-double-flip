#include "./evaluator.hpp"

void evaluate_gtau(const alps::accumulators::result_set &results,
                   const alps::params &parms,
                   alps::hdf5::archive &solver_output)
{

  std::size_t N_t = parms["N"];
  double beta = parms["BETA"];
  std::size_t n_orbitals = parms["FLAVORS"];
  std::size_t n_sites = 1;

  // Imaginary time Green function
  itime_green_function_t G_tau(N_t + 1, n_sites, n_orbitals);
  for (std::size_t i = 0; i < n_orbitals; ++i)
  {
    std::stringstream g_name;
    g_name << "g_" << i;
    std::vector<double> G = results[g_name.str()].mean<std::vector<double>>();
    for (std::size_t t = 0; t < N_t + 1; ++t)
    {
      G_tau(t, 0, 0, i) = G[t];
    }
    G_tau(0, 0, 0, i) *= 2.;   // first and last bin
    G_tau(N_t, 0, 0, i) *= 2.; // have half the size
  }

  for (std::size_t i = 0; i < n_orbitals; ++i)
  { // replace Green function endpoints by corresponding densities
    std::stringstream density_name;
    density_name << "density_" << i;
    double density = results[density_name.str()].mean<double>();
    G_tau(0, 0, 0, i) = -1. * (1 - density);
    G_tau(N_t, 0, 0, i) = -1. * (density);
  }

  // store in hdf5
  G_tau.write_hdf5(solver_output, "/G_tau");

  if (parms.exists("cthyb.DMFT_FRAMEWORK") && parms["cthyb.DMFT_FRAMEWORK"] && parms.exists("solver.OUTFILE_H5GF"))
  {
    alps::gf::itime_sigma_gf_with_tail G_tau_h5gf = translate_Gt_to_h5gf(G_tau, parms);
    alps::hdf5::archive ar(parms["solver.OUTFILE_H5GF"], alps::hdf5::archive::WRITE);
    G_tau_h5gf.save(ar, "/G_tau");
    if (!(parms["cthyb.MEASURE_freq"]))
    {
      throw std::logic_error("we did not measure in frequency. For the ALPS framework you need to add frequency measurement cthyb.MEASURE_freq");
      std::size_t N_w = parms["NMATSUBARA"];
      alps::gf::omega_sigma_gf_with_tail G_omega_h5gf(alps::gf::omega_sigma_gf(alps::gf::matsubara_positive_mesh(beta, N_w), alps::gf::index_mesh(n_orbitals)));
      // alps::gf::fourier_time_to_frequency(G_tau_h5gf, G_omega_h5gf);
      G_omega_h5gf.save(ar, "/G_omega");
    }
  }

  if (parms["cthyb.TEXT_OUTPUT"])
  {
    std::ofstream G_file("Gt.dat");
    for (std::size_t t = 0; t <= N_t; ++t)
    {
      G_file << beta * t / N_t;
      for (std::size_t j = 0; j < n_orbitals; ++j)
      {
        G_file << " " << G_tau(t, 0, 0, j);
      }
      G_file << std::endl;
    }
    G_file.close();
  }
}
