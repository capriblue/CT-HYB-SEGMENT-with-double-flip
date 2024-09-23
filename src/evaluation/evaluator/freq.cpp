#include "./evaluator.hpp"

void evaluate_freq(const alps::accumulators::result_set &results,
                   const alps::params &parms,
                   alps::hdf5::archive &solver_output)
{

  if (!(parms["cthyb.MEASURE_freq"]))
    return;
  // evaluate Matsubara Green's function and self-energy
  double beta = parms["BETA"];
  std::size_t N_w = parms["NMATSUBARA"];
  std::size_t n_orbitals = parms["FLAVORS"];
  std::size_t n_sites = 1;

  matsubara_green_function_t G_omega(N_w, n_sites, n_orbitals);
  matsubara_green_function_t F_omega(N_w, n_sites, n_orbitals);
  matsubara_green_function_t S_omega(N_w, n_sites, n_orbitals);
  for (std::size_t i = 0; i < n_orbitals; ++i)
  {
    std::stringstream gw_re_name;
    gw_re_name << "gw_re_" << i;
    std::stringstream gw_im_name;
    gw_im_name << "gw_im_" << i;
    std::stringstream fw_re_name;
    fw_re_name << "fw_re_" << i;
    std::stringstream fw_im_name;
    fw_im_name << "fw_im_" << i;
    std::vector<double> Gw_re = results[gw_re_name.str()].mean<std::vector<double>>();
    std::vector<double> Gw_im = results[gw_im_name.str()].mean<std::vector<double>>();
    std::vector<double> Fw_re = results[fw_re_name.str()].mean<std::vector<double>>();
    std::vector<double> Fw_im = results[fw_im_name.str()].mean<std::vector<double>>();
    for (std::size_t w = 0; w < N_w; ++w)
    {
      std::complex<double> G(Gw_re[w], Gw_im[w]);
      std::complex<double> F(Fw_re[w], Fw_im[w]);
      G_omega(w, 0, 0, i) = G;
      F_omega(w, 0, 0, i) = F;
      S_omega(w, 0, 0, i) = F / G;
    }
  }

  // store in hdf5
  G_omega.write_hdf5(solver_output, "/G_omega");
  F_omega.write_hdf5(solver_output, "/F_omega");
  S_omega.write_hdf5(solver_output, "/S_omega");

  if (parms.exists("cthyb.DMFT_FRAMEWORK") && parms["cthyb.DMFT_FRAMEWORK"] && parms.exists("solver.OUTFILE_H5GF"))
  {
    int n_matsubara = parms["NMATSUBARA"];
    alps::hdf5::archive ar(parms["solver.OUTFILE_H5GF"], alps::hdf5::archive::WRITE);
    alps::hdf5::archive ar_in(parms["solver.INFILE_H5GF"], alps::hdf5::archive::READ);
    double shift = parms["U"].as<double>() / 2;
    alps::gf::omega_sigma_gf_with_tail G0_omega(alps::gf::omega_sigma_gf(alps::gf::matsubara_positive_mesh(beta, n_matsubara), alps::gf::index_mesh(n_orbitals)));
    G0_omega.load(ar_in, "/G0");
    for (alps::gf::matsubara_index i(0); i < G0_omega.mesh1().extent(); ++i)
    {
      for (alps::gf::index s(0); s < G0_omega.mesh2().extent(); ++s)
      {
        G_omega(i(), 0, 0, s()) = 1.0 / (1.0 / G0_omega(i, s) + shift - S_omega(i(), 0, 0, s()));
      }
    }
    translate_Gw_to_h5gf(G_omega, parms).save(ar, "/G_omega");
  }

  std::ofstream Gw_file("Gw.dat");
  for (std::size_t n = 0; n < N_w; ++n)
  {
    Gw_file << (2. * n + 1) * M_PI / beta;
    for (std::size_t j = 0; j < n_orbitals; ++j)
    {
      Gw_file << " " << G_omega(n, 0, 0, j).real() << " " << G_omega(n, 0, 0, j).imag();
    }
    Gw_file << std::endl;
  }
  Gw_file.close();

  std::ofstream Fw_file("Fw.dat");
  for (std::size_t n = 0; n < N_w; ++n)
  {
    Fw_file << (2. * n + 1) * M_PI / beta;
    for (std::size_t j = 0; j < n_orbitals; ++j)
    {
      Fw_file << " " << F_omega(n, 0, 0, j).real() << " " << F_omega(n, 0, 0, j).imag();
    }
    Fw_file << std::endl;
  }
  Fw_file.close();

  std::ofstream Sw_file("Sw.dat");
  for (std::size_t n = 0; n < N_w; ++n)
  {
    Sw_file << (2. * n + 1) * M_PI / beta;
    for (std::size_t j = 0; j < n_orbitals; ++j)
    {
      Sw_file << " " << S_omega(n, 0, 0, j).real() << " " << S_omega(n, 0, 0, j).imag();
    }
    Sw_file << std::endl;
  }
  Sw_file.close();
}
