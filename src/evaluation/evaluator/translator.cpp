#include "./evaluator.hpp"

alps::gf::itime_sigma_gf_with_tail translate_Gt_to_h5gf(itime_green_function_t Gtau, const alps::params &parms)
{
  double beta = parms["BETA"];
  int n_tau = parms["N"];
  int n_orbitals = parms["FLAVORS"];

  alps::gf::itime_sigma_gf_with_tail Gtau_h5gf(alps::gf::itime_sigma_gf(alps::gf::itime_mesh(beta, n_tau + 1), alps::gf::index_mesh(n_orbitals)));

  typedef alps::gf::one_index_gf<double, alps::gf::index_mesh> density_matrix_type;
  density_matrix_type tail = density_matrix_type(alps::gf::index_mesh(n_orbitals));
  tail.initialize();
  for (alps::gf::index s(0); s < Gtau_h5gf.mesh2().extent(); s++)
  {
    for (alps::gf::itime_index i(0); i < Gtau_h5gf.mesh1().extent(); i++)
    {
      Gtau_h5gf(i, s) = Gtau(i(), 0, 0, s());
    }
    tail(s) = 1.;
  }
  Gtau_h5gf.set_tail(1, tail);
  return Gtau_h5gf;
}

alps::gf::omega_sigma_gf_with_tail translate_Gw_to_h5gf(matsubara_green_function_t Gw, const alps::params &parms)
{
  double beta = parms["BETA"];
  int n_matsubara = parms["N"];
  int n_orbitals = parms["FLAVORS"];

  alps::gf::omega_sigma_gf_with_tail Gw_h5gf(alps::gf::omega_sigma_gf(alps::gf::matsubara_positive_mesh(beta, n_matsubara), alps::gf::index_mesh(n_orbitals)));
  typedef alps::gf::one_index_gf<double, alps::gf::index_mesh> density_matrix_type;
  density_matrix_type tail = density_matrix_type(alps::gf::index_mesh(n_orbitals));
  tail.initialize();
  for (alps::gf::index s(0); s < Gw_h5gf.mesh2().extent(); s++)
  {
    for (alps::gf::matsubara_index w(0); w < Gw_h5gf.mesh1().extent(); w++)
    {
      Gw_h5gf(w, s) = Gw(w(), 0, 0, s());
    }
    tail(s) = 1.;
  }
  Gw_h5gf.set_tail(1, tail);
  return Gw_h5gf;
}
