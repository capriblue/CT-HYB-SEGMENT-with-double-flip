#include "./evaluator.hpp"

void evaluate_nnt(const alps::accumulators::result_set &results,
                  const alps::params &parms,
                  alps::hdf5::archive &solver_output)
{

  if (!(parms["cthyb.MEASURE_nnt"].as<bool>()))
    return;

  std::size_t N_nn = parms["cthyb.N_nn"];
  std::size_t n_orbitals = parms["FLAVORS"];
  double beta = parms["BETA"];

  std::vector<std::vector<double>> nnt(n_orbitals * (n_orbitals + 1) / 2);
  int pos = 0;
  for (std::size_t i = 0; i < n_orbitals; ++i)
  {
    for (std::size_t j = 0; j <= i; ++j)
    {
      std::stringstream nnt_name;
      nnt_name << "nnt_" << i << "_" << j;
      nnt[pos] = results[nnt_name.str()].mean<std::vector<double>>();
      solver_output << alps::make_pvp(nnt_name.str(), nnt[pos++]);
    }
  }

  if (parms["cthyb.TEXT_OUTPUT"].as<bool>())
  {
    std::ofstream nnt_file("nnt.dat");
    nnt_file << "#tau";
    for (std::size_t i = 0; i < n_orbitals; ++i)
      for (std::size_t j = 0; j <= i; ++j)
        nnt_file << " nnt_" << i << j;
    nnt_file << std::endl;
    for (std::size_t n = 0; n <= N_nn; ++n)
    {
      pos = 0;
      double tau = n * beta / N_nn;
      nnt_file << tau;
      for (std::size_t i = 0; i < n_orbitals; ++i)
        for (std::size_t j = 0; j <= i; ++j)
          nnt_file << " " << nnt[pos++][n];
      nnt_file << std::endl;
    }
    nnt_file.close();
  }
}
