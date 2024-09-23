#include "./evaluator.hpp"
void evaluate_basics(const alps::accumulators::result_set &results,
                     const alps::params &parms,
                     alps::hdf5::archive &solver_output)
{

  std::size_t n_orbitals = parms["FLAVORS"];
  double beta = parms["BETA"];

  if (parms["cthyb.TEXT_OUTPUT"])
  {
    std::ofstream sim_file("simulation.dat");
    sim_file << "simulation details:" << std::endl;
    sim_file << "average sign: " << results["Sign"].mean<double>() << std::endl;
    sim_file << "total (effective) number of sweeps, normalized by N_meas: " << results["Sign"].count() * (double)parms["cthyb.N_MEAS"] << std::endl;
    sim_file << "number thermalization sweeps: " << parms["cthyb.THERMALIZATION"] << std::endl;
    sim_file << "inverse temperature: " << beta << std::endl;
    sim_file << "perturbation order:" << std::endl;
    for (std::size_t i = 0; i < n_orbitals; ++i)
    { // replace Green function endpoints by corresponding densities
      std::stringstream order_name;
      order_name << "order_" << i;
      double order = results[order_name.str()].mean<double>();
      sim_file << "orbital " << i << ": " << order << std::endl;
    }
    {
      int tot_acc = 0, cur_prec = sim_file.precision();
      for (int i = 0; i < nacc.size(); i++)
        tot_acc += nacc[i];
      sim_file << std::endl
               << "|------------- Simulation details after " << nsweeps << " sweeps ------------|" << std::endl;
      sim_file << "  Total acceptance rate = " << std::setprecision(2) << std::fixed;
      sim_file << (((double)tot_acc) / nsweeps) * 100 << "%" << std::endl;
      sim_file << "  Individual acceptance rate for update " << std::endl;
      for (int i = 0; i < nacc.size(); i++)
      {
        sim_file << "     " << update_type[i] << " = ";
        sim_file << std::setprecision(2) << std::fixed << (((double)nacc[i]) / nsweeps) * 100 << "%";
        sim_file << " (proposal rate = ";
        sim_file << std::setprecision(2) << std::fixed << (((double)nprop[i]) / nsweeps) * 100 << "%)" << std::endl;
      }
      sim_file << "|-----------------------------------------------------------------|" << std::endl;
    }
    sim_file.close();
    std::ofstream obs_file("observables.dat"); // equal-time correlators
    for (std::size_t i = 0; i < n_orbitals; ++i)
    { // replace Green function endpoints by corresponding densities
      std::stringstream density_name;
      density_name << "density_" << i;
      double density = results[density_name.str()].mean<double>();
      obs_file << "n" << i << "=" << density << ";" << std::endl;
      if (parms["cthyb.MEASURE_nn"])
      {
        for (std::size_t j = 0; j < i; ++j)
        {
          std::stringstream nn_name;
          nn_name << "nn_" << i << "_" << j;
          double nn = results[nn_name.str()].mean<double>();
          obs_file << "n" << i << "n" << j << "=" << nn << ";" << std::endl;
        }
      }
    }
    obs_file.close();

    std::ofstream order_file("orders.dat");
    std::vector<std::vector<double>> order_histogram(n_orbitals);
    std::vector<std::vector<double>> order_histogram_err(n_orbitals);
    for (std::size_t j = 0; j < n_orbitals; ++j)
    {
      std::stringstream order_name;
      order_name << "order_histogram_" << j;
      order_histogram[j] = results[order_name.str()].mean<std::vector<double>>();
      order_histogram_err[j] = results[order_name.str()].error<std::vector<double>>();
    }
    for (std::size_t j = 0; j < n_orbitals; ++j)
    {
      for (std::size_t k = 0; k < order_histogram[0].size(); ++k)
      {
        order_file << k;
        order_file << " " << order_histogram[j][k] << " " << order_histogram_err[j][k] << std::endl;
      }
      order_file << std::endl;
    }
    order_file.close();
  } // text_output
}
