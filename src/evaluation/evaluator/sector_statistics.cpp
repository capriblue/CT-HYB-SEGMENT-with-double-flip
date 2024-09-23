#include "./evaluator.hpp"

void evaluate_sector_statistics(const alps::accumulators::result_set &results,
                                const alps::params &parms,
                                alps::hdf5::archive &solver_output){

  if(!(parms["cthyb.MEASURE_sector_statistics"].as<bool>())) return;

  std::size_t n_orbitals=parms["FLAVORS"];
  std::ofstream stat_file("sector_statistics.dat");
  stat_file << "#state |n_1={0,1} n_2={0,1} ...> n_i={0,1}: orbital i {empty,occupied}" << std::endl;
  stat_file << "#rel weight (in perc)" << std::endl;
  int n_states=pow(2,n_orbitals);
  std::vector<double> sector_statistics=results["sector_statistics"].mean<std::vector<double> >();
  for(int n=0;n<n_states;++n){
    std::stringstream state; state << "|";
    std::size_t i=0; int m=n;
    while(m>=0 && i<n_orbitals){
      int rem=m%2; ++i; m/=2;
      state << rem;
    }
    state << ">";
    stat_file << n << "\t" << sector_statistics[n]*100. << "\t" << state.str() << std::endl;
  }
  stat_file.close();
}
