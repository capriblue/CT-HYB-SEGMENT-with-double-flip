#include "./evaluator.hpp"

void evaluate_nnw(const alps::accumulators::result_set &results,
                  const alps::params &parms,
                  alps::hdf5::archive &solver_output){

  if(!(parms["cthyb.MEASURE_nnw"].as<bool>())) return;

  std::size_t N_W=parms["cthyb.N_W"];
  std::size_t n_orbitals=parms["FLAVORS"];
  double beta=parms["BETA"];

  std::vector<std::vector<double> > nnw_re(n_orbitals*(n_orbitals+1)/2);
  int pos=0;
  for(std::size_t i=0;i<n_orbitals;++i){
    for(std::size_t j=0;j<=i;++j){
      std::stringstream nnw_re_name; nnw_re_name<<"nnw_re_"<<i<<"_"<<j;
      nnw_re[pos]=results[nnw_re_name.str()].mean<std::vector<double> >();
      solver_output<<alps::make_pvp(nnw_re_name.str(), nnw_re[pos++]);
    }
  }

  if(parms["cthyb.TEXT_OUTPUT"].as<bool>()){
    std::ofstream nnw_file("nnw.dat");
    nnw_file << "#w";
    for(std::size_t i=0;i<n_orbitals;++i)
      for(std::size_t j=0;j<=i;++j)
        nnw_file<<" nnw_"<<i<<j;
    nnw_file << std::endl;
    for(std::size_t m=0;m<N_W;++m){
      pos=0;
      double wm=2*m*M_PI/beta;
      nnw_file << wm;
      for(std::size_t i=0;i<n_orbitals;++i)
        for(std::size_t j=0;j<=i;++j)
          nnw_file << " " << nnw_re[pos++][m];
      nnw_file << std::endl;
    }
    nnw_file.close();
  }
}

