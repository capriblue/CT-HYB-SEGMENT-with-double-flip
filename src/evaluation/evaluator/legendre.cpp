#include "./evaluator.hpp"

void evaluate_legendre(const alps::accumulators::result_set &results,
                       const alps::params &parms,
                       alps::hdf5::archive &solver_output){
  if(!(parms["cthyb.MEASURE_legendre"].as<bool>())) return;
  std::cout<<"evaluating legendre polynomial results"<<std::endl;
  double beta = parms["BETA"];
  std::size_t N_l=parms["cthyb.N_LEGENDRE"];
  std::size_t N_w=parms["NMATSUBARA"];
  std::size_t N_t = parms["N"];
  std::size_t n_orbitals = parms["FLAVORS"];
  std::size_t n_sites = 1;

  //Legendre Green function (evaluated in Matsubara)
  matsubara_green_function_t G_l_omega(N_w, n_sites, n_orbitals);
  matsubara_green_function_t F_l_omega(N_w, n_sites, n_orbitals);
  matsubara_green_function_t S_l_omega(N_w, n_sites, n_orbitals);
  //Legendre Green function (evaluated in imaginary time)
  itime_green_function_t G_l_tau(N_t+1, n_sites, n_orbitals);
  itime_green_function_t F_l_tau(N_t+1, n_sites, n_orbitals);
  std::vector<std::vector<std::vector<double> > >gc_conv(2, std::vector<std::vector<double> >(n_orbitals, std::vector<double>(N_l,0.)));
  std::vector<std::vector<std::vector<double> > >fc_conv(2, std::vector<std::vector<double> >(n_orbitals, std::vector<double>(N_l,0.)));
  for(std::size_t i=0;i<n_orbitals;++i){
    std::stringstream gl_name; gl_name<<"gl_"<<i;
    std::stringstream fl_name; fl_name<<"fl_"<<i;
    std::vector<double> Gl=results[gl_name.str()].mean<std::vector<double> >();
    std::vector<double> Fl=results[fl_name.str()].mean<std::vector<double> >();
    for(std::size_t wn=0; wn<N_w; ++wn){
      G_l_omega(wn,0,0,i)=0.;
      F_l_omega(wn,0,0,i)=0.;
      for(std::size_t l=0; l<N_l; ++l){
        G_l_omega(wn,0,0,i)+=t(wn,l)*sqrt(2.*l+1)*Gl[l]; //sqrt(2l+1) has been omitted in the measurement
        F_l_omega(wn,0,0,i)+=t(wn,l)*sqrt(2.*l+1)*Fl[l];
      }
      S_l_omega(wn,0,0,i)=F_l_omega(wn,0,0,i)/G_l_omega(wn,0,0,i);
    }
    //Imaginary time Green function from Legendre
    for(std::size_t t=0;t<N_t+1;++t){
      double tau=t*beta/N_t;
      G_l_tau(t,0,0,i)=0.;
      F_l_tau(t,0,0,i)=0.;
      double x=2.0*tau/beta-1.0;
      double pl_2=1; double pl_1=x; double legendre_p;
      for(std::size_t l=0;l<N_l;++l){
        if(l==0) legendre_p=1;
        else if(l==1) legendre_p=x;
        else{
          legendre_p=((2*l-1)*x*pl_1-(l-1)*pl_2)/static_cast<double>(l);//l
          pl_2=pl_1; //l-2
          pl_1=legendre_p; //l-1
        }
        G_l_tau(t,0,0,i)+=(2.*l+1)*Gl[l]*legendre_p/beta;
        F_l_tau(t,0,0,i)+=(2.*l+1)*Fl[l]*legendre_p/beta;
        if(t==N_t/2){
          gc_conv[0][i][l]=G_l_tau(t,0,0,i);
          fc_conv[0][i][l]=F_l_tau(t,0,0,i);
        }
        if(t==N_t){
          gc_conv[1][i][l]=G_l_tau(t,0,0,i);
          fc_conv[1][i][l]=F_l_tau(t,0,0,i);
        }
      }
    }
  }//i

  //store in hdf5
  G_l_omega.write_hdf5(solver_output, "/G_l_omega");
  F_l_omega.write_hdf5(solver_output, "/F_l_omega");
  S_l_omega.write_hdf5(solver_output, "/S_l_omega");
  G_l_tau.write_hdf5(solver_output, "/G_l_tau");
  F_l_tau.write_hdf5(solver_output, "/F_l_tau");

  if(parms["cthyb.TEXT_OUTPUT"]){
    std::ofstream gc_str("Gl_conv.dat");
    std::ofstream fc_str("Fl_conv.dat");
    gc_str << "#lc";
    fc_str << "#lc";
    for(std::size_t i=0;i<n_orbitals;++i){
      gc_str << " Gt_" << i << "(tau=beta/2,lc) Gt_" << i << "(tau=beta)";
      fc_str << " Ft_" << i << "(tau=beta/2,lc) Ft_" << i << "(tau=beta)";
    }
    gc_str << std::endl;
    fc_str << std::endl;
    for(std::size_t l=0;l<N_l;++l){
      gc_str << l;
      fc_str << l;
      for(std::size_t i=0;i<n_orbitals;++i){
        gc_str << " " << gc_conv[0][i][l] << " " << gc_conv[1][i][l];
        fc_str << " " << fc_conv[0][i][l] << " " << fc_conv[1][i][l];
      }
      gc_str << std::endl;
      fc_str << std::endl;
    }
    gc_str.close();
    fc_str.close();
    std::ofstream Gtl_file("Gtl.dat");
    for(std::size_t t=0;t<=N_t;++t){
      Gtl_file<<beta*t/N_t;
      for(std::size_t j=0;j<n_orbitals;++j){
        Gtl_file<<" " <<G_l_tau(t,0,0,j);
      }
      Gtl_file<<std::endl;
    }
    Gtl_file.close();
    std::ofstream Ftl_file("Ftl.dat");
    for(std::size_t t=0;t<=N_t;++t){
      Ftl_file<<beta*t/N_t;
      for(std::size_t j=0;j<n_orbitals;++j){
        Ftl_file<<" " <<F_l_tau(t,0,0,j);
      }
      Ftl_file<<std::endl;
    }
    Ftl_file.close();
    std::ofstream Gw_file("Gwl.dat");
    for(std::size_t t=0;t<N_w;++t){
      Gw_file<<(2.*t+1)*M_PI/beta;
      for(std::size_t j=0;j<n_orbitals;++j){
        Gw_file<<" "<<G_l_omega(t,0,0,j).real()<<" "<<G_l_omega(t,0,0,j).imag();
      }
      Gw_file<<std::endl;
    }
    Gw_file.close();
    std::ofstream Fw_file("Fwl.dat");
    for(std::size_t t=0;t<N_w;++t){
      Fw_file<<(2.*t+1)*M_PI/beta;
      for(std::size_t j=0;j<n_orbitals;++j){
        Fw_file<<" "<<F_l_omega(t,0,0,j).real()<<" "<<F_l_omega(t,0,0,j).imag();
      }
      Fw_file<<std::endl;
    }
    Fw_file.close();
    std::ofstream Sw_file("Swl.dat");
    for(std::size_t t=0;t<N_w;++t){
      Sw_file<<(2.*t+1)*M_PI/beta;
      for(std::size_t j=0;j<n_orbitals;++j){
        Sw_file<<" "<<S_l_omega(t,0,0,j).real()<<" "<<S_l_omega(t,0,0,j).imag();
      }
      Sw_file<<std::endl;
    }
    Sw_file.close();
  }

  //overwriting G(tau) and G(omega) in the cthyb framework with those from the legendre polys
  if (parms.exists("cthyb.DMFT_FRAMEWORK") && parms["cthyb.DMFT_FRAMEWORK"] && parms.exists("solver.OUTFILE_H5GF")){
    alps::hdf5::archive ar(parms["solver.OUTFILE_H5GF"], alps::hdf5::archive::WRITE);
    translate_Gw_to_h5gf(G_l_omega, parms).save(ar, "/G_omega");
    translate_Gt_to_h5gf(G_l_tau  , parms).save(ar, "/G_tau");
  }

}