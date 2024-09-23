#include "./evaluator.hpp"

void evaluate_2p(const alps::accumulators::result_set &results,
                 const alps::params &parms,
                 alps::hdf5::archive &solver_output){
  //write two-particle functions to text file if desired
  //compute the vertex function if needed;
  bool MEASURE_g2w=parms["cthyb.MEASURE_g2w"];
  bool MEASURE_h2w=parms["cthyb.MEASURE_h2w"];

  if(!(MEASURE_g2w || MEASURE_h2w)) return;

  //int N_w = parms["NMATSUBARA"].as<bool>();
  std::size_t N_W = parms["cthyb.N_W"].as<int>();
  std::size_t N_w2 = 2*parms["cthyb.N_w2"].as<int>();
  std::size_t n_orbitals = parms["FLAVORS"]; //number of orbitals
  bool text_output = parms["cthyb.TEXT_OUTPUT"];
  bool COMPUTE_VERTEX = parms["cthyb.COMPUTE_VERTEX"];
  double beta=parms["BETA"];

  std::ofstream g2w_str;
  std::ofstream h2w_str;
  std::ofstream gam_str;
  if(text_output){
    if(MEASURE_g2w) g2w_str.open("g2w.dat");
    if(MEASURE_h2w) h2w_str.open("h2w.dat");
    if(COMPUTE_VERTEX) gam_str.open("gammaw.dat");
  }
  std::vector<double> g2w_re;
  std::vector<double> g2w_im;
  std::vector<double> h2w_re;
  std::vector<double> h2w_im;
  std::vector<std::vector<double> >gw_re(n_orbitals);
  std::vector<std::vector<double> >gw_im(n_orbitals);
  std::vector<std::vector<double> >fw_re(n_orbitals);
  std::vector<std::vector<double> >fw_im(n_orbitals);
  std::vector<std::complex<double> > vertex;

  if(COMPUTE_VERTEX){
    vertex.resize(N_w2*N_w2*N_W);
    for(std::size_t i=0;i<n_orbitals;++i){
      std::stringstream gw_re_name; gw_re_name<<"gw_re_"<<i;
      std::stringstream gw_im_name; gw_im_name<<"gw_im_"<<i;
      std::stringstream fw_re_name; fw_re_name<<"fw_re_"<<i;
      std::stringstream fw_im_name; fw_im_name<<"fw_im_"<<i;
      gw_re[i]=results[gw_re_name.str()].mean<std::vector<double> >();
      gw_im[i]=results[gw_im_name.str()].mean<std::vector<double> >();

      if(MEASURE_h2w){
        fw_re[i]=results[fw_re_name.str()].mean<std::vector<double> >();
        fw_im[i]=results[fw_im_name.str()].mean<std::vector<double> >();
      }
    }
  }

  std::complex<double> g1,g2,g3,g4,f1;
  std::complex<double> gg,g2w,h2w,g2w_con,gamma;

  for(std::size_t i=0;i<n_orbitals;++i){
    for(std::size_t j=0;j<=i;++j){
      if(MEASURE_g2w){
        std::stringstream g2w_re_name; g2w_re_name<<"g2w_re_"<<i<<"_"<<j;
        std::stringstream g2w_im_name; g2w_im_name<<"g2w_im_"<<i<<"_"<<j;
        g2w_re=results[g2w_re_name.str()].mean<std::vector<double> >();
        g2w_im=results[g2w_im_name.str()].mean<std::vector<double> >();
      }
      if(MEASURE_h2w){
        std::stringstream h2w_re_name; h2w_re_name<<"h2w_re_"<<i<<"_"<<j;
        std::stringstream h2w_im_name; h2w_im_name<<"h2w_im_"<<i<<"_"<<j;
        h2w_re=results[h2w_re_name.str()].mean<std::vector<double> >();
        h2w_im=results[h2w_im_name.str()].mean<std::vector<double> >();
      }
      int N_wh=N_w2/2;
      for(int w2n=-N_wh;w2n<N_wh;++w2n){
        if(COMPUTE_VERTEX) g2=( w2n<0 ? conj(std::complex<double>(gw_re[i][-w2n-1],gw_im[i][-w2n-1])) : std::complex<double>(gw_re[i][w2n],gw_im[i][w2n]) );//i;w2
        for(int w3n=-N_wh;w3n<N_wh;++w3n){
          if(COMPUTE_VERTEX) g3=( w3n<0 ? conj(std::complex<double>(gw_re[j][-w3n-1],gw_im[j][-w3n-1])) : std::complex<double>(gw_re[j][w3n],gw_im[j][w3n]) );//j;w3
          for(std::size_t Wm=0;Wm<N_W;++Wm){
            int index=Wm*N_w2*N_w2 + (w2n+N_wh)*N_w2 + (w3n+N_wh);
            int w1n=w2n+Wm; int w4n=w3n+Wm;

            if(MEASURE_g2w) g2w = std::complex<double>(g2w_re[index],g2w_im[index]);
            if(MEASURE_h2w) h2w = std::complex<double>(h2w_re[index],h2w_im[index]);

            if(COMPUTE_VERTEX){
              g1=( w1n<0 ? conj(std::complex<double>(gw_re[i][-w1n-1],gw_im[i][-w1n-1])) : std::complex<double>(gw_re[i][w1n],gw_im[i][w1n]) );//i;w1
              g4=( w4n<0 ? conj(std::complex<double>(gw_re[j][-w4n-1],gw_im[j][-w4n-1])) : std::complex<double>(gw_re[j][w4n],gw_im[j][w4n]) );//j;w4
              if(MEASURE_h2w)
                f1=( w1n<0 ? conj(std::complex<double>(fw_re[i][-w1n-1],fw_im[i][-w1n-1])) : std::complex<double>(fw_re[i][w1n],fw_im[i][w1n]) );//i;w1

              gg=0.;
              if(w1n==w2n) gg+=beta*g1*g3;
              if(w2n==w3n && i==j) gg-=beta*g1*g3;

              //evaluate the vertex depending on what has been measured
              if(MEASURE_g2w && MEASURE_h2w) g2w_con=g1*h2w - f1*g2w; //usually most accurate
              else if(MEASURE_h2w)           g2w_con=(g1*h2w-f1*gg)/(1.0+f1); //somewhat less accurate; g2w not needed
              else                           g2w_con=g2w-gg; //straightforward evaluation, least accurate; h2w not needed

              gamma=g2w_con/(g1*g2*g3*g4);
              vertex[index]=gamma;
            }
            if(text_output){

              if(MEASURE_g2w)    g2w_str << "w: " << 2*w2n+1 << " wp: " << 2*w3n+1 << " W: " << Wm << " " << " i: " << i << " j: " << j << "  "
                << g2w.real() << " " << g2w.imag() << std::endl;
              if(MEASURE_h2w)    h2w_str << "w: " << 2*w2n+1 << " wp: " << 2*w3n+1 << " W: " << Wm << " " << " i: " << i << " j: " << j << "  "
                << h2w.real() << " " << h2w.imag() << std::endl;
              if(COMPUTE_VERTEX) gam_str << "w: " << 2*w2n+1 << " wp: " << 2*w3n+1 << " W: " << Wm << " " << " i: " << i << " j: " << j << "  "
                << gamma.real() << " " << gamma.imag() << std::endl;
            }
          }//Wm
        }//w3n
      }//w2n
      //write to hdf5
      if(COMPUTE_VERTEX){
        std::stringstream data_path; data_path<<"vertex_"<<i<<"_"<<j;
        solver_output<<alps::make_pvp(data_path.str(), vertex);
      }
    }//j
  }//i
}

