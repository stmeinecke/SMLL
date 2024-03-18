#include "incl/global_v2.cpp"
#include "incl/get_params.cc"
#include "incl/get_specs_sm_v4.cc"
#include "incl/get_correlations.cc"
#include "SML_Laser_parameters.cpp"
#include "SML_Laser_vars.cpp"
#include "incl/vars_vec_v8.cpp"
#include "incl/DDEintegrator_v4.cpp"
#include "incl/evaluateTS_v2.cpp"
#include "SML_Laser_eqs.cpp"


using namespace std;

int main(int argc, char* argv[]){
  
  //noise
  katana::Seed s;
  s.set_seed();
  
  //paramters -- units in ns and 10^11/cm^2 = 10^15/m^2 --> n is in units of 10^15/m^3
  parameters p;
    
  //3D density of states reservoir
  p.h_bulk = katana::getCmdOption(argv, argv+argc, "-h_bulk" , 53.033)*1e-9; // bulkthickness in 10^-9m = nm //used to be 150 with the old Lingnau fit that lacked a factor of 2^1.5 in the DOS prefactor -> 53.033*2^1.5 = 150
  p.confine = katana::getCmdOption(argv, argv+argc, "-confine" , 0.06); //SML confinement in eV
  p.mEff = katana::getCmdOption(argv, argv+argc, "-mEff" , 0.07); // in m_e
  p.temp = katana::getCmdOption(argv, argv+argc, "-temp" , 300.0); // in Kelvin
  p.NC = 2.0 * pow( (sm::m_e * p.mEff * sm::kB * p.temp)/( 2.0*M_PI*sm::hbar*sm::hbar) , 1.5) * 1e-15; // BULK density of states adjusted for 10^11/cm^2 = 10^15/m^2 --> thus in units of 10^15/m^3
  
  p.J = katana::getCmdOption(argv, argv+argc, "-J" , 240); //pump current in 10^11/cm^2
  p.n_SML = katana::getCmdOption(argv, argv+argc, "-n_SML" , 3.3); //SML area density in 10^11/cm^2
  p.gamma_N = katana::getCmdOption(argv, argv+argc, "-gamma_N" , 2.0); // reservoir 1/lifetime
  
  bool setR = katana::getCmdOption_bool(argv, argv+argc, "-R" , false); // set R manually if no dependency on J and Jtr is desired
  p.R = katana::getCmdOption(argv, argv+argc, "-R" , 135.0); //reservoir SML coupling rate
  p.R0 = katana::getCmdOption(argv, argv+argc, "-R0" , 240.0); //reservoir SML coupling rate at transparency
  double Jtr = katana::getCmdOption(argv, argv+argc, "-Jtr" , 95.0); //pump current density to reach transparency - starting value for following calculation. Can be manually set to a fixed value by using -Jtr
  
  p.gamma_rho = katana::getCmdOption(argv, argv+argc, "-gamma_rho" , 15.0); //SMLQD 1/lifetime
  
  p.g = katana::getCmdOption(argv, argv+argc, "-g" , 185.0); //gain
  p.alpha = katana::getCmdOption(argv, argv+argc, "-alpha" , 5.0); 
  p.kappa = katana::getCmdOption(argv, argv+argc, "-kappa" , 110.0); // photon 1/lifetime
  
  p.tau = katana::getCmdOption(argv, argv+argc, "-tau" , 3.8); //feedback time
  p.K = katana::getCmdOption(argv, argv+argc, "-K" , 0.0); //relative feedback strength
  p.C = katana::getCmdOption(argv, argv+argc, "-C" , 0.0); //relative feedback phase
  
  p.fakeSpE = katana::getCmdOption(argv, argv+argc, "-fakeSpe" , 1.0E-50);  
  
  
  //////////////////////////////////////////
  //n at the transparency level rho = 0.5
  //////////////////////////////////////////
  
  auto calc_ntr = [&](){

  double rho_eq_goal = 0.5*(p.gamma_rho/p.R0 + 1.0);
  
  double ntest = 0.0; 
  double nstep = 1e7; // <------ hardcoded stepsize
  double rho_eq_tmp = get_rho_eq(ntest, p.NC, p.confine, p.temp);
  
  while(rho_eq_tmp < rho_eq_goal){
    ntest += nstep;
    rho_eq_tmp = get_rho_eq(ntest, p.NC, p.confine, p.temp);
  }
  while(fabs(rho_eq_tmp - rho_eq_goal) > 1e-10){ // <------ threshold here
    nstep *= 0.5;
    if(rho_eq_tmp > rho_eq_goal){
      ntest -= nstep;
      rho_eq_tmp = get_rho_eq(ntest, p.NC, p.confine, p.temp);
      if(rho_eq_tmp > rho_eq_goal) nstep *= 2.0;
    }
    else{
      ntest += nstep;
      rho_eq_tmp = get_rho_eq(ntest, p.NC, p.confine, p.temp);
      if(rho_eq_tmp < rho_eq_goal) nstep *= 2.0;
    }
//       cout << rho_eq_goal << " " << rho_eq_tmp << " " << ntest << endl;
  }
  
  return ntest;
  };
  
  
  auto calcParas = [&](){
    if(!setR) p.R = p.R0 * p.J * p.J / (Jtr*Jtr);
    p.n_tr = calc_ntr();
    cout << "n_tr: " << p.n_tr << endl;
  };

  
  
  //solver
  DDEintegrator DDEsolver;
  //solver parameters
  unsigned long long int tn = 0;
  double intTime = katana::getCmdOption(argv, argv+argc, "-intTime" , 10.0);
  double dt = katana::getCmdOption(argv, argv+argc, "-dt" , 0.001);
  double outTime = katana::getCmdOption(argv, argv+argc, "-outTime" , 10.0);
  if(outTime < 0.0) outTime = intTime;
  unsigned long long int outTime_ntn = (unsigned long long int)(outTime/dt);
  
  
  //dynamical variables with delay
  vars_vec_wdX Xhist(p.tau,dt);
  
  
  //////////////////////////////////////////
  //noise
  //////////////////////////////////////////
  
  //noise
  p.SqrtNoiseStr = sqrt( katana::getCmdOption(argv, argv+argc, "-noiseStr" , 1.0E-5) / dt );
  
  std::function<void (parameters*)> noise = noise_empty;
  
  std::function<void (parameters*)> GWnoise = [&](parameters* p){ 
      p->noise = p->SqrtNoiseStr*katana::gwNoise();
  };
  //gw noise
  if(katana::getCmdOption_bool(argv, argv+argc, "-noise" , false)) noise = GWnoise;
  
    
  //////////////////////////////////////////
  //calc J as a funktion of rho
  //////////////////////////////////////////
  
  //only below threshold - uses a divide and conquer method with hard coded steps and thresholds...
  
  //watch out - calculation only for a fixed scattering rate R that does NOT depend on J.
  auto calc_JofRho_fixedR = [&](double rholevel, double R){ // <------ R generally depends on J

    double rho_eq_goal = rholevel*(1.0 + p.gamma_rho/R); 
    
    double ntest = 0.0; 
    double nstep = 1e7; // <------ hardcoded stepsize
    double rho_eq_tmp = get_rho_eq(ntest, p.NC, p.confine, p.temp);
    
    while(rho_eq_tmp < rho_eq_goal){
      ntest += nstep;
      rho_eq_tmp = get_rho_eq(ntest, p.NC, p.confine, p.temp);
    }
    while(fabs(rho_eq_tmp - rho_eq_goal) > 1e-10){ // <------ threshold here
      nstep *= 0.5;
      if(rho_eq_tmp > rho_eq_goal){
        ntest -= nstep;
        rho_eq_tmp = get_rho_eq(ntest, p.NC, p.confine, p.temp);
        if(rho_eq_tmp > rho_eq_goal) nstep *= 2.0;
      }
      else{
        ntest += nstep;
        rho_eq_tmp = get_rho_eq(ntest, p.NC, p.confine, p.temp);
        if(rho_eq_tmp < rho_eq_goal) nstep *= 2.0;
      }
//       cout << rho_eq_goal << " " << rho_eq_tmp << " " << ntest << endl;
    }
    
    double Jtmp = p.gamma_N * ntest * p.h_bulk + 2.0*p.n_SML * p.gamma_rho * rholevel;
    
    return Jtmp;
  };
  
  //iterative method to determine J of rho if R depends on J
  auto calc_JofRho = [&](double rholevel){
  
    double Jtest, Jtemp, Rtest;
    Jtest = p.J;
    
    do{
      Jtemp = Jtest;
      Rtest = p.R0 * Jtest*Jtest/(Jtr*Jtr); // R is recalculated in every loop using the last J
      Jtest = calc_JofRho_fixedR(rholevel,Rtest);
//       cout << "Jtest: " << Jtest << endl;
    } while( fabs(Jtemp-Jtest)/Jtest > 1e-8); // <------ threshold here
    
    return Jtest;
  };
  
  
  //calculates Jtr unless the flags -Jtr or -R are used
  if(!setR && !katana::getCmdOption_bool(argv, argv+argc, "-Jtr" , false)){
    //iterative method to converge to the true Jtr as the current value of Jtr affects the calculation of the new Jtr.
    double tmpJtr;
    do{
      tmpJtr = Jtr;
      calcParas();
      Jtr = calc_JofRho(0.5);
//       cout << "Jtr " << Jtr << endl;
    } while( fabs(tmpJtr - Jtr)/Jtr > 1e-6); // <------ threshold here
    
    calcParas();
    cout << "calculated transparency pump current: Jtr=" << Jtr << endl;
  }
  
  if(katana::getCmdOption_bool(argv, argv+argc, "-Jth" , false)){
    double rhothr = 0.5*(p.kappa/p.g + 1.0);
    double Jth = calc_JofRho( rhothr );
    cout << "calculated threshold pump current: Jth=" << Jth << endl;
  }

  if(katana::getCmdOption_bool(argv, argv+argc, "-P" , false)){ // set the pump current J in units of Jth
    double P = katana::getCmdOption(argv, argv+argc, "-P" , 1.5);  
    double rhothr = 0.5*(p.kappa/p.g + 1.0);
    double Jth = calc_JofRho( rhothr );
    cout << "calculated threshold pump current: Jth=" << Jth << endl;
    p.J = P*Jth;
    cout << "calculated pump current: J=" << p.J << endl;
    calcParas();
  }
  
  
  //////////////////////////////////////////
  //outputfunctions
  //////////////////////////////////////////
  
  //for timeseries output
  outputfile.precision(10); //outputfile must be declared in global.hpp
  
  //empty output function
  auto empty = [&](vars* X, unsigned long long int tn, unsigned long long int tn_final){};
  
  //output timeseries of time and all dynamical variables -> very slow and might take lots  of memory
  auto outputAllToFile = [&](vars* X, unsigned long long int tn, unsigned long long int tn_final){
    int k = 0;
    int kmax = 1;
    if(tn >= tn_final - outTime_ntn && k == 0){
      outputfile << tn*dt << "\t";
      outputfile << norm(X->E) << "\t" << X->rho << "\t" << X->N << "\t" << X->N/p.NC;
      outputfile << std::endl;
    }
    k = (k+1)%kmax;
  };
  
  //output to vectors for TS analysis
  std::vector<double> outVec;
  outVec.reserve((int)(outTime/dt));
  
  std::vector<std::complex<double>> outVecComplex;
  outVecComplex.reserve((int)(outTime/dt));
  
  
  auto outputToVector = [&](vars* X, unsigned long long int tn, unsigned long long int tn_final){
    if(tn >= tn_final - outTime_ntn){
      outVec.push_back(norm(X->E));
      outVecComplex.push_back(X->E);
    }
  };
  
  
  //////////////////////////////////////////
  //for sweeps
  //////////////////////////////////////////
  
  std::ofstream out_sweep_State, out_sweep_Ext_GS_max, out_sweep_Ext_GS_min, out_sweep_TS;
  std::ofstream out_sweep_AllExt_GS_max, out_sweep_AllExt_GS_min;
  out_sweep_State.precision(12), out_sweep_TS.precision(12), out_sweep_Ext_GS_min.precision(12);
  out_sweep_AllExt_GS_max.precision(12), out_sweep_AllExt_GS_min.precision(12);
  std::ofstream out_sweep_Ext_freq_max, out_sweep_Ext_freq_min;
  out_sweep_Ext_freq_max.precision(12), out_sweep_Ext_freq_min.precision(12);
  evaluateTS eval_sweep;
  evaluateTS eval_sweep_freq;
  
  double sweep_doubleCountTol = katana::getCmdOption(argv, argv+argc, "-SDoubleCountTol" , 1E-2);
  eval_sweep.doubleCountTol = sweep_doubleCountTol;
  eval_sweep_freq.doubleCountTol = sweep_doubleCountTol;
  
  double* PtrToSweepPar;
  std::string str_sweep_parameter;
  double* PtrToSecPar;
  std::string str_sweep_sec_parameter;
  
  double sweep_IntTime = katana::getCmdOption(argv, argv+argc, "-sIntTime" , 1000.0);
  double sweep_OutTime = katana::getCmdOption(argv, argv+argc, "-sOutTime" , 100.0);
  int nsweep_steps = katana::getCmdOption(argv, argv+argc, "-sSteps" , 100);
  double sweep_start = katana::getCmdOption(argv, argv+argc, "-sStart" , 0);
  double sweep_end = katana::getCmdOption(argv, argv+argc, "-sEnd" , 1);
  double sweep_step = (sweep_end-sweep_start)/(double)nsweep_steps;
  
  std::vector<double> sweep_Phase, sweep_Frequency, sweep_rollingMeanFreq;
  sweep_Phase.reserve((int)(outTime/dt)), sweep_Frequency.reserve((int)(outTime/dt)), sweep_rollingMeanFreq.reserve((int)(outTime/dt));
  
  std::string str_sweep_data = "data/";
  std::string str_sweep_updown = "up";
  
  std::string str_state;
  std::string str_extrema_GS;
  std::string str_extrema_freq;
  std::string str_powerSpec;
  std::string str_opticalSpec;
  std::string str_AC;
  std::string str_TS;
  std::string str_bin;
  
  bool bool_wdownsweep = katana::getCmdOption_bool(argv, argv+argc, "-wdownSweep" , false);
  bool bool_powerSpec = katana::getCmdOption_bool(argv, argv+argc, "-wpowerSpec" , false);
  bool bool_opticalSpec = katana::getCmdOption_bool(argv, argv+argc, "-wopticalSpec" , false);
  bool bool_AC = katana::getCmdOption_bool(argv, argv+argc, "-wAC" , false);
  bool bool_TS = katana::getCmdOption_bool(argv, argv+argc, "-wTS" , false);
  bool bool_Extrema = katana::getCmdOption_bool(argv, argv+argc, "-wExtrema" , false);
  bool bool_AllExtrema = katana::getCmdOption_bool(argv, argv+argc, "-wAllExtrema" , false);
  bool bool_bin = katana::getCmdOption_bool(argv, argv+argc, "-wBin" , false);
  bool bool_LogSweep = katana::getCmdOption_bool(argv, argv+argc, "-wLogSweep" , false);
  


  
  
  auto sweep = [&](){
    
    std::vector<double> SweepVals(nsweep_steps+1);
    SweepVals[0] = sweep_start;

    if(bool_LogSweep){
      double SweepLogMult = pow( (sweep_end/sweep_start), 1.0/((double)(nsweep_steps)) );
      for(unsigned int k = 1; k < SweepVals.size(); k++){
        SweepVals[k] = SweepVals[k-1]*SweepLogMult;
//         cout << k << " " <<  SweepVals[k] << endl;
      }
    }
    else{
      double SweepLinAdd = (sweep_end-sweep_start)/(double)nsweep_steps;
      for(unsigned int k = 1; k < SweepVals.size(); k++){
        SweepVals[k] = SweepVals[k-1]+SweepLinAdd;
//         cout << k << " " <<  SweepVals[k] << endl;
      }
    }
    
    str_state = str_sweep_data + str_sweep_updown + "_state_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_extrema_GS = str_sweep_data + str_sweep_updown + "_extrema_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_extrema_freq = str_sweep_data + str_sweep_updown + "_extrema_freq_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_powerSpec = str_sweep_data + "/powerSpec/powerSpec_"+ str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_opticalSpec = str_sweep_data + "/opticalSpec/opticalSpec_"+ str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_AC = str_sweep_data + "AC/AC_" + str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_TS = str_sweep_data + "TS/TS_" + str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    str_bin = str_sweep_data + "bin/bin_" + str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + to_string(*PtrToSecPar);
    
    if(bool_Extrema){
      out_sweep_Ext_GS_max.open(str_extrema_GS + "_max");
      out_sweep_Ext_GS_min.open(str_extrema_GS + "_min");
      out_sweep_Ext_freq_max.open(str_extrema_freq + "_max");
      out_sweep_Ext_freq_min.open(str_extrema_freq + "_min");
    }
    if(bool_AllExtrema){
      out_sweep_AllExt_GS_max.open(str_extrema_GS + "_All_max");
      out_sweep_AllExt_GS_min.open(str_extrema_GS + "_All_min");
    }
    out_sweep_State.open(str_state);
    
    outTime = sweep_OutTime;
    outTime_ntn = (unsigned long long int)(outTime/dt);
//     *PtrToSweepPar = sweep_start;
    
    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string binhist = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(binhist);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(1E-9);
    }
    
    for(unsigned int k = 0; k < SweepVals.size(); k++){
      *PtrToSweepPar = SweepVals[k];
      
      sm::clck = clock();
      calcParas();
      tn=0;
//       Xhist.setToCnst(1E-6);
      outVec.resize(0), outVecComplex.resize(0);
      DDEsolver.DDE_RK4(SMLL, algs_empty, noise, &Xhist, &p, &tn, sweep_IntTime, dt, outputToVector);
      eval_sweep.newTS(&outVec, dt);
      
      //analize optical frequency
      
      sweep_Phase.resize(0), sweep_Frequency.resize(0), sweep_rollingMeanFreq.resize(0);
      double smfreq =  eval_sweep.opticalFreq(outVecComplex, sweep_Phase, sweep_Frequency, dt);
      int sk_tau = p.tau / dt;
      
      for(int k = sk_tau; k < sweep_Phase.size()-2; k++){
        sweep_rollingMeanFreq.push_back( (sweep_Phase[k] - sweep_Phase[k-sk_tau])/p.tau );
      }
      eval_sweep_freq.newTS(&sweep_rollingMeanFreq,dt);
      
//       out_sweep_TS.open(str_TS+"_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar));
//       for(int k = 0; k < sweep_rollingMeanFreq.size(); k++) out_sweep_TS << dt*k << "\t" << sweep_rollingMeanFreq[k] << std::endl;
//       out_sweep_TS.close();
 
      cout << str_sweep_parameter << ": " << *PtrToSweepPar;
      out_sweep_State << *PtrToSweepPar << "\t" << *PtrToSecPar << "\t"; //2
      out_sweep_State << eval_sweep.average << "\t" << eval_sweep.greatestMax << "\t" << eval_sweep.smallestMin << "\t"; //5
      out_sweep_State << eval_sweep.numberOfMax() << "\t" << eval_sweep.numberOfMin() << "\t"; //7
      
      out_sweep_State << smfreq << "\t"; //8
      out_sweep_State << eval_sweep_freq.average << "\t" << eval_sweep_freq.numberOfMax() << "\t" << eval_sweep_freq.numberOfMin() << "\t" << eval_sweep_freq.loopsRatio() << "\t"; //12
      out_sweep_State << Xhist.at_t0_rPtr()->rho << "\t" << Xhist.at_t0_rPtr()->N << "\t"; //14
      out_sweep_State << endl;
      
      if(bool_Extrema){
        //GS
        out_sweep_Ext_GS_max << *PtrToSweepPar << "\t" << eval_sweep.greatestMax << std::endl;
        for(int k = 0; k < eval_sweep.uniqueMax.size(); k++) out_sweep_Ext_GS_max << *PtrToSweepPar << "\t" << eval_sweep.uniqueMax[k] << std::endl;      
        out_sweep_Ext_GS_min << *PtrToSweepPar << "\t" <<  eval_sweep.smallestMin << std::endl;
        for(int k = 0; k < eval_sweep.uniqueMin.size(); k++) out_sweep_Ext_GS_min << *PtrToSweepPar << "\t" << eval_sweep.uniqueMin[k] << std::endl;
        //Frequency
        out_sweep_Ext_freq_max<< *PtrToSweepPar << "\t" << eval_sweep_freq.greatestMax << std::endl;
        for(int k = 0; k < eval_sweep_freq.uniqueMax.size(); k++) out_sweep_Ext_freq_max << *PtrToSweepPar << "\t" << eval_sweep_freq.uniqueMax[k] << std::endl;      
        for(int k = 0; k < eval_sweep_freq.uniqueMin.size(); k++) out_sweep_Ext_freq_min << *PtrToSweepPar << "\t" << eval_sweep_freq.uniqueMin[k] << std::endl;
        out_sweep_Ext_freq_min << *PtrToSweepPar << "\t" <<  eval_sweep_freq.smallestMin << std::endl;
      }
      if(bool_AllExtrema){
        //GS
        out_sweep_AllExt_GS_max << *PtrToSweepPar << "\t" << eval_sweep.greatestMax << std::endl;
        for(int k = 0; k < eval_sweep.maxima.size(); k++) out_sweep_AllExt_GS_max << *PtrToSweepPar << "\t" << eval_sweep.maxima[k] << std::endl;      
        out_sweep_AllExt_GS_min << *PtrToSweepPar << "\t" <<  eval_sweep.smallestMin << std::endl;
        for(int k = 0; k < eval_sweep.minima.size(); k++) out_sweep_AllExt_GS_min << *PtrToSweepPar << "\t" << eval_sweep.minima[k] << std::endl;
      }
      if(bool_powerSpec){
        std::vector<double> powerSpec;
        std::vector<double> WindowedData;
        //GS
        sm::window_Hann(outVec, WindowedData);
        sm::get_power_spec(WindowedData, powerSpec, dt);
        std::ostringstream sw_num;
        sw_num << std::fixed << std::setprecision(12);
        sw_num << *PtrToSweepPar;
        sm::dump_power_spec(powerSpec, dt, 10.0, str_powerSpec+"_Hann_"+str_sweep_parameter+"_"+sw_num.str());
//         sm::dump_power_spec(powerSpec, dt, 1.0, str_powerSpec+"_GS_windowed_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar));
//         //no window
//         sm::get_power_spec(outVec, powerSpec, dt);
//         sm::dump_power_spec(powerSpec, dt, powerSpec.size()/50, str_powerSpec+"_GS_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar));
      }
      if(bool_opticalSpec){
        std::vector<double> opticalSpec;
        std::ostringstream sw_num;
        sw_num << std::fixed << std::setprecision(12);
        sw_num << *PtrToSweepPar;
        
        std::vector<std::complex<double>> WindowedData;
        sm::window_Hann(outVecComplex, WindowedData);
        sm::get_optical_spec(WindowedData, opticalSpec);
        sm::dump_optical_spec(opticalSpec, dt, 4, str_opticalSpec+"_Hann_"+str_sweep_parameter+"_"+sw_num.str());
//         //no window
//         sm::get_optical_spec(outVecComplex, opticalSpec);
//         sm::dump_optical_spec(opticalSpec, dt, 4, str_opticalSpec+"_Hann_"+str_sweep_parameter+"_"+sw_num.str());
      }
      if(bool_TS){
        out_sweep_TS.open(str_TS+"_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar));
        for(int k = 0; k < outVec.size()/10; k+=2) out_sweep_TS << dt*k << "\t" << outVec[k] << std::endl;
        out_sweep_TS.close();
      }
      if(bool_AC){
        std::vector<double> AC;
        //GS
        katana::get_autocorr_norm_out(outVecComplex, AC);
        katana::print_correlation(AC, dt, str_AC+"_GS_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar), (int)(1000.0/dt), 1, "ps");
      }
      if(bool_bin){
        Xhist.tn = tn;
        Xhist.save(str_bin+"_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar));
        cout << str_bin+"_"+str_sweep_parameter+"_"+to_string(*PtrToSweepPar) << endl;
      }
      
      cout << ", integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
    }
    if(bool_Extrema){
      out_sweep_Ext_GS_max.close();
      out_sweep_Ext_GS_min.close();
      out_sweep_Ext_freq_max.close();
      out_sweep_Ext_freq_min.close();
    }
    if(bool_AllExtrema){
      out_sweep_AllExt_GS_max.close();
      out_sweep_AllExt_GS_min.close();
    }
    out_sweep_State.close();
  };
 
  
  
  //timeseries to file
  if(katana::getCmdOption_bool(argv, argv+argc, "-simpleTS" , false)){
    outputfile.open("data/simpleTS");
//     tn=0;
    
    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string binhist = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(binhist);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(1E-9);
    }
    
    sm::clck = clock();
    DDEsolver.DDE_RK4(SMLL, algs_empty, noise, &Xhist, &p, &tn, intTime, dt, outputAllToFile);
    cout << "integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
    outputfile.close();
    
    //save history to binary file
    if(katana::getCmdOption_bool(argv, argv+argc, "-saveHist" , false)){
      Xhist.tn = tn;
      Xhist.save("data/Xhist.bin");
    }
  }
  
  
  
  //compute and analize time series
  if(katana::getCmdOption_bool(argv, argv+argc, "-TS" , false)){
    
    //evaluateTS objects
    evaluateTS eval;
    evaluateTS eval_avfreq;
    outVec.resize(0), outVecComplex.resize(0);
    
    std::vector<double> Phase, Frequency;
    
    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string binhist = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(binhist);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(0.001);
    }

    sm::clck = clock();
//     DDEsolver.DDE_euler(SMLL, algs_empty, noise, &Xhist, &p, &tn, intTime, dt, outputToVector);
    DDEsolver.DDE_RK4(SMLL, algs_empty, noise, &Xhist, &p, &tn, intTime, dt, outputToVector);
    cout << "integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
    //load new TS
    eval.newTS(&outVec, dt);
    
    cout << "----results----" << endl;
    cout <<"average: " << eval.average << endl;
    cout <<"greates Max: " << eval.greatestMax << ", smallest Min: " << eval.smallestMin << endl;
    cout <<"number of Max: " << eval.numberOfMax() <<", number of Min: " << eval.numberOfMin() << endl;
    cout <<"period: " << eval.period << endl;
    if(katana::getCmdOption_bool(argv, argv+argc, "-TSwOptFreq" , false)){
      double avFreq = eval.opticalFreq(outVecComplex, Phase, Frequency, dt);
      int k_tau = p.tau / dt;
      std::vector<double> avfreq;
      avfreq.reserve(Phase.size());
      for(int k = k_tau; k < Phase.size()-2; k++){
        avfreq.push_back( (Phase[k] - Phase[k-k_tau])/p.tau );
      }
      eval_avfreq.newTS(&avfreq,dt);
      cout <<"average optical frequency: " << avFreq << endl;
      cout <<"Rolling mean frequency: mean: " << eval_avfreq.average << " greatest max: " << eval_avfreq.greatestMax << " smallest min: " << eval_avfreq.smallestMin << " number of max: " << eval_avfreq.numberOfMax() << " number of min: " << eval_avfreq.numberOfMin() << endl;
    }
    cout << "----^^^^^^----" << endl;

//     double Var = 0.0;
//     for(unsigned int k = 0; k < outVec.size(); k++){
//       Var += (outVec[k]-eval.average)*(outVec[k]-eval.average);
//     }
//     Var /= (double)outVec.size();
//     double stdDev = sqrt(Var);
//     cout << Intensity variance << " " << stdDev << " " << stdDev / eval.average <<  endl;
    
    outputfile.open("data/out_TS");
    for(int k = 0; k < outVec.size()-2; k++){
      outputfile << dt*k << "\t" << outVec[k] << "\t";
//       outputfile << Phase[k] << "\t" << Frequency[k];
      outputfile << endl;
    }
    outputfile.close();
    
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-powerSpec" , false)){
      std::vector<double> powerSpec;
      std::vector<double> WindowedData;
      
      sm::window_Hann(outVec, WindowedData);
      sm::get_power_spec(WindowedData, powerSpec, dt);
      sm::dump_power_spec(powerSpec, dt, 10, "data/out_TS_powerSpec_Hann");
      //no window
      sm::get_power_spec(outVec, powerSpec, dt);
      sm::dump_power_spec(powerSpec, dt, 10, "data/out_TS_powerSpec");
    }
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-opticalSpec" , false)){
      std::vector<double> opticalSpec;
      std::vector<std::complex<double>> WindowedData;
      sm::window_Hann(outVecComplex, WindowedData);
      sm::get_optical_spec(WindowedData, opticalSpec);
      sm::dump_optical_spec(opticalSpec, dt, 10, "data/out_TS_opticalSpec_Hann");
      //no window
      sm::get_optical_spec(outVecComplex, opticalSpec);
      sm::dump_optical_spec(opticalSpec, dt, 10, "data/out_TS_opticalSpec");
    }
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-AC" , false)){
      std::vector<double> AC;
      katana::get_autocorr_norm_out(outVecComplex, AC);
      katana::print_correlation(AC, dt, "data/out_TS_AC", AC.size(), 1, "ns");
    }
    if(katana::getCmdOption_bool(argv, argv+argc, "-IntAC" , false)){
      std::vector<double> IntAC;
      katana::get_autocorr(outVec, IntAC);
      katana::print_correlation(IntAC, dt, "data/out_TS_IntAC", IntAC.size(), 1, "ns");
    }
    
    //save history to binary file
    if(katana::getCmdOption_bool(argv, argv+argc, "-saveHist" , false)){
      Xhist.tn = tn;
      Xhist.save("data/Xhist.bin");
    }
    
  }
  
  //K sweep with C
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_K_wC" , false)){
    tn= 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.K;
    str_sweep_parameter = "K";
    PtrToSecPar = &p.C;
    str_sweep_sec_parameter = "C";
    
    str_sweep_updown = "up";
    sweep();

    if(bool_wdownsweep){
      sweep_start = sweep_end;
      sweep_step = -sweep_step;
      
      str_sweep_updown = "down";
      sweep();
    }
  }
  
  
  //K sweep with J
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_K_wJ" , false)){
    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.K;
    str_sweep_parameter = "K";
    PtrToSecPar = &p.J;
    str_sweep_sec_parameter = "J";
    
    str_sweep_updown = "up";
    sweep();

    if(bool_wdownsweep){
      sweep_start = sweep_end;
      sweep_step = -sweep_step;
      
      str_sweep_updown = "down";
      sweep();
    }
  }
  
  //J sweep with g
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_J_wg" , false)){
    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.J;
    str_sweep_parameter = "J";
    PtrToSecPar = &p.g;
    str_sweep_sec_parameter = "g";
    
    str_sweep_updown = "up";
    sweep();

    if(bool_wdownsweep){
      sweep_start = sweep_end;
      sweep_step = -sweep_step;
      
      str_sweep_updown = "down";
      sweep();
    }
  }


 
  return 0;
}



