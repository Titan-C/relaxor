#include "experiment.h"

void Gen_exp(unsigned int L, unsigned int numexps, std::vector<double> rho, std::vector<double>& Tdat,
	     std::vector<double>& Fields, std::vector<double> tau, std::string Exp_ID)
{
  Material relaxor(L);
  unsigned int Equi_iter=350;
  //   Código para bajar (variar) la temperatura a campos fijos
  if (Exp_ID == "cool" || Exp_ID == "heat"){
    for(unsigned int p=0; p<rho.size(); p++){
      std::vector<double> Thermostat;
      step2vec(Tdat[0],Tdat[2],Tdat[1],Thermostat,1);
      for(unsigned int t=0; t< tau.size() ; t++){
	unsigned int Exp_iter = stepEstimator(3000,tau[t],2);

	for(unsigned int E=0; E<Fields.size(); E++){
	  std::ostringstream id_proc;
	  id_proc<<Exp_ID<<"_p"<<rho[p]<<"_E"<<Fields[E]<<"_t"<<tau[t]<<"_L"<<L<<"_n"<<numexps;
	  id_proc<<"_Ti"<<Thermostat[0]<<"Tf"<<Tdat[1]<<"dT"<<Tdat[0]<<"_X"<<Exp_iter<<"_Q"<<Equi_iter;

	  std::cout<<id_proc.str()<<":";
	  clock_t cl_start = clock();
	  unsigned int sim_size = sizeof(double)*Exp_iter*Thermostat.size()*numexps;
	  if (needSimulation(id_proc.str(), sim_size)) {
	    std::vector<double> Equi_field = wavearray(Fields[E],tau[t],Equi_iter, Equi_iter);
	    std::vector<double> Exp_field  = wavearray(Fields[E],tau[t],Exp_iter ,0    );
	    std::cout<<Equi_field.size()<<" "<<Exp_field.size()<<" ";

	    for(unsigned int n=0; n<numexps; n++){
	      relaxor.init(rho[p], id_proc.str(),false);
	      for(unsigned int T=0; T<Thermostat.size(); T++){
		relaxor.state(Thermostat[T], Equi_field ,false);
		relaxor.state(Thermostat[T], Exp_field  ,true);
	      }}
	  }
	  proces_data(Thermostat,Fields[E],tau[t],numexps,L*L*L, rho[p],Exp_iter,id_proc.str());
	  std::cout<<(double) (clock()-cl_start)/CLOCKS_PER_SEC<<"\n";
	}}
    }
  }
}

void proces_data(std::vector< double >& Temps, double Field,
		 unsigned int tau, unsigned int numexps, unsigned int PNR,
		 double p, unsigned int Niter, std::string id_proc){
  //procesar los datos para casos de variación de temperatura
    std::vector<double> pol_stats, pol_int_avg;
    pp_data(pol_stats,pol_int_avg,Temps.size(),numexps,tau,Niter,id_proc);
    std::vector<double> intfield (1,Field);
    eval_pol(pol_stats,numexps,Temps,id_proc,true);
    calc_sus(pol_int_avg,numexps,Temps,intfield,id_proc);
//     eval_frozen(PNR,Niter, Temps, numexps, id_proc);
    intfield.clear();
    pol_int_avg.clear();
    pol_stats.clear();
}

void pp_data(std::vector<double>& pol_stats, std::vector<double>& pol_int_avg, unsigned int data_length,
	      unsigned int numexps, unsigned int tau, unsigned int Niter, std::string id_proc){
  //Generar vectores de Datos
  double * pol_hist = new double[Niter];
  unsigned int dat_vec_size = data_length*numexps*2;
  pol_stats.resize(dat_vec_size);
  pol_int_avg.resize(dat_vec_size);

  /*Generar arreglo de peso sin, cos para la integral*/
  std::vector<double> cos_wave, sin_wave;
  cos_wave = wavearray(1.0,tau,Niter,0);
  sin_wave = wavearray(1.0,tau,Niter,tau/4.0);

  //Abrir Archivo y leer
  std::string name = "log_pol_"+id_proc+".dat";
  std::ifstream file(name.c_str());
  for(unsigned int ind=0; ind<dat_vec_size; ind+=2){
      file.read((char *)&pol_hist[0],Niter*sizeof(double));
      /*Calcular media y desviación estandar de polarización por
       numero y condiciones de experimento*/
      pol_stats[ind]=gsl_stats_mean(pol_hist,1,Niter);
      pol_stats[ind+1]=gsl_stats_sd_m(pol_hist,1,Niter, pol_stats[ind]);
      /*Integración por Simpson, para promedio pesado */
      pol_int_avg[ind]=simpson_int(pol_hist,cos_wave)/Niter;
      pol_int_avg[ind+1]=simpson_int(pol_hist,sin_wave)/Niter;
  }
  //liberar memoria
  cos_wave.clear();
  sin_wave.clear();
  delete[] pol_hist;
  file.close();
}
void eval_frozen(unsigned int PNR, unsigned int Niter, const std::vector<double>& Temps, unsigned int numexps, std::string id_proc){
  int * sigma_hist = new int[PNR];
  std::vector< double > sigmaTemp;

  //Abrir Archivo, leer guardar datos
  std::string name = "log_sigma_"+id_proc+".dat";
  std::ifstream file(name.c_str());
  for(unsigned int n=0; n<numexps ; n++){
    for(unsigned int T=0; T<Temps.size(); T++){
	file.read((char *)&sigma_hist[0],PNR*sizeof(int));
	for(unsigned int s = 0; s<PNR ; s++)
	  sigmaTemp[n*Temps.size()*PNR+T*PNR + s] += sigma_hist[s];
  }}
  array_print(sigmaTemp, "sigmas_"+id_proc+".dat", PNR, Niter);
  delete[] sigma_hist;
  //Evaluar % congelamiento
  std::vector< std::vector<double> > Frozen;
  Frozen.resize(Temps.size());
  for(unsigned int T=0; T<Temps.size(); T++){
    Frozen[T].assign(4,0);
    Frozen[T][0]=Temps[T];
    for(unsigned int s = 0; s<PNR ; s++){
      double avgabssigma= std::abs(sigmaTemp[T*PNR + s]/Niter);
      if (avgabssigma > 0.9){
	Frozen[T][3]++;
	if (avgabssigma > 0.95){
	  Frozen[T][2]++;
	  if (avgabssigma > 0.999999)
	    Frozen[T][1]++;
	}}
    }
    for(unsigned int f=1; f<4;f++)
      Frozen[T][f]/=PNR;
  }
  array_print(Frozen,"fro_"+id_proc+".dat");
  /*Liberar memoria*/
  sigmaTemp.clear();
  for(unsigned int T=0; T<Temps.size(); T++)
    Frozen[T].clear();
  Frozen.clear();
}
void eval_pol(const std::vector<double>& pol_stats, unsigned int numexps, const std::vector<double>& x_array, std::string id_proc, bool absolut) {

  //Polarización, o polarización absoluta y desviación estandar
  double * data_array = new double [numexps];
  double data_length = x_array.size();
  std::vector< std::vector<double> > pol_final;
  pol_final.resize(data_length);
  for(unsigned int x=0;x< data_length; x++){
    pol_final[x].resize(3);
    pol_final[x][0]=x_array[x];

    double pol_std=0;
    for(unsigned int n=0;n<numexps;n++){
      unsigned int ind = 2*(n*data_length+x);
      data_array[n]=(absolut) ? std::abs(pol_stats[ind]) : pol_stats[ind];
      pol_std+=pol_stats[ind+1]*pol_stats[ind+1];
    }
    pol_final[x][1]=gsl_stats_mean(data_array,1,numexps);
    pol_final[x][2]=sqrt(pol_std/numexps);
  }
  delete[] data_array;
  array_print(pol_final,"pol_"+id_proc+".dat");

  /*Liberar memoria*/
  for(unsigned int i=0;i<data_length;i++)
    pol_final[i].clear();
  pol_final.clear();
}

void calc_sus(const std::vector<double>& pol_int_avg, unsigned int numexps,
	      const std::vector<double>& x_array, const std::vector<double>& campo, std::string id_proc){

  /*Calcular susceptibildad más error*/
  double * data_arrayr = new double [numexps];
  double * data_arrayi = new double [numexps];
  bool fieldvec = (campo.size()>1) ? true : false;
  double data_length = x_array.size();
  std::vector< std::vector<double> > X_mat;
  X_mat.resize(data_length);
  for(unsigned int x=0;x<data_length;x++){
    X_mat[x].resize(7);
    X_mat[x][0]=x_array[x];

    for(unsigned int n=0;n<numexps;n++){
      unsigned int ind = 2*(n*data_length+x);
      unsigned int field_ind=(fieldvec) ? x : 0;
      data_arrayr[n]=pol_int_avg[ind]/campo[field_ind];
      data_arrayi[n]=pol_int_avg[ind+1]/campo[field_ind];
    }
    X_mat[x][1]=gsl_stats_mean(data_arrayr,1,numexps);
    X_mat[x][3]=gsl_stats_sd_m(data_arrayr,1,numexps,X_mat[x][1]);
    X_mat[x][2]=0.0001/X_mat[x][3]/X_mat[x][3];
    X_mat[x][4]=gsl_stats_mean(data_arrayi,1,numexps);
    X_mat[x][6]=gsl_stats_sd_m(data_arrayi,1,numexps,X_mat[x][4]);
    X_mat[x][5]=0.0001/X_mat[x][6]/X_mat[x][6];
  }
  array_print(X_mat, "sus_"+id_proc+".dat");

  //liberar memoria
  for(unsigned int i=0; i<X_mat.size();i++)
    X_mat[i].clear();
  X_mat.clear();
  delete[] data_arrayr;
  delete[] data_arrayi;
}

void step2vec(double v_start, double v_end, double dv, std::vector<double>& array, double unidad){
  while(v_start<=v_end) {
    array.push_back(v_start*unidad);
    v_start+=dv;
  }
  while (v_start>v_end) {
    array.push_back(v_start*unidad);
    v_start-=dv;
  }
}

std::vector<double> loop2vec(double max, int divs, double unidad){
  std::vector<double> vec;
  double step = 1.0 / divs;
  unsigned int l=0;
  for(double x=step; x<=1 ;x+=step){
    vec.push_back(x*x*max*unidad);
    l++;
  }
  for(unsigned int i=l-1; i>0; i--)
    vec.push_back(vec[i-1]);

  unsigned int vec_size= vec.size();
  for(unsigned int i=0; i<vec_size; i++)
    vec.push_back(-1*vec[i]);

  for(unsigned int i=0; i<=l;i++)
    vec.push_back(vec[i]);

  return vec;
}

std::vector<double> str2vec(std::string magnitudes, double unidad){
  std::istringstream data(magnitudes);
  std::vector<double> data_array;
  double num;
  while(!data.eof()){
    data >> num;
    data_array.push_back(num*unidad);
  }
  return data_array;
}

unsigned int stepEstimator(unsigned int Niter, unsigned int tau, unsigned int min_periods){
  unsigned int periods = Niter/tau;
  if (periods > min_periods)
    return periods*tau;
  else
    return min_periods*tau;
}

std::vector<double> wavearray(double amplitude, unsigned int tau, unsigned int length, double phase){
  std::vector<double> wave;
  wave.resize(length);
  unsigned int wavetop = (tau>=length)? length : tau;
    for(unsigned int i=0; i<wavetop; i++)
      wave[i]=amplitude*cos(_2pi*(i-phase)/tau);

  unsigned int periods = length/tau;
  for(unsigned int i=1; i<periods;i++){
    for(unsigned int j = 0; j<tau ;j++)
      wave[i*tau+j]=wave[j];
  }
  return wave;
}

double simpson_int(const double f_array[], const std::vector<double>& weight){
  unsigned int length=weight.size();

  double Integral = f_array[0]*weight[0];

  for(unsigned int i=1; i<length-1; i+=2)
    Integral+=4*f_array[i]*weight[i];

  for(unsigned int i=2; i<length; i+=2)
    Integral+=2*f_array[i]*weight[i];

  length--;
  Integral+=f_array[length]*weight[length];

  return Integral/3;
}

//Calcular la desviación estandar de una matriz
double stan_dev(std::vector< std::vector<double> > M){
  unsigned int celdas;
  celdas = M.size() * M[0].size();
  double * Aij = new double [celdas];
  for(unsigned int i = 0 ; i<M.size(); i++){
    for(unsigned int j = 0; j<M[0].size(); j++)
      Aij[i*M[0].size() + j] = M[i][j];
  }
  double sd = gsl_stats_sd (Aij, 1, celdas);
  delete[] Aij;
  return sd;
}

bool needSimulation(std::string id_proc, unsigned int size)
{
  struct stat file;
  id_proc = "log_pol_"+id_proc+".dat";

  if (stat(id_proc.c_str(), &file) == -1)
    return true;

  if (file.st_size != size){
    std::remove(id_proc.c_str());
    return true;
  }

  return false;
}
