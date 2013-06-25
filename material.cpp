#include "material.h"

/*Constructor:
Dimensiona y encera a los vectores del sistema. Llena sus datos iniciales */
Material::Material(unsigned int L,double p, std::string ID,
		 bool polarizar, bool log_H, bool log_S){
  rng = gsl_rng_alloc (gsl_rng_taus);// Start random number generator
  ExpID=ID;
  rho = p;
  logH=log_H;
  logS=log_S;

  // Setup Material
  PNR = L*L*L;
  sigma.resize(PNR);
  mu_E.resize(PNR);

  /*Setup Topological configuration of
    PNRs inside the material */
  G.resize(PNR);
  set_space_config(L);
  J.resize(PNR);
  set_interaction_dipole_config(polarizar);
}

/*Destructor:
Frees Memory*/
Material::~Material(){
  gsl_rng_free (rng);
  for(unsigned int i=0; i<PNR; i++){
    J[i].clear();
    G[i].clear();
  }
  J.clear();
  G.clear();
  mu_E.clear();
  sigma.clear();
}

void Material::set_space_config(unsigned int L){

  /* Generates a simple cubic lattice where each
   * PNRs is assigned to a lattice point,
   * this is only a topological consideration, not
   * a real spatial configuration */

  for(unsigned int i = 0; i < PNR; i++)
    G[i].resize(6);

  unsigned int ind_xy, L2=L*L;
  std::vector< std::vector<unsigned int> > R;
  R.resize(PNR);

  for(unsigned int i=0; i<R.size(); i++){
    // Coeficientes vector posición i-ésima PNR
    ind_xy = i % L2;
    R[i].resize(3);
    R[i][0] = ind_xy % L;
    R[i][1] = ind_xy / L;
    R[i][2] = i / L2;

    /*Encontrar índices de los primeros vecinos.
     *solo existen 6: arriba y abajo(+z, -z), derecha e izquierda(+y, -y), adelante y atraz(+x, -x).
     *También debo aplicar las condiciones de borde en este caso*/
    G[i][0] = (R[i][2] == L-1 )	?i - (L-1)*L2	:i + L2;//arriba
    G[i][1] = (R[i][2] == 0 )	?i + (L-1)*L2	:i - L2;//abajo
    G[i][2] = (R[i][1] == L-1 )	?i - (L-1)*L	:i + L;//derecha
    G[i][3] = (R[i][1] == 0)	?i + (L-1)*L	:i - L;//izquierda
    G[i][4] = (R[i][0] == L-1 )	?i - L+1	:i + 1;//adelante
    G[i][5] = (R[i][0] == 0 )	?i + L-1	:i - 1;//atraz
  }

  //Free memory
  for(unsigned int i=0; i<R.size();i++)
    R[i].clear();
  R.clear();
}

void Material::Jex(){
  // Assigns J with same structure as G
  for(unsigned int i = 0; i < PNR; i++)
    J[i].assign (G[i].size(),-1000);

  // Look for the interaction of the j neighbour with PNR i
  for(unsigned int i=0;i<PNR;i++){
    for(unsigned int j=0; j<J[i].size() ;j++){
      if ( J[i][j] == -1000 ){ /* if not assigned
	* look if the neighbour has already an assigned value*/
        for(unsigned int k=0; k<G[ G[i][j] ].size() ;k++){
	  if (G[ G[i][j] ][k] == i){ // find matching neighbour
	    /* if neighbour hasn't an interaction energy assigned
	     * yet, assign it, then copy it */
	    if ( J[ G[i][j] ][k] == -1000)
		 J[ G[i][j] ][k] = gsl_ran_gaussian(rng,1)+rho;
	    J[i][j] = J[ G[i][j] ][k];
	    k=G[ G[i][j] ].size();//end loop
	  }
	}}}}
}

void Material::set_sigma(bool polarize){
  if (polarize)
    sigma.assign(PNR,1);
  else{
    for(unsigned int i=0; i<PNR; i++)
      sigma[i] = (gsl_rng_uniform(rng)-0.5 > 0)? 1:-1;
  }
}

void Material::set_mu(bool polarize){
  set_sigma(polarize);
  for(unsigned int i=0; i<PNR; i++)
    mu_E[i]  = gsl_rng_uniform(rng);
//   array_print_bin(mu_E,"mu"+ExpID+".dat");
}

void Material::set_interaction_dipole_config(bool polarizar){
  gsl_rng_set(rng, std::time(NULL) );
  // Calculate new values for exchange Energy an dipolar momentum
  Jex();
  set_mu(polarizar);
}

double Material::total_E(double E){
  double Hamil = 0;
  for(unsigned int i = 0; i < PNR; i++){
    for(unsigned int j = 0; j < 6; j++)
      Hamil -= J[i][j]*sigma[i]*sigma[G[i][j]];
    Hamil -= E*mu_E[i]*sigma[i];
  }
  return Hamil;
}

double Material::delta_E(unsigned int idflip, double E){
  double dHamil = 0;
  for(unsigned int i = 0; i<6; i++)
    dHamil += J[idflip][i]*sigma[idflip]*sigma[G[idflip][i]];
  dHamil *=4;
  dHamil += 2*E*mu_E[idflip]*sigma[idflip];

  return dHamil;
}

double Material::norm_pol(){
  double P=0;
  for(unsigned int i=0; i<PNR; i++)
    P += mu_E[i]*sigma[i];

  return (double) P / PNR;
}

void Material::MonteCarloStep(double T, double E_field){
  for(unsigned int idflip = 0; idflip < PNR; idflip++){
    double dH = delta_E(idflip, E_field);
    if ( dH < 0 || exp(-dH/T) >= gsl_rng_uniform(rng) )
      sigma[idflip] *= -1;
  }
}

void Material::state(double T, std::vector< double >& field, unsigned int Equilibration_Iter, bool measure){
  //vector historial de polarización por experimento
  std::vector<double> log_pol,log_H;
  log_pol.resize(field.size());
  log_H.resize(field.size());  
  std::vector<int> log_sigma (PNR,0);
  
  // Evaluate material's state during given amount of steps
  unsigned int start = (Equilibration_Iter>0) ? field.size()-Equilibration_Iter : 0;
  for(unsigned int i = start ; i< field.size(); i++){
    MonteCarloStep(T,field[i]);
    if (measure){
      log_pol[i] = norm_pol();
      if (logS) update_log_sigma(log_sigma);
      if (logH) log_H[i]=total_E(field[i]);
    }
  }

  // Save recorded data in binary format
  if (measure){
    array_print_bin(log_pol,"log_pol_"+ExpID+".dat");
    if (logS) array_print_bin(log_sigma,"log_sigma_"+ExpID+".dat");
    if (logH) array_print_bin<double>(log_H,"log_H"+ExpID+".dat");
  }
  
  log_H.clear();
  log_pol.clear();
  log_sigma.clear();
}

void Material::update_log_sigma(std::vector< int >& log_sigma){
  for(unsigned int s = 0; s<PNR ; s++)
    log_sigma[s] += sigma[s];
}
std::string Material::repr()
{
  std::ostringstream descrip;
  descrip<<"\n Random: "<<gsl_rng_uniform(rng);
  descrip<<"\n Exp ID:"<<ExpID;
  descrip<<"\n rho:"<<rho;
  descrip<<"\n Log H:"<<logH<<" , S:"<<logS;
  descrip<<"\n #PNRs:"<<PNR;
  descrip<<"\n sigma size:"<<sigma.size();
  descrip<<"\n mu size:"<<mu_E.size();
  return descrip.str();

}
