//============================================================================
// Name        : ising2D.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <vector>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include "IsingModel.h"

#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>

#include <fstream>
#include <sstream>

using namespace std;

//Utilities
vec_sz idx(const vec_sz i, const vec_sz j,vec_sz N)
{//Return the correct index for 1d representation of 2d matrices; 
 //See http://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
  return  i*N+j;
}

//Fitting functions
double quadratic(double *x, double *par)
{// f(x)=par[0]+par[1]*x+par[2]*x*x

  double result=par[0]+par[1]* (*x) +par[2]*(*x)*(*x);
  return result;
}
double power_law(double *x,double *par){
  double Tc=2.26;
  double result=par[0]+pow(fabs(*x-Tc),(-1.)*par[1]);
  return result;
}
double line(double *x,double *par){
  return par[0]*(x[0])+par[1];
}

int magnvstime(IsingModel& ising_model, unsigned int max_mcs=2000, double beta=8, const char *outfilename1="magnvstime.root"){
  //this function draws a graph of the magnetization as a function of MCS
  unsigned int max_side_dim=ising_model.getSideDimension(); //declare the size of the lattice, via its side length
  unsigned int lattice_dim= max_side_dim*max_side_dim; //coincide with the size of a mcs
  vector<double> magnetization(max_mcs);
  cout<<"[*]Plot magnetization as a function of time (MCS)"<<endl;
  cout<<"[+]Current parameters: inverse temperature="<<beta<<" Lattice side dimension="<<max_side_dim<<endl;
  //simulate the convergence
  for (unsigned i=0; i <max_mcs; i++)
    {	if(i%100==0)cout<<"[+]MonteCarlo steps performed:"<<(i)<<"/"<<max_mcs<<"..."<<endl;
      ising_model.simulate(beta,lattice_dim); //1 mcs
      magnetization[i]=ising_model.getMagnetization();
    }
  cout<<"[+]MonteCarlo steps performed:"<<max_mcs<<"/"<<max_mcs<<endl;
  cout<<endl;
  cout<<"[+]Creating histogram"<<endl;
   TH1D* magn_vs_time= new TH1D ("MagnvsTime","Magnetization vs MCS",max_mcs,0,max_mcs);
  magn_vs_time->SetDirectory(0); //because we want a persistent output
  for (unsigned i=0; i<max_mcs; i++)
    {
      magn_vs_time->Fill(i, magnetization[i]);
    }
  TFile out_root(outfilename1,"recreate", "Magnetization vs MCS");
  magn_vs_time->Write(); //write histogram to file
  out_root.Close();
  cout<<"[+]The histogram has been saved to a file."<<endl;
  
	ostringstream matrixbuffer;
  unsigned side=ising_model.getSideDimension();
  cout<<side<<endl;
  for(unsigned k=0; k<side*side; k++){
    matrixbuffer<<ising_model.getSpin(k)<<" ";
    if( (k+1)%side==0 ) matrixbuffer<<"\n";
  }
  
  ofstream outmat;
  outmat.open("matrix.dat");
  cout<<"[+]Saving the spin matrix to a file"<<endl;
  outmat<<matrixbuffer.str();
  outmat.close();
  cout<<endl;
cout<<endl;
  return 0;
}

int critical_temperature(IsingModel &ising_model,unsigned max_mcs=1000){
  //This function computes the critical temperature using the Binder cumulant technique and output result on a file
  vector<unsigned> length_list={12,15,16,20};
  vector<const char *> name_list={"12x12","15x15","16x16","20x20"};
  vec_sz list_size = length_list.size();
  vector<const char *> name_hist_list={"hist10","hist12","hist15","hist16","hist20"};
  const char *hfitName;
  cout<<"[*]Computing critical temperature."<<endl;
  cout<<"Simulation performed for the following lattice sizes:"<<endl;
  for (unsigned i=0; i<list_size; i++) cout<<" "<<length_list[i]<<" ";
  cout<<endl;
  double max_beta = 0.5; // min temp T=1
  double min_beta = 0.4; 	//max temp = infinity
  unsigned temp_steps=100;
  vector<double> inv_temperature(temp_steps);
  for (unsigned i = 0; i<temp_steps; i++ ){
    inv_temperature[i]=min_beta+i*(max_beta-min_beta)/temp_steps;
  }
  //we use a 2x2 matrix representation in which binder_cumulant[index(i,j,N)]
  //represents the value of the bc at temp inv_temperature[j] for the i-th lattice in lenght_list
  //and N is the number of elements in a row
  vector<double> binder_cumulants(list_size*temp_steps,0);
  vec_sz bc_vec_sz = binder_cumulants.size();
  vector<double> magnetization(max_mcs,0);
  for(unsigned i=0; i< list_size; i++)
    { 	cout<<"[+]Simulation on the "<<length_list[i]<<" lattice started..."<<endl;
      ising_model.newGraph(length_list[i]);
      ofstream outf;
      string name(name_hist_list[i]);
      outf.open(name.append(".dat"));
      if (!outf.is_open()){
	cout<<"\n[!]Could not open the outfile"<<endl;
	return -1;
      }
      ostringstream buffer;
      for(unsigned j=0; j < temp_steps ; j++)
	{ //vectors of energy and magnetizations of the current model
	  ising_model.resetGraph(); //reset each time!
	  //burn-in time
	  ising_model.simulate(inv_temperature[j],1000); //burn-in time (empirical)
	  for(unsigned k=0; k<max_mcs; k++)
	    {	ising_model.simulate(inv_temperature[j],length_list[i]*length_list[i]);
	      magnetization[k]=ising_model.getMagnetization();
	    }
	  //procedure to compute binder's cumulant and heat capacity
	  double m_2=0;
	  double m_4=0;
	  double m_2k=0;
	  double m_m=0;
	  for (unsigned k=0; k<max_mcs; k++)
	    {
	      m_2k= magnetization[k]*magnetization[k];
	      m_2+= m_2k;
	      m_4+= m_2k*m_2k;
	      m_m+=magnetization[k];
	    }
	  //mean magnetization 
	  m_m=m_m/max_mcs;
	  //second moment of magnetization 
	  m_2=(m_2 / max_mcs);
	  //fourth moment of the magnetization
	  m_4=(m_4 / max_mcs);
	  binder_cumulants[idx(i,j, temp_steps)] = 1 - ( m_4/(3*m_2*m_2)); //compute the binder's cumulant for this T
	  buffer<<binder_cumulants[idx(i,j, temp_steps)]<<endl;
	}
      cout<<"[+]Simulation completed."<<endl;
      outf<<buffer.str();
      outf.close();
    }

  cout<<"[+]Beginning procedure used to find the critical temperature:"<<endl;
  cout<<"Quadratic fit of the Binder's cumulant obtained for beta in the range("<<min_beta<<","<<max_beta<<")"<<endl;
  unsigned n_fitparam=3;
  unsigned M=list_size*n_fitparam;
  vector<double> fitparam_matrix(M);
  TF1 *fit_f = new TF1("fit_f",quadratic,min_beta,max_beta,n_fitparam);

  double step_width=(max_beta-min_beta)/((double)temp_steps);
  for (unsigned i=0; i<list_size; i++){
    TH1D *hfit = new TH1D(name_list[i],name_list[i],temp_steps,min_beta-step_width/2.,max_beta+step_width/2.);
    for(unsigned j=0; j<temp_steps; j++){
      hfit->Fill(inv_temperature[j],binder_cumulants[idx(i,j,temp_steps)]);
    }    hfit->Fit("fit_f","0");
    TF1 *f_i_fit = hfit->GetFunction("fit_f");
    for(unsigned j=0; j<n_fitparam; j++) {
      fitparam_matrix[idx(i,j,n_fitparam)]=f_i_fit->GetParameter(j);
    }
    cout<<endl;
    cout<<"    Fit "<<i+1<<"of "<<list_size<<" completed!"<<endl;
  }
  //once we have a fit we can approximate the  T_C finding the root of fit_i-fit_j;
  //since this should converge with L-> infinity, we compute only fit_i-fit_{i+1}
  //we use TF1::GetX(0,0,5) to find roots in the interval
  vector<double> criticalT(list_size-1);
  for(unsigned i=0; i<list_size-1 ; i++)
    { double params[3] = {0,0,0};
      for (unsigned j=0; j<n_fitparam; j++){
	  params[j]=fitparam_matrix[idx(i,j,M)]-fitparam_matrix[idx(i+1,j,M)];
	}
      fit_f->SetParameters(params);
      criticalT[i]=1./(fit_f->GetX(0,0.44,0.5) ); //revert inverse temp
    }
  cout<<"[+]Computing critical temperature T_c using the Binder's cumulant; "<<endl;
  cout<<"   The following values are the estimates for T_c :"<<endl;
  for(unsigned i=0; i< list_size-1; i++) cout<<"    T_c(L="<<length_list[i]<<")="<<criticalT[i]<<endl;
  //SAVE TO TFILE
  cout<<endl;
  cout<<endl;

  return 0;
}

int magnvstemp(IsingModel& ising_model, unsigned long int max_mcs=6000, unsigned side_dim=16, const char *outfilename="magnvsT.root")
{
  //CREATE A  GRAPH OF ABSOLUTE MAGNETIZATION VS TEMPERATURE
  //declare streams for the output
  cout<<"[*]Magnetization as a function of temperature"<<endl;
  cout<<"[+]Parameters: MCS="<<max_mcs<<" lattice side="<<side_dim<<endl;
  ofstream outf;
  outf.open("magnvsT.dat");
  if (!outf.is_open()){
    cout<<"[!]Could not open the outfile"<<endl;
    return -1;
  }
  ostringstream buffer;
  //parameters of the simulations
  double max_temp = 4;
  double min_temp = 0.1; 	
  unsigned temp_steps=100;
  vector<double> inv_temperature(temp_steps,0);
  for (unsigned i = 0; i<temp_steps; i++ ){
    inv_temperature[i]=1/(min_temp+i*(max_temp-min_temp)/temp_steps );
  }
  //variabl for the results
  double m_m;//average magnetization
  //histogram
  TH1D* magn_vs_temp= new TH1D ("MagnvsT","Magnetization vs T",100*temp_steps,min_temp-0.001,max_temp+0.001);
  magn_vs_temp->SetDirectory(0); //because we want a persistent output
  cout<<"[+]Progress:"<<endl;
  //set the model in an ordered phase
  ising_model.resetGraph();
  for (unsigned j = 0; j<temp_steps; j++ ){
    
    m_m=0;
    if(j%10==0) cout<<"  Step "<<j<<" of "<<temp_steps<<endl;
    ising_model.resetGraph(); //reset each time (in an ordered state)
    ising_model.simulate(inv_temperature[j],1000);//burn-in time
    for(unsigned k=0; k<max_mcs; k++)
      {
	ising_model.simulate(inv_temperature[j],side_dim*side_dim);
	m_m+=ising_model.getMagnetization()/max_mcs;
      }
    buffer<<1./inv_temperature[j]<<"\t"<<fabs(m_m)<<"\n";
    magn_vs_temp->Fill(1./inv_temperature[j], fabs(m_m));
  }
  outf<<buffer.str();
  outf.close();
  TFile out_root(outfilename,"recreate", "Magnetization vs T");
  magn_vs_temp->Write(); //write histogram to file
  out_root.Close();
  cout<<"[+]The histogram has been saved to a file."<<endl;
  cout<<endl;
  cout<<endl;

  return 0;
}


int critical_exponents(IsingModel& ising_model, unsigned long int max_mcs =1000, unsigned  side_dim=16, const char * outfilename1="critical_exponents.root"){
  cout<<"[*]Computing critical exponents"<<endl;
  cout<<"[+]Parameters: MCS="<<max_mcs<<" lattice side="<<side_dim<<endl;
  ising_model.newGraph(side_dim);
  
  //parameters of the simulations
  double max_temp = 4;
  double min_temp = 0.1; 	
  unsigned temp_steps=150;
  vector<double> inv_temperature(temp_steps,0);
  vector<double> heat_capacity(temp_steps,0);
  vector<double> susceptibility(temp_steps,0);
  for (unsigned i = 0; i<temp_steps; i++ ){
    inv_temperature[i]=1/(min_temp+i*(max_temp-min_temp)/temp_steps );
  }

  //variables for the results
  vector<double> energy(max_mcs,0);
  vector<double> magnetization(max_mcs,0);
  cout<<"[+]Progress:"<<endl;
  for (unsigned j = 0; j<temp_steps; j++ ){
    double m_m=0;
    double m_2=0;
    double e_2=0;
    double e_m=0;
   if(j%10==0) cout<<"   Step "<<j<<" of "<<temp_steps<<endl; 
    ising_model.resetGraph(); //reset each time!
    ising_model.simulate(inv_temperature[j],300);
    for(unsigned k=0; k<max_mcs; k++)
      {
	ising_model.simulate(inv_temperature[j],side_dim*side_dim);
	energy[k]=ising_model.getEnergy()/(side_dim*side_dim);
	magnetization[k]=ising_model.getMagnetization();
	m_m+=magnetization[k];
	m_2+=magnetization[k]*magnetization[k];
	e_m+= energy[k]/max_mcs;
	e_2+= energy[k]*energy[k]/max_mcs;
      }
    m_m=m_m/max_mcs;
    m_2=m_2/max_mcs-m_m;
    e_m=e_m;
    e_2=e_2-e_m*e_m;
    heat_capacity[j]=inv_temperature[j]*inv_temperature[j]*(e_2);
    susceptibility[j]=m_2/inv_temperature[j];
  }
  cout<<"[+]Writing the heat capacity and susceptibility vectors"<<endl;
  ofstream outfH;
  outfH.open("heatcap.dat");
  if (!outfH.is_open()){
    cout<<"\n[!]Could not open the outfile"<<endl;
    return -1;
  }
  ostringstream bufferH;

  ofstream outf_chi;
  outf_chi.open("chivsT.dat");
  if (!outf_chi.is_open()){
    cout<<"\n[!]Could not open the outfile"<<endl;
    return -1;
  }
  ostringstream bufferC;
  for(unsigned k=0; k<temp_steps; k++){
    bufferH<<1./inv_temperature[k]<<"\t"<<heat_capacity[k]<<"\n";
    bufferC<<1./inv_temperature[k]<<"\t"<<susceptibility[k]<<"\n";
  }
  outfH<<bufferH.str();
  outfH.close();
  outf_chi<<bufferC.str();
  outf_chi.close();

  //Fit
  TF1 *fit_pl = new TF1("fit_pl",power_law,min_temp,max_temp,2);
  TH1D *heatcap_hist = new TH1D("heat_capacity_hist","Heat capacity vs T",temp_steps,min_temp-0.001,max_temp+0.001);
    for(unsigned j=0; j<temp_steps; j++){
      heatcap_hist->Fill(1./inv_temperature[j],heat_capacity[j]);
    }
    heatcap_hist->Fit("fit_pl","0");
    TF1 *fit_hc = heatcap_hist->GetFunction("fit_pl");
    TFile *f_hc=new TFile("critical_exponents.root","recreate");
    heatcap_hist->Write();
    f_hc->Close();
    cout<<"[+]Fit results for heat capacity"<<endl;
    cout<<"Offset (error):"<<fit_hc->GetParameter(0)<<" ("<<fit_hc->GetParError(0) <<")"<<endl;
    cout<<"Exponent (error):"<<fit_hc->GetParameter(1)<<" ("<<fit_hc->GetParError(1)<<")"<<endl;
    
    cout<<endl;
    cout<<endl;
    return 0;
}


int simulate_ising( unsigned mcs=3000,double beta=10, unsigned  max_side_dim=20)
{ 	
  cout<<"[***]Metropolis Monte Carlo simulation of a 2D Ising model."<<endl;
  cout<<"[+]Creating the Ising model object"<<endl;
  IsingModel ising_model(max_side_dim);
  cout<<endl;
  //Magnetization as a function of MC steps performed
  magnvstime(ising_model,mcs,beta);
  //Magnetization as a function of temperature
  magnvstemp(ising_model);
  //Compute the critical temperature
  critical_temperature(ising_model,mcs);
  //Compute the critical exponents of susceptibility and heat capacity
  critical_exponents(ising_model,10000);
  //Drawings
  //1) Magnetization/ MCS
  new TCanvas();
  TFile *mvst_file=new TFile("magnvstime.root","read");
  TH1D *mvst_hist=(TH1D*) mvst_file->Get("MagnvsTime");
  mvst_hist->SetDirectory(0);
  mvst_hist->Draw("hist");
  //2) Magnetization / T
  new TCanvas();
  TFile *mvsTemp_file=new TFile("magnvsT.root","read");
  TH1D *mvsTemp_hist=(TH1D*) mvsTemp_file->Get("MagnvsT");
  mvsTemp_hist->SetDirectory(0);
  mvsTemp_hist->SetMarkerStyle(7);
  mvsTemp_hist->Draw("PHIST");
  return 0;
}

int estimate_errors(long unsigned side_dim=20, double beta_err=0.5){
  //This is an auxiliary function used for testing the limitation of our statistical analysis
  cout<<"[*]Raw estimate of the Monte Carlo error on magnetization and energy measurements:"<<endl;
  cout<<"[...]the running autocorrelation time is computed at T=2, for 10000 MCS on a 20*20 lattice; the corresponding error is used for all the points."<<endl;
  IsingModel ising_model(side_dim);
  long unsigned max_mcs=10000;

  //variables for the results
  vector<double> energy(max_mcs,0);
  vector<double> magnetization(max_mcs,0);
  double m_m=0;
  double m_2=0;
  double e_2=0;
  double e_m=0;
  double m_acorr_integral=0;
  double m_acorr_time=0;
  double e_acorr_integral=0;
  double e_acorr_time=0;
      
  for(unsigned k=0; k<max_mcs; k++)
    {
      ising_model.simulate(beta_err,side_dim*side_dim);
      energy[k]=ising_model.getEnergy()/(side_dim*side_dim);
      magnetization[k]=ising_model.getMagnetization();
      m_m+=magnetization[k];
      m_2+=magnetization[k]*magnetization[k]/max_mcs;
      e_m+= energy[k]/max_mcs;
      e_2+= energy[k]*energy[k]/max_mcs;
      if(k>0){
	m_acorr_integral+=magnetization[k-1]*magnetization[k]/(max_mcs-1);
	e_acorr_integral+=energy[k-1]*energy[k]/(max_mcs-1);
      }
    }
  m_m=m_m/max_mcs;
  m_2=m_2-m_m*m_m;
  e_m=e_m;
  e_2=e_2-e_m*e_m;
  m_acorr_time=(m_acorr_integral - m_m*m_m)/m_2;
  e_acorr_time=(e_acorr_integral - e_m*e_m)/e_2;
  double m_err=sqrt(m_2/max_mcs*(1+2*m_acorr_time) );
  double e_err=sqrt(e_2/max_mcs*(1+2*e_acorr_time) );
  cout<<"[+]Magn, err"<<m_m<<"\t"<<m_err<<endl;
  cout<<"[+]En, err"<<e_m<<"\t"<<e_err<<endl;
  
  return 0;
}
