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

using namespace std;

//Utilities
vec_sz idx(const vec_sz i, const vec_sz j,vec_sz N)
{	//Return the correct index for 1d representation of 2d matrices; do not change data members
	//See http://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
	return  i*N+j;
}

double quadratic(double *x, double *par)
   {// f(x)=par[0]+par[1]*x+par[2]*x*x

      double result=par[0]+par[1]* (*x) +par[2]*(*x)*(*x);
      return result;
   }

int calibrate_ising(double beta, unsigned int max_mcs, const char * outfilename1){

	unsigned int max_side_dim=20; //declare the size of the lattice, via its side length
	unsigned int lattice_dim= max_side_dim*max_side_dim; //coincide with the size of a mcs
	cout<<"[+]Metropolis Monte Carlo simulation of the 2D Ising model."<<endl;
	cout<<"[+]Creating the run manager for the simulation..."<<endl;
	IsingModel ising_model(max_side_dim );  //initial size represents the maximum size allowed


	TH1D* magn_vs_time= new TH1D ("histCalibration","Magnetization vs MCS",max_mcs,0,max_mcs);
	magn_vs_time->SetDirectory(0); //because we want the output still shown when the file will be closed

	//Convergence time (MCS) computed (UPPER BOUND) near theoretic T_C // see Landau, Binder for motivations
	cout<<"[+]Calibration procedure for a L*L lattice started."<<endl;
	cout<<"   It will compute the convergence time of the Metropolis algorithm for a lattice with maximum L="<<max_side_dim<<endl;
	cout<<"   Using a greater L in the simulation can result in unreliable measurements."<<endl;

	vector<double> magnetization(max_mcs);

	cout<<"[+]Current beta="<<beta<<endl;
	//simulate the convergence
	for (unsigned i=0; i <max_mcs; i++)
	{	if(i%100==0)cout<<"[+]MonteCarlo steps performed:"<<(i)<<"/"<<max_mcs<<"..."<<endl;
		ising_model.simulate(beta,lattice_dim); //1 mcs
		magnetization[i]=ising_model.getMagnetization();
	}


	cout<<"[+]Creating histogram"<<endl;
	new TCanvas() ;
	cout<<"[+]Populating histogram"<<endl;
	for (unsigned i=0; i<max_mcs; i++)
	{
		magn_vs_time->Fill(i, magnetization[i]);
	}
	cout<<"[+]Drawing the histogram..."<<endl;
	new TCanvas();
	magn_vs_time->DrawCopy("hist");
	cout<<"Recreating the file "<<outfilename1<<"..."<<endl;
	TFile out_root(outfilename1,"recreate", "Ising simulation results");
	cout<<"[+]Saving..."<<endl;
	magn_vs_time->Write(); //write histogram to file
	cout<<"[+]Closing the file"<<endl;
	out_root.Close();
	cout<<"[+]The histogram has been saved to a file."<<endl;

	return 1;
}

int simulate_ising(const char * outfilename1, unsigned long int max_mcs =0, unsigned  max_side_dim=20)
{	//User program; writes results to a file resp. in ROOT / raw numerical format

	//Ising model:
	//Create a IsingModel object (inherits the MCMC simulation routine from SimulationModel)
	//Initialize the model at random at any temperature
	//Use MCMC to evaluate convergence toward equilibrium, save the results (in a vector, discrete times)
	//Use the Binder cumulant to identify the critical temperature :
	//To do this, one should vary the size of the lattice and compute the cumulant at various temperatures
	//Then, an interpolation is needed to find the fixed point.
	//When the critical point is found, the critical exponents can be extracted by a fit and compared
	//with the known ones

	//Sheared Ising:
	//Compute the distributions and verify FDT
	/////////////////////////////////////////////////////////////////

	//1. Ising model simulation
	//TODO: add timing!

	cout<<"[+]Metropolis Monte Carlo simulation of a 2D Ising model."<<endl;
	cout<<"[+]Creating the run manager for the simulation..."<<endl;
	IsingModel ising_model(max_side_dim);
	cout<<"[+]Creation completed"<<endl;
	//The aim of the program is to evaluate the
	//T_C; we will compute an array of Binder's cumulant at various temperatures in (0.2,infinity)
	//for a choice of lattice l=(8,16,32,64)
	//then we will perform a quadratic fit using ROOT (see e.g. https://root.cern.ch/root/HowtoFit.html)
	//and take T_C equal to the fixed point

	cout<<"[+]Starting the main simulation procedure."<<endl;
	//Lattice side size should always be tested against max_side_dim
	vector<unsigned> length_list={8}; //,12,14,16,20};
	//vector<const char *> name_list={"L=8","L=12","L=14","L=16","L=20"};
	vec_sz list_size = length_list.size();
	cout<<"[+]Simulating 2D Ising model on a square lattice with L*L spins "<<endl;
	cout<<"   for the following values of L:"<<endl;
	for (unsigned i=0; i<list_size; i++) cout<<" "<<length_list[i]<<" ";
	cout<<endl;
	for (unsigned i =0; i<list_size; i++)
			if(length_list[i]>max_side_dim)
				cout<<"[!]Maximum lattice size ("<<max_side_dim<<") exceeded: results could be unreliable!"<<endl;
	double max_beta = 20.; // min temp T=5e-2
	double min_beta = 0.; 	//max temp = infinity
	unsigned temp_steps=200;
	vector<double> inv_temperature(temp_steps);
	for (unsigned i = 0; i<temp_steps; i++ ){
		inv_temperature[i]=min_beta+i*(max_beta-min_beta)/temp_steps;
	}
	//we used a 2x2 matrix representation in which binder_cumulant[index(i,j,N)]
	//represents the value of the bc at temp inv_temperature[j] for the i-th lattice in lenght_list
	//and N is the number of elements in a row
	vector<double> binder_cumulants(list_size*temp_steps,0);
	vec_sz bc_vec_sz = binder_cumulants.size();
	//same for the energy fluctuations vector
	vector<double> heat_capacity(list_size*temp_steps,0);
	vec_sz hc_vec_sz = heat_capacity.size();
	cout<<"[+]Parameters correctly initialized. "<<endl;

	vector<double> energy(max_mcs,0);
	vector<double> magnetization(max_mcs,0);
	for(unsigned i=0; i< list_size; i++)
	{ 	cout<<"[+]Simulation for L="<<length_list[i]<<" started..."<<endl;
		ising_model.newGraph(length_list[i]);
		cout<<"[+]New graph created"<<endl;
		for(unsigned j=0; j < temp_steps ; j++)
		{	//mcs_i = max_mcs / (max_side_dim/length_list[i]); //integer division, this linear speed up "should not mess up"...
			//unsigned long int max_steps_i = mcs_i*length_list[i]*length_list[i]; //total n of steps
			//vectors of energy and magnetizations of the current model

			//burn-in time
			ising_model.simulate(inv_temperature[j],150); //burn in time (empirical)
			for(unsigned k=0; k<max_mcs; k++)
			{
				ising_model.simulate(inv_temperature[j],length_list[i]*length_list[i]);
				energy[k]=ising_model.getEnergy();
				magnetization[k]=ising_model.getMagnetization();
				if(magnetization[k]<-1) return -3;
			}
			//procedure to compute binder's cumulant and heat capacity
			double m_2=0;
			double m_4=0;
			double m_2k=0;
			double e_2=0;
			double e_m=0;
			for (unsigned k=0; k<max_mcs; k++)
			{
				m_2k= magnetization[k]*magnetization[k];
				m_2+= m_2k;
				m_4+= m_2k*m_2k;
				e_2+= energy[k]*energy[k];
				e_m+= energy[k];
			}
			m_4=m_4 / max_mcs;
			m_2=m_2 / max_mcs;
			e_2=e_2 / max_mcs;
			e_m=e_m / max_mcs;
			cout<<"[D] m_4 "<<m_4<<" "<<"m_2 "<<m_2<<" "<<" e_2 "<<e_2<<" e_m "<<e_m<<endl;
			binder_cumulants[idx(i,j, temp_steps)] = 1 - ( m_4/(3*m_2*m_2)); //compute the binder's cumulant for this T
			heat_capacity[idx(i,j,temp_steps)] =inv_temperature[j]*inv_temperature[j]*(e_2-(e_m*e_m)); // compute the heat capacity(T)=k_b T^2 (var E)
			ising_model.resetGraph(); //reset each time!
		}
		cout<<"   [+]...completed."<<endl;
	}

	cout<<"[+]Beginning procedure used to find the critical temperature:"<<endl;
	cout<<"   quadratic fit of the Binder's cumulant obtained for beta in the range("<<min_beta<<","<<max_beta<<")"<<endl;
	unsigned n_fitparam=3;
	unsigned M=list_size*n_fitparam;
	vector<double> fitparam_matrix(M);
	//cout<<"[D]fitparam_matrix created"<<endl;

	TF1 *fit_f = new TF1("fit_f",quadratic,0.,5.,n_fitparam); //range(0,5), 3 parameters
	cout<<"[D]fit_f TF1 created"<<endl;


	//TODO:fitting procedure for each lattice, update fitparam_matrix
	for (unsigned i=0; i<list_size; i++){
		TH1D *hfit = new TH1D("temp_hist","Binder's cumulant vs Inverse temperature ",temp_steps,min_beta,max_beta);
		cout<<"[D]	Temp hist created"<<endl;
		//populate the hist
		for(unsigned j=0; j<temp_steps; j++){
			//cout<<"[D]		Current binder's cumulant is:"<<binder_cumulants[idx(i,j,temp_steps)]<<endl;
			hfit->Fill(inv_temperature[j],binder_cumulants[idx(i,j,temp_steps)]);
			//cout<<"[D]		hfit populated!"<<endl;
		}
		cout<<"[D]Accessing hfit and fitting fit_f to it"<<endl;
		hfit->Fit("fit_f");
		cout<<"[D] Fit done"<<endl;
		new TCanvas();
		hfit->DrawCopy(); //DEBUG
		//cout<<"Printed"<<endl;
		TF1 *f_i_fit = hfit->GetFunction("fit_f");
		cout<<"[D] Temp fit function obtained from histogram:"<<endl;
		for(unsigned j=0; j<n_fitparam; j++) {
			fitparam_matrix[idx(i,j,n_fitparam)]=f_i_fit->GetParameter(j);
		}
		cout<<endl;
		cout<<"    Fit "<<i+1<<"/"<<list_size<<" completed!"<<endl;
		hfit->SetDirectory(0);

	}

	//once we have a fit we can approximate the  T_C finding the root of fit_i-fit_j;
	//since this should converge with L-> infinity, we compute only fit_i-fit_4
	//we use TF1::GetX(0,0,5) to find roots in the interval
	vector<double> criticalT(list_size-1);
	for(unsigned i=0; i<list_size-1 ; i++)
	{
		double params[3] = {0,0,0};
		for (unsigned j=0; j<n_fitparam; j++)
		{
			//cout<<"[D] param ("<<i<<","<<j<<") = "<<fitparam_matrix[idx(i,j,M)]<<" ";
			//cout<<"[D] param ("<<i<<","<<4<<") = "<<fitparam_matrix[idx(list_size-1,j,M)]<<endl;
			params[j]=fitparam_matrix[idx(i,j,M)]-fitparam_matrix[idx(list_size-1,j,M)];
		}
		fit_f->SetParameters(params);
		criticalT[i]=1./(fit_f->GetX(0,0,5) ); //revert inverse temp
	}
	cout<<"[+]Computing critical temperature T_c using the Binder's cumulant; "<<endl;
	cout<<"   The following values are the estimates for T_c (in growing order of reliability)"<<endl;
	cout<<"   The reference value is the fit of the Binder's cumulant for L="<<length_list[list_size-1]<<endl;
	for(unsigned i=0; i< list_size-1; i++) cout<<"    T_c(L="<<length_list[i]<<")="<<criticalT[i]<<endl;

	//TODO:
	//Now we compute the critical exponent for the heat capacity C \propto (T_C -T)^alpha
	//Heat capacity is computed using the FR theorem for equilibrium inside the SimulationModel class
	//This is done for all the lattices


	//TODO: raw output file close
	return 0;

}




