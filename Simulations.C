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
vec_sz index(const vec_sz x, const vec_sz y,unsigned N)
{	//Return the correct index for 1d representation of 2d matrices; do not change data members
	//See http://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
	return  x + N* y;
}


double quadratic(double *x, double *par)
   {// f(x)=par[0]+par[1]*x+par[2]*x*x
      Double_t arg = 0;
      if (par[3] != 0){
    	  cout<<"[!]A quadratic polinomial needs exactly 3 parameters"<<endl;
    	  return 0;
      }

      double result=par[0]+par[1]* (*x) +par[2]*(*x)*(*x);
      return result;
   }

int simulate_ising( const char * outfilename1, const char * outfilename2, unsigned long int max_mcs = 0)
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

	unsigned int max_side_dim=20; //declare the size of the lattice, via its side length
	unsigned int lattice_dim= max_side_dim*max_side_dim; //coincide with the size of a mcs
	cout<<"[+]Metropolis Monte Carlo simulation of a 2D Ising model."<<endl;
	cout<<"[+]Creating the run manager for the simulation..."<<endl;
	IsingModel *ising_model = new IsingModel(max_side_dim );  //initial size represents the maximum size allowed
	cout<<"Recreating the file "<<outfilename1<<"..."<<endl;
	TFile out_root(outfilename1,"recreate", "Ising simulation results");

	/////////////////////////////////////////////////////////////
	//START CALIBRATION
	if(max_mcs==0)
	{
		//Convergence time (MCS) computed (UPPER BOUND) near theoretic T_C // see Landau, Binder for motivations
		cout<<"[+]Calibration procedure for a L*L lattice started."<<endl;
		cout<<"    It will compute the convergence time of the Metropolis algorithm for a lattice with maximum L="<<max_side_dim<<endl;
		cout<<"    Using a greater L in the simulation can result in unreliable measurements."<<endl;
		max_mcs= 2500; //new max number of mcs
		vector<double> magnetization(max_mcs);
		double beta_c = 0.44; // note that this is computed for k_b=1 and J=1, T_C = 2.269
		int reached_convergence=0;
		//simulate the convergence
		for (unsigned i=0; i <max_mcs; i++)
		{	if(i%100==0)cout<<"[+]MonteCarlo steps performed:"<<(i)<<"/"<<max_mcs<<"..."<<endl;
			ising_model->simulate(beta_c,lattice_dim); //1 mcs
			magnetization[i]=ising_model->getMagnetization();
		}
		cout<<"[+]Creating histogram"<<endl;
		new TCanvas() ;
		TH1D *magn_vs_time = new TH1D("histCalibration","Magnetization vs MCS",max_mcs,0,max_mcs); //create an histogram to display
		cout<<"[+]Populating histogram"<<endl;
		for (unsigned i=0; i<max_mcs; i++)
		{
			magn_vs_time->AddBinContent(i, magnetization[i]);
		}
		cout<<"[+]Drawing the histogram..."<<endl;
		new TCanvas();
		magn_vs_time->Draw();
		cout<<"[+]Saving..."<<endl;
		magn_vs_time->Write(); //write histogram to file
		cout<<"[+]The histogram has been saved to a file.";
		cout<<"[!]Insert the maximum number of MonteCarlo step resulting from calibration:"<<endl;
		cin>>max_mcs;
	}

	//The next step of the programme is to evaluate the
	//T_C; we will compute an array of Binder's cumulant at various temperatures in (0.2,infinity)
	//for a choice of lattice l=(8,16,32,64)
	//then we will perform a quadratic fit using ROOT (see e.g. https://root.cern.ch/root/HowtoFit.html)
	//and take T_C equal to the fixed point

	cout<<"[+]Starting the main simulation procedure."<<endl;

	//Lattice side size should always be tested against max_side_dim
	vector<unsigned> length_list={8,10,14,16,20};
	vec_sz list_size = length_list.size();
	cout<<"[+]Simulating 2D Ising model on a square lattice with L*L spins "<<endl;
	cout<<"   for the following values of L:"<<endl;
	for (unsigned i=0; i<list_size; i++) cout<<" "<<length_list[i]<<" ";
	cout<<endl;
	for (unsigned i =0; i<list_size; i++)
			if(length_list[i]>max_side_dim)
				cout<<"[!]Maximum lattice size ("<<max_side_dim<<") exceeded: results could be unreliable!"<<endl;
	double max_beta = 200; // min temp T=5e-3
	double min_beta = 0; 	//max temp = infinity
	double temp_steps=80;
	vector<double> inv_temperature(temp_steps);
	for (unsigned i = 0; i<temp_steps; i++ ){
		inv_temperature[i]=min_beta+i*(max_beta-min_beta)/temp_steps;
	}
	//we use a 2x2 matrix representation in which binder_cumulant[index(i,j)]
	//represents the value of the bc at temp inv_temperature[j] for the i-th lattice in lenght_list
	vector<double> binder_cumulants(list_size*temp_steps);
	vec_sz bc_vec_sz = binder_cumulants.size();
	//same for the energy fluctuations vector
	vector<double> heat_capacity(list_size*temp_steps);
	vec_sz hc_vec_sz = heat_capacity.size();
	cout<<"[+]Run parameters correctly initialized. "<<endl;
	unsigned long int mcs_i = 0;

	for(unsigned i=0; i< list_size; i++)
	{ 	cout<<"[+]Simulation for L="<<length_list[i]<<" started..."<<endl;
		ising_model->newLattice(length_list[i]);
		for(unsigned j=0; j < temp_steps ; j++)
		{	mcs_i = max_mcs / (max_side_dim/length_list[i]); //integer division, this linear speed up "should not mess up"...
			unsigned long int max_steps_i = mcs_i*length_list[i]*length_list[i]; //total n of steps
			//vectors of energy and magnetizations of the current model
			vector<double> energy(max_steps_i);
			vector<double> magnetization(max_steps_i);
			//burn-in time
			ising_model->simulate(inv_temperature[j],(max_steps_i/100)+1); //empirical...
			for(unsigned k=0; k<max_steps_i; k++)
			{
				ising_model->simulate(inv_temperature[j],1); //advance 1 step at the time
				energy[k]=ising_model->getEnergy();
				magnetization[k]=ising_model->getMagnetization();
			}
			//procedure to compute binder's cumulant and heat capacity
			double m_2=0;
			double m_4=0;
			double m_2k=0;
			double e_2=0;
			double e_m=0;
			for (unsigned k=0; k<max_steps_i; k++)
			{
				m_2k=magnetization[k]*magnetization[k];
				m_2+=m_2k  ;
				m_4+= m_2k*m_2k;
				e_2 = energy[k]*energy[k];
				e_m = energy[k];
			}
			m_4=m_4 / max_steps_i;
			m_2=m_2 / max_steps_i;
			e_2=e_2 / max_steps_i;
			e_m=e_m / max_steps_i;
			binder_cumulants[index(i,j, temp_steps)] = 1- (m_4/(m_2*m_2)); //compute the binder's cumulant for this T
			heat_capacity[index(i,j,temp_steps)] =inv_temperature[j]*inv_temperature[j]*(e_2-(e_m*e_m)); // compute the heat capacity(T)=k_b T^2 (var E)
			ising_model->resetLattice(); //reset each time!
		}
		cout<<"   Completed."<<endl;
	}

	cout<<"[+]Beginning procedure used to find the critical temperature:"<<endl;
	cout<<"   quadratic fit of the Binder's cumulant obtained for beta in the range("<<min_beta<<max_beta<<")"<<endl;
	unsigned n_fitparam=3;
	unsigned M=list_size*n_fitparam;
	vector<double> fitparam_matrix(M);
	TF1 *fit_f = new TF1("fit_f", quadratic,0,5,n_fitparam); //range(0,5), 3 parameters
	//TODO:fitting procedure for each lattice, update fitparam_matrix
	for (unsigned i=0; i<list_size; i++){
		TH1D *hfit = new TH1D("temp_hist","htemporaneo",temp_steps,0,temp_steps);
		cout<<"[D]Temp hist created"<<endl;
		//populate the hist
		for(unsigned j=0; j<temp_steps; j++) hfit->Fill(j,binder_cumulants[index(i,j,temp_steps)]);
		hfit->Fit("fit_f");
		hfit->Draw(); //DEBUG
		TF1 *f_i_fit = hfit->GetFunction("f_i_fit");
		cout<<"[D] Temp fit function obtained from histogram:"<<endl;
		for(unsigned j=0; j<n_fitparam; j++) {
			fitparam_matrix[index(i,j,n_fitparam)]=f_i_fit->GetParameter(j);
			cout<<f_i_fit->GetParameter(j)<<" ";
		}
		cout<<endl;
		cout<<"    Fit "<<i<<"/"<<list_size<<" completed!"<<endl;
		delete hfit;
		delete f_i_fit;
		cout<<"[D]Temporary pointers deleted!"<<endl;
	}

	//once we have a fit we can approximate the  T_C finding the root of fit_i-fit_j;
	//since this should converge with L-> infinity, we compute only fit_i-fit_4
	//we use TF1::GetX(0,0,5) to find roots in the interval

	vector<double> criticalT(list_size-1);
	for(unsigned i=0; i< list_size-1 ; i++){
		double params[3] = {0,0,0};
		for (unsigned j=0; j<n_fitparam; j++) params[j]=fitparam_matrix[index(i,j,M)]-fitparam_matrix[index(list_size-1,j,M)];
		fit_f->SetParameters(params);
		criticalT[i]=fit_f->GetX(0,0,5) ; //revert inverse temp
	}
	cout<<"[+]Computing critical temperature T_c using the Binder's cumulant; "<<endl;
	cout<<"    The following values are the estimates for T_c (in growing order of reliability)"<<endl;
	cout<<"    The reference value is the fit of the Binder's cumulant for L="<<length_list[list_size-1]<<endl;
	for(unsigned i=0; i< list_size-1; i++) cout<<"    T_c(L="<<length_list[i]<<")="<<criticalT[i]<<endl;

	//TODO:
	//Now we compute the critical exponent for the heat capacity C \propto (T_C -T)^alpha
	//Heat capacity is computed using the FR theorem for equilibrium inside the SimulationModel class
	//This is done for all the lattices


	//TODO: raw output file close
	out_root.Close();
	return 0;

}

int simulate_driven_ising( const char* outfilename )
{
	//2. Driven ising model simulation
		//DrivenIsingModel *driven_ising = new DrivenIsingModel(8);

	return 0;
}

/*
int test_lattice() {
	//Test for class lattice.h

	//Tests for constructors
	int total_test=0;
	int failed_test=0;

	cout<<"[!]Testing class Lattice"<<endl;
	cout<<"[!]Testing constructor functions:"<<endl;
	cout<<"[+]Default constructor to create an empty lattice:"<<endl;
	Lattice l1;
	total_test++;
	cout<<"[+]Construction completed"<<endl;
	cout<<"[+]Default constructor and assignment to a new pointer to lattice: "<<endl;
	Lattice *l2= new Lattice();
    total_test++;
	cout<<"[+]Construction completed"<<endl;
	total_test++;
	failed_test++;
	cout<<"[!]No assignment/copy operations are allowed on lattice objects, since they contains const members"<<endl;
	cout<<"[+]Custom constructor with N=0"<<endl;
	try
	{
		Lattice l4(0);
		cout<<"[+]Custom constructor succeeded with N=0"<<endl;
		total_test++;
		failed_test++;
	}
	catch(...){
		cout<<"[!]Can't construct Lattice l4(0) "<<endl;
		total_test++;}
	cout<<"[+]Custom constructor with N=4, q=1"<<endl;
	try
	{
		Lattice l5(4,1);
		total_test++;
		failed_test++;
		cout<<"[+]Custom constructor succeeded with N=4, q=1"<<endl;
	}
	catch(...){
		cout<<"[!]Can't construct Lattice Lattice l5(4,1); "<<endl;
		total_test++;}
	cout<<"[+]Custom constructor with N=4,q=4"<<endl;
	try
	{
		Lattice l6(6);
		cout<<"[+]Custom constructor succeeded with N=6, q=4"<<endl;
		total_test++;
	}
	catch(...){
		cout<<"[!]Can't construct Lattice Lattice l6(4); "<<endl;
		failed_test++;
		total_test++;}

	//Tests for private data members: ONLY WITH DEBUG GETTER FUNCTIONS
	Lattice l6(150);
	cout<<"[+]Test for private data members:"<<endl;

	try	{
		cout<<"[+]The spin matrix is:"<<endl;
		//Best way to print a vector, see e.g. http://stackoverflow.com/questions/10750057/how-to-print-out-the-contents-of-a-vector/11335634#11335634
		//the std::copy function copy the vector to the std output
		vector<int> path1(l6.dbg_get_spin());
		copy(path1.begin(), path1.end(), ostream_iterator<int>(cout, " "));

		cout<<endl;
		cout<<"[+]Size of the spin matrix is:"<<path1.size()<<endl;
		cout<<"\n";
		total_test++;
	}
	catch(...)
	{
		cout<<"[!]Access denied to spin matrix "<<endl;
		total_test++;
		failed_test++;
	}
	try
	{
		cout<<"[+]The adjacency matrix is:";
		vector<int> path2(l6.dbg_get_weight());
		unsigned N = l6.get_dimension();
		//copy(path.begin(), path.end(), ostream_iterator<int>(cout, " "));
		for(unsigned i =0; i<path2.size(); i++) {
			if(i%N==0) cout<<"\n";
			cout<<path2[i]<<" ";
		}
		cout<<"\n";
		cout<<"[+]Size of the adjacency matrix is:"<<path2.size()<<endl;
		total_test++;
		int symmetry=0;
		for(unsigned i=0; i<N; i++)
		{
			for(unsigned j=0; j<N; j++)
					{
						if(path2[index(i,j,N)]!=path2[index(j,i,N)]) symmetry =1;
					}
		}

		if(symmetry==0) cout<<"[+] Adjacency matrix is symmetric"<<endl;
		else cout<<"[!]Adjacency matrix not symmetric"<<endl;

		cout<<"[+]Test of the neighbors() function:"<<endl;
		vector<vec_sz> neighs0 = l6.neighbors(0);
		cout<<"Neighbors of node 0 are:"<<endl;
		copy(neighs0.begin(), neighs0.end(),  ostream_iterator<int>(cout, " ") );
		cout<<endl;
	}
	catch(...)
	{
		cout<<"[!]Access denied to adjacency matrix "<<endl;
		total_test++;
		failed_test++;
	}

	try
	{
		cout<<"[+]The dimension is:";
		cout<<l6.get_dimension()<<endl;
		total_test++;
	}
	catch(...)
	{
		cout<<"[!]Access denied to dimension "<<endl;
		total_test++;
		failed_test++;
	}
	try
	{
		cout<<"[+]The coord number  is:";
		cout<<l6.get_coord_number()<<endl;
		total_test++;
	}
	catch(...)
	{
		cout<<"[!]Access denied to coord number "<<endl;
		total_test++;
		failed_test++;
	}

	return 0;
}
*/
