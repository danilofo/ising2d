/*
 * IsingModel.C
 *
 *  Created on: 07 gen 2017
 *      Author: danilo
 */
#include "IsingModel.h"

//implementazione della classe
ClassImp(IsingModel)

//public constructors
IsingModel::IsingModel(): coupling(0),lattice(),M(0),E(0),Rnd(){}
IsingModel::IsingModel(vec_sz length, double J, const char* type): coupling(J),
		lattice(),M(0.),E(0.),Rnd(){
	//initiallize
	//TODO: CHECK THIS !
	this->Rnd = new TRandom3(123); //TODO: fix the seed choice
	this->lattice = Lattice(length,Rnd);
	this->E=this->hamiltonian();
	vd_sz N = length*length; 
	//initial value of magnetization
	for (vec_sz i=1; i<N; i++) (this->M)+=(this->lattice.spin_value(i)/N);
}
IsingModel::~IsingModel(){
	delete Rnd;
}
//Utilities
void IsingModel::resetLattice(){
	this->lattice.reset();
}
void IsingModel::newLattice(vec_sz new_length){
	this->lattice = Lattice(new_length, Rnd);
}
//functions needed to simulate
const double IsingModel::hamiltonian()
{
	double H=0;
	const vec_sz N = this->lattice.get_dimension();
	for(vec_sz i=0; i<N ; i++)
	{
		const vector<vec_sz> neighs_i=this->lattice.neighbors(i);
		const unsigned q = this->lattice.get_coord_number();
		for(unsigned j=0; j<q; j++)
		{
			//compute the hamiltonian; since the graph is not directed, each node appear twice
			//in the neighbors list; to speed up the computation we use only nodes with j>i instead
			//of dividing the result by 2
			if(neighs_i[j]>i) H+= this->coupling * this->lattice.spin_value(neighs_i[j]);
		}
	}
	return H;
}

const double IsingModel::getMagnetization() const {
	return this->M;
}

const double IsingModel::getEnergy() const {
	return this->E;
}

void IsingModel::simulate(double beta, //inverse temperature
		unsigned n_iterations)
//Metropolis Monte Carlo simulation of a model defined by an hamiltonian on a generic graph;
//Randomly choose a spin to flip; accept the change in the configuration according to the Metropolis choice
{
	double magnetization = this->M ;
	double H0 = this->E; //energy of the current distribution
	if(beta < 0){
		cout<<"[!]Invalid inverse temperature, must be greater than zero."<<endl;
	}
	else
	{
		//initialization
		vd_sz N = this->lattice.get_dimension(); //total number of spin 
		double NewH=0; //energy of the proposed distribution
		vec_sz index =0;
		double delta=0;
		double acceptance=0;
		double boltz=0;
		double spin_value_i0=0;

		//loop until n_iterations is reached
		for(unsigned i=1; i<n_iterations; i++)
		{
			index= Rnd->Integer(N);
			//proposal
			spin_value_i0=this->lattice.spin_value(index);
			this->lattice.flip(index);
			NewH=this->hamiltonian();
			delta=H0-NewH;//negative if NewH>H0;energy improvements should always be accepted!
			boltz=std::exp( beta * delta );
			acceptance = std::min(1.0, boltz);
			//metropolis choice
			if (acceptance==1)
			{
				H0+=delta;
				magnetization-= 2.0*spin_value_i0/N; //THIS DEPENDS ON THE HAMILTONIAN !
				//Maybe can be handled by a function of IsingModel newMagnetization(oldSpin, newSpin)
			}
			else
			{
				if(Rnd->Rndm()<acceptance )
				{
					H0+=delta;
					magnetization-= 2*spin_value_i0/N;
				}
				else this->lattice.flip(index);
			}
		}
	}
	//assignment of final energy / magn
	this->M=magnetization;
	this->E=H0;
	return;
}



