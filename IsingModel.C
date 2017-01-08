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
IsingModel::IsingModel(): coupling(0),lattice(){}
IsingModel::IsingModel(vec_sz lenght, double J, const char* type): coupling(J), lattice(lenght){}

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

const double IsingModel::simulate(double beta, //inverse temperature
		unsigned n_iterations)
//Metropolis Monte Carlo simulation of a model defined by an hamiltonian on a generic graph;
//Randomly choose a spin to flip; accept the change in the configuration according to the Metropolis choice
{
	double magnetization = 0 ;
if(beta < 0){
		cout<<"[!]Invalid inverse temperature, must be greater than zero."<<endl;
	}
else
	{
		//initialization
		vec_sz N = this->lattice.get_dimension();
		double H0 = this->hamiltonian(); //energy of the current distribution
		double NewH=0; //energy of the proposed distribution
		TRandom *Rnd = new TRandom3(1750); //random number generator; one for each instance of IsingModel!
		vec_sz index =0;
		double delta=0;
		double acceptance=0;
		double boltz=0;
		int spin_value_i0=0;
		//initial value of magnetization
		for (vec_sz i=1; i<N; i++) magnetization+= this->lattice.spin_value(i)/N;
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
				magnetization-=  2*spin_value_i0/N;
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
return magnetization;
}



