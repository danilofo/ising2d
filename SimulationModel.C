/*
 * SimulationModel.C
 *
 *  Created on: 08 gen 2017
 *      Author: danilo
 */
#include "SimulationModel.h"

//implementazione della classe
ClassImp(SimulationModel)


//public constructors
SimulationModel::SimulationModel():graph(NULL),M(0),E(0),Rnd(NULL){}

SimulationModel::~SimulationModel(){
    if(Rnd!=NULL) delete Rnd;
    if(graph!=NULL) delete graph;
}


const double SimulationModel::getMagnetization() const {
	return this->M;
}

const double SimulationModel::getEnergy() const {
	return this->E;
}

void SimulationModel::simulate(double beta, //inverse temperature
		unsigned n_iterations){
//Metropolis Monte Carlo simulation of a model defined by an hamiltonian on a generic graph;
//Randomly choose a spin to flip; accept the change in the configuration according to the Metropolis choice
	//cout<<"[D]IsingModel:simulate started"<<endl;
	double magnetization = this->M ;
	double H0 = this->E; //energy of the current distribution
	if(beta < 0){
		cout<<"[!]IsingModel:Invalid inverse temperature, must be greater than zero."<<endl;
	}
	else
	{
		//initialization
		vd_sz N = this->graph->getDimension(); //total number of spin
		double NewH=0; //energy of the proposed distribution
		vec_sz index =0;
		double delta=0;
		double acceptance=0;
		double boltz=0;
		double spin_value_i0=0;
		double spin_value_i1 =0;

		//loop until n_iterations is reached
		for(unsigned i=1; i<n_iterations; i++)
		{
			//cout<<"[D]	IsingModel:Rnd called...";
			index= Rnd->Integer(N);
			//proposal
			spin_value_i0=this->graph->getNodeValue(index);
			//cout<<"[D]	IsingModel:.flip() called..."<<endl;
			this->graph->flip(index);
			//cout<<"[D]	IsingModel:.flip ok"<<endl;
			spin_value_i1=this->graph->getNodeValue(index);
			NewH=this->hamiltonian();
			delta=H0-NewH;//negative if NewH>H0;energy improvements should always be accepted!
			boltz=std::exp( beta * delta );
			acceptance = std::min(1.0, boltz);
			//metropolis choice
			if (acceptance==1)
			{
				H0+=delta;
				magnetization-= this->energyVar(spin_value_i0,spin_value_i1); //THIS DEPENDS ON THE HAMILTONIAN !
				//Maybe can be handled by a function of IsingModel newMagnetization(oldSpin, newSpin)
			}
			else
			{
				if(Rnd->Rndm()<acceptance )
				{
					H0+=delta;
					magnetization-= 2*spin_value_i0/N;
				}
				else this->graph->flip(index);
			}
		}
	}
	//assignment of final energy / magn
	this->M=magnetization;
	this->E=H0;

	return;
}
