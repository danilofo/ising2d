/*
 * SimulationModel.C
 *
 *  Created on: 08 gen 2017
 *      Author: danilo
 */
#include "SimulationModel.h"
#include "Riostream.h"
using namespace std; 

//implementazione della classe
ClassImp(SimulationModel)
//public constructors
SimulationModel::SimulationModel():graph(NULL),M(0),E(0),Rnd(NULL){}
SimulationModel::~SimulationModel(){
	//TRandom3* Rnd is created in this class, so the destructor deletes it.
	//The pointer Graph* graph, instead, is assigned to an object only in the
	//derived class, so it should not be deleted.
    if(Rnd!=NULL) delete Rnd;
}
void SimulationModel::setEnergy(double e) {this->E=e;}
void SimulationModel::setMagnetization(double m ) {this->M=m;}
double SimulationModel::getMagnetization() const { 	return this->M;}
double SimulationModel::getEnergy() const {	return this->E;}
void SimulationModel::flipSpin(vd_sz i){
	double old_v = this->graph->getNodeValue(i);
	this->graph->flip(i);
	double new_v = this->graph->getNodeValue(i);
	this->setMagnetization(this->M+magnVar(old_v,new_v));
}
void SimulationModel::simulate(double beta, //inverse temperature
		unsigned n_iterations){
	//Metropolis Monte Carlo simulation of a model defined by an hamiltonian on a generic graph;
	//Randomly choose a spin to flip; accept the change in the configuration according to the Metropolis choice
	//cout<<"[D]IsingModel:simulate started"<<endl;
	if(graph==NULL){
		cout<<"[!]SimulationModel:simulate:invalid pointer graph received"<<endl;
		return;
	}
	double magnetization = this->getMagnetization();
	cout<<"[D]Current magn:"<<magnetization<<endl;
	double H0 = this->getEnergy(); //energy of the current distribution
	if(beta < 0){
		cout<<"[!]IsingModel:Invalid inverse temperature, must be greater than zero."<<endl;
	}
	else
	{
		//initialize loop variables
		vd_sz dim = this->graph->getDimension();
		vd_sz node_i= 0;
		double spin_i0 = 0;
		double spin_i1 = 0;
		double E0 = this->hamiltonian();
		double E1= 0;
		double boltz_w =0;

		for(unsigned i=0; i<n_iterations; i++){
			node_i = Rnd->Integer(dim);
			//save current state of node i
			spin_i0 = this->graph->getNodeValue(node_i);
			//flip current spin
			this->flipSpin(node_i);
			//save new value (general)
			spin_i1 = this->graph->getNodeValue(node_i);
			//update E1
			E1=this->hamiltonian();
			if(E0<E1){ // we are accepting with prob boltz_w
				boltz_w = std::exp(beta*(E0-E1));
				if(Rnd->Rndm() < boltz_w) { //we accept an higher energy
					E0=E1;
				}
				else //we are rejecting the increase in energy
					this->flipSpin(node_i);
			}
			else{ //we are accepting with probability one
				E0=E1; //update current minimum
			}
		}
	this->setEnergy(E0);
	}

	return ;
}















		/*


		//initialization
		vd_sz N = this->graph->getDimension(); //total number of spin
		double NewH=0; //energy of the proposed distribution
		vec_sz index =0;
		double delta=0;
		double acceptance=0;
		double boltz=0;
		double spin_value_i0=0;
		double spin_value_i1=0;
		//loop until n_iterations is reached
		for(unsigned i=1; i<n_iterations; i++)
		{
			//index=(index+1)%N;
			index= Rnd->Integer(N);
			//proposal
			//spin_value_i0=this->graph->getNodeValue(index);
			this->graph->flip(index);
			//spin_value_i1=this->graph->getNodeValue(index);
			//magnetization+= this->magnVar(spin_value_i0,spin_value_i1);
			NewH=this->hamiltonian();
			delta=H0-NewH;//negative if NewH>H0; improvements should always be accepted!
			boltz=std::exp( beta * delta );
			acceptance = std::min(1.0, boltz);
			//metropolis choice
			if (acceptance==1)
			{
				H0=NewH;
				magnetization=this->magnetization();
				//magnetization+=this->magnVar(spin_value_i1,spin_value_i0);
			}
			else
			{
				if(Rnd->Rndm()<acceptance )
				{
					H0=NewH;
					magnetization=this->magnetization();
					//magnetization+= this->magnVar(spin_value_i1,spin_value_i0);
				}
				else this->graph->flip(index);
			}
		}
	}
	//assignment of final energy / magn
	this->setMagnetization(magnetization);
	this->setEnergy(H0);
	return;
}*/
