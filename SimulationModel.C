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
double SimulationModel::getMagnetization() const { return this->M;}
double SimulationModel::getEnergy() const {	return this->E;}
void SimulationModel::flipSpin(vd_sz i){
	double old_v = this->graph->getNodeValue(i);
	this->graph->flip(i);
	double new_v = this->graph->getNodeValue(i);
	this->M+= mVar(old_v,new_v);
	this->E+=energyVar(i,old_v,new_v);
}
void SimulationModel::simulate(double beta, //inverse temperature
		vd_sz n_iterations){
	//Metropolis Monte Carlo simulation of a model defined by an hamiltonian on a generic graph;
	//Randomly choose a spin to flip; accept the change in the configuration according to the Metropolis choice
	//cout<<"[D]IsingModel:simulate started"<<endl;
	if(graph==NULL){
		cout<<"[!]SimulationModel:simulate:invalid pointer graph received"<<endl;
		return;
	}
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
		double E0 = this->getEnergy();
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
			//energy after update
			E1=this->getEnergy();
			if(E0-E1<0){ // we are accepting with prob boltz_w
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
	}
	return ;
}
