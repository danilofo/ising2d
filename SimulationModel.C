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
SimulationModel::SimulationModel():graph(NULL),M(0.),E(0.),Rnd(NULL){}

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
