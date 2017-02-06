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
SimulationModel::SimulationModel(): coupling(0),graph(),M(0),E(0),Rnd(){}
SimulationModel::SimulationModel(vec_sz length, double J, const char* type): coupling(J),
graph(),M(0.),E(0.),Rnd(){
    //initiallize
    //TODO: CHECK THIS !
    this->Rnd = new TRandom3(123); //TODO: fix the seed choice
    this->graph = Graph(length,Rnd);
    this->E=this->hamiltonian();
    vd_sz N = length*length;
    //initial value of magnetization
    for (vec_sz i=1; i<N; i++) (this->M)+=(this->lattice.node_value(i)/N);
}
SimulationModel::~SimulationModel(){
    delete Rnd;
}
//Utilities
void SimulationModel::resetGraph(){
    this->graph.reset();
}
void SimulationModel::newGraph(vec_sz new_length){
    this->graph = Graph(new_length, Rnd);
}
