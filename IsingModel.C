/*
 * IsingModel.C
 *
 *  Created on: 07 gen 2017
 *      Author: danilo
 */
#include "IsingModel.h"
#include "SimulationModel.h"


//implementazione della classe
ClassImp(IsingModel)

//public constructors
IsingModel::IsingModel(vec_sz length, double J): SimulationModel(),lattice(NULL),coupling(J)
{
	//initialize
	//cout<<"[D]IsingModel: constructing TRandom3 Rnd...";
	this->Rnd = new TRandom3(823); //TODO: fix the seed choice
	//cout<<"[D]IsingModel: Rnd...completed"<<endl;
	//cout<<"[D]IsingModel: constructing Lattice lattice...";
	this->lattice = new Lattice(length,Rnd);
	this->graph = this->lattice;
	//cout<<"[D]IsingModel: lattice...completed"<<endl;
	//cout<<"[D]IsingModel: computing initial energy...";
	this->E=this->hamiltonian();
	//cout<<"[D]IsingModel: hamiltonian...completed"<<endl;
	vd_sz N = length*length;
	//initial value of magnetization
	for (vec_sz i=1; i<N; i++)
		(this->M)+=(this->lattice->getNodeValue(i)/N);

	//cout<<"[D]IsingModel:Construction completed"<<endl;
}
//functions needed to simulate
const double IsingModel::hamiltonian()
{
    double H=0;
    const vec_sz N = this->lattice->getDimension();
    for(vec_sz i=0; i<N ; i++)
    {
        const vector<vec_sz> neighs_i=this->lattice->neighbors(i);
        const unsigned q = this->lattice->getCoordN() ;
        for(unsigned j=0; j<q; j++)
        {
            //compute the hamiltonian; since the graph is not directed, each node appear twice
            //in the neighbors list; to speed up the computation we use only nodes with j>i instead
            //of dividing the result by 2
            if(neighs_i[j]>i) H+= this->coupling * this->lattice->getNodeValue(neighs_i[j]);
        }
    }
    return H;
}

const double energyVar(double old_val, double new_val){
	return (-2.)*(old_val);
}

void resetGraph(){
	this->lattice->reset();
}
virtual void newGraph(vec_sz N){
	this->lattice = new Lattice(N,Rnd);
	this->graph = this->lattice;
}



