/*
 * SimulationModel.h
 *
 *  Created on: 08 gen 2017
 *      Author: danilo
 */
#if !defined(__CINT__) || defined(__MAKECINT__) //This directive is used by ROOT
#endif

#ifndef SIMULATIONMODEL_H_
#define SIMULATIONMODEL_H_

#include <algorithm>
#include "Graph.h"

typedef std::vector<int>::size_type vec_sz;
typedef vector<Spin>::size_type vd_sz;

class SimulationModel//: public TObject
{
    // classe in cui definisco il grafo
    
    public:
	SimulationModel();
	virtual ~SimulationModel();

	//
    void simulate(double beta, unsigned n_iterations=1000);
    virtual double hamiltonian()=0;
    virtual double magnVar(double old_val, double new_val)=0;

    //Actions on graph
    virtual void resetGraph()=0; //
    virtual void newGraph(vd_sz ,const char*)=0; //

    double getEnergy() const ;
    double getMagnetization() const;
    


    protected:
    Graph* graph;
	TRandom3* Rnd;

	//debug only
    void setEnergy(double ) ;
    void setMagnetization(double ) ;

	double E;//current energy
	double M; //current magnetization

	ClassDef(SimulationModel,1); //Used by to define a class ROOT
};

#endif /* SIMULATIONMODEL_H_ */
