
/*
 * IsingModel.h
 *
 *  Created on: 07 gen 2017
 *      Author: danilo
 */
#if !defined(__CINT__) || defined(__MAKECINT__) //This directive is used by ROOT
#endif


#ifndef ISINGMODEL_H_
#define ISINGMODEL_H_
#include "Lattice.h"
#include "SimulationModel.h"
#include <algorithm>


class IsingModel: public SimulationModel {
	//An Ising Model defined on a lattice (class Lattice)

	public:
	//public constructors
	IsingModel(vec_sz length, double J=1.);
	~IsingModel();

	//simulation
	const double hamiltonian();
	const double energyVar(double old_val, double new_val);

	//actions on the lattice
    void resetGraph(); //
    void newGraph(vec_sz N); //

	private:
    Lattice* lattice;
	const double coupling;

	ClassDef(IsingModel,1); //Used by to define a class ROOT

};



#endif /* ISINGMODEL_H_ */
