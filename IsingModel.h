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
#include "lattice.h"
#include <algorithm>

typedef std::vector<int>::size_type vec_sz;

class IsingModel: public TObject {
	public:
	//public constructors
	IsingModel();
	IsingModel(vec_sz lenght, double J=1., const char* type="lattice" );

	//simulation
	//ALL FUNCTIONS HERE SHOULD IGNORE THE TYPE OF GRAPH
	//HOW TO DO THE SAME FOR THE HAMILTONIAN?
	const double simulate(double beta, unsigned n_iterations=1000);
	const double hamiltonian();

	private:
	//The default construction of a Lattice object should be replaced with a Graph object!
	//At runtime the constructor IsingModel should assign the correct object;
	Lattice lattice;
	const double coupling;

	ClassDef(IsingModel,1); //Used by to define a class ROOT

};



#endif /* ISINGMODEL_H_ */
