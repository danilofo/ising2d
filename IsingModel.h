
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
#include <algorithm>

typedef std::vector<int>::size_type vec_sz;

class IsingModel: public TObject {
	//An Ising Model defined on a lattice (class Lattice)

	public:
	//public constructors
	IsingModel();
	IsingModel(vec_sz length, double J=1., const char* type="lattice" );
	~IsingModel();

	//simulation
	//ALL FUNCTIONS HERE SHOULD IGNORE THE TYPE OF GRAPH
	//HOW TO DO THE SAME FOR THE HAMILTONIAN?
	void simulate(double beta, unsigned n_iterations=1000);
	const double hamiltonian();
	void resetLattice(); //reset lattice function
	void newLattice(vec_sz new_length); // new lattice;


	const double getEnergy() const;
	const double getMagnetization() const;

	private:
	//The default construction of a Lattice object should be replaced with a Graph object!
	//At runtime the constructor IsingModel should assign the correct object;
	Lattice lattice;
	const double coupling;

	double E;//current energy
	double M; //current magnetization

	TRandom* Rnd;
	ClassDef(IsingModel,1); //Used by to define a class ROOT

};



#endif /* ISINGMODEL_H_ */
