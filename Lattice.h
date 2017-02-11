/*
 * lattice.h
 *
 *  Created on: 04 gen 2017
 *      Author: danilo
 */
#if !defined(__CINT__) || defined(__MAKECINT__) //This directive is used by ROOT
#endif

#ifndef LATTICE_H_
#define LATTICE_H_


#include "Graph.h"


class Lattice: public Graph {
	//DESCRIPTION:
	//This class instantiates a lattice in which each node (spin) can switch between two states;
	//the class gives access to each node in the lattice and provides a neighbors' list for each of them

	public:
		//Member functions:

		// public constructors
		Lattice(int ,TRandom3* , int q=4,const char * flag="random");
		~Lattice();
		
		//Functionalities

		virtual const vector<vd_sz> neighbors( const vd_sz ) ; //returns an array with neighbors' index
		virtual void reset(const char * flag="random");
		unsigned getCoordN();

		virtual void flip(const vd_sz ); //flip the i-th node

	private:
		//Data members
		vd_sz length;
		unsigned int coord_number;


		//Private member functions:
		//initialization procedures
		void initSpins(const char * );
		void initRectangularLattice();

	ClassDef(Lattice,1); //Used by to define a class ROOT
	}; //class Lattice


#endif /* LATTICE_H_ */
