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


#include <TObject.h> //The primitive library for ROOT objects
#include <TRandom3.h> //A ROOT library for pseudo-random number
#include "Riostream.h"
#include <vector>

using namespace std;
typedef vector<int>::size_type vec_sz; //define a name for the type of the size of vector<int>
typedef vector<double>::size_type vd_sz; //TODO: substitute vec_sz -> vd_sz for vector<double>


class Lattice: public TObject {
	//DESCRIPTION:
	//This class instantiates a lattice in which each node (spin) can switch between two states;
	//the class gives access to each node in the lattice and provides a neighbors' list for each of them
	//
	public:
		//NO PUBLIC DATA MEMBER
		//Member functions:

		// public constructors
		Lattice();
		Lattice(int ,TRandom* , int q=4,const char * flag="random");
		//copy and assignment operators
		//Lattice(const Lattice& other);
		//Lattice operator=(const Lattice rhs); //(const &) should not be used here! see e.g http://www.cplusplus.com/articles/y8hv0pDG/

		//Functionalities
		void flip(const vec_sz ); //flip the i-th node
		const vector<vec_sz> neighbors( const vec_sz ) const; //returns an array with neighbors' index
		const double spin_value(const vd_sz); //return the state of an individual spin

		void reset(const char * flag="random");
		//TODO: implement a reset(const char * flag ) function (identical to initSpins;
		//Needed to simulate an ising model several times with the same lattice

		//DEBUG ONLY GETTER FUNCTIONS

		vector<double> dbg_get_spin();
		vector<int> dbg_get_weight();

		//Public getter functions:
		vec_sz get_dimension();
		unsigned int get_coord_number();

	private:
		//Data members
		vector<double> spin; //std::vector instead of C-like array! contains matrix of spins, allocated dynamically
					      //for this choice see: http://stackoverflow.com/questions/381621/using-arrays-or-stdvectors-in-c-whats-the-performance-gap
		vector<int> weight; // contains the constant adjacency matrix of the lattice, allocated dynamically
		vec_sz dimension; // For the use of size_t instead of unsigned int
								//see e.g. http://stackoverflow.com/questions/1951519/when-should-i-use-stdsize-t
		vec_sz lenght;
		unsigned int coord_number;
		TRandom* Rnd;
		//Private member functions:
		//initialization procedures
		void initSpins(const char * );
		void initRectangularLattice();

		//Utilities
		const vec_sz index(const vec_sz ,const vec_sz) const; //return the correct index for the 1d representation of matrices

	ClassDef(Lattice,1); //Used by to define a class ROOT
	}; //class Lattice


#endif /* LATTICE_H_ */
