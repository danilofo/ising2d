/*
 * Graph.h
 *
 *  Created on: 08 gen 2017
 *      Author: danilo
 */
#if !defined(__CINT__) || defined(__MAKECINT__) //This directive is used by ROOT
#endif



#ifndef GRAPH_H_
#define GRAPH_H_

#include <TObject.h> //The primitive library for ROOT objects
#include <TRandom3.h> //A ROOT library for pseudo-random number
#include "Riostream.h"
#include <vector>
#include "Spin.h"

using namespace std;
typedef vector<int>::size_type vec_sz; //define a name for the type of the size of vector<int>
typedef vector<Spin>::size_type vd_sz; //TODO: substitute vec_sz -> vd_sz for vector<double>

class Graph: public TObject
{

	public:
	Graph(vd_sz =0, TRandom3* =NULL);
	~Graph();

	//Functionalities
	virtual void flip(const vec_sz )=0; //flip the i-th node
	virtual const vector<vec_sz> neighbors( const vec_sz )=0; //returns an array with neighbors' index
	const double getNodeValue(const vd_sz); //return the state of an individual spin

	virtual void reset(const char * flag="random")=0;
	//TODO: implement a reset(const char * flag ) function (identical to initSpins;
	//Needed to simulate an ising model several times with the same lattice

	//Public getter functions:
	vd_sz getDimension();

	//Public setter functions:

	//protected members can be accessed by the derived class
	protected:
	vector<Spin> spin; //std::vector instead of C-like array! contains matrix of spins, allocated dynamically
					  //for this choice see: http://stackoverflow.com/questions/381621/using-arrays-or-stdvectors-in-c-whats-the-performance-gap
	vector<int> weight; // contains the constant adjacency matrix of the lattice, allocated dynamically

	//Utilities
	const vec_sz index(const vec_sz ,const vec_sz) const; //return the correct index for the 1d representation of matrices

	TRandom3* Rnd;


	vd_sz dimension; // For the use of size_t instead of unsigned int
					//see e.g. http://stackoverflow.com/questions/1951519/when-should-i-use-stdsize-t



	ClassDef(Graph,1);
};

#endif /* GRAPH_H_ */
