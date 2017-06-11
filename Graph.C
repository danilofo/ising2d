/*
 * Graph.C
 *
 *  Created on: 08 gen 2017
 *      Author: danilo
 */
#include "Graph.h"

Graph::Graph(vd_sz d, TRandom3* rnd): dimension(d), Rnd(rnd),spin(d,0),weight(d*d,0){}

Graph::~Graph(){}

vd_sz Graph::getDimension(){
	return this->dimension;
}
vec_sz Graph::index(const vec_sz i, const vec_sz j) const
{	//Return the correct index for the matrix element A_{ij}
	//Return the correct index for 1d representation of 2d matrices; const because does not change data members
	//See http://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
	const vec_sz  &n_rows = this->dimension;
	return  i*n_rows +  j;
}
double Graph::getNodeValue(const vd_sz i){ //return the state of an individual spin
	double spin_i=0;
	const vd_sz  &N=this->dimension;
	if( i<N )spin_i= this->spin[i].getSpinValue();
	else{
		cout<<"[!]Graph:Node "<<i<<" is not in the graph"<<endl;
		return -100.;
	}
	return spin_i;
}

void Graph::flip(const vd_sz i){
	const vd_sz N = this->getDimension();
	if(i>N-1){
		cout<<"[!]Graph:Node "<<i<<" is not in the graph (dimension="<<N<<")"<<endl;
	}
	else{
		this->spin[i].flipSpin();
	}
}
