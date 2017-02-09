/*
 * Graph.C
 *
 *  Created on: 08 gen 2017
 *      Author: danilo
 */

#include "Graph.h"

Graph::Graph(vd_sz d, TRandom3* rnd): dimension(d), Rnd(rnd),spin(d,0),weight(d*d,0){}

Graph::~Graph(){if(Rnd!=NULL) delete Rnd;}

vd_sz Graph::getDimension(){
	return this->dimension;
}

const vec_sz Graph::index(const vec_sz x, const vec_sz y) const
{	//Return the correct index for 1d representation of 2d matrices; do not change data members
	//See http://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
	const vec_sz  &n_rows = this->dimension;
	return  x + n_rows* y;
}

const double Graph::getNodeValue(const vd_sz i){ //return the state of an individual spin
	double spin_i=0;
	const vd_sz  &N=this->dimension;
	if( i<N )spin_i= this->spin[i].getSpinValue();
	else cout<<"[!]Graph:Node "<<i<<" is not in the graph"<<endl;
	return spin_i;
}

