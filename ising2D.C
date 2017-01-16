//============================================================================
// Name        : ising2D.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "lattice.h"

#include <vector>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator

using namespace std;

vec_sz index(const vec_sz x, const vec_sz y,unsigned N)
{	//Return the correct index for 1d representation of 2d matrices; do not change data members
	//See http://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
	return  x + N* y;
}

int test_lattice() {
	//Test for class lattice.h

	//Tests for constructors
	int total_test=0;
	int failed_test=0;

	cout<<"[!]Testing class Lattice"<<endl;
	cout<<"[!]Testing constructor functions:"<<endl;
	cout<<"[+]Default constructor to create an empty lattice:"<<endl;
	Lattice l1;
	total_test++;
	cout<<"[+]Construction completed"<<endl;
	cout<<"[+]Default constructor and assignment to a new pointer to lattice: "<<endl;
	Lattice *l2= new Lattice();
    total_test++;
	cout<<"[+]Construction completed"<<endl;
	total_test++;
	failed_test++;
	cout<<"[!]No assignment/copy operations are allowed on lattice objects, since they contains const members"<<endl;
	cout<<"[+]Custom constructor with N=0"<<endl;
	try
	{
		Lattice l4(0);
		cout<<"[+]Custom constructor succeeded with N=0"<<endl;
		total_test++;
		failed_test++;
	}
	catch(...){
		cout<<"[!]Can't construct Lattice l4(0) "<<endl;
		total_test++;}
	cout<<"[+]Custom constructor with N=4, q=1"<<endl;
	try
	{
		Lattice l5(4,1);
		total_test++;
		failed_test++;
		cout<<"[+]Custom constructor succeeded with N=4, q=1"<<endl;
	}
	catch(...){
		cout<<"[!]Can't construct Lattice Lattice l5(4,1); "<<endl;
		total_test++;}
	cout<<"[+]Custom constructor with N=4,q=4"<<endl;
	try
	{
		Lattice l6(6);
		cout<<"[+]Custom constructor succeeded with N=6, q=4"<<endl;
		total_test++;
	}
	catch(...){
		cout<<"[!]Can't construct Lattice Lattice l6(4); "<<endl;
		failed_test++;
		total_test++;}

	//Tests for private data members: ONLY WITH DEBUG GETTER FUNCTIONS
	Lattice l6(150);
	cout<<"[+]Test for private data members:"<<endl;

	try	{
		cout<<"[+]The spin matrix is:"<<endl;
		//Best way to print a vector, see e.g. http://stackoverflow.com/questions/10750057/how-to-print-out-the-contents-of-a-vector/11335634#11335634
		//the std::copy function copy the vector to the std output
		vector<int> path1(l6.dbg_get_spin());
		copy(path1.begin(), path1.end(), ostream_iterator<int>(cout, " "));

		cout<<endl;
		cout<<"[+]Size of the spin matrix is:"<<path1.size()<<endl;
		cout<<"\n";
		total_test++;
	}
	catch(...)
	{
		cout<<"[!]Access denied to spin matrix "<<endl;
		total_test++;
		failed_test++;
	}
	try
	{
		cout<<"[+]The adjacency matrix is:";
		vector<int> path2(l6.dbg_get_weight());
		unsigned N = l6.get_dimension();
		//copy(path.begin(), path.end(), ostream_iterator<int>(cout, " "));
		for(unsigned i =0; i<path2.size(); i++) {
			if(i%N==0) cout<<"\n";
			cout<<path2[i]<<" ";
		}
		cout<<"\n";
		cout<<"[+]Size of the adjacency matrix is:"<<path2.size()<<endl;
		total_test++;
		int symmetry=0;
		for(unsigned i=0; i<N; i++)
		{
			for(unsigned j=0; j<N; j++)
					{
						if(path2[index(i,j,N)]!=path2[index(j,i,N)]) symmetry =1;
					}
		}

		if(symmetry==0) cout<<"[+] Adjacency matrix is symmetric"<<endl;
		else cout<<"[!]Adjacency matrix not symmetric"<<endl;

		cout<<"[+]Test of the neighbors() function:"<<endl;
		vector<vec_sz> neighs0 = l6.neighbors(0);
		cout<<"Neighbors of node 0 are:"<<endl;
		copy(neighs0.begin(), neighs0.end(),  ostream_iterator<int>(cout, " ") );
		cout<<endl;
	}
	catch(...)
	{
		cout<<"[!]Access denied to adjacency matrix "<<endl;
		total_test++;
		failed_test++;
	}

	try
	{
		cout<<"[+]The dimension is:";
		cout<<l6.get_dimension()<<endl;
		total_test++;
	}
	catch(...)
	{
		cout<<"[!]Access denied to dimension "<<endl;
		total_test++;
		failed_test++;
	}
	try
	{
		cout<<"[+]The coord number  is:";
		cout<<l6.get_coord_number()<<endl;
		total_test++;
	}
	catch(...)
	{
		cout<<"[!]Access denied to coord number "<<endl;
		total_test++;
		failed_test++;
	}
	return 0;
}
