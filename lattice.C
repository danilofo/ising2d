#include "/home/danilo/workspace/ising2D/src/lattice.h"

//implementazione della classe
ClassImp(Lattice)

//constructors
Lattice::Lattice(): lenght(0),dimension(0),coord_number(0),spin(),weight(){}

Lattice::Lattice(int N, int q,const char *flag):  lenght(N),dimension(N*N), coord_number(q),spin(N*N,0),weight(N*N*N*N,0)
{	//Construct a lattice with N spin IN EACH ROW/COLUMN
	if (N<1 )
	{
		cout<<"[!]Initialization failed. Dimension must be greater than 0."<<endl;
		throw N ;
	}
	else if(q<2)
	{
		cout<<"[!]Initialization failed. Coordination number must be greater than 1."<<endl;
		throw q ;
	}
	if(coord_number==4)
	{
		this->initRectangularLattice(); //fill  weight
		this->initSpins(flag); //fill spin
	}
	else
	{
		cout<<"[!]A lattice with coordination number q="<<coord_number<<" is not yet available."<<endl;
		cout<<"[+]Admitted choices: q=4"<<endl;
		throw q;
	}
}


//utilities
const vec_sz Lattice::index(const vec_sz x, const vec_sz y) const
{	//Return the correct index for 1d representation of 2d matrices; do not change data members
	//See http://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
	const vec_sz &n_rows = this->dimension;
	return  x + n_rows* y;
}

//Public getter functions
vec_sz Lattice::get_dimension()
{	vec_sz N = this->dimension;
	return N;
	}
unsigned int Lattice::get_coord_number()
{	unsigned int q = this->coord_number;
	return q;
	}

// DEBUG ONLY FUNCTIONS
vector<int> Lattice::dbg_get_spin()
	{
	vector<int> result(this-> spin);
	return result;}
vector<int> Lattice::dbg_get_weight()
{
	vector<int> result(this->weight);
	return result;
	}


//initialization procedures:
void Lattice::initSpins(const char * flag)
//Initialize the value of each spin ; choice allowed are
//"random" or "ordered". The latter means that all spins are assigned to +1
{	const vector<int> test_vc(this->dimension,0);
	if (test_vc==spin ) //protect against re-initializations AND prevent assignment of spin
	{
		//const vec_sz mat_dim = this->dimension*this->dimension; //Read dimension
		//vector<int> Lattice::spin(this->dimension); //Declare and allocate the spin array

		if (strncmp(flag,"random",10)==0)
		{	delete gRandom;
			TRandom* gRandom = new TRandom3(123);
			for(unsigned i=0; i<this->dimension; i++) //initialize a matrix using a single for loop
			{
				this->spin[i]= (gRandom->Rndm() < 0.5) ? -1 : 1 ; //assign the value +1 with probability 0.5
			}
		}

		else if (strncmp(flag,"ordered",10)==0)
		{	for(unsigned i=0; i<this->dimension; i++)
			{
				this->spin[i]=1;
			}
		}
		else
		{
			cout<<"[!]Initialization failed. Admitted choices : \"random\", \"ordered\" "<<endl;
		}
	}
	else{cout<<"[!]Re-initialization attempt ignored"<<endl;}
}

void Lattice::initRectangularLattice()
{

	const vec_sz &N = this->dimension;
	const vec_sz &L = this->lenght;
	const vector<int> test_vc(N*N,0);
	if(test_vc==this->weight)
	{
		for (unsigned i=0; i<N; i++)
		{
			//this is the correct way to assign the neighborhood relation:
			//node i is assigned neighbor of i+1 mod N and i+(node in a row) mod N,
			//that is, the nodes respectively on the right and at the bottom
			//of the currently selected node. Since the assignment is symmetric,
			//we only assign 2 neighbors at each iteration.
			//The other 2 are assigned when the top or left nodes are selected.
			this->weight[index((i+1)%N,i)] = 1;
			this->weight[index(i,(i+1)%N)] = 1;
			this->weight[index((i+L)%N,i)] = 1;
			this->weight[index(i,(i+L)%N)] = 1;
		}
	}
	else{cout<<"[!]Re-initialization attempt ignored"<<endl;}
}


//class functionalities:
void Lattice::flip(const vec_sz i)
	{
		const vec_sz &N = this->dimension;
		if(i>N-1)
		{
			cout<<"Node "<<i<<" is not in the lattice (dimension="<<N<<")"<<endl;
		}
		else{
			this->spin[i]-=2*this->spin[i];
		}
	}

const vector<vec_sz> Lattice::neighbors( const vec_sz i) const
{
	const unsigned int &q = this->coord_number;
	const vec_sz &N = this->dimension;

	if (i>N-1) //check if i is in the lattice
	{
		cout<<"Node "<<i<<" is not in the lattice (dimension="<<N<<")"<<endl;
		vector<vec_sz> neighs;
		return neighs; //two distinct returns are needed to avoid vector neighs to go out-of-scope
	}
	else
	{
		vector<vec_sz> neighs(q) ;
		unsigned int count=0;
		//given an adjacency matrix A , neighbors of node i are simply the non-zero entries of
		//the vector A[i]
		for(unsigned j=0; j<N && count<q; j++)
		{
			if(this->weight[index(i,j)]==1)
			{
				neighs[count]=j;
				count++;
			}
		}
		return neighs;
	}

}




