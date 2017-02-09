#include "Lattice.h"

//implementazione della classe
ClassImp(Lattice)

//constructors


Lattice::Lattice(int N, TRandom3* rnd,int q, const char *flag):
Graph(N*N,rnd),lenght(N),coord_number(q)
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

const unsigned Lattice::getCoordN(){
	return this->coord_number;
}

//initialization procedures:
void Lattice::initSpins(const char * flag)
//Initialize the value of each spin ; choice allowed are
//"random" or "ordered". The latter means that all spins are assigned to +1
{	const vec_sz N = this->getDimension();
	if (strncmp(flag,"random",10)==0)
	{
		for(vec_sz i=0; i<N; i++)
		{
			if(Rnd->Rndm() < 0.5) this->spin[i].setSpinUp();
			else this->spin[i].setSpinDown();
		}
	}

	else if (strncmp(flag,"ordered",10)==0)
	{	for(vec_sz i=0; i<N; i++)
		{
			this->spin[i].setSpinUp();
		}
	}
	else
	{
		cout<<"[!]Initialization failed. Admitted choices : \"random\", \"ordered\" "<<endl;
	}

}

void Lattice::reset(const char* flag)
{	const vd_sz N = this->getDimension();
	if (strncmp(flag,"random",10)==0)
	{
		for(vec_sz i=0; i<N; i++)
				{
					if(Rnd->Rndm() < 0.5) this->spin[i].setSpinUp();
					else this->spin[i].setSpinDown();
				}
	}

	else if (strncmp(flag,"ordered",10)==0)
	{	for(vec_sz i=0; i<N; i++)
	{
		this->spin[i].setSpinUp();
	}
	}
	else
	{
		cout<<"[!]Initialization failed. Admitted choices : \"random\", \"ordered\" "<<endl;
	}
}



void Lattice::initRectangularLattice()
{

	const vd_sz N = this->getDimension();
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
			if( (i+1)%L ==0)  //end of row
			{
				this->weight[Lattice::index((i+1-L) ,i)] = 1;
				this->weight[Lattice::index(i,(i+1-L))] = 1;
				this->weight[Lattice::index((i+L)%N,i)] = 1;
				this->weight[Lattice::index(i,(i+L)%N)] = 1;
			}
			else{
				this->weight[Lattice::index((i+1),i)] = 1;
				this->weight[Lattice::index(i,(i+1))] = 1;
				this->weight[Lattice::index((i+L)%N,i)] = 1;
				this->weight[Lattice::index(i,(i+L)%N)] = 1;
			}

		}
	}
	else{cout<<"[!]Re-initialization attempt ignored"<<endl;}
}


//class functionalities:
void Lattice::flip(const vec_sz i)
{
	const vd_sz N = this->getDimension();
	if(i>N-1)
	{
		cout<<"Node "<<i<<" is not in the lattice (dimension="<<N<<")"<<endl;
	}
	else{
		this->spin[i].flipSpin();
	}
}

const vector<vec_sz> Lattice::neighbors( const vec_sz i)
{
	const unsigned int &q = this->coord_number;
	const vd_sz N = getDimension();

	if (i>N-1) //check if i is in the lattice
	{
		cout<<"[!]Node "<<i<<" is not in the lattice (dimension="<<N<<")"<<endl;
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



