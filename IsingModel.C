/*
 * IsingModel.C
 *
 *  Created on: 07 gen 2017
 *      Author: danilo
 */
#include "IsingModel.h"
#include "SimulationModel.h"


//implementazione della classe
ClassImp(IsingModel)

//Counter for instances of IsingModel
unsigned IsingModel::n_instances=0;
//public constructors
IsingModel::IsingModel(vec_sz length, double J): SimulationModel(),lattice(NULL),coupling(J)
{
	if(n_instances>0){
		cout<<"[!]IsingModel: only one instance of type IsingModel is allowed!"<<endl;
		return;
	}
	//initialize
	//cout<<"[D]IsingModel: constructing TRandom3 Rnd...";
	this->Rnd = new TRandom3(823); //TODO: fix the seed choice
	//cout<<"[D]IsingModel: Rnd...completed"<<endl;
	//cout<<"[D]IsingModel: constructing Lattice lattice...";
	this->lattice = new Lattice(length,Rnd);
	this->graph = this->lattice;
	if(Rnd==NULL || lattice==NULL || graph == NULL)
	{//Check the previous allocations
		cout<<"[!]IsingModel: Invalid memory allocation !"<<endl;
		return;
	}
	//cout<<"[D]IsingModel: lattice...completed"<<endl;
	//cout<<"[D]IsingModel: computing initial energy...";
	this->E=this->hamiltonian();
	//cout<<"[D]IsingModel: hamiltonian...completed"<<endl;
	vd_sz N = length*length;
	this->M=this->magnetization();
	n_instances+=1;
}

IsingModel::~IsingModel()
{
if(this->lattice!=NULL) delete lattice; 
IsingModel::n_instances=0;
}

//
double IsingModel::magnetization(){
	if(lattice==NULL){
		cout<<"[!]IsingModel:invalid lattice pointer received!"<<endl;
		return -100;
	}
	vd_sz N = this->lattice->getDimension();
	double result =0;
	for (vd_sz i=1; i<N; i++){
		result+=(this->lattice->getNodeValue(i)/N);
	}
	return result;
}


//functions needed to simulate
double IsingModel::hamiltonian()
{	if(lattice==NULL){
	cout<<"[!]IsingModel:invalid lattice pointer received!"<<endl;
	return 0;
	}
    double H=0;
    double deltaE = 0;
    const vd_sz N = this->lattice->getDimension();
    for(vd_sz i=0; i<N ; i++)
    {
        const vector<vd_sz> neighs_i=this->lattice->neighbors(i);
        const unsigned q = this->lattice->getCoordN() ;
        for(vd_sz j=0; j<q; j++)
        {
            //compute the hamiltonian; since the graph is not directed, each node appear twice
            //in the neighbors list; to speed up the computation we use only nodes with j>i instead
            //of dividing the result by 2
            if(neighs_i[j]>i) {
            	deltaE = (this->lattice->getNodeValue(i))*(this->lattice->getNodeValue(neighs_i[j]));
            	H+=  deltaE;
            }
        }
    }
    H=this->coupling*H; //multiply by the coupling constant
    return H;
}

double IsingModel::magnVar(double old_val, double new_val){
	vd_sz N = this->lattice->getDimension();
	return (-2.)*(old_val/N);
}

void IsingModel::resetGraph(){
	if(lattice==NULL){
		cout<<"[!]IsingModel:resetGraph: invalid lattice pointer"<<endl;
		return;
	}
	//cout<<"[D]resetGraph:old energy:"<<this->E<<endl;
	this->lattice->reset();
	//reset the values of E,M
	this->E = this->hamiltonian();
	this->M = this->magnetization();
	//cout<<"[D]newGraph:old energy:"<<this->E<<endl;

}
void IsingModel::newGraph(vec_sz N,const char* flag){
	//cout<<"[D]newGraph:old energy:"<<this->E<<endl;
	if(this->lattice!=NULL)
	{delete lattice; //prevent memory leak
	//cout<<"[D]IsingModel:lattice deleted"<<endl;
	this->lattice = new Lattice(N,this->Rnd,4,flag);
	if(lattice==NULL){
		//cout<<"[!]IsingModel:newGraph:Invalid allocation!"<<endl;
		return;
	}
	//cout<<"[D]IsingModel: new lattice created"<<endl;
	this->graph = this->lattice;
	this->E=this->hamiltonian();
	//cout<<"[D]IsingModel: hamiltonian...completed"<<endl;
	this->M = this->magnetization();
	//cout<<"[D]IsingModel: new lattice initializated"<<endl;
	//cout<<"[D]newGraph:new energy:"<<this->E<<endl;
	}
	else cout<<"[!]IsingModel:newGraph:Null pointer received"<<endl;
}



