/*
 * IsingModel.C
 *
 *  Created on: 07 gen 2017
 *      Author: danilo
 */
#include "IsingModel.h"
#include "SimulationModel.h"

ClassImp(IsingModel)

IsingModel::IsingModel(vec_sz length, double J): SimulationModel(),lattice(NULL),coupling(J)
{ //initialize
  //cout<<"[D]IsingModel: constructing TRandom3 Rnd...";
  this->Rnd = new TRandom3(); //Default choice is a UUID based on machine time
  //cout<<"[D]IsingModel: Rnd...completed"<<endl;
  //cout<<"[D]IsingModel: constructing Lattice lattice...";
  this->lattice = new Lattice(length,Rnd);
  this->graph = this->lattice;
  if(Rnd==NULL || lattice==NULL || graph == NULL)
    {//Check the previous allocations
      cout<<"[!]IsingModel: Invalid memory allocation !"<<endl;
      return;
    }
  this->E=this->hamiltonian();
  this->M=this->magnetization();
}

IsingModel::~IsingModel()
{
  if(this->lattice!=NULL) delete lattice; 
}

vd_sz IsingModel::getSideDimension(){
  const vd_sz&  raw_dim = this->lattice->getDimension();
  return static_cast<vd_sz>(sqrt(raw_dim));
}

double IsingModel::getSpin( vec_sz i ) const{
  return this->lattice->getNodeValue(i);
}
double IsingModel::magnetization(){
  if(lattice==NULL){
    cout<<"[!]IsingModel:invalid lattice pointer received!"<<endl;
    return -100;
  }
  const vd_sz& N = this->lattice->getDimension();
  double result=0;
  for (vd_sz i=0; i<N; i++) {

    result+=(this->lattice->getNodeValue(i))/N;
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
	    H+= deltaE;
	  }
        }
    }
  H=this->coupling*H; //multiply by the coupling constant
  return H;
}

double IsingModel::energyVar(vd_sz node_i, double old_val, double new_val){
  double result =0;
  vector<vd_sz> neighs_i= this->lattice->neighbors(node_i);
  unsigned q = neighs_i.size();
  // each product in the hamiltonian changes sign once-->  (-2)*old_i
  for(vd_sz j=0; j<q; j++){
    result-= old_val*this->lattice->getNodeValue(neighs_i[j]);
    result+= new_val*this->lattice->getNodeValue(neighs_i[j]);
  }
  result = this->coupling*result;
  return result;
}

double IsingModel::mVar(double old_val, double new_val){
  const vd_sz &N = this->lattice->getDimension();
  return (-2.0)*(old_val)/N;
}

void IsingModel::resetGraph(const char* mod){
  if(lattice==NULL){
    cout<<"[!]IsingModel:resetGraph: invalid lattice pointer"<<endl;
    return;
  }
  //cout<<"[D]resetGraph:old energy:"<<this->E<<endl;
  this->lattice->reset();
  //reset the values of E,M
  this->E = this->hamiltonian();
  this->M = this->magnetization();
}
void IsingModel::newGraph(vec_sz N,const char* flag){
  if(this->lattice!=NULL)
    {delete lattice; //prevent memory leak
      this->lattice = new Lattice(N,this->Rnd,4,flag);
      if(lattice==NULL){
	cout<<"[!]IsingModel:newGraph:Invalid allocation!"<<endl;
	return;
      }
      this->graph = this->lattice;
      this->E=this->hamiltonian();
      this->M=this->magnetization();
    }
  else cout<<"[!]IsingModel:newGraph:Null pointer received"<<endl;
}



