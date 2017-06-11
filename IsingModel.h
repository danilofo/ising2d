
/*
 * IsingModel.h
 *
 *  Created on: 07 gen 2017
 *      Author: danilo
 */
#if !defined(__CINT__) || defined(__MAKECINT__) //This directive is used by ROOT
#endif

#ifndef ISINGMODEL_H_
#define ISINGMODEL_H_
#include "Lattice.h"
#include "SimulationModel.h"

class IsingModel: public SimulationModel {
  //An Ising Model defined on a lattice (class Lattice)
  private:
  Lattice* lattice;
  const double coupling;
  
 public:
   IsingModel(vec_sz length, double J=-1.); //constructor
  ~IsingModel();
  //
  vec_sz getSideDimension();
  //simulation
  virtual double magnetization();
  virtual double hamiltonian(); //virtual qualifier here is ignored, but improves readability
  virtual double mVar(double old_val, double new_val);
  virtual double energyVar(vd_sz node_i, double old_val, double new_val);
  //actions on the lattice
  virtual void resetGraph(); //
  virtual void newGraph(vec_sz N,const char* flag="random"); 

 ClassDef(IsingModel,1); //Used by to define a class ROOT

};
#endif /* ISINGMODEL_H_ */
