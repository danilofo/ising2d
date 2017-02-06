/*
 * SimulationModel.h
 *
 *  Created on: 08 gen 2017
 *      Author: danilo
 */

#ifndef SIMULATIONMODEL_H_
#define SIMULATIONMODEL_H_
#include <algorithm>
#include "Graph.h"


typedef std::vector<int>::size_type vec_sz;

class SimulationModel: public TObject {
    // classe in cui definisco il grafo
    
    public:
    //public constructors
    void simulate(double beta, unsigned n_iterations=1000);
    void resetLattice(); //reset lattice function
    void newLattice(vec_sz new_length); // new lattice;
    
    
    const double getEnergy() const;
    const double getMagnetization() const;
    
    private:
    
    Graph graph;
    const double coupling;

    
    
    
};




#endif /* SIMULATIONMODEL_H_ */
