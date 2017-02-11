#include "IsingModel.h"
#include "Riostream.h"
#include "TRandom3.h"

int test_ising(int N)
{	//cout<<"Creating is1"<<endl;
	IsingModel *is1=new IsingModel(1);
	TRandom3* rd = new TRandom3(10*98921);
 	cout<<"Using newGraph(100) to change dimension:"<<endl;
 	is1->newGraph(10,"random");
 	cout<<"new graph"<<endl;
 	cout<<is1->getMagnetization()<<endl;
 	cout<<"magn"<<endl;
 	double en0=is1->hamiltonian();
 	cout<<"hamil"<<endl;
 	double en1=0;
 	for(unsigned i=0; i<N; i++) {
 		unsigned j = rd->Integer(10);
 		is1->flipSpin(j);
 		en1=is1->hamiltonian();
 		if(en0>en1) is1->flipSpin(j);
 		else en0=en1;
 		cout<<is1->getMagnetization()<<endl;
 	}
 	delete is1;
return 1; 
}

vec_sz idx(const vec_sz x, const vec_sz y,unsigned N)
{	//Return the correct index for 1d representation of 2d matrices; do not change data members
	//See http://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
	return  x*N + y;
}

int testing_simulate(int max_mcs ,unsigned long max_side_dim=10)
{	cout<<"[+]TESTING UTILITY"<<endl;
	IsingModel ising_model(max_side_dim);
	cout<<"[+]Starting the simulation procedure."<<endl;

	vector<unsigned> length_list={10};
	vec_sz list_size = length_list.size();
	cout<<"[+]Simulating 2D Ising model on a square lattice with L*L spins "<<endl;
	cout<<"   for the following values of L:"<<endl;
	for (unsigned i=0; i<list_size; i++) cout<<" "<<length_list[i]<<" ";
	cout<<endl;
	for (unsigned i =0; i<list_size; i++)
			if(length_list[i]>max_side_dim)
				cout<<"[!]Maximum lattice size ("<<max_side_dim<<") exceeded: results could be unreliable!"<<endl;
	double max_beta = 200.; // min temp T=5e-3
	double min_beta = 0.; 	//max temp = infinity
	vec_sz temp_steps=80;
	vector<double> inv_temperature(temp_steps);
	for (int i = 0; i<temp_steps; i++ ){
		inv_temperature[i]=min_beta+i*(max_beta-min_beta)/temp_steps;
	}
	//we used a 2x2 matrix representation in which binder_cumulant[index(i,j,N)]
	//represents the value of the bc at temp inv_temperature[j] for the i-th lattice in lenght_list
	//and N is the number of elements in a row
	vector<double> binder_cumulants(list_size*temp_steps,0);
	//vector<double> binder_cumulants(80,0);
	vec_sz bc_vec_sz = binder_cumulants.size();
	//same for the energy fluctuations vector
	vector<double> heat_capacity(list_size*temp_steps,0);
	vec_sz hc_vec_sz = heat_capacity.size();
	cout<<"[+]Run parameters correctly initialized. "<<endl;
	//unsigned long int mcs_i = 0;
	//unsigned long int max_steps = max_mcs* max_side_dim*max_side_dim;
	vector<double> energy(max_mcs,0);
	vector<double> magnetization(max_mcs,0);
	for(unsigned i=0; i< list_size; i++)
	{ 	cout<<"[+]Simulation for L="<<length_list[i]<<" started..."<<endl;
		ising_model.newGraph(length_list[i]);
		cout<<"[+]New graph created"<<endl;
		for(unsigned j=0; j < temp_steps ; j++)
		{	//mcs_i = max_mcs / (max_side_dim/length_list[i]); //integer division, this linear speed up "should not mess up"...
			//unsigned long int max_steps_i = mcs_i*length_list[i]*length_list[i]; //total n of steps
			//vectors of energy and magnetizations of the current model

			//burn-in time
			ising_model.simulate(inv_temperature[j],(max_mcs/10)+1); //empirical...
			for(unsigned k=0; k<max_mcs; k++)
			{
				ising_model.simulate(inv_temperature[j],length_list[i]*length_list[i]);
				energy[k]=ising_model.getEnergy();
				magnetization[k]=ising_model.getMagnetization();
				cout<<"[D] Current values of E,M ="<<energy[k]<<","<<magnetization[k]<<endl;
			}
			//procedure to compute binder's cumulant and heat capacity
			double m_2=0;
			double m_4=0;
			double m_2k=0;
			double e_2=0;
			double e_m=0;
			for (unsigned k=0; k<max_mcs; k++)
			{
				m_2k= magnetization[k]*magnetization[k];
				m_2+= m_2k;
				m_4+= m_2k*m_2k;
				e_2 = energy[k]*energy[k];
				e_m = energy[k];
			}
			m_4=m_4 / max_mcs;
			m_2=m_2 / max_mcs;
			e_2=e_2 / max_mcs;
			e_m=e_m / max_mcs;
			//cout<<"m_4 "<<m_4<<" "<<"m_2 "<<m_2<<" "<<" e_2 "<<e_2<<" e_m "<<e_m<<endl;
			cout<<"index i,j | real index | max index = "<<i<<","<<j<<"|"<<idx(i,j,temp_steps)<<"|"<<binder_cumulants.size()<<endl;
			binder_cumulants[idx(i,j, temp_steps)] = 1 - ( m_4/(3*m_2*m_2)); //compute the binder's cumulant for this T
			//heat_capacity[idx(i,j,temp_steps)] =inv_temperature[j]*inv_temperature[j]*(e_2-(e_m*e_m)); // compute the heat capacity(T)=k_b T^2 (var E)
			ising_model.resetGraph(); //reset each time!
		}
		cout<<"   [D]...completed."<<endl;
	}
	return -1;
}
