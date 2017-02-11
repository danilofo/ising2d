/*
 * Spin.h
 *
 *  Created on: 08 gen 2017
 *      Author: danilo
 */
#if !defined(__CINT__) || defined(__MAKECINT__) //This directive is used by ROOT
#endif

#ifndef SPIN_H_
#define SPIN_H_

class Spin{
//A basic data structure for spin variables
//
	public:
	Spin(double=0); //default value at construction time

	//public methods
	void setSpinUp();
	void setSpinDown();
	void flipSpin();
	double getSpinValue();

	private:
	double spin;


};

#endif /* SPIN_H_ */
