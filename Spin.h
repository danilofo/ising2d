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

	public:
	Spin(double=0);

	void setSpinUp();
	void setSpinDown();
	void flipSpin();

	const double getSpinValue();

	private:
	double spin;


};

#endif /* SPIN_H_ */
