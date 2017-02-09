/*
 * Spin.C
 *
 *  Created on: 08 gen 2017
 *      Author: danilo
 */
#include "Spin.h"
#include "iostream"


Spin::Spin(double s){

	if(s==0) this->spin=0;
	else if(s==1)  this->setSpinUp();

	else if(s==-1) this->setSpinDown();
	else std::cout<<"[!]Spin:"<<s<<"is an invalid value for a spin variable"<<std::endl;
}

void Spin::setSpinUp(){
	this->spin=1;
}
void Spin::setSpinDown(){
	this->spin=-1;
}
void Spin::flipSpin(){
	this->spin-= 2*this->spin;
}

const double Spin::getSpinValue(){
	return this->spin;
}
