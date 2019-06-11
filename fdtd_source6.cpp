//    PCB FDTD Simulation Engine
//    Copyright (C) 2003 Fabian Kung Wai Lee
//
//    This file is part of PCB FDTD Simulation Engine.
//
//    PCB FDTD Simulation Engine is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    PCB FDTD Simulation Engine is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with PCB FDTD Simulation Engine.  
//    If not, see <http://www.gnu.org/licenses/>.


// Filename		 : fdtd_source6.cpp 
// Remark		 : Code for simulating a PN junction or diode
// Programming Language  : C++
// Programmer	 : Fabian Kung Wai Lee


#include <math.h>  // Standard ANSI C/C++ math library header
#include "fdtd_header.h"
#include "fdtd_class.h"
#include "fdtd_component_class.h"


///////////////////////////////////////////////////////////
// --- Definition of global variables --- /////////////////
///////////////////////////////////////////////////////////
const double eo = 8.85418E-12; // Permittivity of free space
const double uo = 4.0E-7*3.1415926535898; // Permeability of free space
const double c = 2.998E8;	// Speed of light in free space

///////////////////////////////////////////////////////////
// --- Definition of external variables --- ///////////////
///////////////////////////////////////////////////////////
extern double fdel_x;
extern double fdel_y;
extern double fdel_z;
extern double fdel_t;

///////////////////////////////////////////////////////////
// --- Functions of Class tEDiode --- /////////////////////
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
// --- Default Constructor --- ////////////////////////////
///////////////////////////////////////////////////////////
// The constructor performs the following basic functions:
// (1) Initialize the diode model parameters
// (2) precompute the computation coefficients
// arguments:
// fdel_x = discretization in x, default = 1.0mm.
// fdel_y = discretization in y, default = 1.0mm.
// fdel_z = discretization in z, default = 1.0mm.
// fdel_t = discretization of time, default = 1.0psec.
// fe = dielectric constant, default = 1.0.
// ucor = orientation specification. 0=positive x, 1=negative x
//									 2=positive y, 3=negative y
//									 4=positive z, otherwise=negative z.
// plist = point to parameter linked-list, default = NULL.
tEDiode::tEDiode(double fdel_x, double fdel_y, double fdel_z,
				 double fdel_t, double fe, unsigned char ucor, lpDOUBLE plist)
{
	double fC; // temporary variable

	ucxyz = ucor; // set orientation of device
	if ((ucor==0) || (ucor==2) || (ucor==4))
		fn = 1;	// positive orientation
	else
		fn = -1;	// negative orientation
	if (plist != NULL) // skip O
		plist = plist->next; 
	AssignVariable(plist, fIs, 1.0E-10); // read Is
	if (plist == NULL) // read n
		fC = 1/38.6957; // Vt=kT/q, assume room temparature of 25 degrees Celcius.
	else
	{
		fC = (plist->fP)/38.6957; // nVt
		plist = plist->next; // point to next item
	}
	fnVt = fC; //  nVt
	AssignVariable(plist,fVj,0.6); // read Vj
	AssignVariable(plist,fM,0.5);	// read M
	AssignVariable(plist,fCj0,0);	// read Cj0
	AssignVariable(plist,fTT,0);	// read TT
	AssignVariable(plist,fFC,0.5);	// read FC
	if (fFC > 1)					// limit the value of fFC
		fFC = 0.98;		
	AssignVariable(plist,fR,0);		// read R

	if ((ucxyz==0) || (ucxyz==1)) // x oriented diode
	{
		fC0 = fdel_t/(fe*eo*fdel_y); // C0x
		fC1 = fdel_t/(fe*eo*fdel_z); // C1x
		fC2 = (fn*fdel_t)/(fe*eo*fdel_z*fdel_y); // C2x
		fC3 = fC2*fdel_x/fn; // C3x
		fC4 = (fn*fC3*fdel_x)/2; // C4x
		fdel_d = fdel_x;
	}
	else if ((ucxyz==2) || (ucxyz==3)) // y oriented diode
	{
		fC0 = fdel_t/(fe*eo*fdel_z); // C0y
		fC1 = fdel_t/(fe*eo*fdel_x); // C1y
		fC2 = (fn*fdel_t)/(fe*eo*fdel_x*fdel_z); // C2y
		fC3 = fC2*fdel_y/fn; // C3y
		fC4 = (fn*fC3*fdel_y)/2; // C4y
		fdel_d = fdel_y;
	}
	else // z oriented diode
	{
		fC0 = fdel_t/(fe*eo*fdel_x); // C0z
		fC1 = fdel_t/(fe*eo*fdel_y); // C1z
		fC2 = (fn*fdel_t)/(fe*eo*fdel_x*fdel_y); // C2z
		fC3 = fC2*fdel_z/fn; // C3z
		fC4 = (fn*fC3*fdel_z)/2; // C4z
		fdel_d = fdel_z;
	}

	fC5 = 1/fnVt;	// C5x, C5y or C5z.
	fC6 = 1/fdel_t; // C6x, C6y or C6z.
}

///////////////////////////////////////////////////////////
// --- Public functions --- ///////////////////////////////
///////////////////////////////////////////////////////////

// function to update the electric field component within the PN junction
// arguments:
// fE = current E field value
// fH01 = adjacent H field component
// fH00 = adjacent H field component
// fH11 = adjacent H field component
// fH10 = adjacent H field component
double tEDiode::updateEfield(double fE, double fH01, double fH00, double fH11,
						 double fH10)
{
	double fA, fB, fC, fV;

	fV = fn*fE*fdel_d;
	fdi = PN_junc_di(fV);
	ftemp2 = (fC0)*(fH11 - fH10) - (fC1)*(fH01 - fH00);
	fA = fC4*(0.25*PN_junc_ddi(fV)+fC6*PN_junc_dcap(fV));
	fB = 1.0 + fC3*(0.5*fdi+fC6*PN_junc_cap(fV));
	fC = fC2*PN_junc_i(fV)-ftemp2;
	ftemp2 = (fB*fB)-4*fA*fC;
	if (fA == 0) // to prevent divide by zero error
	{
		ftemp1 = -fC/fB;
	}
	else
	{	
		if (ftemp2 >= 0) // prevent square root of imaginery number
		{
			ftemp1 = (-fB+sqrt(ftemp2))/(2*fA);
		}
		else // if ftemp2 is zero, the following steps 1.) ignore higher differentiation of Id
			 // 2.) if step 1 is not succesfully, 1st order solution is used.
		{
			if ((ucxyz==0)||(ucxyz==2)||(ucxyz==4)) // check orientation, whether 
													// +ve or -ve direction.
			{
				fA = fC4*(fC6*PN_junc_dcap(fE*fdel_d));
			}
			else
			{
				fA = fC4*(-fC6*PN_junc_dcap(-fE*fdel_d));
			}

			ftemp2 = (fB*fB)-4*fA*fC;
			if (fA == 0)
			{
				ftemp1 = -fC/fB;
			}
			else
			{
				ftemp1 = (-fB+sqrt(ftemp2))/(2*fA);
			}
		}
	}
	
	// Sometime due to the finite precision of the computer where when a large number
	// and a very small number (<< 1) are added together, the smaller number is 
	// truncated and the sum equal the large number. This can result in -B+sqrt(B**2-4AC)
	// =0.  And the electric across the diode artificially locked to a fixed value.
	// To prevent such condition, the 1st order solution is used.
	if (ftemp1 == 0)
		ftemp1 = -fC/fB;
		
	// update E field
	fE = fE + ftemp1;

	return fE;
}

// function to calculate the total dynamic junction capacitance of the PN junction
// arguments: fV = biasing voltage
// return value: the capacitance
double tEDiode::PN_junc_cap(double fV)
{
	double ftemp; // temporary variable
	double fF2,fF3;

	if (fV < (fFC*fVj))
	{
		ftemp = fCj0/pow((1-fV/fVj),fM);
		return (ftemp + ((fTT*fIs)/fnVt)*exp(fV/fnVt));
	}
	else
	{
		ftemp = ((fTT*fIs)/fnVt)*exp(fV/fnVt);
		fF2 = pow(1-fFC,1+fM);
		fF3 = 1-fFC*(1+fM);
		return (ftemp + (fCj0/fF2)*(fF3+fM*fV/fVj));
	}
}

// function to calculate the differentiation of the total dynamic junction 
// capacitance of the PN junction
// arguments: fV = biasing voltage
// return value: the capacitance
double tEDiode::PN_junc_dcap(double fV)
{
	double ftemp1; // temporary variable
	double ftemp2;
	double fF2;

	if (fV < (fFC*fVj))
	{
		ftemp1 = (fM/fVj)*fCj0/pow((1-fV/fVj),fM+1);
		ftemp2 = ((fTT*fIs)*exp(fV/fnVt))/(fnVt*fnVt);
		return ftemp1+ftemp2;
	}
	else
	{
		fF2 = pow(1-fFC,1+fM);
		ftemp1 = exp(fV/fnVt);
		ftemp1 = ftemp1/(fnVt*fnVt);
		return (fTT*fIs)*ftemp1 + ((fCj0*fM)/(fF2*fVj));
	}
}

// function to calculate the static PN junction current as a function of biasing 
// voltage
// argument: fV = bias voltage
// return: static diode current Id
double tEDiode::PN_junc_i(double fV)
{
	return fIs*(exp(fV*fC5)-1);
}

// function to calculate the differentiation of static PN junction current as a 
// function of biasing voltage
// argument: fV = bias voltage
// return: dI/dV
double tEDiode::PN_junc_di(double fV)
{
	return fIs*fC5*exp(fV*fC5);
}

// function to calculate the double differentiation of static PN junction current  
// as a function of biasing voltage
// argument: fV = bias voltage
// return: (d/dV)(d/dV)I
double tEDiode::PN_junc_ddi(double fV)
{
	return fIs*fC5*fC5*exp(fV*fC5);
}

// function to assign a value in the parameter list of the component to the coresponding
// variable in the object.
// arguments:
//	plist = reference to parameter list pointer
//	fP	  = reference to the variable (a double datatype)
//  fdefault = the
void tEDiode::AssignVariable(lpDOUBLE &plist, double &fP, double fdefault)
{
	if (plist == NULL)	// if no valid data
		fP = fdefault;	// assign default value
	else
	{
		fP = plist->fP;
		plist = plist->next; // point to next item
	}
}