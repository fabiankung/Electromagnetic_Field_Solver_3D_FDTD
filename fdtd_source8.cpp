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


// Filename		 : fdtd_source8.cpp
// Remark		 : Code for simulating a PN junction or diode employing Ciampolini 
//                         and Piket-May approach
// Programming Language	 : C++
// Author		 : Fabian Kung Wai Lee


#include <math.h> // Standard ANSI C/C++ math library header
#include "fdtd_header.h"
#include "fdtd_class.h"
#include "fdtd_component_class.h"


///////////////////////////////////////////////////////////
// --- Definition of global variables --- /////////////////
///////////////////////////////////////////////////////////
const double eo = 8.854E-12; // Permittivity of free space
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
// --- Functions of Class tEDiodec --- /////////////////////
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
// plist = point to parameter linked-list, default = NULL.
tEDiodec::tEDiodec(double fdel_x, double fdel_y, double fdel_z,
				  double fdel_t, double fe, lpDOUBLE plist)
{
	double fC; // temporary variable

	AssignVariable(plist, fIs, 1.0E-10); // read Is
	if (plist == NULL) // read n
		fC = 1/38.6957; // Vt=kT/q, assume roome temparature of 25 degrees Celcius.
	else
	{
		fC = (plist->fP)/38.6957; // nVt
		plist = plist->next; // point to next item
	}
	fnVt = fC; //  nVt
	AssignVariable(plist,fVj,0.56); // read Vj
	AssignVariable(plist,fM,0.5);	// read M
	AssignVariable(plist,fCj0,0);	// read Cj0
	AssignVariable(plist,fTT,0);	// read TT
	AssignVariable(plist,fFC,0.5);	// read FC
	if (fFC > 1)					// limit the value of fFC
		fFC = 0.98;		
	AssignVariable(plist,fR,0);		// read R

	fC0 = fdel_t/(fe*eo*fdel_x);
	fC1 = fdel_t/(fe*eo*fdel_y);
	fC2 = fdel_t/(fe*eo*fdel_x*fdel_y);
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
// NOTE : This is an experimental routine, it is only intended for -ve oriented
//        diode.
double tEDiodec::updateEfield(double fE, double fH01, double fH00, double fH11,
						 double fH10)
{
	double ftemp1,ftemp2,ftemp3,fEnew;
	double fEold;
	double ftol = 1;
	int nindex;

	fEnew = fE ; // the 1st estimate of the new E field is the current value
	fEold = fE;
	ftemp3 = (fC0)*(fH11 - fH10) - (fC1)*(fH01 - fH00);
	for (nindex=0;nindex<5;nindex++)
	{
		// This is the function f(En)
		ftemp1 = fEnew - fE - ftemp3 + fC2*Idiode(fEnew,fE);
		// This is the function df(En)/dEn
		ftemp2 = 1+fC2*dIdiode(fEnew,fE);
		if (ftemp2 != 0) // prevent divide by 0
		{
			fEnew = fEnew - ftemp1/ftemp2;
			if (((fEnew-fEold)/fEold*100) < ftol) // check for convergence
				break;
			else
				fEold = fEnew;
		}
	}
	return fEnew;
}


// function to calculate the diode instantaneous current
// arguments:
// fEn = the estimated E field at E(n+1)
// fE = E field at E(n)
double tEDiodec::Idiode(double fEn, double fE)
{
	double fV,fC0i,fC1i,fC2i,ftemp1,fF2,fF3;
	double ftemp2;

	fV = -fE*fdel_z;
	fC0i = 1/fnVt;
	fC1i = 0.5*(fEn+fE)*fdel_z;
	fC2i = (fEn-fE)*fdel_z/fdel_t;

	if (fV < fFC*fVj)
	{
		ftemp1 = -fIs*(exp(-fC0i*fC1i)-1);
		ftemp2 = ftemp1 + fC0i*fTT*fIs*exp(-fC0i*fC1i)*fC2i;
		ftemp2 = ftemp2 + (fCj0/(pow(1+fC1i/fVj,fM)))*fC2i;
	}
	else
	{
		fF2 = pow(1-fFC,1+fM);
		fF3 = 1-fFC*(1+fM);
		ftemp1 = -fIs*(exp(-fC0i*fC1i)-1);
		ftemp2 = ftemp1 + fC0i*fTT*fIs*exp(-fC0i*fC1i)*fC2i;
		ftemp2 = ftemp2 + (fCj0/fF2)*(fF3-fM*fC1i/fVj)*fC2i;
	}
	
	return ftemp2;
}

// function to calculate the diode instantaneous current
// arguments:
// fEn = the estimated E field at E(n+1)
// fE = E field at E(n)
double tEDiodec::dIdiode(double fEn, double fE)
{	
	// approximate differentiation with finite difference.
	return (Idiode(fEn+0.0005,fE) - Idiode(fEn,fE))/0.0005;
}

// function to assign a value in the parameter list of the component to the coresponding
// variable in the object.
// arguments:
//	plist = reference to parameter list pointer
//	fP	  = reference to the variable (a double datatype)
//  fdefault = the
void tEDiodec::AssignVariable(lpDOUBLE &plist, double &fP, double fdefault)
{
	if (plist == NULL)	// if no valid data
		fP = fdefault;	// assign default value
	else
	{
		fP = plist->fP;
		plist = plist->next; // point to next item
	}
}