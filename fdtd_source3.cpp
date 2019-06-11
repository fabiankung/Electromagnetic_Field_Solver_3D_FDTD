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


// Filename		 : fdtd_source3.cpp
// Remark		 : Functions for class tEHfdtd
// Programming Language  : C++
// Author	 : Fabian Kung Wai Lee

#include <math.h>	// Standard ANSI C/C++ math library header
#include "fdtd_header.h"
#include "fdtd_class.h"
#include "fdtd_component_class.h"

///////////////////////////////////////////////////////////
// --- Definition of external variables --- ///////////////
///////////////////////////////////////////////////////////
extern double fdel_x;
extern double fdel_y;
extern double fdel_z;
extern double fdel_t;

///////////////////////////////////////////////////////////
// --- Functions of Class tEHfdtd --- /////////////////////
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
// --- Default Constructor --- ////////////////////////////
///////////////////////////////////////////////////////////
// The constructor performs 1 basic function:
// (1) Invoke the base class constructor.
tEHfdtd::tEHfdtd() : tEfield(), tHfield()
{
}


///////////////////////////////////////////////////////////
// --- Public functions --- ///////////////////////////////
///////////////////////////////////////////////////////////

// function to update the Ex electric field component
// umode = 0: default update mode.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEx_type0(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamex[i][j][k];
	fEx[i][j][k] = fEx[i][j][k] +
				   (stCCex[uparam].fC0)*(fHz[i][j][k] - fHz[i][j-1][k]) -
				   (stCCex[uparam].fC1)*(fHy[i][j][k] - fHy[i][j][k-1]);
}

// function to update the Ey electric field component
// umode = 0: default update mode.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEy_type0(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamey[i][j][k];
	fEy[i][j][k] = fEy[i][j][k] +
				   (stCCey[uparam].fC0)*(fHx[i][j][k] - fHx[i][j][k-1]) -
				   (stCCey[uparam].fC1)*(fHz[i][j][k] - fHz[i-1][j][k]);
}

// function to update the Ez electric field component
// umode = 0: default update mode.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEz_type0(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamez[i][j][k];
	fEz[i][j][k] = fEz[i][j][k] +
				   (stCCez[uparam].fC0)*(fHy[i][j][k] - fHy[i-1][j][k]) -
				   (stCCez[uparam].fC1)*(fHx[i][j][k] - fHx[i][j-1][k]);
}

// function to update the Ex electric field component
// umode = 1: field in PEC.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEx_type1(int i, int j, int k)
{
	fEx[i][j][k] = 0;
}

// function to update the Ey electric field component
// umode = 1: field in PEC.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEy_type1(int i, int j, int k)
{
	fEy[i][j][k] = 0;
}

// function to update the Ez electric field component
// umode = 1: field in PEC.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEz_type1(int i, int j, int k)
{
	fEz[i][j][k] = 0;
}

// function to update the Ex electric field component
// umode = 2: lossy dielectric update mode.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEx_type2(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamex[i][j][k];
	fEx[i][j][k] = (stCCex[uparam].fC0)*fEx[i][j][k] +
				   (stCCex[uparam].fC1)*(fHz[i][j][k] - fHz[i][j-1][k]) -
				   (stCCex[uparam].fC2)*(fHy[i][j][k] - fHy[i][j][k-1]);
}

// function to update the Ey electric field component
// umode = 2: lossy dielectric update mode.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEy_type2(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamey[i][j][k];
	fEy[i][j][k] = (stCCey[uparam].fC0)*fEy[i][j][k] +
				   (stCCey[uparam].fC1)*(fHx[i][j][k] - fHx[i][j][k-1]) -
				   (stCCey[uparam].fC2)*(fHz[i][j][k] - fHz[i-1][j][k]);
}

// function to update the Ez electric field component
// umode = 2: lossy dielectric update mode.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEz_type2(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamez[i][j][k];
	fEz[i][j][k] = (stCCez[uparam].fC0)*fEz[i][j][k] +
				   (stCCez[uparam].fC1)*(fHy[i][j][k] - fHy[i-1][j][k]) -
				   (stCCez[uparam].fC2)*(fHx[i][j][k] - fHx[i][j-1][k]);
}

// function to update the Ex electric field component
// umode = 50: Resistive voltage source, trapezoidal.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// t = current time step
// return: none
void tEHfdtd::UpdateEx_type50(int i, int j, int k, double t)
{
	unsigned char uparam;
	double fdelay, frise, fhigh, flow, fperiod, fsource, fVlow, fVhigh;

	uparam = uparamex[i][j][k];
	
	fdelay = stCCex[uparam].fC6;
	frise = stCCex[uparam].fC7;
	fhigh = stCCex[uparam].fC8;
	flow = stCCex[uparam].fC9;
	fVlow = stCCex[uparam].fC4;
	fVhigh = stCCex[uparam].fC5;

	// compute source output
	t = t-fdelay;
	fperiod = 2*frise+fhigh+flow;

	while (t>=fperiod)
	{
		t = t-fperiod;
	}

	if (t<0)
		fsource = 0;
	else if (t<frise)
		fsource = t/frise;
	else if ((t-frise)<fhigh)
		fsource = 1;
	else if ((t-frise-fhigh)<frise)
		fsource = 1 - (t-frise-fhigh)/frise;
	else 
		fsource = 0;

	fsource = fsource*(fVhigh-fVlow) + fVlow;
	fEx[i][j][k] = (stCCex[uparam].fC0)*fEx[i][j][k] +
				   (stCCex[uparam].fC1)*(fHz[i][j][k] - fHz[i][j-1][k]) -
				   (stCCex[uparam].fC2)*(fHy[i][j][k] - fHy[i][j][k-1]) -
				   (stCCex[uparam].fC3)*(fsource);
}

// function to update the Ey electric field component
// umode = 50: Resistive voltage source, trapezoidal.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// t = current time step
// return: none
void tEHfdtd::UpdateEy_type50(int i, int j, int k, double t)
{
	unsigned char uparam;
	double fdelay, frise, fhigh, flow, fperiod, fsource, fVlow, fVhigh;

	uparam = uparamey[i][j][k];
	
	fdelay = stCCey[uparam].fC6;
	frise = stCCey[uparam].fC7;
	fhigh = stCCey[uparam].fC8;
	flow = stCCey[uparam].fC9;
	fVlow = stCCey[uparam].fC4;
	fVhigh = stCCey[uparam].fC5;

	// compute source output
	t = t-fdelay;
	fperiod = 2*frise+fhigh+flow;

	while (t>=fperiod)
	{
		t = t-fperiod;
	}

	if (t<0)
		fsource = 0;
	else if (t<frise)
		fsource = t/frise;
	else if ((t-frise)<fhigh)
		fsource = 1;
	else if ((t-frise-fhigh)<frise)
		fsource = 1 - (t-frise-fhigh)/frise;
	else 
		fsource = 0;

	fsource = fsource*(fVhigh-fVlow) + fVlow;
	fEy[i][j][k] = (stCCey[uparam].fC0)*fEy[i][j][k] +
				   (stCCey[uparam].fC1)*(fHx[i][j][k] - fHx[i][j][k-1]) -
				   (stCCey[uparam].fC2)*(fHz[i][j][k] - fHz[i-1][j][k]) -
				   (stCCey[uparam].fC3)*(fsource);
}

// function to update the Ez electric field component
// umode = 50: Resistive voltage source, trapezoidal.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// t = current time step
// return: none
void tEHfdtd::UpdateEz_type50(int i, int j, int k, double t)
{
	unsigned char uparam;
	double fdelay, frise, fhigh, flow, fperiod, fsource, fVlow, fVhigh;
	double fyted;

	uparam = uparamez[i][j][k];
	
	fdelay = stCCez[uparam].fC6;
	frise = stCCez[uparam].fC7;
	fhigh = stCCez[uparam].fC8;
	flow = stCCez[uparam].fC9;
	fVlow = stCCez[uparam].fC4;
	fVhigh = stCCez[uparam].fC5;

	// compute source output
	t = t-fdelay;
	fperiod = 2*frise+fhigh+flow;

	while (t>=fperiod)
	{
		t = t-fperiod;
	}

	if (t<0)
		fsource = 0;
	else if (t<frise)
		fsource = t/frise;
	else if ((t-frise)<fhigh)
		fsource = 1;
	else if ((t-frise-fhigh)<frise)
		fsource = 1 - (t-frise-fhigh)/frise;
	else 
		fsource = 0;

	fsource = fsource*(fVhigh-fVlow) + fVlow;

	fyted = (stCCez[uparam].fC1)*(fHy[i][j][k] - fHy[i-1][j][k]) -
		   (stCCez[uparam].fC2)*(fHx[i][j][k] - fHx[i][j-1][k]);

//	if (fsource > 0) // limit the magnitude of curl(H), so that the source is proper
//	{
//		if (fyted > (2*fVhigh/fdel_z))
//			fyted = (2*fVhigh/fdel_z);
//	}
//	else if (fsource < 0) 
//	{ 
//		if (fyted < (2*fVlow/fdel_z))
//			fyted = (2*fVlow/fdel_z); 
//	}

	fEz[i][j][k] = (stCCez[uparam].fC0)*fEz[i][j][k] + fyted -
				   (stCCez[uparam].fC3)*(fsource);
}

// function to update the Ex electric field component
// umode = 51: Resistive voltage source, sinusoidal.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// t = current time step
// nSimMode = Simulation mode, to differentiate between REUSE and NORMAL mode.
// In REUSE mode, the field value at t=0 will be added to the previous E field.
// In NORMAL mode, the previous E field at t=0 is set to 0.
// return: none
void tEHfdtd::UpdateEx_type51(int i, int j, int k, double t, int SimMode)
{
	unsigned char uparam;
	double fphase, ffreq, fatten, famp, fsource;		
	double fEx_prev;

	uparam = uparamex[i][j][k];
	
	if (t == 0)
	{
		if (SimMode != REUSE)
			fEx_prev = 0;
		else
			fEx_prev = fEx[i][j][k];
	}
	else
		fEx_prev = 0;
	fphase = stCCex[uparam].fC4;
	ffreq = stCCex[uparam].fC5;
	fatten = stCCex[uparam].fC6;
	famp = stCCex[uparam].fC7;

	// compute source output
	fsource = exp(-fatten*t)*sin(6.283185307*ffreq*t-fphase);


	fEx[i][j][k] = (stCCex[uparam].fC0)*fEx[i][j][k] +
				   (stCCex[uparam].fC1)*(fHz[i][j][k] - fHz[i][j-1][k]) -
				   (stCCex[uparam].fC2)*(fHy[i][j][k] - fHy[i][j][k-1]) -
				   (stCCex[uparam].fC3)*(fsource*famp)+
				   fEx_prev; // this is the previous field
}

// function to update the Ey electric field component
// umode = 51: Resistive voltage source, sinusoidal.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// t = current time step
// nSimMode = Simulation mode, to differentiate between REUSE and NORMAL mode.
// In REUSE mode, the field value at t=0 will be added to the previous E field.
// In NORMAL mode, the previous E field at t=0 is set to 0.
// return: none
void tEHfdtd::UpdateEy_type51(int i, int j, int k, double t, int SimMode)
{
	unsigned char uparam;
	double fphase, ffreq, fatten, famp, fsource;		
	double fEy_prev;

	uparam = uparamey[i][j][k];
	
	if (t == 0)
	{
		if (SimMode != REUSE)
			fEy_prev = 0;
		else
			fEy_prev = fEy[i][j][k];
	}
	else
		fEy_prev = 0;
	fphase = stCCey[uparam].fC4;
	ffreq = stCCey[uparam].fC5;
	fatten = stCCey[uparam].fC6;
	famp = stCCey[uparam].fC7;

	// compute source output
	fsource = exp(-fatten*t)*sin(6.283185307*ffreq*t-fphase);


	fEy[i][j][k] = (stCCez[uparam].fC0)*fEy[i][j][k] +
				   (stCCez[uparam].fC1)*(fHx[i][j][k] - fHx[i][j][k-1]) -
				   (stCCez[uparam].fC2)*(fHz[i][j][k] - fHz[i-1][j][k]) -
				   (stCCez[uparam].fC3)*(fsource*famp)+
				   fEy_prev; // this is the previous field
}

// function to update the Ez electric field component
// umode = 51: Resistive voltage source, sinusoidal.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// t = current time step
// nSimMode = Simulation mode, to differentiate between REUSE and NORMAL mode.
// In REUSE mode, the field value at t=0 will be added to the previous E field.
// In NORMAL mode, the previous E field at t=0 is set to 0.
// return: none
void tEHfdtd::UpdateEz_type51(int i, int j, int k, double t, int SimMode)
{
	unsigned char uparam;
	double fphase, ffreq, fatten, famp, fsource;		
	double fEz_prev, fyted;

	uparam = uparamez[i][j][k];
	
	if (t == 0)
	{
		if (SimMode != REUSE)
			fEz_prev = 0;
		else
			fEz_prev = fEz[i][j][k];
	}
	else
		fEz_prev = 0;
	fphase = stCCez[uparam].fC4;
	ffreq = stCCez[uparam].fC5;
	fatten = stCCez[uparam].fC6;
	famp = stCCez[uparam].fC7;

	// compute source output
	fsource = exp(-fatten*t)*sin(6.283185307*ffreq*t-fphase);

	fyted = (stCCez[uparam].fC1)*(fHy[i][j][k] - fHy[i-1][j][k]) -
		   (stCCez[uparam].fC2)*(fHx[i][j][k] - fHx[i][j-1][k]);

//	if (fsource > 0) // limit the magnitude of curl(H), so that the source is proper
//	{
//		if (fyted > (2*famp/fdel_z))
//			fyted = (2*famp/fdel_z);
//	}
//	else if (fsource < 0) 
//	{ 
//		if (fyted < (-2*famp/fdel_z))
//			fyted = (-2*famp/fdel_z); 
//	}

	fEz[i][j][k] = (stCCez[uparam].fC0)*fEz[i][j][k] + fyted -
				   (stCCez[uparam].fC3)*(fsource*famp)+
				   fEz_prev; // this is the previous field
}

// function to update the Ex electric field component
// umode = 52: Resistive dc voltage source.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// t = current time step
// return: none
void tEHfdtd::UpdateEx_type52(int i, int j, int k, double t)
{
	unsigned char uparam;
	double fVdc;

	uparam = uparamex[i][j][k];
	
	fVdc = stCCex[uparam].fC4;

	if (t==0) // we want to set the initialized value to the required dc voltage.
		fEx[i][j][k] = -fVdc/fdel_x;

	fEx[i][j][k] = (stCCex[uparam].fC0)*fEx[i][j][k] +
				   (stCCex[uparam].fC1)*(fHz[i][j][k] - fHz[i][j-1][k]) -
				   (stCCex[uparam].fC2)*(fHy[i][j][k] - fHy[i][j][k-1]) -
				   (stCCex[uparam].fC3)*fVdc;
}

// function to update the Ey electric field component
// umode = 52: Resistive dc voltage source.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// t = current time step
// return: none
void tEHfdtd::UpdateEy_type52(int i, int j, int k, double t)
{
	unsigned char uparam;
	double fVdc;

	uparam = uparamey[i][j][k];
	
	fVdc = stCCey[uparam].fC4;

	if (t==0) // we want to set the initialized value to the required dc voltage.
		fEy[i][j][k] = -fVdc/fdel_y;

	fEy[i][j][k] = (stCCey[uparam].fC0)*fEy[i][j][k] +
				   (stCCey[uparam].fC1)*(fHx[i][j][k] - fHx[i][j][k-1]) -
				   (stCCey[uparam].fC2)*(fHz[i][j][k] - fHz[i-1][j][k]) -
				   (stCCey[uparam].fC3)*fVdc;
}

// function to update the Ez electric field component
// umode = 52: Resistive dc voltage source.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// t = current time step
// return: none
void tEHfdtd::UpdateEz_type52(int i, int j, int k, double t)
{
	unsigned char uparam;
	double fVdc, fyted;

	uparam = uparamez[i][j][k];
	
	fVdc = stCCez[uparam].fC4;

	if (t==0) // we want to set the initialized value to the required dc voltage.
		fEz[i][j][k] = -fVdc/fdel_z;

	fyted = (stCCez[uparam].fC1)*(fHy[i][j][k] - fHy[i-1][j][k]) -
		   (stCCez[uparam].fC2)*(fHx[i][j][k] - fHx[i][j-1][k]);

//	if (fVdc > 0) // limit the magnitude of curl(H), so that the source is proper
//	{
//		if (fyted > (2*fVdc/fdel_z))
//			fyted = (2*fVdc/fdel_z);
//	}
//	else if (fVdc < 0) 
//	{ 
//		if (fyted < (-2*fVdc/fdel_z))
//			fyted = (-2*fVdc/fdel_z); 
//	}

	fEz[i][j][k] = (stCCez[uparam].fC0)*fEz[i][j][k] + fyted -
				   (stCCez[uparam].fC3)*fVdc;
}

// function to update the Ex electric field component
// umode = 100: Resistor.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEx_type100(int i, int j, int k)
{
	unsigned char uparam;
	double fEold, fpower, fJ;

	uparam = uparamex[i][j][k];
	fEold = fEx[i][j][k];
	fJ = (stCCex[uparam].fC3)*fEold;
	
	fEx[i][j][k] = (stCCex[uparam].fC0)*fEx[i][j][k] +
				   (stCCex[uparam].fC1)*(fHz[i][j][k] - fHz[i][j-1][k]) -
				   (stCCex[uparam].fC2)*(fHy[i][j][k] - fHy[i][j][k-1]);

	fpower = 0.5*(fEx[i][j][k]+fEold)*fJ; // compute numerical power 
	if (fpower < 0)	// perform compensation if numerical power is -ve
	{
		fEx[i][j][k] = -fEold;
	}
}

// function to update the Ey electric field component
// umode = 100: Resistor.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEy_type100(int i, int j, int k)
{
	unsigned char uparam;
	double fEold, fpower, fJ;

	uparam = uparamey[i][j][k];
	fEold = fEy[i][j][k];
	fJ = (stCCey[uparam].fC3)*fEold;

	fEy[i][j][k] = (stCCey[uparam].fC0)*fEy[i][j][k] +
				   (stCCey[uparam].fC1)*(fHx[i][j][k] - fHx[i][j][k-1]) -
				   (stCCey[uparam].fC2)*(fHz[i][j][k] - fHz[i-1][j][k]);

	fpower = 0.5*(fEy[i][j][k]+fEold)*fJ;  // compute numerical power
	if (fpower < 0)	// perform compensation if numerical power is -ve
	{
		fEy[i][j][k] = -fEold;
	}
}

// function to update the Ez electric field component
// umode = 100: Resistor.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEz_type100(int i, int j, int k)
{
	unsigned char uparam;
	double fEold, fpower, fJ;

	uparam = uparamez[i][j][k];
	fEold = fEz[i][j][k];
	fJ = (stCCez[uparam].fC3)*fEold;

	fEz[i][j][k] = (stCCez[uparam].fC0)*fEz[i][j][k] +
				   (stCCez[uparam].fC1)*(fHy[i][j][k] - fHy[i-1][j][k]) -
				   (stCCez[uparam].fC2)*(fHx[i][j][k] - fHx[i][j-1][k]);
	fpower = 0.5*(fEz[i][j][k]+fEold)*fJ;  // compute numerical power
	if (fpower < 0)	// perform compensation if numerical power is -ve
	{
		fEz[i][j][k] = -fEold;
	}
}

// function to update the Ex electric field component
// umode = 101: Capacitor.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEx_type101(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamex[i][j][k];
	fEx[i][j][k] = fEx[i][j][k] +
				   (stCCex[uparam].fC0)*(fHz[i][j][k] - fHz[i][j-1][k]) -
				   (stCCex[uparam].fC1)*(fHy[i][j][k] - fHy[i][j][k-1]);
}

// function to update the Ey electric field component
// umode = 101: Capacitor.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEy_type101(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamey[i][j][k];
	fEy[i][j][k] = fEy[i][j][k] +
				   (stCCey[uparam].fC0)*(fHx[i][j][k] - fHx[i][j][k-1]) -
				   (stCCey[uparam].fC1)*(fHz[i][j][k] - fHz[i-1][j][k]);
}

// function to update the Ez electric field component
// umode = 101: Capacitor.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEz_type101(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamez[i][j][k];
	fEz[i][j][k] = fEz[i][j][k] +
				   (stCCez[uparam].fC0)*(fHy[i][j][k] - fHy[i-1][j][k]) -
				   (stCCez[uparam].fC1)*(fHx[i][j][k] - fHx[i][j-1][k]);
}

// function to update the Ex electric field component
// umode = 150: Diode (2nd order).
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEx_type150(int i, int j, int k)
{
	unsigned char uparam;
	tEDiode *pDiode;

	uparam = uparamex[i][j][k];
	pDiode = (tEDiode *) stCCex[uparam].plist; // cast the void pointer into a 
											   // pointer to tEDiode class
	fEx[i][j][k] = pDiode->updateEfield(fEx[i][j][k],fHy[i][j][k],fHy[i][j][k-1],
										 fHz[i][j][k],fHz[i][j-1][k-1]);
}

// function to update the Ey electric field component
// umode = 150: Diode (2nd order).
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEy_type150(int i, int j, int k)
{
	unsigned char uparam;
	tEDiode *pDiode;

	uparam = uparamey[i][j][k];
	pDiode = (tEDiode *) stCCey[uparam].plist; // cast the void pointer into a 
											   // pointer to tEDiode class
	fEy[i][j][k] = pDiode->updateEfield(fEy[i][j][k],fHz[i][j][k],fHz[i-1][j][k],
										 fHx[i][j][k],fHx[i][j][k-1]);
}

// function to update the Ez electric field component
// umode = 150: Diode (2nd order).
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEz_type150(int i, int j, int k)
{
	unsigned char uparam;
	tEDiode *pDiode;

	uparam = uparamez[i][j][k];
	pDiode = (tEDiode *) stCCez[uparam].plist; // cast the void pointer into a 
											   // pointer to tEDiode class
	fEz[i][j][k] = pDiode->updateEfield(fEz[i][j][k],fHx[i][j][k],fHx[i][j-1][k],
										 fHy[i][j][k],fHy[i-1][j][k]);
}

// function to update the Ez electric field component
// umode = 151: Diode (Ciampolini & Piket-May approach).
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEz_type151(int i, int j, int k)
{
	unsigned char uparam;
	tEDiodec *pDiodec;

	uparam = uparamez[i][j][k];
	pDiodec = (tEDiodec *) stCCez[uparam].plist; // cast the void pointer into a 
											   // pointer to tEDiodec class
	fEz[i][j][k] = pDiodec->updateEfield(fEz[i][j][k],fHx[i][j][k],fHx[i][j-1][k],
										 fHy[i][j][k],fHy[i-1][j][k]);
}

// function to update the Ez electric field component
// umode = 160: PWL two terminal device.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEz_type160(int i, int j, int k)
{
	unsigned char uparam;
	int r=-1; // r to keep track of range
	double I=0, I0=0, I2=0, I_1=0, V=0, V0=0, V2=0, V_1=0;
	double ftemp1, M1, M, Ic;
	lpDOUBLE plist;

	uparam = uparamez[i][j][k]; // assign coefficient table index
	plist = (lpDOUBLE) stCCez[uparam].plist;

	// check the range r of the 2-port voltage
	// The parameter list is in the format:
	// I0 V0 I1 V1 I2 V2 ... IN VN
	while (plist != NULL)
	{
		I = plist->fP; // get In
		plist = plist->next; // point to next item
		if (plist != NULL)
		{
			V = plist->fP; // get Vn
			if ((-fEz[i][j][k]*fdel_z) < V) // Ez falls within the range
			{	
				plist = plist->next; // point to next item
				break;
			}
			else
			{
				r++; // increase range by 1
				plist = plist->next; // point to next item
				V_1 = V0;
				I_1 = I0;
				V0 = V;
				I0 = I;
			}
		}
	}
	
	// compute C and M here
	if (r < 0)
	{
		M = 0;
		Ic = I;
	}
	else 
	{
		if (V != V0)
		{
			M = (I-I0)/(V-V0);
			Ic = I0+M*((-fEz[i][j][k]*fdel_z)-V0);
		}
		else // to prevent divide by zero
		{
			M = 0;
			Ic = I;
		}
	}

	// Update Ez here
	stCCez[uparam].fC4 = fEz[i][j][k]; // backup Ezn
	ftemp1 = (stCCez[uparam].fC1)*(fHy[i][j][k] - fHy[i-1][j][k]) -
    		 (stCCez[uparam].fC2)*(fHx[i][j][k] - fHx[i][j-1][k]);
	ftemp1 = ftemp1/(1+(stCCez[uparam].fC0)*M);
	fEz[i][j][k] = fEz[i][j][k] + ftemp1 +
				   ((stCCez[uparam].fC3)*Ic)/(1+(stCCez[uparam].fC0)*M);

	// To smooth out abrupt changes at junction when V>=V1
	if ((-fEz[i][j][k]*fdel_z)>= V)
	{
		if (plist != NULL)
		{
			I2 = plist->fP; // get In
			plist = plist->next; // point to next item
			if (plist != NULL)
			{
				V2 = plist->fP; // get Vn
				if (V2 != V) // to prevent divide by 0
				{
					M1 = (I2-I)/(V2-V);
					M = (M+M1)/2; // taking the average
				}
				ftemp1 = (stCCez[uparam].fC1)*(fHy[i][j][k] - fHy[i-1][j][k]) -
    					(stCCez[uparam].fC2)*(fHx[i][j][k] - fHx[i][j-1][k]);
				ftemp1 = ftemp1/(1+(stCCez[uparam].fC0)*M);
				fEz[i][j][k] = stCCez[uparam].fC4 + ftemp1 +
							  ((stCCez[uparam].fC3)*Ic)/(1+(stCCez[uparam].fC0)*M);
			}
		}
	}

	// To smooth out abrupt changes at junction when V<V0
	if ((-fEz[i][j][k]*fdel_z) < V0)
	{
		if (V_1 != V0) // to prevent divide by 0
		{
			M1 = (I0-I_1)/(V0-V_1);
			M = (M+M1)/2; // taking the average
		}
		ftemp1 = (stCCez[uparam].fC1)*(fHy[i][j][k] - fHy[i-1][j][k]) -
    			(stCCez[uparam].fC2)*(fHx[i][j][k] - fHx[i][j-1][k]);
		ftemp1 = ftemp1/(1+(stCCez[uparam].fC0)*M);
		fEz[i][j][k] = stCCez[uparam].fC4 + ftemp1 +
					  ((stCCez[uparam].fC3)*Ic)/(1+(stCCez[uparam].fC0)*M);
	}
}

// function to update the Ex electric field component
// umode = 180: Bjt BE junction.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEx_type180(int i, int j, int k)
{
	unsigned char uparam, uor;
	tEBjt *pBjt;

	uparam = uparamex[i][j][k];
	pBjt = (tEBjt *) stCCex[uparam].plist; // cast the void pointer into a pointer to
										   // tEBjt class
	uor = (unsigned char) stCCex[uparam].fC0; // get the orientation of the component.
	pBjt->Update_Efield(this,i,j,k,uor,fHy[i][j][k],fHy[i][j][k-1],fHz[i][j][k],fHz[i][j-1][k]);
}

// function to update the Ey electric field component
// umode = 180: Bjt.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEy_type180(int i, int j, int k)
{
	unsigned char uparam, uor;
	tEBjt *pBjt;

	uparam = uparamey[i][j][k];
	pBjt = (tEBjt *) stCCey[uparam].plist; // cast the void pointer into a pointer to
										   // tEBjt class
	uor = (unsigned char) stCCey[uparam].fC0; // get the orientation of the device.
	pBjt->Update_Efield(this,i,j,k,uor,fHz[i][j][k],fHz[i-1][j][k],fHx[i][j][k],fHx[i][j][k-1]);
}

// function to update the Ez electric field component
// umode = 180: Bjt.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEz_type180(int i, int j, int k)
{
	unsigned char uparam, uor;
	tEBjt *pBjt;

	uparam = uparamez[i][j][k];
	pBjt = (tEBjt *) stCCez[uparam].plist; // cast the void pointer into a pointer to
										   // tEBjt class
	uor = (unsigned char) stCCez[uparam].fC0; // get the orientation of the device.
	pBjt->Update_Efield(this, i,j,k,uor,fHx[i][j][k],fHx[i][j-1][k],fHy[i][j][k],fHy[i-1][j][k]);
}

// function to update the Ex electric field component
// umode = 200: Inductor.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEx_type200(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamex[i][j][k];
	stCCex[uparam].fC3 = stCCex[uparam].fC3 + stCCex[uparam].fC2*fEx[i][j][k];
	fEx[i][j][k] = fEx[i][j][k] +
				   (stCCex[uparam].fC0)*(fHz[i][j][k] - fHz[i][j-1][k]) -
				   (stCCex[uparam].fC1)*(fHy[i][j][k] - fHy[i][j][k-1]) -
				   stCCex[uparam].fC3;
			
}

// function to update the Ey electric field component
// umode = 200: Inductor.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEy_type200(int i, int j, int k)
{
	unsigned char uparam;
	
	uparam = uparamey[i][j][k];
	stCCey[uparam].fC3 = stCCey[uparam].fC3 + stCCey[uparam].fC2*fEy[i][j][k];
	fEy[i][j][k] = fEy[i][j][k] +
				   (stCCey[uparam].fC0)*(fHx[i][j][k] - fHx[i][j][k-1]) -
				   (stCCey[uparam].fC1)*(fHz[i][j][k] - fHz[i-1][j][k]) -
				   stCCey[uparam].fC3;
			

}

// function to update the Ez electric field component
// umode = 200: Inductor.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateEz_type200(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamez[i][j][k];
	stCCez[uparam].fC3 = stCCez[uparam].fC3 + stCCez[uparam].fC2*fEz[i][j][k];
	fEz[i][j][k] = fEz[i][j][k] +
				   (stCCez[uparam].fC0)*(fHy[i][j][k] - fHy[i-1][j][k]) -
				   (stCCez[uparam].fC1)*(fHx[i][j][k] - fHx[i][j-1][k]) -
				   stCCez[uparam].fC3;
			
}

// function to update the Hx magnetic field component
// umode = 0: default update mode.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHx_type0(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamhx[i][j][k];
	fHx[i][j][k] = fHx[i][j][k] -
				   (fChx[uparam][0])*(fEz[i][j+1][k] - fEz[i][j][k]) +
				   (fChx[uparam][1])*(fEy[i][j][k+1] - fEy[i][j][k]);
}

// function to update the Hy magnetic field component
// umode = 0: default update mode.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHy_type0(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamhy[i][j][k];
	fHy[i][j][k] = fHy[i][j][k] -
				   (fChy[uparam][0])*(fEx[i][j][k+1] - fEx[i][j][k]) +
				   (fChy[uparam][1])*(fEz[i+1][j][k] - fEz[i][j][k]);
}

// function to update the Hz magnetic field component
// umode = 0: default update mode.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHz_type0(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamhz[i][j][k];
	fHz[i][j][k] = fHz[i][j][k] -
				   (fChz[uparam][0])*(fEy[i+1][j][k] - fEy[i][j][k]) +
				   (fChz[uparam][1])*(fEx[i][j+1][k] - fEx[i][j][k]);
}

// function to update the Hx magnetic field component
// umode = 1: field normal to PEC.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHx_type1(int i, int j, int k)
{
	fHx[i][j][k] = 0;
}

// function to update the Hy magnetic field component
// umode = 1: field normal PEC.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHy_type1(int i, int j, int k)
{
	fHy[i][j][k] = 0;
}

// function to update the Hz magnetic field component
// umode = 1: field normal to PEC.
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHz_type1(int i, int j, int k)
{
	fHz[i][j][k] = 0;
}

// function to update the Hx magnetic field component
// umode = 2: partial field normal to PEC, upper triangular PEC plane
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHx_type2(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamhx[i][j][k];
	fHx[i][j][k] = fHx[i][j][k] - (fChx[uparam][0])*fEy[i][j][k] -
				   (fChx[uparam][1])*fEz[i][j+1][k];
}

// function to update the Hy magnetic field component
// umode = 2: partial field normal PEC, upper triangular PEC plane
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHy_type2(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamhy[i][j][k];
	fHy[i][j][k] = fHy[i][j][k] + (fChy[uparam][0])*fEx[i][j][k] -
				   (fChy[uparam][1])*fEz[i][j][k];
}

// function to update the Hz magnetic field component
// umode = 2: partial field normal to PEC, upper triangular PEC plane
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHz_type2(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamhz[i][j][k];
	fHz[i][j][k] = fHz[i][j][k] - (fChz[uparam][0])*fEy[i+1][j][k] +
				   (fChz[uparam][1])*fEx[i][j+1][k];
}

// function to update the Hx magnetic field component
// umode = 3: partial field normal to PEC, lower triangular PEC plane
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHx_type3(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamhx[i][j][k];
	fHx[i][j][k] = fHx[i][j][k] + (fChx[uparam][0])*fEy[i][j][k] +
				   (fChx[uparam][1])*fEz[i][j+1][k];
}

// function to update the Hy magnetic field component
// umode = 3: partial field normal PEC, lower triangular PEC plane
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHy_type3(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamhy[i][j][k];
	fHy[i][j][k] = fHy[i][j][k] - (fChy[uparam][0])*fEx[i][j][k] +
				   (fChy[uparam][1])*fEz[i][j][k];
}

// function to update the Hz magnetic field component
// umode = 3: partial field normal to PEC, lower triangular PEC plane
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHz_type3(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamhz[i][j][k];
	fHz[i][j][k] = fHz[i][j][k] + (fChz[uparam][0])*fEy[i][j][k] -
				   (fChz[uparam][1])*fEx[i][j][k];
}

// function to update the Hz magnetic field component
// umode = 4: partial field normal to PEC, cross lower triangular PEC plane
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHz_type4(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamhz[i][j][k];
	fHz[i][j][k] = fHz[i][j][k] + (fChz[uparam][0])*fEy[i][j][k] +
				   (fChz[uparam][1])*fEx[i][j+1][k];
}

// function to update the Hz magnetic field component
// umode = 5: partial field normal to PEC, cross upper triangular PEC plane
// Arguments:
// i = x index of x field component
// j = y index of y field component
// k = z index of z field component
// return: none
void tEHfdtd::UpdateHz_type5(int i, int j, int k)
{
	unsigned char uparam;

	uparam = uparamhz[i][j][k];
	fHz[i][j][k] = fHz[i][j][k] - (fChz[uparam][0])*fEy[i+1][j][k] -
				   (fChz[uparam][1])*fEx[i][j][k];
}