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

// Filename		 : fdtd_source4.cpp 
// Remark		 : Functions for class tABC
// Programming Language  : C++
// Author	 : Fabian Kung Wai Lee	


#include <math.h>  // standard ANSI C/C++ math library header
#include "fdtd_header.h"
#include "fdtd_class_abc.h"

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
// --- Functions of Class tABC --- ////////////////////////
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
// --- Default Constructor --- ////////////////////////////
///////////////////////////////////////////////////////////
// The constructor performs 4 basic functions:
// (1) Allocate dynamic memory for boundary field components.
// (2) Check out for memory error.
// (3) Initialize all the contents of the dynamic memory.
// (4) Initialize the limit of the model, nNx, nNy and nNz.
tABC::tABC(int nx = MAXX-1, int ny = MAXY-1, int nz = MAXZ-1)
{
	// Allocate dynamic memories
	// For storing previous field components at x=0 and x=Max_x boundary
	fEy_x011 = new double[MAXY][MAXZ+1];
	fEy_xm11 = new double[MAXY][MAXZ+1];
	fEy_x012 = new double[MAXY][MAXZ+1];
	fEy_xm12 = new double[MAXY][MAXZ+1];
	fEy_x001 = new double[MAXY][MAXZ+1];
	fEy_xm01 = new double[MAXY][MAXZ+1];
	fEy_x002 = new double[MAXY][MAXZ+1];
	fEy_xm02 = new double[MAXY][MAXZ+1];

	fEz_x011 = new double[MAXY+1][MAXZ];
	fEz_xm11 = new double[MAXY+1][MAXZ];
	fEz_x012 = new double[MAXY+1][MAXZ];
	fEz_xm12 = new double[MAXY+1][MAXZ];
	fEz_x001 = new double[MAXY+1][MAXZ];
	fEz_xm01 = new double[MAXY+1][MAXZ];
	fEz_x002 = new double[MAXY+1][MAXZ];
	fEz_xm02 = new double[MAXY+1][MAXZ];

	// For storing previous field components at y=0 and y=Max_y boundary
	fEx_y011 = new double[MAXX][MAXZ+1];
	fEx_ym11 = new double[MAXX][MAXZ+1];
	fEx_y012 = new double[MAXX][MAXZ+1];
	fEx_ym12 = new double[MAXX][MAXZ+1];
	fEx_y001 = new double[MAXX][MAXZ+1];
	fEx_ym01 = new double[MAXX][MAXZ+1];
	fEx_y002 = new double[MAXX][MAXZ+1];
	fEx_ym02 = new double[MAXX][MAXZ+1];

	fEz_y011 = new double[MAXX+1][MAXZ];
	fEz_ym11 = new double[MAXX+1][MAXZ];
	fEz_y012 = new double[MAXX+1][MAXZ];
	fEz_ym12 = new double[MAXX+1][MAXZ];
	fEz_y001 = new double[MAXX+1][MAXZ];
	fEz_ym01 = new double[MAXX+1][MAXZ];
	fEz_y002 = new double[MAXX+1][MAXZ];
	fEz_ym02 = new double[MAXX+1][MAXZ];

	// For storing previous field components at z=0 and z=Max_z boundary
	fEx_z011 = new double[MAXX][MAXY+1];
	fEx_zm11 = new double[MAXX][MAXY+1];
	fEx_z012 = new double[MAXX][MAXY+1];
	fEx_zm12 = new double[MAXX][MAXY+1];
	fEx_z001 = new double[MAXX][MAXY+1];
	fEx_zm01 = new double[MAXX][MAXY+1];
	fEx_z002 = new double[MAXX][MAXY+1];
	fEx_zm02 = new double[MAXX][MAXY+1];

	fEy_z011 = new double[MAXX+1][MAXY];
	fEy_zm11 = new double[MAXX+1][MAXY];
	fEy_z012 = new double[MAXX+1][MAXY];
	fEy_zm12 = new double[MAXX+1][MAXY];
	fEy_z001 = new double[MAXX+1][MAXY];
	fEy_zm01 = new double[MAXX+1][MAXY];
	fEy_z002 = new double[MAXX+1][MAXY];
	fEy_zm02 = new double[MAXX+1][MAXY];

	// carry out memory allocation check
	if ((fEy_x011 == NULL) || (fEy_x012 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEy_x001 == NULL) || (fEy_x002 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEy_xm11 == NULL) || (fEy_xm12 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEy_xm01 == NULL) || (fEy_xm02 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEz_x011 == NULL) || (fEz_x012 == NULL)) 
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEz_x001 == NULL) || (fEz_x002 == NULL)) 
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEz_xm11 == NULL) || (fEz_xm12 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEz_xm01 == NULL) || (fEz_xm02 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEx_y011 == NULL) || (fEx_y012 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEx_y001 == NULL) || (fEx_y002 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEx_ym11 == NULL) || (fEx_ym12 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEx_ym01 == NULL) || (fEx_ym02 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEz_y011 == NULL) || (fEz_y012 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEz_y001 == NULL)  || (fEz_y002 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEz_ym11 == NULL) || (fEz_ym12 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEz_ym01 == NULL) || (fEz_ym02 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEx_z011 == NULL) || (fEx_z012 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEx_z001 == NULL) || (fEx_z002 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEx_zm11 == NULL) || (fEx_zm12 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEx_zm01 == NULL) || (fEx_zm02 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEy_z011 == NULL) || (fEy_z012 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEy_z001 == NULL) || (fEy_z002 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEy_zm11 == NULL) || (fEy_zm12 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if ((fEy_zm01 == NULL)  || (fEy_zm02 == NULL))
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	fEvalue = 0; // reset E field monitor value

	// Indicate normal termination
	nStatus = NORMAL;
}

///////////////////////////////////////////////////////////
// --- Default Destructor --- /////////////////////////////
///////////////////////////////////////////////////////////
tABC::~tABC()
{
	// remove all dynamic allocated memories
	delete[] fEy_x011;
	delete[] fEy_xm11;
	delete[] fEy_x012;
	delete[] fEy_xm12;
	delete[] fEy_x001;
	delete[] fEy_xm01;
	delete[] fEy_x002;
	delete[] fEy_xm02;

	delete[] fEz_x011;
	delete[] fEz_xm11;
	delete[] fEz_x012;
	delete[] fEz_xm12;
	delete[] fEz_x001;
	delete[] fEz_xm01;
	delete[] fEz_x002;
	delete[] fEz_xm02;

	delete[] fEx_y011;
	delete[] fEx_ym11;
	delete[] fEx_y012;
	delete[] fEx_ym12;
	delete[] fEx_y001;
	delete[] fEx_ym01;
	delete[] fEx_y002;
	delete[] fEx_ym02;

	delete[] fEz_y011;
	delete[] fEz_ym11;
	delete[] fEz_y012;
	delete[] fEz_ym12;
	delete[] fEz_y001;
	delete[] fEz_ym01;
	delete[] fEz_y002;
	delete[] fEz_ym02;

	delete[] fEx_z011;
	delete[] fEx_zm11;
	delete[] fEx_z012;
	delete[] fEx_zm12;
	delete[] fEx_z001;
	delete[] fEx_zm01;
	delete[] fEx_z002;
	delete[] fEx_zm02;

	delete[] fEy_z011;
	delete[] fEy_zm11;
	delete[] fEy_z012;
	delete[] fEy_zm12;
	delete[] fEy_z001;
	delete[] fEy_zm01;
	delete[] fEy_z002;
	delete[] fEy_zm02;
}

///////////////////////////////////////////////////////////
// --- Public functions --- ///////////////////////////////
///////////////////////////////////////////////////////////

// function to initialize the absorbing boundary condition variables
// argument:
// nx = number of cells along x axis in model
// ny = number of cells along y axis in model
// nz = number of cells along z axis in model
// return: none
void tABC::InitABC(int nx, int ny, int nz)
{
	int i,j,k;

// --- Step 1 for Mur's ABC ---
// Initialize all field components for ABC to 0.
	// For x=0 and x = max_x
	for (j=0; j<=MAXY-1; j++)
	{
		for (k=0; k<=MAXZ; k++)
		{
			fEy_x011[j][k] = 0;
			fEy_xm11[j][k] = 0;
			fEy_x012[j][k] = 0;
			fEy_xm12[j][k] = 0;
			fEy_x001[j][k] = 0;
			fEy_xm01[j][k] = 0;
			fEy_x002[j][k] = 0;
			fEy_xm02[j][k] = 0;
		}
	}

	for (j=0; j<=MAXY; j++)
	{
		for (k=0; k<=MAXZ-1; k++)
		{
			fEz_x011[j][k] = 0;
			fEz_xm11[j][k] = 0;
			fEz_x012[j][k] = 0;
			fEz_xm12[j][k] = 0;
			fEz_x001[j][k] = 0;
			fEz_xm01[j][k] = 0;
			fEz_x002[j][k] = 0;
			fEz_xm02[j][k] = 0;
		}
	}

	// For y=0 and y = max_y
	for (i=0; i<=MAXX-1; i++)
	{
		for (k=0; k<=MAXZ; k++)
		{
			fEx_y011[i][k] = 0;
			fEx_ym11[i][k] = 0;
			fEx_y012[i][k] = 0;
			fEx_ym12[i][k] = 0;
			fEx_y001[i][k] = 0;
			fEx_ym01[i][k] = 0;
			fEx_y002[i][k] = 0;
			fEx_ym02[i][k] = 0;
		}
	}

	for (i=0; i<=MAXX; i++)
	{
		for (k=0; k<=MAXZ-1; k++)
		{
			fEz_y011[i][k] = 0;
			fEz_ym11[i][k] = 0;
			fEz_y012[i][k] = 0;
			fEz_ym12[i][k] = 0;
			fEz_y001[i][k] = 0;
			fEz_ym01[i][k] = 0;
			fEz_y002[i][k] = 0;
			fEz_ym02[i][k] = 0;
		}
	}

	// For z=0 and z = max_z
	for (i=0; i<=MAXX-1; i++)
	{
		for (j=0; j<=MAXY; j++)
		{
			fEx_z011[i][j] = 0;
			fEx_zm11[i][j] = 0;
			fEx_z012[i][j] = 0;
			fEx_zm12[i][j] = 0;
			fEx_z001[i][j] = 0;
			fEx_zm01[i][j] = 0;
			fEx_z002[i][j] = 0;
			fEx_zm02[i][j] = 0;
		}
	}

	for (i=0; i<=MAXX; i++)
	{
		for (j=0; j<=MAXY-1; j++)
		{
			fEy_z011[i][j] = 0;
			fEy_zm11[i][j] = 0;
			fEy_z012[i][j] = 0;
			fEy_zm12[i][j] = 0;
			fEy_z001[i][j] = 0;
			fEy_zm01[i][j] = 0;
			fEy_z002[i][j] = 0;
			fEy_zm02[i][j] = 0;
		}
	}

	// Initialize model limit
	nNx = nx;
	nNy = ny;
	nNz = nz;
}

// function to shift the boundary electric field to implement 
// Mur's ABC
// Arguments: 
// &tEH - The address to class tEHfdtd.
// return: none
void tABC::ShiftEfield(tEHfdtd &tEH)
{
	int ni, nj, nk;

	// --- Step 2 of Mur's ABC ---
	// for x=0 and x=nNx+1 for Ey
	for (nj=1; nj<=nNy; nj++)
	{
		for (nk=1; nk<=(nNz+1); nk++)
		{
		 // shift Ey(n-1)[2,j,k]	
		 fEy_x011[nj][nk] = tEH.fEy[2][nj][nk];
		 fEy_x012[nj][nk] = fEy_x011[nj][nk];
		 // shift Ey(n-1)[nNx,j,k]	
		 fEy_xm11[nj][nk] = tEH.fEy[nNx][nj][nk];
		 fEy_xm12[nj][nk] = fEy_xm11[nj][nk];
		 // shift Ey(n-1)[1,j,k]
		 fEy_x001[nj][nk] = tEH.fEy[1][nj][nk];
		 fEy_x002[nj][nk] = fEy_x001[nj][nk];
		 // shift Ey(n-1)[nNx+1,j,k]	
		 fEy_xm01[nj][nk] = tEH.fEy[nNx+1][nj][nk];
		 fEy_xm02[nj][nk] = fEy_xm01[nj][nk];
		}
	}

	// for x=0 and x=nNx+1 for Ez
	for (nj=1; nj<=(nNy+1); nj++)
	{
		for (nk=1; nk<=nNz; nk++)
		{
		 // shift Ez(n-1)[2,j,k]	 
		 fEz_x011[nj][nk] = tEH.fEz[2][nj][nk];
		 fEz_x012[nj][nk] = fEz_x011[nj][nk];
		 // shift Ez(n-1)[nNx,j,k]	
		 fEz_xm11[nj][nk] = tEH.fEz[nNx][nj][nk];
		 fEz_xm12[nj][nk] = fEz_xm11[nj][nk];
		 // shift Ez(n-1)[1,j,k]	 
		 fEz_x001[nj][nk] = tEH.fEz[1][nj][nk];
		 fEz_x002[nj][nk] = fEz_x001[nj][nk];
		 // shift Ez(n-1)[nNx+1,j,k]	
		 fEz_xm01[nj][nk] = tEH.fEz[nNx+1][nj][nk];
		 fEz_xm02[nj][nk] = fEz_xm01[nj][nk];
		}
	}

// for y=0 and y=nNy+1 for Ex
	for (ni=1; ni<=nNx; ni++)
	{
		for (nk=1; nk<=(nNz+1); nk++)
		{
		 // shift Ex(n-1)[i,2,k]	
		 fEx_y011[ni][nk] = tEH.fEx[ni][2][nk];
		 fEx_y012[ni][nk] = fEx_y011[ni][nk];
		 // shift Ex(n-1)[i,nNy,k]	
		 fEx_ym11[ni][nk] = tEH.fEx[ni][nNy][nk];
		 fEx_ym12[ni][nk] = fEx_ym11[ni][nk];
		 // shift Ex(n-1)[i,1,k]	
		 fEx_y001[ni][nk] = tEH.fEx[ni][1][nk];
		 fEx_y002[ni][nk] = fEx_y001[ni][nk];
		 // shift Ex(n-1)[i,nNy+1,k]	
		 fEx_ym01[ni][nk] = tEH.fEx[ni][nNy+1][nk];
		 fEx_ym02[ni][nk] = fEx_ym01[ni][nk];
		}
	}

// for y=0 and y=nNy+1 for Ez
	for (ni=1; ni<=(nNx+1); ni++)
	{
		for (nk=1; nk<=nNz; nk++)
		{	
		 // shift Ez(n-1)[i,2,k]
		 fEz_y011[ni][nk] = tEH.fEz[ni][2][nk];
		 fEz_y012[ni][nk] = fEz_y011[ni][nk];
		 // shift Ez(n-1)[i,nNy,k]	
		 fEz_ym11[ni][nk] = tEH.fEz[ni][nNy][nk];
		 fEz_ym12[ni][nk] = fEz_ym11[ni][nk];
		 // shift Ez(n-1)[i,1,k]
		 fEz_y001[ni][nk] = tEH.fEz[ni][1][nk];
		 fEz_y002[ni][nk] = fEz_y001[ni][nk];
		 // shift Ez(n-1)[i,nNy+1,k]	
		 fEz_ym01[ni][nk] = tEH.fEz[ni][nNy+1][nk];
		 fEz_ym02[ni][nk] = fEz_ym01[ni][nk];
		}
	}

// for z=0 and z=nNz+1 for Ex
	for (ni=1; ni<=nNx; ni++)
	{
		for (nj=1; nj<=(nNy+1); nj++)
		{
		 // shift Ex(n-1)[i,j,2]	
		 fEx_z011[ni][nj] = tEH.fEx[ni][nj][2];
		 fEx_z012[ni][nj] = fEx_z011[ni][nj];
		 // shift Ex(n-1)[i,j,nNz]	
		 fEx_zm11[ni][nj] = tEH.fEx[ni][nj][nNz];
		 fEx_zm12[ni][nj] = fEx_zm11[ni][nj];
		 // shift Ex(n-1)[i,j,1]	
		 fEx_z001[ni][nj] = tEH.fEx[ni][nj][1];
		 fEx_z002[ni][nj] = fEx_z001[ni][nj];
		 // shift Ex(n-1)[i,j,nNz+1]	
		 fEx_zm01[ni][nj] = tEH.fEx[ni][nj][nNz+1];
		 fEx_zm02[ni][nj] = fEx_zm01[ni][nj];
		}
	}

// for z=0 and z=nNz+1 for Ey
	for (ni=1; ni<=(nNx+1); ni++)
	{
		for (nj=1; nj<=nNy; nj++)
		{
		 // shift Ey(n-1)[i,j,2]	
		 fEy_z011[ni][nj] = tEH.fEy[ni][nj][2];
		 fEy_z012[ni][nj] = fEy_z011[ni][nj];
		 // shift Ey(n-1)[i,j,nNz]	
		 fEy_zm11[ni][nj] = tEH.fEy[ni][nj][nNz];
		 fEy_zm12[ni][nj] = fEy_zm11[ni][nj];
		 // shift Ey(n-1)[i,j,1]	
		 fEy_z001[ni][nj] = tEH.fEy[ni][nj][1];
		 fEy_z002[ni][nj] = fEy_z001[ni][nj];
		 // shift Ey(n-1)[i,j,nNz+1]	
		 fEy_zm01[ni][nj] = tEH.fEy[ni][nj][nNz+1];
		 fEy_zm02[ni][nj] = fEy_zm01[ni][nj];
		}
	}
}

// function to enforce the boundary condition on electric field 
// to implement Mur's ABC
// Arguments: 
// &tEH - The address to class tEHfdtd.
// return: none
void tABC::EnforceBoundary(tEHfdtd &tEH)
{
	int ni,nj,nk;

	double Cx;
	double Cy;
	double Cz;

	Cx = (c*fdel_t - fdel_x) / (c*fdel_t + fdel_x);
	Cy = (c*fdel_t - fdel_y) / (c*fdel_t + fdel_y);
	Cz = (c*fdel_t - fdel_z) / (c*fdel_t + fdel_z);
// Enforce ABC for x=1 and x=nNx+1
	// For Ey
	for (nj=1; nj<=nNy; nj++)
	{
		// E(n+1)[0] = E(n)[1] + C[E(n+1)[1] - E(n)[0])
		for (nk=1; nk<=(nNz+1); nk++)
		{
			tEH.fEy[1][nj][nk] = fEy_x011[nj][nk] + Cx*(tEH.fEy[2][nj][nk] - tEH.fEy[1][nj][nk]);
			tEH.fEy[nNx+1][nj][nk] = fEy_xm11[nj][nk] + Cx*(tEH.fEy[nNx][nj][nk] - tEH.fEy[nNx+1][nj][nk]); 
		}
	}
	// For Ez
	for (nj=1; nj<=(nNy+1); nj++)
	{
		for (nk=1; nk<=nNz; nk++)
		{
			tEH.fEz[1][nj][nk] = fEz_x011[nj][nk] + Cx*(tEH.fEz[2][nj][nk] - tEH.fEz[1][nj][nk]);
			tEH.fEz[nNx+1][nj][nk] = fEz_xm11[nj][nk] + Cx*(tEH.fEz[nNx][nj][nk] - tEH.fEz[nNx+1][nj][nk]); 
		}
	}

// Enforce ABC for y=1 and y=nNy+1
	// For Ex
	for (ni=1; ni<=nNx; ni++)
	{
		for (nk=1; nk<=(nNz+1); nk++)
		{
			tEH.fEx[ni][1][nk] = fEx_y011[ni][nk] + Cy*(tEH.fEx[ni][2][nk] - tEH.fEx[ni][1][nk]);
			tEH.fEx[ni][nNy+1][nk] = fEx_ym11[ni][nk] + Cy*(tEH.fEx[ni][nNy][nk] - tEH.fEx[ni][nNy+1][nk]);
		}
	}
	// For Ez
	for (ni=1; ni<=(nNx+1); ni++)
	{
		for (nk=1; nk<=nNz; nk++)
		{
			tEH.fEz[ni][1][nk] = fEz_y011[ni][nk] + Cy*(tEH.fEz[ni][2][nk] - tEH.fEz[ni][1][nk]);
			tEH.fEz[ni][nNy+1][nk] = fEz_ym11[ni][nk] + Cy*(tEH.fEz[ni][nNy][nk] - tEH.fEz[ni][nNy+1][nk]);
		}
	}

// Enforce ABC for z=1 and z=nNz+1
	// For Ex
	for (ni=1; ni<=nNx; ni++)
	{
		for (nj=1; nj<=(nNy+1); nj++)
		{
			tEH.fEx[ni][nj][1] = fEx_z011[ni][nj] + Cz*(tEH.fEx[ni][nj][2] - tEH.fEx[ni][nj][1]);
			tEH.fEx[ni][nj][nNz+1] = fEx_zm11[ni][nj] + Cz*(tEH.fEx[ni][nj][nNz] - tEH.fEx[ni][nj][nNz+1]);
		}
	}
	// For Ey
	for (ni=1; ni<=(nNx+1); ni++)
	{
		for (nj=1; nj<=nNy; nj++)
		{
			tEH.fEy[ni][nj][1] = fEy_z011[ni][nj] + Cz*(tEH.fEy[ni][nj][2] - tEH.fEy[ni][nj][1]);
			tEH.fEy[ni][nj][nNz+1] = fEy_zm11[ni][nj] + Cz*(tEH.fEy[ni][nj][nNz] - tEH.fEy[ni][nj][nNz+1]);
		}
	}
}

// function to check the stability of absorbing boundary condition 
// The magnitude of Ex field component on the center of top, bottom, right and 
// left of the model space will be summed up.  Similarly the magnitude of Ey on
// the center of front and back of the model space will be summed up.  The total
// of the two values will be constantly checked against the INSTABILITY_LIMIT 
// constant when this function is called.
// argument: tEH = pointer to object tEHfdtd
// return: TRUE if check value fEvalue exceed the limit INSTABILITY_LIMIT
// as defined in the fdtd_header.h header file.
BOOL tABC::CheckABCStability(tEHfdtd &tEH)
{
	fEvalue = tEH.fEx[(int)ceil(nNx/2)][1][(int)ceil(nNz/2)]*
				tEH.fEx[(int)ceil(nNx/2)][1][(int)ceil(nNz/2)]; 
	fEvalue += tEH.fEx[(int)ceil(nNx/2)][nNy+1][(int)ceil(nNz/2)]*
				tEH.fEx[(int)ceil(nNx/2)][nNy+1][(int)ceil(nNz/2)];
	fEvalue += tEH.fEx[(int)ceil(nNx/2)][(int)ceil(nNy/2)][1]*
				tEH.fEx[(int)ceil(nNx/2)][(int)ceil(nNy/2)][1];
	fEvalue += tEH.fEx[(int)ceil(nNx/2)][(int)ceil(nNy/2)][nNz+1]*
				tEH.fEx[(int)ceil(nNx/2)][(int)ceil(nNy/2)][nNz+1];
	fEvalue += tEH.fEy[1][(int)ceil(nNy/2)][(int)ceil(nNz/2)]*
				tEH.fEy[nNx+1][(int)ceil(nNy/2)][(int)ceil(nNz/2)];
	if (sqrt(fEvalue) > INSTABILITY_LIMIT)
	{	
		fEvalue = 0; // reset monitor value
		return TRUE;
	}
	else
	{
		fEvalue = 0; // reset monitor value
		return FALSE;
	}
}
