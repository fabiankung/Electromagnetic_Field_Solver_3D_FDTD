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

// Filename		 : fdtd_source5.cpp
// Remark		 : Functions for class tEHfdtd
// Programming Language  : C++
// Author	 : Fabian Kung Wai Lee



#include "fdtd_header.h"
#include "fdtd_class.h"

///////////////////////////////////////////////////////////
// --- Functions of Class tHfield --- /////////////////////
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
// --- Default Constructor --- ////////////////////////////
///////////////////////////////////////////////////////////
// The constructor performs 1 basic function:
// (1) Allocate dynamic memories.
// (2) Check for memory error.
// (3) Intialize all data in dynamic memories.
tHfield::tHfield()
{
	// allocate dynamic memories
	fHx	= new double[MAXX+1][MAXY][MAXZ];
	fHy	= new double[MAXX][MAXY+1][MAXZ];
	fHz	= new double[MAXX][MAXY][MAXZ+1];
	umodhx = new unsigned char[MAXX+1][MAXY][MAXZ];
	umodhy = new unsigned char[MAXX][MAXY+1][MAXZ];
	umodhz = new unsigned char[MAXX][MAXY][MAXZ+1];
	uparamhx = new unsigned char[MAXX+1][MAXY][MAXZ];
	uparamhy = new unsigned char[MAXX][MAXY+1][MAXZ];
	uparamhz = new unsigned char[MAXX][MAXY][MAXZ+1];
	fChx = new double[MAXECONST][4]; 
	fChy = new double[MAXECONST][4];
	fChz = new double[MAXECONST][4];

	// carry out memory allocation check
	if (fHx == NULL)
	{
		nHStatus = OUT_OF_MEMORY;
		return;
	}
	if (fHy == NULL)
	{
		nHStatus = OUT_OF_MEMORY;
		return;
	}		
	if (fHz == NULL)
	{
		nHStatus = OUT_OF_MEMORY;
		return;
	}	
	if (umodhx == NULL)
	{
		nHStatus = OUT_OF_MEMORY;
		return;
	}
	if (umodhy == NULL)
	{
		nHStatus = OUT_OF_MEMORY;
		return;
	}
	if (umodhz == NULL)
	{
		nHStatus = OUT_OF_MEMORY;
		return;
	}
	if (uparamhx == NULL)
	{
		nHStatus = OUT_OF_MEMORY;
		return;
	}
	if (uparamhy == NULL)
	{
		nHStatus = OUT_OF_MEMORY;
		return;
	}
	if (uparamhz == NULL)
	{
		nHStatus = OUT_OF_MEMORY;
		return;
	}
	if (fChx == NULL)
	{
		nHStatus = OUT_OF_MEMORY;
		return;
	}
	if (fChy == NULL)
	{
		nHStatus = OUT_OF_MEMORY;
		return;
	}
	if (fChz == NULL)
	{
		nHStatus = OUT_OF_MEMORY;
		return;
	}

	nHStatus = NORMAL;
}


///////////////////////////////////////////////////////////
// --- Default Destructor --- /////////////////////////////
///////////////////////////////////////////////////////////
tHfield::~tHfield()
{
	// remove all dynamic allocated memories
	delete[] fHx;
	delete[] fHy;
	delete[] fHz;
	delete[] umodhx;
	delete[] umodhy;
	delete[] umodhz;
	delete[] uparamhx;
	delete[] uparamhy;
	delete[] uparamhz;
	delete[] fChx;
	delete[] fChy;
	delete[] fChz;

}

///////////////////////////////////////////////////////////
// --- Public & Protected Function Definition --- /////////
///////////////////////////////////////////////////////////

// function to initialize H field components, including computation coefficient
// index and mode
// Argument:
// SimMode = 0 for NORMAL mode, otherwise for REUSE mode.
void tHfield::InitHField(int SimMode)
{
	int ni, nj, nk;

	// Initialize all field components and parameters
	// For x component:
	for (ni=0; ni<=MAXX; ni++)
	{
		for (nj=0; nj<=MAXY-1; nj++)
		{	
			for (nk=0; nk<=MAXZ-1; nk++)
			{
				if (SimMode == NORMAL)	// only clear the field for NORMAL mode
					fHx[ni][nj][nk] = 0;		// Initialize field component
	  		    umodhx[ni][nj][nk] = 0; // Initialize update mode
				uparamhx[ni][nj][nk] = 0;	// Initialize index to computation parameters 
			}
		}
	}
	// For y component:
	for (ni=0; ni<=MAXX-1; ni++)
	{
		for (nj=0; nj<=MAXY; nj++)
		{	
			for (nk=0; nk<=MAXZ-1; nk++)
			{
				if (SimMode == NORMAL)	// only clear the field for NORMAL mode				
					fHy[ni][nj][nk] = 0; // Initialize field component
	  		    umodhy[ni][nj][nk] = 0; // Initialize update mode
				uparamhy[ni][nj][nk] = 0;	// Initialize index to computation parameters 
			}
		}
	}
	// For z component:
	for (ni=0; ni<=MAXX-1; ni++)
	{
		for (nj=0; nj<=MAXY-1; nj++)
		{	
			for (nk=0; nk<=MAXZ; nk++)
			{
				if (SimMode == NORMAL)	// only clear the field for NORMAL mode		
					fHz[ni][nj][nk] = 0; // Initialize field component
	  		    umodhz[ni][nj][nk] = 0; // Initialize update mode
				uparamhz[ni][nj][nk] = 0;	// Initialize index to computation parameters 
			}
		}
	}
}

// function to get the computation mode index for 
// magnetic field components in x direction
// Arguments:
// umode = computation mode
// nx = x index of x field component
// ny = y index of x field component
// nz = z index of x field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tHfield::GetHxMode(int nx, int ny, int nz)
{
	return umodhx[nx][ny][nz]; 
}

// function to get the computation mode index for 
// magnetic field components in y direction
// Arguments:
// umode = computation mode
// nx = x index of y field component
// ny = y index of y field component
// nz = z index of y field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tHfield::GetHyMode(int nx, int ny, int nz)
{
	return umodhy[nx][ny][nz]; 
}

// function to get the computation mode index for 
// magnetic field components in z direction
// Arguments:
// umode = computation mode
// nx = x index of y field component
// ny = y index of y field component
// nz = z index of y field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tHfield::GetHzMode(int nx, int ny, int nz)
{
	return umodhz[nx][ny][nz]; 
}

// function to set the computation parameters index for 
// magnetic field components in x direction
// Arguments:
// uindex = computation index
// nx = x index of y field component
// ny = y index of y field component
// nz = z index of y field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tHfield::SetHxParam(unsigned char uindex, int nx, int ny, int nz)
{
	uparamhx[nx][ny][nz] = uindex;
	return uindex;
}

// function to set the computation parameters index for 
// magnetic field components in y direction
// Arguments:
// uindex = computation index
// nx = x index of y field component
// ny = y index of y field component
// nz = z index of y field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tHfield::SetHyParam(unsigned char uindex, int nx, int ny, int nz)
{
	uparamhy[nx][ny][nz] = uindex;
	return uindex;
}

// function to set the computation parameters index for 
// magnetic field components in z direction
// Arguments:
// uindex = computation index
// nx = x index of y field component
// ny = y index of y field component
// nz = z index of y field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tHfield::SetHzParam(unsigned char uindex, int nx, int ny, int nz)
{
	uparamhz[nx][ny][nz] = uindex;
	return uindex;
}

// function to set the computation mode for magnetic field
// components in x direction
// Arguments:
// umode = computation mode
// nx = x index of x field component
// ny = y index of x field component
// nz = z index of x field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tHfield::SetHxMode(unsigned char umode, int nx, int ny, int nz)
{
	umodhx[nx][ny][nz] = umode;
	return umode;
}

// function to set the computation mode for electric field
// components in y direction
// Arguments:
// umode = computation mode
// nx = x index of y field component
// ny = y index of y field component
// nz = z index of y field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tHfield::SetHyMode(unsigned char umode, int nx, int ny, int nz)
{
	umodhy[nx][ny][nz] = umode;
	return umode;
}

// function to set the computation mode for magnetic field
// components in z direction
// Arguments:
// umode = computation mode
// nx = x index of z field component
// ny = y index of z field component
// nz = z index of z field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tHfield::SetHzMode(unsigned char umode, int nx, int ny, int nz)
{
	umodhz[nx][ny][nz] = umode;
	return umode;
}
