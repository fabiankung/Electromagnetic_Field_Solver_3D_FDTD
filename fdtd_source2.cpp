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

// Filename		 : fdtd_source2.cpp
// Remark		 : Functions for class tEfield
// Programming Language  : C++
// Author		 : Fabian Kung Wai Lee

#include "fdtd_header.h"
#include "fdtd_class.h"

///////////////////////////////////////////////////////////
// --- Default Constructor --- ////////////////////////////
///////////////////////////////////////////////////////////
// The constructor performs 5 basic functions:
// (1) Allocate dynamic memories.
// (2) Check for out of memory error.
// (3) set all electric field components to 0.
// (4) set all parameters to 0.
// (5) set the field update mode to 0 - normal update.
tEfield::tEfield()
{
	// Allocate dynamic memories
	umodex = new unsigned char[MAXX][MAXY+1][MAXZ+1];
	umodey = new unsigned char[MAXX+1][MAXY][MAXZ+1];
	umodez = new unsigned char[MAXX+1][MAXY+1][MAXZ];
	uparamex = new unsigned char[MAXX][MAXY+1][MAXZ+1];
	uparamey = new unsigned char[MAXX+1][MAXY][MAXZ+1];
	uparamez = new unsigned char[MAXX+1][MAXY+1][MAXZ];
	fEx	= new double[MAXX][MAXY+1][MAXZ+1];
	fEy	= new double[MAXX+1][MAXY][MAXZ+1];
	fEz	= new double[MAXX+1][MAXY+1][MAXZ];

	// carry out memory allocation check
	if (umodex == NULL)
	{
		nEStatus = OUT_OF_MEMORY;
		return;
	}
	if (umodey == NULL)
	{
		nEStatus = OUT_OF_MEMORY;
		return;
	}		
	if (umodez == NULL)
	{
		nEStatus = OUT_OF_MEMORY;
		return;
	}		
	if (uparamex == NULL)
	{
		nEStatus = OUT_OF_MEMORY;
		return;
	}		
	if (uparamey == NULL)
	{
		nEStatus = OUT_OF_MEMORY;
		return;
	}		
	if (uparamez == NULL)
	{
		nEStatus = OUT_OF_MEMORY;
		return;
	}		
	if (fEx == NULL)
	{
		nEStatus = OUT_OF_MEMORY;
		return;
	}		
	if (fEy == NULL)
	{
		nEStatus = OUT_OF_MEMORY;
		return;
	}	
	if (fEz == NULL)
	{
		nEStatus = OUT_OF_MEMORY;
		return;
	}		
	if (stCCex == NULL)
	{
		nEStatus = OUT_OF_MEMORY;
		return;
	}
	if (stCCey == NULL)
	{
		nEStatus = OUT_OF_MEMORY;
		return;
	}
	if (stCCez == NULL)
	{
		nEStatus = OUT_OF_MEMORY;
		return;
	}

	nEStatus = NORMAL;
}

///////////////////////////////////////////////////////////
// --- Default Destructor --- /////////////////////////////
///////////////////////////////////////////////////////////
tEfield::~tEfield()
{
	// remove all dynamic allocated memories
	delete[] umodex;
	delete[] umodey;
	delete[] umodez;
	delete[] uparamex;
	delete[] uparamey;
	delete[] uparamez;
	delete[] fEx;
	delete[] fEy;
	delete[] fEz;
}


///////////////////////////////////////////////////////////
// --- Public functions --- ///////////////////////////////
///////////////////////////////////////////////////////////

// function to initialize all E field components, computation coefficient
// index and mode.
// Argument:
// SimMode = 0 for NORMAL, otherwise for REUSE
void tEfield::InitEField(int SimMode)
{
	int ni,nj,nk;

	// Initialize all field components and parameters
	// For x component:
	for (ni=0; ni<=MAXX-1; ni++)
	{
		for (nj=0; nj<=MAXY; nj++)
		{
			for (nk=0; nk<=MAXZ; nk++)
			{
			  if (SimMode == NORMAL) // Only clear the field on NORMAL mode
				fEx[ni][nj][nk] = 0; // Initialize field component
			  umodex[ni][nj][nk] = 0; // Initialize update mode
			  uparamex[ni][nj][nk] = 0; // Initialize parameter
			}
		}
	}

	// For y component:
	for (ni=0; ni<=MAXX; ni++)
	{
		for (nj=0; nj<=MAXY-1; nj++)
		{
			for (nk=0; nk<=MAXZ; nk++)
			{
			  if (SimMode == NORMAL) // Only clear the field on NORMAL mode
				fEy[ni][nj][nk] = 0; // Initialize field component
			  umodey[ni][nj][nk] = 0; // Initialize update mode
			  uparamey[ni][nj][nk] = 0; // Initialize parameter
			}
		}
	}

	// For z component:
	for (ni=0; ni<=MAXX; ni++)
	{
		for (nj=0; nj<=MAXY; nj++)
		{
			for (nk=0; nk<=MAXZ-1; nk++)
			{
			 if (SimMode == NORMAL) // Only clear the field on NORMAL mode
				fEz[ni][nj][nk] = 0; // Initialize field component
			 umodez[ni][nj][nk] = 0; // Initialize update mode
			 uparamez[ni][nj][nk] = 0; // Initialize parameter
			}
		}
	}
}

// function to set the computation mode for electric field
// components in x direction
// Arguments:
// umode = computation mode
// nx = x index of x field component
// ny = y index of x field component
// nz = z index of x field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tEfield::SetExMode(unsigned char umode, int nx, int ny, int nz)
{
	umodex[nx][ny][nz] = umode;
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
unsigned char tEfield::SetEyMode(unsigned char umode, int nx, int ny, int nz)
{
	umodey[nx][ny][nz] = umode;
	return umode;
}

// function to set the computation mode for electric field
// components in z direction
// Arguments:
// umode = computation mode
// nx = x index of z field component
// ny = y index of z field component
// nz = z index of z field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tEfield::SetEzMode(unsigned char umode, int nx, int ny, int nz)
{
	umodez[nx][ny][nz] = umode;
	return umode;
}

// function to set the computation parameters index for 
// electric field components in x direction
// Arguments:
// uindex = computation index
// nx = x index of x field component
// ny = y index of x field component
// nz = z index of x field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tEfield::SetExParam(unsigned char uindex, int nx, int ny, int nz)
{
	uparamex[nx][ny][nz] = uindex;
	return uindex;
}

// function to set the computation parameters index for 
// electric field components in y direction
// Arguments:
// uindex = computation index
// nx = x index of y field component
// ny = y index of y field component
// nz = z index of y field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tEfield::SetEyParam(unsigned char uindex, int nx, int ny, int nz)
{
	uparamey[nx][ny][nz] = uindex;
	return uindex;
}

// function to set the computation parameters index for 
// electric field components in z direction
// Arguments:
// uindex = computation index
// nx = x index of z field component
// ny = y index of z field component
// nz = z index of z field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tEfield::SetEzParam(unsigned char uindex, int nx, int ny, int nz)
{
	uparamez[nx][ny][nz] = uindex;
	return uindex;
}

// function to get the computation mode index for 
// electric field components in x direction
// Arguments:
// umode = computation mode
// nx = x index of x field component
// ny = y index of x field component
// nz = z index of x field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tEfield::GetExMode(int nx, int ny, int nz)
{
	return umodex[nx][ny][nz]; 
}

// function to get the computation mode index for 
// electric field components in y direction
// Arguments:
// umode = computation mode
// nx = x index of y field component
// ny = y index of y field component
// nz = z index of y field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tEfield::GetEyMode(int nx, int ny, int nz)
{
	return umodey[nx][ny][nz]; 
}

// function to get the computation mode index for 
// electric field components in z direction
// Arguments:
// umode = computation mode
// nx = x index of z field component
// ny = y index of z field component
// nz = z index of z field component
// return: 
// If succesfull the computation mode for the field  
unsigned char tEfield::GetEzMode(int nx, int ny, int nz)
{
	return umodez[nx][ny][nz]; 
}