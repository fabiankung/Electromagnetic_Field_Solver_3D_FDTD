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


// Filename		 : fdtd_class_abc.h
// Remark		 : Definition for absorbing boundary condition class
// Progamming Language	 : C++
// Programmer	 : Fabian Kung Wai Lee

#include "fdtd_class.h"

#ifndef _FDTD_ABC_CLASS_
#define _FDTD_ABC_CLASS_

///////////////////////////////////////////////////////////
// --- Class object for Absorbing Boundary Condition --- //
///////////////////////////////////////////////////////////
class tABC // This class is tailored for using Mur's ABC
{
protected:
// data prototype
double fEvalue;

// For storing previous field components at x=0 and x=Max_x boundary
double (*fEy_x011)[MAXZ+1];
double (*fEy_xm11)[MAXZ+1];
double (*fEy_x012)[MAXZ+1];
double (*fEy_xm12)[MAXZ+1];
double (*fEy_x001)[MAXZ+1];
double (*fEy_xm01)[MAXZ+1];
double (*fEy_x002)[MAXZ+1];
double (*fEy_xm02)[MAXZ+1];

double (*fEz_x011)[MAXZ];
double (*fEz_xm11)[MAXZ];
double (*fEz_x012)[MAXZ];
double (*fEz_xm12)[MAXZ];
double (*fEz_x001)[MAXZ];
double (*fEz_xm01)[MAXZ];
double (*fEz_x002)[MAXZ];
double (*fEz_xm02)[MAXZ];

// For storing previous field components at y=0 and y=Max_y boundary
double (*fEx_y011)[MAXZ+1];
double (*fEx_ym11)[MAXZ+1];
double (*fEx_y012)[MAXZ+1];
double (*fEx_ym12)[MAXZ+1];
double (*fEx_y001)[MAXZ+1];
double (*fEx_ym01)[MAXZ+1];
double (*fEx_y002)[MAXZ+1];
double (*fEx_ym02)[MAXZ+1];

double (*fEz_y011)[MAXZ];
double (*fEz_ym11)[MAXZ];
double (*fEz_y012)[MAXZ];
double (*fEz_ym12)[MAXZ];
double (*fEz_y001)[MAXZ];
double (*fEz_ym01)[MAXZ];
double (*fEz_y002)[MAXZ];
double (*fEz_ym02)[MAXZ];

// For storing previous field components at z=0 and z=Max_z boundary
double (*fEx_z011)[MAXY+1];
double (*fEx_zm11)[MAXY+1];
double (*fEx_z012)[MAXY+1];
double (*fEx_zm12)[MAXY+1];
double (*fEx_z001)[MAXY+1];
double (*fEx_zm01)[MAXY+1];
double (*fEx_z002)[MAXY+1];
double (*fEx_zm02)[MAXY+1];

double (*fEy_z011)[MAXY];
double (*fEy_zm11)[MAXY];
double (*fEy_z012)[MAXY];
double (*fEy_zm12)[MAXY];
double (*fEy_z001)[MAXY];
double (*fEy_zm01)[MAXY];
double (*fEy_z002)[MAXY];
double (*fEy_zm02)[MAXY];

// function prototype

public:
// data prototype
int nStatus;
int nNx;
int nNy;
int nNz;

// function prototype
tABC(int, int, int); // constructor
~tABC(); // destructor
void InitABC(int, int, int);
void ShiftEfield(tEHfdtd &);
void EnforceBoundary(tEHfdtd &);
BOOL CheckABCStability(tEHfdtd &);
};

#endif
