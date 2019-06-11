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

// Filename		 : fdtd_source1.cpp
// Remark		 : Functions for class tCube
// Programming Language	 : C++
// Author	 : Fabian Kung wai Lee

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
extern lpCOMP lpcomp1;
extern int nSimMode;

///////////////////////////////////////////////////////////
// --- Functions of Class tCube --- ///////////////////////
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
// --- Default Constructor --- ////////////////////////////
///////////////////////////////////////////////////////////
// The constructor performs 4 basic functions:
// (1) Allocate dynamic memory.
// (2) Check for out of memory error.
// (3) Initialize the dielectric constant of all E field
//     components to 1.
// (4) Initialize the limit of the model.
tCube::tCube(int nx=10, int ny=10, int nz=10)
{
	int ni, nj, nk;

	// initializing the array pointers
//	Px = new stDOUBLE*[MAXX][MAXY+1][MAXZ+1];
//	Py = new stDOUBLE*[MAXX+1][MAXY][MAXZ+1];
//	Pz = new stDOUBLE*[MAXX+1][MAXY+1][MAXZ];
	stCube = new stCUBE[MAXX][MAXY][MAXZ];

	if (stCube == NULL)
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}

// initialize the dielectric constant, characteristics and
// optional parameter of each cube.
for (ni=0; ni<=MAXX-1; ni++)
{
	for (nj=0; nj<=MAXY-1; nj++)
	{
		for (nk=0; nk<=MAXZ-1; nk++)
		{
		  stCube[ni][nj][nk].e = 1;
		  stCube[ni][nj][nk].sg = 0;
		  stCube[ni][nj][nk].n = 0;
		  stCube[ni][nj][nk].plist = NULL;
		}
	}
}

	// Initialize the limit of the model
	nNx = nx;
	nNy = ny;
	nNz = nz;
	// Initialize computational constant look-up table index
	uEimaxx = 0;
	uEimaxy = 0;
	uEimaxz = 0;
	uHimaxx = 0;
	uHimaxy = 0;
	uHimaxz = 0;
	// Initialize the object tracking pointer and index
	ptObjtrack = NULL;
	uIndtrack = 0;
	// Indicate normal termination
	nStatus = NORMAL;
}

///////////////////////////////////////////////////////////
// --- Default Destructor --- /////////////////////////////
///////////////////////////////////////////////////////////
tCube::~tCube()
{
	int ni,nj,nk,ni2,nj2,nk2;
	lpDOUBLE lpdtemp, lpdtemp1;

	// remove all linked list
for (ni=1; ni<=nNx; ni++)
{
	for (nj=1; nj<=nNy; nj++)
	{
		for (nk=1; nk<=nNz; nk++)
		{
			if (stCube[ni][nj][nk].plist != NULL)
			{	
				lpdtemp1 = stCube[ni][nj][nk].plist; // store the head pointer of the list
				while (stCube[ni][nj][nk].plist->next != NULL)
				{
					lpdtemp = stCube[ni][nj][nk].plist;
					stCube[ni][nj][nk].plist = (stCube[ni][nj][nk].plist)->next;
					delete lpdtemp;
				}
				delete stCube[ni][nj][nk].plist;
				stCube[ni][nj][nk].plist = NULL;
				// recursively compares with all other cubes, if cube has similar pointer
				// values set the pointer to NULL.  This would prevent deleting memory 
				// location which has been disaollocated.  Since current indexes are at ni,
				// nj, nk, the recursive check begins at the current indexes (pointer 
				// plist of previous stCubes is already NULL).
				for (ni2=ni; ni2<=nNx; ni2++)
				{
					for (nj2=nj; nj2<=nNy; nj2++)
					{
						for (nk2=nk; nk2<=nNz; nk2++)
						{
							if (stCube[ni2][nj2][nk2].plist == lpdtemp1)
								stCube[ni2][nj2][nk2].plist = NULL;
						}
					}
				}
			}
		}
	}
}

	// remove all dynamic allocated memories
	delete[] stCube;
}

///////////////////////////////////////////////////////////
// --- Public functions --- ///////////////////////////////
///////////////////////////////////////////////////////////

// function to update the content of stCube structure within the object
// arguments:
// nx = x position of cube
// ny = y position of cube
// nz = z position of cube
// fEpsilon = dielectric constant of cube
// nAttribute = characteristic of cube
// fParam = additional parameter
// return:
// TRUE if success
// FALSE if fail
BOOL tCube::UpdateCube(int nx, int ny, int nz, float fEpsilon, float fSigma, int nDieAttrib,
					   int nSurAttrib, lpDOUBLE lpdouble)
{
	// check for size error
	if (nx > nNx)
	{
		return FALSE;
	}
	if (ny > nNy) 
	{
		return FALSE;
	}
	if (nz > nNz)
	{
		return FALSE;
	}
	stCube[nx][ny][nz].e = fEpsilon;
	stCube[nx][ny][nz].sg = fSigma;
	stCube[nx][ny][nz].d = nDieAttrib;
	stCube[nx][ny][nz].n = nSurAttrib;
	stCube[nx][ny][nz].plist = lpdouble;
	return TRUE;
}

// function to assign dielectric constant to each electric field component
// arguments:
// utype = field type, x, y or z.
// ni = index along x axis
// nj = index along y axis
// nk = index along z axis
// fep = reference to the dielectric constant
// fsigma = reference to the conductivity
// return: none
// If successful: nStatus = NORMAL.
// If not successfull: nStatus = OUT_OF_BOUND.
void tCube::AssignEpsilonSigma(unsigned char utype, int ni, int nj, int nk, float &fe,
								float &fsigma)
{
	float C1e,C2e,C3e,C4e; // temporary variables for dielectric constant
	float C1s,C2s,C3s,C4s; // temporary variables for conductivity constant
// Note: The author is aware that index for array begins at 0 in C++.  However to 
// be consistent with the algorithm in MATLAB, we begin the index at 1.  Thus 
// some memory is wasted.


	// check for size error
	if (nNx > MAXX-1)
	{
		nStatus = OUT_OF_BOUND;
	}
	if (nNy > MAXY-1) 
	{
		nStatus = OUT_OF_BOUND;
	}
	if (nNz > MAXZ-1)
	{
		nStatus = OUT_OF_BOUND;
	}

	switch(utype)
	{
	case 'x':
	// For Ex: ex = (C1e+C2e+C3e+C4e)/4
	//     Ex: sigmax = (C1s+C2s+C3s+C4s)/4
		 // for C1e and C1s:
         if (nk == (nNz+1))
         {
			 // nk = nNz+1 
               if ((nj-1)==nNy)
			   {
				   C1e = stCube[ni][nNy][nNz].e;
				   C1s = stCube[ni][nNy][nNz].sg;
               }
			   else
               {
				   C1e = stCube[ni][nj][nNz].e;
				   C1s = stCube[ni][nj][nNz].sg;
			   }
		 }
		 else
		 {               
			// nk < nNz+1
               if ((nj-1)==nNy)
			   {
               	   C1e = stCube[ni][nNy][nk].e;
	           	   C1s = stCube[ni][nNy][nk].sg;
			   }
               else
			   {
				   C1e = stCube[ni][nj][nk].e;
				   C1s = stCube[ni][nj][nk].sg;
			   }
		 }                 
               
         // for C2e and C2s:
         if (nk==(nNz+1))
		 {
			 // nk = nNz+1 
               if ((nj-1)==0)
			   {
               	C2e = stCube[ni][nj][nNz].e;
              	C2s = stCube[ni][nj][nNz].sg;
			   }
               else
			   {
                C2e = stCube[ni][nj-1][nNz].e;
				C2s = stCube[ni][nj-1][nNz].sg;
			   }
		 } 
         else
		 {
			 // nk < nNz+1
               if ((nj-1)==0)
			   {
               	C2e = stCube[ni][nj][nk].e;
               	C2s = stCube[ni][nj][nk].sg;
			   }
               else
			   {
                C2e = stCube[ni][nj-1][nk].e;
                C2s = stCube[ni][nj-1][nk].sg;
			   }
		 }         
         
         // for C3e and C3s:
         if (nk==1)
		 {
			 // nk = 1 
               if ((nj-1)==0)
			   {
               	C3e = stCube[ni][nj][1].e;
               	C3s = stCube[ni][nj][1].sg;
			   }
               else
			   {
                C3e = stCube[ni][nj-1][1].e;
                C3s = stCube[ni][nj-1][1].sg;
			   }
		 }
		 else
		 {
			 // nk > 1
               if ((nj-1)==0)
			   {
               	C3e = stCube[ni][nj][nk-1].e;
               	C3s = stCube[ni][nj][nk-1].sg;
               }
			   else
			   {
                C3e = stCube[ni][nj-1][nk-1].e;
                C3s = stCube[ni][nj-1][nk-1].sg;
			   }
		 }
         
         // for C4e and C4s:        
         if (nk==1)
		 {
			 // k = 1 
               if ((nj-1)==nNy)
			   {
               	C4e = stCube[ni][nNy][1].e;
              	C4s = stCube[ni][nNy][1].sg;
			   }
               else
			   {
                C4e = stCube[ni][nj][1].e;
                C4s = stCube[ni][nj][1].sg;
			   }
		 }
		 else
		 {
			 // k > 1
               if ((nj-1)==nNy)
			   {
               	C4e = stCube[ni][nNy][nk-1].e;
               	C4s = stCube[ni][nNy][nk-1].sg;
			   }
               else
			   {
                C4e = stCube[ni][nj][nk-1].e;
                C4s = stCube[ni][nj][nk-1].sg;
			   }
		 }
		 // Compute the average conductivity
		 fsigma = (C1s+C2s+C3s+C4s)/4;
         // Compute the average dielectric constant 
		 fe = (C1e+C2e+C3e+C4e)/4;		
		 break;	

	case 'y':
	// For Ey: ey = (C1e+C2e+C3e+C4e)/4
	//     Ey: sigmay = (C1s+C2s+C3s+C4)/4
		// for C1e and C1s:
         if (nk == (nNz+1))
		 {
			 // nk = nNz+1 
               if ((ni-1)==0)
			   {
               	C1e = stCube[1][nj][nNz].e;
              	C1s = stCube[1][nj][nNz].sg;
			   }
               else
			   {
                C1e = stCube[ni-1][nj][nNz].e;
                C1s = stCube[ni-1][nj][nNz].sg;
			   }
		 }
		 else
		 {
			 // nk < Nnz+1
               if ((ni-1)==0)
			   {
               	C1e = stCube[1][nj][nk].e;
               	C1s = stCube[1][nj][nk].sg;
			   }
               else
			   {
                C1e = stCube[ni-1][nj][nk].e;
                C1s = stCube[ni-1][nj][nk].sg;
			   }
		 }                 
            
         
         // for C2e and C2s:
         if (nk == (nNz+1))
		 {
			 // k = nNz+1 
               if ((ni-1)==nNx)
			   {
               	C2e = stCube[nNx][nj][nNz].e;
               	C2s = stCube[nNx][nj][nNz].sg;
			   }
               else
			   {
                C2e = stCube[ni][nj][nNz].e;
                C2s = stCube[ni][nj][nNz].sg;
			   }
		 }
         else
		 {
			// nk < nNz+1
               if ((ni-1)==nNx)
			   {
               	C2e = stCube[nNx][nj][nk].e;
               	C2s = stCube[nNx][nj][nk].sg;
			   }
               else
			   {
                C2e = stCube[ni][nj][nk].e;
                C2s = stCube[ni][nj][nk].sg;
			   }
		 }        
         
         // for C3e and C3s:
         if (nk ==1)
		 {
			 // nk = 1 
               if ((ni-1)==nNx)
			   {
               	C3e = stCube[nNx][nj][1].e;
              	C3s = stCube[nNx][nj][1].sg;
			   }
			   else
			   {
                C3e = stCube[ni][nj][1].e;
                C3s = stCube[ni][nj][1].sg;
			   }     
		 }
		 else
		 {
			 // nk > 1
               if ((ni-1)==nNx)
			   {
               	C3e = stCube[nNx][nj][nk-1].e;
               	C3s = stCube[nNx][nj][nk-1].sg;
			   }
			   else
			   {
                C3e = stCube[ni][nj][nk-1].e;    
                C3s = stCube[ni][nj][nk-1].sg; 
			   }
         }

         // for C4e and C4s:        
         if (nk == 1)
		 {
			 // nk = 1 
               if ((ni-1)==0)
			   {
               	C4e = stCube[1][nj][1].e;
               	C4s = stCube[1][nj][1].sg;
			   }
               else
			   {
                C4e = stCube[ni-1][nj][1].e;
                C4s = stCube[ni-1][nj][1].sg;
			   }
		 }
		 else
		 {
             // nk > 1
               if ((ni-1)==0)
			   {
               	C4e = stCube[1][nj][nk-1].e;
               	C4s = stCube[1][nj][nk-1].sg;
			   }
               else
			   {
                C4e = stCube[ni-1][nj][nk-1].e;
                C4s = stCube[ni-1][nj][nk-1].sg;
			   }
		 }
		 // Compute the average conductivity
		 fsigma = (C1s+C2s+C3s+C4s)/4;
         // Compute the average dielectric constant
		 fe = (C1e+C2e+C3e+C4e)/4;
		 break;
		 
	case 'z':
	// For Ez: ez = (C1e+C2e+C3e+C4e)/4
	//     Ez: sigmaz = (C1s+C2s+C3s+C4s)/4
		 // for C1e and C1s:
         if (ni == (nNx+1))
		 {
			 // ni = nNx+1 
               if ((nj-1)==nNy)
			   {
               	C1e = stCube[nNx][nNy][nk].e;
               	C1s = stCube[nNx][nNy][nk].sg;
			   }
               else
			   {
                C1e = stCube[nNx][nj][nk].e;
                C1s = stCube[nNx][nj][nk].sg;
			   }
		 }
		 else
		 {               
			 // ni < nNx+1
               if ((nj-1)==nNy)
			   {
               	C1e = stCube[ni][nNy][nk].e;
               	C1s = stCube[ni][nNy][nk].sg;
			   }
               else
			   {
                C1e = stCube[ni][nj][nk].e;
                C1s = stCube[ni][nj][nk].sg;
			   }
		 }   
         
         // for C2e and C2s:
         if (ni == 1)
		 {
			 // ni = 1 
               if ((nj-1)==nNy)
			   {
               	C2e = stCube[1][nNy][nk].e;
              	C2s = stCube[1][nNy][nk].sg;
			   }
               else
			   {
                C2e = stCube[1][nj][nk].e;
                C2s = stCube[1][nj][nk].sg;
			   }
		 }
		 else
		 {
             // ni > 1
               if ((nj-1)==nNy)
			   {
               	C2e = stCube[ni-1][nNy][nk].e;
               	C2s = stCube[ni-1][nNy][nk].sg;
			   }
               else
			   {
                C2e = stCube[ni-1][nj][nk].e;
                C2s = stCube[ni-1][nj][nk].sg;
			   }
		 }
         
         // for C3e and C3s:
         if (ni == 1)
		 {
			 // ni = 1 
               if ((nj-1)==0)
			   {
               	C3e = stCube[1][1][nk].e;
               	C3s = stCube[1][1][nk].sg;
			   }
               else
			   {
                C3e = stCube[1][nj-1][nk].e;
                C3s = stCube[1][nj-1][nk].sg;
			   }
		 }
         else
		 {
             // ni > 1
               if ((nj-1)==0)
			   {
               	C3e = stCube[ni-1][1][nk].e;
               	C3s = stCube[ni-1][1][nk].sg;
			   }
               else
			   {
                C3e = stCube[ni-1][nj-1][nk].e;
                C3s = stCube[ni-1][nj-1][nk].sg;
			   }
		 }
         
         // for C4e and C4s:        
         if (ni == (nNx+1))
		 {
			 // ni = nNx+1 
               if ((nj-1)==0)
			   {
               	C4e = stCube[nNx][1][nk].e;
               	C4s = stCube[nNx][1][nk].sg;
			   }
               else
			   {
                C4e = stCube[nNx][nj-1][nk].e;
                C4s = stCube[nNx][nj-1][nk].sg;
			   }
		 }
		 else
		 {
			 // ni < nNx+1
               if ((nj-1)==0)
			   {
               	C4e = stCube[ni][1][nk].e;
               	C4s = stCube[ni][1][nk].sg;
			   }
               else
			   {
                C4e = stCube[ni][nj-1][nk].e;
                C4s = stCube[ni][nj-1][nk].sg;
			   }
		 }
		 
		 // Compute the average conductivity
		 fsigma = (C1s+C2s+C3s+C4s)/4;
         // Compute the average dielectric constant
		 fe = (C1e+C2e+C3e+C4e)/4;
		 break;

	default:
		fe = 1;
		fsigma = 0;
	}
}

// function to assign computation mode to each electric field component based on cube
// characteristics.
// Local arguments:
// nNx = no. of cubes along x axis
// nNy = no. of cubes along y axis
// nNz = no. of cubes along z axis
// return: void
// If successful: nStatus = NORMAL.
// If not successfull: nStatus = OUT_OF_BOUND.
void tCube::AssignComMode(tEHfdtd &teh)
{
	int ni, nj, nk;

	// check for size error
	if (nNx > MAXX-1)
	{
		nStatus = OUT_OF_BOUND;
		return;
	}
	if (nNy > MAXY-1) 
	{
		nStatus = OUT_OF_BOUND;
		return;
	}
	if (nNz > MAXZ-1)
	{
		nStatus = OUT_OF_BOUND;
		return;
	}

for (ni=1; ni<=nNx; ni++)
{
	for (nj=1; nj<=nNy; nj++)
	{
		for (nk=1; nk<=nNz; nk++)
		{
			// Set the electric field and magnetic field computation mode based on cube
			// dielectric characteristic.
			switch (stCube[ni][nj][nk].d)
			{
			case (1): // Linear lossy dielectric, no dispersion.
			// Here we need to check whether the E field component is already a conductor, this is to
			// prevent the assignment of lossy dielectric update mode to overwrite the PEC update mode
			// which is set by adjacent cell.  Discovered this error on 11th June 2002.
				if (teh.GetExMode(ni,nj,nk)!=1) // if not PEC
					teh.SetExMode(2,ni,nj,nk);
				if (teh.GetExMode(ni,nj+1,nk)!=1) // if not PEC
					teh.SetExMode(2,ni,nj+1,nk);
				if (teh.GetExMode(ni,nj,nk+1)!=1) // if not PEC
					teh.SetExMode(2,ni,nj,nk+1);
				if (teh.GetExMode(ni,nj+1,nk+1)!=1) // if not PEC
					teh.SetExMode(2,ni,nj+1,nk+1);
				if (teh.GetEyMode(ni,nj,nk)!=1) // if not PEC
					teh.SetEyMode(2,ni,nj,nk);
				if (teh.GetEyMode(ni+1,nj,nk)!=1) // if not PEC
					teh.SetEyMode(2,ni+1,nj,nk);
				if (teh.GetEyMode(ni,nj,nk+1)!=1) // if not PEC
					teh.SetEyMode(2,ni,nj,nk+1);
				if (teh.GetEyMode(ni+1,nj,nk+1)!=1) // if not PEC	
					teh.SetEyMode(2,ni+1,nj,nk+1);
				if (teh.GetEzMode(ni+1,nj,nk)!=1) // if not PEC
					teh.SetEzMode(2,ni+1,nj,nk);
				if (teh.GetEzMode(ni+1,nj+1,nk)!=1) // if not PEC
					teh.SetEzMode(2,ni+1,nj+1,nk);
				if (teh.GetEzMode(ni,nj+1,nk)!=1) // if not PEC
					teh.SetEzMode(2,ni,nj+1,nk);
				if (teh.GetEzMode(ni,nj,nk)!=1) // if not PEC
					teh.SetEzMode(2,ni,nj,nk);				
				break;
			default: // Cube is linear lossless dielectric, no dispersion.
					 // (computation mode = 0)
				break;
			}
			// Set the electric field and magnetic field computation mode based on cube
			// surface characteristics.
			switch (stCube[ni][nj][nk].n)
			{
			case (1): // PEC on face 1 of cube (i,j,k)
			// When mode = 1, the e field component will be
			// automatically set to 0.
				teh.SetEyMode(1,ni+1,nj,nk);
				teh.SetEzMode(1,ni+1,nj,nk);
				teh.SetEzMode(1,ni+1,nj+1,nk);
				teh.SetEyMode(1,ni+1,nj,nk+1);
				teh.SetHxMode(1,ni+1,nj,nk);
				break;

			case (2): // PEC on face 2 of cube (i,j,k)
				teh.SetExMode(1,ni,nj+1,nk);
				teh.SetEzMode(1,ni+1,nj+1,nk);
				teh.SetEzMode(1,ni,nj+1,nk);
				teh.SetExMode(1,ni,nj+1,nk+1);
				teh.SetHyParam(1,ni,nj+1,nk);
				break;
			
			case (3): // PEC on face 3 of cube (i,j,k)
				teh.SetEyMode(1,ni,nj,nk);
				teh.SetEzMode(1,ni,nj+1,nk);
				teh.SetEzMode(1,ni,nj,nk);
				teh.SetEyMode(1,ni,nj,nk+1);
				teh.SetHxMode(1,ni,nj,nk);
				break;

			case (4): // PEC on face 4 of cube (i,j,k)
				teh.SetExMode(1,ni,nj,nk);
				teh.SetEzMode(1,ni,nj,nk);
				teh.SetEzMode(1,ni+1,nj,nk);
				teh.SetExMode(1,ni,nj,nk+1);
				teh.SetHyMode(1,ni,nj,nk);
				break;

			case (5): // PEC on face 5 of cube (i,j,k)
				teh.SetEyMode(1,ni+1,nj,nk);
				teh.SetExMode(1,ni,nj+1,nk);
				teh.SetEyMode(1,ni,nj,nk);
				teh.SetExMode(1,ni,nj,nk);
				teh.SetHzMode(1,ni,nj,nk);
				break;

			case (6): // PEC on face 6 of cube (i,j,k)
				teh.SetEyMode(1,ni+1,nj,nk+1);
				teh.SetExMode(1,ni,nj+1,nk+1);
				teh.SetEyMode(1,ni,nj,nk+1);
				teh.SetExMode(1,ni,nj,nk+1);
				teh.SetHzMode(1,ni,nj,nk+1);
				break;

			case (7): // Upper triangular PEC on face 1 
					  // of cube (i,j,k)
				teh.SetEyMode(1,ni+1,nj,nk+1);
				teh.SetEzMode(1,ni+1,nj,nk);
				teh.SetHxMode(2,ni+1,nj,nk);
				break;

			case (8): // Lower triangular PEC on face 1 
					  // of cube (i,j,k)
				teh.SetEyMode(1,ni+1,nj,nk);
				teh.SetEzMode(1,ni+1,nj+1,nk);
				teh.SetHxMode(3,ni+1,nj,nk);
				break;

			case (9): // Upper triangular PEC on face 2 
					  // of cube (i,j,k)
				teh.SetExMode(1,ni,nj+1,nk+1);
				teh.SetEzMode(1,ni+1,nj+1,nk);
				teh.SetHyMode(2,ni,nj+1,nk);
				break;

			case (10): // Lower triangular PEC on face 2 
					  // of cube (i,j,k)
				teh.SetExMode(1,ni,nj+1,nk);
				teh.SetEzMode(1,ni,nj+1,nk);
				teh.SetHyMode(3,ni,nj+1,nk);
				break;

			case (11): // Upper triangular PEC on face 3 
					  // of cube (i,j,k)
				teh.SetEyMode(1,ni,nj,nk);
				teh.SetEzMode(1,ni,nj+1,nk);
				teh.SetHxMode(2,ni,nj,nk);
				break;

			case (12): // Lower triangular PEC on face 3 
					  // of cube (i,j,k)
				teh.SetEyMode(1,ni,nj,nk+1);
				teh.SetEzMode(1,ni,nj,nk);
				teh.SetHxMode(3,ni,nj,nk);
				break;

			case (13): // Upper triangular PEC on face 4 
					  // of cube (i,j,k)
				teh.SetExMode(1,ni,nj,nk+1);
				teh.SetEzMode(1,ni,nj,nk);
				teh.SetHyMode(2,ni,nj,nk);
				break;

			case (14): // Lower triangular PEC on face 4 
					  // of cube (i,j,k)
				teh.SetExMode(1,ni,nj,nk);
				teh.SetEzMode(1,ni+1,nj,nk);
				teh.SetHyMode(3,ni,nj,nk);
				break;

			case (15): // Upper triangular PEC on face 5 
					  // of cube (i,j,k)
				teh.SetExMode(1,ni,nj,nk);
				teh.SetEyMode(1,ni+1,nj,nk);
				teh.SetHzMode(2,ni,nj,nk);
				break;

			case (16): // Lower triangular PEC on face 5 
					  // of cube (i,j,k)
				teh.SetExMode(1,ni,nj+1,nk);
				teh.SetEyMode(1,ni,nj,nk);
				teh.SetHzMode(3,ni,nj,nk);
				break;

			case (17): // Upper triangular PEC on face 6 
					  // of cube (i,j,k)
				teh.SetExMode(1,ni,nj,nk+1);
				teh.SetEyMode(1,ni,nj,nk+1);
				teh.SetHzMode(2,ni,nj,nk+1);
				break;

			case (18): // Lower triangular PEC on face 6 
					  // of cube (i,j,k)
				teh.SetExMode(1,ni,nj+1,nk+1);
				teh.SetEyMode(1,ni+1,nj,nk+1);
				teh.SetHzMode(3,ni,nj,nk+1);
				break;
			
			case (19): // Cross lower triangular PEC on face 6
					   // of cube (i,j,k)
				teh.SetExMode(1,ni,nj,nk+1);
				teh.SetEyMode(1,ni+1,nj,nk+1);
				teh.SetHzMode(4,ni,nj,nk+1);				
				break;

			case (20): // Cross upper triangular PEC on face 6 
					   // of cube (i,j,k)
				teh.SetExMode(1,ni,nj+1,nk+1);
				teh.SetEyMode(1,ni,nj,nk+1);
				teh.SetHzMode(5,ni,nj,nk+1);				
				break;

			case (21): // whole cube is covered with PEC
				teh.SetExMode(1,ni,nj,nk);
				teh.SetExMode(1,ni,nj+1,nk);
				teh.SetExMode(1,ni,nj,nk+1);
				teh.SetExMode(1,ni,nj+1,nk+1);
				teh.SetEyMode(1,ni,nj,nk);
				teh.SetEyMode(1,ni+1,nj,nk);
				teh.SetEyMode(1,ni,nj,nk+1);
				teh.SetEyMode(1,ni+1,nj,nk+1);
				teh.SetEzMode(1,ni+1,nj,nk);
				teh.SetEzMode(1,ni+1,nj+1,nk);
				teh.SetEzMode(1,ni,nj+1,nk);
				teh.SetEzMode(1,ni,nj,nk);
				teh.SetHxMode(1,ni,nj,nk);
				teh.SetHxMode(1,ni+1,nj,nk);
				teh.SetHyMode(1,ni,nj,nk);
				teh.SetHyMode(1,ni,nj+1,nk);
				teh.SetHzMode(1,ni,nj,nk);
				teh.SetHzMode(1,ni,nj,nk+1);
				break;

			default: // Cube is dielectric, dielectric constant = 1.
				break;
			}
		}
	}
}
}

// function to assign computation parameter index to each electric field component
// based on lumped component specification.
// Local arguments:
// teh = pointer the class object tEHfdtd, this object contains the look-up
// table.
// return: void
// If successful: nStatus = NORMAL.
// If not successfull: nStatus = OUT_OF_BOUND.
void tCube::AssignComParam(tEHfdtd &teh)
{
	int ni, nj, nk;
	int nmatch = 0;
	unsigned char uEi=0, uHi=0;
	float fep, fsg;

	struct stEauxiliary // Auxiliary look-up table row structure
	{					// for electric field
		float fe;		// dielectric constant
		float fs;		// conductivity
		unsigned char  umode;	// update mode
	};

	struct stHauxiliary // Auxiliary look-up table row structure
	{					// for magnetic field
		unsigned char umode;	// update mode
	};

	stEauxiliary *stCeax;
	stEauxiliary *stCeay;
	stEauxiliary *stCeaz;
	stHauxiliary *stChax;
	stHauxiliary *stChay;
	stHauxiliary *stChaz;

	// Initialize auxiliary look-up table for E field
	stCeax = new stEauxiliary[MAXECONST];
	if (stCeax == NULL)
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}

	stCeay = new stEauxiliary[MAXECONST];
	if (stCeay == NULL)
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}

	stCeaz = new stEauxiliary[MAXECONST];
	if (stCeaz == NULL)
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}

	// Initialize auxiliary look-up table for H field
	stChax = new stHauxiliary[MAXHCONST];
	if (stChax == NULL)
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}

	stChay = new stHauxiliary[MAXHCONST];
	if (stChay == NULL)
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}

	stChaz = new stHauxiliary[MAXHCONST];
	if (stChaz == NULL)
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}

	// check for size error
	if (nNx > MAXX-1)
	{
		nStatus = OUT_OF_BOUND;
		return;
	}
	if (nNy > MAXY-1) 
	{
		nStatus = OUT_OF_BOUND;
		return;
	}
	if (nNz > MAXZ-1)
	{
		nStatus = OUT_OF_BOUND;
		return;
	}


// --- For Ex ---: 
// Begin to create constants for 1st valid cell
	AssignEpsilonSigma('x',1,1,1,fep,fsg); // assign dielectric constant and conductivity
	CreateExConstEntry(teh,uEimaxx,teh.GetExMode(1,1,1),fep,fsg,1,1,1,NULL);
	// update auxiliary look-up table for E field
	stCeax[uEimaxx-1].fe = fep;
	stCeax[uEimaxx-1].fs = fsg;
	stCeax[uEimaxx-1].umode = teh.GetExMode(1,1,1);

// Create entry when field component in PEC
	CreateExConstEntry(teh,uEimaxx,1,fep,fsg,0,0,0,NULL); 
					// the values fep and fsg are just a dummies
	// update auxiliary look-up table for E field
	stCeax[uEimaxx-1].fe = fep;
	stCeax[uEimaxx-1].fs = fsg;
	stCeax[uEimaxx-1].umode = 1;

for (ni=1; ni<=nNx; ni++)
{
	for (nj=1; nj<=nNy+1; nj++)
	{
		for (nk=1; nk<=nNz+1; nk++)
		{
			nmatch = 0;
			AssignEpsilonSigma('x',ni,nj,nk,fep,fsg); // get dielectric constant, conductivity
			// search for match entry
			for (uEi=0; uEi<uEimaxx; uEi++)
			{
				if ((fep==stCeax[uEi].fe) &&
					(fsg==stCeax[uEi].fs) &&
					(teh.GetExMode(ni,nj,nk)==stCeax[uEi].umode))
				{	// dielectric constant and mode match
					teh.SetExParam(uEi,ni,nj,nk);
					nmatch = 1; // indicate find match
					break;
				}
				else if (teh.GetExMode(ni,nj,nk) == 1) // 1 is mode for PEC
				{ // field in PEC, dielectric is not important
					teh.SetExParam(1,ni,nj,nk);
					nmatch = 1; // indicate find match
					break;
				}
			}

			if (nmatch == 0) // no match, create new entry. 
			{
				if (CreateExConstEntry(teh,uEimaxx,teh.GetExMode(ni,nj,nk),fep,fsg,
					ni,nj,nk,NULL)==TRUE)
				{
					teh.SetExParam(uEimaxx-1,ni,nj,nk); // assign Ex to this entry.
					// update auxiliary look-up table E field
					stCeax[uEimaxx-1].fe = fep;
					stCeax[uEimaxx-1].fs = fsg;
					stCeax[uEimaxx-1].umode = teh.GetExMode(ni,nj,nk);
				}
				else
					nStatus = OUT_OF_BOUND;
			}
		}
	}
}

// --- For Ey ---: 
// Begin to create constants for 1st valid cell
	AssignEpsilonSigma('y',1,1,1,fep,fsg); // assign dielectric constant and conductivity 
	CreateEyConstEntry(teh,uEimaxy,teh.GetEyMode(1,1,1),fep,fsg,1,1,1,NULL);
	// update auxiliary look-up table for E field
	stCeay[uEimaxy-1].fe = fep;
	stCeay[uEimaxy-1].fs = fsg;
	stCeay[uEimaxy-1].umode = teh.GetEyMode(1,1,1);

// Create entry when field component in PEC
	CreateEyConstEntry(teh,uEimaxy,1,fep,fsg,0,0,0,NULL); 
					// Note: fep and fsg are just dummy values.
	// update auxiliary look-up table for E field
	stCeay[uEimaxy-1].fe = fep;
	stCeay[uEimaxy-1].fs = fsg;
	stCeay[uEimaxy-1].umode = 1;

for (ni=1; ni<=nNx+1; ni++)
{
	for (nj=1; nj<=nNy; nj++)
	{
		for (nk=1; nk<=nNz+1; nk++)
		{
			nmatch = 0;
			AssignEpsilonSigma('y',ni,nj,nk,fep,fsg); // get dielectric constant, conductivity
			// search for match entry
			for (uEi=0; uEi<uEimaxy; uEi++)
			{
				if ((fep==stCeay[uEi].fe) &&
					(fsg==stCeay[uEi].fs) &&
					(teh.GetEyMode(ni,nj,nk)==stCeay[uEi].umode))
				{	// dielectric constant and mode match
					teh.SetEyParam(uEi,ni,nj,nk);
					nmatch = 1; // indicate find match
					break;
				}
				else if (teh.GetEyMode(ni,nj,nk) == 1) // 1 is mode for PEC
				{ // field in PEC, dielectric is not important
					teh.SetEyParam(1,ni,nj,nk);
					nmatch = 1; // indicate find match
					break;
				}
			}

			if (nmatch == 0) // no match, create new entry. 
			{
				if (CreateEyConstEntry(teh,uEimaxy,teh.GetEyMode(ni,nj,nk),
					fep,fsg,ni,nj,nk,NULL)==TRUE)
				{
					teh.SetEyParam(uEimaxy-1,ni,nj,nk); // assign Ey to this entry.
					// update auxiliary look-up table E field
					stCeay[uEimaxy-1].fe = fep;
					stCeay[uEimaxy-1].fs = fsg;
					stCeay[uEimaxy-1].umode = teh.GetEyMode(ni,nj,nk);
				}
				else
					nStatus = OUT_OF_BOUND;
			}
		}
	}
}

// --- For Ez ---:
// Begin to create constants for 1st valid cell
	AssignEpsilonSigma('z',1,1,1,fep,fsg); // assign dielectric constant and conductivity
	CreateEzConstEntry(teh,uEimaxz,teh.GetEzMode(1,1,1),fep,fsg,1,1,1,NULL);
	// update auxiliary look-up table for E field
	stCeaz[uEimaxz-1].fe = fep;
	stCeaz[uEimaxz-1].fs = fsg;
	stCeaz[uEimaxz-1].umode = teh.GetEzMode(1,1,1);

// Create entry when field component in PEC
	CreateEzConstEntry(teh,uEimaxz,1,fep,fsg,0,0,0,NULL); 
					// Note: fep and fsg are just a dummy values.
	// update auxiliary look-up table for E field
	stCeaz[uEimaxz-1].fe = fep;
	stCeaz[uEimaxz-1].fs = fsg;
	stCeaz[uEimaxz-1].umode = 1;

for (ni=1; ni<=nNx+1; ni++)
{
	for (nj=1; nj<=nNy+1; nj++)
	{
		for (nk=1; nk<=nNz; nk++)
		{
			nmatch = 0;
			AssignEpsilonSigma('z',ni,nj,nk,fep,fsg); // get dielectric constant, conductivity
			// search for match entry
			for (uEi=0; uEi<uEimaxz; uEi++)
			{
				if ((fep==stCeaz[uEi].fe) &&
					(fsg==stCeaz[uEi].fs) &&
					(teh.GetEzMode(ni,nj,nk)==stCeaz[uEi].umode))
				{	// dielectric constant and mode match
					teh.SetEzParam(uEi,ni,nj,nk);
					nmatch = 1; // indicate find match
					break;
				}
				else if (teh.GetEzMode(ni,nj,nk) == 1) // 1 is mode for PEC
				{ // field in PEC, dielectric is not important
					teh.SetEzParam(1,ni,nj,nk);
					nmatch = 1; // indicate find match
					break;
				}
			}

			if (nmatch == 0) // no match, create new entry. 
			{			
				if (CreateEzConstEntry(teh, uEimaxz,teh.GetEzMode(ni,nj,nk),fep,fsg,
					ni,nj,nk,NULL)==TRUE)
				{
					teh.SetEzParam(uEimaxz-1,ni,nj,nk); // assign Ex to this entry.
					// update auxiliary look-up table E field
					stCeaz[uEimaxz-1].fe = fep;
					stCeaz[uEimaxz-1].fs = fsg;
					stCeaz[uEimaxz-1].umode = teh.GetEzMode(ni,nj,nk);
				}
				else
					nStatus = OUT_OF_BOUND;
			}
		}
	}
}

// --- for Hx ---:
// Create the 1st entry
	CreateHxConstEntry(teh, uHimaxx, teh.GetHxMode(1,1,1),
						ni, nj, nk); // Create a new entry
// Update the auxiliary look-up table for H field
	stChax[uHimaxx-1].umode = teh.GetHxMode(1,1,1);

// Create the 2nd entry, when field component normal to PEC
	CreateHxConstEntry(teh, uHimaxx, 1, 0, 0, 0); // Create a new entry
// Update the auxiliary look-up table for H field
	stChax[uHimaxx-1].umode = 1;

for (ni=1; ni<=nNx+1; ni++)
{
	for (nj=1; nj<=nNy; nj++)
	{
		for (nk=1; nk<=nNz; nk++)
		{
			nmatch = 0;

			// search for match entry
			for (uHi=0; uHi<uHimaxx; uHi++)
			{
				if (teh.GetHxMode(ni,nj,nk)==stChax[uHi].umode)
				{ // computation mode match
					teh.SetHxParam(uHi,ni,nj,nk); // assign Hx to this entry.
					nmatch = 1; // indicate find match
					break;
				}
			}

			if (nmatch == 0) // no match, create new entry
			{
				if (CreateHxConstEntry(teh,uHimaxx,teh.GetHxMode(ni,nj,nk),
					ni,nj,nk)==TRUE)
				{
					teh.SetHxParam(uHimaxx-1,ni,nj,nk); // assign Hx to this entry.
					// update auxiliary look-up table H field
					stChax[uHimaxx-1].umode = teh.GetHxMode(ni,nj,nk);
				}
				else
					nStatus = OUT_OF_BOUND;
			}
		}
	}

}

// --- for Hy ---:
// Create the 1st entry
	CreateHyConstEntry(teh, uHimaxy, teh.GetHyMode(1,1,1),
						ni, nj, nk); // Create a new entry
// Update the auxiliary look-up table for H field
	stChay[uHimaxy-1].umode = teh.GetHyMode(1,1,1);

// Create the 2nd entry, when field component normal to PEC
	CreateHyConstEntry(teh, uHimaxy, 1, 0, 0, 0); // Create a new entry
// Update the auxiliary look-up table for H field
	stChay[uHimaxy-1].umode = 1;

for (ni=1; ni<=nNx; ni++)
{
	for (nj=1; nj<=nNy+1; nj++)
	{
		for (nk=1; nk<=nNz; nk++)
		{
			nmatch = 0;

			// search for match entry
			for (uHi=0; uHi<uHimaxy; uHi++)
			{
				if (teh.GetHyMode(ni,nj,nk)==stChay[uHi].umode)
				{ // computation mode match
					teh.SetHyParam(uHi,ni,nj,nk); // assign Hy to this entry.
					nmatch = 1; // indicate find match
					break;
				}
			}

			if (nmatch == 0) // no match, create new entry
			{
				if (CreateHyConstEntry(teh,uHimaxy,teh.GetHyMode(ni,nj,nk),
					ni,nj,nk)==TRUE)
				{
					teh.SetHyParam(uHimaxy-1,ni,nj,nk); // assign Hy to this entry.
					// update auxiliary look-up table H field
					stChay[uHimaxy-1].umode = teh.GetHyMode(ni,nj,nk);
				}
				else
					nStatus = OUT_OF_BOUND;
			}				
		}
	}
}

// --- for Hz ---:
// Create the 1st entry
	CreateHzConstEntry(teh, uHimaxz, teh.GetHzMode(1,1,1),
						ni, nj, nk); // Create a new entry
// Update the auxiliary look-up table for H field
	stChaz[uHimaxz-1].umode = teh.GetHzMode(1,1,1);

// Create the 2nd entry, when field component normal to PEC
	CreateHzConstEntry(teh, uHimaxz, 1, 0, 0, 0); // Create a new entry
// Update the auxiliary look-up table for H field
	stChaz[uHimaxz-1].umode = 1;

for (ni=1; ni<=nNx; ni++)
{
	for (nj=1; nj<=nNy; nj++)
	{
		for (nk=1; nk<=nNz+1; nk++)
		{
			nmatch = 0;

			// search for match entry
			for (uHi=0; uHi<uHimaxz; uHi++)
			{
				if (teh.GetHzMode(ni,nj,nk)==stChaz[uHi].umode)
				{ // computation mode match
					teh.SetHzParam(uHi,ni,nj,nk); // assign Hz to this entry.
					nmatch = 1; // indicate find match
					break;
				}
			}

			if (nmatch == 0) // no match, create new entry
			{
				if (CreateHzConstEntry(teh,uHimaxz,teh.GetHzMode(ni,nj,nk),
					ni,nj,nk)==TRUE)
				{
					teh.SetHzParam(uHimaxz-1,ni,nj,nk); // assign Hz to this entry.
					// update auxiliary look-up table H field
					stChaz[uHimaxz-1].umode = teh.GetHzMode(ni,nj,nk);
				}
				else
					nStatus = OUT_OF_BOUND;
			}	
		}
	}
}

	delete[] stCeax; // remove all auxilliary look-up table.
	delete[] stCeay;
	delete[] stCeaz;
	delete[] stChax;
	delete[] stChay;
	delete[] stChaz;
}


// function to assign computation parameter index to each electric 
// components
// arguments:
// teh = pointer the class object tEHfdtd, this object contains the look-up table.
// return: void
// If successful: nStatus = NORMAL.
// If not successfull: nStatus = OUT_OF_BOUND.
void tCube::AssignCompParam(tEHfdtd &teh)
{
	unsigned char uEi, uHi=0;
	unsigned char uEmaxx=0, uEmaxy=0, uEmaxz=0;
	unsigned char uEstx, uHstx; // temporary variables to store the starting
	unsigned char uEsty, uHsty;	// index in the computational constants look-
	unsigned char uEstz, uHstz;	// up table for components
	lpCOMP lpcompcur=NULL,*refCaux,*refCauy,*refCauz;
	uEstx = uEimaxx;
	uEsty = uEimaxy;
	uEstz = uEimaxz;
	int nmatch = 0;
	float fep, fsg; // temporary variable for storing dielectric constant and conductivity

	refCaux = new lpCOMP[MAXECONST];
	refCauy = new lpCOMP[MAXECONST];
	refCauz = new lpCOMP[MAXECONST];
	if (refCaux == NULL)
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if (refCauy == NULL)
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}
	if (refCauz == NULL)
	{
		nStatus = OUT_OF_MEMORY;
		return;
	}

	if (lpcomp1 != NULL)
	{
		lpcompcur = lpcomp1;
		refCaux[0] = NULL;	// initialize auxiliary table
		refCauy[0] = NULL;  
		refCauz[0] = NULL; 
	}

	while (lpcompcur != NULL)
	{	// compare field type
		if (strcmp("ex",lpcompcur->cType)==0) // Ex field
		{	
			nmatch = 0; // reset match indicator
			if (lpcompcur->uType < 180) // only execute comparison routines if the
										// component is memoryless type.  Otherwise
										// immediately create a new entry in the 
										// computation coefficient table.
			{
				// check for match entries, memoryless Ex component
				for (uEi=0; uEi<uEmaxx; uEi++)
				{	
					// check if the auxiliary table is empty
					if (uEi == 0)
					{
						if (refCaux[uEi] == NULL)
							break;
					}
					// compare component type
					if ((refCaux[uEi]->uType)==(lpcompcur->uType))
					{
						// compare parameter list pointer
						if (CompareParam(refCaux[uEi]->plist,lpcompcur->plist)==TRUE)
						{	// update parameter index for x,y or z component
							teh.SetExParam(uEi+uEstx,lpcompcur->nX,lpcompcur->nY,
										lpcompcur->nZ);
							teh.SetExMode(lpcompcur->uType,lpcompcur->nX,
										  lpcompcur->nY,lpcompcur->nZ);
							nmatch = 1; // indicate found match
						}
					}
				}
			}
			if (nmatch == 0) // no match
			{	
				AssignEpsilonSigma('x',lpcompcur->nX,lpcompcur->nY,lpcompcur->nZ,fep,fsg); 
												// get dielectric constant and conductivity
				if (CreateExConstEntry(teh,uEimaxx,lpcompcur->uType,fep,fsg,
					lpcompcur->nX,lpcompcur->nY,lpcompcur->nZ,
					lpcompcur->plist)==TRUE) // create new entry, if fail do not 
				{						 // proceed
					if (uEmaxx < uEimaxx-uEstx) // only increment local index for x and
					{							// assign new pointer for auxiliary compo- 
						uEmaxx++;				// nent pointer when new entry is inserted	
												// Note: uEimaxx = uEstx+uEmaxx, when 
												// uEimaxx is incremented within the method,
						refCaux[uEmaxx-1] = lpcompcur;	// the equality will not hold.
					}
				}
				else
					nStatus = OUT_OF_BOUND;
			}

		}
		else if (strcmp("ey",lpcompcur->cType)==0) // Ey field
		{	
			nmatch = 0; // reset match indicator
			if (lpcompcur->uType < 180) // only execute comparison routines if the
										// component is memoryless type.  Otherwise
										// immediately create a new entry in the 
										// computation coefficient table.
			{ 
			// check for match entries, memoryless Ey component
				for (uEi=0; uEi<uEmaxy; uEi++)
				{	
					// check if the auxiliary table is empty
					if (uEi == 0)
					{
						if (refCauy[uEi] == NULL)
							break;
					}				
					// compare component type
					if ((refCauy[uEi]->uType)==(lpcompcur->uType))
					{
						// compare parameter list 
						if (CompareParam(refCauy[uEi]->plist,lpcompcur->plist)==TRUE)
						{	// update parameter index for x,y or z component
							teh.SetEyParam(uEi+uEsty,lpcompcur->nX,lpcompcur->nY,
										lpcompcur->nZ);
							teh.SetEyMode(lpcompcur->uType,lpcompcur->nX,
										  lpcompcur->nY,lpcompcur->nZ);
							nmatch = 1; // indicate found match
						}
					}
				}
			}
			if (nmatch == 0) // no match
			{	
				AssignEpsilonSigma('y',lpcompcur->nX,lpcompcur->nY,lpcompcur->nZ,fep,fsg); 
												// get dielectric constant and conductivity
				if (CreateEyConstEntry(teh,uEimaxy,lpcompcur->uType,fep,fsg,
					lpcompcur->nX,lpcompcur->nY,lpcompcur->nZ,
					lpcompcur->plist)==TRUE) // create new entry, if fail do not 
				{						 // proceed
					if (uEmaxy < uEimaxy-uEsty) // only increment local index for y and
					{							// assign new pointer for auxiliary compo- 
						uEmaxy++;				// nent pointer when new entry is inserted	
												// Note: uEimaxy = uEsty+uEmaxy, when 
												// uEimaxy is incremented within the method,
						refCauy[uEmaxy-1] = lpcompcur;	// the equality will not hold.
					}
				}
				else
					nStatus = OUT_OF_BOUND;
			}

		}
		else if (strcmp("ez",lpcompcur->cType)==0) // Ez field
		{	
			nmatch = 0; // reset match indicator
			if (lpcompcur->uType < 180) // only execute comparison routines if the
										// component is memoryless type.  Otherwise
										// immediately create a new entry in the 
										// computation coefficient table.
			{ 
				// check for match entries, memoryless Ez component.
				for (uEi=0; uEi<uEmaxz; uEi++)
				{	
					// check if the auxiliary table is empty
					if (uEi == 0)
					{
						if (refCauz[uEi] == NULL)
							break;
					}				
					// compare component type
					if ((refCauz[uEi]->uType)==(lpcompcur->uType))
					{
						// compare parameter list
						if (CompareParam(refCauz[uEi]->plist,lpcompcur->plist)
							==TRUE)
						{	// update parameter index for x,y or z component
							teh.SetEzParam(uEi+uEstz,lpcompcur->nX,lpcompcur->nY,
										lpcompcur->nZ);
							teh.SetEzMode(lpcompcur->uType,lpcompcur->nX,
										  lpcompcur->nY,lpcompcur->nZ);
							nmatch = 1; // indicate found match
						}
					}
				}
			}
			if (nmatch == 0) // no match
			{	
				AssignEpsilonSigma('z',lpcompcur->nX,lpcompcur->nY,lpcompcur->nZ,fep,fsg); 
												// get dielectric constant and conductivity
				if (CreateEzConstEntry(teh,uEimaxz,lpcompcur->uType,fep,fsg,
					lpcompcur->nX,lpcompcur->nY,lpcompcur->nZ,
					lpcompcur->plist)==TRUE) // create new entry, if fail do not proceed
				{							
					if (uEmaxz < uEimaxz-uEstz) // only increment local index for z and
					{							// assign new pointer for auxiliary compo- 
						uEmaxz++;				// nent pointer when new entry is inserted	
												// Note: uEimaxz = uEstz+uEmaxz, when 
												// uEimaxz is incremented within the method,
						refCauz[uEmaxz-1] = lpcompcur;	// the equality will not hold.
					}
				}
				else
					nStatus = OUT_OF_BOUND;
			}

		}
		// point to next item in component list
		lpcompcur = lpcompcur->next;
	}

	delete[] refCaux; // remove all auxialary look up table
	delete[] refCauy;
	delete[] refCauz;
}


// Function to compare the parameters of two components
// arguments:
// plist1 = pointer to parameter list 1
// plist2 = pointer to parameter list 2
// return:
// TRUE if match and FALSE if not match.
BOOL tCube::CompareParam(lpDOUBLE plist1, lpDOUBLE plist2)
{
	int nmatch = 1; // match indicator, 0 = no match
					//					1 = match
	while ((plist1!=NULL) && (plist2!=NULL))
	{
		if ((plist1->fP) == (plist2->fP))
		{
			plist1 = plist1->next; // next item on plist1
			plist2 = plist2->next; // next item on plist2
		}
		else // on item not match, terminate comparison.
		{
			nmatch = 0;
			break;			
		}
	}

	if (nmatch == 1)
		return TRUE;
	else
		return FALSE;
}


// Function to create new entry in Ex field computation constants array.
// Each row of the array can store 10 constants: C0x, C1ex to C9x.
// These constants depends on the dielectric constant and characteristic
// of the cube.
// arguments:
// teh = Reference to class object tEHfdtd
// uExm = The last valid row index to array
// umode = The update mode for the field
// fe = The associated dielectric constant
// fsg = The associated conductivity
// i = x coordinate
// j = y coordinate
// k = z coordinate
// plist = pointer to parameter list
// return value: TRUE if success FALSE otherwise 
BOOL tCube::CreateExConstEntry(tEHfdtd &teh, unsigned char &uExm, unsigned char umode,
							   float fe, float fs, int i, int j, int k, lpDOUBLE plist)
{	
	double fC; // temporary variable
	unsigned char uor; // temporary variable, for storing orientation information

	// check for out of boundary 
	if (uExm > MAXECONST-1)
	{
		nStatus = OUT_OF_BOUND;
		return FALSE;
	}

	switch (umode)
	{
	case (0): // Normal update mode
		teh.stCCex[uExm].fC0 = fdel_t/(eo*fe*fdel_y); // C0x
		teh.stCCex[uExm].fC1 = fdel_t/(eo*fe*fdel_z); // C1x	
		break;
	case (1): // PEC
		teh.stCCex[uExm].fC0 = 0; // C0x
		teh.stCCex[uExm].fC1 = 0; // C1x	
		break;
	case (2): // Lossy dielectric, no dispersion.	
		fC = fs*fdel_t/2;
		teh.stCCex[uExm].fC0 = (eo*fe-fC)/(eo*fe+fC); // C0x
		teh.stCCex[uExm].fC1 = fdel_t/((eo*fe+fC)*fdel_y); // C1x
		teh.stCCex[uExm].fC2 = fdel_t/((eo*fe+fC)*fdel_z); // C2x
		break;
	case (50): // Resistive voltage source in x, periodic trapezoidal pulse
			   // syntax: comp X Y Z ez 50 R Vlow Vhigh tdly trise thigh tlow 
			   // note: all timings are in picoseconds
			   // default: R=50, Vlow=0, Vhigh=1, tdly=0ps, trise=1ns, thigh=3ns,
			   // tlow=5ns, amp=1.	
		fC = (fdel_t*fdel_x)/(2*fe*eo*fdel_y*fdel_z);
		if (plist == NULL)
			fC = fC/50.0;
		else
			fC = fC/(plist->fP);
		teh.stCCex[uExm].fC0 = (1-fC)/(1+fC);
		teh.stCCex[uExm].fC1 = fdel_t/(fdel_y*fe*eo*(1+fC));
		teh.stCCex[uExm].fC2 = fdel_t/(fdel_z*fe*eo*(1+fC));
		teh.stCCex[uExm].fC3 = (2*fC)/(fdel_x*(1+fC));
		
		// Assign default values
		teh.stCCex[uExm].fC4 = 0.0; // 0V
		teh.stCCex[uExm].fC5 = 1.0; // 1V
		teh.stCCex[uExm].fC6 = 0;   // tdly=0
		teh.stCCex[uExm].fC7 = 1000E-12; // trise/tfall=1ns
		teh.stCCex[uExm].fC8 = 3000E-12; // thigh=3ns
		teh.stCCex[uExm].fC9 = 5000E-12; // tlow=5ns

		// get other parameters
		if (plist == NULL)
			break;
		else
			plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCex[uExm].fC4 = plist->fP; // get Vlow
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCex[uExm].fC5 = plist->fP; // get Vhigh
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCex[uExm].fC6 = (plist->fP)*(1E-12); // get delay
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCex[uExm].fC7 = (plist->fP)*(1E-12); // get rise/fall time
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCex[uExm].fC8 = (plist->fP)*(1E-12); // get high time
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCex[uExm].fC9 = (plist->fP)*(1E-12); // get low time
		break;
	case (51): // Resistive voltage source in x, periodic trapezoidal pulse
			   // syntax: comp X Y Z ez 51 R phase freq atten amplitude
			   // note: all timings are in picoseconds
			   // default: R=50, phase=0rad, freq=1E06, atten=0, amp=1.	
		fC = (fdel_t*fdel_x)/(2*fe*eo*fdel_y*fdel_z);
		if (plist == NULL)
			fC = fC/50.0;
		else
			fC = fC/(plist->fP);
		teh.stCCex[uExm].fC0 = (1-fC)/(1+fC);
		teh.stCCex[uExm].fC1 = fdel_t/(fdel_y*fe*eo*(1+fC));
		teh.stCCex[uExm].fC2 = fdel_t/(fdel_z*fe*eo*(1+fC));
		teh.stCCex[uExm].fC3 = (2*fC)/(fdel_x*(1+fC));
		
		// Assign default values
		teh.stCCex[uExm].fC4 = 0;	// 0 radian
		teh.stCCex[uExm].fC5 = 1E06; // 1MHz
		teh.stCCex[uExm].fC6 = 0; // 0 Neper
		teh.stCCex[uExm].fC7 = 1; // 1 Volt

		// get other parameters
		if (plist == NULL)
			break;
		else
			plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCex[uExm].fC4 = plist->fP; // get phase
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCex[uExm].fC5 = plist->fP; // get frequency
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCex[uExm].fC6 = plist->fP; // get attenuation
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCex[uExm].fC7 = plist->fP; // get amplitude
		break;
	case (52): // x oriented DC source
			   // syntax: comp X Y Z ez 52 R Vdc
			   // default: R=50, Vdc=1.0
		fC = (fdel_t*fdel_x)/(2*fe*eo*fdel_y*fdel_z);
		if (plist == NULL)
			fC = fC/50.0;
		else
			fC = fC/(plist->fP);
		teh.stCCex[uExm].fC0 = (1-fC)/(1+fC);
		teh.stCCex[uExm].fC1 = fdel_t/(fdel_y*fe*eo*(1+fC));
		teh.stCCex[uExm].fC2 = fdel_t/(fdel_z*fe*eo*(1+fC));
		teh.stCCex[uExm].fC3 = (2*fC)/(fdel_x*(1+fC));
		// Assign default values
		teh.stCCex[uExm].fC4 = 1;	// 1V

		// Get other parameters
		if (plist == NULL)
			break;
		else
		{
			plist = plist->next;
			teh.stCCex[uExm].fC4 = plist->fP; // get dc voltage
		}
		break;
	case (100): // X oriented resistor
				// syntax: comp X Y Z ex 100 R
				// default: R=50.
		if (plist == NULL)
			fC = 50;
		else
			fC = plist->fP;
		fC = (fdel_t*fdel_x)/(2*fe*eo*fC*fdel_y*fdel_z);
		teh.stCCex[uExm].fC0 = (1-fC)/(1+fC);
		teh.stCCex[uExm].fC1 = fdel_t/(fdel_y*fe*eo*(1+fC));
		teh.stCCex[uExm].fC2 = fdel_t/(fdel_z*fe*eo*(1+fC));
		teh.stCCex[uExm].fC3 = fC*(2*fe*eo/fdel_t); // store effective conductance
		break;
	case (101): // X oriented capacitor
				// syntax: comp X Y Z ex 101 C
				// default: C = 1pF
		if (plist == NULL)
			fC = (1.0E-12*fdel_x)/(fe*eo*fdel_y*fdel_z);
		else
			fC = ((plist->fP)*fdel_x)/(fe*eo*fdel_y*fdel_z);
		teh.stCCex[uExm].fC0 = fdel_t/(fe*eo*(1+fC)*fdel_y);
		teh.stCCex[uExm].fC1 = fdel_t/(fe*eo*(1+fC)*fdel_z);
		break;
	case (150): // X oriented diode.
				// syntax: comp X Y Z ex 150 Is n Vj M Cj0 TT FC R
				// default: Is = 0.1nA, n = 1, Vj=0.6, M=0.5, Cj0=0
				//			TT=0, FC=0.5, R=0.
				// FC is also limited to a maximum of 0.98, to prevent
				// divide by 0 error.
		unsigned char ucor; // temporary variable

		if (plist == NULL) // assign direction/orientation
			ucor = 0;	   // 0 for default	
		else if (plist->fP == 0)
			ucor = 0;	   // +ve x direction
		else
			ucor = 1;	   // -ve x direction
		
		teh.stCCex[uExm].plist = new tEDiode(fdel_x,fdel_y,fdel_z,fdel_t,
									fe,ucor,plist); // create a new object for diode
		teh.stCCex[uExm].umode = 150; // set computation mode
		if (teh.stCCex[uExm].plist == NULL)
		{
			teh.nEStatus = OUT_OF_MEMORY; // error
		}
		break;
	case (180): // X oriented NPN transistor, assignment of BE junction
				// syntax: comp X Y Z ez 180 O Is Bf Br nF nR Va Ikf Ise ne Vb Ikr Isc nc 
			    // Vje Me Cje TTf Vjc Mc Cjc TTr FC
				// default: O=0, Is=0.1nA, Bf=100, Br=1, nF=1, nR=1, Va=500, Ikf=100, 
			    // Ise=0, ne=1.5, Vb=500, Ikr=100, Isc=0, nc=2, Vje=Vjc=0.6, Me=Mc=0.5, 
				// Cje=Cjc=0, TTf=TTr=0, FC=0.5
		tEBjt * ptrBjt1; // temporary pointer to BJT object

		if (plist == NULL) // assign direction/orientation
			uor = 0;	   // 0 for default	
		else if (plist->fP == 0)
			uor = 0;	   // +ve x direction
		else
			uor = 1;	   // -ve x direction

		teh.stCCex[uExm].plist = new tEBjt;
		teh.stCCex[uExm].fC0 = uor;	// 0 for +ve x, 1 for -ve x. 
		teh.stCCex[uExm].umode = 180; 
		ptObjtrack = teh.stCCex[uExm].plist;	// track current object
		uIndtrack = uExm;						// track current index
		ptrBjt1 = (tEBjt *) ptObjtrack;
		ptrBjt1->AssignBEjunction(i,j,k,fdel_x,fdel_y,fdel_z,fdel_t,fe,uor,plist); 
												// assign BE junction
		if (teh.stCCex[uExm].plist == NULL)
		{
			teh.nEStatus = OUT_OF_MEMORY; // error
			break;
		}
		if (nSimMode == REUSE) // if simulation mode is REUSE initialize internal BE field
		{
			ptrBjt1->Set_BEfield(teh.fEx[i][j][k]);
		}
		break;
	case (181):	// X oriented NPN transistor, assignment of BC junction
				// syntax: comp X Y Z ex 181 O
				// note: this statement must always follows the comp X Y Z e 180....
				//  	 statement.
		tEBjt * ptrBjt2; // temporary pointer to BJT object

		if (plist == NULL) // assign direction/orientation
			uor = 0;	   // 0 for default	
		else if (plist->fP == 0)
			uor = 0;	   // +ve x direction
		else
			uor = 1;	   // -ve x direction

		if (ptObjtrack != NULL) // if the object tracking pointer is properly assigned
		{
			ptrBjt2 = (tEBjt *) ptObjtrack;
			ptrBjt2->AssignBCjunction(i,j,k,fdel_x,fdel_y,fdel_z,fdel_t,fe,uor);
										// assign BC junction
			teh.stCCex[uExm].plist = ptrBjt2;
			teh.stCCex[uExm].umode = 181;
			teh.stCCex[uExm].fC0 = uor;	// 0 for +ve x, 1 for -ve x. 
			teh.SetExMode(180,i,j,k);	// fix computation mode the same as 180
			teh.SetExParam(uExm,i,j,k);	// fix index into computation coefficient structure
			uExm++;
			ptObjtrack = NULL; // reset object tracking pointer
			uIndtrack = 0; // reset index tracking variable
			if (nSimMode == REUSE) // if simulation mode is REUSE initialize internal BE field
			{
				ptrBjt2->Set_BCfield(teh.fEx[i][j][k]);
			}
			return TRUE;
		}
		else
			return FALSE;
	case (200): // X oriented inductor, a memory element.
				// syntax: comp X Y Z ex 200 L
				// default: L = 1nH
		if (plist == NULL)
			fC = 1.0E-9;
		else
			fC = plist->fP;
		teh.stCCex[uExm].fC2 = (fdel_x*fdel_t*fdel_t)/(fC*fe*eo*fdel_y*fdel_z); // C2x
		teh.stCCex[uExm].fC0 = fdel_t/(eo*fe*fdel_y); // C0x
		teh.stCCex[uExm].fC1 = fdel_t/(eo*fe*fdel_z); // C1x	
		teh.stCCex[uExm].fC3 = 0; // Initialize the memory element.
		break;
	default:
		return FALSE;
	}

	teh.SetExMode(umode,i,j,k);	// fix computation mode
	teh.SetExParam(uExm,i,j,k);	// fix index into computation coefficient structure
	uExm++; //increment max index
	return TRUE;
}

// Function to create new entry in Ey field computation constants array.
// Each row of the array can store 10 constants: C0y, C1y to C9y.
// These constants depends on the dielectric constant and characteristic
// of the cube.
// arguments:
// teh = Reference to class object tEHfdtd
// uEym = The last valid row index to array
// umode = The update mode for the field
// fe = The associated dielectric constant
// fsg = The associated conductivity
// i = x coordinate
// j = y coordinate
// k = z coordinate
// plist = pointer to parameter list
// return value: TRUE if success FALSE otherwise 
BOOL tCube::CreateEyConstEntry(tEHfdtd &teh, unsigned char &uEym, unsigned char umode,
							   float fe, float fs, int i, int j, int k, lpDOUBLE plist)
{	
	double fC; // temporary variable
	unsigned char uor; // temporary variable, for storing orientation information

	// check for out of boundary 
	if (uEym > MAXECONST-1)
	{
		nStatus = OUT_OF_BOUND;
		return FALSE;
	}

	switch (umode)
	{
	case (0): // Normal update mode	
		teh.stCCey[uEym].fC0 = fdel_t/(eo*fe*fdel_z); // C0y
		teh.stCCey[uEym].fC1 = fdel_t/(eo*fe*fdel_x); // C1y	
		break;
	case (1): // PEC
		teh.stCCey[uEym].fC0 = 0; // C0y
		teh.stCCey[uEym].fC1 = 0; // C1y
		break;
	case (2): // Lossy dielectric, no dispersion.
		fC = fs*fdel_t/2;
		teh.stCCey[uEym].fC0 = (eo*fe-fC)/(eo*fe+fC); // C0y
		teh.stCCey[uEym].fC1 = fdel_t/((eo*fe+fC)*fdel_z); // C1y
		teh.stCCey[uEym].fC2 = fdel_t/((eo*fe+fC)*fdel_x); // C2y	
		break;
	case (50): // Resistive voltage source in y, periodic trapezoidal pulse
			   // syntax: comp X Y Z ez 50 R Vlow Vhigh tdly trise thigh tlow 
			   // note: all timings are in picoseconds
			   // default: R=50, Vlow=0, Vhigh=1, tdly=0ps, trise=1ns, thigh=3ns,
			   // tlow=5ns, amp=1.	
		fC = (fdel_t*fdel_y)/(2*fe*eo*fdel_x*fdel_z);
		if (plist == NULL)
			fC = fC/50.0;
		else
			fC = fC/(plist->fP);
		teh.stCCey[uEym].fC0 = (1-fC)/(1+fC);
		teh.stCCey[uEym].fC1 = fdel_t/(fdel_z*fe*eo*(1+fC));
		teh.stCCey[uEym].fC2 = fdel_t/(fdel_x*fe*eo*(1+fC));
		teh.stCCey[uEym].fC3 = (2*fC)/(fdel_y*(1+fC));
		
		// Assign default values
		teh.stCCey[uEym].fC4 = 0.0; // 0V
		teh.stCCey[uEym].fC5 = 1.0; // 1V
		teh.stCCey[uEym].fC6 = 0;   // tdly=0
		teh.stCCey[uEym].fC7 = 1000E-12; // trise/tfall=1ns
		teh.stCCey[uEym].fC8 = 3000E-12; // thigh=3ns
		teh.stCCey[uEym].fC9 = 5000E-12; // tlow=5ns

		// get other parameters
		if (plist == NULL)
			break;
		else
			plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCey[uEym].fC4 = plist->fP; // get Vlow
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCey[uEym].fC5 = plist->fP; // get Vhigh
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCey[uEym].fC6 = (plist->fP)*(1E-12); // get delay
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCey[uEym].fC7 = (plist->fP)*(1E-12); // get rise/fall time
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCey[uEym].fC8 = (plist->fP)*(1E-12); // get high time
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCey[uEym].fC9 = (plist->fP)*(1E-12); // get low time
		break;
	case (51): // Resistive voltage source in y, periodic trapezoidal pulse
			   // syntax: comp X Y Z ez 51 R phase freq atten amplitude
			   // note: all timings are in picoseconds
			   // default: R=50, phase=0rad, freq=1E06, atten=0, amp=1.	
		fC = (fdel_t*fdel_y)/(2*fe*eo*fdel_x*fdel_z);
		if (plist == NULL)
			fC = fC/50.0;
		else
			fC = fC/(plist->fP);
		teh.stCCey[uEym].fC0 = (1-fC)/(1+fC);
		teh.stCCey[uEym].fC1 = fdel_t/(fdel_z*fe*eo*(1+fC));
		teh.stCCey[uEym].fC2 = fdel_t/(fdel_x*fe*eo*(1+fC));
		teh.stCCey[uEym].fC3 = (2*fC)/(fdel_y*(1+fC));
		
		// Assign default values
		teh.stCCey[uEym].fC4 = 0;	// 0 radian
		teh.stCCey[uEym].fC5 = 1E06; // 1MHz
		teh.stCCey[uEym].fC6 = 0; // 0 Neper
		teh.stCCey[uEym].fC7 = 1; // 1 Volt

		// get other parameters
		if (plist == NULL)
			break;
		else
			plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCey[uEym].fC4 = plist->fP; // get phase
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCey[uEym].fC5 = plist->fP; // get frequency
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCey[uEym].fC6 = plist->fP; // get attenuation
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCey[uEym].fC7 = plist->fP; // get amplitude
		break;
	case (52): // y oriented DC source
			   // syntax: comp X Y Z ez 52 R Vdc
			   // default: R=50, Vdc=1.0
		fC = (fdel_t*fdel_y)/(2*fe*eo*fdel_y*fdel_z);
		if (plist == NULL)
			fC = fC/50.0;
		else
			fC = fC/(plist->fP);
		teh.stCCey[uEym].fC0 = (1-fC)/(1+fC);
		teh.stCCey[uEym].fC1 = fdel_t/(fdel_z*fe*eo*(1+fC));
		teh.stCCey[uEym].fC2 = fdel_t/(fdel_x*fe*eo*(1+fC));
		teh.stCCey[uEym].fC3 = (2*fC)/(fdel_y*(1+fC));
		// Assign default values
		teh.stCCey[uEym].fC4 = 1;	// 1V

		// Get other parameters
		if (plist == NULL)
			break;
		else
		{
			plist = plist->next;
			teh.stCCey[uEym].fC4 = plist->fP; // get dc voltage
		}
		break;
	case (100): // Y oriented resistor
				// syntax: comp X Y Z ey 100 R
				// default: R=50.
		if (plist == NULL)
			fC = 50;
		else
			fC = plist->fP;
		fC = (fdel_t*fdel_y)/(2*fe*eo*fC*fdel_x*fdel_z);
		teh.stCCey[uEym].fC0 = (1-fC)/(1+fC);
		teh.stCCey[uEym].fC1 = fdel_t/(fdel_z*fe*eo*(1+fC));
		teh.stCCey[uEym].fC2 = fdel_t/(fdel_x*fe*eo*(1+fC));
		teh.stCCey[uEym].fC3 = fC*(2*fe*eo/fdel_t);	// store effective conductance
		break;
	case (101): // Y oriented capacitor
				// syntax: comp X Y Z ey 101 C
				// default: C = 1pF
		if (plist == NULL)
			fC = (1.0E-12*fdel_y)/(fe*eo*fdel_x*fdel_z);
		else
			fC = ((plist->fP)*fdel_y)/(fe*eo*fdel_x*fdel_z);
		teh.stCCey[uEym].fC0 = fdel_t/(fe*eo*(1+fC)*fdel_z);
		teh.stCCey[uEym].fC1 = fdel_t/(fe*eo*(1+fC)*fdel_x);
		break;
	case (150): // Y oriented diode.
				// syntax: comp X Y Z ey 150 Is n Vj M Cj0 TT FC R
				// default: Is = 0.1nA, n = 1, Vj=0.6, M=0.5, Cj0=0
				//			TT=0, FC=0.5, R=0.
				// FC is also limited to a maximum of 0.98, to prevent
				// divide by 0 error.
		unsigned char ucor; // temporary variable

		if (plist == NULL) // assign direction/orientation
			ucor = 2;	   // 2 for default	
		else if (plist->fP == 0)
			ucor = 2;	   // +ve y direction
		else
			ucor = 3;	   // -ve y direction
		
		teh.stCCey[uEym].plist = new tEDiode(fdel_x,fdel_y,fdel_z,fdel_t,
									fe,ucor,plist); // create a new object for diode
		teh.stCCey[uEym].umode = 150; // set computation mode
		if (teh.stCCey[uEym].plist == NULL)
		{
			teh.nEStatus = OUT_OF_MEMORY; // error
		}
		break;
	case (180): // Y oriented NPN transistor, assignment of BE junction
				// syntax: comp X Y Z ez 180 O Is Bf Br nF nR Va Ikf Ise ne Vb Ikr Isc nc 
			    // Vje Me Cje TTf Vjc Mc Cjc TTr FC
				// default: O=0, Is=0.1nA, Bf=100, Br=1, nF=1, nR=1, Va=500, Ikf=100, 
			    // Ise=0, ne=1.5, Vb=500, Ikr=100, Isc=0, nc=2, Vje=Vjc=0.6, Me=Mc=0.5, 
				// Cje=Cjc=0, TTf=TTr=0, FC=0.5
		tEBjt * ptrBjt1; // temporary pointer to BJT object

		if (plist == NULL) // assign direction/orientation
			uor = 2;	   // 2 for default	
		else if (plist->fP == 0)
			uor = 2;	   // +ve y direction
		else
			uor = 3;	   // -ve y direction

		teh.stCCey[uEym].plist = new tEBjt;
		teh.stCCey[uEym].fC0 = uor;	// 2 for +ve y, 3 for -ve y.
		teh.stCCey[uEym].umode = 180; 
		ptObjtrack = teh.stCCey[uEym].plist;	// track current object
		uIndtrack = uEym;						// track current index
		ptrBjt1 = (tEBjt *) ptObjtrack;
		ptrBjt1->AssignBEjunction(i,j,k,fdel_x,fdel_y,fdel_z,fdel_t,fe,uor,
								plist); // assign BE junction
		if (teh.stCCey[uEym].plist == NULL)
		{
			teh.nEStatus = OUT_OF_MEMORY; // error
			break;
		}
		if (nSimMode == REUSE) // if simulation mode is REUSE initialize internal BE field
		{
			ptrBjt1->Set_BEfield(teh.fEy[i][j][k]);
		}
		break;
	case (181):	// Y oriented NPN transistor, assignment of BC junction
				// syntax: comp X Y Z ey 181 O
				// note: this statement must always follows the comp X Y Z e 180....
				//  	 statement.
		tEBjt * ptrBjt2; // temporary pointer to BJT object

		if (plist == NULL) // assign direction/orientation
			uor = 2;	   // 2 for default	
		else if (plist->fP == 0)
			uor = 2;	   // +ve y direction
		else
			uor = 3;	   // -ve y direction

		if (ptObjtrack != NULL) // if the object tracking pointer is properly assigned
		{
			ptrBjt2 = (tEBjt *) ptObjtrack;
			ptrBjt2->AssignBCjunction(i,j,k,fdel_x,fdel_y,fdel_z,fdel_t,fe,
									  uor); // assign BC junction
			teh.SetEyMode(180,i,j,k);	// fix computation mode the same as 180
			teh.SetEyParam(uEym,i,j,k);	// fix index into computation coefficient structure
			teh.stCCey[uEym].plist = ptrBjt2;
			teh.stCCey[uEym].fC0 = uor;	// 2 for +ve y, 3 for -ve y.
			teh.stCCey[uEym].umode = 181;  
			uEym++;
			ptObjtrack = NULL; // reset object tracking pointer
			uIndtrack = 0; // reset index tracking variable
			if (nSimMode == REUSE) // if simulation mode is REUSE initialize internal BE field
			{
				ptrBjt2->Set_BCfield(teh.fEy[i][j][k]);
			}
			return TRUE;
		}
		else
			return FALSE;
	case (200): // Y oriented inductor, a memory element.
				// syntax: comp X Y Z ey 200 L
				// default: L = 1nH
		if (plist == NULL)
			fC = 1.0E-9;
		else
			fC = plist->fP;
		teh.stCCey[uEym].fC2 = (fdel_y*fdel_t*fdel_t)/(fC*fe*eo*fdel_x*fdel_z); // C2y
		teh.stCCey[uEym].fC0 = fdel_t/(eo*fe*fdel_z); // C0y
		teh.stCCey[uEym].fC1 = fdel_t/(eo*fe*fdel_x); // C1y	
		teh.stCCey[uEym].fC3 = 0; // Initialize the memory element.
		break;
	default:
		return FALSE;
	}

	teh.SetEyMode(umode,i,j,k);	// fix computation mode
	teh.SetEyParam(uEym,i,j,k);	// fix index into computation coefficient structure
	uEym++; //increment max index
	return TRUE;
}

// Function to create new entry in Ez field computation constants array.
// Each row of the array can store 10 constants: C0z, C1z to C9z.
// These constants depends on the dielectric constant and characteristic
// of the cube.
// arguments:
// teh = reference to class object tEHfdtd
// uEzm = The last valid row index to array
// umode = The update mode for the field
// fe = The associated dielectric constant
// fsg = The associated conductivity 
// i = x coordinate
// j = y coordinate
// k = z coordinate
// plist = pointer to parameter list
// return value: TRUE if success FALSE otherwise 
BOOL tCube::CreateEzConstEntry(tEHfdtd &teh, unsigned char &uEzm, unsigned char umode,
							   float fe, float fs, int i, int j, int k, lpDOUBLE plist)
{	
	double fC; // general purpose variable
	unsigned char uor; // temporary variable, for storing orientation information

	// check for out of boundary 
	if (uEzm > MAXECONST-1)
	{
		nStatus = OUT_OF_BOUND;
		return FALSE;
	}

	switch (umode)
	{
	case (0): // Normal update mode	
		teh.stCCez[uEzm].fC0 = fdel_t/(eo*fe*fdel_x); // C0z
		teh.stCCez[uEzm].fC1 = fdel_t/(eo*fe*fdel_y); // C1z
		break;
	case (1): // PEC
		teh.stCCez[uEzm].fC0 = 0; // C0z
		teh.stCCez[uEzm].fC1 = 0; // C1z	
		break;
	case (2): // Lossy dielectric, no dispersion.
		fC = fs*fdel_t/2;
		teh.stCCez[uEzm].fC0 = (eo*fe-fC)/(eo*fe+fC); // C0z
		teh.stCCez[uEzm].fC1 = fdel_t/((eo*fe+fC)*fdel_x); // C1z
		teh.stCCez[uEzm].fC2 = fdel_t/((eo*fe+fC)*fdel_y); // C2z	
		break;
	case (50): // Resistive voltage source in z, periodic trapezoidal pulse
			   // syntax: comp X Y Z ez 50 R Vlow Vhigh tdly trise thigh tlow 
			   // note: all timings are in picoseconds
			   // default: R=50, Vlow=0, Vhigh=1, tdly=0ps, trise=1ns, thigh=3ns,
			   // tlow=5ns, amp=1.	
		fC = (fdel_t*fdel_z)/(2*fe*eo*fdel_x*fdel_y);
		if (plist == NULL)
			fC = fC/50.0;
		else
			fC = fC/(plist->fP);
		teh.stCCez[uEzm].fC0 = (1-fC)/(1+fC);
		teh.stCCez[uEzm].fC1 = fdel_t/(fdel_x*fe*eo*(1+fC));
		teh.stCCez[uEzm].fC2 = fdel_t/(fdel_y*fe*eo*(1+fC));
		teh.stCCez[uEzm].fC3 = (2*fC)/(fdel_z*(1+fC));
		
		// Assign default values
		teh.stCCez[uEzm].fC4 = 0.0; // 0V
		teh.stCCez[uEzm].fC5 = 1.0; // 1V
		teh.stCCez[uEzm].fC6 = 0;   // tdly=0
		teh.stCCez[uEzm].fC7 = 1000E-12; // trise/tfall=1ns
		teh.stCCez[uEzm].fC8 = 3000E-12; // thigh=3ns
		teh.stCCez[uEzm].fC9 = 5000E-12; // tlow=5ns

		// get other parameters
		if (plist == NULL)
			break;
		else
			plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCez[uEzm].fC4 = plist->fP; // get Vlow
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCez[uEzm].fC5 = plist->fP; // get Vhigh
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCez[uEzm].fC6 = (plist->fP)*(1E-12); // get delay
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCez[uEzm].fC7 = (plist->fP)*(1E-12); // get rise/fall time
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCez[uEzm].fC8 = (plist->fP)*(1E-12); // get high time
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCez[uEzm].fC9 = (plist->fP)*(1E-12); // get low time
		break;
	case (51): // Resistive voltage source in z, periodic trapezoidal pulse
			   // syntax: comp X Y Z ez 51 R phase freq atten amplitude
			   // note: all timings are in picoseconds
			   // default: R=50, phase=0rad, freq=1E06, atten=0, amp=1.	
		fC = (fdel_t*fdel_z)/(2*fe*eo*fdel_x*fdel_y);
		if (plist == NULL)
			fC = fC/50.0;
		else
			fC = fC/(plist->fP);
		teh.stCCez[uEzm].fC0 = (1-fC)/(1+fC);
		teh.stCCez[uEzm].fC1 = fdel_t/(fdel_x*fe*eo*(1+fC));
		teh.stCCez[uEzm].fC2 = fdel_t/(fdel_y*fe*eo*(1+fC));
		teh.stCCez[uEzm].fC3 = (2*fC)/(fdel_z*(1+fC));
		
		// Assign default values
		teh.stCCez[uEzm].fC4 = 0;	// 0 radian
		teh.stCCez[uEzm].fC5 = 1E06; // 1MHz
		teh.stCCez[uEzm].fC6 = 0; // 0 Neper
		teh.stCCez[uEzm].fC7 = 1; // 1 Volt

		// get other parameters
		if (plist == NULL)
			break;
		else
			plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCez[uEzm].fC4 = plist->fP; // get phase
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCez[uEzm].fC5 = plist->fP; // get frequency
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCez[uEzm].fC6 = plist->fP; // get attenuation
		plist = plist->next;
		if (plist == NULL)
			break;
		teh.stCCez[uEzm].fC7 = plist->fP; // get amplitude
		break;
	case (52): // Z oriented DC source
			   // syntax: comp X Y Z ez 52 R Vdc
			   // default: R=50, Vdc=1.0
		fC = (fdel_t*fdel_z)/(2*fe*eo*fdel_x*fdel_y);
		if (plist == NULL)
			fC = fC/50.0;
		else
			fC = fC/(plist->fP);
		teh.stCCez[uEzm].fC0 = (1-fC)/(1+fC);
		teh.stCCez[uEzm].fC1 = fdel_t/(fdel_x*fe*eo*(1+fC));
		teh.stCCez[uEzm].fC2 = fdel_t/(fdel_y*fe*eo*(1+fC));
		teh.stCCez[uEzm].fC3 = (2*fC)/(fdel_z*(1+fC));
		// Assign default values
		teh.stCCez[uEzm].fC4 = 1;	// 1V

		// Get other parameters
		if (plist == NULL)
			break;
		else
		{
			plist = plist->next;
			teh.stCCez[uEzm].fC4 = plist->fP; // get dc voltage
		}
		break;
	case (100): // Z oriented resistor
				// syntax: comp X Y Z ez 100 R
				// default: R=50.
		if (plist == NULL)
			fC = 50;
		else
			fC = plist->fP;
		fC = (fdel_t*fdel_z)/(2*fe*eo*fC*fdel_x*fdel_y);
		teh.stCCez[uEzm].fC0 = (1-fC)/(1+fC);
		teh.stCCez[uEzm].fC1 = fdel_t/(fdel_x*fe*eo*(1+fC));
		teh.stCCez[uEzm].fC2 = fdel_t/(fdel_y*fe*eo*(1+fC));
		teh.stCCez[uEzm].fC3 = fC*(2*fe*eo/fdel_t);  // store effective conductance
		break;
	case (101): // Z oriented capacitor
				// syntax: comp X Y Z ez 101 C
				// default: C = 1pF
		if (plist == NULL)
			fC = (1.0E-12*fdel_z)/(fe*eo*fdel_x*fdel_y);
		else
			fC = ((plist->fP)*fdel_z)/(fe*eo*fdel_x*fdel_y);
		teh.stCCez[uEzm].fC0 = fdel_t/(fe*eo*(1+fC)*fdel_x);
		teh.stCCez[uEzm].fC1 = fdel_t/(fe*eo*(1+fC)*fdel_y);
		break;
	case (150): // Z oriented diode.
				// syntax: comp X Y Z ez 150 O Is n Vj M Cj0 TT FC R
				// default: O=0, Is=0.1nA, n=1, Vj=0.6, M=0.5, Cj0=0
				//			TT=0, FC=0.5, R=0.
				// FC is also limited to a maximum of 0.98, to prevent
				// divide by 0 error.
		if (plist == NULL) // assign direction/orientation
			uor = 4;	   // 4 for default	
		else if (plist->fP == 0)
			uor = 4;	   // +ve z direction
		else
			uor = 5;	   // -ve z direction	

		teh.stCCez[uEzm].plist = new tEDiode(fdel_x,fdel_y,fdel_z,fdel_t,
									fe,uor,plist); // create a new object for diode
		teh.stCCez[uEzm].umode = 150; // set computation mode
		if (teh.stCCez[uEzm].plist == NULL)
		{
			teh.nEStatus = OUT_OF_MEMORY; // error
		}
		break;
	case (151): // Z oriented diode.
				// syntax: comp X Y Z ez 150 Is n Vj M Cj0 TT FC R
				// default: O=0, Is=0.1nA, n=1, Vj=0.6, M=0.5, Cj0=0
				//			TT=0, FC=0.5, R=0.
				// FC is also limited to a maximum of 0.98, to prevent
				// divide by 0 error.
		teh.stCCez[uEzm].plist = new tEDiodec(fdel_x,fdel_y,fdel_z,fdel_t,
									fe,plist); // create a new object for diode
		teh.stCCez[uEzm].umode = 151; // set computation mode
		if (teh.stCCez[uEzm].plist == NULL)
		{
			teh.nEStatus = OUT_OF_MEMORY; // error
		}
		break;
	case (160): // Z oriented PWL two terminal device.
				// syntax: comp X Y Z ez 160 I0 V0 I1 V1 ... In Vn
				// default: I0 = 0, V0 = 0
		teh.stCCez[uEzm].fC0=(fdel_t*fdel_z)/(2*fe*eo*fdel_x*fdel_y); // C0z
		teh.stCCez[uEzm].fC1 = fdel_t/(fe*eo*fdel_x); // C1z
		teh.stCCez[uEzm].fC2 = fdel_t/(fe*eo*fdel_y);	 // C2z
		teh.stCCez[uEzm].fC3 = fdel_t/(fe*eo*fdel_x*fdel_y); // C3z
		teh.stCCez[uEzm].plist = plist; // assign parameter list
		teh.stCCez[uEzm].umode = 160; // set computation mode
		break;
	case (180): // Z oriented NPN transistor, assignment of BE junction
				// syntax: comp X Y Z ez 180 O Is Bf Br nF nR Va Ikf Ise ne Vb Ikr Isc nc 
			    // Vje Me Cje TTf Vjc Mc Cjc TTr FC
				// default: O=0, Is=0.1nA, Bf=100, Br=1, nF=1, nR=1, Va=500, Ikf=100, 
			    // Ise=0, ne=1.5, Vb=500, Ikr=100, Isc=0, nc=2, Vje=Vjc=0.6, Me=Mc=0.5, 
				// Cje=Cjc=0, TTf=TTr=0, FC=0.5
		tEBjt * ptrBjt1; // temporary pointer to BJT object

		if (plist == NULL) // assign direction/orientation
			uor = 4;	   // 4 for default	
		else if (plist->fP == 0)
			uor = 4;	   // +ve z direction
		else
			uor = 5;	   // -ve z direction

		teh.stCCez[uEzm].plist = new tEBjt;
		teh.stCCez[uEzm].fC0 = uor;	// 4 for +ve z orientation, 5 for -ve z orientation.
		teh.stCCez[uEzm].umode = 180; 
		ptObjtrack = teh.stCCez[uEzm].plist;	// track current object
		uIndtrack = uEzm;						// track current index
		ptrBjt1 = (tEBjt *) ptObjtrack;
		ptrBjt1->AssignBEjunction(i,j,k,fdel_x,fdel_y,fdel_z,fdel_t,fe,uor,
								plist); // assign BE junction
		if (teh.stCCez[uEzm].plist == NULL)
		{
			teh.nEStatus = OUT_OF_MEMORY; // error
			break;
		}
		if (nSimMode == REUSE) // if simulation mode is REUSE initialize internal BE field
		{
			ptrBjt1->Set_BEfield(teh.fEz[i][j][k]);
		}
		break;
	case (181):	// Z oriented NPN transistor, assignment of BC junction
				// syntax: comp X Y Z ez 181 O
				// note: this statement must always follows the comp X Y Z e 180....
				//  	 statement.
		tEBjt * ptrBjt2; // temporary pointer to BJT object

		if (plist == NULL) // assign direction/orientation
			uor = 4;	   // 4 for default	
		else if (plist->fP == 0)
			uor = 4;	   // +ve z direction
		else
			uor = 5;	   // -ve z direction

		if (ptObjtrack != NULL) // if the object tracking pointer is properly assigned
		{
			ptrBjt2 = (tEBjt *) ptObjtrack;
			ptrBjt2->AssignBCjunction(i,j,k,fdel_x,fdel_y,fdel_z,fdel_t,fe,
									  uor); // assign BC junction
			teh.SetEzMode(180,i,j,k);	// fix computation mode the same as 180
			teh.SetEzParam(uEzm,i,j,k);	// fix index into computation coefficient structure
			teh.stCCez[uEzm].fC0 = uor;	// 4 for +ve z orientation, 5 for -ve z orientation.
			teh.stCCez[uEzm].plist = ptrBjt2;
			teh.stCCez[uEzm].umode = 181;
			uEzm++;
			ptObjtrack = NULL; // reset object tracking pointer
			uIndtrack = 0; // reset index tracking variable
			if (nSimMode == REUSE) // if simulation mode is REUSE initialize internal BE field
			{
				ptrBjt2->Set_BCfield(teh.fEz[i][j][k]);
			}
			return TRUE;
		}
		else
			return FALSE;
	case (200): // Z oriented inductor, a memory element.
				// syntax: comp X Y Z ez 200 L
				// default: L = 1nH
		if (plist == NULL)
			fC = 1.0E-9;
		else
			fC = plist->fP;
		teh.stCCez[uEzm].fC2 = (fdel_z*fdel_t*fdel_t)/(fC*fe*eo*fdel_x*fdel_y); // C2z
		teh.stCCez[uEzm].fC0 = fdel_t/(eo*fe*fdel_x); // C0z
		teh.stCCez[uEzm].fC1 = fdel_t/(eo*fe*fdel_y); // C1z	
		teh.stCCez[uEzm].fC3 = 0; // Initialize the memory element.
		break;
	default:
		return FALSE;
	}

	teh.SetEzMode(umode,i,j,k);	// fix computation mode
	teh.SetEzParam(uEzm,i,j,k);	// fix index into computation coefficient structure
	uEzm++; //increment max index
	return TRUE;
}

// Function to create new entry in Hx field computation constants array.
// Each row of the array can store 4 constants: C0x, C1x to C3x.
// These constants depends on the relative permeability and characteristic
// of the cube.
// arguments:
// teh = Reference to class object tEHfdtd
// uHxm = The last valid row index to array
// umode = The update mode for the field
// i = x coordinate
// j = y coordinate
// k = z coordinate
// return: TRUE if success FALSE otherwise
BOOL tCube::CreateHxConstEntry(tEHfdtd &teh,unsigned char &uHxm, unsigned char umode,
							   int i, int j, int k)
{
	// check for out of boundary 
	if (uHxm > MAXECONST-1)
	{
		nStatus = OUT_OF_BOUND;
		return FALSE;
	}

	switch (umode)
	{
	case (0): // Normal update mode
		(teh.fChx)[uHxm][0] = fdel_t/(uo*fdel_y); // C0x
		(teh.fChx)[uHxm][1] = fdel_t/(uo*fdel_z); // C1x	
		break;
	case (1): // Field normal to PEC
		(teh.fChx)[uHxm][0] = 0;
		(teh.fChx)[uHxm][1] = 0;
		break;
	case (2): // Partial field normal to PEC, upper diagonal triangle is PEC
		(teh.fChx)[uHxm][0] = (2*fdel_t)/(uo*fdel_z); // C0x
		(teh.fChx)[uHxm][1] = (2*fdel_t)/(uo*fdel_y); // C1x
		break;
	case (3): // Partial field normal to PEC, lower diagonal triangle is PEC
		(teh.fChx)[uHxm][0] = (2*fdel_t)/(uo*fdel_z); // C0x
		(teh.fChx)[uHxm][1] = (2*fdel_t)/(uo*fdel_y); // C1x
		break;
	default:
		return FALSE;
	}
	uHxm++; // increment max index
	return TRUE;
}

// Function to create new entry in Hy field computation constants array.
// Each row of the array can store 4 constants: C0y, C1y to C3y.
// These constants depends on the relative permeability and characteristic
// of the cube.
// arguments:
// teh = Reference to class object tEHfdtd
// uHym = The last valid row index to array
// umode = The update mode for the field
// i = x coordinate
// j = y coordinate
// k = z coordinate
// return: TRUE if success FALSE otherwise
BOOL tCube::CreateHyConstEntry(tEHfdtd &teh,unsigned char &uHym, unsigned char umode,
							   int i, int j, int k)
{
	// check for out of boundary 
	if (uHym > MAXECONST-1)
	{
		nStatus = OUT_OF_BOUND;
		return FALSE;
	}

	switch (umode)
	{
	case (0): // Normal update mode
		(teh.fChy)[uHym][0] = fdel_t/(uo*fdel_z); // C0y
		(teh.fChy)[uHym][1] = fdel_t/(uo*fdel_x); // C1y	
		break;
	case (1): // Field normal to PEC
		(teh.fChy)[uHym][0] = 0;
		(teh.fChy)[uHym][1] = 0;
		break;
	case (2): // Partial field normal to PEC, upper diagonal triangle is PEC
		(teh.fChy)[uHym][0] = (2*fdel_t)/(uo*fdel_z); // C0y
		(teh.fChy)[uHym][1] = (2*fdel_t)/(uo*fdel_x); // C1y
		break;
	case (3): // Partial field normal to PEC, lower diagonal triangle is PEC
		(teh.fChy)[uHym][0] = (2*fdel_t)/(uo*fdel_z); // C0y
		(teh.fChy)[uHym][1] = (2*fdel_t)/(uo*fdel_x); // C1y
		break;
	default:
		return FALSE;
	}
	uHym++; 	// increment max index
	return TRUE;
}

// Function to create new entry in Hz field computation constants array.
// Each row of the array can store 4 constants: C0z, C1z to C3z.
// These constants depends on the relative permeability and characteristic
// of the cube.
// arguments:
// teh = Reference to class object tEHfdtd
// uHzm = The last valid row index to array
// umode = The update mode for the field
// i = x coordinate
// j = y coordinate
// k = z coordinate
// return: TRUE if success FALSE otherwise
BOOL tCube::CreateHzConstEntry(tEHfdtd &teh,unsigned char &uHzm, unsigned char umode,
							   int i, int j, int k)
{
	// check for out of boundary 
	if (uHzm > MAXECONST-1)
	{
		nStatus = OUT_OF_BOUND;
		return FALSE;
	}

	switch(umode)
	{
	case (0): // Normal update mode
		(teh.fChz)[uHzm][0] = fdel_t/(uo*fdel_x); // C0z
		(teh.fChz)[uHzm][1] = fdel_t/(uo*fdel_y); // C1z	
		break;
	case (1): // Field normal to PEC
		(teh.fChz)[uHzm][0] = 0;
		(teh.fChz)[uHzm][1] = 0;
		break;
	case (2): // Partial field normal to PEC, upper diagonal triangle is PEC
		(teh.fChz)[uHzm][0] = (2*fdel_t)/(uo*fdel_x); // C0z
		(teh.fChz)[uHzm][1] = (2*fdel_t)/(uo*fdel_y); // C1z
		break;
	case (3): // Partial field normal to PEC, lower diagonal triangle is PEC
		(teh.fChz)[uHzm][0] = (2*fdel_t)/(uo*fdel_x); // C0z
		(teh.fChz)[uHzm][1] = (2*fdel_t)/(uo*fdel_y); // C1z
		break;
	case (4): // Partial field normal to PEC, cross lower diagonal triangle is PEC
		(teh.fChz)[uHzm][0] = (2*fdel_t)/(uo*fdel_x); // C0z
		(teh.fChz)[uHzm][1] = (2*fdel_t)/(uo*fdel_y); // C1z	
		break;
	case (5): // Partial field normal to PEC, cross upper diagonal triangle is PEC
		(teh.fChz)[uHzm][0] = (2*fdel_t)/(uo*fdel_x); // C0z
		(teh.fChz)[uHzm][1] = (2*fdel_t)/(uo*fdel_y); // C1z	
		break;
	default:
		return FALSE;
	}	
	uHzm++; // increment max index
	return TRUE;
}