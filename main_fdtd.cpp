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

// Filename		 : main_fdtd.cpp
// Remark		 : This file contains example of subroutine to use the 
//                         objects in "fdtd_sourceX.cpp" files, where x=1,2,3...8
//                         The subroutine below illustrate the calling of all FDTD functions in 
//                         correct sequence.
// Programming Language  : C++
// Author	 : Fabian Kung Wai Lee


#include "fdtd_header.h"
#include "fdtd_class.h"
#include "fdtd_class_abc.h"
#include "fdtd_component_class.h"


// Global variable declaration for FDTD
// Spatial and time discretization
double fdel_x = 0.00075;		// delta x
double fdel_y = 0.0008;		// delta y
double fdel_z = 0.00055;		// delta z
double fdel_t = 1.0E-12;		// delta t


// Local functions prototype
int FDTD_Run(int, int);

// Global datatype declaration
// The object for storing E and H field components and update methods
tEHfdtd tEH;
tABC	tabc(MAXX-1,MAXY-1,MAXZ-1); 


///////////////////////////////////////////////////////////
// --- User Routine Definition --- ////////////////////////
///////////////////////////////////////////////////////////



// Function to perform FDTD simulation
// It is assumed that the global objects 'tEH', the object for storing 
// E and H field components and their respective methods, and 'tABC', the object for
// storing boundary E and H fields and the respective methods, are properly initialized.
// This means the initial values of all E and H field components must be set (typically
// 0), and each field components update mode (e.g. 'umode') must be set properly.
// You need to write your own routines to perform this.  For instance you can read a 
// text file that describe the model and update all the properties in tEH and tABC
// properly before calling this subroutine.
// argument:
// nTstop - the no. of time step to run the FDTD update routines from current time step
// nCurstep - current time step
// return - the actual no. of time update routines are performed
int FDTD_Run(int nTstop, int nCurstep)
{
	int ntindex;
	int i, j, k;	

	//Check initialization of object
	if (tabc.nStatus != NORMAL)
		return 0;

	// assign limit of model space to ABC
	tabc.nNx = nMaxx;
	tabc.nNy = nMaxy;
	tabc.nNz = nMaxz;

	for (ntindex=0; ntindex <= nTstop; ntindex++)
	{
		// Shift and update boundary E fields 
		tabc.ShiftEfield(tEH);

		// Update interior E & H fields
		// Ex field
		for (i=1; i<=nMaxx; i++)
		{
			for (j=2; j<=nMaxy; j++)
			{
				for (k=2; k<=nMaxz; k++)
				{
					switch(tEH.GetExMode(i,j,k))
					{
						case 0:
							tEH.UpdateEx_type0(i,j,k);
							break;
						case 1:
							tEH.UpdateEx_type1(i,j,k);
							break;
						case 2:
							tEH.UpdateEx_type2(i,j,k);
							break;
						case 50:
							tEH.UpdateEx_type50(i,j,k,(nCurstep + ntindex)*fdel_t);
							break;
						case 51:
							tEH.UpdateEx_type51(i,j,k,(nCurstep + ntindex)*fdel_t,nSimMode);
							break;
						case 52:
							tEH.UpdateEx_type52(i,j,k,(nCurstep + ntindex)*fdel_t);
							break;
						case 100:
							tEH.UpdateEx_type100(i,j,k);
							break;
						case 101:
							tEH.UpdateEx_type101(i,j,k);
							break;
						case 150:
							tEH.UpdateEx_type150(i,j,k);
							break;
						case 180:
							tEH.UpdateEx_type180(i,j,k);
							break;
						case 200:
							tEH.UpdateEx_type200(i,j,k);
						default:
							break;
					}
				}
			}
		}

		// Ey field
		for (i=2; i<=nMaxx; i++)
		{
			for (j=1; j<=nMaxy; j++)
			{
				for (k=2; k<=nMaxz; k++)
				{
					switch(tEH.GetEyMode(i,j,k))
					{
						case 0:
							tEH.UpdateEy_type0(i,j,k);
							break;
						case 1:
							tEH.UpdateEy_type1(i,j,k);
							break;
						case 2:
							tEH.UpdateEy_type2(i,j,k);
							break;
						case 50:
							tEH.UpdateEy_type50(i,j,k,(nCurstep + ntindex)*fdel_t);
							break;
						case 51:
							tEH.UpdateEy_type51(i,j,k,(nCurstep + ntindex)*fdel_t,nSimMode);
							break;
						case 52:
							tEH.UpdateEy_type52(i,j,k,(nCurstep + ntindex)*fdel_t);
							break;
						case 100:
							tEH.UpdateEy_type100(i,j,k);
							break;
						case 101:
							tEH.UpdateEy_type101(i,j,k);
							break;
						case 150:
							tEH.UpdateEy_type150(i,j,k);
							break;
						case 180:
							tEH.UpdateEy_type180(i,j,k);
							break;
						case 200:
							tEH.UpdateEy_type200(i,j,k);
						default:
							break;
					}
				}
			}
		}

		// Ez field
		for (i=2; i<=nMaxx; i++)
		{
			for (j=2; j<=nMaxy; j++)
			{
				for (k=1; k<=nMaxz; k++)
				{
					switch(tEH.GetEzMode(i,j,k))
					{
						case 0:
							tEH.UpdateEz_type0(i,j,k);
							break;
						case 1:
							tEH.UpdateEz_type1(i,j,k);
							break;
						case 2:
							tEH.UpdateEz_type2(i,j,k);
							break;
						case 50:
							tEH.UpdateEz_type50(i,j,k,(nCurstep + ntindex)*fdel_t);
							break;
						case 51:
							tEH.UpdateEz_type51(i,j,k,(nCurstep + ntindex)*fdel_t,nSimMode);
							break;
						case 52:
							tEH.UpdateEz_type52(i,j,k,(nCurstep + ntindex)*fdel_t);
							break;
						case 100:
							tEH.UpdateEz_type100(i,j,k);
							break;
						case 101:
							tEH.UpdateEz_type101(i,j,k);
							break;
						case 150:
							tEH.UpdateEz_type150(i,j,k);
							break;
						case 151:
							tEH.UpdateEz_type151(i,j,k);
							break;
						case 160:
							tEH.UpdateEz_type160(i,j,k);
							break;
						case 180:
							tEH.UpdateEz_type180(i,j,k);
							break;
						case 200:
							tEH.UpdateEz_type200(i,j,k);
							break;
						default:
							break;
					}
				}
			}
		}

		// Enforce absorbing boundary condition
		tabc.EnforceBoundary(tEH); 

		// Hx field
		for (i=2; i<=nMaxx; i++)
		//for (i=1; i<=nMaxx+1; i++)
		{
			for (j=1; j<=nMaxy; j++)
			{
				for (k=1; k<=nMaxz; k++)
				{
					switch(tEH.GetHxMode(i,j,k))
					{
					case 0:				
						tEH.UpdateHx_type0(i,j,k);
						break;
					case 1:
						tEH.UpdateHx_type1(i,j,k);
						break;
					case 2:				
						tEH.UpdateHx_type2(i,j,k);
						break;
					case 3:
						tEH.UpdateHx_type3(i,j,k);
						break;
					default:
						break;
					}
				}
			}
		}

		// Hy field
		for (i=1; i<=nMaxx; i++)
		{
			for (j=2; j<=nMaxy; j++)
			//for (j=1; j<=nMaxy+1; j++)
			{
				for (k=1; k<=nMaxz; k++)
				{
					switch(tEH.GetHyMode(i,j,k))
					{
					case 0:				
						tEH.UpdateHy_type0(i,j,k);
						break;
					case 1:
						tEH.UpdateHy_type1(i,j,k);
						break;
					case 2:				
						tEH.UpdateHy_type2(i,j,k);
						break;
					case 3:
						tEH.UpdateHy_type3(i,j,k);
						break;
					default:
						break;
					}
				}
			}
		}

		// Hz field
		for (i=1; i<=nMaxx; i++)
		{
			for (j=1; j<=nMaxy; j++)
			{
				for (k=2; k<=nMaxz; k++)
				//for (k=1; k<=nMaxz+1; k++)
				{
					switch(tEH.GetHzMode(i,j,k))
					{
					case 0:				
						tEH.UpdateHz_type0(i,j,k);
						break;
					case 1:
						tEH.UpdateHz_type1(i,j,k);
						break;
					case 2:				
						tEH.UpdateHz_type2(i,j,k);
						break;
					case 3:
						tEH.UpdateHz_type3(i,j,k);
						break;
					case 4:
						tEH.UpdateHz_type4(i,j,k);
						break;
					case 5:
						tEH.UpdateHz_type5(i,j,k);
						break;
					default:
						break;
					}
				}
			}
		}

		// retrieve the data in the probe list and write to output file
		// The write interval is determined by global variable nOutint.
		if ((nCurstep + ntindex) > nOutint_off)
		{
			if (((nCurstep + ntindex)%nOutint) == 0)
			{	
				nOutCount++;
				WriteOutData(hdatafile, lpprobe1, nCurstep + ntindex, tEH);
			}
		}
	
		// check for instability - terminate if instability occur
		//if (tabc.CheckABCStability(tEH) == TRUE)
		//	return ntindex;
	}

	return ntindex;
}

