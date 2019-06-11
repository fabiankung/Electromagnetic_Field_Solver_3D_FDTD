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


// Filename		 : fdtd_header.h
// Remark		 : Constants definition for FDTD datatypes
// Programming Language  : C++
// Author		 : Fabian Kung Wai Lee

// --- Definition of macro ---
#define MAXX 100
#define MAXY 100
#define MAXZ 36

#define MAXECONST 84 
#define MAXHCONST 12

#define DATAFILENAME data.out	// data file name

#define INSTABILITY_LIMIT 1E+64


// Enumerated Constants for variable nProgramStatus, nStatus and nSimMode
enum {NORMAL,OUT_OF_MEMORY,OUT_OF_BOUND,FILE_ERROR,SIM_ERROR};
//#define NORMAL			0
#define SIMULATION		1
#define REUSE			2



