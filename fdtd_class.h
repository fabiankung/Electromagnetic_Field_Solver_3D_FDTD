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


// Filename		 : fdtd_class.h
// Remark		 : Declaration for structures and class for E and H fields objects
// Programming Language	 : C++
// Author		 : Fabian Kung Wai Lee


#ifndef _FDTD_CLASS_
#define _FDTD_CLASS_

///////////////////////////////////////////////////////////
// --- Cube object definition --- /////////////////////////
///////////////////////////////////////////////////////////
// This structure is used to store the parameter list of a field
// component
typedef struct stDOUBLE
{
	double fP;		// a double variable 
	struct stDOUBLE *next;	// pointer to stDOUBLE data structure
} stDouble;

typedef stDOUBLE *lpDOUBLE;

// Each cube is equivalent to a single Yee's cell
struct stCUBE 
{
	float e;	// Permittivity of the cube
	float sg;	// Conductivity of the cube 
	short d;	// Characteristics of the cube dielectric
	short n;	// Characteristics of the cube surface
	lpDOUBLE plist; // pointer to optional parameter list
};

///////////////////////////////////////////////////////////
// --- Component object definition --- ////////////////////
///////////////////////////////////////////////////////////

typedef struct stCOMP
{
	char cType[4]; // determine the field type to be implemented as a component
	unsigned char  uType;	   // determine the component type
	int	 nX;	   // X position
	int  nY;       // Y position 
	int  nZ;       // Z position
	lpDOUBLE plist;	// parameter list
	int  nR;	   // reference integer
	stCOMP *next;  // pointer to next item
} stComp;

typedef stCOMP *lpCOMP;

///////////////////////////////////////////////////////////
// --- Probe object definition --- ////////////////////////
///////////////////////////////////////////////////////////
typedef struct stPROBE		
{
	int	nX;	// x coordinate 
	int	nY;	// y coordinate
	int	nZ;	// z coordinate
	int nType;	// type: 0=Hx, 1=Hy, 2=Hz, 3=Ex, 4=Ey, 5=Ez.
				// 6=Line integral (V)
	struct stPROBE *next; // pointer to its own datatype
} stProbe;

typedef stPROBE *lpPROBE; 

///////////////////////////////////////////////////////////
// --- Parameter object definition --- ////////////////////
///////////////////////////////////////////////////////////
typedef struct stPARAM
{
	int nR; // reference integer
	lpDOUBLE plist; // pointer to parameter list
	struct stPARAM *next; // pointer to its own datatype
} stParam;

typedef stPARAM *lpPARAM;

///////////////////////////////////////////////////////////
// --- Computation coefficient object definition --- //////
///////////////////////////////////////////////////////////
struct stCCTAB 
{
	double fC0;	
	double fC1;
	double fC2;
	double fC3;
	double fC4;
	double fC5;
	double fC6;
	double fC7;
	double fC8;
	double fC9;
	unsigned char umode;	
	void	*plist;
};

///////////////////////////////////////////////////////////
// --- Class object for electric field component --- //////
///////////////////////////////////////////////////////////

// --- Base class for electric field ---
class tEfield 
{
protected:
// Data prototype
	unsigned char (*umodex)[MAXY+1][MAXZ+1];	// pointer to array for storing
	unsigned char (*umodey)[MAXY][MAXZ+1];		// computational mode and parameter
	unsigned char (*umodez)[MAXY+1][MAXZ];		// corresponding the field component.
	unsigned char (*uparamex)[MAXY+1][MAXZ+1];
	unsigned char (*uparamey)[MAXY][MAXZ+1];
	unsigned char (*uparamez)[MAXY+1][MAXZ];

// function prototype

public:
// Data prototype:
	int		nEStatus;
	double	(*fEx)[MAXY+1][MAXZ+1];		// pointers to field component arrays
	double	(*fEy)[MAXY][MAXZ+1];
	double	(*fEz)[MAXY+1][MAXZ];
	stCCTAB	stCCex[MAXECONST];	// structure for storing computational coefficients				
	stCCTAB stCCey[MAXECONST];
	stCCTAB stCCez[MAXECONST];

// function prototype:
	void InitEField(int SimMode);
	unsigned char SetExMode(unsigned char umode, int nx,
							   int ny, int nz);
	unsigned char SetEyMode(unsigned char umode, int nx,
							   int ny, int nz);
	unsigned char SetEzMode(unsigned char umode, int nx,
							   int ny, int nz);
	unsigned char SetExParam(unsigned char uindex, int nx,
							   int ny, int nz);
	unsigned char SetEyParam(unsigned char uindex, int nx,
							   int ny, int nz);
	unsigned char SetEzParam(unsigned char uindex, int nx,
							   int ny, int nz);
	unsigned char GetExMode(int nx, int ny, int nz);
	unsigned char GetEyMode(int nx, int ny, int nz);
	unsigned char GetEzMode(int nx, int ny, int nz);
	tEfield();	// constructor
	~tEfield();	// destructor
};


///////////////////////////////////////////////////////////
// --- Class object for magnetic field component --- //////
///////////////////////////////////////////////////////////

// --- Base class for magnetic field ---
class tHfield
{
protected:
// data prototype
	unsigned char (*umodhx)[MAXY][MAXZ];
	unsigned char (*umodhy)[MAXY+1][MAXZ];
	unsigned char (*umodhz)[MAXY][MAXZ+1];
	unsigned char (*uparamhx)[MAXY][MAXZ];
	unsigned char (*uparamhy)[MAXY+1][MAXZ];
	unsigned char (*uparamhz)[MAXY][MAXZ+1];

// function prototype

public:
// data prototype
	int		nHStatus;
	double	(*fHx)[MAXY][MAXZ];		// pointers to field component arrays
	double	(*fHy)[MAXY+1][MAXZ];
	double	(*fHz)[MAXY][MAXZ+1];
	double	(*fChx)[4];				// pointer to array for storing computation constants	
	double  (*fChy)[4];
	double	(*fChz)[4];

// function prototype
	tHfield(); // constructor
	~tHfield(); // destructor
	void InitHField(int SimMode);
	unsigned char SetHxMode(unsigned char umode, int nx,
							   int ny, int nz);
	unsigned char SetHyMode(unsigned char umode, int nx,
							   int ny, int nz);
	unsigned char SetHzMode(unsigned char umode, int nx,
							   int ny, int nz);
	unsigned char SetHxParam(unsigned char uindex, int nx,
							   int ny, int nz);
	unsigned char SetHyParam(unsigned char uindex, int nx,
							   int ny, int nz);
	unsigned char SetHzParam(unsigned char uindex, int nx,
							   int ny, int nz);
	unsigned char GetHxMode(int nx, int ny, int nz);
	unsigned char GetHyMode(int nx, int ny, int nz);
	unsigned char GetHzMode(int nx, int ny, int nz);
};


///////////////////////////////////////////////////////////
// --- Derived class for electric and magnetic field --- //
///////////////////////////////////////////////////////////

class tEHfdtd: public tEfield, public tHfield
{
private:
// data prototype

// function prototype


public:
//data prototype

// function prototype
tEHfdtd();	// constructor
~tEHfdtd() {}; // destructor

// These function updates the electric field components
void UpdateEx_type0(int i, int j, int k);
void UpdateEy_type0(int i, int j, int k);
void UpdateEz_type0(int i, int j, int k);

void UpdateEx_type1(int i, int j, int k);
void UpdateEy_type1(int i, int j, int k);
void UpdateEz_type1(int i, int j, int k);

void UpdateEx_type2(int i, int j, int k);
void UpdateEy_type2(int i, int j, int k);
void UpdateEz_type2(int i, int j, int k);

void UpdateEx_type50(int i, int j, int k, double t);
void UpdateEy_type50(int i, int j, int k, double t);
void UpdateEz_type50(int i, int j, int k, double t);

void UpdateEx_type51(int i, int j, int k, double t, int nSimMode);
void UpdateEy_type51(int i, int j, int k, double t, int nSimMode);
void UpdateEz_type51(int i, int j, int k, double t, int nSimMode);

void UpdateEx_type52(int i, int j, int k, double t);
void UpdateEy_type52(int i, int j, int k, double t);
void UpdateEz_type52(int i, int j, int k, double t);

void UpdateEx_type100(int i, int j, int k);
void UpdateEy_type100(int i, int j, int k);
void UpdateEz_type100(int i, int j, int k);

void UpdateEx_type101(int i, int j, int k);
void UpdateEy_type101(int i, int j, int k);
void UpdateEz_type101(int i, int j, int k);

void UpdateEx_type150(int i, int j, int k);
void UpdateEy_type150(int i, int j, int k);
void UpdateEz_type150(int i, int j, int k);

void UpdateEz_type151(int i, int j, int k);

void UpdateEz_type160(int i, int j, int k);

void UpdateEx_type180(int i, int j, int k);
void UpdateEy_type180(int i, int j, int k);
void UpdateEz_type180(int i, int j, int k);

void UpdateEx_type200(int i, int j, int k);
void UpdateEy_type200(int i, int j, int k);
void UpdateEz_type200(int i, int j, int k);

// These function updates the magnetic field components
void UpdateHx_type0(int i, int j, int k);
void UpdateHy_type0(int i, int j, int k);
void UpdateHz_type0(int i, int j, int k);

void UpdateHx_type1(int i, int j, int k);
void UpdateHy_type1(int i, int j, int k);
void UpdateHz_type1(int i, int j, int k);

void UpdateHx_type2(int i, int j, int k);
void UpdateHy_type2(int i, int j, int k);
void UpdateHz_type2(int i, int j, int k);

void UpdateHx_type3(int i, int j, int k);
void UpdateHy_type3(int i, int j, int k);
void UpdateHz_type3(int i, int j, int k);

void UpdateHz_type4(int i, int j, int k);

void UpdateHz_type5(int i, int j, int k);
};


///////////////////////////////////////////////////////////
// --- Class object for cube --- //////////////////////////
///////////////////////////////////////////////////////////
class tCube	
{
private:
// data prototype:
//	stDOUBLE *(*Px)[MAXY+1][MAXZ+1]; // a pointer to an array of stDOUBLE structure !
//	stDOUBLE *(*Py)[MAXY][MAXZ+1];	// The stDOUBLE structures form a linked list to 
//	stDOUBLE *(*Pz)[MAXY+1][MAXZ];	// store the parameters of a field component
	stCUBE (*stCube)[MAXY][MAXZ]; //pointer to array of stCube structures
	unsigned char uEimaxx;		// index to last valid entry in the computational
	unsigned char uEimaxy;		// constants look-up table for E field.
	unsigned char uEimaxz;
	unsigned char uHimaxx;		// index to last valid entry in the computatiobal
	unsigned char uHimaxy;		// constants look-up table for H field.
	unsigned char uHimaxz;
	void *ptObjtrack;			// a void pointer, for tracking object
	unsigned char uIndtrack;	// track the index to computation coefficient table

// function prototype:
	void AssignEpsilonSigma(unsigned char, int, int, int, float &, float &);
	BOOL CreateExConstEntry(tEHfdtd &, unsigned  char &,unsigned char,
				            float, float, int, int, int, lpDOUBLE);
	BOOL CreateEyConstEntry(tEHfdtd &, unsigned  char &,unsigned char,
				            float, float, int, int, int, lpDOUBLE);
	BOOL CreateEzConstEntry(tEHfdtd &, unsigned  char &,unsigned char,
				            float, float, int, int, int, lpDOUBLE);
	BOOL CreateHxConstEntry(tEHfdtd &, unsigned char &, unsigned char,
							int, int, int);
	BOOL CreateHyConstEntry(tEHfdtd &, unsigned char &, unsigned char,
							int, int, int);
	BOOL CreateHzConstEntry(tEHfdtd &, unsigned char &, unsigned char,
							int, int, int);
	BOOL CompareParam(lpDOUBLE, lpDOUBLE);

public:
// data prototype:
	int nStatus;
	int nNx;	// Max number of cubes along x axis
	int nNy;	// Max number of cubes along y axis
	int nNz;	// Max number of cubes along z axis

// function prototype:
	tCube(int, int, int);	// constructor
	~tCube();	// destructor
	BOOL UpdateCube(int, int ,int, float, float,int, int, lpDOUBLE);
	void AssignComMode(tEHfdtd &);
	void AssignComParam(tEHfdtd &);
	void AssignCompParam(tEHfdtd &);
};

#endif