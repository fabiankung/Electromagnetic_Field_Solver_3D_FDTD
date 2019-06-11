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

// Filename		 : fdtd_component_class.h
// Remark		 : Definition for lump component class
// Programming Language	 : C++
// Programmer	 : Fabian Kung Wai Lee


#ifndef _COMPONENT_CLASS_
#define _COMPONENT_CLASS_

///////////////////////////////////////////////////////////
// --- Class object for component --- /////////////////////
///////////////////////////////////////////////////////////

// --- Base class for diode ---
class tEDiode
{
	// data prototype
public:
	unsigned char ucxyz; // orientation indicator
						 // 0=positive x, 1=negative x
						 // 2=positive y, 3=negative y
						 // 4=positive z, 5=negative z.
private:
	// Diode models parameters
	double fIs; // reverse saturation current
	double fnVt; // product of emission coefficient and kT/q.
	double fVj; // junction potential
	double fCj0; // space charge capacitance when no bias
	double fM;	// grading coefficient
	double fTT; // sum of transit time for hole and electron
	double fFC; // coefficient for forward-bias depletion capacitance formula
	double fR;  // junction resistance, also to account for high level injection
	// Computation coefficients
	double fC0;
	double fC1;
	double fC2;
	double fC3;
	double fC4;
	double fC5;
	double fC6;
	double fC7;
	float	fn; // orientation constant
	double fdel_d; // length of junction
	double fdi; // Register to store current di/dt
	// general purpose variables
	double ftemp1;
	double ftemp2;

	// function prototype
private:
	void AssignVariable(lpDOUBLE &plist, double &fP, double fdefault);

protected:
	
	double PN_junc_cap(double fV);
	double PN_junc_dcap(double fV);
	double PN_junc_i(double fV);
	double PN_junc_di(double fV);
	double PN_junc_ddi(double fV);

public:
	double updateEfield(double fE, double fH01, double fH00, double fH11,
						 double fH10);
	tEDiode(double fdel_x=0.001, double fdel_y=0.001, double fdel_z=0.001, // constructor
		 double fdel_t = 1.0E-12, double fe=1.0, unsigned char ucor=0, lpDOUBLE plist=NULL);  
	~tEDiode() { }; // destructor
};

// --- Base class for transistor ---
class tEBjt
{
	// data prototype
private:
	// transistor model parameters
	double fIs;	// saturation current
	double fBf; // ideal maximum forward current gain
	double fBr; // ideal maximum reverse current gain
	double fVa;	// forward Early voltage
	double fIkf; // forward beta high-current roll-off 'knee' current
	double fIse; // BE leakage saturation current
	double fVb; // reverse Early voltage
	double fIkr; // corner for reverse beta high-current roll-off
	double fIsc; // BC leakage saturation current
	double fnVtc; // product of emission coefficient and kT/q for collector junction
	double fnVte; // product of emission coefficient and kT/q for emitter junction
	double fnVtcl; // product of emission coefficient and kT/q for BC leakage
	double fnVtel; // product of emission coefficient and kT/q for BE leakage
	double frc; // collector resistance
	double fre; // emitter resistance
	double frb; // base resistance
	double fVje; // junction potential (BE)
	double fCje; // space charge capacitance when no bias (BE)
	double fMe;	// grading coefficient (BE)
	double fTTf; // sum of transit time for hole and electron (BE)
	double fVjc; // junction potential (BC)
	double fCjc; // space charge capacitance when no bias (BC)
	double fMc;	// grading coefficient (BC)
	double fTTr; // sum of transit time for hole and electron (BC)
	double fFC; // coefficient for forward-bias depletion capacitance formula (BE&BC)

	// Computation coefficients for BC junction
	double fC0c;
	double fC1c;
	double fC2c;
	double fC3c;
	double fC4c;
	double fC5c;
	double fC6c;
	float fnc; // BC junction orientation coefficient
	// Computation coefficients for BE junction
	double fC0e;
	double fC1e;
	double fC2e;
	double fC3e;
	double fC4e;
	double fC5e;
	double fC6e;
	float fne; // BE junction orientation coefficient

	// Device state variables
	double fQb; // majority charge carrier in Base
	double fIcc; // 
	double fIec; //

	// Orientation coefficients

	unsigned char uobe; // orientation indicator for BE junction
	unsigned char uobc; // orientation indicator for BC junction
						// 0=positive x, 1=negative x
						// 2=positive y, 3=negative y
						// 4=positive z, otherwise=negative z.
	// Internal fields
	double fEbe;	// BE junction E field
	double fEbc;	// BC junction E field
	double fH00c,fH01c,fH10c,fH11c; // BC junction H fields
	double fH00e,fH01e,fH10e,fH11e;	// BE junction H fields
	double fdel_de;	// BE junction length
	double fdel_dc;	// BC junction length

	double fEben;	// previous BE junction E field
	double fEbcn;	// previous BC junction E field 
	double febe;	// dielectric constant of BE junction
	double febc;    // dielectric constant of BC junction

	double fEbc1;
	double fEbc2;
	double fEbc3;
	double fEbc4;
	// index of internal E fields
	int ie,je,ke;	// BE junction E field index
	int ic,jc,kc;	// BC junction E field index

	// General purpose coefficients
	unsigned char uState;	// state variable
							// 0 = no field update
							// 1 = update either BE or BC field
							// 2 = update both BE and BC field
	// function prototype
private:
	void AssignVariable(lpDOUBLE &plist, double &fP, double fdefault);

protected:
	void Update_BJTEfield(void);
	double Ce(double fVe, double fVc);
	double Cc(double fVe,double fVc);
	double Qb(double fVe, double fVc);
	double dQb_Vbe(double fVe, double fVc);
	double dQb_Vbc(double fVe, double fVc);
	double Icc(double fVe, double fVc);
	double dIcc_Vbe(double fVe, double fVc);
	double dIcc_Vbc(double fVe, double fVc);
	double Iec(double fVe, double fVc);
	double dIec_Vbe(double fVe, double fVc);
	double dIec_Vbc(double fVe, double fVc);
	double Ic1(double fVe, double fVc);
	double dIc1_Vbe(double fVe, double fVc);
	double dIc1_Vbc(double fVe, double fVc);
	double Ic2(double fVe, double fVc);
	double dIc2_Vbe(double fVe, double fVc);
	double dIc2_Vbc(double fVe, double fVc);
	double ILe(double fVe);
	double ILc(double fVc);
	double dILe_Vbe(double fVe);
	double dILc_Vbc(double fVc);

public:
	void AssignBEjunction(int i=1, int j=1, int k=1, double fdel_x=0.001, 
			double fdel_y=0.001, double fdel_z=0.001, double fdel_t=1.0E-12, 
			double fe=1, unsigned char ucor=0, lpDOUBLE plist=NULL);
	void AssignBCjunction(int i=1, int j=1, int k=1, double fdel_x=0.001,
			double fdel_y=0.001, double fdel_z=0.001, double fdel_t=1.0E-12,
			double fe=1, unsigned char ucor=0);
	void Update_Efield(tEHfdtd *pteh, int i, int j, int k, unsigned char uor, double fH01,
					double fH00, double fH11, double fH10);
	void Set_BEfield(double fE);
	void Set_BCfield(double fE);
	tEBjt(double fEe=0, double fEc=0);  // constructor
	~tEBjt() {};  // destructor
};


// --- Base class for diode (Ciampolini and Piket-May's approach) ---
class tEDiodec
{
	// data prototype
public:

private:
	// Diode models parameters
	double fIs; // reverse saturation current
	double fnVt; // product of emission coefficient and kT/q.
	double fVj; // junction potential
	double fCj0; // space charge capacitance when no bias
	double fM;	// grading coefficient
	double fTT; // sum of transit time for hole and electron
	double fFC; // coefficient for forward-bias depletion capacitance formula
	double fR;  // junction resistance, also to account for high level injection
	// Computation coefficients
	double fC0;
	double fC1;
	double fC2;
	
	// function prototype
private:
	void AssignVariable(lpDOUBLE &plist, double &fP, double fdefault);

protected:
	double Idiode(double fEn, double fE);
	double dIdiode(double fEn, double fE);

public:
	double updateEfield(double fE, double fH01, double fH00, double fH11,
						 double fH10);
	tEDiodec(double fdel_x=0.001, double fdel_y=0.001, double fdel_z=0.001, // constructor
		 double fdel_t = 1.0E-12, double fe=1.0, lpDOUBLE plist=NULL);  
	~tEDiodec() { }; // destructor
};

#endif