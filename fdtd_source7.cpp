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

// Filename		 : fdtd_source7.cpp 
// Remark		 : Code for simulating a bipolar junction transistor
// Programming Language  : C++
// Author	 : Fabian Kung Wai Lee

#include <math.h> // standard ANSI C/C++ math library header
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

///////////////////////////////////////////////////////////
// --- Functions of Class tEBjt --- ///////////////////////
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
// --- Default Constructor --- ////////////////////////////
///////////////////////////////////////////////////////////
// The constructor performs the following basic functions:
// (1) Initialize the object state.
// (2) Initialize internal E fields to 0.
// arguments:
// fEe = initial BE junction E field
// fEc = initial BC junction E field
tEBjt::tEBjt(double fEe, double fEc)
{
	fEbe = fEe;
	fEbc = fEc;
	
	// initialize state variable
	uState = 0;

	fEbc1 = 0;
	fEbc2 = 0;
	fEbc3 = 0;
	fEbc4 = 0;
}

///////////////////////////////////////////////////////////
// --- Public functions --- ///////////////////////////////
///////////////////////////////////////////////////////////

// function to update the information for BE junction
// arguments:
// i = x coordinate of corresponding E field, default 1.
// j = y coordinate of corresponding E field, default 1.
// k = z coordinate of corresponding E field, default 1.
// fdel_x = x discretization, default 1.0mm.
// fdel_y = y discretization, default 1.0mm.
// fdel_z = z discretization, default 1.0mm.
// fdel_t = time discretization, default 1.0E-12.
// fe = dielectric constant of the E field, default 1.
// ucor = orientation specification. 0=positive x, 1=negative x
//									 2=positive y, 3=negative y
//									 4=positive z, otherwise=negative z.
// plist = point to parameter linked-list, default = NULL.
// return: none.
void tEBjt::AssignBEjunction(int i, int j, int k, double fdel_x, double fdel_y,
							 double fdel_z, double fdel_t, double fe, unsigned char ucor,
							 lpDOUBLE plist)
{
	if (plist != NULL)
		plist = plist->next;	// skip the 1st parameter, the orientation indicator.
	AssignVariable(plist,fIs,1.0e-10); // read Is into fIs
	AssignVariable(plist,fBf,100);	// read Bf
	AssignVariable(plist,fBr,1);	// read Br
	AssignVariable(plist,fnVte,1);  // read ne
	fnVte = fnVte/38.695652; // fnVte
	AssignVariable(plist,fnVtc,1);	// read nc
	fnVtc = fnVtc/38.695652; // fnVtc
	AssignVariable(plist,fVa,5000);	// read Va
	AssignVariable(plist,fIkf,10000);  // read Ikf
	AssignVariable(plist,fIse,0);   // read Ise
	AssignVariable(plist,fnVtel,1.5); // read ne
	fnVtel = fnVtel/38.695652; // fnVtel
	AssignVariable(plist,fVb,5000);  // read Vb
	AssignVariable(plist,fIkr,10000);	// read Ikr
	AssignVariable(plist,fIsc,0);	// read Isc
	AssignVariable(plist,fnVtcl,2.0); // read nc
	fnVtcl = fnVtcl/38.695652; // fnVtcl	
	AssignVariable(plist,fVje,0.6); // read Vje
	AssignVariable(plist,fMe,0.5);	// read Me
	AssignVariable(plist,fCje,0);	// read Cje
	AssignVariable(plist,fTTf,0);	// read TTf
	AssignVariable(plist,fVjc,0.6); // read Vjc
	AssignVariable(plist,fMc,0.5);	// read Mc
	AssignVariable(plist,fCjc,0);	// read Cjc
	AssignVariable(plist,fTTr,0);	// read TTr
	AssignVariable(plist,fFC,0.5);	// read FC
	uobe = ucor;	// assign orientation indicator
	if ((ucor==0) || (ucor==2) || (ucor==4))
		fne = 1;	// positive orientation
	else
		fne = -1;	// negative orientation

	if ((ucor==0) || (ucor==1)) // x oriented diode
	{
		fC0e = fdel_t/(fe*eo*fdel_y); // C0x
		fC1e = fdel_t/(fe*eo*fdel_z); // C1x
		fC2e = (fdel_t*fne)/(fe*eo*fdel_z*fdel_y); // C2x
		fC3e = fC2e*fdel_x/fne; // C3x
		fC4e = (fC3e*fdel_x*fne)/2; // C4x
		fdel_de = fdel_x;
	}
	else if ((ucor==2) || (ucor==3)) // y oriented diode
	{
		fC0e = fdel_t/(fe*eo*fdel_z); // C0y
		fC1e = fdel_t/(fe*eo*fdel_x); // C1y
		fC2e = (fdel_t*fne)/(fe*eo*fdel_x*fdel_z); // C2y
		fC3e = fC2e*fdel_y/fne; // C3y
		fC4e = (fC3e*fdel_y*fne)/2; // C4y
		fdel_de = fdel_y;
	}
	else // z oriented diode
	{
		fC0e = fdel_t/(fe*eo*fdel_x); // C0z
		fC1e = fdel_t/(fe*eo*fdel_y); // C1z
		fC2e = (fdel_t*fne)/(fe*eo*fdel_x*fdel_y); // C2z
		fC3e = fC2e*fdel_z/fne; // C3z
		fC4e = (fC3e*fdel_z*fne)/2; // C4z
		fdel_de = fdel_z;
	}

	fC5e = 1/fnVte;		// C5
	fC6e = 1/fdel_t;	// C6

	// assign indexes
	ie = i;
	je = j;
	ke = k;

	// assign dielectric constant
	febe = fe*eo;
}

// function to update the information for BC junction
// arguments:
// i = x coordinate of corresponding E field, default 1.
// j = y coordinate of corresponding E field, default 1.
// k = z coordinate of corresponding E field, default 1.
// fdel_x = x discretization, default 1.0mm.
// fdel_y = y discretization, default 1.0mm.
// fdel_z = z discretization, default 1.0mm.
// fdel_t = time discretization, default 1.0E-12.
// fe = dielectric constant of the E field, default 1.
// ucor = orientation specification. 0=positive x, 1=negative x
//									 2=positive y, 3=negative y
//									 4=positive z, otherwise=negative z.
// return: none.
void tEBjt::AssignBCjunction(int i, int j, int k, double fdel_x, double fdel_y,
							 double fdel_z, double fdel_t, double fe, unsigned char ucor)
{
	uobc = ucor;	// assign orientation indicator
	if ((ucor==0) || (ucor==2) || (ucor==4))
		fnc = 1;	// positive orientation
	else
		fnc = -1;	// negative orientation

	if ((ucor==0) || (ucor==1)) // x oriented diode
	{
		fC0c = fdel_t/(fe*eo*fdel_y); // C0x
		fC1c = fdel_t/(fe*eo*fdel_z); // C1x
		fC2c = (fdel_t*fnc)/(fe*eo*fdel_z*fdel_y); // C2x
		fC3c = fC2c*fdel_x/fnc; // C3x
		fC4c = (fC3c*fdel_x*fnc)/2; // C4x
		fdel_dc = fdel_x;
	}
	else if ((ucor==2) || (ucor==3)) // y oriented diode
	{
		fC0c = fdel_t/(fe*eo*fdel_z); // C0y
		fC1c = fdel_t/(fe*eo*fdel_x); // C1y
		fC2c = (fdel_t*fnc)/(fe*eo*fdel_x*fdel_z); // C2y
		fC3c = fC2c*fdel_y/fnc; // C3y
		fC4c = (fC3c*fdel_y*fnc)/2; // C4y
		fdel_dc = fdel_y;
	}
	else // z oriented diode
	{
		fC0c = fdel_t/(fe*eo*fdel_x); // C0z
		fC1c = fdel_t/(fe*eo*fdel_y); // C1z
		fC2c = (fdel_t*fnc)/(fe*eo*fdel_x*fdel_y); // C2z
		fC3c = fC2c*fdel_z/fnc; // C3z
		fC4c = (fC3c*fdel_z*fnc)/2; // C4z
		fdel_dc = fdel_z;
	}

	fC5c = 1/fnVtc;		// C5
	fC6c = 1/fdel_t;	// C6

	// assign indexes
	ic = i;
	jc = j;
	kc = k;

	// assign dielectric constant
	febc = fe*eo;
}

// function to update the BC junction E field
void tEBjt::Update_BJTEfield(void)
{
	double fA1, fB1, fC1, fD1, fE1;	// BC junction coefficients
	double fA2, fB2, fC2, fD2, fE2; // BE junction coefficients
	double Vbe, Vbc;
	double ftemp1, ftemp2;
	static double fden, fnome, fnomc;
	static int nCnt = 0;

	nCnt++;
	Vbe = fne*fdel_de*fEbe;
	Vbc = fnc*fdel_dc*fEbc;
	fQb = Qb(Vbe,Vbc);
	fIcc = Icc(Vbe,Vbc);
	fIec = Iec(Vbe,Vbc);
	// compute coefficients for BE junction
	ftemp1 = (fC0e)*(fH11e - fH10e) - (fC1e)*(fH01e - fH00e);
	fA2 = 0;
	fB2 = 0.5*fC3e*fne*fnc*(dIc2_Vbc(Vbe,Vbc)-dIec_Vbc(Vbe,Vbc));
	fC2 = fC2e*(ILe(Vbe)+Ic2(Vbe,Vbc)-Iec(Vbe,Vbc))-ftemp1;
	fD2 = 0;
	fE2 = 1+fC3e*(0.5*(dILe_Vbe(Vbe)+dIc2_Vbe(Vbe,Vbc)-dIec_Vbe(Vbe,Vbc))+fC6e*Ce(Vbe,Vbc));

	// compute coefficients for BC junction
	ftemp2 = (fC0c)*(fH11c - fH10c) - (fC1c)*(fH01c - fH00c);
	fA1 = 0;
	fB1 = 1+fC3c*(0.5*(dILc_Vbc(Vbc)+dIc1_Vbc(Vbe,Vbc)-dIcc_Vbc(Vbe,Vbc))+fC6c*Cc(Vbe,Vbc));
	fC1 = fC2c*(ILc(Vbc)+Ic1(Vbe,Vbc)-Icc(Vbe,Vbc))-ftemp2;
	fD1 = 0;
	fE1 = 0.5*fC3c*fne*fnc*(dIc1_Vbe(Vbe,Vbc)-dIcc_Vbe(Vbe,Vbc));

	// compute Ebc and Ebe
	fden = fB2*fE1-fB1*fE2;
	fnome = fB1*fC2-fB2*fC1;
	fnomc = fC1*fE2-fC2*fE1;

	// update E fields
	fEbe = fEbe+(fnome/fden);
	fEbc = fEbc+(fnomc/fden); 	
}

// function to update both the BE and BC junctions E field
void tEBjt::Update_Efield(tEHfdtd *pteh, int i, int j, int k, unsigned char uor, double fH01,
						  double fH00, double fH11, double fH10)
{
	if ((i==ie)&&(j==je)&&(k==ke)&&(uor==uobe)) // its BE junction
	{
		fH00e = fH00;
		fH01e = fH01;
		fH10e = fH10;
		fH11e = fH11;
	}
	else	// its BC junction
	{
		fH00c = fH00;
		fH01c = fH01;
		fH10c = fH10;
		fH11c = fH11;
	}
	
	uState++;	// increment object state
	if (uState == 2)
	{
		Update_BJTEfield();	// compute BC and BE junction fields simultaneously
		uState = 0; // reset object state
		switch (uobe) // check whether BE field is x, y or z oriented
		{
			case 0:
				pteh->fEx[ie][je][ke] = fEbe;	// update BE field
				break;
			case 1:
				pteh->fEx[ie][je][ke] = fEbe;	// update BE field
				break;
			case 2:
				pteh->fEy[ie][je][ke] = fEbe;	// update BE field
				break;
			case 3:
				pteh->fEy[ie][je][ke] = fEbe;	// update BE field
				break;
			case 4:
				pteh->fEz[ie][je][ke] = fEbe;	// update BE field
				break;
			default:
				pteh->fEz[ie][je][ke] = fEbe;	// update BE field
		}
		switch (uobc)	// check whether BC field is x, y or z oriented
		{
			case 0:
				pteh->fEx[ic][jc][kc] = fEbc;	// update BC field
				break;
			case 1:
				pteh->fEx[ic][jc][kc] = fEbc; // update BC field
				break;
			case 2:
				pteh->fEy[ic][jc][kc] = fEbc; // update BC field
				break;
			case 3:
				pteh->fEy[ic][jc][kc] = fEbc; // update BC field
				break;
			case 4:
				pteh->fEz[ic][jc][kc] = fEbc; // update BC field
				break;
			default:
				pteh->fEz[ic][jc][kc] = fEbc; // update BC field
		}
	}
}

// function to set the BE junction field 
// argument: fE = the field value
// return: none
void tEBjt::Set_BEfield(double fE)
{
	fEbe = fE;
}

// function to set the BC junction field
// argument: fE = the field value
// return: none
void tEBjt::Set_BCfield(double fE)
{
	fEbc = fE;
}

// function to calculate the total dynamic junction capacitance of the BC junction
// arguments: fV = biasing voltage
// return value: the capacitance
double tEBjt::Cc(double fVe, double fVc)
{
	double ftemp1,ftemp2;
	double fF;

	ftemp1 = fTTr*(dIec_Vbc(fVe,fVc));
	if (fVc < (fFC*fVjc))
	{
		ftemp2 = fCjc/(pow(1-(fVc/fVjc),fMc));
	}
	else
	{
		fF = 1-fFC*(1+fMc);
		ftemp2 = fF+(fMc*fVc/fVjc);
		fF = pow(1-fFC,1+fMc);
		ftemp2 = ftemp2*fCjc/fF;
	}
	return ftemp1+ftemp2;
}


// function to calculate the total dynamic junction capacitance of the BE junction
// arguments: fV = biasing voltage
// return value: the capacitance
double tEBjt::Ce(double fVe, double fVc)
{
	double ftemp1,ftemp2;
	double fF;

	ftemp1 = fTTf*(dIcc_Vbe(fVe,fVc));
	if (fVe < (fFC*fVje))
	{
		ftemp2 = fCje/(pow(1-(fVe/fVje),fMe));
	}
	else
	{
		fF = 1-fFC*(1+fMe);
		ftemp2 = fF+(fMe*fVe/fVje);
		fF = pow(1-fFC,1+fMe);
		ftemp2 = ftemp2*fCje/fF;
	}
	return ftemp1+ftemp2;
}

// function to calculate the majority carrier charge in Base (Qb) 
// capacitance of the BE junction
// arguments: fV = biasing voltage
// return value: Qb
double tEBjt::Qb(double fVe, double fVc)
{
	double fq1, fq2;
	double ftemp;

	fq1 = 1 + fVe/fVb + fVc/fVa;
	fq2 = (fIs/fIkf)*(exp(fVe/fnVte)-1) + (fIs/fIkr)*(exp(fVc/fnVtc)-1);
	ftemp = (0.25*fq1*fq1) + fq2;
	return ((0.5*fq1) + sqrt(ftemp));
}

// function to calculate dQb/dVbc
// capacitance of the BE junction
// arguments: fV = biasing voltage
// return value: dQb/dVbc
double tEBjt::dQb_Vbc(double fVe, double fVc)
{
	double fq1, ftemp;

	fq1 = 1 + fVe/fVb + fVc/fVa;
	ftemp = 0.5*(fq1/fVa) + (fIs/(fIkr*fnVtc))*exp(fVc/fnVtc);
	ftemp = ftemp/(fQb -(0.5*fq1));
	return 0.5*((1/fVa) + ftemp);
}

// function to calculate dQb/dVbe
// capacitance of the BE junction
// arguments: fV = biasing voltage
// return value: dQb/dVbe
double tEBjt::dQb_Vbe(double fVe, double fVc)
{
	double fq1, ftemp;

	fq1 = 1 + fVe/fVb + fVc/fVa;
	ftemp = 0.5*(fq1/fVb) + (fIs/(fIkf*fnVte))*exp(fVe/fnVte);
	ftemp = ftemp/(fQb-(0.5*fq1));
	return 0.5*((1/fVb) + ftemp);
}


// function to calculate the static junction current Icc as a function of biasing 
// voltage
// argument: fVe = bias voltage Vbe
//			 fVc = bias voltage Vbc
// return: static current Icc
double tEBjt::Icc(double fVe, double fVc)
{
	return (fIs/fQb)*(exp(fVe*fC5e)-1);
}

// function to calculate the differentiation of static junction current Icc as a 
// function of biasing voltage
// argument: fVe = bias voltage Vbe
//			 fVc = bias voltage Vbc
// return: dIcc/dVbe
double tEBjt::dIcc_Vbe(double fVe, double fVc)
{
	double ftemp;

	ftemp = fIs*fC5e*(exp(fVe*fC5e));
	ftemp = ftemp - dQb_Vbe(fVe,fVc)*fIcc;
	return ftemp/fQb;
}

// function to calculate the differentiation of static junction current Icc as a 
// function of biasing voltage
// argument: fVe = bias voltage Vbe
//			 fVc = bias voltage Vbc
// return: dIcc/dVbc
double tEBjt::dIcc_Vbc(double fVe, double fVc)
{
	return -dQb_Vbc(fVe,fVc)*fIcc/fQb;
}

// function to calculate the static junction current Iec as a function of biasing 
// voltage
// argument: fVe = bias voltage Vbe
//			 fVc = bias voltage Vbc
// return: static current Iec
double tEBjt::Iec(double fVe, double fVc)
{
	return (fIs/fQb)*(exp(fVc*fC5c)-1);
}

// function to calculate the differentiation of static junction current Iec as a 
// function of biasing voltage
// argument: fVe = bias voltage Vbe
//			 fVc = bias voltage Vbc
// return: dIec/dVbe
double tEBjt::dIec_Vbe(double fVe, double fVc)
{
	return (-1/fQb)*dQb_Vbe(fVe,fVc)*fIec;
}

// function to calculate the differentiation of static junction current Iec as a 
// function of biasing voltage
// argument: fVe = bias voltage Vbe
//			 fVc = bias voltage Vbc
// return: dIec/dVbc
double tEBjt::dIec_Vbc(double fVe, double fVc)
{
	double ftemp;

	ftemp = (fIs*fC5c)*exp(fVc*fC5c);
	ftemp = ftemp - dQb_Vbc(fVe,fVc)*fIec;
	return ftemp/fQb;	
}


// function to calculate the static junction current Ic1 as a function of biasing 
// voltage
// argument: fVe = bias voltage Vbe
//			 fVc = bias voltage Vbc
// return: static current Ic1
double tEBjt::Ic1(double fVe, double fVc)
{
	return (1+1/fBr)*Iec(fVe,fVc);
}

// function to calculate the differentiation of static junction current Ic1 as a 
// function of biasing voltage
// argument: fVe = bias voltage Vbe
//			 fVc = bias voltage Vbc
// return: dIc1/dVbc
double tEBjt::dIc1_Vbc(double fVe, double fVc)
{
	return (1+1/fBr)*dIec_Vbc(fVe,fVc);
}

// function to calculate the differentiation of static junction current Ic1 as a 
// function of biasing voltage
// argument: fVe = bias voltage Vbe
//			 fVc = bias voltage Vbc
// return: dIc1/dVbe
double tEBjt::dIc1_Vbe(double fVe, double fVc)
{
	return (1+1/fBr)*dIec_Vbe(fVe,fVc);
}

// function to calculate the static junction current Ic2 as a function of biasing 
// voltage
// argument: fVe = bias voltage Vbe
//			 fVc = bias voltage Vbc
// return: static current Ic2
double tEBjt::Ic2(double fVe, double fVc)
{
	return	(1+1/fBf)*Icc(fVe,fVc);
}

// function to calculate the differentiation of static junction current Ic2 as a 
// function of biasing voltage
// argument: fVe = bias voltage Vbe
//			 fVc = bias voltage Vbc
// return: dIc2/dVbe
double tEBjt::dIc2_Vbe(double fVe, double fVc)
{
	return (1+1/fBf)*dIcc_Vbe(fVe,fVc);
}

// function to calculate the differentiation of static junction current Ic2 as a 
// function of biasing voltage
// argument: fVe = bias voltage Vbe
//			 fVc = bias voltage Vbc
// return: dIc2/dVbc
double tEBjt::dIc2_Vbc(double fVe, double fVc)
{
	return (1+1/fBf)*dIcc_Vbc(fVe,fVc);
}

// function to calculate the ILc 
// function of biasing voltage
// argument: fVc = bias voltage Vbc
// return: ILc
double tEBjt::ILc(double fVc)
{
	return fIsc*(exp(fVc/fnVtcl)-1);
}

// function to calculate the differentiation of static junction current ILc as a 
// function of biasing voltage
// argument: fVc = bias voltage Vbc
// return: dILc/dVbc
double tEBjt::dILc_Vbc(double fVc)
{
	return (fIsc/fnVtcl)*(exp(fVc/fnVtcl));
}

// function to calculate the ILe 
// function of biasing voltage
// argument: fVe = bias voltage Vbe
// return: ILe
double tEBjt::ILe(double fVe)
{
	return fIse*(exp(fVe/fnVtel)-1);
}

// function to calculate the differentiation of static junction current ILe as a 
// function of biasing voltage
// argument: fVe = bias voltage Vbe
// return: dILe/dVbe
double tEBjt::dILe_Vbe(double fVe)
{
	return (fIse/fnVtel)*(exp(fVe/fnVtel));
}

// function to assign a value in the parameter list of the component to the coresponding
// variable in the object.
// arguments:
//	plist = reference to parameter list pointer
//	fP	  = reference to the variable (a double datatype)
//  fdefault = the default value
// return: none
void tEBjt::AssignVariable(lpDOUBLE &plist, double &fP, double fdefault)
{
	if (plist == NULL)	// if no valid data
		fP = fdefault;	// assign default value
	else
	{
		fP = plist->fP;
		plist = plist->next; // point to next item
	}
}
