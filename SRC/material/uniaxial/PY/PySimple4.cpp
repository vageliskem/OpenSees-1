/* *********************************************************************
**    Module:	PySimple4.cpp 
**
**    Purpose:	Provide a simple p-y spring for OpenSees.
**
**    Developed by Ross W. Boulanger
**    
** Copyright @ 2002 The Regents of the University of California (The Regents). All Rights Reserved.
**
** The Regents grants permission, without fee and without a written license agreement, for (a) use, 
** reproduction, modification, and distribution of this software and its documentation by educational, 
** research, and non-profit entities for noncommercial purposes only; and (b) use, reproduction and 
** modification of this software by other entities for internal purposes only. The above copyright 
** notice, this paragraph and the following three paragraphs must appear in all copies and modifications 
** of the software and/or documentation.
**
** Permission to incorporate this software into products for commercial distribution may be obtained 
** by contacting the University of California 
** Office of Technology Licensing 
** 2150 Shattuck Avenue #510, 
** Berkeley, CA 94720-1620, 
** (510) 643-7201.
**
** This software program and documentation are copyrighted by The Regents of the University of California. 
** The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The 
** end-user understands that the program was developed for research purposes and is advised not to rely 
** exclusively on the program for any reason.
**
** IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR 
** CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS 
** DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. REGENTS GRANTS 
** NO EXPRESS OR IMPLIED LICENSE IN ANY PATENT RIGHTS OF REGENTS BUT HAS IMPLEMENTED AN INDIVIDUAL 
** CONTRIBUTOR LICENSE AGREEMENT FOR THE OPENSEES PROJECT AT THE UNIVERISTY OF CALIFORNIA, BERKELEY 
** TO BENEFIT THE END USER.
**
** REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
** OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION,
** IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, 
** SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
**
** ****************************************************************** */

// $Revision: 1.0
// $Date: 2001/12/15
// $Source: /OpenSees/SRC/material/uniaxial/PySimple4.cpp

// Written: RWB
// Created: Dec 2001
// Revision: A
// tested and checked: Boris Jeremic (jeremic@ucdavis.edu) Spring 2002
//
// Description: This file contains the class implementation for PySimple4

#include <stdlib.h>
#include <math.h>

#include "PySimple4.h"
#include <Vector.h>
#include <Channel.h>
#include <elementAPI.h>

// Controls on internal iteration between spring components
const int PYmaxIterations = 20;
const double PYtolerance = 1.0e-12;

void* OPS_PySimple4()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 5) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial PySimple4 tag? soilType? pult? y50? drag? dashpot?\n";
	return 0;
    }
    
    int idata[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid int inputs\n";
	return 0;
    }
    
    double ddata[8] = {0,0,0,0,0,0,0,0};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 8) numdata = 8;
    if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }
    
    UniaxialMaterial *theMaterial = 0;
    theMaterial = new PySimple4(idata[0], MAT_TAG_PySimple4, idata[1], ddata[0], ddata[1],
    				ddata[2],ddata[3],ddata[4],ddata[5],ddata[6],ddata[7]);
    
    return theMaterial;
}

/////////////////////////////////////////////////////////////////////
//	Constructor with data

PySimple4::PySimple4(int tag, int classtag, int soil, double p_ult, double y_50,
				 double dragratio, double dash_pot, double s1or3 , double pyield_3 , double ke_3 , double C_3)
:UniaxialMaterial(tag,classtag),
 soilType(soil), pult(p_ult), y50(y_50), drag(dragratio), dashpot(dash_pot), simple1or3(s1or3), pyield(pyield_3), ke(ke_3), C(C_3)
{
  // Initialize PySimple variables and history variables
  //
  this->revertToStart();
  initialTangent = Ttangent;
}

/////////////////////////////////////////////////////////////////////
//	Default constructor

PySimple4::PySimple4()
:UniaxialMaterial(0,0),
 soilType(0), pult(0.0), y50(0.0), drag(0.0), dashpot(0.0), simple1or3(0.0), pyield(0.0), ke(0.0), C(0.0)
{
}

/////////////////////////////////////////////////////////////////////
//	Default destructor
PySimple4::~PySimple4()
{
    // Does nothing
}














/////////////////////////////////////////////////////////////////////
void PySimple4::getGap(double ylast, double dy, double dy_old)
{
	// For stability in Closure spring, may limit "dy" step size to avoid
	// overshooting on the closing of this gap.
	//
	TGap_y = ylast + dy;
	if(TGap_y > TClose_yright) {dy = 0.75*(TClose_yright - ylast);}
	if(TGap_y < TClose_yleft)  {dy = 0.75*(TClose_yleft  - ylast);}

	// Limit "dy" step size if it is oscillating in sign and not shrinking
	//
	if(dy*dy_old < 0.0 && fabs(dy/dy_old) > 0.5) dy = -dy_old/2.0;
	
	// Combine the Drag and Closure elements in parallel, starting by
	// resetting TGap_y in case the step size was limited.
	//
	TGap_y   = ylast + dy;
	getClosure(ylast,dy);
	getDrag(ylast,dy);
	TGap_p = TDrag_p + TClose_p;
	TGap_tang = TDrag_tang + TClose_tang;

	// Ensure that |p|<pmax.
	//
	if(fabs(TGap_p)>=pult) TGap_p =(TGap_p/fabs(TGap_p))*(1.0-PYtolerance)*pult;

	return;
}

/////////////////////////////////////////////////////////////////////
void PySimple4::getFarField(double y)
{
	TFar_y   = y;
	TFar_tang= TFar_tang;
	TFar_p   = TFar_tang * TFar_y;

	return;
}

/////////////////////////////////////////////////////////////////////
void PySimple4::getClosure(double ylast, double dy)
{
	// Reset the history terms to the last Committed values, and let them
	// reset if the reversal of loading persists in this step.
	//
	if(TClose_yleft != CClose_yleft)  TClose_yleft = CClose_yleft;
	if(TClose_yright!= CClose_yright) TClose_yright= CClose_yright;

	// Check if plastic deformation in Near Field should cause gap expansion
	//
	TClose_y = ylast + dy;
	double yrebound=1.5*y50;
	if(TNF_y+TClose_y > -TClose_yleft + yrebound)
		TClose_yleft=-(TNF_y+TClose_y) + yrebound;
	if(TNF_y+TClose_y < -TClose_yright - yrebound)
		TClose_yright=-(TNF_y+TClose_y) - yrebound;

	// Spring force and tangent stiffness
	//
	TClose_p=1.8*pult*(y50/50.0)*(pow(y50/50.0 + TClose_yright - TClose_y,-1.0)
		-pow(y50/50.0 + TClose_y - TClose_yleft,-1.0));
	TClose_tang=1.8*pult*(y50/50.0)*(pow(y50/50.0+ TClose_yright - TClose_y,-2.0)
		+pow(y50/50.0 + TClose_y - TClose_yleft,-2.0));

	// Ensure that tangent not zero or negative.
	//	
	if(TClose_tang <= 1.0e-2*pult/y50) {TClose_tang = 1.0e-2*pult/y50;}

	return;
}

/////////////////////////////////////////////////////////////////////
void PySimple4::getDrag(double ylast, double dy)
{
	TDrag_y = ylast + dy;
	double pmax=drag*pult;
	double dyTotal=TDrag_y - CDrag_y;

	// Treat as elastic if dyTotal is below PYtolerance
	//
	if(fabs(dyTotal*TDrag_tang/pult) < 10.0*PYtolerance) 
	{
		TDrag_p = TDrag_p + dy*TDrag_tang;
		if(fabs(TDrag_p) >=pmax) TDrag_p =(TDrag_p/fabs(TDrag_p))*(1.0-1.0e-8)*pmax;
		return;
	}
	// Reset the history terms to the last Committed values, and let them
	// reset if the reversal of loading persists in this step.
	//
	if(TDrag_pin != CDrag_pin)
	{
		TDrag_pin = CDrag_pin;
		TDrag_yin = CDrag_yin;
	}

	// Change from positive to negative direction
	//
	if(CDrag_y > CDrag_yin && dyTotal < 0.0)
	{
		TDrag_pin = CDrag_p;
		TDrag_yin = CDrag_y;
	}
	// Change from negative to positive direction
	//
	if(CDrag_y < CDrag_yin && dyTotal > 0.0)
	{
		TDrag_pin = CDrag_p;
		TDrag_yin = CDrag_y;
	}
	
	// Positive loading
	//
	if(dyTotal >= 0.0)
	{
		TDrag_p=pmax-(pmax-TDrag_pin)*pow(y50/2.0,nd)
					*pow(y50/2.0 + TDrag_y - TDrag_yin,-nd);
		TDrag_tang=nd*(pmax-TDrag_pin)*pow(y50/2.0,nd)
					*pow(y50/2.0 + TDrag_y - TDrag_yin,-nd-1.0);
	}
	// Negative loading
	//
	if(dyTotal < 0.0)
	{
		TDrag_p=-pmax+(pmax+TDrag_pin)*pow(y50/2.0,nd)
					*pow(y50/2.0 - TDrag_y + TDrag_yin,-nd);
		TDrag_tang=nd*(pmax+TDrag_pin)*pow(y50/2.0,nd)
					*pow(y50/2.0 - TDrag_y + TDrag_yin,-nd-1.0);
	}
	// Ensure that |p|<pmax and tangent not zero or negative.
	//
	if(fabs(TDrag_p) >=pmax) {
		TDrag_p =(TDrag_p/fabs(TDrag_p))*(1.0-PYtolerance)*pmax;}
	if(TDrag_tang <=1.0e-2*pult/y50) TDrag_tang = 1.0e-2*pult/y50;

	return;
}




///////////////////////////////////////////////////////////////////
void PySimple4::getNearField(double ylast, double dy, double dy_old) // Proper pysimple from pysimple1
{
	// Limit "dy" step size if it is oscillating in sign and not shrinking
	//
	if(dy*dy_old < 0.0 && fabs(dy/dy_old) > 0.5) dy = -dy_old/2.0;

	// Set "dy" so "y" is at middle of elastic zone if oscillation is large.
	// Note that this criteria is based on the min step size in setTrialStrain.
	//
	if(dy*dy_old < -y50*y50) dy = (TNFyinr + TNFyinl)/2.0 - ylast;
	
	// Establish trial "y" and direction of loading (with NFdy) for entire step
	//
	TNF_y = ylast + dy;
	double NFdy = TNF_y - CNF_y;

	// Treat as elastic if NFdy is below PYtolerance
	//
	if(fabs(NFdy*TNF_tang/pult) < 10.0*PYtolerance) 
	{
		TNF_p = TNF_p + dy*TNF_tang;
		if(fabs(TNF_p) >=pult) TNF_p=(TNF_p/fabs(TNF_p))*(1.0-PYtolerance)*pult;
		return;
	}

	// Reset the history terms to the last Committed values, and let them
	// reset if the reversal of loading persists in this step.
	//
	if(TNFpinr != CNFpinr || TNFpinl != CNFpinl)
	{
		TNFpinr = CNFpinr;
		TNFpinl = CNFpinl;
		TNFyinr = CNFyinr;
		TNFyinl = CNFyinl;
	}

	// For stability, may have to limit "dy" step size if direction changed.
	//
	bool changeDirection = false;
	
	// Direction change from a yield point triggers new Elastic range
	//
	double minE = 0.25;		// The min Elastic range on +/- side of p=0
	if(CNF_p > CNFpinr && NFdy <0.0){				// from pos to neg
		changeDirection = true;
		TNFpinr = CNF_p;
		if(fabs(TNFpinr)>=(1.0-PYtolerance)*pult){TNFpinr=(1.0-2.0*PYtolerance)*pult;}
		TNFpinl = TNFpinr - 2.0*pult*Elast;
		if (TNFpinl > -minE*pult) {TNFpinl = -minE*pult;}
		TNFyinr = CNF_y;
		TNFyinl = TNFyinr - (TNFpinr-TNFpinl)/NFkrig; 
	}
	if(CNF_p < CNFpinl && NFdy > 0.0){				// from neg to pos
		changeDirection = true;
		TNFpinl = CNF_p;
		if(fabs(TNFpinl)>=(1.0-PYtolerance)*pult){TNFpinl=(-1.0+2.0*PYtolerance)*pult;}
		TNFpinr = TNFpinl + 2.0*pult*Elast;
		if (TNFpinr < minE*pult) {TNFpinr = minE*pult;}
		TNFyinl = CNF_y;
		TNFyinr = TNFyinl + (TNFpinr-TNFpinl)/NFkrig; 
	}
	// Now if there was a change in direction, limit the step size "dy"
	//
	if(changeDirection == true) {
		double maxdy = 0.25*pult/NFkrig;
		if(fabs(dy) > maxdy) dy = (dy/fabs(dy))*maxdy;
	}

	// Now, establish the trial value of "y" for use in this function call.
	//
	TNF_y = ylast + dy;

	// Postive loading
	//
	if(NFdy >= 0.0){
		// Check if elastic using y < yinr
		if(TNF_y <= TNFyinr){							// stays elastic
			TNF_tang = NFkrig;
			TNF_p = TNFpinl + (TNF_y - TNFyinl)*NFkrig;
		}
		else {
			TNF_tang = np * (pult-TNFpinr) * pow(yref,np) 
				* pow(yref - TNFyinr + TNF_y, -np-1.0);
			TNF_p = pult - (pult-TNFpinr)* pow(yref/(yref-TNFyinr+TNF_y),np);
		}
	}

	// Negative loading
	//
	if(NFdy < 0.0){
		// Check if elastic using y < yinl
		if(TNF_y >= TNFyinl){							// stays elastic
			TNF_tang = NFkrig;
			TNF_p = TNFpinr + (TNF_y - TNFyinr)*NFkrig;
		}
		else {
			TNF_tang = np * (pult+TNFpinl) * pow(yref,np) 
				* pow(yref + TNFyinl - TNF_y, -np-1.0);
			TNF_p = -pult + (pult+TNFpinl)* pow(yref/(yref+TNFyinl-TNF_y),np);
		}
	}

	// Ensure that |p|<pult and tangent not zero or negative.
	//
	if(fabs(TNF_p) >=pult) TNF_p=(TNF_p/fabs(TNF_p))*(1.0-PYtolerance)*pult;
	if(TNF_tang <= 1.0e-2*pult/y50) TNF_tang = 1.0e-2*pult/y50;

    return;
}
























// /////////////////////////////////////////////////////////////////////
// void PySimple4::getNearField_3(double ylast, double dy, double dy_old) //get Near Field according to PYsimple1
// {
// 	// Limit "dy" step size if it is oscillating in sign and not shrinking
// 	//
// 	if(dy*dy_old < 0.0 && fabs(dy/dy_old) > 0.5) dy = -dy_old/2.0;

// 	// Set "dy" so "y" is at middle of elastic zone if oscillation is large.
// 	// Note that this criteria is based on the min step size in setTrialStrain.
// 	//
// 	if(dy*dy_old < -y50*y50) dy = (TNFyinr + TNFyinl)/2.0 - ylast;
	
// 	// Establish trial "y" and direction of loading (with NFdy) for entire step
// 	//
// 	TNF_y = ylast + dy;  //ylast commited y and dy _V
// 	double NFdy = TNF_y - CNF_y;

// 	// Treat as elastic if NFdy is below PYtolerance
// 	//
// 	if(fabs(NFdy*TNF_tang/pult) < 10.0*PYtolerance) 
// 	{
// 		TNF_p = TNF_p + dy*TNF_tang;
// 		if(fabs(TNF_p) >=pult) TNF_p=(TNF_p/fabs(TNF_p))*(1.0-PYtolerance)*pult;
// 		return;
// 	}

// 	// Reset the history terms to the last Committed values, and let them
// 	// reset if the reversal of loading persists in this step.
// 	//
// 	if(TNFpinr != CNFpinr || TNFpinl != CNFpinl)
// 	{
// 		TNFpinr = CNFpinr;
// 		TNFpinl = CNFpinl;
// 		TNFyinr = CNFyinr;
// 		TNFyinl = CNFyinl;
// 	}

// 	// For stability, may have to limit "dy" step size if direction changed.
// 	//
// 	bool changeDirection = false;
	
// 	// Direction change from a yield point triggers new Elastic range
// 	//
// 	double minE = 0.25;		// The min Elastic range on +/- side of p=0
// 	if(CNF_p > CNFpinr && NFdy <0.0){				// from pos to neg direction
// 		changeDirection = true;
// 		TNFpinr = CNF_p;
// 		if(fabs(TNFpinr)>=(1.0-PYtolerance)*pult){TNFpinr=(1.0-2.0*PYtolerance)*pult;}
// 		TNFpinl = TNFpinr - 2.0*pult*Elast;
// 		if (TNFpinl > -minE*pult) {TNFpinl = -minE*pult;}
// 		TNFyinr = CNF_y;
// 		TNFyinl = TNFyinr - (TNFpinr-TNFpinl)/NFkrig; 
// 		// TNFyinl = TNFyinr - (TNFpinr-TNFpinl)/ke; \\ _v
// 	}
// 	if(CNF_p < CNFpinl && NFdy > 0.0){				// from neg to pos direction 
// 		changeDirection = true;
// 		TNFpinl = CNF_p;
// 		if(fabs(TNFpinl)>=(1.0-PYtolerance)*pult){TNFpinl=(-1.0+2.0*PYtolerance)*pult;}
// 		TNFpinr = TNFpinl + 2.0*pult*Elast;
// 		if (TNFpinr < minE*pult) {TNFpinr = minE*pult;}
// 		TNFyinl = CNF_y;
// 		TNFyinr = TNFyinl + (TNFpinr-TNFpinl)/NFkrig; 
// 		// TNFyinr = TNFyinl + (TNFpinr-TNFpinl)/ke;  \\ _v
// 	}
// 	// Now if there was a change in direction, limit the step size "dy"
// 	//
// 	if(changeDirection == true) {
// 		double maxdy = 0.25*pult/NFkrig;
// 		// double maxdy = 0.25*pult/ke; // _v

// 		if(fabs(dy) > maxdy) dy = (dy/fabs(dy))*maxdy;
// 	}

// 	// Now, establish the trial value of "y" for use in this function call.
// 	//
// 	TNF_y = ylast + dy;

// 	// Postive loading
// 	//
// 	if(NFdy >= 0.0){
// 		// Check if elastic using y < yinr
// 		if(TNF_y <= TNFyinr){							// stays elastic
// 			// NFkrig  = 100.0 * (0.5 * pult) / y50
// 			// TNF_tang = NFkrig;
// 			// TNF_p = TNFpinl + (TNF_y - TNFyinl)*NFkrig;
// 			TNF_tang = ke;
// 			TNF_p = TNFpinl + (TNF_y - TNFyinl)*ke;

// 		}
// 		else {

// 			TNF_tang = 1.0*np * (pult-TNFpinr) * pow(yref,np) 
// 				* pow(yref - TNFyinr + TNF_y, -np-1.0);
// 			TNF_p = pult - (pult-TNFpinr)* pow(yref/(yref-TNFyinr+TNF_y),np); // This is the first equation in the 
// 			// paper  appendix. according to this yref=c*y50 with c=10 for clay and 0.5 for sand; _V
			
// 			// opserr << "WARNING" << simple1or3 << endln;

// 		}
// 	}

// 	// Negative loading
// 	//
// 	if(NFdy < 0.0){
// 		// Check if elastic using y < yinl
// 		if(TNF_y >= TNFyinl){							// stays elastic
// 			// TNF_tang = NFkrig;
// 			// TNF_p = TNFpinr + (TNF_y - TNFyinr)*NFkrig;
// 			TNF_tang = ke;
// 			TNF_p = TNFpinr + (TNF_y - TNFyinr)*ke;
// 		}
// 		else {
// 			// TNF_tang = NFkrig;
// 			// TNF_p = TNFpinr + (TNF_y - TNFyinr)*NFkrig;
// 			np=np*1.0;
// 			TNF_tang = np * (pult+TNFpinl) * pow(yref,np) 
// 				* pow(yref + TNFyinl - TNF_y, -np-1.0);
// 			TNF_p = -pult + (pult+TNFpinl)* pow(yref/(yref+TNFyinl-TNF_y),np);

// 		}
// 	}

// 	// Ensure that |p|<pult and tangent not zero or negative.
// 	//
// 	if(fabs(TNF_p) >=pult) TNF_p=(TNF_p/fabs(TNF_p))*(1.0-PYtolerance)*pult;
// 	if(TNF_tang <= 1.0e-2*pult/y50) TNF_tang = 1.0e-2*pult/y50;

//     return;
// }















/////////////////////////////////////////////////////////////////////
void PySimple4::getNearField_3(double ylast, double dy, double dy_old) //get Near Field according to PYsimple1
{
	// Limit "dy" step size if it is oscillating in sign and not shrinking
	//
	if(dy*dy_old < 0.0 && fabs(dy/dy_old) > 0.5) dy = -dy_old/2.0;

	// Set "dy" so "y" is at middle of elastic zone if oscillation is large.
	// Note that this criteria is based on the min step size in setTrialStrain.
	//
	if(dy*dy_old < -y50*y50) dy = (TNFyinr + TNFyinl)/2.0 - ylast;
	
	// Establish trial "y" and direction of loading (with NFdy) for entire step
	//
	TNF_y = ylast + dy;  //ylast commited y and dy _V
	double NFdy = TNF_y - CNF_y;

	// Treat as elastic if NFdy is below PYtolerance
	//
	if(fabs(NFdy*TNF_tang/pult) < 10.0*PYtolerance) 
	{
		TNF_p = TNF_p + dy*TNF_tang;
		if(fabs(TNF_p) >=pult) TNF_p=(TNF_p/fabs(TNF_p))*(1.0-PYtolerance)*pult;
		      // opserr << "WARNING P1=" << P1 << endln;

		return;
	}
	// else
	// {	
	// 	TNF_p = TNF_p + dy*TNF_tang; 
	// 	return;
	// }



		
		
			// HEREEE
	
	
		signdy = sign(dy);	//note sign(dy) gives same value as sign(yRate)
		signdyLast = sign(dyLast);
		TpinF = CpinF;
		TpinB = CpinB;
		double CpinFB4check = CpinF;
		double CpinBB4check = CpinB;
		TpinUse = CpinLast;
		Tyin = Cyin;
		TLastYieldDir = CLastYieldDir;
		int yRate3 = 0.0; //set previously defined yRate tp yRate3 and set it to zero;





  if(yRate3 == 0.0)
    {
      // P1 = Cp+ke*dy;	//this is the viscoelastic force predictor which assumes all the displacement is elastic. If yielding occurs, a revised version of this equation will be used which considers only the viscoeslatic portion of the deformation.
      TNF_p = TNF_p+ke*dy;	//this is the viscoelastic force predictor which assumes all the displacement is elastic. If yielding occurs, a revised version of this equation will be used which considers only the viscoeslatic portion of the deformation.
      // opserr << "WARNING TNF_p=" << TNF_p << endln;
      // return;
    }
  else
    {
      TNF_p = TNF_p+ke*dy+((dashpot/tstep1)*dy)-((dashpot/tstep2)*dyELast);	
    }
  if (signdy!=0 && signdy != TLastYieldDir && sign(TNF_p-Cp)!=signdy) {
    TNF_p = TNF_p+(0.000001*signdy);	//this ensures that when the displacement reverses, the load does not continue in the previous direction. This would only occur in unusual situations where deceleration yielding or perhaps a load reversal had also occured during the previous step.	
  }	
  Tp = TNF_p;
  // 	Ttangent
  if(yRate3 == 0.0)
    {
      Ttangent = ke;
    }
  else if (dy != 0.0 && dashpot == 0.0) {
    Ttangent = ke;	
  }
  else if (dy != 0.0) {
    Ttangent = ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*dy)));
  }
  else if (dy == 0.0 && dashpot != 0.0 && sign(Tp-Cp) != 0.0){
    // for case where dy is zero, if previous dashpot force was nonzero, the tangent will either be pos. or neg. infinity, since the force will
    // change without y changing. This is handled as a special case by the code below, but for now the tangent is saved as ke
    Ttangent = Ctangent;
  }
  else {
    Ttangent = ke;
  }
  Tpalpha = Cpalpha;
  f = fabs(Tp-Cpalpha)-pyield;	//yield function
  if (f <= 0.0)
    {
      TyeTotal =  (Tp-Cp+ke*CyeTotal+dashpot*CyeTotal/tstep1+dashpot*dyELast/tstep2)/(ke+dashpot/tstep1);
      TdyE    = TyeTotal-CyeTotal;
      return;	//no yielding
    }
	






  // Return mapping if yielding occurs
  if(f>0)
    {
      double CpUse = Cp;
      int signdyLast = sign(dyLast);
      int signP = sign(Tp-Cp);
      double dyELastUse = dyELast;
      double CyeTotalUse = CyeTotal;
      double bump = 0.0;
      /////If the current P is very close to Pult and displacement increment is in same dir. as last step,
      // // or the displacement increment is zero but the viscoelastic guess still lands above Pult,
      // // there can be convergence issues because the residual for the correct value of next P gets very
      // // small. Instead just assign Pult minus tol. as next P.
      if(  (fabs(Cp) > 0.99*pult && fabs(TNF_p) > fabs(pult) && signdy == signdyLast) || (fabs(Cp) > 0.99*pult && fabs(P1) > fabs(pult) && signdy == 0)  )
	{
	  Tp = (pult-PYtolerance)*signdy;
	  if (signP == 1) {TpinUse = TpinF;}
	  else if (signP == -1) {TpinUse = TpinB;}
	  else {TpinUse = CpinLast;}				
	  TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
	  TdyE    = TyeTotal-CyeTotalUse;
	  //for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occured)
	  if(dy==0.0 && TdyE == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else if (dy==0.0 && dashpot != 0.0){
	    Ttangent = Ctangent;
	  }
	  else if (dy == 0.0 && dashpot == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else if (TdyE == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else{
	    Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
	  }
	  Tpalpha = Tp - pyield*signdy;
	  // Compute viscoelastic displacement increment
	  TLastYieldDir = signP;
	  return;	
	  }

	


























// Test for deceleration yielding
      //Special case... if displacement increment is in same dir. as last converged step but velocity is much less
      // than the velocity in the previous step (or zero), the decreased force in the damper component could cause the
      // viscoelastic force to decrease relative to Cp (this can only occur if the decrease in dashpot force is large
      // relative to the increase in force in the elastic component, so with a relatively large dashpot coeff. and large 
      // deceleration). This is somewhat counterituitive, because you would expect continuing to displace in the same
      // dir. would keep increasing the load, but a large loss of damper force could cause the cummulative load to drop.
      // If this drop is large enough, yielding could occur in the opposite direction of displacement, so the previous 
      // if statements for re-assigning Tpin and Tyin won't work because they rely on sign(dy). Instead, need to use
      // the term sign(P1-Cp), where P1 is the viscoelastic predictor guess.
      //	If sign(P1-Cp) and signdy have opposite signs, this condition has occured...
	
      if(signP != signdy && signP != 0 && dy == 0.0 && signP != CLastYieldDir)
	// special case for dy = zero
	{
	  //compute new Tpin based on signP instead of signdy...
	  TpinUse = Cpalpha + pyield*signP;
	  Tyin = Cy;	
	  //Compute the change in viscoelastic displacement that occurs when the decrease in dashpot force causes P
	  //	to drop from the current P to the yield surface. First compute time this takes by proportioning forces:
	  TLastYieldDir = signP;
	  pn1_a = TpinUse;									//lower bracket is current commited value of p
	  pn1_b = (1.0-0.5*PYtolerance)*signP*pult;			//upper bracket is pult, but with opposite sign
	  //...then don't change anything, just send in the regular values...
	  CpUse = TpinUse;
	  tstep1 = tstep;
	  tstep2 = tstep;
	  dyELastUse = dyELast;
	  CyeTotalUse = CyeTotal;
	  bump=Cp-TpinUse;
	}
      else if(signP != signdy && signP != 0 && signP != CLastYieldDir)
	//	Decleration yielding for a non-zero dy case, still flip the search direction:
	{
	  TpinUse = Cpalpha - pyield*signdy;
	  Tyin = Cy + (TpinUse-Cp)/(-1.0*((Tp-Cp)/dy));
	  CpUse = TpinUse;
	  TLastYieldDir = (-1)*signdy;
	  pn1_a = TpinUse;
	  pn1_b = (1.0-0.5*PYtolerance)*(-1.0)*signdy*pult;	//upper bracket is pult, but with opposite sign
	  bump = Cp-TpinUse;
	}
      else
	{
	  // if no deceleration yielding, only re-assign Tpin if previous yield was not in the current direction
	  if(TLastYieldDir != signdy && signdy!= 0)
	    {
	      TpinUse = Cpalpha + pyield*signdy;
	      Tyin = Cy + (TpinUse-Cp)/((Tp-Cp)/dy);
	    }
	  pn1_a = Cp;											//lower bracket is current commited value of p
	  if(signP == 0)
	    {
	      pn1_b = (1.0-0.5*PYtolerance)*CLastYieldDir*pult;			//upper bracket is pult
	      TLastYieldDir = CLastYieldDir;
	    }
	  else{
	    pn1_b = (1.0-0.5*PYtolerance)*signP*pult;			//upper bracket is pult
	    TLastYieldDir = signP;
	  }
	}
	
      // Bumping routine and update backstress
      //If decleration yielding occured, bumped values have already been computed and this loop will not be entered.
      if((pn1_a < TpinUse && pn1_b > TpinUse) || (pn1_a > TpinUse && pn1_b < TpinUse)) // test to see if brackets are on opposite sides of yield surface
	{
	  pn1_a = TpinUse;
	  CpUse = TpinUse;
	  bump = TpinUse-Cp;	
	}	
      //Determine if backstress should be updated
      double TpinUseB4check = TpinUse;
      if (signP == 1 && TpinUse < CpinF) {
	TpinF = TpinUse;
      }
      else if (signP ==1 && TpinUse >= CpinF && CLastYieldDir != 0){
	TpinUse = CpinF;
	TpinF = CpinF;
      }
      else if (signP ==1 && TpinUse >= CpinF && CLastYieldDir == 0){
	TpinF = TpinUse;	
      }
      if (signP == -1 && TpinUse > CpinB) {
	TpinB = TpinUse;
      }
      else if (signP ==-1 && TpinUse <= CpinB && CLastYieldDir != 0){
	TpinUse = CpinB;
	TpinB = CpinB;
      }
      else if (signP ==-1 && TpinUse <= CpinB && CLastYieldDir == 0){
	TpinB = TpinUse;
      }		
      if (signP == 0){
	TpinUse = CpinLast;
      }

      //Compute residuals 	
      R1 = getResidual(ke,CpUse,pn1_a,dy,pult,C,TpinUse,dashpot,tstep1,dyELastUse,CyeTotalUse,tstep2,TNF_p,bump);
      R2 = getResidual(ke,CpUse,pn1_b,dy,pult,C,TpinUse,dashpot,tstep1,dyELastUse,CyeTotalUse,tstep2,P1,bump);
		
      // Iterate to find P next			
      //first check if either of the residuals computed at the starting brackets satisfies the tolerance			
      if(fabs(R1)<PYtolerance*pult)
	{
	  Tp = pn1_a;
	  TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
	  TdyE    = TyeTotal-CyeTotalUse;
	  //for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occured)
	  if(dy==0.0 && TdyE == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else if (dy==0.0 && dashpot != 0.0){
	    Ttangent = Ctangent;
	  }
	  else if ((dy == 0.0 && dashpot == 0.0) || TdyE == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else{
	    Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
	  }
	  //check if decleration cause a decrease in force
	  signP = sign(Tp-Cp);
	  if(signP != signdy)
	    {
	      signPalphaNew = (-1)*signdy;
	    }
	  else 
	    {
	      signPalphaNew = signdy;
	    }
	  Tpalpha = Tp - pyield*signPalphaNew;
	  return;
	}
      if(fabs(R2)<PYtolerance*pult)
	{
	  Tp = pn1_b;
	  TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
	  TdyE    = TyeTotal-CyeTotalUse;
	  //for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occured)
	  if(dy==0.0 && TdyE == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else if (dy==0.0 && dashpot != 0.0){
	    Ttangent = Ctangent;
	  }
	  else if ((dy == 0.0 && dashpot == 0.0) || TdyE == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else{
	    Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
	  }
	  //check if decleration cause in decrease in force
	  signP = sign(Tp-Cp);
	  if(signP != signdy)
	    {
	      signPalphaNew = (-1)*signdy;
	    }
	  else 
	    {
	      signPalphaNew = signdy;
	    }
	  Tpalpha = Tp - pyield*signPalphaNew;
	  return;
	}
			
      // If R1 and R2 weren't within tolerance, check to see if they have the same sign...
      if(sign(R1)==sign(R2))
	{
	  // 	Check for a condition in which, for a very large trial displacement step (keep in mind the solver may send a trial displacement step that is significantly larger than the eventual converged disp. step), the state goes from being either elastic or mostly elastic all the way to the yield surface. If this jump happens from far away from the yield surface, the residual for the upper bound guess can be very large and won't satisfy tolerance, even though it is essentially the correct guess. Test for this condition by nudging the upper bound guess even closer to pult and seeing if the residual continues to decrease.
	  pn1_b = (1.0-0.0000000001*PYtolerance)*TLastYieldDir*pult;
	  R2 = getResidual(ke,CpUse,pn1_b,dy,pult,C,TpinUse,dashpot,tstep1,dyELastUse,CyeTotalUse,tstep2,TNF_p,bump);
	  if(fabs(R2)<PYtolerance*pult)
	    {
	      Tp = (1.0-0.5*PYtolerance)*TLastYieldDir*pult;	//keep consistent with other cases to avoid tiny fluctuations near the yield surface
	      TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
	      TdyE    = TyeTotal-CyeTotalUse;
	      //for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occured)
	      if(dy==0.0 && TdyE == 0.0){
		Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	      }
	      else if (dy==0.0 && dashpot != 0.0){
		Ttangent = Ctangent;
	      }
	      else if ((dy == 0.0 && dashpot == 0.0) || TdyE == 0.0){
		Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	      }
	      else{
		Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
	      }
	      //check if decleration cause in decrease in force
	      signP = sign(Tp-Cp);
	      if(signP != signdy)
		{
		  signPalphaNew = (-1)*signdy;
		}
	      else 
		{
		  signPalphaNew = signdy;
		}
	      Tpalpha = Tp - pyield*signPalphaNew;
	      return;
	    }	
	}

      //start Ridder's loop, initialize variables
      double pn1_3 = 0.0;
      double pn1_4 = 0.0;
      double R3 = 0.0;
      double R4 = 0.0;
      double S = 0.0;
      int i=0;
      while((fabs(R1)>PYtolerance*pult) && (fabs(R2)>PYtolerance*pult))
	{
	  pn1_3 = (pn1_a+pn1_b)/2.0;
	  R3    = getResidual(ke,CpUse,pn1_3,dy,pult,C,TpinUse,dashpot,tstep1,dyELastUse,CyeTotalUse,tstep2,TNF_p,bump);
	  if(fabs(R3)<PYtolerance*pult)
	    {
	      Tp = pn1_3;
	      TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
	      TdyE    = TyeTotal-CyeTotalUse;
	      //for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occured)
	      if(dy==0.0 && TdyE == 0.0){
		Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	      }
	      else if (dy==0.0 && dashpot != 0.0){
		Ttangent = Ctangent;
	      }
	      else if ((dy == 0.0 && dashpot == 0.0) || TdyE == 0.0){
		Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	      }
	      else{
		Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
	      }
	      //check if decleration causes in decrease in force
	      signP = sign(Tp-Cp);
	      if(signP != signdy)
		{
		  signPalphaNew = (-1)*signdy;
		}
	      else 
		{
		  signPalphaNew = signdy;
		}
	      Tpalpha = Tp - pyield*signPalphaNew;
	      return ;
	    }
	  S = (sqrt(R3*R3 - R1*R2));
	  pn1_4 = pn1_3+(pn1_3-pn1_a)*((sign(R1-R2)*R3)/(sqrt(R3*R3 - R1*R2)));
	  R4    = getResidual(ke,CpUse,pn1_4,dy,pult,C,TpinUse,dashpot,tstep1,dyELastUse,CyeTotalUse,tstep2,TNF_p,bump);	
	  if(fabs(R4)<PYtolerance*pult)
	    {
	      Tp = pn1_4;
	      TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
	      TdyE    = TyeTotal-CyeTotalUse;
	      //for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occured)
	      if(dy==0.0 && TdyE == 0.0){
		Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	      }
	      else if (dy==0.0 && dashpot != 0.0){
		Ttangent = Ctangent;
	      }
	      else if ((dy == 0.0 && dashpot == 0.0) || TdyE == 0.0){
		Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	      }
	      else{
		Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
	      }
	      //check if decleration cause in decrease in force
	      signP = sign(Tp-Cp);
	      if(signP != signdy)
		{
		  signPalphaNew = (-1)*signdy;
		}
	      else 
		{
		  signPalphaNew = signdy;
		}
	      Tpalpha = Tp - pyield*signPalphaNew;
	      return ;
	    }
	  if(sign(R3) != sign(R4))
	    {
	      pn1_a = pn1_3;
	      pn1_b = pn1_4;
	      R1 = R3;
	      R2 = R4;
	    }
	  else
	    {
	      if(sign(R1) != sign(R4))
		{
		  pn1_a = pn1_a;
		  pn1_b = pn1_4;
		  R1 = R1;
		  R2 = R4;
		}
	      else
		{
		  if(sign(R2) != sign(R4))
		    {
		      pn1_a = pn1_4;
		      pn1_b = pn1_b;
		      R1 = R4;
		      R2 = R2;
		    }
		  else
		    {
		      //printf ("none of the Ridder's values had opposite signs-- error! \n");*/
		    }
		}
	    }
	  i++;
	  if(i==PYmaxIterations){
					
	    if(fabs(pn1_a) > 0.995*pult && fabs(pn1_b) > 0.995*pult)
	      {
		Tp = (pn1_a+pn1_b)/2.0;
		if(fabs(Tp) >= pult){
		  Tp = (pult-PYtolerance)*signdy;
		}
		TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
		TdyE    = TyeTotal-CyeTotalUse;
		//for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occured)
		if(dy==0.0 && TdyE == 0.0){
		  Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
		}
		else if (dy==0.0 && dashpot != 0.0){
		  Ttangent = Ctangent;
		}
		else if ((dy == 0.0 && dashpot == 0.0) || TdyE == 0.0){
		  Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
		}
		else{
		  Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
		}
		//check if decleration cause in decrease in force
		signP = sign(Tp-Cp);
		if(signP != signdy)
		  {
		    signPalphaNew = (-1)*signdy;
		  }
		else 
		  {
		    signPalphaNew = signdy;
		  }
		Tpalpha = Tp - pyield*signPalphaNew;
		// Compute viscoelastic displacement increment
		return ;	
	      }
	    else {
	      opserr << "Ridder's method for material tag " << this->getTag() << " failed to find a working value for P in " << PYmaxIterations << " iterations." << endln;
	      Tp = Cp;
	      Ttangent = Ctangent;
	      Tpalpha = Cpalpha;
	      TdyE = dy;
	      return ;
	    }
	  }
	}	
			
			
    }	//end return mapping algorithm
  return ;


}
















/////////////////////////////////////////////////////////////////////
int 
PySimple4::setTrialStrain (double newy, double yRate)
{
	// Set trial values for displacement and load in the material
	// based on the last Tangent modulus.
	//
	double dy = newy - Ty;
	double dp = Ttangent * dy;
	TyRate    = yRate;

	// Limit the size of step (dy or dp) that can be imposed. Prevents
	// numerical difficulties upon load reversal at high loads
	// where a soft loading modulus becomes a stiff unloading modulus.
	//
	int numSteps = 1;
	double stepSize = 1.0;
	if(fabs(dp/pult) > 0.5) numSteps = 1 + int(fabs(dp/(0.5*pult)));
	if(fabs(dy/y50)  > 1.0 ) numSteps = 1 + int(fabs(dy/(1.0*y50)));
	stepSize = 1.0/float(numSteps);
	if(numSteps > 100) numSteps = 100;

	dy = stepSize * dy;

	




	// Main loop over the required number of substeps
	
	for(int istep=1; istep <= numSteps; istep++)
	{
		Ty = Ty + dy;
		dp = Ttangent * dy;
		
	// May substep within Gap or NearField element if oscillating, which can happen
	// when they jump from soft to stiff.
	//
		double dy_gap_old = ((Tp + dp) - TGap_p)/TGap_tang;
		double dy_nf_old  = ((Tp + dp) - TNF_p) /TNF_tang;
		// double dy_nf_old  = ((Tp + dp) - TNF_p) /ke; // Vag

	// Iterate to distribute displacement among the series components.
	// Use the incremental iterative strain & iterate at this strain.
	//
	for (int j=1; j < PYmaxIterations; j++)
	{
		Tp = Tp + dp;

		// Stress & strain update in Near Field element
		double dy_nf = (Tp - TNF_p)/TNF_tang;
		// double dy_nf = (Tp - TNF_p)/ke; // Vag


		
		if(simple1or3==3) 
			{
		    opserr << "PySimple3 mat" << endln;
			getNearField_3(TNF_y,dy_nf,dy_nf_old); // for pysimple3
								// getNearField_3(TNF_y,dy_nf,dy_nf_old); // for pysimple3

			}
		else
			{
		    opserr << "PySimple1 mat" << endln;
			getNearField(TNF_y,dy_nf,dy_nf_old);
			} // for pysimple1



		// Residuals in Near Field element
		double p_unbalance = Tp - TNF_p;
		double yres_nf = (Tp - TNF_p)/TNF_tang;
		dy_nf_old = dy_nf;

		// Stress & strain update in Gap element
		double dy_gap = (Tp - TGap_p)/TGap_tang;
		getGap(TGap_y,dy_gap,dy_gap_old);

		// Residuals in Gap element
		double p_unbalance2 = Tp - TGap_p;
		double yres_gap = (Tp - TGap_p)/TGap_tang;
		dy_gap_old = dy_gap;

		// Stress & strain update in Far Field element
		double dy_far = (Tp - TFar_p)/TFar_tang;
		TFar_y = TFar_y + dy_far;
		getFarField(TFar_y);

		// Residuals in Far Field element
		double p_unbalance3 = Tp - TFar_p;
		double yres_far = (Tp - TFar_p)/TFar_tang;

		// Update the combined tangent modulus
		Ttangent = pow(1.0/TGap_tang + 1.0/TNF_tang + 1.0/TFar_tang, -1.0);

		// Residual deformation across combined element
		double dv = Ty - (TGap_y + yres_gap)
			- (TNF_y + yres_nf) - (TFar_y + yres_far);

		// Residual "p" increment 
		dp = Ttangent * dv;

		// Test for convergence
		double psum = fabs(p_unbalance) + fabs(p_unbalance2) + fabs(p_unbalance3);
		if(psum/pult < PYtolerance) break;
	}
	}

	return 0;
}









/////////////////////////////////////////////////////////////////////
double 
PySimple4::getStress(void)
{
	// Dashpot force is only due to velocity in the far field.
	// If converged, proportion by Tangents.
	// If not converged, proportion by ratio of displacements in components.
	//
	double ratio_disp =(1.0/TFar_tang)/(1.0/TFar_tang + 1.0/TNF_tang + 1.0/TGap_tang);
	if(Ty != Cy) {
		ratio_disp = (TFar_y - CFar_y)/(Ty - Cy);
		if(ratio_disp > 1.0) ratio_disp = 1.0;
		if(ratio_disp < 0.0) ratio_disp = 0.0;
	}
	double dashForce = dashpot * TyRate * ratio_disp;

	// Limit the combined force to pult.
	//
	if(fabs(Tp + dashForce) >= (1.0-PYtolerance)*pult)
		return (1.0-PYtolerance)*pult*(Tp+dashForce)/fabs(Tp+dashForce);
	else return Tp + dashForce;
}
/////////////////////////////////////////////////////////////////////
double 
PySimple4::getTangent(void)
{
    return this->Ttangent;
}
/////////////////////////////////////////////////////////////////////
double 
PySimple4::getInitialTangent(void)
{
    return this->initialTangent;
}
/////////////////////////////////////////////////////////////////////
double 
PySimple4::getDampTangent(void)
{
	// Damping tangent is produced only by the far field component.
	// If converged, proportion by Tangents.
	// If not converged, proportion by ratio of displacements in components.
	//
	double ratio_disp =(1.0/TFar_tang)/(1.0/TFar_tang + 1.0/TNF_tang + 1.0/TGap_tang);
	if(Ty != Cy) {
		ratio_disp = (TFar_y - CFar_y)/(Ty - Cy);
		if(ratio_disp > 1.0) ratio_disp = 1.0;
		if(ratio_disp < 0.0) ratio_disp = 0.0;
	}

	double DampTangent = dashpot * ratio_disp;

	// Minimum damping tangent referenced against Farfield spring
	//
	if(DampTangent < TFar_tang * 1.0e-12) DampTangent = TFar_tang * 1.0e-12;

	// Check if damping force is being limited
	//
	double totalForce = Tp + dashpot * TyRate * ratio_disp;
	if(fabs(totalForce) >= (1.0-PYtolerance)*pult) DampTangent = 0.0;

	return DampTangent;
}
/////////////////////////////////////////////////////////////////////
double 
PySimple4::getStrain(void)
{
    return this->Ty;
}
/////////////////////////////////////////////////////////////////////
double 
PySimple4::getStrainRate(void)
{
    return this->TyRate;
}
/////////////////////////////////////////////////////////////////////
int 
PySimple4::commitState(void)
{
	// Commit trial history variable -- Combined element
    Cy       = Ty;
    Cp       = Tp;
    Ctangent = Ttangent;
    
	// Commit trial history variables for Near Field component
	CNFpinr   = TNFpinr;
	CNFpinl   = TNFpinl; 
	CNFyinr   = TNFyinr;
	CNFyinl   = TNFyinl;	
	CNF_p     = TNF_p;
	CNF_y     = TNF_y;
	CNF_tang  = TNF_tang;

	// Commit trial history variables for Drag component
	CDrag_pin = TDrag_pin;
	CDrag_yin = TDrag_yin;
	CDrag_p   = TDrag_p;
	CDrag_y   = TDrag_y;
	CDrag_tang= TDrag_tang;

	// Commit trial history variables for Closure component
	CClose_yleft  = TClose_yleft;
	CClose_yright = TClose_yright;
	CClose_p      = TClose_p;
	CClose_y      = TClose_y;
	CClose_tang   = TClose_tang;

	// Commit trial history variables for the Gap
	CGap_y    = TGap_y;
	CGap_p    = TGap_p;
	CGap_tang = TGap_tang;
    
	// Commit trial history variables for the Far Field
	CFar_y    = TFar_y;
	CFar_p    = TFar_p;
	CFar_tang = TFar_tang;
    
    return 0;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple4::revertToLastCommit(void)
{
  // Reset to committed values
  
  Ty       = Cy;
  Tp       = Cp;
  Ttangent = Ctangent;
  
  
  TNFpinr   = CNFpinr;
  TNFpinl   = CNFpinl; 
  TNFyinr   = CNFyinr;
  TNFyinl   = CNFyinl;	
  TNF_p     = CNF_p;
  TNF_y     = CNF_y;
  TNF_tang  = CNF_tang;
  
  TDrag_pin = CDrag_pin;
  TDrag_yin = CDrag_yin;
  TDrag_p   = CDrag_p;
  TDrag_y   = CDrag_y;
  TDrag_tang= CDrag_tang;
  
  TClose_yleft  = CClose_yleft;
  TClose_yright = CClose_yright;
  TClose_p      = CClose_p;
  TClose_y      = CClose_y;
  TClose_tang   = CClose_tang;
  
  TGap_y    = CGap_y;
  TGap_p    = CGap_p;
  TGap_tang = CGap_tang;
  
  TFar_y    = CFar_y;
  TFar_p    = CFar_p;
  TFar_tang = CFar_tang;
  
  return 0;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple4::revertToStart(void)
{

	// If soilType = 0, then it is entering with the default constructor.
	// To avoid division by zero, set small nonzero values for terms.
	//
	if(soilType == 0){
		pult = 1.0e-12;
		y50  = 1.0e12;
	}

	// Reset gap "drag" if zero (or negative).
	//
	if(drag <= PYtolerance) drag = PYtolerance;

	// Only allow zero or positive dashpot values
	//
	if(dashpot < 0.0) dashpot = 0.0;

	// Do not allow zero or negative values for y50 or pult.
	//
	if(pult <= 0.0 || y50 <= 0.0) {
		opserr << "WARNING -- only accepts positive nonzero pult and y50" << endln;
		opserr << "PyLiq1: " << endln;
		opserr << "pult: " << pult << "   y50: " << y50 << endln;
		exit(-1);
	}

	// Initialize variables for Near Field rigid-plastic spring
	//
	if(soilType ==0) {	// This will happen with default constructor
		yref  = 10.0*y50;
		np    = 5.0;
		Elast = 0.35;
		nd    = 1.0;
		TFar_tang   = pult/(8.0*pow(Elast,2.0)*y50);
	}
	else if(soilType ==1) {
		yref  = 10.0*y50;
		np    = 5.0;
		Elast = 0.35;
		nd    = 1.0;
		TFar_tang   = pult/(8.0*pow(Elast,2.0)*y50);
	}
	else if (soilType == 2){
		yref  = 0.5*y50;
		np    = 2.0;
		Elast = 0.2;
		nd    = 1.0;
		// This TFar_tang assumes Elast=0.2, but changes very little for other
		// reasonable Elast values. i.e., API curves are quite linear initially.
		TFar_tang   = 0.542*pult/y50;
	}
	else{
		opserr << "WARNING -- only accepts soilType of 1 or 2" << endln;
		opserr << "PyLiq1: " << endln;
		opserr << "soilType: " << soilType << endln;
		exit(-1);
	}

	// Far Field components: TFar_tang was set under "soil type" statements.
	//
	TFar_p  = 0.0;
	TFar_y  = 0.0;

	// Near Field components
	//
	NFkrig  = 100.0 * (0.5 * pult) / y50;
	TNFpinr = Elast*pult;
	TNFpinl = -TNFpinr;
	TNFyinr = TNFpinr / NFkrig;
	TNFyinl = -TNFyinr;
	TNF_p   = 0.0;
	TNF_y   = 0.0;
	TNF_tang= NFkrig;
	// TNF_tang= ke; // _vag

	// Drag components
	//
	TDrag_pin = 0.0;
	TDrag_yin = 0.0;
	TDrag_p   = 0.0;
	TDrag_y   = 0.0;
	TDrag_tang= nd*(pult*drag-TDrag_p)*pow(y50/2.0,nd)
					*pow(y50/2.0 - TDrag_y + TDrag_yin,-nd-1.0);

	// Closure components
	//
	TClose_yleft = -y50/100.0;
	TClose_yright=  y50/100.0;
	TClose_p     = 0.0; 
	TClose_y     = 0.0;
	TClose_tang  = 1.8*pult*(y50/50.0)*(pow(y50/50.0+ TClose_yright - TClose_y,-2.0)
		+pow(y50/50.0 + TClose_y - TClose_yleft,-2.0));

	// Gap (Drag + Closure in parallel)
	//
	TGap_y   = 0.0;
	TGap_p   = 0.0;
	TGap_tang= TClose_tang + TDrag_tang;

	// Entire element (Far field + Near field + Gap in series)
	//
	Ty       = 0.0;
	Tp       = 0.0;
	Ttangent = pow(1.0/TGap_tang + 1.0/TNF_tang + 1.0/TFar_tang, -1.0);
	TyRate   = 0.0;




	//Vagelis entries

	

	dy = 0.0;
  CLastYieldDir = 0;
  TLastYieldDir = 0;
  signdy = 0;
  Cp = 0.0;
  Cy = 0.0;
  Cyin = 0.0;
  CpinF = 0.0;
  CpinB = 0.0;
  CpinLast = 0.0;
  Ty = 0.0;
  Tyin = 0.0;
  Tp = 0.0;
  TyRate = 0.0;
  // Ttangent = ke;
  TpinF = 0.0;
  TpinB = 0.0;
  TpinUse = 0.0;
  Tpalpha = 0.0;
  dyLast = 0.0;
  signdyLast = 0.0;
  dyELast = 0.0;
  lam = 0;		//lambda, plastic multiplier
  lamLB = 0;
  lamUB = 0;		//lower- and upper-bound guesses for lambda
  ypDot = 0;		//plastic displacement rate
  P1 = 0;			//Force (dP1/dt) in viscoelastic spring component
  P2 = 0;			//Force (dP2/dt) in plastic spring component
  tstep = 0;		//time step, computed from yRate
  yLast = 0;		//y from the timestep before the currently committed timestep, i.e. y(i)=y, y(i-1)=Cy, y(i-2)=yLast
  ypRate = 0;		//plastic deformation rate
  yeRate = 0;		//elastic deformation rate
  TdyP = 0;		//plastic deformation increment
  TdyE = 0;		//elastic deformation increment
  signdy = 0;
  residualLam = 0;	//residual when iterating to determine correct value of plastic multiplier
  Rlam1 = 0;		//residual for lambda using lower- and upper-bound guesses
  Rlam2 = 0;
  sysTimeStep = 0;	//system time
  bumped = 0;
  pP1 = 0;
  dyELastUse = 0;
  signBump = 0.0;
  P1veGuess = 0.0;
  CyeTotal = 0.0;
  TyeTotal = 0.0;
  CtstepLast = 0.0;
  tstep1 = 0.0;
  tstep2 = 0.0;




	// Now get all the committed variables initiated
	//
	this->commitState();

    return 0;
}

/////////////////////////////////////////////////////////////////////
UniaxialMaterial *
PySimple4::getCopy(void)
{
    PySimple4 *theCopy;			// pointer to a PySimple4 class
	theCopy = new PySimple4();	// new instance of this class
	*theCopy= *this;			// theCopy (dereferenced) = this (dereferenced pointer)
	return theCopy;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple4::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(39);
  
  data(0) = this->getTag();
  data(1) = soilType;
  data(2) = pult;
  data(3) = y50;
  data(4) = drag;
  data(5) = dashpot;
  data(6) = yref;
  data(7) = np;
  data(8) = Elast;
  data(9) = nd;
  data(10)= NFkrig;

  data(11) = CNFpinr;
  data(12) = CNFpinl;
  data(13) = CNFyinr;
  data(14) = CNFyinl;
  data(15) = CNF_p;
  data(16) = CNF_y;
  data(17) = CNF_tang;

  data(18) = CDrag_pin;
  data(19) = CDrag_yin;
  data(20) = CDrag_p;
  data(21) = CDrag_y;
  data(22) = CDrag_tang;

  data(23) = CClose_yleft;
  data(24) = CClose_yright;
  data(25) = CClose_p;
  data(26) = CClose_y;
  data(27) = CClose_tang;

  data(28) = CGap_y;
  data(29) = CGap_p;
  data(30) = CGap_tang;

  data(31) = CFar_y;
  data(32) = CFar_p;
  data(33) = CFar_tang;

  data(34) = Cy;
  data(35) = Cp;
  data(36) = Ctangent;
  data(37) = TyRate;

  data(38) = initialTangent;

  data(39) = simple1or3; // Trial input by Vag

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "PySimple4::sendSelf() - failed to send data\n";

  return res;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple4::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(39);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "PySimple4::recvSelf() - failed to receive data\n";
      CNF_tang = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
	soilType = (int)data(1);
	pult     = data(2);
	y50      = data(3);
	drag     = data(4);
	dashpot  = data(5);
	yref     = data(6);
	np       = data(7);
	Elast    = data(8);
	nd       = data(9);
	NFkrig   = data(10);

	CNFpinr  = data(11);
	CNFpinl  = data(12);
	CNFyinr  = data(13);
	CNFyinl  = data(14);
	CNF_p    = data(15);
	CNF_y    = data(16);
	CNF_tang = data(17);

	CDrag_pin = data(18);
	CDrag_yin = data(19);
	CDrag_p   = data(20);
	CDrag_y   = data(21);
	CDrag_tang= data(22);

	CClose_yleft = data(23);
	CClose_yright= data(24);
	CClose_p     = data(25);
	CClose_y     = data(26);
	CClose_tang  = data(27);

	CGap_y    = data(28);
	CGap_p    = data(29);
	CGap_tang = data(30);

	CFar_y    = data(31);
	CFar_p    = data(32);
	CFar_tang = data(33);

	Cy        = data(34);
	Cp        = data(35);
	Ctangent  = data(36);
	TyRate    = data(37);
	
	initialTangent = data(38);

	simple1or3 = data(39); // Trial input by Vag


	// set the trial quantities
	this->revertToLastCommit();
  }

  return res;
}

/////////////////////////////////////////////////////////////////////
void 
PySimple4::Print(OPS_Stream &s, int flag)
{
    s << "PySimple4, tag: " << this->getTag() << endln;
    s << "  soilType: " << soilType << endln;
    s << "  pult: " << pult << endln;
    s << "  y50: " << y50 << endln;
    s << "  drag: " << drag << endln;
	s << "  dashpot: " << dashpot << endln;
}

/////////////////////////////////////////////////////////////////////









//  Below that point complementing functions from PYsimple3

/////////////////////////////////////////////////////////////////////

int
PySimple4::sign(double val)
{
  if (val > 0) return 1;
  if (val < 0) return -1;
  return 0;
}


/////////////////////////////////////////////////////////////////////
//	Residual function for plastic component
double
PySimple4::getResidual(double ke, double Cp, double Tp, double dy, double pu, double C, double Tpin, double dashpot, double tstepCurrent, double dyELast, double CyeTotal, double tstepLast, double Pveguess, double bump)
{
  signdy = sign(dy);
  double signVEguess = sign(Pveguess-Cp);
  if(pu>Tp*signdy)
    {	
      if(tstep != 0.0)
	{		
	  return    (C*ke*(dy-((Tp-Cp+(dashpot*(dyELast/tstepLast)-bump))/(ke+(dashpot/tstepCurrent)))))+(Tp-Cp)+(Tpin-pu*signVEguess)*(log(pu-Cp*signVEguess)-log(pu-Tp*signVEguess));
	}
      else{
	// for static analysis with tstep = 0 the previous form would have blown up...
	return ((Tp-Cp)*(1.0-1.0/C)+((Tpin-pu*signdy)*(log(pu-Tp*signdy)-log(pu-Cp*signdy)))/C - ke*dy);	
      }
    }	
  else
    return 0.0;
}