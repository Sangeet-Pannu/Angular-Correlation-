#include "Riostream.h"
#include <iostream>
#include <assert.h>
#include <cmath>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "RooMath.h"


/*
SOFTWARE: 112 COULEX DATA FITTING FUNCTION.
Author: Sangeet-Pal Pannu

FUNCTIONALITY: This code implements two skewed gaussian fitting,
in which each peak is fit using two skewed gaussians in order to
yield a high degree of variablity when fitting the skewed tail feature.

- Fits two peaks together on top a quadratic background.
*/

//Linear background
//par[0] = Slope
//par[1] = Constant

 // additional peaks to be added to system besides groundstate.




//W(theta) = 1 + a_2*P_2(cos(theta)) + a_4*P_4(cos(theta))

double AngCorFunction(double *x,double* par)
{
	double rad = x[0];   
	double Pl2 = (1.0/2.0)*(3*pow(rad,2)-1);
	double Pl4 = (1.0/8)*(35*pow(rad,4)-30*pow(rad,2)+3);	
	return (1 + par[2]*par[0]*Pl2+par[1]*par[3]*Pl4)+par[4];

}




/*
TH1* OpenFile() Function's functionality is opening ascii files and writing them to
.root files via TFile *file = new TFile().

To change ascii file, simply remove pre-existing file in "in.open()" and place the specific ascii 
file to be read in the brackets.

The naming of the root file can be changed by altering the content of line 97.
*/

// WANT TO CREATE A GRAPH OF SCATTER POINTS NOT A 1D HISTOGRAM.
   
#include <TROOT.h>
#include <TFile.h>
#include <string>
#include <stdio.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include "TGraphErrors.h"
#include "TFrame.h"
#include "TF1.h"
#include "TTree.h"
/*
    data.root is the output file for which the histograms will be written to
    for the next setup of the analysis of angular correlations.
    
    exp.txt is the experimental coinicidence fits on the respective cascade
    one is looking at, with the implemented EVENT MIXING.
    
    comb??? is the simulation files for the respective cascade.
*/
//angs "String" changes between "Fangs" (FIPPS ANGLE ID) and "Iangs" (IFIN ANGLE ID).

TGraph* CreateHisto(
   std::string data = "./data.root",
   std::string exp = "./exp.txt",std::string angs = "./Fang.txt"){

   TGraph2D graphexp(exp.c_str());
   TGraph2D graphangs(angs.c_str());

   auto binEdgesexp = graphangs.GetX();
   
   auto valuesexp = graphexp.GetY();
   
   auto errorexp = graphexp.GetZ();

   auto histexp = new TGraphErrors(graphexp.GetN()-1,binEdgesexp,valuesexp,nullptr,errorexp);

   return histexp;

}




void FitCode(){
    //Uses the exp file to create a graph which will be used to fit the angular correlation results on.
    //auto c1 = new TCanvas("c1", "Fit spectrum with Residual plot");
    TGraph* fhist = CreateHisto();
    fhist->SetMarkerColor(kBlue);
    fhist->SetMarkerStyle(kFullCircle);
    

    TF1 *fitfunc = new TF1("fitfunc",AngCorFunction,-1,1,5); 

    fitfunc->SetParName(2,"Q2");
    fitfunc->SetParName(3,"Q4");
   

    fitfunc->FixParameter(0,0.3571428571428571); // Background slope
    fitfunc->FixParameter(1,1.1428571428571415); // Background offset 
    fitfunc->FixParameter(2,0.8); // Background slope
    fitfunc->FixParameter(3,0.6);
   // fitfunc->SetParLimits(0,0,10);
    // fitfunc->SetParLimits(1,0,10);
    fhist->Fit("fitfunc","REM");
    fhist->Draw("AP");
    
   
 
    

}
