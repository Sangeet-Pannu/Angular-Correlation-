#include "Riostream.h"
#include <iostream>
#include <assert.h>
#include <cmath>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "RooMath.h"
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
FitCode.C Will be used to evaluate the Attenuation Factors for 
Each of the Detector Pair Configurations.
*/

float Wtheo[29];
float y[29];
float yr[29];

int J1=2; // change this for whatever spin it is.
int J2=2;
int J3=0;

float delta1 = 3.7;
float delta2 = 0;

float A2;
float A4;

float R2LLJ2J1 = 0;
float R2LMJ2J1 = 0;
float R2MMJ2J1 = 0;
float R2LLJ2J3 = 0;
float R2LMJ2J3 = 0;
float R2MMJ2J3 = 0;
float R4LLJ2J1 = 0;
float R4LMJ2J1 = 0;
float R4MMJ2J1 = 0;
float R4LLJ2J3 = 0;
float R4LMJ2J3 = 0;
float R4MMJ2J3 = 0;

//=================================================================+|
void Coefficient_Set(){
if (J1 == 1 && J2 == 2 && J3 == 0) {
        R2LLJ2J1 = 0.4183;
        R2LMJ2J1 = 0.9354;
        R2MMJ2J1 = -0.2988;
        R2LLJ2J3 = -0.5976;
        R2LMJ2J3 = 0.0;
        R2MMJ2J3 = 0.0;
        R4LLJ2J1 = 0.0;
        R4LMJ2J1 = 0.0;
        R4MMJ2J1 = 0.7127;
        R4LLJ2J3 = -1.0690;
        R4LMJ2J3 = 0.0;
        R4MMJ2J3 = 0.0;
    }else if (J1 == 2 && J2 == 2 && J3 == 0) {
        R2LLJ2J1 = -0.4183;
        R2LMJ2J1 = 0.6124;
        R2MMJ2J1 = 0.1281;
        R2LLJ2J3 = -0.5976;
        R2LMJ2J3 = 0.0;
        R2MMJ2J3 = 0.0;
        R4LLJ2J1 = 0.0;
        R4LMJ2J1 = 0.0;
        R4MMJ2J1 = -0.3064;
        R4LLJ2J3 = -1.0690;
        R4LMJ2J3 = 0.0;
        R4MMJ2J3 = 0.0;
        }else if (J1 == 3 && J2 == 2 && J3 == 0) {
            R2LLJ2J1 = 0.1195;
            R2LMJ2J1 = -0.6547;
            R2MMJ2J1 = 0.3415;
            R2LLJ2J3 = -0.5976;
            R2LMJ2J3 = 0.0;
            R2MMJ2J3 = 0.0;
            R4LLJ2J1 = 0.0;
            R4LMJ2J1 = 0.0;
            R4MMJ2J1 = 0.0764;
            R4LLJ2J3 = -1.0690;
            R4LMJ2J3 = 0.0;
            R4MMJ2J3 = 0.0;
            }else if (J1 == 4 && J2 == 2 && J3 == 0) {
                R2LLJ2J1 = -0.1707;
                R2LMJ2J1 = -0.5051;
                R2MMJ2J1 = 0.4482;
                R2LLJ2J3 = -0.5976;
                R2LMJ2J3 = 0.0;
                R2MMJ2J3 = 0.0;
                R4LLJ2J1 = -0.0085;
                R4LMJ2J1 = 0.0627;
                R4MMJ2J1 = -0.0297;
                R4LLJ2J3 = -1.0690;
                R4LMJ2J3 = 0.0;
                R4MMJ2J3 = 0.0;
            }else if (J1 == 0 && J2 == 2 && J3 == 0){
                R2LLJ2J1 = -0.5976;
                R2LMJ2J1 = 0.0;
                R2MMJ2J1 = 0.0;
                R2LLJ2J3 = -0.5976;
                R2LMJ2J3 = 0.0;
                R2MMJ2J3 = 0.0;
                R4LLJ2J1 = -1.0690;
                R4LMJ2J1 = -0.0;
                R4MMJ2J1 = -0.0;
                R4LLJ2J3 = -1.0690;
                R4LMJ2J3 = 0.0;
                R4MMJ2J3 = 0.0;
            
            }else{
                cout << "Some Thing Went Wrong!" << endl;
            }                
}
//=================================================================+|

/*
    exp.txt is the experimental coinicidence fits on the respective cascade
    one is looking at, with the implemented EVENT MIXING.
*/
//angs "String" changes between "Fangs" (FIPPS ANGLE ID) and "Iangs" (IFIN ANGLE ID).

TGraph* CreateHisto(std::string exp = "./Normexp.txt",std::string angs = "/home/Sangeet/Master_project/100Ru-Data/Angular_Correlations/Fangs.txt"){
   TGraph2D graphexp(exp.c_str());
   TGraph2D graphangs(angs.c_str());

   auto binEdgesexp = graphangs.GetX();
   auto valuesexp = graphexp.GetY();
   auto errorexp = graphexp.GetZ();
   auto histexp = new TGraphErrors(graphexp.GetN()-1,binEdgesexp,valuesexp,nullptr,errorexp);

   for(int j =0; j<graphexp.GetN();j++){
      y[j] = valuesexp[j];
      yr[j] = errorexp[j];
   }

 
   return histexp;

}

// implementing a simple normalization.
double AngCorFunction_M1(double *x,double* par)
{
   double rad = x[0];   
   double Pl2 = (1.0/2.0)*(3.*pow(rad,2)-1.);
   double Pl4 = (1.0/8.0)*(35.*pow(rad,4)-30.*pow(rad,2)+3.);

   return (par[4])*(1 + par[2]*par[0]*Pl2+ par[3]*par[1]*Pl4);

}

//implementing Sambu's method
double AngCorFunction_M2(double *x,double* par)
{   
    double rad = x[0];   
    double Pl2 = (1.0/2.0)*(3*pow(rad,2)-1);
    double Pl4 = (1.0/8)*(35*pow(rad,4)-30*pow(rad,2)+3);

    float a2 = ((R2LLJ2J1 + 2 * delta1 * R2LMJ2J1 + pow(delta1, 2) * R2MMJ2J1) * (R2LLJ2J3 + 2 * delta2 * R2LMJ2J3 + pow(delta2, 2) * R2MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));
    float a4 = ((R4LLJ2J1 + 2 * delta1 * R4LMJ2J1 + pow(delta1, 2) * R4MMJ2J1) * (R4LLJ2J3 + 2 * delta2 * R4LMJ2J3 + pow(delta2, 2) * R4MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));

    double theo_W = (1 + par[0]*a2*Pl2+ par[1]*a4*Pl4);
    
    return theo_W;

}

//All cascades other than 0-2-0
double AngCorFunction_M3(double *x,double* par)
{
    

    double rad = x[0];   
    double Pl2 = (1.0/2.0)*(3*pow(rad,2)-1);
    double Pl4 = (1.0/8)*(35*pow(rad,4)-30*pow(rad,2)+3);

    //?-2-0 Any with mixing ratios
  
    float a2 = ((R2LLJ2J1 + 2 * delta1 * R2LMJ2J1 + pow(delta1, 2) * R2MMJ2J1) * (R2LLJ2J3 + 2 * delta2 * R2LMJ2J3 + pow(delta2, 2) * R2MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));
    float a4 = ((R4LLJ2J1 + 2 * delta1 * R4LMJ2J1 + pow(delta1, 2) * R4MMJ2J1) * (R4LLJ2J3 + 2 * delta2 * R4LMJ2J3 + pow(delta2, 2) * R4MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));

    double theo_W = (1 + a2*Pl2+a4*Pl4);
    
    return theo_W;

}



void FitCode(){
    //Uses the exp file to create a graph which will be used to fit the angular correlation results on.

    Coefficient_Set();

    auto c1 = new TCanvas("c1", "Angular Distribution With Using SAMBU");
    c1->Divide(2,1);
    c1->cd(1);

    TGraph* fhist = CreateHisto();
    TGraph* fhist1 = CreateHisto();
    
    fhist->SetMarkerColor(kBlue);
    fhist->SetMarkerStyle(kFullCircle);

    fhist1->SetMarkerColor(kGreen);
    fhist1->SetMarkerStyle(kFullCircle);

    
    TF1 *fitfunc = new TF1("fitfunc",AngCorFunction_M2,-1,1,2);

    //Below for when we are looking at cascades other than 0-2-0.
    //TF1 *fitfunc = new TF1("fitfunc",AngCorFunction_M3,-1,1,2);

    fitfunc->SetParName(0,"Q2");
    fitfunc->SetParName(1,"Q4");

    fitfunc->SetParameter(0,0.8); // Background slope
    fitfunc->SetParameter(1,0.8);
    fitfunc->SetParLimits(0,0.1,1);
    fitfunc->SetParLimits(1,0.1,1);

    fhist->Fit("fitfunc","REM");
    cout << "|**----------FIT FINISHED for SAMBU METHOD----------**|"<< endl;
    cout << "Ch1_2/NDF= " << fitfunc->GetChisquare()/fitfunc->GetNDF() << endl;

    cout << "\nQ2 = " << fitfunc->GetParameter(0) <<" +/- "<< fitfunc->GetParError(0) << endl;
    cout << "Q4 = " << fitfunc->GetParameter(1) <<" +/- "<< fitfunc->GetParError(1) << endl;
    fhist->SetTitle("Angular Distribution With Using SAMBU");
    fhist->Draw("AP");
    

    c1->cd(2);
    cout << "\n" << endl;

    TF1 *fitfunc1 = new TF1("fitfunc1",AngCorFunction_M1,-1,1,5);

    fitfunc1->SetParName(0,"a2");
    fitfunc1->SetParName(1,"a4");
    fitfunc1->SetParName(2,"Q2");
    fitfunc1->SetParName(3,"Q4");
    fitfunc1->SetParName(4,"Norm C");   

    fitfunc1->FixParameter(0,0.1020408163265306); //A2 Value Gotten From GRIFFIN TOOLS Angular Distribution Calculator
    fitfunc1->FixParameter(1,0.009070294784580489); //A4 Value Gotten From GRIFFIN TOOLS Angular Distribution Calculator
    
    fitfunc1->SetParameter(2,0.9); 
    fitfunc1->SetParameter(3,0.6);
    fitfunc1->SetParLimits(3,0.1,1);
    fitfunc1->SetParLimits(2,0.1,1);
    
    fitfunc1->SetParameter(4,1);

    fhist1->Fit("fitfunc1","REM");
    cout << "|**----------FIT FINISHED for NORM C METHOD----------**|"<< endl;
    cout << "Ch1_2/NDF= " << fitfunc1->GetChisquare()/fitfunc1->GetNDF() << endl;

    cout << "\nQ2 = " << fitfunc1->GetParameter(2) <<" +/- "<< fitfunc1->GetParError(2) << endl;
    cout << "Q4 = " << fitfunc1->GetParameter(3) <<" +/- "<< fitfunc1->GetParError(3) << endl;

    cout << "\nNORM C = " << fitfunc1->GetParameter(4) <<" +/- "<< fitfunc1->GetParError(4) << endl;
    fhist1->SetTitle("Angular Distribution of NORM C");
    fhist1->Draw("AP");
    
   
 
    

}
