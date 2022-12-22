#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"


float y[55];
float yr[55];

// definition of shared parameter
// background function 
int iparB[4] = { 0,      // A2
                 1,      // A4
                 2,	//Respective Q2 factor
                 3,	//Respective Q4 factor
};

// signal + background function 
int iparSB[4] = { 0, // A2
                  1, // A4
                  4, //Respective Q2 factor
                  5, //Respective Q4 factor

};

int iparSBL[4] = { 0, // A2
                  1, // A4
                  6, //Respective Q2 factor
                  7, //Respective Q4 factor

};

struct GlobalChi2 { 
   GlobalChi2(  vector< ROOT::Fit::Chi2Function * > & fVec)  
             {
	            for(unsigned int i = 0; i < fVec.size(); i++){
		           fChi2_Vec.push_back( fVec[i] );
				}
			 }

   // parameter vector is first background (in common 1 and 2) 
   // and then is signal (only in 2)
   double operator() (const double *par) const {
      vector< vector<double> > pVec;
      vector< double > dummyVec;
      
      for (int i = 0; i < 4; ++i) dummyVec.push_back(par[iparB[i] ]);
      pVec.push_back(dummyVec);

      dummyVec.clear();
      for (int i = 0; i < 4; ++i) dummyVec.push_back(par[iparSB[i] ]);
      pVec.push_back(dummyVec);
      
      dummyVec.clear();
      for (int i = 0; i < 4; ++i) dummyVec.push_back(par[iparSBL[i] ]);
      pVec.push_back(dummyVec);
      
      double val = 0;
      for( size_t i = 0; i < fChi2_Vec.size(); i++ ){
        val += (*fChi2_Vec[i])(&(pVec[i][0]));
      }

      return val;
   } 

   vector< const ROOT::Math::IMultiGenFunction * > fChi2_Vec;
};

double AngCorFunction_M1(double *x,double* par)
{
	double rad = x[0];   
	double Pl2 = (1.0/2.0)*(3*pow(rad,2)-1);
	double Pl4 = (1.0/8)*(35*pow(rad,4)-30*pow(rad,2)+3);

	return (1 + par[2]*par[0]*Pl2+par[1]*par[3]*Pl4);

}


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

void combinedFit() { 
  //Remember to change the 
  TGraph* fhist = CreateHisto(); //FIPPS-FIPPS Distribution
  TGraph* fhist1 = CreateHisto("./Normexp1.txt","/home/Sangeet/Master_project/100Ru-Data/Angular_Correlations/Iangs.txt"); //Ifin-Ifin Distribution
  TGraph* fhist2 = CreateHisto("./Normexp2.txt","/home/Sangeet/Master_project/100Ru-Data/Angular_Correlations/FIangs.txt"); //Fipps-Ifin Distribution

  //|----------------FIPPS-FIPPS---------------------------------------|
  TF1 *fitfunc1 = new TF1("fitfunc1",AngCorFunction_M1,-1,1,4);
  fitfunc1->SetParName(0,"a2");
  fitfunc1->SetParName(1,"a4");
  fitfunc1->SetParName(2,"Q2");
  fitfunc1->SetParName(3,"Q4"); 

  fitfunc1->SetParameter(0,-0.07142857142857141); 
  fitfunc1->SetParameter(1,0);
  fitfunc1->FixParameter(2,0.8833); //Q2 for FIPPS-FIPPS
  fitfunc1->FixParameter(3,0.6834); //Q2 for FIPPS-FIPPS
  
  //|----------------FIPPS-FIPPS---------------------------------------|

  //|----------------IFIN-IFIN-----------------------------------------|
  TF1 *fitfunc = new TF1("fitfunc",AngCorFunction_M1,-1,1,4);
  fitfunc->SetParName(0,"a2");
  fitfunc->SetParName(1,"a4");
  fitfunc->SetParName(2,"Q2");
  fitfunc->SetParName(3,"Q4"); 

  fitfunc->SetParameter(0,-0.07142857142857141); 
  fitfunc->SetParameter(1,0); 
  fitfunc->FixParameter(2,0.9390); //Q2 for IFIN-IFINS
  fitfunc->FixParameter(3,0.7359); //Q4 for IFIN-IFINS
  
  //|----------------IFIN-IFIN-----------------------------------------|
  
  //|----------------FIPP-IFIN-----------------------------------------|
  TF1 *fitfunc2 = new TF1("fitfunc2",AngCorFunction_M1,-1,1,4);
  fitfunc->SetParName(0,"a2");
  fitfunc->SetParName(1,"a4");
  fitfunc->SetParName(2,"Q2");
  fitfunc->SetParName(3,"Q4"); 

  fitfunc->SetParameter(0,-0.07142857142857141); 
  fitfunc->SetParameter(1,0); 
  fitfunc->FixParameter(2,0.9390); //Q2 for FIPPS-IFINS
  fitfunc->FixParameter(3,0.7359); //Q4 for FIPPS-IFINS
  
  //|----------------FIPP-IFIN-----------------------------------------|


  ROOT::Math::WrappedMultiTF1 wfB(*fitfunc1,1); //Fipps-Fipps
  ROOT::Math::WrappedMultiTF1 wfSB(*fitfunc,1); //Ifin-Ifin
  ROOT::Math::WrappedMultiTF1 wfSBL(*fitfunc2,1); //Fipps-Ifin
  	
  ROOT::Fit::DataOptions opt; 
  
  
  // set the data range
  //-------|FIPPS-FIPPS|----------------|
  ROOT::Fit::DataRange rangeB; 
  rangeB.SetRange(-1,1);
  ROOT::Fit::BinData dataB(opt,rangeB); 
  ROOT::Fit::FillData(dataB, fhist);
  //-------|FIPPS-FIPPS|----------------|
  
  //-------|IFIN-IFIN|----------------|
  ROOT::Fit::DataRange rangeSB; 
  rangeSB.SetRange(-1,1);
  ROOT::Fit::BinData dataSB(opt,rangeSB); 
  ROOT::Fit::FillData(dataSB, fhist1);
  //-------|IFIN-IFIN|----------------|
  
  //-------|FIPPS-IFIN|----------------|
  ROOT::Fit::DataRange rangeSBL; 
  rangeSBL.SetRange(-1,1);
  ROOT::Fit::BinData dataSBL(opt,rangeSBL); 
  ROOT::Fit::FillData(dataSBL, fhist2);
  //-------|FIPPS-IFIN|----------------|
  
  vector< ROOT::Fit::Chi2Function * > chi2_vec;
  chi2_vec.push_back( new ROOT::Fit::Chi2Function(dataB, wfB) );
  chi2_vec.push_back( new ROOT::Fit::Chi2Function(dataSB, wfSB) );
  chi2_vec.push_back( new ROOT::Fit::Chi2Function(dataSBL, wfSBL) );
  
  GlobalChi2 globalChi2(chi2_vec);
  
  ROOT::Fit::Fitter fitter;
  
  const int Npar = 8; 
  double par0[Npar] = { 0.5,0.5,0.8833,0.6834,0.9390,0.7359,0.9390,0.7359}; //{ A2 , A4 , Q2_FF , Q4_FF , Q2_II , Q4_II , Q2_FI , Q4_FI }
  fitter.Config().SetParamsSettings(8,par0);
  
  fitter.Config().ParSettings(2).Fix();
  fitter.Config().ParSettings(3).Fix();
  
  fitter.Config().ParSettings(4).Fix();
  fitter.Config().ParSettings(5).Fix();
  
  fitter.Config().ParSettings(6).Fix();
  fitter.Config().ParSettings(7).Fix();
  
  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");
  
  fitter.FitFCN(8,globalChi2,0,dataB.Size()+dataSB.Size()+dataSBL.Size(),true);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);
 
  fhist->SetMarkerColor(kBlue);
  fhist->SetMarkerStyle(kFullCircle);

  fhist1->SetMarkerColor(kGreen);
  fhist1->SetMarkerStyle(kFullCircle);
  
  fhist2->SetMarkerColor(kBlack);
  fhist2->SetMarkerStyle(kFullCircle);
    
  TCanvas * c1 = new TCanvas("Simfit","Simultaneous fit of Three histograms",10,10,1030,1030);
  c1->Divide(1,3);
  gStyle->SetOptFit(1111);
  //-------|FIPPS-FIPPS|----------------|
  c1->cd(1);
  fitfunc1->SetFitResult( result, iparB);
  fitfunc1->SetRange(rangeB().first, rangeB().second);   
  fitfunc1->SetLineColor(kBlue);
  fhist->GetListOfFunctions()->Add(fitfunc1);
  fhist->Draw("AP"); 
  //-------|FIPPS-FIPPS|----------------|

  //-------|IFIN-IFIN|----------------|
  c1->cd(2);
  fitfunc->SetFitResult( result, iparSB);
  fitfunc->SetRange(rangeSB().first, rangeSB().second);   
  fitfunc->SetLineColor(kGreen);
  fhist1->GetListOfFunctions()->Add(fitfunc);
  fhist1->Draw("AP");
  //-------|IFIN-IFIN|----------------|
  
  //-------|FIPPS-IFIN|----------------|
  c1->cd(3);
  fitfunc2->SetFitResult( result, iparSBL);
  fitfunc2->SetRange(rangeSBL().first, rangeSBL().second);   
  fitfunc2->SetLineColor(kBlack);
  fhist2->GetListOfFunctions()->Add(fitfunc2);
  fhist2->Draw("AP");
  //-------|FIPPS-IFIN|----------------|
  }


