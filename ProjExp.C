#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TGraphErrors.h"
#include "TPeak.h"
#include "TPeakFitter.h"
#include "TRWPeak.h"
#include "TSinglePeak.h"

//name is the common name between the other files in the .root file from the selector.
/*
ProjExp:

This Function takes a 2D matrix and preforms a projection of the Y-axis whilst 
subtracting the background from the Timerandom matrix (background).

name: The common name across all the histograms, if you use AngularCorrelationSelectorAddback.C then the common name is "addbackAddback"
lowProj: The low edge of projection gate
highProj: The high edge of projection gate

lowBGProj: The low edge of background projection gate
higBGhProj: The high edge of background projection gate

background_normalization: if the background window and projection window are not the same size then a normalization needs to be placed on the background gate

nang: number of angles
*/
void ProjExp(std::string name, float lowProj, float highProj, float lowBGProj = -1., float highBGProj = -1., float lowBG2Proj = -1., float highBG2Proj = -1., float background_normalization = 1, int nang = 50) {
   
   TH1::SetDefaultSumw2();
   TFile *finput = (TFile*)gDirectory->GetFile();
   TFile output("Spectra.root","recreate");
  
   bool subtractBG = (lowBGProj != -1.) || (highBGProj != -1.);
   bool subtractBG2 = (lowBG2Proj != -1.) || (highBG2Proj != -1.);

   std::string name_tmp, filename;
  
   TH2F *matrix = (TH2F*)finput->Get(name.c_str());
   //std::string name_tmpBG = name+"BG";
   //TH2F *matrix2 = (TH2F*)finput->Get(name_tmpBG.c_str());
   
   //TH1D *proj_F =new TH1D("proj_F","proj_F",matrix->GetXaxis()->GetNbins(),matrix->GetXaxis()->GetXmin(),matrix->GetXaxis()->GetXmax());
   //TH1D *proj_FBG =new TH1D("proj_FBG","proj_FBG",matrix->GetXaxis()->GetNbins(),matrix->GetXaxis()->GetXmin(),matrix->GetXaxis()->GetXmax());
   
   TH1D *proj =new TH1D("proj","proj",matrix->GetXaxis()->GetNbins(),matrix->GetXaxis()->GetXmin(),matrix->GetXaxis()->GetXmax()); // defines the histogram: name, name, Number of total bins, min bin value, max bin value.
   TH1D *projBG =new TH1D("projBG","projBG",matrix->GetXaxis()->GetNbins(),matrix->GetXaxis()->GetXmin(),matrix->GetXaxis()->GetXmax());
   TH1D *projBG2 =new TH1D("projBG2","projBG2",matrix->GetXaxis()->GetNbins(),matrix->GetXaxis()->GetXmin(),matrix->GetXaxis()->GetXmax());
   TH1D *projMixed =new TH1D("projMixed","projMixed",matrix->GetXaxis()->GetNbins(),matrix->GetXaxis()->GetXmin(),matrix->GetXaxis()->GetXmax());
   projBG->Sumw2();
   projMixed->Sumw2();
   matrix->Sumw2();
   std::cout << std::endl << "Start to project and export " << nang << " " << name <<"-like matrixes:" << std::endl << std::endl;
   std::cout << "   Low projection limit = " << lowProj << std::endl;
   std::cout << "   High projection limit = " << highProj << std::endl;
   
   //----------------------------------------------------------------------------------
   if(subtractBG) {
     std::cout << "   Low background projection limit = " << lowBGProj << std::endl;
     std::cout << "   High background projection limit = " << highBGProj << std::endl;
   } else {
     std::cout << "   No background projection" << std::endl;   
   }

   if(subtractBG2) {
     std::cout << "   Second Low background projection limit = " << lowBG2Proj << std::endl;
     std::cout << "   Second High background projection limit = " << highBG2Proj << std::endl;
   } else {
     std::cout << "   No background Second projection" << std::endl;   
   }

   //----------------------------------------------------------------------------------
   std::cout << "   Time coincidence normalization = " << background_normalization << std::endl << std::endl;
   //----------------------------------------------------------------------------------
    //matrix->Add(matrix2,-1*background_normalization);
    //matrix->ProjectionX("proj_F",matrix->GetXaxis()->FindBin(lowProj),matrix->GetXaxis()->FindBin(highProj));
    //matrix->ProjectionX("proj_FBG",matrix->GetXaxis()->FindBin(lowProj),matrix->GetXaxis()->FindBin(highProj));
    //proj_F->Sumw2();
    //proj_F->Add(proj_FBG,-1*((highProj-lowProj)/(highBG2Proj-lowBG2Proj)));
   // proj_F->Write(Form("%s",name.c_str()));
    
   //----------------------------------------------------------------------------------
   for(int i=0;i<nang;i++){ // nang are the total number of angles in experiment.
      std::cout << "   Doing matrix " << name <<  i << std::endl;
      //----------------------------------------------------------------------------------
      name_tmp=name+std::to_string(i);
      TH2F *matrix = (TH2F*)finput->Get(name_tmp.c_str()); // Gets the name+i histogram.
      //----------------------------------------------------------------------------------
      name_tmp=name+"BG"+std::to_string(i); 
      TH2F *matrixBG = (TH2F*)finput->Get(name_tmp.c_str());// Gets the name+i background histogram.
      
      matrix->Add(matrixBG,-1*background_normalization); // Some form of background normalization.
      //----------------------------------------------------------------------------------
      
      matrix->ProjectionX("proj",matrix->GetXaxis()->FindBin(lowProj),matrix->GetXaxis()->FindBin(highProj)); //Makes a projection around the lower limit and higher limit set by user.
      
       //----------------------------------------------------------------------------------
       if(subtractBG) {
         matrix->ProjectionX("projBG",matrix->GetXaxis()->FindBin(lowBGProj),matrix->GetXaxis()->FindBin(highBGProj));//Makes a projection of the background around the lower limit and higher limit set by user.
         proj->Sumw2();
         proj->Add(projBG,-1*((highProj-lowProj)/(highBGProj-lowBGProj)));
      }else if(subtractBG2) { 
         matrix->ProjectionX("projBG2",matrix->GetXaxis()->FindBin(lowBG2Proj),matrix->GetXaxis()->FindBin(highBG2Proj));//Makes a projection of the background around the lower limit and higher limit set by user.
         proj->Sumw2();
         proj->Add(projBG2,-1*((highProj-lowProj)/(highBG2Proj-lowBG2Proj)));
      }
      
       //----------------------------------------------------------------------------------
      proj->Write(Form("%s%d",name.c_str(),i));
      
      name_tmp=name+"Mixed"+std::to_string(i);
      TH2F *matrixMixed = (TH2F*)finput->Get(name_tmp.c_str());
      matrixMixed->ProjectionX("projMixed",matrixMixed->GetXaxis()->FindBin(lowProj),matrixMixed->GetXaxis()->FindBin(highProj));
      projMixed->Write(Form("%sMixed%d",name.c_str(),i));
    
      
      
   }
   output.Close();
   std::cout << std::endl << "All done, root spectra saved" << std::endl << std::endl;

}

Double_t DrawingFunction(double *dim, double *par)
{
   Double_t x      = dim[0]; // channel number used for fitting
   Double_t height = par[0]; // height of photopeak
   Double_t c      = par[1]; // Peak Centroid of non skew gaus
   Double_t sigma  = par[2]; // standard deviation of gaussian
   Double_t step   = par[5]; // Size of the step function;
   
   Double_t beta   = par[3]; // Skewness parameter
   Double_t R      = par[4]; // relative height of skewed gaussian

   Double_t gauss      = height * (1.0 - R / 100.0) * TMath::Gaus(x, c, sigma);
   Double_t step_func  = TMath::Abs(step) * height / 100.0 * TMath::Erfc((x - c) / (TMath::Sqrt(2.0) * sigma));
   
   if(beta == 0.0)
      return gauss;
   else
      return (par[6]+x*par[7])+ step_func + gauss + R * height / 100.0 * (TMath::Exp((x - c) / beta)) *
          (TMath::Erfc(((x - c) / (TMath::Sqrt(2.0) * sigma)) + sigma / (TMath::Sqrt(2.0) * beta)));

 
}



void DrawFits(TFitResult* p,double lowlimit, double highlimit,TH1* hist1)
{
   
   int par_t = size(p->Parameters());
   
   double par[par_t];
   int num_f = (par_t-4)/6;      
   for(int i=0;i<par_t;i++) par[i] = p->Parameters()[i]; 
   
   double bk_Cons = par[par_t-4];
   double bk_slope = par[par_t-3];
   
   //Drawing of first function
   TF1 *F = new TF1("F",DrawingFunction,lowlimit,highlimit,8);
   for(int j=0;j<6;j++) F->SetParameter(j,par[j]); //Grabs the first 6 parameters which relate to the first peak
   F->SetParameter(7,0);
   F->SetParameter(8,0);
   F->SetLineColor(kGreen);
   F->Draw();
   
   if(num_f>1){
      TF1 *F1 = new TF1("F1",DrawingFunction,lowlimit,highlimit,8);
      for(int j=0;j<6;j++) F1->SetParameter(j,par[j+6]); //Grabs the first 6 parameters which relate to the first peak
      F1->SetParameter(7,0);
      F1->SetParameter(8,0);
      F1->SetLineColor(kBlue);
      F1->Draw("same");
   
   }
   
   hist1->Draw("same");
   
   
}

 TH1* CreateSliceHist(TH1* hist1,double lowlimit, double highlimit) {

     
    TH1F*hist = new TH1F("hist","Histo",highlimit-lowlimit,lowlimit,highlimit);
    
    double limL = lowlimit - hist1->GetXaxis()->GetXmin();
    double limH = highlimit - hist1->GetXaxis()->GetXmin();
    
    for(int j =0;j<=(limH-limL);++j) {
        
      hist->SetBinContent(hist->FindBin(j+lowlimit),hist1->GetBinContent(j+limL+1)); 
  
    }
    
    return hist;

}

void Testing(float lowProj=686., float highProj=720.){
   std::string name_tmp1 = "addbackAddbackMixed41";
   TFile *finput = (TFile*)gDirectory->GetFile();
   auto Histo = (TH1D*)finput->Get(name_tmp1.c_str()); // Opens each Histogram in the Spectra.root file given.

   auto h = CreateSliceHist(Histo,lowProj-10,highProj+10);
   Histo->SetAxisRange(lowProj-10,highProj+10);
   Histo->SetLineColor(kRed);
   h->Draw();
   Histo->Draw("hist same");


}

void ProjFit(std::string name, float cent = 1., float lowProj=1., float highProj=1.,bool errin = false,int nang = 50, int start = 1,int stop =50) {// Input parameter: "histogram common name", centroidal value, lower projection, high projection, number of angles
   std::vector<std::string> lines;
   std::vector<std::string> lines_area;
   TFile *finput = (TFile*)gDirectory->GetFile();
   std::string name_tmp1;
   std::string name_tmp2;
   bool Badfit = false;
   bool BadfitMixed = false;
   bool determine = false;
   bool filecreated = false;
   
   bool fitplots = false;
   bool viewratio = false;
   
   ofstream exp;

   std::ifstream input_file;
   std::ifstream input_file2;
   
   //|===============================================================================[1]

   //---------------------------------------------------------------------|[Normal AC Histograms]
      //Fit function definitions of the Normal (Non-eventmixed histograms)
      TRWPeak *p1 = new TRWPeak(cent); //Centroidal value defined by user.
      TRWPeak *pM1 = new TRWPeak(1343);
      TRWPeak *pL1 = new TRWPeak(1322.61 );
      TRWPeak *pK1 = new TRWPeak(598.19 );
      
      //TGauss *p1 = new TGauss(cent,-1);
      //TGauss *pM1 = new TGauss(688.89,-1);
      
      // Define fitter
      TPeakFitter *pf = new TPeakFitter(lowProj,highProj);
      // Adding peaks
      pf->AddPeak(p1);
      pf->AddPeak(pM1);
     // pf->AddPeak(pL1);
     // pf->AddPeak(pK1);
       
      // p1->GetFitFunction()->FixParameter(6,5);
     // p1->GetFitFunction()->SetParLimits(1,1341.2,1341.6);
      // p1->GetFitFunction()->SetParameter(0,200);
      // pM1->GetFitFunction()->FixParameter(5,0);
      p1->GetFitFunction()->FixParameter(5,0);    
      pM1->GetFitFunction()->SetParLimits(1,1343.3,1343.6);
      
      pM1->GetFitFunction()->FixParameter(5,0);
      pL1->GetFitFunction()->SetParLimits(1,1322.1,1322.9);
      pL1->GetFitFunction()->FixParameter(5,0);
      pK1->GetFitFunction()->FixParameter(1,598.19);
      pK1->GetFitFunction()->FixParameter(5,0);
      
      // pM1->GetFitFunction()->FixParameter(1,688.89);
      //pM1->GetFitFunction()->FixParameter(5,0);
      //p1->GetFitFunction()->SetParameter(0,500);
      //p1->GetFitFunction()->SetParLimits(5,538.5,540);
      //---------------------------------------------------------------------|[Normal AC Histograms]
   
      //---------------------------------------------------------------------|[Event Mixed AC Histograms]
      //Fit function definitions of the Normal (Event-Mixed histograms)
      TRWPeak *p2 = new TRWPeak(cent); //Centroidal value defined by user.    
      TRWPeak *pM2 = new TRWPeak(1343.49);
      TRWPeak *pM3 = new TRWPeak(1322.61 );
      TRWPeak *pM4 = new TRWPeak(598.19 );
     // TGauss *p2 = new TGauss(cent,-1);
    //  TGauss *pM2 = new TGauss(688.89,-1);
    
      // Define fitter
      TPeakFitter *pf2 = new TPeakFitter(lowProj,highProj);

      // Adding peaks
      pf2->AddPeak(p2);
      pf2->AddPeak(pM2);
     // pf2->AddPeak(pM3);
     // pf2->AddPeak(pM4);
      //pf2->AddPeak(pM5);
      //pf2->AddPeak(pM6);
      
      p2->GetFitFunction()->SetParLimits(1,1341.2,1341.6);
      p2->GetFitFunction()->FixParameter(5,0);    
      pM2->GetFitFunction()->SetParLimits(1,1343.3,1343.6);
      pM2->GetFitFunction()->FixParameter(5,0);
      pM3->GetFitFunction()->SetParLimits(1,1322.1,1322.9);
      pM3->GetFitFunction()->FixParameter(5,0);
      pM4->GetFitFunction()->FixParameter(1,598.19);
      pM4->GetFitFunction()->FixParameter(5,0);
       //p2->GetFitFunction()->FixParameter(5,0);
       //p2->GetFitFunction()->SetParameter(2,7519366.042613);
      // pM2->GetFitFunction()->FixParameter(5,0);
           // pM2->GetFitFunction()->FixParameter(1,688.89);
        //    pM2->GetFitFunction()->SetParLimits(2,1,1.5);
   //---------------------------------------------------------------------|[Event Mixed AC Histograms]

   //|===============================================================================[1]
   

 /*
   This block of code below is here to recognize whether there is an "exp.txt" file 
   that will be updated (from an earlier run of this ProjFit() code).      
 */
 //|===============================================================================[2]
   if (fopen("exp.txt", "r")) { 
   std::cout << "exp.txt exists! Updating of the previous will occur" << std::endl;
   filecreated = true;
   
   //---------------------------------------------------------------------|[Reading in exp.txt lines]
   input_file.open("exp.txt");
   std::string input;
   while (std::getline(input_file, input)) lines.push_back(input);
   //---------------------------------------------------------------------|[Reading in exp.txt lines]
    
   //---------------------------------------------------------------------|[Reading in Area.txt lines]
   std::string input2;
   input_file2.open("Areas.txt");
   while (std::getline(input_file2, input2)) lines_area.push_back(input2);
   //---------------------------------------------------------------------|[Reading in Area.txt lines]
   
   }else{
   exp.open("exp.txt");
   }
 //|===============================================================================[2] 
 
   ofstream fout ("Fitted_Spectra_Details.txt", std::ofstream::app);
   ofstream foutup ("Areas.txt");
   if(filecreated){fout << "|--------------------------------------------------------------------------------------------------------------------------------------------|"<<endl;}
   
  
   if(!filecreated)fout << "Fitted Results Below" << endl;


   //----------------------------------------------------------------------------------
   Double_t w = 1000;
   Double_t h = 700;
   

   if(!filecreated)
   {
         exp << 1 << endl;
         exp << 2 << endl;
   }

   
   auto c1 = new TCanvas("c1", "CANVAS", w, h); 
  // auto c2 = new TCanvas("c2", "CANVAS", w, h); 
   for(int i=start;i<stop;i++){ // nang are the total number of angles in experiment.  
   
         
            std::cout << "\nHistogram: " << name <<  i <<"\n"<< std::endl;
            
         
        //  |======================================================================================|[3][Fitting of the NORMAL AC]
        
            name_tmp1=name+std::to_string(i);
            

            c1->cd();
            auto Histo = (TH1D*)finput->Get(name_tmp1.c_str()); // Opens each Histogram in the Spectra.root file given.
 
          
            name_tmp1 = name_tmp1 + ".root";
            auto h = CreateSliceHist(Histo,lowProj-10,highProj+10);
            //pf->GetBackground()->FixParameter(0,0);
            pf->GetBackground()->FixParameter(1,0);
            Histo->SetTitle(Form("Angle #:%d",i));
            Histo->SetAxisRange(lowProj-1,highProj+1);
            
            if(true){ //to algin the sigma values of all peaks.
               p1->GetFitFunction()->SetParLimits(2,0.5,3);
               pf->Fit(Histo,"NEQM");
               p1->GetFitFunction()->FixParameter(2,p1->GetFitFunction()->GetParameter(2));
               pM1->GetFitFunction()->FixParameter(2,p1->GetFitFunction()->GetParameter(2));
               pL1->GetFitFunction()->FixParameter(2,p1->GetFitFunction()->GetParameter(2)); //p1->GetFitFunction()->GetParameter(2)
               pK1->GetFitFunction()->FixParameter(2,p1->GetFitFunction()->GetParameter(2));
               
            }
            //
            TFitResultPtr Normal1= pf->Fit(Histo,"SEM");//fits the addbackAddback projected histogram
            c1->Update();
          if(!fitplots && !viewratio) c1->WaitPrimitive(); 
             
          if(viewratio)
          {
          //c1->Clear(); // Fit does not draw into correct pad
          auto rp1 = new TRatioPlot(h,"",Normal1.Get());
          rp1->Draw("noconfint");
          c1->Update();
          c1->Modified();
          
          rp1->GetLowerRefYaxis()->SetTitle("Residual");
          rp1->GetLowerRefGraph()->SetMinimum(-5);
          rp1->GetLowerRefGraph()->SetMaximum(5);
          
          c1->Update();
          if(!fitplots) c1->WaitPrimitive(); 
            }
            
            if(fitplots)
            {
          //=================Drawing out the Fits=================|
          auto c2 = new TCanvas("c2", "Fitted Histograms",1000, 600); 
          //h->Draw();
          
          DrawFits(Normal1.Get(),lowProj,highProj,Histo);
          
          c2->Update();
          c2->WaitPrimitive(); 
             c2->Close();
       }        
          //=================Drawing out the Fits=================|
           
           // if(p1->GetChi2()/p1->GetNDF()>30)
            //{
               //c1->WaitPrimitive(); 
              // Badfit = true;
               c1->SaveAs(name_tmp1.c_str());
            //}
   
   //  |======================================================================================|[3] [Fitting of the NORMAL AC]
       
   //  |======================================================================================|[3][Fitting of the EVENT MIXED AC]
       name_tmp2 = name+"Mixed"+std::to_string(i);
            TH1D *Histomixed = (TH1D*)finput->Get(name_tmp2.c_str());
            
            name_tmp2 = name_tmp2 + ".root";

       auto h1_mixed = CreateSliceHist(Histomixed,lowProj-1,highProj+1);           
       Histomixed->SetTitle(Form("Mixed Angle #:%d",i));
       Histomixed->SetAxisRange(lowProj-5,highProj+5);
       pf2->GetBackground()->FixParameter(1,0);
       if(true){ //to algin the sigma values of all peaks.
       	       p2->GetFitFunction()->SetParLimits(2,0.5,3);
               pf2->Fit(Histomixed,"NEQM");
               p2->GetFitFunction()->FixParameter(2,p2->GetFitFunction()->GetParameter(2));
               pM2->GetFitFunction()->FixParameter(2,p2->GetFitFunction()->GetParameter(2));
              pM3->GetFitFunction()->FixParameter(2,p2->GetFitFunction()->GetParameter(2));
               pM4->GetFitFunction()->FixParameter(2,p2->GetFitFunction()->GetParameter(2));
               
            }
            //p2->GetBackgroundFunction()->FixParameter(1,0);
           
            TFitResultPtr Normal2 = pf2->Fit(Histomixed,"SEM"); //fits the addbackAddbackMIXED projected histograms
            c1->Update();
            if(!fitplots && !viewratio) c1->WaitPrimitive(); 
            
       if(viewratio)
       {     
          // TRATIO PLOT 
          c1->Clear(); // Fit does not draw into correct pad
          auto rp2 = new TRatioPlot(h1_mixed,"",Normal2.Get());
          rp2->Draw("noconfint");
          
          rp2->GetLowerRefYaxis()->SetTitle("Residual");
          rp2->GetLowerRefGraph()->SetMinimum(-10);
          rp2->GetLowerRefGraph()->SetMaximum(10);
          c1->Update();
          if(!fitplots) c1->WaitPrimitive(); 
       }
            
            if(fitplots)
            {
          //=================Drawing out the Fits=================|
          auto c3 = new TCanvas("c3", "Fitted Histograms",1000, 600); 
          //h->Draw();
          
          DrawFits(Normal2.Get(),lowProj,highProj,h1_mixed);
          
          c3->Update();
          c3->WaitPrimitive(); 
             c3->Close();
          //=================Drawing out the Fits=================|
            }
            
            if(filecreated)  c1->WaitPrimitive();
      
           // if(p2->GetChi2()/p2->GetNDF()>30)
           // {
      
               //BadfitMixed = true;
               c1->SaveAs(name_tmp2.c_str());
           // }
         //  |======================================================================================|[3][Fitting of the EVENT MIXED AC]   

         
      //this will be for drawing the secondary or more peaks that are used, for better visualization of whats going on. 
   



         //  |======================================================================================|[4][WRITING TO FILES]  
            float sqrtred1 = TMath::Sqrt(p1->GetChi2()/p1->GetNDF());
            float sqrtred2 = TMath::Sqrt(p2->GetChi2()/p2->GetNDF());
            if(!filecreated)
            {
               fout << "\t\taddbackAddback_Projected:" << i <<"\t\t" << "addbackAddbackMixed_Projected:" << i << endl;
               if(Badfit){fout <<"---------BAD FIT---------" << endl;}
               if(BadfitMixed){fout <<"---------BAD MIXED FIT---------" << endl;}
               fout << "\t| Normal Centroid |" << " +/- " << "| Normal Centroid Error\t|"  << "\t|" << " Mixed Centroid |" << " +/- " << "| Mixed Centroid Error|"  << "\t|" << endl;
               fout << "\t| " <<p1->Centroid() <<  "\t\t\t" <<p1->CentroidErr() << "\t|" << "\t| " <<p2->Centroid() <<  "\t\t\t" <<p2->CentroidErr() << "\t|" << endl;
               fout << "\t| Area:\t"<<p1->Area()<<" +/- "<<p1->AreaErr() <<"\t\t\t\t"<< "   Mixed Area:\t"<<p2->Area()<<" +/- "<<p2->AreaErr() << "\t\t|" << endl;
               
               //Stores results of the areas.
               //important for the normalization of the data points the code does.
               foutup << i+3 << "\t" <<p1->Area()<<"\t"<<p1->AreaErr()<<"\t"<< p2->Area() <<"\t"<<p2->AreaErr()<< endl;
            }

            if(filecreated)
            {
               fout << "\t\taddbackAddback_Projected:" << i <<"\t\t" << "addbackAddbackMixed_Projected:" << i << endl;
               if(Badfit){fout <<"---------BAD FIT---------" << endl;}
               if(BadfitMixed){fout <<"---------BAD MIXED FIT---------" << endl;}
                fout << "\t| Normal Centroid |" << " +/- " << "| Normal Centroid Error\t|"  << "\t|" << " Mixed Centroid |" << " +/- " << "| Mixed Centroid Error|"  << "\t|" << endl;
               fout << "\t| " <<p1->Centroid() <<  "\t\t" <<p1->CentroidErr() << "\t|" << "\t| " <<p2->Centroid() <<  "\t\t" <<p2->CentroidErr() << "\t|" << endl;
               fout << "\t| Area:\t"<<p1->Area()<<" +/- "<<p1->AreaErr() <<"\t\t\t\t"<< "   Mixed Area:\t"<<p2->Area()<<" +/- "<<p2->AreaErr() << "\t\t|" << endl;
            }   
               
   //  |=================================================================================|[5][Updating Results]
            double value_area = p1->Area()/p2->Area();      
            if(filecreated)
            {
           
               for (auto& line : lines)
               {
                  
                   if((i+3) == atoi((line.substr(0,2).c_str())))
                  {
                     double term1 = value_area*TMath::Sqrt( TMath::Power( (p2->AreaErr()/p2->Area()) ,2) + TMath::Power( (p1->AreaErr()/p1->Area()) ,2) );
                     line = to_string(i+3) + "\t" + to_string(value_area) + "\t" + to_string(term1);

                  }
               }  
            }else{
               
               exp << i+3 <<"\t"<< value_area <<"\t"<< value_area*TMath::Sqrt( TMath::Power( (p2->AreaErr()/p2->Area()) ,2) + TMath::Power( (p1->AreaErr()/p1->Area()) ,2) ) << endl;
               
            }
            
            if(filecreated)
            {
           
               for (auto& line1 : lines_area)
               {
                  
                   if((i+3) == atoi((line1.substr(0,2).c_str())))
                  {
           
                     line1 = to_string(i+3) + "\t" + p1->Area() + "\t" + p1->AreaErr() + "\t" + p2->Area() + "\t" + p2->AreaErr();

                  }
               }  
            }
   //  |=================================================================================|[5][Updating Results]    
            Badfit = false;
            BadfitMixed = false;
   }
   //  |======================================================================================|[4][WRITING TO FILES]
   
   //  |==============================================|[6][Writing Updated FILES]
         //The writing of the updated file
         
            if(filecreated)
            {
              input_file.close();

               std::ofstream output_file("exp.txt");
               for (auto const& line : lines)
               output_file << line << '\n';
            
               output_file.close();          
            }
      
            if(filecreated)
            {
              input_file.close();
          
               std::ofstream output_file1("Areas.txt");
               for (auto const& line1 : lines_area)
               output_file1 << line1 << '\n';
            
               output_file1.close();          
            } 
    //  |==============================================|[6][Writing Updated FILES]
    
           
   fout.close();


}

