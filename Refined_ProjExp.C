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
  
   TH1D *proj =new TH1D("proj","proj",matrix->GetXaxis()->GetNbins(),matrix->GetXaxis()->GetXmin(),matrix->GetXaxis()->GetXmax()); // defines the histogram: name, name, Number of total bins, min bin value, max bin value.
  
   TH1D *projBG =new TH1D("projBG","projBG",matrix->GetXaxis()->GetNbins(),matrix->GetXaxis()->GetXmin(),matrix->GetXaxis()->GetXmax());

   TH1D *projBG2 =new TH1D("projBG2","projBG2",matrix->GetXaxis()->GetNbins(),matrix->GetXaxis()->GetXmin(),matrix->GetXaxis()->GetXmax());

   TH1D *projMixed =new TH1D("projMixed","projMixed",matrix->GetXaxis()->GetNbins(),matrix->GetXaxis()->GetXmin(),matrix->GetXaxis()->GetXmax());
	projBG->Sumw2();
	projMixed->Sumw2();
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
   
   //----------------------------------------------------------------------------------
   for(int i=0;i<nang;i++){ // nang are the total number of angles in experiment.
      std::cout << "   Doing matrix " << name <<  i << std::endl;
      //----------------------------------------------------------------------------------
      name_tmp=name+std::to_string(i);
      TH2F *matrix = (TH2F*)finput->Get(name_tmp.c_str()); // Gets the name+i histogram.
      //----------------------------------------------------------------------------------
      name_tmp=name+"BG"+std::to_string(i); 
      TH2F *matrixBG = (TH2F*)finput->Get(name_tmp.c_str());// Gets the name+i background histogram.
      matrix->Sumw2();
      matrix->Add(matrixBG,-1*background_normalization); // Some form of background normalization.
      //----------------------------------------------------------------------------------
      
      matrix->ProjectionX("proj",matrix->GetXaxis()->FindBin(lowProj),matrix->GetXaxis()->FindBin(highProj)); //Makes a projection around the lower limit and higher limit set by user.
      
       //----------------------------------------------------------------------------------
      if(subtractBG) {
         matrix->ProjectionX("projBG",matrix->GetXaxis()->FindBin(lowBGProj),matrix->GetXaxis()->FindBin(highBGProj));//Makes a projection of the background around the lower limit and higher limit set by user.
         proj->Sumw2();
         if(subtractBG2){
            proj->Add(projBG,-(1/2)*((highProj-lowProj)/(highBGProj-lowBGProj)));   
         }else{
            proj->Add(projBG,-1*((highProj-lowProj)/(highBGProj-lowBGProj)));   
         }
         
      }
      if(subtractBG2) {
         matrix->ProjectionX("projBG2",matrix->GetXaxis()->FindBin(lowBG2Proj),matrix->GetXaxis()->FindBin(highBG2Proj));//Makes a projection of the background around the lower limit and higher limit set by user.
         proj->Sumw2();
         proj->Add(projBG2,-(1/2)*((highProj-lowProj)/(highBG2Proj-lowBG2Proj)));
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

   ofstream exp;

   std::ifstream input_file;
   std::ifstream input_file2;
   
   //|===============================================================================[1]

	//---------------------------------------------------------------------|[Normal AC Histograms]
   	//Fit function definitions of the Normal (Non-eventmixed histograms)
	   TRWPeak *p1 = new TRWPeak(cent); //Centroidal value defined by user.
	   //TRWPeak *pM1 = new TRWPeak(677);
	   //TRWPeak *pL1 = new TRWPeak(699);
	   //TRWPeak *pK1 = new TRWPeak(702);
	   
	   // Define fitter
	   TPeakFitter *pf = new TPeakFitter(lowProj,highProj);
	   // Adding peaks
	   pf->AddPeak(p1);
	   //pf->AddPeak(pM1);
	   //pf->AddPeak(pL1);
	   //pf->AddPeak(pK1);
	   
   	//---------------------------------------------------------------------|[Normal AC Histograms]
   
   	//---------------------------------------------------------------------|[Event Mixed AC Histograms]
   	//Fit function definitions of the Normal (Event-Mixed histograms)
	   TRWPeak *p2 = new TRWPeak(cent); //Centroidal value defined by user.
	   //TRWPeak *pM2 = new TRWPeak(677); //Centroidal value defined by user.
	   //TRWPeak *pM3 = new TRWPeak(699); //Centroidal value defined by user.
	   //TRWPeak *pM4 = new TRWPeak(702); //Centroidal value defined by user.
	   //TRWPeak *pM5 = new TRWPeak(835.904); //Centroidal value defined by user.
	   //TRWPeak *pM6 = new TRWPeak(711.356); //Centroidal value defined by user.
	   
	   // Define fitter
	   TPeakFitter *pf2 = new TPeakFitter(lowProj,highProj);

	   // Adding peaks
	   pf2->AddPeak(p2);
	   //pf2->AddPeak(pM2);
	   //pf2->AddPeak(pM3);
	   //pf2->AddPeak(pM4);
	   //pf2->AddPeak(pM5);
	   //pf2->AddPeak(pM6);
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
   Double_t w = 600;
   Double_t h = 600;
   auto c1 = new TCanvas("c", "c", w, h); 

   if(!filecreated)
   {
         exp << 1 << endl;
         exp << 2 << endl;
   }

   
   
   for(int i=start;i<stop;i++){ // nang are the total number of angles in experiment.  //SKIPS THE 0 DEGREE ANGLE!!!
   
         
            std::cout << "\nHistogram: " << name <<  i <<"\n"<< std::endl;
      
         
        //  |======================================================================================|[3][Fitting of the NORMAL AC]
        
            name_tmp1=name+std::to_string(i);
            TH1D *Histo = (TH1D*)finput->Get(name_tmp1.c_str()); // Opens each Histogram in the Spectra.root file given.
               
            name_tmp1 = name_tmp1 + ".root";
           
            Histo->SetAxisRange(lowProj-5,highProj+5);
            
            p1->GetFitFunction()->FixParameter(5,0);
            //p1->GetFitFunction()->SetParLimits(5,538.5,540);
            
            pf->Fit(Histo,"REM");//fits the addbackAddback projected histogram
      
            c1->Modified();
            c1->Update();
	
      	   if(filecreated)  c1->WaitPrimitive(); 
            //Will wait for user to update the canvas (via the click of enter) as there is a 
            //probability that the fit has initially failed.
            
            if(p1->GetChi2()/p1->GetNDF()>30)
            {
               c1->WaitPrimitive(); 
   	         Badfit = true;
               c1->SaveAs(name_tmp1.c_str());
            }
	
	//  |======================================================================================|[3] [Fitting of the NORMAL AC]
	
	//  |======================================================================================|[3][Fitting of the EVENT MIXED AC]
	    name_tmp2 = name+"Mixed"+std::to_string(i);
            TH1D *Histomixed = (TH1D*)finput->Get(name_tmp2.c_str());
            
            name_tmp2 = name_tmp2 + ".root";
      	
      	    p2->GetFitFunction()->FixParameter(5,0);
   	
           
            Histomixed->SetAxisRange(lowProj-5,highProj+5);
            pf2->Fit(Histomixed,"REM"); //fits the addbackAddbackMIXED projected histograms
            c1->Modified();
            c1->Update();
          
            if(filecreated)  c1->WaitPrimitive();
      
            if(p2->GetChi2()/p2->GetNDF()>30)
            {
               c1->WaitPrimitive(); 
   		      BadfitMixed = true;
               c1->SaveAs(name_tmp2.c_str());
            }
         //  |======================================================================================|[3][Fitting of the EVENT MIXED AC]   


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
               foutup << i+3 << "\t" <<p1->Area()<<"\t"<< p2->Area() << endl;
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
           
                     line1 = to_string(i+3) + "\t" + p1->Area() + "\t" + p2->Area();

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

