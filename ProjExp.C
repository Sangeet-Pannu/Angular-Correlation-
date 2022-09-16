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
;
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
         proj->Add(projBG,-1*((highProj-lowProj)/(highBGProj-lowBGProj)));
      }
      if(subtractBG2) {
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

/*
ProjFit:

This Function takes a 2D matrix and preforms a projection of the Y-axis whilst 
subtracting the background from the Timerandom matrix (background).

name: The common name across all the histograms, if you use AngularCorrelationSelectorAddback.C then the common name is "addbackAddback"
lowProj: The low edge of projection gate
highProj: The high edge of projection gate

errin: ignore

start/stop: the start and stop of the which angle to begin and end the fitting process.

nang: number of angles

after each fit of the matrix, the fit will stop and to continue simply press "Enter". THis is to view the fit.
*/

void ProjFit(std::string name, float cent = 1., float lowProj=1., float highProj=1.,bool errin = false,int nang = 50, int start = 1,int stop =50) {// Input parameter: "histogram common name", centroidal value, lower projection, high projection, number of angles

   std::vector<std::string> lines;
   TFile *finput = (TFile*)gDirectory->GetFile();
   std::string name_tmp1;
   std::string name_tmp2;
   bool Badfit = false;
   bool BadfitMixed = false;
   bool determine = false;
   bool filecreated = false;

   ofstream exp;
   std::ifstream input_file;
   

   //Gaussian Peak, with the argument being the centroid of the plot
   
   TRWPeak *p1 = new TRWPeak(cent); //Centroidal value defined by user.
   TRWPeak *pM1 = new TRWPeak(581.064);
   //TRWPeak *pL1 = new TRWPeak(631);
   
   // Define fitter
   TPeakFitter *pf = new TPeakFitter(lowProj,highProj);
   // Adding peaks
   pf->AddPeak(p1);
   //pf->AddPeak(pM1);
   //pf->AddPeak(pL1);
   
   TRWPeak *p2 = new TRWPeak(cent); //Centroidal value defined by user.
   TRWPeak *pM2 = new TRWPeak(815.63); //Centroidal value defined by user.
   TRWPeak *pM3 = new TRWPeak(831.592); //Centroidal value defined by user.
   TRWPeak *pM4 = new TRWPeak(835.88); //Centroidal value defined by user.
   //TRWPeak *pM5 = new TRWPeak(836.509); //Centroidal value defined by user.
   //TRWPeak *pM6 = new TRWPeak(711.356); //Centroidal value defined by user.
   
   // Define fitter
   TPeakFitter *pf2 = new TPeakFitter(810,845);

   // Adding peaks
   pf2->AddPeak(p2);
   pf2->AddPeak(pM2);
   pf2->AddPeak(pM3);
   pf2->AddPeak(pM4);
   //pf2->AddPeak(pM5);
   //pf2->AddPeak(pM6);
   

   //I want to make it check whether the exp file already exists to be able to update those.
   ofstream fout ("Fitted_Spectra.txt", std::ofstream::app);

  
   if (fopen("exp.txt", "r")) { 
   std::cout << "exp.txt exists! Updating of the previous will occur" << std::endl;
   filecreated = true;
   input_file.open("exp.txt");
   std::string input;
   while (std::getline(input_file, input)) lines.push_back(input);  
   }else{
   exp.open("exp.txt");
   }

   ofstream foutup ("Updated_Fitted_spectra.root");

   
   if(!filecreated)fout << "Fitted Results Below" << endl;
   if(filecreated) foutup << "Updated Fitted Results Below" << endl;

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
   
         
            std::cout << "Histogram: " << name <<  i << std::endl;
            //----------------------------------------------------------------------------------
               
            name_tmp1=name+std::to_string(i);
            TH1D *Histo = (TH1D*)finput->Get(name_tmp1.c_str()); // Opens each Histogram in the Spectra.root file given.
               
            name_tmp2 = name+"Mixed"+std::to_string(i);
            TH1D *Histomixed = (TH1D*)finput->Get(name_tmp2.c_str());
               
            //----------------------------------------------------------------------------------
            
         
            //----------------------------------------------------------------------------------
            name_tmp1 = name_tmp1 + ".root";
            Histo->SetAxisRange(lowProj-50,highProj+50);
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

            name_tmp2 = name_tmp2 + ".root";
      
      
            //if(determine){ //runs once 
            //pf2->Fit(Histomixed,"QN0"); //fits the addbackAddbackMIXED projected histograms
            //cout << "Fixing Centroid Parameters" << endl;
            //p2->GetFitFunction()->SetParLimits(1,615.,616.5);//For TAB3Peak
            //cout << "Fixed 611. Centroid Parameter" << endl;
            //pM2->GetFitFunction()->SetParLimits(1,580,582);//For TAB3Peak
            pM2->GetFitFunction()->FixParameter(1,815.68);
            //pM3->GetFitFunction()->FixParameter(1,598.333);
            //pM4->GetFitFunction()->FixParameter(1,832.149);
            //pM5->GetFitFunction()->FixParameter(1,836.304);
            //cout << "Fixed 621. Centroid Parameter" << endl;
            //pM3->GetFitFunction()->SetParLimits(1,815,816.2);//For TAB3Peak
            //pM4->GetFitFunction()->SetParLimits(1,831,832);//For TAB3Peak
            //pM5->GetFitFunction()->SetParLimits(1,835,836);//For TAB3Peak
            //pM6->GetFitFunction()->FixParameter(1,711.356);//For TAB3Peak
            //cout << "Fixed 611. Centroid Parameter" << endl;
            //pM1->GetFitFunction()->SetParLimits(1,610.,612.);//For TAB3Peak
            //cout << "Fixed 611. Centroid Parameter" << endl;
            //pL1->GetFitFunction()->SetParLimits(1,630.,632.);//For TAB3Peak
            //cout << "Fixed 611. Centroid Parameter" << endl;
            //p1->GetFitFunction()->SetParLimits(1,615.,617.);//For TAB3Peak
            //cout << "Fixed 621. Centroid Parameter" << endl;
            //determine = false;
            //}
            
            Histomixed->SetAxisRange(lowProj-50,highProj+50);
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
            


            //----------------------------------------------------------------------------------
            
               
            //----------------------------------------------------------------------------------
            float sqrtred1 = TMath::Sqrt(p1->GetChi2()/p1->GetNDF());
            float sqrtred2 = TMath::Sqrt(p2->GetChi2()/p2->GetNDF());
            if(!filecreated)
            {
               fout << "addbackAddback_Projected:" << i <<"\t\t" << "addbackAddbackMixed_Projected:" << i << endl;
               if(Badfit){fout <<"---------BAD FIT---------" << endl;}
               if(BadfitMixed){fout <<"---------BAD MIXED FIT---------" << endl;}
               if(!errin){fout << " Area:\t"<<p1->Area()<<"\t+/-\t"<<p1->AreaErr() <<"\t\t"<< " Mixed Area:\t"<<p2->Area()<<"\t+/-\t"<<p2->AreaErr()<< endl;}
            }

            if(filecreated)
            {
               foutup << "addbackAddback_Projected:" << i <<"\t\t" << "addbackAddbackMixed_Projected:" << i << endl;
               if(Badfit){foutup <<"---------BAD FIT---------" << endl;}
               if(BadfitMixed){foutup <<"---------BAD MIXED FIT---------" << endl;}
               if(!errin){foutup << " Area:\t"<<p1->Area()<<"\t+/-\t"<<p1->AreaErr() <<"\t\t"<< " Mixed Area:\t"<<p2->Area()<<"\t+/-\t"<<p2->AreaErr()<< endl;} 
            }   
               
            //----------------------------------------------------------------------------------
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
   }
         //The writing of the updated file
            if(filecreated)
            {
              input_file.close();

               std::ofstream output_file("exp.txt");
               for (auto const& line : lines)
               output_file << line << '\n';
            
               exp.close();          
            }

            Badfit = false;
            BadfitMixed = false;  
         
   fout.close();


}


void ProjSimFit(std::string name,std::string name2,float cent = 1., float lowProj=1., float highProj=1., int nang = 52, int start = 1) {// Input parameter: "histogram common name", centroidal value, lower projection, high projection, number of angles
   
   TFile *finput = (TFile*)gDirectory->GetFile();
   std::string name_tmp1;
   std::string name_tmp2 = name2 + ".txt";
   
   //Gaussian Peak, with the argument being the centroid of the plot
   
   TABPeak *p1 = new TABPeak(cent); //Centroidal value defined by user.
   // Define fitter
   TPeakFitter *pf = new TPeakFitter(lowProj,highProj);
   // Adding peaks
   pf->AddPeak(p1);
   
   
   

      
   
   //----------------------------------------------------------------------------------
   
   ofstream fout ("Fitted_Spectra.txt");
   ofstream exp (name_tmp2);
   ofstream maps ("ANGLE_MAPP.txt");
   
   fout << "Fitted Results Below" << endl;
   //----------------------------------------------------------------------------------
   Double_t w = 600;
      Double_t h = 600;
      auto c1 = new TCanvas("c", "c", w, h); 
      
      
      float Areas[nang];
      float Areaserr[nang];
   for(int i=start;i<nang;i++)
   { // nang are the total number of angles in experiment.  //SKIPS THE 0 DEGREE ANGLE!!!
   
         
            std::cout << "Histogram: " << name <<  i << std::endl;
            //----------------------------------------------------------------------------------
            
            name_tmp1=name+std::to_string(i);
            TH1D *Histo = (TH1D*)finput->Get(name_tmp1.c_str()); // Opens each Histogram in the Spectra.root file given.
               
               
            //----------------------------------------------------------------------------------
               
               
            //----------------------------------------------------------------------------------
         name_tmp1 = name_tmp1 + ".root";
            pf->Fit(Histo,"REM");//fits the simulated projected histograms
            //c1->SaveAs(name_tmp1.c_str());
   
            //----------------------------------------------------------------------------------
            
               
            //----------------------------------------------------------------------------------
            fout << "Histogram" << i << endl;
            fout << " Centroidal Value: "<<p1->Centroid()<<"+/-"<<p1->CentroidErr() << endl;
      fout << " Area: "<<p1->Area()<<"+/-"<<p1->AreaErr() << endl;
      fout << " Chi2/NDF: "<<p1->GetChi2()/p1->GetNDF()<< endl;
      fout <<"--------------------------------------|"<<endl;
         //----------------------------------------------------------------------------------
               
               
            //----------------------------------------------------------------------------------
            Areas[i-1] = p1->Area();
            Areaserr[i-1] = p1->AreaErr();
            
            
            
            
            //----------------------------------------------------------------------------------
            
            
   }
   
   //----------------------------------------------------------------------------------
   std::map<float, float, std::less<float>> angles;
   std::map<float, float, std::less<float>> angles1;
   //std::map<float, float, std::less<float>> anglemap;
   
   int x; // variable for input value
   int i = 0;
   int y;
         
   float angle_val;        
   double_t value;
   double_t value1;
      
   string name1 = "/home/Sangeet/Desktop/Master_project/100Zr/Sim_files/Real_Files/id_to_angle.txt";
            
   ifstream indata; // indata is like cin 
   indata.open(name1); // opens the file
      
   if(!indata) 
   { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
         exit(1);
   }
      
      
   while ( !indata.eof() ) // keep reading until end-of-file.   
   { 
      //Below is a descripiton of the content inside id_to_angle.txt
      /*
      | spectrum id | angle [deg] | #pairs |
      */
      
      indata >> x >> angle_val >> y; // takes in the values from the ascii file id_to_angle.txt 
            
      if(angle_val == 0) break;// the eof places a value of 0 to "angle_val" after reading everything which places\
               a false entry in final results. this is to counter such situations.
               
      /*
      angles:
         First entry is the angle value and the second is the area divided by number of pairs.
         
      angles1:
         First entry is the angle value and the second is the area ERROR divided by number of pairs.
         
      Maps are used as they have the inherent ability to sort in descending order the "key" value by which 
      the order of the angles are corrected!!
      */
      //anglemap.insert({angle_val,x});
      angles.insert({angle_val,Areas[i]/y}); 
      angles1.insert({angle_val,Areaserr[i]/y});
      i++;
   
   }
      
         //----------------------------------------------------------------------------------
         
         int k = 1; //input counter for the output file ___.txt.
         float a; // variable for the area error of the respective area.
         exp.precision(4);    
         for (auto itr = angles.begin(); itr != angles.end(); ++itr) {
         a = angles1.at(itr->first);// finds the 
         
         //the Below populates the exp file (named: ___.txt) with index | Area | Area_error
         exp << k <<"\t"<< itr->second << "\t"<< a << endl;
         k++;
      }
      /*
      int count = 1;
      for (auto itr = anglemap.begin(); itr != anglemap.end(); ++itr) {
         
         maps << count <<"\t"<< itr->first <<"\t"<< itr->second << endl;
         count++;
      }
      */
   //----------------------------------------------------------------------------------
   fout.close();
   exp.close();
   maps.close();
   //----------------------------------------------------------------------------------
}


