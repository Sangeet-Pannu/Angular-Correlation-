#include "Riostream.h"
#include <iostream>
#include <assert.h>
#include "TH1.h"
#include "TFile.h"
//#include "TMath.h"
#include "RooMath.h"
#include <vector>


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

int j=9; // additional peaks to be added to system besides groundstate.

/*
 * Below block holds the Fit model definition of the code.
 * Based off of the Electron_Model peak shape.
 * 
 * */

//------------------------------------------
//------------------------------------------
double ElectronModel_NFit(double *x,double* par)
{
    
    bool reject=false;
    if (reject && x[0] > 1580 && x[0] < 1750) {
      TF1::RejectPoint();
      return 0;
    }
   
 
    
    // Parameters:
    // par[0]  = Centroid 1st peak
    // par[1]  = Amplitude 1st peak
    
    // par[9]  = Centroid 2nd peak
    // par[10] = Amplitude 2nd peak
    
    // par[2]  = Sigma
    // par[3]  = Scale factor 1st left tail
    // par[4]  = Decay constant 1st left tail
    // par[5]  = Scale factor 2nd left tail
    // par[6]  = Decay constant 2nd left tail
    // par[7]  = Scale factor right tail
    // par[8]  = Decay constant right tail

    // First peak
    // First term, Gaussian
    double f1 = ( 1/(sqrt(2*TMath::Pi())* par[2]) ) * exp( -0.5*( TMath::Power((x[0]-par[0])/par[2],2) ));
    // First left tail
    double skew1_1 = 1-erf( ( (par[2]*par[2])+par[4]*(x[0]-par[0]) ) / (sqrt(2)*par[2]*par[4]) );
    double f1_1 = ( par[3]/(2* par[4]) ) * exp( ( par[2]*par[2]+2*par[4]*x[0]-2*par[4]*par[0] ) / ( 2*par[4]*par[4] ) );
    // Second left tail
    double skew1_2 = 1-erf( ( (par[2]*par[2])+par[6]*(x[0]-par[0]) ) / (sqrt(2)*par[2]*par[6]) );
    double f1_2 = ( par[5]/(2* par[6]) ) * exp( ( par[2]*par[2]+2*par[6]*x[0]-2*par[6]*par[0] ) / ( 2*par[6]*par[6] ) );
    // Right tail    
    double skew1_3 = 1-erf( ( (par[2]*par[2])-par[8]*(x[0]-par[0]) ) / (sqrt(2)*par[2]*par[8]) );
    double f1_3 = ( par[7]/(2* par[8]) ) * exp( ( par[2]*par[2]-2*par[8]*x[0]+2*par[8]*par[0] ) / ( 2*par[8]*par[8] ) );    
    // Total peak
    double gaus1 = par[1]*( f1*(1-par[3]-par[5]-par[7])+f1_1*skew1_1+f1_2*skew1_2+f1_3*skew1_3);

    double sumfunc = gaus1;
    

    for(int i = 9; i<=(2*j+9); i+=2)
    {
    	f1 = ( 1/(sqrt(2*TMath::Pi())* par[2]) ) * exp( -0.5*( TMath::Power((x[0]-par[i])/par[2],2) ));
    	// First left tail
    	skew1_1 = 1-erf( ( (par[2]*par[2])+par[4]*(x[0]-par[i]) ) / (sqrt(2)*par[2]*par[4]) );
    	f1_1 = ( par[3]/(2* par[4]) ) * exp( ( par[2]*par[2]+2*par[4]*x[0]-2*par[4]*par[i] ) / ( 2*par[4]*par[4] ) );
    	// Second left tail
    	skew1_2 = 1-erf( ( (par[2]*par[2])+par[6]*(x[0]-par[i]) ) / (sqrt(2)*par[2]*par[6]) );
    	f1_2 = ( par[5]/(2* par[6]) ) * exp( ( par[2]*par[2]+2*par[6]*x[0]-2*par[6]*par[i] ) / ( 2*par[6]*par[6] ) );
   	 // Right tail    
    	skew1_3 = 1-erf( ( (par[2]*par[2])-par[8]*(x[0]-par[i]) ) / (sqrt(2)*par[2]*par[8]) );
    	f1_3 = ( par[7]/(2* par[8]) ) * exp( ( par[2]*par[2]-2*par[8]*x[0]+2*par[8]*par[i] ) / ( 2*par[8]*par[8] ) );    
    	// Total peak
    	gaus1 = par[i+1]*( f1*(1-par[3]-par[5]-par[7])+f1_1*skew1_1+f1_2*skew1_2+f1_3*skew1_3);
    
    	sumfunc += gaus1;
    
    }
	
    // Sum
    return sumfunc;
}

double ElectronModel_singles(double *x,double* par)
{
    double f1 = ( 1/(sqrt(2*TMath::Pi())* par[2]) ) * exp( -0.5*( TMath::Power((x[0]-par[0])/par[2],2) ));
    // First left tail
    double skew1_1 = 1-erf( ( (par[2]*par[2])+par[4]*(x[0]-par[0]) ) / (sqrt(2)*par[2]*par[4]) );
    double f1_1 = ( par[3]/(2* par[4]) ) * exp( ( par[2]*par[2]+2*par[4]*x[0]-2*par[4]*par[0] ) / ( 2*par[4]*par[4] ) );
    // Second left tail
    double skew1_2 = 1-erf( ( (par[2]*par[2])+par[6]*(x[0]-par[0]) ) / (sqrt(2)*par[2]*par[6]) );
    double f1_2 = ( par[5]/(2* par[6]) ) * exp( ( par[2]*par[2]+2*par[6]*x[0]-2*par[6]*par[0] ) / ( 2*par[6]*par[6] ) );
    // Right tail    
    double skew1_3 = 1-erf( ( (par[2]*par[2])-par[8]*(x[0]-par[0]) ) / (sqrt(2)*par[2]*par[8]) );
    double f1_3 = ( par[7]/(2* par[8]) ) * exp( ( par[2]*par[2]-2*par[8]*x[0]+2*par[8]*par[0] ) / ( 2*par[8]*par[8] ) );    
    // Total peak
    double gaus1 = par[1]*( f1*(1-par[3]-par[5]-par[7])+f1_1*skew1_1+f1_2*skew1_2+f1_3*skew1_3);

    return gaus1;
}

double background(double *x,double* par){
    return (par[0]*x[0]) + par[1];
}

double EModelwBK(double *x,double* par){
    return background(x,par) + ElectronModel_singles(x,&par[2]);
}

double FitFunction(double *x,double* par){

    return background(x,par) + ElectronModel_NFit(x,&par[2]);
    
}

//------------------------------------------
//------------------------------------------

void DrawFits(TF1* fitfunc,double lowlimit, double highlimit)
{
	
    Double_t pars[11+(2*j)];

    double par1 = fitfunc->GetParameter(0);
    double par2 = fitfunc->GetParameter(1);

    Double_t par[9];
     
    TF1 *bkfunc = new TF1("bkfunc",background,lowlimit,highlimit,2);
    bkfunc->SetParameter(0,par1);
    bkfunc->SetParameter(1,par2);
    bkfunc->SetLineColor(kBlack);
    bkfunc->Draw("SAME");

    par[0]=fitfunc->GetParameter(2);
    par[1]=fitfunc->GetParameter(3);
    par[2]=fitfunc->GetParameter(4);
    par[3]=fitfunc->GetParameter(5);
    par[4]=fitfunc->GetParameter(6);
    par[5]=fitfunc->GetParameter(7);
    par[6]=fitfunc->GetParameter(8);
    par[7]=fitfunc->GetParameter(9);
    par[8]=fitfunc->GetParameter(10);
    
    TF1 *SG1 = new TF1("SG1",EModelwBK,lowlimit,highlimit,11);
    SG1->SetParameter(0,par1);
    SG1->SetParameter(1,par2);
    SG1->SetParameter(2,par[0]);
    SG1->SetParameter(3,par[1]);
    SG1->SetParameter(4,par[2]);
    SG1->SetParameter(5,par[3]);
    SG1->SetParameter(6,par[4]);
    SG1->SetParameter(7,par[5]);
    SG1->SetParameter(8,par[6]);
    SG1->SetParameter(9,par[7]);
    SG1->SetParameter(10,par[8]);
    
    SG1->SetLineColor(kGreen);
    SG1->Draw("SAME");
    
  
    int l = 1;   
    for(int i=11; i < (2*j+11); i+=2)
    {
    	par[0]=fitfunc->GetParameter(i);
		par[1]=fitfunc->GetParameter(i+1);
    	
    	TF1 * func_draw = new TF1("func_draw",EModelwBK,lowlimit,highlimit,11);
    	func_draw ->SetParameter(0,par1);
    	func_draw ->SetParameter(1,par2);
    	func_draw ->SetParameter(2,par[0]);
    	func_draw ->SetParameter(3,par[1]);
    	func_draw ->SetParameter(4,par[2]);
    	func_draw ->SetParameter(5,par[3]);
		func_draw ->SetParameter(6,par[4]);
		func_draw ->SetParameter(7,par[5]);
		func_draw ->SetParameter(8,par[6]);
		func_draw ->SetParameter(9,par[7]);
    	func_draw ->SetParameter(10,par[8]);
		func_draw ->SetLineColor(l);
    	func_draw ->Draw("SAME");
    	l++;
    	
    
    }
    
    

  
  
 }

/*
TH1* OpenFile() Function's functionality is opening ascii files and writing them to
.root files via TFile *file = new TFile().

To change ascii file, simply remove pre-existing file in "in.open()" and place the specific ascii 
file to be read in the brackets.

The naming of the root file can be changed by altering the content of line 97.
*/
   // Opening the ASCII file Function
    TH1* OpenFile(string name = "/home/sangeet/Desktop/Zr_experiment/92Zr_ascii_files/92Zr_34/replay_run162_163_92Zr_34_hXew_1-inverted-1col.ascii") {

    ifstream in;

    in.open(name.c_str());
    
    float y;
    TFile *file = new TFile("96Zr_Deg45.root","RECREATE");

    TH1F*hist = new TH1F("hist","94Zr_Deg41",2000,0,2000);
   
    int i=0;
    while (1) {
      in >> y ;
      if (!in.good()) break;    
      hist->SetBinContent(hist->GetXaxis()->FindBin(i),y); 
      i++;
    }
    
    in.close();
 
    return hist;

}

void DisplayArea(TF1* fitfunc,double rebin, int status[10]){

    cout << "---------------- Area from Amplitudes ------------------------" <<endl;

    cout << endl << "Rebin = " << rebin << endl << endl;  
    double cen1 = fitfunc->GetParameter(2);
    double area_gs = fitfunc->GetParameter(3)/rebin;
    double errarea_gs = fitfunc->GetParError(3)/rebin;
    cout << "Centroidal Value: " << cen1 << endl;
    cout << "Area of Ground-State = " << area_gs << " +- " << errarea_gs << " (err pc = " << 100*errarea_gs/area_gs << ")" << endl << endl;
    
	double area_21;
	double errarea_21;
	string Areastrn = "Area of ";
    
    	
    int l = 1;
    int p = 1;
    int p1 = 1;
    for(int i = 12; i<=(2*j+10); i+=2)
    {
    	double ceni = fitfunc->GetParameter(i-1);
    	area_21 = fitfunc->GetParameter(i)/rebin;
    	errarea_21 = fitfunc->GetParError(i)/rebin;
    	if(status[l] == 1)
    	{
    	cout << "Centroidal Value: " << ceni << endl;
    	cout << Areastrn+" Excited State "+p <<" = "<< area_21 << " +- " << errarea_21 << " (err pc = " << 100*errarea_21/area_21 << ")" << endl << endl;
    	p++;
    	}else{cout << "Centroidal Value: " << ceni << endl; cout << Areastrn+" impurity "+p1 <<" = "<< area_21 << " +- " << errarea_21 << " (err pc = " << 100*errarea_21/area_21 << ")" << endl << endl;p1++;}
	l++;
    }
}
 

void Integration(TF1* fitfunc,TH1* fhist, double lowlimit, double highlimit, double rebin, int status[10]){
    
    
    TFitResultPtr fr = fhist->Fit("fitfunc","RSEM");
    
    fitfunc->SetParameters(fitfunc->GetParameters());
   	
    //residual plot
    auto c1 = new TCanvas("c1", "Fit spectrum with Residual plot");
    fhist->Fit("fitfunc","RSEM");  
    c1->Clear(); // Fit does not draw into correct pad
    auto rp1 = new TRatioPlot(fhist);
    rp1->Draw();

    rp1->GetLowerRefYaxis()->SetTitle("Residual");
    c1->Update();
    
    
    auto c2 = new TCanvas("c2", "Fit spectrum with individual fits");
    fhist->Draw();
    DrawFits(fitfunc,lowlimit,highlimit);
    DisplayArea(fitfunc,rebin,status);
   
   

}

/*
 The block of code below pretains to the click event functionality in
 * the fit software.
 */
//---------------------------
float x = 0, y = 0;
bool check = true;
void Clicked() {
 	int event = gPad->GetEvent();
	int px = gPad->GetEventX();
	int py = gPad->GetEventY();	
 	//cout << "event # = " << event << endl;
   	if (event==61) { // Upon double clicking
   		//if (!gPad->GetSelected()) return;
   		
      		Float_t xx = gPad->AbsPixeltoX(px);
     		Float_t yy = gPad->AbsPixeltoY(py);
			x  = gPad->PadtoX(xx);
     		y  = gPad->PadtoY(yy);
		return;
	}
	else if(event ==24){ //Any key will cause system to exit 
		check = false;
		
		return;
	}
}

void trig(TCanvas *c){
	
	while((x == 0 or y ==0)){
		c->AddExec("ex","Clicked()");
		gPad->WaitPrimitive();	
		c->DeleteExec("Clicked()");
		
		
		
	}
}
//---------------------------

void FitCode(){
 
 	   TCanvas *c1 = new TCanvas("c1","c1",200,10,800,600);
 	   //cout<< "Enter File (ascii) directory location+name: "<<endl;
 	   //cin>>filename;
 
	    TH1* fhist = OpenFile();
	    fhist->Draw(); //Position changed from ***
	    gPad->Update();
	    
	    double limits[2];
	    //Click positions for integration domain
	    for(int i = 0; i<2; ++i){
			trig(c1); 
			limits[i] = x;
			TLine *line1 = new TLine(x,0,x,10000);
			line1->SetLineColor(kGreen);
			line1->Draw();
			x = 0;
			y = 0;
	    }
	    
	    cout << "Limit position 1: "<<limits[0]<<endl;
	    cout << "Limit position 2: "<<limits[1]<<endl;
	    
	    double limitsbkx[2];
	    double limitsbky[2];
	    //Click positions for background subtraction domain.
	    string BKsel = "y";
	    //To ask user whether the slope is be flat or with some slope.
	    
	    
	    for(int i = 0; i<2; ++i){
			trig(c1); 
			limitsbkx[i] = x;
			limitsbky[i] = y;
			TLine *line1 = new TLine(x,0,x,10000);	
			line1->SetLineColor(kBlack);
			line1->Draw();
			x = 0;
			y = 0;
	    }
	    
	    TLine *line2 = new TLine(limitsbkx[0],limitsbky[0],limitsbkx[1],limitsbky[0]);
	    line2->SetLineColor(kBlack);
		line2->Draw();
		
	    cout << "BK Limit position 1: "<<limitsbkx[0]<<endl;
	    cout << "BK Limit position 2: "<<limitsbkx[1]<<endl;
	    
	    
	    float slope;
	    if(BKsel == "n"){
			slope = (limitsbky[1]-limitsbky[0])/(limitsbkx[1]-limitsbkx[0]);
			
		};//Computes the slope if user requires a slope on BK.
	    
	    
	    cout << "NOTE: that the first peak selected must be the ground state" << endl;
	    cout << "furthermore, the second peak should be the 1st excited state." << endl;
	    //Centrodial peak value selector via clicks
	    std::vector< float > centval; //Plays as a dynamic array; (however takes more memory)
	    check = true;
	    /*
	     * While statement below allows for the click method for Cent.
	     * */
	    while(check){
			trig(c1); 
			if(!check){break;};
			centval.push_back(x);
			TLine *line1 = new TLine(x,0,x,1e5);
			line1->SetLineColor(kRed);
			line1->Draw();
			x = 0;
			y = 0;
	    }
	    
	   
		

	    // "lowlimit" ~ Lower boundary for fitting domain
	    // "Highlimit" ~ Higher boundary for fitting domain
	    // Changing these parameters, changes the fitting domain.
	    
	    j = centval.size()-1; //minus 1 coming from the fact that GS peak and 1st excited is included.
	    
	    int status[11+(2*j)];
	    int par = 11+(2*j);
	
	    double lowlimit =limits[0];
	    double highlimit = limits[1];
	   
	    double rebin =2;
	    fhist->Rebin(rebin);  


	    TF1 *fitfunc = new TF1("fitfunc",FitFunction,lowlimit,highlimit,par); 

	    fitfunc->SetParName(0,"BK Slope");
	    fitfunc->SetParName(1,"BK Y-Intercept");
	    fitfunc->SetParName(2,"Ground-State Cent. #1");
	    fitfunc->SetParName(3,"Area 1");
	    fitfunc->SetParName(4,"Sigma");
	    fitfunc->SetParName(5,"Scale 1l");
	    fitfunc->SetParName(6,"Decay 1l");
	    fitfunc->SetParName(7,"Scale 2l");
	    fitfunc->SetParName(8,"Decay 2l");
	    fitfunc->SetParName(9,"Scale r");
	    fitfunc->SetParName(10,"Decay r");
	    
	    /*
	     * This below block of code is to discern the impurity peaks
	     * from the "real" peaks.
	     * 
	     * UNDER CONSTRUCTION
	     * */
	     
	    //--------------------------------------------------------------
	    for(int i=11; i < (2*j+11); i+=2)
	    {
	    	string cent_name = "Centroid ";
	    	string area_name = "Area ";

	    	fitfunc->SetParName(i,cent_name+(i-9));  	
	    	fitfunc->SetParName(i+1,area_name+(i-9));
	    	
	    }
	    //--------------------------------------------------------------
	    
	    
	    // Background Fitting parameters--------------------------------	
	    
		fitfunc->FixParameter(0,0);
	    fitfunc->FixParameter(1,limitsbky[0]); // Background offset
	     
		// Background Fitting parameters--------------------------------	
		
		/*
		 * 
		 * This block distinishes the ground-state peak parameter limits
		 * as the remaining peaks of the spectrum should have the exact 
		 * peak shape parameters. 
		 * 
		 */
		
		//-----------------------------------------------------------------------
	    fitfunc->FixParameter(2,centval[0]); // Centroid 1st peak
	    fitfunc->SetParLimits(2,centval[0]-15,centval[0]+15); //Parameter limiter on ground state centroid.
	    fitfunc->SetParameter(3,500000); // Amplitude 1st peak
	    fitfunc->SetParameter(4,9.61468e+00); // Sigma
	    fitfunc->FixParameter(5,7.85495e-02); // Scale factor 1st left tail
	    fitfunc->SetParLimits(5,0,1);
	    fitfunc->FixParameter(6,6.76023e+01); // Decay constant 1st left tail
	    fitfunc->SetParLimits(6,0,1e+03);
	    fitfunc->FixParameter(7,6.57981e-01); // Scale factor 2nd left tail
	    fitfunc->SetParLimits(7,0,1);
	    fitfunc->FixParameter(8,1.82368e+01); // Decay constant 2nd left tail
	    fitfunc->SetParLimits(8,0,1e+03);
	    fitfunc->FixParameter(9,3.42822e-02); // Scale factor right tail
	    fitfunc->SetParLimits(9,0,1);
	    fitfunc->FixParameter(10,8.52217e+00); // Decay constant right tail
	    fitfunc->SetParLimits(10,0,1e+03);
	    //-----------------------------------------------------------------------
	    
	    fitfunc->FixParameter(11,centval[1]);
	    fitfunc->SetParLimits(11,centval[1]-5,centval[1]+5);
	    fitfunc->SetParameter(12,5000); // 2nd Peak
	    
	    for(int i=1;i<j;i++){
			fitfunc->FixParameter(11+(2*i),centval[i+1]); // Input Parameters: Parameter index, Centroidal Value.
			fitfunc->SetParLimits(11+(2*i),centval[i+1]-10,centval[i+1]+10); //Input Parameters: Parameter index, Centroidal Value Min, Centroidal Value Max.
			fitfunc->SetParameter(12+(2*i),500); // Input Parameters: Parameter index, Amplitude Value.
		}
		
	    gStyle->SetOptFit(0000);
	    fhist->SetStats(0);
	    
	    Integration(fitfunc,fhist,lowlimit,highlimit,rebin,status);// Input Parameters: Fit function, Histogram to be fit, Low fit limit, High fit limit, rebin value, centoridal Status (repair)
	    
	    
	    
	    
}




















