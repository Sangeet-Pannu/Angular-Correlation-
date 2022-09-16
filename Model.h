#include "TMath.h"

int l = 9;


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
	    

	    for(int i = 9; i<=(2*l+9); i+=2)
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
