//
//  Chi2_Mix_Xtract.cxx
//  Program to extracting spins and Mixing ratios
//  of x-2-0 transition cascades.
//  ** Modification of ChiXe126_x-2-0 which was Created by Sambuu on 2019-07-10. **
//  Created by Sangeet-Pal on 2022-09-19.
//

#include <array>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <limits.h>

#define PI 3.14159265



int main(){

    //---------------------
    std::ifstream infile("exp.txt");
    std::string lline;
    float y[29];
    float yr[29];
    int count1 = 0;
    if (infile.is_open())
        {
     
        while (std::getline(infile, lline))
        {
            std::stringstream ss(lline);
            float b, c;
            int a;
                        
            if (ss >> a >> b >> c)
            {
                y[count1] = b;
            yr[count1] = c;
            }
            count1++;
        }
        
        infile.close();
        
        }
        else
        {
        std::cout << "error: infile is not open" << std::endl;
        }   
    //---------------------

    for(int k=2;k<(count1-1);k++) std::cout << "Normalized Peak areas: " << y[k] << "\tCorresponding errors: " << yr[k] << std::endl;

    bool isFipps = true;
    bool isIfin = false;

    

    //Fipps: All the possible angles
    float F_angs[21] = {0.579857,0.615795,0.7812,0.809525,0.982443,1.00657,1.36004,
        1.38068,1.56068,1.58091,1.761,1.78155,2.13503,2.15915,2.33207,2.36039,2.5258,
        2.56174,2.85618,2.94012,3.14159};

    //Ifin: All the possible angles.
    float I_angs[29] = {0.795613, 0.909671,0.945734,0.995199,1.048,1.0742,1.081,1.15188,
        1.20067,1.29824,1.38323,1.39967,1.56271,1.57888,1.7419,1.75836,1.8433,1.94093,1.98971,
        2.06059,2.06739,2.09359,2.1463,2.19586,2.23192,2.34598,2.8865,2.96151,3.14159};

    float Wtheo[29];
    

    int J2 = 2, J3 = 0;
    float ac_min = 0;
    float ac_max = 0;
    int count = 0;
    float delta_max = 0;
    float delta_min = 0;
    float del_max = 0;
    float del_min = 0;
    float chi_min = 0;
    float chi_max = 0;
    float ch_min = 0;
    float ch_max = 0;
    float c_min = 0;
    float c_max = 0;
    float chi_min0 = 0;
    float chi_max0 = 0;
    
    float ChiSquare2 = 0;
    float ChiSquare_0 = 0;
    int countchi = 0;
    float delmin120 = 0;
    float delmin220 = 0;
    float delmin320 = 0;
    float c_min020 = 0;
    float c_min120 = 0;
    float c_min220 = 0;
    float c_min320 = 0;
    float c_min420 = 0;
    float c_max120 = 0;
    float c_max220 = 0;
    float c_max320 = 0;
    float delmax120 = 0;
    float delmax220 = 0;
    float delmax320 = 0;

    float Q2 = 1;//quenching factors
    float Q4 = 1;

    float chi_rp120 = 0; //chi^2+1
    float chi_rp220 = 0;
    float chi_rp320 = 0;
    float chi_rm120 = 0; //chi^1-1
    float chi_rm220 = 0;
    float chi_rm320 = 0;
    int cc = 0;
    float drp120 = 0;
    float drp220 = 0;
    float drp320 = 0;
    float drm120 = 0;
    float drm220 = 0;
    float drm320 = 0;

    float CV95 =2.37;
    int loopc;

    //==========================================================================================================|J|
    for (int J1 = 0; J1 < 1; J1++)
    {
	
	
	
    	std::stringstream chi0filename;
        chi0filename << "chi_square" << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Ru100.dat";
        std::string chi0name = chi0filename.str();

	    std::stringstream chifilename;
       
        chifilename << "_"<< J1 << "-" << J2 << "-" << J3 << "_Xe126.dat";
        std::string chiname = chifilename.str();

        //resets the wcoeff
        float wcoeff = 0;
        float wcoeff1 = 0;
        float wcoeff2 = 0;
        double ChiSquare = 0;
        int cchi = 0;
        int cochi = 0;
        int countchi0 = 0;

        //starting to find chi square distribution
        
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
        } else
            if (J1 == 2 && J2 == 2 && J3 == 0) {
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
            } else
                if (J1 == 3 && J2 == 2 && J3 == 0) {
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
                } else
                    if (J1 == 4 && J2 == 2 && J3 == 0) {
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
                        std::cout << "you entered wrong value for J" << std::endl;
                        return 0;
                    }
        
        float a22 = Q2 * R2LLJ2J1 * R2LLJ2J3;
        float a44 = Q4 * R4LLJ2J1 * R4LLJ2J3;
        
        
        std::ofstream chifile (chiname);
        //==================================================================================|2|
        if (chifile.is_open())
        {
            double delta1;
            for (delta1 = -50.000; delta1 < 50.000; delta1 += 0.010)
            {
                //Set to zero as the second gamma goes from 2_1+ to 0_1+ which is a pure E2 transition.
                int delta2 = 0;
                
                float a2 = Q2 * ((R2LLJ2J1 + 2 * delta1 * R2LMJ2J1 + pow(delta1, 2) * R2MMJ2J1) * (R2LLJ2J3 + 2 * delta2 * R2LMJ2J3 + pow(delta2, 2) * R2MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * delta1 * R4LMJ2J1 + pow(delta1, 2) * R4MMJ2J1) * (R4LLJ2J3 + 2 * delta2 * R4LMJ2J3 + pow(delta2, 2) * R4MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));

                //Definition of the theoritcal angular correlation results.
                //---------------------------------------------------------**
                if(isIfin){
                    loopc = 0;
                    for(const auto& ang : I_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }

                if(isFipps)
                {
                    loopc = 0;
                    for(const auto& ang : F_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }
                //========================================================**
            
        
           

             
                //=============================================================================***
                // "Normalization Coefficent"
                //wcoeff1 = SUM_{Y*W_theo/Yr^2}
                //wcoeff2 = SUM_{W_theo^2/Yr^2}
                //woceff = wcoeff1/wcoeff2 = UM_{Y*W_theo/Yr^2}/SUM_{W_theo^2/Yr^2}

                for(int i = 0; i<loopc; i++) wcoeff1 +=  (y[i+2]*Wtheo[i]/pow(yr[i+2],2));
                for(int i = 0; i<loopc; i++) wcoeff2 +=  pow(Wtheo[i]/yr[i+2],2);

                wcoeff = wcoeff1/wcoeff2;

                // "Calculation of Chi2/NDF value"
                for(int j = 0; j<loopc; j++) ChiSquare += pow((y[j+2]-wcoeff*Wtheo[j])/yr[j+2],2);

                ChiSquare = ChiSquare/(loopc-1); //chi2/NDF, where NDF = DOF - 1; 
                //=============================================================================***     
                
                //=============================================================================****              
                //starting to find min and max chi square values
                if (J1 == 1 && J2 == 2 && J3 == 0) {
                    if (countchi==0) {
                        chi_min = ChiSquare;
                        chi_max = ChiSquare;
                        delta_min = delta1;
                        delta_max = delta1;
                    }else
                        if (chi_min > ChiSquare){
                            chi_min = ChiSquare;
                            delta_min = delta1;
                        }else if(chi_max < ChiSquare){
                            chi_max = ChiSquare;
                            delta_max = delta1;
                        }
                    countchi++;
                }else if (J1 == 2 && J2 == 2 && J3 == 0) {
                    if (countchi==0) {
                        chi_min = ChiSquare;
                        chi_max = ChiSquare;
                        delta_min = delta1;
                        delta_max = delta1;
                    }else
                        if (chi_min > ChiSquare){
                            chi_min = ChiSquare;
                            delta_min = delta1;
                        }else if(chi_max < ChiSquare){
                            chi_max = ChiSquare;
                            delta_max = delta1;
                        }
                    countchi++;
                }else if (J1 == 3 && J2 == 2 && J3 == 0) {
                    if (countchi==0) {
                        chi_min = ChiSquare;
                        chi_max = ChiSquare;
                        delta_min = delta1;
                        delta_max = delta1;
                    }else
                        if (chi_min > ChiSquare){
                            chi_min = ChiSquare;
                            delta_min = delta1;
                        }else if(chi_max < ChiSquare){
                            chi_max = ChiSquare;
                            delta_max = delta1;
                        }
                    countchi++;
                }
                if (J1 == 1 && J2 == 2 && J3 == 0) {
                    if (cochi==0) {
                        c_min = ChiSquare;
                        c_max = ChiSquare;
                        delmax120 = delta1;
                        delmin120 = delta1;
                    }else
                        if (c_min > ChiSquare){
                            c_min = ChiSquare;
                            delmin120 = delta1;
                        }else if(c_max < ChiSquare){
                            c_max = ChiSquare;
                            delmax120= delta1;
                        }
                    c_min120 = c_min;
                    c_max120 = c_max;
                    cochi++;
                }else if (J1 == 2 && J2 == 2 && J3 == 0) {
                    if (cochi==0) {
                        c_min = ChiSquare;
                        c_max = ChiSquare;
                        delmax220 = delta1;
                        delmin220 = delta1;
                    }else
                        if (c_min > ChiSquare){
                            c_min = ChiSquare;
                            delmin220 = delta1;
                        }else if(c_max < ChiSquare){
                            c_max = ChiSquare;
                            delmax220 = delta1;
                        }
                    c_min220 = c_min;
                    c_max220 = c_max;
                    cochi++;
                }else if (J1 == 3 && J2 == 2 && J3 == 0) {
                    if (cochi==0) {
                        c_min = ChiSquare;
                        c_max = ChiSquare;
                        delmax320 = delta1;
                        delmin320 = delta1;
                    }else
                        if (c_min > ChiSquare){
                            c_min = ChiSquare;
                            delmin320 = delta1;
                        }else if(c_max < ChiSquare){
                            c_max = ChiSquare;
                            delmax320 = delta1;
                        }
                    c_min320 = c_min;
                    c_max320 = c_max;
                    cochi++;
                }
                if (cchi==0) {
                    ch_min = ChiSquare;
                    ch_max = ChiSquare;
                    del_min = delta1;
                    del_max = delta1;
                }else
                    if (ch_min > ChiSquare){
                        ch_min = ChiSquare;
                        del_min = delta1;
                    }else if(chi_max < ChiSquare){
                        ch_max = ChiSquare;
                        del_max = delta1;
                    }
                
                cchi++;
                chifile << atan (delta1) << "  " << ChiSquare << std::endl;
                
                //=============================================================================****
            } 
            
            chifile.close();        
        }
        //==================================================================================|2|

        //find uncetainty on mixing ratio based on 1 sigma uncertainty on chi^2
        double d1;

        if (J1 == 1 && J2 == 2 && J3 == 0) {
            for (d1 = delmin120; d1 <= 51.000; d1 += 0.010){
                ChiSquare2 = 0;
                float a2 = Q2 * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                
                if(isIfin){
                    loopc = 0;
                    for(const auto& ang : I_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }

                if(isFipps)
                {
                    loopc = 0;
                    for(const auto& ang : F_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }
                
               // "Normalization Coefficent"
                wcoeff1 = 0; 
                wcoeff2 = 0;

                for(int i = 0; i<loopc; i++) wcoeff1 +=  (y[i+2]*Wtheo[i]/pow(yr[i+2],2));
                for(int i = 0; i<loopc; i++) wcoeff2 +=  pow(Wtheo[i]/yr[i+2],2);

                wcoeff = wcoeff1/wcoeff2;

                // "Calculation of Chi2/NDF value"
                for(int j = 0; j<loopc; j++) ChiSquare2 += pow((y[j+2]-wcoeff*Wtheo[j])/yr[j+2],2);


                ChiSquare2 = ChiSquare2/(loopc-1); //chi2/NDF, where NDF = DOF - 1; 

                //The chi2 with a plus one for the uncertainty.
                float chim = (c_min120 * (loopc-1) + 1) / (loopc-1);
                if (cc == 0 ) {
                    if (chim <= ChiSquare2 ) {
                        chi_rp120 = ChiSquare2;
                        drp120 = d1;
                        cc++;
                    }else if (chim >= ChiSquare2 && d1 >= 50.000) {  //this if statement prevents from if chi square curve is flatten and give too big wrong uncertainty
                        chi_rp120 = ChiSquare2;
                        drp120 = d1;
                        cc++;
                    }
                }
            }
            for (d1 = delmin120; d1 >= -51.000; d1 += -0.010){
                ChiSquare2 = 0;
                float a2 = Q2 * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
    
                if(isIfin){
                    loopc = 0;
                    for(const auto& ang : I_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }

                if(isFipps)
                {
                    loopc = 0;
                    for(const auto& ang : F_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }
                
               // "Normalization Coefficent"
                wcoeff1 = 0; 
                wcoeff2 = 0;

                for(int i = 0; i<loopc; i++) wcoeff1 +=  (y[i+2]*Wtheo[i]/pow(yr[i+2],2));
                for(int i = 0; i<loopc; i++) wcoeff2 +=  pow(Wtheo[i]/yr[i+2],2);

                wcoeff = wcoeff1/wcoeff2;

                // "Calculation of Chi2/NDF value"
                for(int j = 0; j<loopc; j++) ChiSquare2 += pow((y[j+2]-wcoeff*Wtheo[j])/yr[j+2],2);


                ChiSquare2 = ChiSquare2/(loopc-1); //chi2/NDF, where NDF = DOF - 1; 

                float chim = (c_min120 * (loopc-1) + 1) /  (loopc-1);
                if (cc == 1 ) {
                    if (chim <= ChiSquare2 ) {
                        chi_rm120 = ChiSquare2;
                        drm120 = d1;
                        cc++;
                    }else if (chim >= ChiSquare2 && d1 <= -50.000) {  //this if statement prevents from if chi square curve is flatten and give too big wrong uncertainty
                        chi_rm120 = ChiSquare2;
                        drm120 = d1;
                        cc++;
                    }
                }
            }
        }else if (J1 == 2 && J2 == 2 && J3 == 0) {
            for (d1 = delmin220; d1 <= 52.000; d1 += 0.010){
                ChiSquare2 = 0;
                float a2 = Q2 * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
            
                if(isIfin){
                    loopc = 0;
                    for(const auto& ang : I_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }

                if(isFipps)
                {
                    loopc = 0;
                    for(const auto& ang : F_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }
                
               // "Normalization Coefficent"
                wcoeff1 = 0; 
                wcoeff2 = 0;

                for(int i = 0; i<loopc; i++) wcoeff1 +=  (y[i+2]*Wtheo[i]/pow(yr[i+2],2));
                for(int i = 0; i<loopc; i++) wcoeff2 +=  pow(Wtheo[i]/yr[i+2],2);

                wcoeff = wcoeff1/wcoeff2;

                // "Calculation of Chi2/NDF value"
                for(int j = 0; j<loopc; j++) ChiSquare2 += pow((y[j+2]-wcoeff*Wtheo[j])/yr[j+2],2);


                ChiSquare2 = ChiSquare2/(loopc-1); //chi2/NDF, where NDF = DOF - 1;
                float chim = (c_min220 * (loopc-1) + 1) /  (loopc-1);
                if (cc == 2 ) {
                    if (chim <= ChiSquare2 ) {
                        chi_rp220 = ChiSquare2;
                        drp220 = d1;
                        cc++;
                    }else if (chim >= ChiSquare2 && d1 >= 51.000) {  //this if statement prevents from if chi square curve is flatten and give too big wrong uncertainty
                        chi_rp220 = ChiSquare2;
                        drp220 = d1;
                        cc++;
                    }
                }
            }
            for (d1 = delmin220; d1 >= -52.00; d1 += -0.010){
                ChiSquare2 = 0;
                float a2 = Q2 * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));


                if(isIfin){
                    loopc = 0;
                    for(const auto& ang : I_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }

                if(isFipps)
                {
                    loopc = 0;
                    for(const auto& ang : F_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }
                
               // "Normalization Coefficent"
                wcoeff1 = 0; 
                wcoeff2 = 0;

                for(int i = 0; i<loopc; i++) wcoeff1 +=  (y[i+2]*Wtheo[i]/pow(yr[i+2],2));
                for(int i = 0; i<loopc; i++) wcoeff2 +=  pow(Wtheo[i]/yr[i+2],2);

                wcoeff = wcoeff1/wcoeff2;

                // "Calculation of Chi2/NDF value"
                for(int j = 0; j<loopc; j++) ChiSquare2 += pow((y[j+2]-wcoeff*Wtheo[j])/yr[j+2],2);


                ChiSquare2 = ChiSquare2/(loopc-1); //chi2/NDF, where NDF = DOF - 1;

                float chim = (c_min220 *  (loopc-1) + 1) / (loopc-1);
                if (cc == 3 ) {
                    if (chim <= ChiSquare2 ) {
                        chi_rm220 = ChiSquare2;
                        drm220 = d1;
                        cc++;
                    }else if (chim >= ChiSquare2 && d1 <= -50.000) {  //this if statement prevents from if chi square curve is flatten and give too big wrong uncertainty
                        chi_rm220 = ChiSquare2;
                        drm220 = d1;
                        cc++;
                    }
                }
            }
        }else if (J1 == 3 && J2 == 2 && J3 == 0) {
            for (d1 = delmin320; d1 <= 51.000; d1 += 0.010){
                ChiSquare2 = 0;
                float a2 = Q2 * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                
                if(isIfin){
                    loopc = 0;
                    for(const auto& ang : I_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }

                if(isFipps)
                {
                    loopc = 0;
                    for(const auto& ang : F_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }
                    
              // "Normalization Coefficent"
                wcoeff1 = 0; 
                wcoeff2 = 0;

                for(int i = 0; i<loopc; i++) wcoeff1 +=  (y[i+2]*Wtheo[i]/pow(yr[i+2],2));
                for(int i = 0; i<loopc; i++) wcoeff2 +=  pow(Wtheo[i]/yr[i+2],2);

                wcoeff = wcoeff1/wcoeff2;

                // "Calculation of Chi2/NDF value"
                for(int j = 0; j<loopc; j++) ChiSquare2 += pow((y[j+2]-wcoeff*Wtheo[j])/yr[j+2],2);


                ChiSquare2 = ChiSquare2/(loopc-1); //chi2/NDF, where NDF = DOF - 1;
              
                
                float chim = (c_min320 * (loopc-1) + 1) / (loopc-1);
                if (cc == 4 ) {
                    if (chim <= ChiSquare2 ) {
                        chi_rp320 = ChiSquare2;
                        drp320 = d1;
                        cc++;
                    }else if (chim >= ChiSquare2 && d1 >= 50.000) {  //this if statement prevents from if chi square curve is flatten and give too big wrong uncertainty
                        chi_rp320 = ChiSquare2;
                        drp320 = d1;
                        cc++;
                    }
                }
            }
            for (d1 = delmin320; d1 >= -51.000; d1 += -0.010){
                ChiSquare2 = 0;
                float a2 = Q2 * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                
               

                if(isIfin){
                    loopc = 0;
                    for(const auto& ang : I_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }

                if(isFipps)
                {
                    loopc = 0;
                    for(const auto& ang : F_angs){
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }
                
               // "Normalization Coefficent"
                wcoeff1 = 0; 
                wcoeff2 = 0;

                for(int i = 0; i<loopc; i++) wcoeff1 +=  (y[i+2]*Wtheo[i]/pow(yr[i+2],2));
                for(int i = 0; i<loopc; i++) wcoeff2 +=  pow(Wtheo[i]/yr[i+2],2);

                wcoeff = wcoeff1/wcoeff2;

                // "Calculation of Chi2/NDF value"
                for(int j = 0; j<loopc; j++) ChiSquare2 += pow((y[j+2]-wcoeff*Wtheo[j])/yr[j+2],2);


                ChiSquare2 = ChiSquare2/(loopc-1); //chi2/NDF, where NDF = DOF - 1;

                float chim = (c_min320 * (loopc-1) + 1) / (loopc-1);
                if (cc == 5 ) {
                    if (chim <= ChiSquare2 ) {
                        chi_rm320 = ChiSquare2;
                        drm320 = d1;
                        cc++;
                    }else if (chim >= ChiSquare2 && d1 <= -50.000) {  //this if statement prevents from if chi square curve is flatten and give too big wrong uncertainty
                        chi_rm320 = ChiSquare2;
                        drm320 = d1;
                        cc++;
                    }
                }
            }
            
        }
        // finishing to find uncetainty on mixing ratio based on 1 sigma uncertainty on chi^2
    
        /*
        //===================================================================================[4]
        std::ofstream chi0file (chi0name);
        if (chi0file.is_open())
        {
            
            
            
            
            if(isIfin){
                    loopc = 0;
                    for(const auto& ang : I_angs){
                        Wtheo[loopc] =  1 + a22 * (3 * pow(cos (ang), 2) - 1) / 2 + a44 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }

                if(isFipps)
                {
                    loopc = 0;
                    for(const auto& ang : F_angs){
                        Wtheo[loopc] =  1 + a22 * (3 * pow(cos (ang), 2) - 1) / 2 + a44 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }

                }
                wcoeff = 0; 
               // "Normalization Coefficent"
                for(int i = 0; i<=loopc; i++) wcoeff +=  (y[i+2]*Wtheo[i]/pow(yr[i+2],2));

                // "Calculation of Chi2/NDF value"
                for(int j = 0; j<=loopc; j++) ChiSquare_0 += pow((y[j+2]-wcoeff*Wtheo[j])/yr[j+2],2);


                ChiSquare_0 = ChiSquare_0/(loopc-1); //chi2/NDF, where NDF = DOF - 1;
                
                
            if (countchi0==0) {
                chi_min0 = ChiSquare_0;
                chi_max0 = ChiSquare_0;
            }else if (chi_min0 > ChiSquare_0){
                chi_min0 = ChiSquare_0;
            }else if(chi_max0 < ChiSquare_0){
                chi_max0 = ChiSquare_0;
            }
            countchi0++;
            
            if (J1 == 0 && J2 == 2 && J3 == 0) {
                if (countchi==0) {
                    chi_min = ChiSquare_0;
                    chi_max = ChiSquare_0;
                }else
                    if (chi_min > ChiSquare_0){
                        chi_min = ChiSquare_0;
                    }else if(chi_max < ChiSquare_0){
                        chi_max = ChiSquare_0;
                    }
                countchi++;
            }else if (J1 == 4 && J2 == 2 && J3 == 0) {
                if (countchi==0) {
                    chi_min = ChiSquare_0;
                    chi_max = ChiSquare_0;
                }else
                    if (chi_min > ChiSquare_0){
                        chi_min = ChiSquare_0;
                    }else if(chi_max < ChiSquare_0){
                        chi_max = ChiSquare_0;
                    }
                countchi++;
            }
            if (J1 == 0 && J2 == 2 && J3 == 0) {
                c_min020 = ChiSquare_0;
            }else if (J1 == 4 && J2 == 2 && J3 == 0) {
                c_min420 = ChiSquare_0;
            }
            chi0file << atan (0) * 180 / PI << "  " << ChiSquare_0 << std::endl;
           
        }//===================================================================================[4]
        chi0file.close();   
        */
        /*
        // finishing to find uncetainty on mixing ratio based on 1 sigma uncertainty on chi^2

        std::cout << J1 << "-" << J2 << "-" << J3 << std::endl;
        if (J1 == 1 && J2 == 2 && J3 == 0) {
            std::cout << "Minimum chi square value is " << c_min << std::endl;
            std::cout << "Delta value at minimum chi square is " << delmin120 << std::endl;
            std::cout << "Maximum chi square value is " << c_max << std::endl;
            std::cout << "Delta value at maximum chi square is " << delmax120 << std::endl;
            //std::cout << "Uncertiant of delta^2 is +" << pow(drp120, 2) << std::endl;
            //std::cout << "Minimum delta^2 value is " << pow(delmin120, 2) << std::endl;
            //std::cout << "Uncertiant of delta^2 is -" << pow(drm120, 2) << std::endl;
            //std::cout << "Uncertiant of delta is + (" << drp120 - delmin120 << ")" << std::endl;
            //std::cout << "Minimum delta value is (" << delmin120 << ")" << std::endl;
            //std::cout << "Uncertiant of delta is - (" << drm120 - delmin120 << ")" << std::endl;
        }else if (J1 == 2 && J2 == 2 && J3 == 0) {
            std::cout << "Minimum chi square value is " << c_min << std::endl;
            std::cout << "Delta value at minimum chi square is " << delmin220 << std::endl;
            std::cout << "Maximum chi square value is " << c_max << std::endl;
            //std::cout << "Delta value at maximum chi square is " << delmax220 << std::endl;
            //std::cout << "Uncertiant of delta^2 is +" << pow(drp220, 2) << std::endl;
            //std::cout << "Minimum delta^2 value is " << pow(delmin220, 2) << std::endl;
            //std::cout << "Uncertiant of delta^2 is -" << pow(drm220, 2) << std::endl;
            //std::cout << "Uncertiant of delta is + (" << drp220 - delmin220 << ")" << std::endl;
            //std::cout << "Minimum delta value is (" << delmin220 << ")" << std::endl;
            //std::cout << "Uncertiant of delta is - (" << trunc((drm220 - delmin220)*100)/100  << ")" << std::endl;
        }else if (J1 == 3 && J2 == 2 && J3 == 0) {
            std::cout << "Minimum chi square value is " << c_min << std::endl;
            std::cout << "Delta value at minimum chi square is " << delmin320 << std::endl;
            std::cout << "Maximum chi square value is " << c_max << std::endl;
            //std::cout << "Delta value at maximum chi square is " << delmax320 << std::endl;
            //std::cout << "Uncertiant of delta^2 is +" << pow(drp320, 2) << std::endl;
            //std::cout << "Minimum delta^2 value is " << pow(delmin320, 2) << std::endl;
            //std::cout << "Uncertiant of delta^2 is -" << pow(drm320, 2) << std::endl;
            //std::cout << "Uncertiant of delta is + (" << drp320 - delmin320 << ")" << std::endl;
            //std::cout << "Minimum delta value is (" << delmin320 << ")" << std::endl;
            //std::cout << "Uncertiant of delta is - (" << drm320 - delmin320 << ")" << std::endl;
        }else {
            std::cout << "Minimum chi square value at 0 mixing is " << chi_min0 << std::endl;
            std::cout << "Delta value at minimum chi square is " << 0 << std::endl;
            std::cout << "Maximum chi square value at 0 mixing is " << chi_max0 << std::endl;
            std::cout << "Delta value at maximum chi square is " << 0 << std::endl;
        }*/

    }//=========================================================================================================|J|
            
   return 0;    
 }    
    
