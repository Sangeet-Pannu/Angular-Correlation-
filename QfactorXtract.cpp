//
//  QfactorXtract.cxx
//  Program to extracting spins and Mixing ratios
//  of x-2-0 transition cascades.
//  ** Modification of Q_factor.cpp which was Created by Sambuu on 2019-07-10. **
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

    //Fipps: All the possible angles
    float F_angs[21] = {0.579857,0.615795,0.7812,0.809525,0.982443,1.00657,1.36004,
        1.38068,1.56068,1.58091,1.761,1.78155,2.13503,2.15915,2.33207,2.36039,2.5258,
        2.56174,2.85618,2.94012,3.14159};

    //Ifin: All the possible angles.
    float I_angs[29] = {0.795613, 0.909671,0.945734,0.995199,1.048,1.0742,1.081,1.15188,
        1.20067,1.29824,1.38323,1.39967,1.56271,1.57888,1.7419,1.75836,1.8433,1.94093,1.98971,
        2.06059,2.06739,2.09359,2.1463,2.19586,2.23192,2.34598,2.8865,2.96151,3.14159};

    float Wtheo[29];
    float Wtheo_0[29];

    int J2 = 2, J3 = 0;
   
    float chi_min0 = 0;
    float chi_max0 = 0;
    float ChiSquare_0 = 0;

    float exp_max = 0;

    float exp_min = 0;
   

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


    for(int num =2;num<(count1-1);num++){

        if(num==2){
            exp_max = y[num];
            exp_min = y[num];
            continue;
        }

        if(y[num]>exp_max){
            exp_max = y[num];
        }else if(y[num]<exp_min){
            exp_min = y[num];

        }
    

    }

    std::cout << exp_max << " " << exp_min << std::endl;

   

    
    //==========================================================================================================|J|
        int J1 = 0;

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
    
    double Q2,Q4;

    double Q2min,Q4min = 0;

    double Q2max,Q4max = 0;        


    for(Q4=1.00; Q4>= 0.10; Q4 -= 0.01)
    {    
        for(Q2=1.00; Q2>= 0.10; Q2 -= 0.01)
        { 

            float a22 = Q2 * R2LLJ2J1 * R2LLJ2J3;
            float a44 = Q4 * R4LLJ2J1 * R4LLJ2J3;


            std::stringstream chi0filename;
            chi0filename << "chi_square" << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Qxtract.dat";
            std::string chi0name = chi0filename.str();

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
                wcoeff1 = 0; 
                wcoeff2 = 0;

                for(int i = 0; i<loopc; i++) wcoeff1 +=  (y[i+2]*Wtheo[i]/pow(yr[i+2],2));
                for(int i = 0; i<loopc; i++) wcoeff2 +=  pow(Wtheo[i]/yr[i+2],2);

                wcoeff = wcoeff1/wcoeff2;

                // "Calculation of Chi2/NDF value"
   
                for(int j = 0; j<loopc; j++) ChiSquare_0 += pow((y[j+2]-wcoeff*Wtheo[j])/yr[j+2],2);

                ChiSquare_0 = ChiSquare_0/(loopc-1); //chi2/NDF, where NDF = DOF - 1;

                if(countchi0==0){
                    chi_min0 = ChiSquare_0;
                    chi_max0 = ChiSquare_0;

                    Q2min = Q2;
                    Q4min = Q4;
                }else if(chi_min0 > ChiSquare_0){
                    chi_min0 = ChiSquare_0;
                    Q2min = Q2;
                    Q4min = Q4;
                }else if(chi_max0 < ChiSquare_0){
                    chi_max0 = ChiSquare_0;
                    Q2max = Q2;
                    Q4max = Q4;
                }

                countchi0++;

                chi0file << atan(0)*180/PI << "  " << ChiSquare_0 << std::endl;
            }

            chi0file.close();
    
        }   
     
    }
   
    std::cout << "Minimum chi_square value is " << chi_min0 << " at Q2 = " << Q2min << " Q4 = " << Q4min << std::endl;
    std::cout << "Maximum chi_square value is " << chi_max0 << " at Q2 = " << Q2max << " Q4 = " << Q4max << std::endl;
        
   return 0;    
 }    
    
