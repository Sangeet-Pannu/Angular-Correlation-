//
//  Chi2_Mix_Xtract.cxx
//  Program to extracting spins and Mixing ratios
//  of x-2-0 transition cascades.
//  ** Modification of ChiRu100_x-2-0 which was Created by Sambuu on 2019-07-10. **
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
    int peak = 2445;

    //Fipps: All the possible angles
    float F_angs[21] = {0.579857,0.615795,0.7812,0.809525,0.982443,1.00657,1.36004,
        1.38068,1.56068,1.58091,1.761,1.78155,2.13503,2.15915,2.33207,2.36039,2.5258,
        2.56174,2.85618,2.94012,3.14159};

    //Ifin: All the possible angles.
    float I_angs[29] = {0.795613, 0.909671,0.945734,0.995199,1.048,1.0742,1.081,1.15188,
        1.20067,1.29824,1.38323,1.39967,1.56271,1.57888,1.7419,1.75836,1.8433,1.94093,1.98971,
        2.06059,2.06739,2.09359,2.1463,2.19586,2.23192,2.34598,2.8865,2.96151,3.14159};

    //Ifin: All the possible angles.
    float M_angs[56] = {0.5910522605,0.6227631477,0.7709206573,0.7806159612,0.7924842002,
        0.7965403453,0.8175750535,0.8953835769,0.9090739395,0.9305484706,0.9493613746,0.9722827837,
        0.9916594291,0.9940138782,1.033130197,1.063060849,1.074906398,1.079959127,1.082505562,1.149833383,
        1.166417502,1.204371432,1.215853953,1.313170021,1.338644847,1.468825465,1.494504495,1.518713957,
        1.622878697,1.647088159,1.672767188,1.802942571,1.828424378,1.925743937,1.937228203,1.975171661,
        1.991752289,2.059087092,2.061635272,2.066679274,2.078530059,2.108462456,2.147575285,2.149931479,
        2.169304634,2.19223826,2.211035456,2.232520459,2.246203841,2.324010619,2.34505929,2.349108453,2.360976692,
        2.370680723,2.51882427,2.550536902};


    float Wtheo[106];
    float Wtheo_0[106];

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
    float exp_min = 0;
    float exp_max = 0;

    //Fipps-Fipps
    float Q2F = 0.8833;// Has a 1.03% error
    float Q4F = 0.6834;// Has a 0.6% Error
    //Ifin-Ifin
    float Q2I = 0.9390;// Has a 3.91% Error
    float Q4I = 0.7359;//Has a 0.85% Error
    //Mixing
    float Q2M = 0.9006;//Has a 1.65% Error
    float Q4M = 0.7109;//Has a 0.72% Error

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

    //Calculated Confidence value (99%)
    float CV95 =1.34; 
    
    int loopc;

    std::ifstream infile("exp.txt");
    std::string lline;
    float y[106];
    float yr[106];
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

    std::ifstream infileI("exp.txt");
    std::string llineI;
  
    if (infileI.is_open())
        {
     
        while (std::getline(infileI, llineI))
        {

            std::stringstream ss(llineI);
            float b, c;
            int a;
            
                    
            if (ss >> a >> b >> c)
            {
                y[count1] = b;
                yr[count1] = c;
            }
            count1++;
        }
        
        infileI.close();
        
        }
    else
    {
        std::cout << "error: infile IFIN is not open" << std::endl;
    } 

    std::ifstream infileM("exp.txt");
    std::string llineM;
    if (infileM.is_open())
        {
     
        while (std::getline(infileM, llineM))
        {

            std::stringstream ss(llineM);
            float b, c;
            int a;
            
                    
            if (ss >> a >> b >> c)
            {
                y[count1] = b;
                yr[count1] = c;
            }
            count1++;
        }
        
        infileM.close();
        
        }
    else
    {
        std::cout << "error: infile Mixing is not open" << std::endl;
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

    std::stringstream cfilename1;
    std::stringstream cfilename2;
    std::stringstream cfilename3;
        /*
        cfilename << "exp.dat";
        std::string cname = cfilename.str();
        
        std::ofstream cfile (cname);
        if (cfile.is_open())
        {
    
                for(int i = 0; i<106; i++){
                    if(i<21)cfile << cos(F_angs[i]) << " " << y[i+2] << " " << yr[i+2] << std::endl;
                    if(i>=21 && i<50 )cfile << cos(I_angs[i]) << " " << y[i+2] << " " << yr[i+2] << std::endl;
                    if(i>=50)cfile << cos(M_angs[i]) << " " << y[i+2] << " " << yr[i+2] << std::endl;
                }
     
 
        }
        cfile.close();
*/
        cfilename1 << "exp_Fipps.dat";
        std::string cname1 = cfilename1.str();
        
        std::ofstream cfile1 (cname1);
        if (cfile1.is_open())
        {
    
                for(int i = 0; i<21; i++){
                    cfile1 << cos(F_angs[i]) << " " << y[i+2] << " " << yr[i+2] << std::endl;
                   
                }
     
 
        }
        cfile1.close();

        cfilename2 << "exp_Ifin.dat";
        std::string cname2 = cfilename2.str();
        
        std::ofstream cfile2 (cname2);
        if (cfile2.is_open())
        {
    
                for(int i = 21; i<50; i++){
                   
                    cfile2 << cos(I_angs[i-21]) << " " << y[i+2] << " " << yr[i+2] << std::endl;
                    
                }
     
 
        }
        cfile2.close();

        cfilename3 << "exp_Mix.dat";
        std::string cname3 = cfilename3.str();
        
        std::ofstream cfile3 (cname3);
        if (cfile3.is_open())
        {
    
                for(int i = 50; i<106; i++){
                 
                    cfile3 << cos(M_angs[i-50]) << " " << y[i+2] << " " << yr[i+2] << std::endl;
                }
     
 
        }
        cfile3.close();
    
    //==========================================================================================================|J|
    for (int J1 = 0; J1 < 5; J1++)
    {
    
        std::stringstream chi0filename;
        chi0filename << "chi_square" << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Ru100.dat";
        std::string chi0name = chi0filename.str();

        std::stringstream chifilename;
       
        chifilename << "chi_square" << "_"<< J1 << "-" << J2 << "-" << J3 << "_Ru100.dat";
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
        
        
        
       
        std::ofstream chifile (chiname);
        //==================================================================================|2|
        if (chifile.is_open())
        {
            double delta1;
            for (delta1 = -50.000; delta1 < 50.000; delta1 += 0.010)
            {
                //Set to zero as the second gamma goes from 2_1+ to 0_1+ which is a pure E2 transition.
                int delta2 = 0;
                
                float a2=0;
                float a4=0;
                loopc = 0;

                //Definition of the theoritcal angular correlation results.
                //---------------------------------------------------------**
                if(loopc<21)
                {  
                    for(const auto& ang : F_angs){
                        a2= Q2F * ((R2LLJ2J1 + 2 * delta1 * R2LMJ2J1 + pow(delta1, 2) * R2MMJ2J1) * (R2LLJ2J3 + 2 * delta2 * R2LMJ2J3 + pow(delta2, 2) * R2MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));
                        a4= Q4F * ((R4LLJ2J1 + 2 * delta1 * R4LMJ2J1 + pow(delta1, 2) * R4MMJ2J1) * (R4LLJ2J3 + 2 * delta2 * R4LMJ2J3 + pow(delta2, 2) * R4MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else if(loopc>=21 && loopc<50)
                {
                    for(const auto& ang : I_angs){
                        a2= Q2I * ((R2LLJ2J1 + 2 * delta1 * R2LMJ2J1 + pow(delta1, 2) * R2MMJ2J1) * (R2LLJ2J3 + 2 * delta2 * R2LMJ2J3 + pow(delta2, 2) * R2MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));
                        a4= Q4I * ((R4LLJ2J1 + 2 * delta1 * R4LMJ2J1 + pow(delta1, 2) * R4MMJ2J1) * (R4LLJ2J3 + 2 * delta2 * R4LMJ2J3 + pow(delta2, 2) * R4MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else
                {
                    for(const auto& ang : M_angs){
                         a2= Q2M * ((R2LLJ2J1 + 2 * delta1 * R2LMJ2J1 + pow(delta1, 2) * R2MMJ2J1) * (R2LLJ2J3 + 2 * delta2 * R2LMJ2J3 + pow(delta2, 2) * R2MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));
                        a4= Q4M * ((R4LLJ2J1 + 2 * delta1 * R4LMJ2J1 + pow(delta1, 2) * R4MMJ2J1) * (R4LLJ2J3 + 2 * delta2 * R4LMJ2J3 + pow(delta2, 2) * R4MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }
                
                //========================================================**
            
        
           

             
                //=============================================================================***
            
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
                chifile << atan (delta1)* 180 / PI << "  " << ChiSquare << std::endl;
                
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
                
                float a2=0;
                float a4=0;
                loopc = 0;

                //Definition of the theoritcal angular correlation results.
                //---------------------------------------------------------**
                if(loopc<21)
                {  
                    for(const auto& ang : F_angs){
                        a2 = Q2F * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4F * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else if(loopc>=21 && loopc<50)
                {
                    for(const auto& ang : I_angs){
                        a2 = Q2I * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4I * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else
                {
                    for(const auto& ang : M_angs){
                        a2 = Q2M * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4M * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
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
                float a2=0;
                float a4=0;
                loopc = 0;

                //Definition of the theoritcal angular correlation results.
                //---------------------------------------------------------**
                if(loopc<21)
                {  
                    for(const auto& ang : F_angs){
                        a2 = Q2F * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4F * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else if(loopc>=21 && loopc<50)
                {
                    for(const auto& ang : I_angs){
                        a2 = Q2I * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4I * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else
                {
                    for(const auto& ang : M_angs){
                        a2 = Q2M * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4M * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
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
                float a2=0;
                float a4=0;
                loopc = 0;

                //Definition of the theoritcal angular correlation results.
                //---------------------------------------------------------**
                if(loopc<21)
                {  
                    for(const auto& ang : F_angs){
                        a2 = Q2F * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4F * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else if(loopc>=21 && loopc<50)
                {
                    for(const auto& ang : I_angs){
                        a2 = Q2I * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4I * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else
                {
                    for(const auto& ang : M_angs){
                        a2 = Q2M * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4M * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
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
                float a2=0;
                float a4=0;
                loopc = 0;

                //Definition of the theoritcal angular correlation results.
                //---------------------------------------------------------**
                if(loopc<21)
                {  
                    for(const auto& ang : F_angs){
                        a2 = Q2F * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4F * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else if(loopc>=21 && loopc<50)
                {
                    for(const auto& ang : I_angs){
                        a2 = Q2I * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4I * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else
                {
                    for(const auto& ang : M_angs){
                        a2 = Q2M * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4M * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
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
                float a2=0;
                float a4=0;
                loopc = 0;

                //Definition of the theoritcal angular correlation results.
                //---------------------------------------------------------**
                if(loopc<21)
                {  
                    for(const auto& ang : F_angs){
                        a2 = Q2F * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4F * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else if(loopc>=21 && loopc<50)
                {
                    for(const auto& ang : I_angs){
                        a2 = Q2I * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4I * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else
                {
                    for(const auto& ang : M_angs){
                        a2 = Q2M * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4M * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
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
                float a2=0;
                float a4=0;
                loopc = 0;

                //Definition of the theoritcal angular correlation results.
                //---------------------------------------------------------**
                if(loopc<21)
                {  
                    for(const auto& ang : F_angs){
                        a2 = Q2F * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4F * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else if(loopc>=21 && loopc<50)
                {
                    for(const auto& ang : I_angs){
                        a2 = Q2I * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4I * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else
                {
                    for(const auto& ang : M_angs){
                        a2 = Q2M * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                        a4 = Q4M * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
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
    
        
        //===================================================================================[4]
        std::ofstream chi0file (chi0name);
        if (chi0file.is_open())
        {
                float a22 = 0;
                float a44 = 0;
                loopc = 0;

                //Definition of the theoritcal angular correlation results.
                //---------------------------------------------------------**
                if(loopc<21)
                {  
                    for(const auto& ang : F_angs){
                        a22=Q2F * R2LLJ2J1 * R2LLJ2J3;
                        a44=Q4F * R4LLJ2J1 * R4LLJ2J3;
                        Wtheo_0[loopc] =  1 + a22 * (3 * pow(cos (ang), 2) - 1) / 2 + a44 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else if(loopc>=21 && loopc<50)
                {
                    for(const auto& ang : I_angs){
                        a22=Q2I * R2LLJ2J1 * R2LLJ2J3;
                        a44=Q4I * R4LLJ2J1 * R4LLJ2J3;
                        Wtheo_0[loopc] =  1 + a22 * (3 * pow(cos (ang), 2) - 1) / 2 + a44 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else
                {
                    for(const auto& ang : M_angs){
                        a22=Q2M * R2LLJ2J1 * R2LLJ2J3;
                        a44=Q4M * R4LLJ2J1 * R4LLJ2J3;
                        Wtheo_0[loopc] =  1 + a22 * (3 * pow(cos (ang), 2) - 1) / 2 + a44 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }
                wcoeff1 = 0; 
                wcoeff2 = 0;

                for(int i = 0; i<loopc; i++) wcoeff1 +=  (y[i+2]*Wtheo_0[i]/pow(yr[i+2],2));
                for(int i = 0; i<loopc; i++) wcoeff2 +=  pow(Wtheo_0[i]/yr[i+2],2);

                wcoeff = wcoeff1/wcoeff2;

                // "Calculation of Chi2/NDF value"
   
                for(int j = 0; j<loopc; j++) ChiSquare_0 += pow((y[j+2]-wcoeff*Wtheo_0[j])/yr[j+2],2);

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
        
        
        // finishing to find uncetainty on mixing ratio based on 1 sigma uncertainty on chi^2

        std::cout << J1 << "-" << J2 << "-" << J3 << std::endl;
        if (J1 == 1 && J2 == 2 && J3 == 0) {
            std::cout << "Minimum chi square value is " << c_min << std::endl;
            std::cout << "Delta value at minimum chi square is " << delmin120 << std::endl;
            std::cout << "Maximum chi square value is " << c_max << std::endl;
            std::cout << "Delta value at maximum chi square is " << delmax120 << std::endl;
            std::cout << "Uncertiant of delta^2 is +" << pow(drp120, 2) << std::endl;
            std::cout << "Minimum delta^2 value is " << pow(delmin120, 2) << std::endl;
            std::cout << "Uncertiant of delta^2 is -" << pow(drm120, 2) << std::endl;
            std::cout << "Uncertiant of delta is + (" << drp120 - delmin120 << ")" << std::endl;
            std::cout << "Minimum delta value is (" << delmin120 << ")" << std::endl;
            std::cout << "Uncertiant of delta is - (" << drm120 - delmin120 << ")" << std::endl;
        }else if (J1 == 2 && J2 == 2 && J3 == 0) {
            std::cout << "Minimum chi square value is " << c_min << std::endl;
            std::cout << "Delta value at minimum chi square is " << delmin220 << std::endl;
            std::cout << "Maximum chi square value is " << c_max << std::endl;
            std::cout << "Delta value at maximum chi square is " << delmax220 << std::endl;
            std::cout << "Uncertiant of delta^2 is +" << pow(drp220, 2) << std::endl;
            std::cout << "Minimum delta^2 value is " << pow(delmin220, 2) << std::endl;
            std::cout << "Uncertiant of delta^2 is -" << pow(drm220, 2) << std::endl;
            std::cout << "Uncertiant of delta is + (" << drp220 - delmin220 << ")" << std::endl;
            std::cout << "Minimum delta value is (" << delmin220 << ")" << std::endl;
            std::cout << "Uncertiant of delta is - (" << trunc((drm220 - delmin220)*100)/100  << ")" << std::endl;
        }else if (J1 == 3 && J2 == 2 && J3 == 0) {
            std::cout << "Minimum chi square value is " << c_min << std::endl;
            std::cout << "Delta value at minimum chi square is " << delmin320 << std::endl;
            std::cout << "Maximum chi square value is " << c_max << std::endl;
            std::cout << "Delta value at maximum chi square is " << delmax320 << std::endl;
            std::cout << "Uncertiant of delta^2 is +" << pow(drp320, 2) << std::endl;
            std::cout << "Minimum delta^2 value is " << pow(delmin320, 2) << std::endl;
            std::cout << "Uncertiant of delta^2 is -" << pow(drm320, 2) << std::endl;
            std::cout << "Uncertiant of delta is + (" << drp320 - delmin320 << ")" << std::endl;
            std::cout << "Minimum delta value is (" << delmin320 << ")" << std::endl;
            std::cout << "Uncertiant of delta is - (" << drm320 - delmin320 << ")" << std::endl;
        }else {
            std::cout << "Minimum chi square value at 0 mixing is " << chi_min0 << std::endl;
            std::cout << "Delta value at minimum chi square is " << 0 << std::endl;
            std::cout << "Maximum chi square value at 0 mixing is " << chi_max0 << std::endl;
            std::cout << "Delta value at maximum chi square is " << 0 << std::endl;
       }

         std::stringstream conflev95filename;
        conflev95filename << "confidencelev95_Ru100.dat";
        std::string conflevname = conflev95filename.str();
        
        std::ofstream cvfile (conflevname);
        if (cvfile.is_open())
        {
            double delta1;
            for (delta1 = -50.000; delta1 < 50.000; delta1 += 0.010){
                cvfile << atan (delta1) * 180 / PI << "  " << CV95 << std::endl;
            }
        }
    
        /*
         beyond this point is all about to draw expected AC at min mixing ratio
         */
        
        int delta2 = 0;
        float a22 = 0;
        float a44 = 0;
        float a2 = 0;
        float a4 =0;
        
        loopc = 0;

                //Definition of the theoritcal angular correlation results.
                //---------------------------------------------------------**
                if(loopc<21)
                {  
                    for(const auto& ang : F_angs){
                        a2 = Q2F * ((R2LLJ2J1 + 2 * del_min * R2LMJ2J1 + pow(del_min, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(del_min, 2)));
                        a4 = Q4F * ((R4LLJ2J1 + 2 * del_min * R4LMJ2J1 + pow(del_min, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(del_min, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else if(loopc>=21 && loopc<50)
                {
                    for(const auto& ang : I_angs){
                        a2 = Q2I * ((R2LLJ2J1 + 2 * del_min * R2LMJ2J1 + pow(del_min, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(del_min, 2)));
                        a4 = Q4I * ((R4LLJ2J1 + 2 * del_min * R4LMJ2J1 + pow(del_min, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(del_min, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }else
                {
                    for(const auto& ang : M_angs){
                        a2 = Q2M * ((R2LLJ2J1 + 2 * del_min * R2LMJ2J1 + pow(del_min, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(del_min, 2)));
                        a4 = Q4M * ((R4LLJ2J1 + 2 * del_min * R4LMJ2J1 + pow(del_min, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(del_min, 2)));
                        Wtheo[loopc] =  1 + a2 * (3 * pow(cos (ang), 2) - 1) / 2 + a4 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc++;
                    }
                }

                int loopc1 = 0;

                //Definition of the theoritcal angular correlation results.
                //---------------------------------------------------------**
                if(loopc1<21)
                {  
                    for(const auto& ang : F_angs){
                        a22=Q2F * R2LLJ2J1 * R2LLJ2J3;
                        a44=Q4F * R4LLJ2J1 * R4LLJ2J3;
                        Wtheo_0[loopc1] =  1 + a22 * (3 * pow(cos (ang), 2) - 1) / 2 + a44 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc1++;
                    }
                }else if(loopc1>=21 && loopc1<50)
                {
                    for(const auto& ang : I_angs){
                        a22=Q2I * R2LLJ2J1 * R2LLJ2J3;
                        a44=Q4I * R4LLJ2J1 * R4LLJ2J3;
                        Wtheo_0[loopc1] =  1 + a22 * (3 * pow(cos (ang), 2) - 1) / 2 + a44 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc1++;
                    }
                }else
                {
                    for(const auto& ang : M_angs){
                        a22=Q2M * R2LLJ2J1 * R2LLJ2J3;
                        a44=Q4M * R4LLJ2J1 * R4LLJ2J3;
                        Wtheo_0[loopc1] =  1 + a22 * (3 * pow(cos (ang), 2) - 1) / 2 + a44 * (35 * pow(cos (ang), 4) - 30 * pow(cos (ang), 2) + 3) / 8;
                        loopc1++;
                    }
                }
                  
        
                double wcoeff1F = 0; 
                double wcoeff2F = 0;

                double wcoeff1I = 0; 
                double wcoeff2I = 0;

                double wcoeff1M = 0; 
                double wcoeff2M = 0;

                double wcoeffF = 0;
                double wcoeffI = 0;
                double wcoeffM = 0;

                for(int i = 0; i<21; i++) wcoeff1F +=  (y[i+2]*Wtheo[i]/pow(yr[i+2],2));
                for(int i = 0; i<21; i++) wcoeff2F +=  pow(Wtheo[i]/yr[i+2],2);

                for(int i = 21; i<50; i++) wcoeff1I +=  (y[i+2]*Wtheo[i]/pow(yr[i+2],2));
                for(int i = 21; i<50; i++) wcoeff2I +=  pow(Wtheo[i]/yr[i+2],2);

                for(int i = 50; i<loopc; i++) wcoeff1M +=  (y[i+2]*Wtheo[i]/pow(yr[i+2],2));
                for(int i = 50; i<loopc; i++) wcoeff2M +=  pow(Wtheo[i]/yr[i+2],2);

                wcoeffF = wcoeff1F/wcoeff2F;
                wcoeffI = wcoeff1I/wcoeff2I;
                wcoeffM = wcoeff1M/wcoeff2M;

                double wcoeff1F_0 = 0; 
                double wcoeff2F_0 = 0;

                double wcoeff1I_0 = 0; 
                double wcoeff2I_0 = 0;

                double wcoeff1M_0 = 0; 
                double wcoeff2M_0 = 0;

                double wcoeffF_0 = 0;
                double wcoeffI_0 = 0;
                double wcoeffM_0 = 0;

                for(int i = 0; i<21; i++) wcoeff1F_0 +=  (y[i+2]*Wtheo_0[i]/pow(yr[i+2],2));
                for(int i = 0; i<21; i++) wcoeff2F_0 +=  pow(Wtheo_0[i]/yr[i+2],2);

                for(int i = 21; i<50; i++) wcoeff1I_0 +=  (y[i+2]*Wtheo_0[i]/pow(yr[i+2],2));
                for(int i = 21; i<50; i++) wcoeff2I_0 +=  pow(Wtheo_0[i]/yr[i+2],2);

                for(int i = 50; i<loopc; i++) wcoeff1M_0 +=  (y[i+2]*Wtheo_0[i]/pow(yr[i+2],2));
                for(int i = 50; i<loopc; i++) wcoeff2M_0 +=  pow(Wtheo_0[i]/yr[i+2],2);
                

                wcoeffF_0 = wcoeff1F_0/wcoeff2F_0;
                wcoeffI_0 = wcoeff1I_0/wcoeff2I_0;
                wcoeffM_0 = wcoeff1M_0/wcoeff2M_0;
        
        //--------------------------------------------------------------------------------------------------------------------------------------------------------------------|
        
        //Fipps Distribution
        std::stringstream mixedfilename;
        //mixedfilename << Output2 << "AC" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_min_Ru100.dat";
        mixedfilename << "AC_Fipps" << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_min_Ru100.dat";
        std::string mixedname = mixedfilename.str();
        std::ofstream mixedfile (mixedname);
        
        std::stringstream zerofilename;
        //zerofilename << Output2 << "AC" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Ru100.dat";
        zerofilename << "AC_Fipps" << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Ru100.dat";
        std::string zeroname = zerofilename.str();
        std::ofstream zerofile (zeroname);
        //-------------------------------|

        //Ifin Distribution
        std::stringstream mixedfilenameI;
        //mixedfilename << Output2 << "AC" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_min_Ru100.dat";
        mixedfilenameI << "AC_Ifin" << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_min_Ru100.dat";
        std::string mixednameI = mixedfilenameI.str();
        std::ofstream mixedfileI (mixednameI);
        
        std::stringstream zerofilenameI;
        //zerofilename << Output2 << "AC" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Ru100.dat";
        zerofilenameI << "AC_Ifin" << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Ru100.dat";
        std::string zeronameI = zerofilenameI.str();
        std::ofstream zerofileI (zeronameI);
        //-------------------------------|


        //Mixing Distribution
        std::stringstream mixedfilenameM;
        //mixedfilename << Output2 << "AC" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_min_Ru100.dat";
        mixedfilenameM << "AC_Mix" << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_min_Ru100.dat";
        std::string mixednameM = mixedfilenameM.str();
        std::ofstream mixedfileM (mixednameM);
        
        std::stringstream zerofilenameM;
        //zerofilename << Output2 << "AC" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Ru100.dat";
        zerofilenameM << "AC_Mix" << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Ru100.dat";
        std::string zeronameM = zerofilenameM.str();
        std::ofstream zerofileM (zeronameM);

        //--------------------------------------------------------------------------------------------------------------------------------------------------------------------|

        float Wtheomixed;
        float Wtheomixed_0;

        //===========================================================================|FIPPS AC SAVING|

        a2 = Q2F * ((R2LLJ2J1 + 2 * del_min * R2LMJ2J1 + pow(del_min, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(del_min, 2)));
        a4 = Q4F * ((R4LLJ2J1 + 2 * del_min * R4LMJ2J1 + pow(del_min, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(del_min, 2)));
        a22=Q2F * R2LLJ2J1 * R2LLJ2J3;
        a44=Q4F * R4LLJ2J1 * R4LLJ2J3;

        if (zerofile.is_open())
        {
            if (mixedfile.is_open())
            {
                for (int tt = 0; tt < 181; tt++)
                {
                    Wtheomixed = (1 + a2 * (3 * pow(cos (tt*PI/180), 2) - 1) / 2 + a4 * (35 * pow(cos (tt*PI/180), 4) - 30 * pow(cos (tt*PI/180), 2) + 3) / 8) * wcoeffF;
                    Wtheomixed_0 = (1 + a22 * (3 * pow(cos (tt*PI/180), 2) - 1) / 2 + a44 * (35 * pow(cos (tt*PI/180), 4) - 30 * pow(cos (tt*PI/180), 2) + 3) / 8) * wcoeffF_0;
                    
                    if (count==0) {
                        if (Wtheomixed > Wtheomixed_0){
                            ac_min = Wtheomixed_0;
                            ac_max = Wtheomixed;
                        }else{
                            ac_min = Wtheomixed;
                            ac_max = Wtheomixed_0;
                        }
                    }else if (ac_min > Wtheomixed_0){
                        ac_min = Wtheomixed_0;
                    }else if(ac_min > Wtheomixed){
                        ac_min = Wtheomixed;
                        
                    }else if(ac_max < Wtheomixed_0){
                        ac_max = Wtheomixed_0;
                    }else if(ac_max < Wtheomixed){
                        ac_max = Wtheomixed;
                        
                    }
                    count++;
                    
                    zerofile << cos(tt*PI/180) << " " << Wtheomixed_0 << std::endl;
                    mixedfile << cos(tt*PI/180) << " " << Wtheomixed << std::endl;
                }
            }
            zerofile.close();
        }
        mixedfile.close();
        //===========================================================================|FIPPS AC SAVING|
        
        
        //===========================================================================|IFIN AC SAVING|

        a2 = Q2I * ((R2LLJ2J1 + 2 * del_min * R2LMJ2J1 + pow(del_min, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(del_min, 2)));
        a4 = Q4I * ((R4LLJ2J1 + 2 * del_min * R4LMJ2J1 + pow(del_min, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(del_min, 2)));
        a22=Q2I * R2LLJ2J1 * R2LLJ2J3;
        a44=Q4I * R4LLJ2J1 * R4LLJ2J3;

        if (zerofileI.is_open())
        {
            if (mixedfileI.is_open())
            {
                for (int tt = 0; tt < 181; tt++)
                {
                    Wtheomixed = (1 + a2 * (3 * pow(cos (tt*PI/180), 2) - 1) / 2 + a4 * (35 * pow(cos (tt*PI/180), 4) - 30 * pow(cos (tt*PI/180), 2) + 3) / 8) * wcoeffI;
                    Wtheomixed_0 = (1 + a22 * (3 * pow(cos (tt*PI/180), 2) - 1) / 2 + a44 * (35 * pow(cos (tt*PI/180), 4) - 30 * pow(cos (tt*PI/180), 2) + 3) / 8) * wcoeffI_0;
                    
                    if (count==0) {
                        if (Wtheomixed > Wtheomixed_0){
                            ac_min = Wtheomixed_0;
                            ac_max = Wtheomixed;
                        }else{
                            ac_min = Wtheomixed;
                            ac_max = Wtheomixed_0;
                        }
                    }else if (ac_min > Wtheomixed_0){
                        ac_min = Wtheomixed_0;
                    }else if(ac_min > Wtheomixed){
                        ac_min = Wtheomixed;
                        
                    }else if(ac_max < Wtheomixed_0){
                        ac_max = Wtheomixed_0;
                    }else if(ac_max < Wtheomixed){
                        ac_max = Wtheomixed;
                        
                    }
                    count++;
                    
                    zerofileI << cos(tt*PI/180) << " " << Wtheomixed_0 << std::endl;
                    mixedfileI << cos(tt*PI/180) << " " << Wtheomixed << std::endl;
                }
            }
            zerofileI.close();
        }
        mixedfileI.close();
        //===========================================================================|IFIN AC SAVING|
       
        //===========================================================================|Mixed AC SAVING|
        a2 = Q2M * ((R2LLJ2J1 + 2 * del_min * R2LMJ2J1 + pow(del_min, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(del_min, 2)));
        a4 = Q4M * ((R4LLJ2J1 + 2 * del_min * R4LMJ2J1 + pow(del_min, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(del_min, 2)));
        a22=Q2M * R2LLJ2J1 * R2LLJ2J3;
        a44=Q4M * R4LLJ2J1 * R4LLJ2J3;

        if (zerofileM.is_open())
        {
            if (mixedfileM.is_open())
            {
                for (int tt = 0; tt < 181; tt++)
                {
                    Wtheomixed = (1 + a2 * (3 * pow(cos (tt*PI/180), 2) - 1) / 2 + a4 * (35 * pow(cos (tt*PI/180), 4) - 30 * pow(cos (tt*PI/180), 2) + 3) / 8) * wcoeffM;
                    Wtheomixed_0 = (1 + a22 * (3 * pow(cos (tt*PI/180), 2) - 1) / 2 + a44 * (35 * pow(cos (tt*PI/180), 4) - 30 * pow(cos (tt*PI/180), 2) + 3) / 8) * wcoeffM_0;
                    
                    if (count==0) {
                        if (Wtheomixed > Wtheomixed_0){
                            ac_min = Wtheomixed_0;
                            ac_max = Wtheomixed;
                        }else{
                            ac_min = Wtheomixed;
                            ac_max = Wtheomixed_0;
                        }
                    }else if (ac_min > Wtheomixed_0){
                        ac_min = Wtheomixed_0;
                    }else if(ac_min > Wtheomixed){
                        ac_min = Wtheomixed;
                        
                    }else if(ac_max < Wtheomixed_0){
                        ac_max = Wtheomixed_0;
                    }else if(ac_max < Wtheomixed){
                        ac_max = Wtheomixed;
                        
                    }
                    count++;
                    
                    zerofileM << cos(tt*PI/180) << " " << Wtheomixed_0 << std::endl;
                    mixedfileM << cos(tt*PI/180) << " " << Wtheomixed << std::endl;
                }
            }
            zerofileM.close();
        }
        mixedfileM.close();
        //===========================================================================|Mixed AC SAVING|
    }//=========================================================================================================|J|

    float ACmax = 0, ACmin =0;
    if (exp_max > ac_max) {
        ACmax = exp_max;
    }else{
        ACmax = ac_max;
    }
    if (exp_min < ac_min) {
        ACmin = exp_min;
    }else{
        ACmin = ac_min;
    }
    
    float ymin=0;
    float ymax=0;
    
    if (peak == 1073){
        ymin= 0.75 * ACmin;
        ymax= 1.12 * ACmax;
    }else{
        ymin= 0.75 * ACmin;
        ymax= 1.07 * ACmax;
    }
    

    //---------------------------------------------------------------------------------------------------------------------------|
    std::stringstream bfile;
    //bfile << Output2 << "ACscript" << peak << ".bfile"; //script file includes xmgrace Angular Correlation plots
    bfile << "ACscript_Fipps" << ".bfile"; //script file includes xmgrace Angular Correlation plots
    std::string bfilename = bfile.str();
    std::ofstream bfilebash(bfilename);
    
    if (bfilebash.is_open())
    {

        bfilebash << "read xy \"AC_Fipps" << "_0-2-0_delta_0_Ru100.dat\"" << std::endl;
        bfilebash << "read xy \"AC_Fipps" << "_1-2-0_delta_min_Ru100.dat\"" << std::endl;
        bfilebash << "read xy \"AC_Fipps" << "_2-2-0_delta_min_Ru100.dat\"" << std::endl;
        bfilebash << "read xy \"AC_Fipps" << "_3-2-0_delta_min_Ru100.dat\"" << std::endl;
        bfilebash << "read xy \"AC_Fipps" << "_4-2-0_delta_0_Ru100.dat\"" << std::endl;
        bfilebash << "read xydy \"" << "exp_Fipps.dat\"" << std::endl;
        
        bfilebash << "s0 line color 1" << std::endl;
        bfilebash << "s0 linestyle 4" << std::endl;
        bfilebash << "s0 linewidth 2" << std::endl;
        bfilebash << "s1 line color 2" << std::endl;
        bfilebash << "s1 linewidth 2" << std::endl;
        bfilebash << "s2 line color 3" << std::endl;
        bfilebash << "s2 linewidth 2" << std::endl;
        bfilebash << "s3 line color 4" << std::endl;
        bfilebash << "s3 linewidth 2" << std::endl;
        bfilebash << "s4 line color 1" << std::endl;
        bfilebash << "s4 linewidth 2" << std::endl;
        
        bfilebash << "s5 symbol 1" << std::endl;
        bfilebash << "s5 linestyle 0" << std::endl;
        bfilebash << "s5 symbol color 1" << std::endl;
        bfilebash << "s5 symbol size 1.2" << std::endl;
        bfilebash << "s5 symbol linewidth 1.5" << std::endl;
        bfilebash << "s5 errorbar linewidth 1.5" << std::endl;
        bfilebash << "s5 errorbar color 1" << std::endl;
        
        bfilebash << "s0 legend \"0-2-0 (\\xd\\0=0)\"" << std::endl;

        if (c_min120 <= chi_min) {
            bfilebash << "s1 legend \"1-2-0 (\\xd\\0=" << delmin120 << "\\S+" << drp120-delmin120 << "\\N" << "\\s" << round((drm120-delmin120) * pow(10, 2))/pow(10, 2) << "\\N)\"" << std::endl;
        }else{
            bfilebash << "s1 legend \"1-2-0 (\\xd\\0=" << delmin120 << ")\"" << std::endl;
        }
        if (c_min220 <= chi_min) {
            bfilebash << "s2 legend \"2-2-0 (\\xd\\0=" << delmin220 << "\\S+" << drp220-delmin220 << "\\N" << "\\s" << round((drm220-delmin220) * 100) / 100 << "\\N)\"" << std::endl;
        }else{
            bfilebash << "s2 legend \"2-2-0 (\\xd\\0=" << delmin220 << ")\"" << std::endl;
        }
        if (c_min320 <= chi_min) {
            bfilebash << "s3 legend \"3-2-0 (\\xd\\0=" << round(delmin320 * pow(10, 2))/pow(10, 2) << "\\S+" << round((drp320-delmin320)*100)/100 << "\\N" << "\\s" << round((drm320-delmin320)*100)/100 << "\\N)\"" << std::endl;
        }else{
            bfilebash << "s3 legend \"3-2-0 (\\xd\\0=" << round(delmin320 * pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        }
        bfilebash << "s4 legend \"4-2-0 (\\xd\\0=0)\"" << std::endl;

        bfilebash << "s5 legend \"exp " << peak << " keV\"" << std::endl; //EXPERIMENT PEAK keV

        

        //------------LEGEND CONFIGURATIONS---------------|
        bfilebash << "legend 0.49, 0.77" << std::endl;
        bfilebash << "legend char size 1.2" << std::endl;

        //------------------------------------------------|
        
        
        bfilebash << "with string" << std::endl;
        bfilebash << "string on" << std::endl;
        bfilebash << "string loctype view" << std::endl;
        bfilebash << "string 0.89, 0.65" << std::endl;
        bfilebash << "string char size 1.6" << std::endl;
        bfilebash << "string def \"(1-3)\\S+\\N(0\\S+\\N)\"" << std::endl;
        
        bfilebash << "with string" << std::endl;
        bfilebash << "string on" << std::endl;
        bfilebash << "string loctype view" << std::endl;
        bfilebash << "string 0.48, 0.44" << std::endl;
        bfilebash << "string char size 1.3" << std::endl;
        bfilebash << "string def \"gated on 658 keV \\xg\\0 ray (Ru-100)\"" << std::endl;
       
        bfilebash << "xaxis label \"Cos(\\xq\\0) \"" << std::endl;
        bfilebash << "xaxis label char size 1.6" << std::endl;
        bfilebash << "xaxis ticklabel char size 1.6" << std::endl;
        bfilebash << "xaxis tick major 0.5" << std::endl;
        bfilebash << "xaxis tick minor 0.25" << std::endl;
        bfilebash << "world xmin -1.2" << std::endl;
        bfilebash << "world xmax 1.2" << std::endl;
        
        bfilebash << "yaxis label \"Counts\"" << std::endl;
        bfilebash << "yaxis label char size 1.6" << std::endl;
        bfilebash << "yaxis ticklabel char size 1.6" << std::endl;
        
        bfilebash << "world ymin " << ymin << std::endl;
        bfilebash << "world ymax " << ymax << std::endl;
        
        if (exp_max >= 2 && exp_max < 10){
            bfilebash << "yaxis tick major 1" << std::endl;
            bfilebash << "yaxis tick minor 0.5" << std::endl;
            bfilebash << "view 0.11, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max >= 0 && exp_max < 2){
            bfilebash << "yaxis tick major 0.5" << std::endl;
            bfilebash << "yaxis tick minor 0.1" << std::endl;
            bfilebash << "view 0.11, 0.13, 1.26, 0.8" << std::endl;
        }else{
            std::cout << "Statistic for angular correlation analysis is too low or too high" << std::endl;
            return 0;
        }
        bfilebash << "saveall \"AC_Fipps"<< ".arg\"" << std::endl;
        bfilebash << "print to \"AC_Fipps"<< ".eps\"" << std::endl;
        bfilebash << "device \"EPS\" OP \"level2\"" << std::endl;
        bfilebash << "print" << std::endl;
        
        bfilebash.close();
    }
    else
    {
        std::cout << "error: the Fipps file is not open";
    }

    std::stringstream bfileI;
    //bfile << Output2 << "ACscript" << peak << ".bfile"; //script file includes xmgrace Angular Correlation plots
    bfileI << "ACscript_Ifin" << ".bfile"; //script file includes xmgrace Angular Correlation plots
    std::string bfilenameI = bfileI.str();
    std::ofstream bfilebashI(bfilenameI);
    
    if (bfilebashI.is_open())
    {

        bfilebashI << "read xy \"AC_Ifin" << "_0-2-0_delta_0_Ru100.dat\"" << std::endl;
        bfilebashI << "read xy \"AC_Ifin" << "_1-2-0_delta_min_Ru100.dat\"" << std::endl;
        bfilebashI << "read xy \"AC_Ifin" << "_2-2-0_delta_min_Ru100.dat\"" << std::endl;
        bfilebashI << "read xy \"AC_Ifin" << "_3-2-0_delta_min_Ru100.dat\"" << std::endl;
        bfilebashI << "read xy \"AC_Ifin" << "_4-2-0_delta_0_Ru100.dat\"" << std::endl;
        bfilebashI << "read xydy \"" << "exp_Ifin.dat\"" << std::endl;
        
        bfilebashI << "s0 line color 1" << std::endl;
        bfilebashI << "s0 linestyle 4" << std::endl;
        bfilebashI << "s0 linewidth 2" << std::endl;
        bfilebashI << "s1 line color 2" << std::endl;
        bfilebashI << "s1 linewidth 2" << std::endl;
        bfilebashI << "s2 line color 3" << std::endl;
        bfilebashI << "s2 linewidth 2" << std::endl;
        bfilebashI << "s3 line color 4" << std::endl;
        bfilebashI << "s3 linewidth 2" << std::endl;
        bfilebashI << "s4 line color 1" << std::endl;
        bfilebashI << "s4 linewidth 2" << std::endl;
        
        bfilebashI << "s5 symbol 1" << std::endl;
        bfilebashI << "s5 linestyle 0" << std::endl;
        bfilebashI << "s5 symbol color 1" << std::endl;
        bfilebashI << "s5 symbol size 1.2" << std::endl;
        bfilebashI << "s5 symbol linewidth 1.5" << std::endl;
        bfilebashI << "s5 errorbar linewidth 1.5" << std::endl;
        bfilebashI << "s5 errorbar color 1" << std::endl;
        
        bfilebashI << "s0 legend \"0-2-0 (\\xd\\0=0)\"" << std::endl;

        if (c_min120 <= chi_min) {
            bfilebashI << "s1 legend \"1-2-0 (\\xd\\0=" << delmin120 << "\\S+" << drp120-delmin120 << "\\N" << "\\s" << round((drm120-delmin120) * pow(10, 2))/pow(10, 2) << "\\N)\"" << std::endl;
        }else{
            bfilebashI << "s1 legend \"1-2-0 (\\xd\\0=" << delmin120 << ")\"" << std::endl;
        }
        if (c_min220 <= chi_min) {
            bfilebashI << "s2 legend \"2-2-0 (\\xd\\0=" << delmin220 << "\\S+" << drp220-delmin220 << "\\N" << "\\s" << round((drm220-delmin220) * 100) / 100 << "\\N)\"" << std::endl;
        }else{
            bfilebashI << "s2 legend \"2-2-0 (\\xd\\0=" << delmin220 << ")\"" << std::endl;
        }
        if (c_min320 <= chi_min) {
            bfilebashI << "s3 legend \"3-2-0 (\\xd\\0=" << round(delmin320 * pow(10, 2))/pow(10, 2) << "\\S+" << round((drp320-delmin320)*100)/100 << "\\N" << "\\s" << round((drm320-delmin320)*100)/100 << "\\N)\"" << std::endl;
        }else{
            bfilebashI << "s3 legend \"3-2-0 (\\xd\\0=" << round(delmin320 * pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        }
        bfilebashI << "s4 legend \"4-2-0 (\\xd\\0=0)\"" << std::endl;

        bfilebashI << "s5 legend \"exp " << peak << " keV\"" << std::endl; //EXPERIMENT PEAK keV

        

        //------------LEGEND CONFIGURATIONS---------------|
        bfilebashI << "legend 0.49, 0.77" << std::endl;
        bfilebashI << "legend char size 1.2" << std::endl;

        //------------------------------------------------|
        
        
        bfilebashI << "with string" << std::endl;
        bfilebashI << "string on" << std::endl;
        bfilebashI << "string loctype view" << std::endl;
        bfilebashI << "string 0.89, 0.65" << std::endl;
        bfilebashI << "string char size 1.6" << std::endl;
        bfilebashI << "string def \"(1-3)\\S+\\N(0\\S+\\N)\"" << std::endl;
        
        bfilebashI << "with string" << std::endl;
        bfilebashI << "string on" << std::endl;
        bfilebashI << "string loctype view" << std::endl;
        bfilebashI << "string 0.48, 0.44" << std::endl;
        bfilebashI << "string char size 1.3" << std::endl;
        bfilebashI << "string def \"gated on 658 keV \\xg\\0 ray (Ru-100)\"" << std::endl;
       
        bfilebashI << "xaxis label \"Cos(\\xq\\0) \"" << std::endl;
        bfilebashI << "xaxis label char size 1.6" << std::endl;
        bfilebashI << "xaxis ticklabel char size 1.6" << std::endl;
        bfilebashI << "xaxis tick major 0.5" << std::endl;
        bfilebashI << "xaxis tick minor 0.25" << std::endl;
        bfilebashI << "world xmin -1.2" << std::endl;
        bfilebashI << "world xmax 1.2" << std::endl;
        
        bfilebashI << "yaxis label \"Counts\"" << std::endl;
        bfilebashI << "yaxis label char size 1.6" << std::endl;
        bfilebashI << "yaxis ticklabel char size 1.6" << std::endl;
        
        bfilebashI << "world ymin " << ymin << std::endl;
        bfilebashI << "world ymax " << ymax << std::endl;
        
        if (exp_max >= 2 && exp_max < 10){
            bfilebashI << "yaxis tick major 1" << std::endl;
            bfilebashI << "yaxis tick minor 0.5" << std::endl;
            bfilebashI << "view 0.11, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max >= 0 && exp_max < 2){
            bfilebashI << "yaxis tick major 0.5" << std::endl;
            bfilebashI << "yaxis tick minor 0.1" << std::endl;
            bfilebashI << "view 0.11, 0.13, 1.26, 0.8" << std::endl;
        }else{
            std::cout << "Statistic for angular correlation analysis is too low or too high" << std::endl;
            return 0;
        }
        bfilebashI << "saveall \"AC_Ifin"<< ".arg\"" << std::endl;
        bfilebashI << "print to \"AC_Ifin"<< ".eps\"" << std::endl;
        bfilebashI << "device \"EPS\" OP \"level2\"" << std::endl;
        bfilebashI << "print" << std::endl;
        
        bfilebashI.close();
    }
    else
    {
        std::cout << "error: the Ifin file is not open";
    }

    std::stringstream bfileM;
    //bfile << Output2 << "ACscript" << peak << ".bfile"; //script file includes xmgrace Angular Correlation plots
    bfileM << "ACscript_Mix" << ".bfile"; //script file includes xmgrace Angular Correlation plots
    std::string bfilenameM = bfileM.str();
    std::ofstream bfilebashM(bfilenameM);
    
    if (bfilebashM.is_open())
    {

        bfilebashM << "read xy \"AC_Mix" << "_0-2-0_delta_0_Ru100.dat\"" << std::endl;
        bfilebashM << "read xy \"AC_Mix" << "_1-2-0_delta_min_Ru100.dat\"" << std::endl;
        bfilebashM << "read xy \"AC_Mix" << "_2-2-0_delta_min_Ru100.dat\"" << std::endl;
        bfilebashM << "read xy \"AC_Mix" << "_3-2-0_delta_min_Ru100.dat\"" << std::endl;
        bfilebashM << "read xy \"AC_Mix" << "_4-2-0_delta_0_Ru100.dat\"" << std::endl;
        bfilebashM << "read xydy \"" << "exp_Mix.dat\"" << std::endl;
        
        bfilebashM << "s0 line color 1" << std::endl;
        bfilebashM << "s0 linestyle 4" << std::endl;
        bfilebashM << "s0 linewidth 2" << std::endl;
        bfilebashM << "s1 line color 2" << std::endl;
        bfilebashM << "s1 linewidth 2" << std::endl;
        bfilebashM << "s2 line color 3" << std::endl;
        bfilebashM << "s2 linewidth 2" << std::endl;
        bfilebashM << "s3 line color 4" << std::endl;
        bfilebashM << "s3 linewidth 2" << std::endl;
        bfilebashM << "s4 line color 1" << std::endl;
        bfilebashM << "s4 linewidth 2" << std::endl;
        
        bfilebashM << "s5 symbol 1" << std::endl;
        bfilebashM << "s5 linestyle 0" << std::endl;
        bfilebashM << "s5 symbol color 1" << std::endl;
        bfilebashM << "s5 symbol size 1.2" << std::endl;
        bfilebashM << "s5 symbol linewidth 1.5" << std::endl;
        bfilebashM << "s5 errorbar linewidth 1.5" << std::endl;
        bfilebashM << "s5 errorbar color 1" << std::endl;
        
        bfilebashM << "s0 legend \"0-2-0 (\\xd\\0=0)\"" << std::endl;

        if (c_min120 <= chi_min) {
            bfilebashM << "s1 legend \"1-2-0 (\\xd\\0=" << delmin120 << "\\S+" << drp120-delmin120 << "\\N" << "\\s" << round((drm120-delmin120) * pow(10, 2))/pow(10, 2) << "\\N)\"" << std::endl;
        }else{
            bfilebashM << "s1 legend \"1-2-0 (\\xd\\0=" << delmin120 << ")\"" << std::endl;
        }
        if (c_min220 <= chi_min) {
            bfilebashM << "s2 legend \"2-2-0 (\\xd\\0=" << delmin220 << "\\S+" << drp220-delmin220 << "\\N" << "\\s" << round((drm220-delmin220) * 100) / 100 << "\\N)\"" << std::endl;
        }else{
            bfilebashM << "s2 legend \"2-2-0 (\\xd\\0=" << delmin220 << ")\"" << std::endl;
        }
        if (c_min320 <= chi_min) {
            bfilebashM << "s3 legend \"3-2-0 (\\xd\\0=" << round(delmin320 * pow(10, 2))/pow(10, 2) << "\\S+" << round((drp320-delmin320)*100)/100 << "\\N" << "\\s" << round((drm320-delmin320)*100)/100 << "\\N)\"" << std::endl;
        }else{
            bfilebashM << "s3 legend \"3-2-0 (\\xd\\0=" << round(delmin320 * pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        }
        bfilebashM << "s4 legend \"4-2-0 (\\xd\\0=0)\"" << std::endl;

        bfilebashM << "s5 legend \"exp " << peak << " keV\"" << std::endl; //EXPERIMENT PEAK keV

        

        //------------LEGEND CONFIGURATIONS---------------|
        bfilebashM << "legend 0.49, 0.77" << std::endl;
        bfilebashM << "legend char size 1.2" << std::endl;

        //------------------------------------------------|
        
        
        bfilebashM << "with string" << std::endl;
        bfilebashM << "string on" << std::endl;
        bfilebashM << "string loctype view" << std::endl;
        bfilebashM << "string 0.89, 0.65" << std::endl;
        bfilebashM << "string char size 1.6" << std::endl;
        bfilebashM << "string def \"(1-3)\\S+\\N(0\\S+\\N)\"" << std::endl;
        
        bfilebashM << "with string" << std::endl;
        bfilebashM << "string on" << std::endl;
        bfilebashM << "string loctype view" << std::endl;
        bfilebashM << "string 0.48, 0.44" << std::endl;
        bfilebashM << "string char size 1.3" << std::endl;
        bfilebashM << "string def \"gated on 658 keV \\xg\\0 ray (Ru-100)\"" << std::endl;
       
        bfilebashM << "xaxis label \"Cos(\\xq\\0) \"" << std::endl;
        bfilebashM << "xaxis label char size 1.6" << std::endl;
        bfilebashM << "xaxis ticklabel char size 1.6" << std::endl;
        bfilebashM << "xaxis tick major 0.5" << std::endl;
        bfilebashM << "xaxis tick minor 0.25" << std::endl;
        bfilebashM << "world xmin -1.2" << std::endl;
        bfilebashM << "world xmax 1.2" << std::endl;
        
        bfilebashM << "yaxis label \"Counts\"" << std::endl;
        bfilebashM << "yaxis label char size 1.6" << std::endl;
        bfilebashM << "yaxis ticklabel char size 1.6" << std::endl;
        
        bfilebashM << "world ymin " << ymin << std::endl;
        bfilebashM << "world ymax " << ymax << std::endl;
        
        if (exp_max >= 2 && exp_max < 10){
            bfilebashM << "yaxis tick major 1" << std::endl;
            bfilebashM << "yaxis tick minor 0.5" << std::endl;
            bfilebashM << "view 0.11, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max >= 0 && exp_max < 2){
            bfilebashM << "yaxis tick major 0.5" << std::endl;
            bfilebashM << "yaxis tick minor 0.1" << std::endl;
            bfilebashM << "view 0.11, 0.13, 1.26, 0.8" << std::endl;
        }else{
            std::cout << "Statistic for angular correlation analysis is too low or too high" << std::endl;
            return 0;
        }
        bfilebashM << "saveall \"AC_Mix"<< ".arg\"" << std::endl;
        bfilebashM << "print to \"AC_Mix"<< ".eps\"" << std::endl;
        bfilebashM << "device \"EPS\" OP \"level2\"" << std::endl;
        bfilebashM << "print" << std::endl;
        
        bfilebashM.close();
    }
    else
    {
        std::cout << "error: the Mixed file is not open";
    }

    //-------------------------------------------------------------------------------------------------------------------|
    float yminchi=0;
    float ymaxchi=0;

    yminchi= 0.3 * chi_min;
    ymaxchi= 4.0 * chi_max;
   
    
    std::stringstream Chi;
    Chi << "Chiscript" << ".bfile";//script file includes xmgrace Chi Square plots
    std::string bfileChi = Chi.str();
    std::ofstream bfilebashChi(bfileChi);
    
    if (bfilebashChi.is_open())
    {
        bfilebashChi << "read xy \"chi_square" << "_0-2-0_delta_0_Ru100.dat\"" << std::endl;
        bfilebashChi << "read xy \"chi_square" << "_1-2-0_Ru100.dat\"" << std::endl;
        bfilebashChi << "read xy \"chi_square" << "_2-2-0_Ru100.dat\"" << std::endl;
        bfilebashChi << "read xy \"chi_square" << "_3-2-0_Ru100.dat\"" << std::endl;
        bfilebashChi << "read xy \"chi_square" << "_4-2-0_delta_0_Ru100.dat\"" << std::endl;
        bfilebashChi << "read xy \"confidencelev95_Ru100.dat\"" << std::endl;
        
        bfilebashChi << "s0 symbol 1" << std::endl;
        bfilebashChi << "s0 linestyle 0" << std::endl;
        bfilebashChi << "s0 symbol color 1" << std::endl;
        bfilebashChi << "s0 symbol fill 1" << std::endl;
        bfilebashChi << "s0 symbol size 1.2" << std::endl;
        bfilebashChi << "s0 symbol linewidth 1.5" << std::endl;
        
        bfilebashChi << "s1 line color 2" << std::endl;
        bfilebashChi << "s1 linewidth 2" << std::endl;
        bfilebashChi << "s2 line color 3" << std::endl;
        bfilebashChi << "s2 linewidth 2" << std::endl;
        bfilebashChi << "s3 line color 4" << std::endl;
        bfilebashChi << "s3 linewidth 2" << std::endl;
        
        bfilebashChi << "s4 symbol 2" << std::endl;
        bfilebashChi << "s4 linestyle 0" << std::endl;
        bfilebashChi << "s4 symbol color 1" << std::endl;
        bfilebashChi << "s4 symbol fill 1" << std::endl;
        bfilebashChi << "s4 symbol size 1.2" << std::endl;
        bfilebashChi << "s4 symbol linewidth 1.5" << std::endl;
        
        bfilebashChi << "s5 linestyle 4" << std::endl;
        bfilebashChi << "s5 line color 1" << std::endl;
        bfilebashChi << "s5 line linewidth 1.5" << std::endl;
        
        
        bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
        bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        
        
        bfilebashChi << "legend char size 1.2" << std::endl;
        
        bfilebashChi << "with string" << std::endl;
        bfilebashChi << "string on" << std::endl;
        bfilebashChi << "string loctype view" << std::endl;
        

        bfilebashChi << "string 0.75, 0.2" << std::endl;
        bfilebashChi << "string char size 1.3" << std::endl;
        bfilebashChi << "string def \"" << peak << " keV \\xg\\0 ray\"" << std::endl; //PEAK GAMMA
        
        bfilebashChi << "with string" << std::endl;
        bfilebashChi << "string on" << std::endl;
        bfilebashChi << "string loctype view" << std::endl;
        bfilebashChi << "string 0.1, 0.23" << std::endl;
        bfilebashChi << "string char size 1.3" << std::endl;
        bfilebashChi << "string def \"95\%\"" << std::endl;
        
        bfilebashChi << "with string" << std::endl;
        bfilebashChi << "string on" << std::endl;
        bfilebashChi << "string loctype view" << std::endl;
        bfilebashChi << "string 0.5, 0.74" << std::endl;
        bfilebashChi << "string char size 1.3" << std::endl;
        bfilebashChi << "string def \"gated on 539 keV \\xg\\0 ray (Ru-100)\"" << std::endl; //CHANGE EXPERIMENTAL
       
        
        bfilebashChi << "xaxis label \"arctan(\\xd\\0), deg.\"" << std::endl;
        bfilebashChi << "xaxis label char size 1.6" << std::endl;
        bfilebashChi << "xaxis ticklabel char size 1.6" << std::endl;
        bfilebashChi << "xaxis tick major 30" << std::endl;
        bfilebashChi << "xaxis tick minor 10" << std::endl;
        bfilebashChi << "world xmin -90" << std::endl;
        bfilebashChi << "world xmax 90" << std::endl;
        
        bfilebashChi << "yaxis label \"\\xc\\N\\S2\\N/\\xn\"" << std::endl;
        bfilebashChi << "yaxis label char size 1.6" << std::endl;
        bfilebashChi << "yaxis ticklabel char size 1.6" << std::endl;
        bfilebashChi << "yaxes scale logarithmic" << std::endl;
        bfilebashChi << "yaxis ticklabel format power" << std::endl;
        bfilebashChi << "yaxis ticklabel prec 0" << std::endl;
        bfilebashChi << "yaxis tick major 10" << std::endl;
        bfilebashChi << "yaxis tick minor 1" << std::endl;
        
        
        bfilebashChi << "world ymin " << yminchi << std::endl;
        bfilebashChi << "world ymax " << ymaxchi << std::endl;
        
        bfilebashChi << "view 0.17, 0.13, 1.26, 0.8" << std::endl;
        bfilebashChi << "saveall \"Chi_"<< ".arg\"" << std::endl;
        bfilebashChi << "print to \"Chi_"<< ".eps\"" << std::endl;
        bfilebashChi << "device \"EPS\" OP \"level2\"" << std::endl;
        bfilebashChi << "print" << std::endl;
        
        bfilebashChi.close();
    }
    else
    {
        std::cout << "error: the file is not open";
    }
    
    std::stringstream chib;
    chib <<  "AChi" << "bash.txt";
    std::string bfileb = chib.str();
    std::ofstream bbashAChi(bfileb);
    
    if (bbashAChi.is_open())
    {
        bbashAChi << "xmgrace -batch Chiscript" << ".bfile -nosafe -hardcopy" << std::endl;
        bbashAChi << "xmgrace -batch ACscript_Fipps" << ".bfile -nosafe -hardcopy" << std::endl;
        bbashAChi << "xmgrace -batch ACscript_Ifin" << ".bfile -nosafe -hardcopy" << std::endl;
        bbashAChi << "xmgrace -batch ACscript_Mix" << ".bfile -nosafe -hardcopy" << std::endl;

        bbashAChi << "ps2pdf" << " AC_Fipps" << ".eps" << std::endl;
        bbashAChi << "pdfcrop" << " AC_Fipps" << ".pdf" << std::endl;

        bbashAChi << "ps2pdf" << " AC_Ifin" << ".eps" << std::endl;
        bbashAChi << "pdfcrop" << " AC_Ifin" << ".pdf" << std::endl;

        bbashAChi << "ps2pdf" << " AC_Mix" << ".eps" << std::endl;
        bbashAChi << "pdfcrop" << " AC_Mix" << ".pdf" << std::endl;

        bbashAChi << "ps2pdf" << " Chi_" << ".eps" << std::endl;
        bbashAChi << "pdfcrop" << " Chi_" << ".pdf" << std::endl;

        bbashAChi << "xpdf" << " Chi_" << "-crop.pdf" << std::endl;
        bbashAChi << "xpdf" << " AC_" << "-crop.pdf" << std::endl;
        bbashAChi.close();
    }
    else
    {
        std::cout << "error: the AChi" << "bash.txt file is not open";
    }
        
   return 0;    
 }    
    





