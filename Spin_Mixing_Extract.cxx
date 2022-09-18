//
//  main.cpp
//  ChiXe126_x-2-0
//
//  Created by Sambuu on 2019-07-10.
//  Copyright © 2019 Sambuu. All rights reserved.
//  Modifeid for 126Xe by Farnaz Ghazi Moradi January 2020
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

int main()
{
    float peak = 0;
    std::cout << "Enter the energy of the transition for peak fit: ";
    std::cin >> peak;

    std::string line;
    //std::string Output0 = "/home/fghazimo/My_Analysis/8Pi/126Xe/2020/Angular_correlation_for_8pi/matrices/x-2-0/";
    std::string Output0 = "/home/fghazimo/My_Analysis/8Pi/126Xe/2020/Angular_correlation_for_8pi/x-2-0/";
    std::stringstream transdir;
    transdir << Output0 << peak << std::endl;
    std::string transdirfilename = transdir.str();
    
    std::stringstream engdir;
    engdir << Output0 << peak << "/eng/";
    std::string engdirfilename = engdir.str();
    std::stringstream datadir;
    datadir << Output0 << peak << "/data/";
    std::string datadirfilename = datadir.str();
    std::stringstream peakareadir;
    peakareadir << Output0 << peak << "/peakarea/";
    std::string peakareadirfilename = peakareadir.str();
//    std::cout << "Creating directories" << " " << "\033[1;34m\t " << peak << " \033[0m" << std::endl;
  //  std::cout << "The output of " << peak<< ".sh "<< " (" << peak << ".eng and " << peak << "gf3.sto) are saved in:" << " " << "\033[1;32m\t " << engdirfilename << " \033[0m" << std::endl;
    std::cout << ".eng and .sto files will be saved in:" << "\033[1;32m " << engdirfilename << " \033[0m" << std::endl;
    std::cout << ".dat files will be saved in:" << "\033[1;32m " << datadirfilename << " \033[0m" << std::endl;
    std::cout << "Peakareas will be saved in:" << "\033[1;32m " << peakareadirfilename << " \033[0m" << std::endl;

    //std::string Output1 = "/home/fghazimo/My_Analysis/8Pi/126Xe/2020/Angular_correlation_for_8pi/x-2-0/1970/Xe126eng/";
    //std::string Output2 = "/home/fghazimo/My_Analysis/8Pi/126Xe/2020/Angular_correlation_for_8pi/x-2-0/1970/Xe126data/";
    //std::string Output3 = "/home/fghazimo/My_Analysis/8Pi/126Xe/2020/Angular_correlation_for_8pi/x-2-0/1970/Xe126peakarea/";

//    std::string Output1 = "/Users/bjigmedd/Sambuu/PhD_project/Xe126/Xe126sto/";
//    std::string Output2 = "/Users/bjigmedd/Sambuu/PhD_project/Xe126/Xe126data/";
//    std::string Output3 = "/Users/bjigmedd/Sambuu/PhD_project/Xe126/Xe126peakarea/";
    
    
    int J2 = 2, J3 = 0;
    float ac_min = 0;
    float ac_max = 0;
    float Wtheomixed = 0;
    float Wtheomixed_0 = 0;
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
    double ChiSquare = 0;
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
    
    float Q2 = 0.96;//quenching factors
    float Q4 = 0.89;

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
    
    std::stringstream engout;
    engout << peakareadirfilename << peak << "engout.txt";
    std::string engoutfilename = engout.str();
    
    std::stringstream filename;
    filename << engdirfilename << peak << ".eng";
    std::string engname = filename.str();
    
    std::ifstream file(engname); //opening eng file which was converted to file.txt
    if (file.is_open())
    {
        
        std::ofstream out(engoutfilename); //creating newfile in which only numerical data lines will be saved
        
        while (getline(file, line))
        {
            //taking out texts from eng (txt) file
            if (line != " No.  Centroid +- error      Area +- error      Energy +- error    Sp.name      Date     Time")
                
                out << line << "\n";
            
        }
        out.close();
        file.close();
    }
    else
    {
        std::cout << "error: file \""<<engname<<"\" is not open" << std::endl;
    }
    std::string lline;
    
    std::ifstream infile(engoutfilename); //opening file which includes only peak area and their uncertainties and date and time and energy
    
    std::stringstream peakarea;
    //peakarea << Output3 << peak << "peakarea_Xe126.txt";
    peakarea << peakareadirfilename << peak << "peakarea.txt";
    std::string peakareafilename = peakarea.str();
    
    float peakmax = peak + 1.2;
    float peakmin = peak - 1.2;
    
    if (infile.is_open())
    {
        std::ofstream outpeak(peakareafilename); //opening file in which only peak area and its uncertainties will be saved
        
        while (std::getline(infile, lline))
        {
            std::stringstream ss(lline);
            int a, d, e;
            float b, c, f, g;
            std::string h, i, j;
            
            if (ss >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j)
            {
                if (f >= peakmin && f <= peakmax) {
                    
                    outpeak << d << " " << e << std::endl;
                }
                // Add a, b, and c to their respective arrays
            }
        }
        outpeak.close();
        infile.close();
        
    }
    else
    {
        std::cout << "error: infile is not open" << std::endl;
    }
    
    //float p42 = 24, p71 = 49, p109 = 48, p138 = 24, p180 = 8; //
    float p42 = 27, p71 = 54, p109 = 54, p138 = 27, p180 = 9; //Farnaz-pair number for S1292 126Xe Ge 14 was not considered
    //float p42 = 30, p71 = 60, p109 = 60, p138 = 30, p180 = 10; //
    //float p42 = 21, p71 = 45, p109 = 42, p138 = 21, p180 = 7; //pair number
    float data[10], a42, a42r, a71, a71r, a109, a109r, a138, a138r, a180, a180r; //r for uncertainty
    float NorToPair42, NorToPair71, NorToPair109, NorToPair138, NorToPair180;
    float NorToPair42r, NorToPair71r, NorToPair109r, NorToPair138r, NorToPair180r; // r for uncertainty
    float Sys42r, Sys71r, Sys109r, Sys138r, Sys180r; // r for uncertainty
    //float R42 = 0.803416, R71 = 0.884586, R109 = 0.901429, R138 = 0.824117, R180 = 1.000; //relative renormalization of 818 transition for each angle
    //float R42r = 0.00583759, R71r = 0.00540157, R109r = 0.00553977, R138r = 0.00603736, R180r = 0.00726115; // r for uncertainty
    
    float CorPeakArea42, CorPeakArea71, CorPeakArea109, CorPeakArea138, CorPeakArea180;
    float CorPeakArea42r, CorPeakArea71r, CorPeakArea109r, CorPeakArea138r, CorPeakArea180r;
    
    std::ifstream peakfile;
    
    peakfile.open(peakareafilename);
    
    for (int i = 0; i <= 9; i++)
    {
        peakfile >> data[i];
	std::cout << "\033[1;31m\t " << data[i] <<" \033[0m" << std::endl;
        
    }
    peakfile.close();
    
    a42 = data[0];
    a42r = data[1];
    a71 = data[2];
    a71r = data[3];
    a109 = data[4];
    a109r = data[5];
    a138 = data[6];
    a138r = data[7];
    a180 = data[8];
    a180r = data[9];
    
    int AngleList [5] = {42, 71, 109, 138, 180};
    NorToPair42 = a42 / p42;
    NorToPair71 = a71 / p71;
    NorToPair109 = a109 / p109;
    NorToPair138 = a138 / p138;
    NorToPair180 = a180 / p180;
    
    NorToPair42r = a42r / p42;
    NorToPair71r = a71r / p71;
    NorToPair109r = a109r / p109;
    NorToPair138r = a138r / p138;
    NorToPair180r = a180r / p180;
    
    Sys42r = 0.01 * NorToPair42; // 2.5% systematic uncertainty adopted to make chi^2 close to 1 for 818 keV transition 0_2_0
    Sys71r = 0.01 * NorToPair71;
    Sys109r = 0.01 * NorToPair109;
    Sys138r = 0.01 * NorToPair138;
    Sys180r = 0.01 * NorToPair180;
    
    CorPeakArea42 = NorToPair42;
    CorPeakArea71 = NorToPair71;
    CorPeakArea109 = NorToPair109;
    CorPeakArea138 = NorToPair138;
    CorPeakArea180 = NorToPair180;
    
    CorPeakArea42r = sqrt(pow(NorToPair42r, 2)+pow(Sys42r, 2));
    CorPeakArea71r = sqrt(pow(NorToPair71r, 2)+pow(Sys71r, 2));
    CorPeakArea109r = sqrt(pow(NorToPair109r, 2)+pow(Sys109r, 2));
    CorPeakArea138r = sqrt(pow(NorToPair138r, 2)+pow(Sys138r, 2));
    CorPeakArea180r = sqrt(pow(NorToPair180r, 2)+pow(Sys180r, 2));
    
    if (CorPeakArea42 > exp_max) {
        exp_max = CorPeakArea42;
        exp_min = CorPeakArea42;
    }
    if (CorPeakArea71 > exp_max) {
        exp_max = CorPeakArea71;
    }else if (CorPeakArea71 < exp_min) {
        exp_min = CorPeakArea71;
    }
    if (CorPeakArea109 > exp_max){
        exp_max = CorPeakArea109;
    }else if (CorPeakArea109 < exp_min){
        exp_min = CorPeakArea109;
    }
    if(CorPeakArea138 > exp_max){
        exp_max = CorPeakArea138;
    }else if (CorPeakArea138 < exp_min){
        exp_min = CorPeakArea138;
    }
    if(CorPeakArea180 > exp_max){
        exp_max = CorPeakArea180;
    }else if (CorPeakArea180 < exp_min){
        exp_min = CorPeakArea180;
    }
    
    std::cout << exp_max << " " << exp_min << std::endl;
    
    
    for (int J1 = 0; J1 < 5; J1++){
        
        int cchi = 0;
        int cochi = 0;
        int countchi0 = 0;
        
        std::stringstream cfilename;
        
        //cfilename << Output2 << peak << "exp.dat";
        cfilename << datadirfilename << peak << "exp.dat";
        std::string cname = cfilename.str();
        
        std::ofstream cfile (cname);
        if (cfile.is_open())
        {
            
            cfile << cos(AngleList [0]*PI/180) << " " << CorPeakArea42 << " " << CorPeakArea42r << std::endl;
            cfile << cos(AngleList [1]*PI/180) << " " << CorPeakArea71 << " " << CorPeakArea71r << std::endl;
            cfile << cos(AngleList [2]*PI/180) << " " << CorPeakArea109 << " " << CorPeakArea109r << std::endl;
            cfile << cos(AngleList [3]*PI/180) << " " << CorPeakArea138 << " " << CorPeakArea138r << std::endl;
            cfile << cos(AngleList [4]*PI/180) << " " << CorPeakArea180 << " " << CorPeakArea180r << std::endl;
        }
        cfile.close();
        
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
        
        /*std::stringstream a2a4filename;
         a2a4filename << Output1 << "a2a4_" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_Xe126.dat";
         std::string a2a4name = a2a4filename.str();
         
         std::ofstream a2a4file (a2a4name);
         if (a2a4file.is_open())
         {
         std::stringstream a4filename;
         a4filename << Output1 << "a4_" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_Xe126.dat";
         std::string a4name = a4filename.str();
         
         std::ofstream a4file (a4name);
         if (a4file.is_open())
         {
         
         
         std::stringstream a2filename;
         a2filename << Output1 << "a2_" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_Xe126.dat";
         std::string a2name = a2filename.str();
         
         std::ofstream a2file (a2name);
         if (a2file.is_open())
         {
         */
        std::stringstream chifilename;
        //chifilename << Output2 << "chi_square" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_Xe126.dat";
        chifilename << datadirfilename << "chi_square" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_Xe126.dat";
        std::string chiname = chifilename.str();
        
        std::ofstream chifile (chiname);
        if (chifile.is_open())
        {
            double delta1;
            for (delta1 = -50.000; delta1 < 50.000; delta1 += 0.010){
                
                int delta2 = 0;
                
                float a2 = Q2 * ((R2LLJ2J1 + 2 * delta1 * R2LMJ2J1 + pow(delta1, 2) * R2MMJ2J1) * (R2LLJ2J3 + 2 * delta2 * R2LMJ2J3 + pow(delta2, 2) * R2MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * delta1 * R4LMJ2J1 + pow(delta1, 2) * R4MMJ2J1) * (R4LLJ2J3 + 2 * delta2 * R4LMJ2J3 + pow(delta2, 2) * R4MMJ2J3) / (1 + pow(delta1, 2)) / (1 + pow(delta2, 2)));
                
                float Wtheo = 0;
                float t42 = 0, t71 = 0, t109 = 0, t138 = 0, t180 = 0;
                
                for (int theta = 0; theta < 181; theta++)
                {
                    
                    Wtheo = 1 + a2 * (3 * pow(cos (theta*PI/180), 2) - 1) / 2 + a4 * (35 * pow(cos (theta*PI/180), 4) - 30 * pow(cos (theta*PI/180), 2) + 3) / 8;
                    
                    
                    if (theta == 42)
                    {
                        t42 = Wtheo;
                    }
                    if (theta == 71)
                    {
                        t71 = Wtheo;
                    }
                    if (theta == 109)
                    {
                        t109 = Wtheo;
                    }
                    if (theta == 138)
                    {
                        t138 = Wtheo;
                    }
                    if (theta == 180)
                    {
                        t180 = Wtheo;
                    }
                }
                
                float wcoeff = 0;
                
                wcoeff = (CorPeakArea42 * t42 / pow(CorPeakArea42r, 2) + CorPeakArea71 * t71 / pow(CorPeakArea71r, 2) + CorPeakArea109 * t109 / pow(CorPeakArea109r, 2) + CorPeakArea138 * t138 / pow(CorPeakArea138r, 2) + CorPeakArea180 * t180 / pow(CorPeakArea180r, 2)) / (pow(t42, 2) / pow(CorPeakArea42r, 2) + pow(t71, 2) / pow(CorPeakArea71r, 2) + pow(t109, 2) / pow(CorPeakArea109r, 2) + pow(t138, 2) / pow(CorPeakArea138r, 2) + pow(t180, 2) / pow(CorPeakArea180r, 2)); // wcoeff is weighted coefficient
                
                ChiSquare = (pow((CorPeakArea42 - wcoeff * t42), 2) / pow(CorPeakArea42r, 2) + pow((CorPeakArea71 - wcoeff * t71), 2) / pow(CorPeakArea71r, 2) + pow((CorPeakArea109 - wcoeff * t109), 2) / pow(CorPeakArea109r, 2) + pow((CorPeakArea138 - wcoeff * t138), 2) / pow(CorPeakArea138r, 2) + pow((CorPeakArea180 - wcoeff * t180), 2) / pow(CorPeakArea180r, 2))/4;
                
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
                
                chifile << atan (delta1) * 180 / PI << "  " << ChiSquare << std::endl;
                /*
                 a2file << delta1 << "  " << a2 << std::endl;
                 a4file << delta1 << "  " << a4 << std::endl;
                 a2a4file << a2 << "  " << a4 << std::endl;
                 
                 }
                 a2file.close();
                 }
                 a4file.close();
                 }
                 a2a4file.close();
                 */
            }
            chifile.close();
        }
        
        // starting to find uncetainty on mixing ratio based on 1 sigma uncertainty on chi^2
        double d1;
        if (J1 == 1 && J2 == 2 && J3 == 0) {
            for (d1 = delmin120; d1 <= 51.000; d1 += 0.010){
                
                float a2 = Q2 * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                
                float Wtheo = 0;
                float t42 = 0, t71 = 0, t109 = 0, t138 = 0, t180 = 0;
                for (int theta = 0; theta < 181; theta++)
                {
                    
                    Wtheo = 1 + a2 * (3 * pow(cos (theta*PI/180), 2) - 1) / 2 + a4 * (35 * pow(cos (theta*PI/180), 4) - 30 * pow(cos (theta*PI/180), 2) + 3) / 8;
                    
                    if (theta == 42)
                    {
                        t42 = Wtheo;
                    }
                    if (theta == 71)
                    {
                        t71 = Wtheo;
                    }
                    if (theta == 109)
                    {
                        t109 = Wtheo;
                    }
                    if (theta == 138)
                    {
                        t138 = Wtheo;
                    }
                    if (theta == 180)
                    {
                        t180 = Wtheo;
                    }
                }
                
                float wcoeff = 0;
                
                wcoeff = (CorPeakArea42 * t42 / pow(CorPeakArea42r, 2) + CorPeakArea71 * t71 / pow(CorPeakArea71r, 2) + CorPeakArea109 * t109 / pow(CorPeakArea109r, 2) + CorPeakArea138 * t138 / pow(CorPeakArea138r, 2) + CorPeakArea180 * t180 / pow(CorPeakArea180r, 2)) / (pow(t42, 2) / pow(CorPeakArea42r, 2) + pow(t71, 2) / pow(CorPeakArea71r, 2) + pow(t109, 2) / pow(CorPeakArea109r, 2) + pow(t138, 2) / pow(CorPeakArea138r, 2) + pow(t180, 2) / pow(CorPeakArea180r, 2)); // wcoeff is weighted coefficient
                
                ChiSquare2 = (pow((CorPeakArea42 - wcoeff * t42), 2) / pow(CorPeakArea42r, 2) + pow((CorPeakArea71 - wcoeff * t71), 2) / pow(CorPeakArea71r, 2) + pow((CorPeakArea109 - wcoeff * t109), 2) / pow(CorPeakArea109r, 2) + pow((CorPeakArea138 - wcoeff * t138), 2) / pow(CorPeakArea138r, 2) + pow((CorPeakArea180 - wcoeff * t180), 2) / pow(CorPeakArea180r, 2))/4;
                
                float chim = (c_min120 * 4 + 1) / 4;
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
                
                float a2 = Q2 * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                
                float Wtheo = 0;
                float t42 = 0, t71 = 0, t109 = 0, t138 = 0, t180 = 0;
                for (int theta = 0; theta < 181; theta++)
                {
                    
                    Wtheo = 1 + a2 * (3 * pow(cos (theta*PI/180), 2) - 1) / 2 + a4 * (35 * pow(cos (theta*PI/180), 4) - 30 * pow(cos (theta*PI/180), 2) + 3) / 8;
                    
                    if (theta == 42)
                    {
                        t42 = Wtheo;
                    }
                    if (theta == 71)
                    {
                        t71 = Wtheo;
                    }
                    if (theta == 109)
                    {
                        t109 = Wtheo;
                    }
                    if (theta == 138)
                    {
                        t138 = Wtheo;
                    }
                    if (theta == 180)
                    {
                        t180 = Wtheo;
                    }
                }
                
                float wcoeff = 0;
                
                wcoeff = (CorPeakArea42 * t42 / pow(CorPeakArea42r, 2) + CorPeakArea71 * t71 / pow(CorPeakArea71r, 2) + CorPeakArea109 * t109 / pow(CorPeakArea109r, 2) + CorPeakArea138 * t138 / pow(CorPeakArea138r, 2) + CorPeakArea180 * t180 / pow(CorPeakArea180r, 2)) / (pow(t42, 2) / pow(CorPeakArea42r, 2) + pow(t71, 2) / pow(CorPeakArea71r, 2) + pow(t109, 2) / pow(CorPeakArea109r, 2) + pow(t138, 2) / pow(CorPeakArea138r, 2) + pow(t180, 2) / pow(CorPeakArea180r, 2)); // wcoeff is weighted coefficient
                
                ChiSquare2 = (pow((CorPeakArea42 - wcoeff * t42), 2) / pow(CorPeakArea42r, 2) + pow((CorPeakArea71 - wcoeff * t71), 2) / pow(CorPeakArea71r, 2) + pow((CorPeakArea109 - wcoeff * t109), 2) / pow(CorPeakArea109r, 2) + pow((CorPeakArea138 - wcoeff * t138), 2) / pow(CorPeakArea138r, 2) + pow((CorPeakArea180 - wcoeff * t180), 2) / pow(CorPeakArea180r, 2))/4;
                
                float chim = (c_min120 * 4 + 1) / 4;
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
                
                float a2 = Q2 * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                
                float Wtheo = 0;
                float t42 = 0, t71 = 0, t109 = 0, t138 = 0, t180 = 0;
                for (int theta = 0; theta < 181; theta++)
                {
                    
                    Wtheo = 1 + a2 * (3 * pow(cos (theta*PI/180), 2) - 1) / 2 + a4 * (35 * pow(cos (theta*PI/180), 4) - 30 * pow(cos (theta*PI/180), 2) + 3) / 8;
                    
                    if (theta == 42)
                    {
                        t42 = Wtheo;
                    }
                    if (theta == 71)
                    {
                        t71 = Wtheo;
                    }
                    if (theta == 109)
                    {
                        t109 = Wtheo;
                    }
                    if (theta == 138)
                    {
                        t138 = Wtheo;
                    }
                    if (theta == 180)
                    {
                        t180 = Wtheo;
                    }
                }
                
                float wcoeff = 0;
                
                wcoeff = (CorPeakArea42 * t42 / pow(CorPeakArea42r, 2) + CorPeakArea71 * t71 / pow(CorPeakArea71r, 2) + CorPeakArea109 * t109 / pow(CorPeakArea109r, 2) + CorPeakArea138 * t138 / pow(CorPeakArea138r, 2) + CorPeakArea180 * t180 / pow(CorPeakArea180r, 2)) / (pow(t42, 2) / pow(CorPeakArea42r, 2) + pow(t71, 2) / pow(CorPeakArea71r, 2) + pow(t109, 2) / pow(CorPeakArea109r, 2) + pow(t138, 2) / pow(CorPeakArea138r, 2) + pow(t180, 2) / pow(CorPeakArea180r, 2)); // wcoeff is weighted coefficient
                
                ChiSquare2 = (pow((CorPeakArea42 - wcoeff * t42), 2) / pow(CorPeakArea42r, 2) + pow((CorPeakArea71 - wcoeff * t71), 2) / pow(CorPeakArea71r, 2) + pow((CorPeakArea109 - wcoeff * t109), 2) / pow(CorPeakArea109r, 2) + pow((CorPeakArea138 - wcoeff * t138), 2) / pow(CorPeakArea138r, 2) + pow((CorPeakArea180 - wcoeff * t180), 2) / pow(CorPeakArea180r, 2))/4;
                
                float chim = (c_min220 * 4 + 1) / 4;
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
                
                float a2 = Q2 * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                
                float Wtheo = 0;
                float t42 = 0, t71 = 0, t109 = 0, t138 = 0, t180 = 0;
                for (int theta = 0; theta < 181; theta++)
                {
                    
                    Wtheo = 1 + a2 * (3 * pow(cos (theta*PI/180), 2) - 1) / 2 + a4 * (35 * pow(cos (theta*PI/180), 4) - 30 * pow(cos (theta*PI/180), 2) + 3) / 8;
                    
                    if (theta == 42)
                    {
                        t42 = Wtheo;
                    }
                    if (theta == 71)
                    {
                        t71 = Wtheo;
                    }
                    if (theta == 109)
                    {
                        t109 = Wtheo;
                    }
                    if (theta == 138)
                    {
                        t138 = Wtheo;
                    }
                    if (theta == 180)
                    {
                        t180 = Wtheo;
                    }
                }
                
                float wcoeff = 0;
                
                wcoeff = (CorPeakArea42 * t42 / pow(CorPeakArea42r, 2) + CorPeakArea71 * t71 / pow(CorPeakArea71r, 2) + CorPeakArea109 * t109 / pow(CorPeakArea109r, 2) + CorPeakArea138 * t138 / pow(CorPeakArea138r, 2) + CorPeakArea180 * t180 / pow(CorPeakArea180r, 2)) / (pow(t42, 2) / pow(CorPeakArea42r, 2) + pow(t71, 2) / pow(CorPeakArea71r, 2) + pow(t109, 2) / pow(CorPeakArea109r, 2) + pow(t138, 2) / pow(CorPeakArea138r, 2) + pow(t180, 2) / pow(CorPeakArea180r, 2)); // wcoeff is weighted coefficient
                
                ChiSquare2 = (pow((CorPeakArea42 - wcoeff * t42), 2) / pow(CorPeakArea42r, 2) + pow((CorPeakArea71 - wcoeff * t71), 2) / pow(CorPeakArea71r, 2) + pow((CorPeakArea109 - wcoeff * t109), 2) / pow(CorPeakArea109r, 2) + pow((CorPeakArea138 - wcoeff * t138), 2) / pow(CorPeakArea138r, 2) + pow((CorPeakArea180 - wcoeff * t180), 2) / pow(CorPeakArea180r, 2))/4;
                
                float chim = (c_min220 * 4 + 1) / 4;
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
                
                float a2 = Q2 * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                
                float Wtheo = 0;
                float t42 = 0, t71 = 0, t109 = 0, t138 = 0, t180 = 0;
                for (int theta = 0; theta < 181; theta++)
                {
                    
                    Wtheo = 1 + a2 * (3 * pow(cos (theta*PI/180), 2) - 1) / 2 + a4 * (35 * pow(cos (theta*PI/180), 4) - 30 * pow(cos (theta*PI/180), 2) + 3) / 8;
                    
                    if (theta == 42)
                    {
                        t42 = Wtheo;
                    }
                    if (theta == 71)
                    {
                        t71 = Wtheo;
                    }
                    if (theta == 109)
                    {
                        t109 = Wtheo;
                    }
                    if (theta == 138)
                    {
                        t138 = Wtheo;
                    }
                    if (theta == 180)
                    {
                        t180 = Wtheo;
                    }
                }
                
                float wcoeff = 0;
                
                wcoeff = (CorPeakArea42 * t42 / pow(CorPeakArea42r, 2) + CorPeakArea71 * t71 / pow(CorPeakArea71r, 2) + CorPeakArea109 * t109 / pow(CorPeakArea109r, 2) + CorPeakArea138 * t138 / pow(CorPeakArea138r, 2) + CorPeakArea180 * t180 / pow(CorPeakArea180r, 2)) / (pow(t42, 2) / pow(CorPeakArea42r, 2) + pow(t71, 2) / pow(CorPeakArea71r, 2) + pow(t109, 2) / pow(CorPeakArea109r, 2) + pow(t138, 2) / pow(CorPeakArea138r, 2) + pow(t180, 2) / pow(CorPeakArea180r, 2)); // wcoeff is weighted coefficient
                
                ChiSquare2 = (pow((CorPeakArea42 - wcoeff * t42), 2) / pow(CorPeakArea42r, 2) + pow((CorPeakArea71 - wcoeff * t71), 2) / pow(CorPeakArea71r, 2) + pow((CorPeakArea109 - wcoeff * t109), 2) / pow(CorPeakArea109r, 2) + pow((CorPeakArea138 - wcoeff * t138), 2) / pow(CorPeakArea138r, 2) + pow((CorPeakArea180 - wcoeff * t180), 2) / pow(CorPeakArea180r, 2))/4;
                
                float chim = (c_min320 * 4 + 1) / 4;
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
                
                float a2 = Q2 * ((R2LLJ2J1 + 2 * d1 * R2LMJ2J1 + pow(d1, 2) * R2MMJ2J1) * (R2LLJ2J3) / (1 + pow(d1, 2)));
                float a4 = Q4 * ((R4LLJ2J1 + 2 * d1 * R4LMJ2J1 + pow(d1, 2) * R4MMJ2J1) * (R4LLJ2J3) / (1 + pow(d1, 2)));
                
                float Wtheo = 0;
                float t42 = 0, t71 = 0, t109 = 0, t138 = 0, t180 = 0;
                for (int theta = 0; theta < 181; theta++)
                {
                    
                    Wtheo = 1 + a2 * (3 * pow(cos (theta*PI/180), 2) - 1) / 2 + a4 * (35 * pow(cos (theta*PI/180), 4) - 30 * pow(cos (theta*PI/180), 2) + 3) / 8;
                    
                    if (theta == 42)
                    {
                        t42 = Wtheo;
                    }
                    if (theta == 71)
                    {
                        t71 = Wtheo;
                    }
                    if (theta == 109)
                    {
                        t109 = Wtheo;
                    }
                    if (theta == 138)
                    {
                        t138 = Wtheo;
                    }
                    if (theta == 180)
                    {
                        t180 = Wtheo;
                    }
                }
                
                float wcoeff = 0;
                
                wcoeff = (CorPeakArea42 * t42 / pow(CorPeakArea42r, 2) + CorPeakArea71 * t71 / pow(CorPeakArea71r, 2) + CorPeakArea109 * t109 / pow(CorPeakArea109r, 2) + CorPeakArea138 * t138 / pow(CorPeakArea138r, 2) + CorPeakArea180 * t180 / pow(CorPeakArea180r, 2)) / (pow(t42, 2) / pow(CorPeakArea42r, 2) + pow(t71, 2) / pow(CorPeakArea71r, 2) + pow(t109, 2) / pow(CorPeakArea109r, 2) + pow(t138, 2) / pow(CorPeakArea138r, 2) + pow(t180, 2) / pow(CorPeakArea180r, 2)); // wcoeff is weighted coefficient
                
                ChiSquare2 = (pow((CorPeakArea42 - wcoeff * t42), 2) / pow(CorPeakArea42r, 2) + pow((CorPeakArea71 - wcoeff * t71), 2) / pow(CorPeakArea71r, 2) + pow((CorPeakArea109 - wcoeff * t109), 2) / pow(CorPeakArea109r, 2) + pow((CorPeakArea138 - wcoeff * t138), 2) / pow(CorPeakArea138r, 2) + pow((CorPeakArea180 - wcoeff * t180), 2) / pow(CorPeakArea180r, 2))/4;
                
                float chim = (c_min320 * 4 + 1) / 4;
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
        
        //starting to find chi square at 0 mixing which will be used to draw angular correlation of 0-2-0 cascade
        
        std::stringstream chi0filename;
        //chi0filename << Output2 << "chi_square" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Xe126.dat";
        chi0filename << datadirfilename << "chi_square" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Xe126.dat";
        std::string chi0name = chi0filename.str();
        
        std::ofstream chi0file (chi0name);
        if (chi0file.is_open())
        {
            /*
             std::stringstream a2a4filename0;
             a2a4filename0 << Output1 << "a2a4_" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Xe126.dat";
             std::string a2a4name0 = a2a4filename0.str();
             
             std::ofstream a2a4file0 (a2a4name0);
             if (a2a4file0.is_open())
             {
             std::stringstream a4filename0;
             a4filename0 << Output1 << "a4_" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Xe126.dat";
             std::string a4name0 = a4filename0.str();
             
             std::ofstream a4file0 (a4name0);
             if (a4file0.is_open())
             {
             std::stringstream a2filename0;
             a2filename0 << Output1 << "a2_" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Xe126.dat";
             std::string a2name0 = a2filename0.str();
             
             std::ofstream a2file0 (a2name0);
             if (a2file0.is_open())
             {
             */
            
            float Wtheo_0 = 0;
            float t42_0 = 0, t71_0 = 0, t109_0 = 0, t138_0 = 0, t180_0 = 0;
            
            for (int theta = 0; theta < 181; theta++)
            {
                
                Wtheo_0 = 1 + a22 * (3 * pow(cos (theta*PI/180), 2) - 1) / 2 + a44 * (35 * pow(cos (theta*PI/180), 4) - 30 * pow(cos (theta*PI/180), 2) + 3) / 8;
                
                if (theta == 42)
                {
                    t42_0 = Wtheo_0;
                }
                if (theta == 71)
                {
                    t71_0 = Wtheo_0;
                }
                if (theta == 109)
                {
                    t109_0 = Wtheo_0;
                }
                if (theta == 138)
                {
                    t138_0 = Wtheo_0;
                }
                if (theta == 180)
                {
                    t180_0 = Wtheo_0;
                }
            }
            
            float wcoeff_0 = 0;
            
            wcoeff_0 = (CorPeakArea42 * t42_0 / pow(CorPeakArea42r, 2) + CorPeakArea71 * t71_0 / pow(CorPeakArea71r, 2) + CorPeakArea109 * t109_0 / pow(CorPeakArea109r, 2) + CorPeakArea138 * t138_0 / pow(CorPeakArea138r, 2) + CorPeakArea180 * t180_0 / pow(CorPeakArea180r, 2)) / (pow(t42_0, 2) / pow(CorPeakArea42r, 2) + pow(t71_0, 2) / pow(CorPeakArea71r, 2) + pow(t109_0, 2) / pow(CorPeakArea109r, 2) + pow(t138_0, 2) / pow(CorPeakArea138r, 2) + pow(t180_0, 2) / pow(CorPeakArea180r, 2));
            
            ChiSquare_0 = (pow((CorPeakArea42 - wcoeff_0 * t42_0), 2) / pow(CorPeakArea42r, 2) + pow((CorPeakArea71 - wcoeff_0 * t71_0), 2) / pow(CorPeakArea71r, 2) + pow((CorPeakArea109 - wcoeff_0 * t109_0), 2) / pow(CorPeakArea109r, 2) + pow((CorPeakArea138 - wcoeff_0 * t138_0), 2) / pow(CorPeakArea138r, 2) + pow((CorPeakArea180 - wcoeff_0 * t180_0), 2) / pow(CorPeakArea180r, 2))/4;
            
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
            /*
             a2file0 << "0" << "  " << a22 << std::endl;
             a4file0 << "0" << "  " << a44 << std::endl;
             a2a4file0 << a22 << "  " << a44 << std::endl;
             
             }
             a2file0.close();
             }
             a4file0.close();
             }
             a2a4file0.close();
             */
        }
        chi0file.close();
        
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
        conflev95filename << datadirfilename << "confidencelev95_8pi_Xe126.dat";
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
        
        float a2 = Q2 * ((R2LLJ2J1 + 2 * del_min * R2LMJ2J1 + pow(del_min, 2) * R2MMJ2J1) * (R2LLJ2J3 + 2 * delta2 * R2LMJ2J3 + pow(delta2, 2) * R2MMJ2J3) * 1 / (1 + pow(del_min, 2)) * 1 / (1 + pow(delta2, 2)));
        float a4 = Q4 * ((R4LLJ2J1 + 2 * del_min * R4LMJ2J1 + pow(del_min, 2) * R4MMJ2J1) * (R4LLJ2J3 + 2 * delta2 * R4LMJ2J3 + pow(delta2, 2) * R4MMJ2J3) * 1 / (1 + pow(del_min, 2)) * 1 / (1 + pow(delta2, 2)));
        
        
        float Wtheo = 0;
        float t42 = 0, t71 = 0, t109 = 0, t138 = 0, t180 = 0;
        float Wtheo_0 = 0;
        float t42_0 = 0, t71_0 = 0, t109_0 = 0, t138_0 = 0, t180_0 = 0;
        
        for (int theta = 0; theta < 181; theta++)
        {
            
            Wtheo = 1 + a2 * (3 * pow(cos (theta*PI/180), 2) - 1) / 2 + a4 * (35 * pow(cos (theta*PI/180), 4) - 30 * pow(cos (theta*PI/180), 2) + 3) / 8;
            
            Wtheo_0 = 1 + a22 * (3 * pow(cos (theta*PI/180), 2) - 1) / 2 + a44 * (35 * pow(cos (theta*PI/180), 4) - 30 * pow(cos (theta*PI/180), 2) + 3) / 8;
            
            if (theta == 42)
            {
                t42 = Wtheo;
                t42_0 = Wtheo_0;
            }
            if (theta == 71)
            {
                t71 = Wtheo;
                t71_0 = Wtheo_0;
            }
            if (theta == 109)
            {
                t109 = Wtheo;
                t109_0 = Wtheo_0;
            }
            if (theta == 138)
            {
                t138 = Wtheo;
                t138_0 = Wtheo_0;
            }
            if (theta == 180)
            {
                t180 = Wtheo;
                t180_0 = Wtheo_0;
            }
        }
        
        float wcoeff = 0;
        
        wcoeff = (CorPeakArea42 * t42 / pow(CorPeakArea42r, 2) + CorPeakArea71 * t71 / pow(CorPeakArea71r, 2) + CorPeakArea109 * t109 / pow(CorPeakArea109r, 2) + CorPeakArea138 * t138 / pow(CorPeakArea138r, 2) + CorPeakArea180 * t180 / pow(CorPeakArea180r, 2)) / (pow(t42, 2) / pow(CorPeakArea42r, 2) + pow(t71, 2) / pow(CorPeakArea71r, 2) + pow(t109, 2) / pow(CorPeakArea109r, 2) + pow(t138, 2) / pow(CorPeakArea138r, 2) + pow(t180, 2) / pow(CorPeakArea180r, 2));
        
        float wcoeff_0 = 0;
        
        wcoeff_0 = (CorPeakArea42 * t42_0 / pow(CorPeakArea42r, 2) + CorPeakArea71 * t71_0 / pow(CorPeakArea71r, 2) + CorPeakArea109 * t109_0 / pow(CorPeakArea109r, 2) + CorPeakArea138 * t138_0 / pow(CorPeakArea138r, 2) + CorPeakArea180 * t180_0 / pow(CorPeakArea180r, 2)) / (pow(t42_0, 2) / pow(CorPeakArea42r, 2) + pow(t71_0, 2) / pow(CorPeakArea71r, 2) + pow(t109_0, 2) / pow(CorPeakArea109r, 2) + pow(t138_0, 2) / pow(CorPeakArea138r, 2) + pow(t180_0, 2) / pow(CorPeakArea180r, 2));
        
        std::stringstream mixedfilename;
        //mixedfilename << Output2 << "AC" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_min_Xe126.dat";
        mixedfilename << datadirfilename << "AC" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_min_Xe126.dat";
        std::string mixedname = mixedfilename.str();
        std::ofstream mixedfile (mixedname);
        
        std::stringstream zerofilename;
        //zerofilename << Output2 << "AC" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Xe126.dat";
        zerofilename << datadirfilename << "AC" << peak << "_"<< J1 << "-" << J2 << "-" << J3 << "_delta_0_Xe126.dat";
        std::string zeroname = zerofilename.str();
        std::ofstream zerofile (zeroname);
        if (zerofile.is_open())
        {
            if (mixedfile.is_open())
            {
                for (int tt = 0; tt < 181; tt++)
                {
                    Wtheomixed = (1 + a2 * (3 * pow(cos (tt*PI/180), 2) - 1) / 2 + a4 * (35 * pow(cos (tt*PI/180), 4) - 30 * pow(cos (tt*PI/180), 2) + 3) / 8) * wcoeff;
                    Wtheomixed_0 = (1 + a22 * (3 * pow(cos (tt*PI/180), 2) - 1) / 2 + a44 * (35 * pow(cos (tt*PI/180), 4) - 30 * pow(cos (tt*PI/180), 2) + 3) / 8) * wcoeff_0;
                    
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
        
    }
    
    std::cout << "MAIN" << std::endl;
    std::cout << "Numerical data will be saved in " << peak << "engout.txt." << std::endl;
    std::cout << "Peak area in " << peak << "peakarea.txt." << std::endl;
    std::cout << "ChiSq in chi_square" << peak << "_*_Xe126.dat" << std::endl;
    std::cout << "Theo AC in ACJ1-J2-J3_delta_min_Xe126.dat" << std::endl;
    std::cout << "maximum AC is " << ac_max << std::endl;
    std::cout << "minimum AC is " << ac_min << std::endl;
    std::cout << "maximum ChiSq is " << chi_max << std::endl;
    std::cout << "minimum ChiSq is " << chi_min << std::endl;
    std::cout << "Delta value at minimum chi square is " << delta_min << std::endl;
    std::cout << "Delta value at maximum chi square is " << delta_max << std::endl;
    std::cout << "deltamin120 is " << delmin120 << std::endl;
    std::cout << "deltamin220 is " << delmin220 << std::endl;
    std::cout << "deltamin320 is " << delmin320 << std::endl;
    std::cout << "c_min020 is " << c_min020 << std::endl;
    std::cout << "c_min120 is " << c_min120 << std::endl;
    std::cout << "c_min220 is " << c_min220 << std::endl;
    std::cout << "c_min320 is " << c_min320 << std::endl;
    std::cout << "c_min420 is " << c_min420 << std::endl;
    
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
    }else if (peak == 884){
        ymin= 1.95 * ACmin;
        ymax= 0.8 * ACmax;
    }else if (peak == 1126){
        ymin= 1.8 * ACmin;
        ymax= 0.65 * ACmax;
    }else if (peak == 1421){
        ymin= 0.97 * exp_min;
        ymax= 1.2 * exp_max;
    }else if (peak == 1505){
        ymin= 0.65 * exp_min;
        ymax= 2.2 * exp_max;
    }else if (peak == 1630){
        ymin= 0.8 * exp_min;
        ymax= 1.5 * exp_max;
    }else if (peak == 1698){
        ymin= 0.9 * exp_min;
        ymax= 1.4 * exp_max;
    }else if (peak == 1822){
        ymin= 0.9 * exp_min;
        ymax= 1.4 * exp_max;
    }else if (peak == 1904){
        ymin= 0.92 * exp_min;
        ymax= 1.3 * exp_max;
    }else if (peak == 2004){
        ymin= 0.84 * exp_min;
        ymax= 1.5 * exp_max;
    }else if (peak == 2130){
        ymin= 0.95 * exp_min;
        ymax= 1.2 * exp_max;
    }else if (peak == 2212){
        ymin= 0.9 * exp_min;
        ymax= 1.3 * exp_max;
    }else if (peak == 2318){
        ymin= 0.9 * exp_min;
        ymax= 1.3 * exp_max;
    }else if (peak == 2421){
        ymin= 0.75 * exp_min;
        ymax= 2.1 * exp_max;
    }else if (peak == 2445){
        ymin= 0.85 * exp_min;
        ymax= 1.4 * exp_max;
    }else if (peak == 2537){
        ymin= 0.85 * exp_min;
        ymax= 1.4 * exp_max;
    }else if (peak == 2658){
        ymin= 0.85 * exp_min;
        ymax= 1.4 * exp_max;
    }else if (peak == 2747){
        ymin= 0.78 * exp_min;
        ymax= 1.7 * exp_max;
    }else if (peak == 2809){
        ymin= 0.92 * exp_min;
        ymax= 1.38 * exp_max;
    }else if (peak == 2892){
        ymin= 0.85 * exp_min;
        ymax= 1.4 * exp_max;
    }else if (peak == 603){
        ymin= 0.85 * exp_min;
        ymax= 1.6 * exp_max;
    }else if (peak == 686){
        ymin= 0.65 * exp_min;
        ymax= 2.1 * exp_max;
    }else if (peak == 744){
        ymin= 0.93 * exp_min;
        ymax= 1.2 * exp_max;
    }else if (peak == 774){
        ymin= 0.7 * exp_min;
        ymax= 1.8 * exp_max;
    }else if (peak == 1085){
        ymin= 0.93 * exp_min;
        ymax= 1.3 * exp_max;
    }else if (peak == 1186){
        ymin= 0.2 * exp_min;
        ymax= 1.55 * exp_max;
    }else{
        ymin= 0.75 * ACmin;
        ymax= 1.07 * ACmax;
    }
    
    std::stringstream bfile;
    //bfile << Output2 << "ACscript" << peak << ".bfile"; //script file includes xmgrace Angular Correlation plots
    bfile << datadirfilename << "ACscript" << peak << ".bfile"; //script file includes xmgrace Angular Correlation plots
    std::string bfilename = bfile.str();
    std::ofstream bfilebash(bfilename);
    
    if (bfilebash.is_open())
    {
        bfilebash << "read xy \"AC" << peak << "_0-2-0_delta_0_Xe126.dat\"" << std::endl;
        bfilebash << "read xy \"AC" << peak << "_1-2-0_delta_min_Xe126.dat\"" << std::endl;
        bfilebash << "read xy \"AC" << peak << "_2-2-0_delta_min_Xe126.dat\"" << std::endl;
        bfilebash << "read xy \"AC" << peak << "_3-2-0_delta_min_Xe126.dat\"" << std::endl;
        bfilebash << "read xy \"AC" << peak << "_4-2-0_delta_0_Xe126.dat\"" << std::endl;
        bfilebash << "read xydy \"" << peak << "exp.dat\"" << std::endl;
        
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
        bfilebash << "s5 legend \"exp " << peak << " keV\"" << std::endl;
        if (peak == 1073){
            bfilebash << "legend 0.55, 0.77" << std::endl;
        }else if (peak == 1421||peak==1630||peak==1698||peak == 2445){
            bfilebash << "legend 0.50, 0.77" << std::endl;
        }else if (peak==1822||peak == 1904||peak==2004||peak==2130||peak==2421||peak == 2537||peak == 2658||peak == 2747||peak == 2809||peak == 2892||peak == 603||peak == 686){
            bfilebash << "legend 0.49, 0.77" << std::endl;
        }else{
            bfilebash << "legend 0.51, 0.77" << std::endl;
        }
        bfilebash << "legend char size 1.2" << std::endl;
        
        if (peak == 818||peak==306){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 1.00, 0.65" << std::endl;
            bfilebash << "string char size 1.8" << std::endl;
            bfilebash << "string def \"2\\s2\\N\\S+\"" << std::endl;
        }else if (peak == 1073||peak==873){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 1.00, 0.65" << std::endl;
            bfilebash << "string char size 1.8" << std::endl;
            bfilebash << "string def \"0\\s3\\N\\S+\"" << std::endl;
        }else if (peak == 1126||peak==873){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.95, 0.65" << std::endl;
            bfilebash << "string char size 1.8" << std::endl;
            bfilebash << "string def \"2\\s3\\N\\S+\"" << std::endl;
        }else if (peak == 1421||peak==873){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.90, 0.65" << std::endl;
            bfilebash << "string char size 1.8" << std::endl;
            bfilebash << "string def \"2\\S+\\N(3\\S-\\N)\"" << std::endl;
        }else if (peak == 1505||peak==873){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 1.00, 0.65" << std::endl;
            bfilebash << "string char size 1.8" << std::endl;
            bfilebash << "string def \"(1,3)\\S+\\N(3\\S+\\N)\"" << std::endl;
        }else if (peak == 884){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.9, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(2-4)\\S+\\N(4\\S+\\N)\"" << std::endl;
        }else if (peak == 1630){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.9, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(1-4)\\S+\\N(2\\S+\\N)\"" << std::endl;
        }else if (peak == 1698){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.9, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(1-4)\\S+\\N(2\\S+\\N)\"" << std::endl;
        }else if (peak==1822){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.9, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(1,3)\\S+\\N(2\\S+\\N)\"" << std::endl;
        }else if (peak == 1904){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.9, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(2-4)\\S+\\N(4\\S+\\N)\"" << std::endl;
        }else if (peak == 2004){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(1-3)\\S+\\N(0\\S+\\N)\"" << std::endl;
        }else if (peak == 2130){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(2-4)\\S+\\N(2\\S+\\N)\"" << std::endl;
        }else if (peak == 2212){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(2,3)\\S+\\N(2\\S+\\N)\"" << std::endl;
        }else if (peak == 2318){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(2,3)\\S+\\N(2\\S+\\N)\"" << std::endl;
        }else if (peak == 2421){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(1,3)\\S+\\N(1\\S+\\N)\"" << std::endl;
        }else if (peak == 2445){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.95, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"2\\S+\"" << std::endl;
        }else if (peak == 2537){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(1,2)\\S+\\N(2\\S+\\N)\"" << std::endl;
        }else if (peak == 2658){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(2,3)\\S+\\N(2\\S+\\N)\"" << std::endl;
        }else if (peak == 2747){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(1-3)\\S+\\N(2\\S+\\N)\"" << std::endl;
        }else if (peak == 2809){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(2,3)\\S+\\N(2\\S+\\N)\"" << std::endl;
        }else if (peak == 2892){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.95, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(2)\\S+\\N(4\\S+\\N)\"" << std::endl;
        }else if (peak == 603){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(1-4)\\S+\\N(0\\S+\\N)\"" << std::endl;
        }else if (peak == 686){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(1,3)\\S+\\N(3\\S+\\N)\"" << std::endl;
        }else if (peak == 744){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.91, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(2-4)\\S+\\N(4\\S+\\N)\"" << std::endl;
        }else if (peak == 774){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(1-4)\\S+\\N(4\\S+\\N)\"" << std::endl;
        }else if (peak == 1085){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.89, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"(1-4)\\S+\\N(4\\S+\\N)\"" << std::endl;
        }else if (peak == 1186){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.95, 0.65" << std::endl;
            bfilebash << "string char size 1.6" << std::endl;
            bfilebash << "string def \"0\\S+\"" << std::endl;
        }
        if (peak ==1073||peak==884||peak == 1126||peak==818||peak==1421||peak==1505||peak==1630||peak==1698||peak==1822||peak== 1904||peak==2004||peak==2130||peak == 2212||peak == 2318||peak == 2421||peak == 2445||peak == 2537||peak == 2658||peak == 2747||peak == 2809||peak == 2892){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.48, 0.44" << std::endl;
            bfilebash << "string char size 1.3" << std::endl;
            bfilebash << "string def \"gated on 658 keV \\xg\\0 ray (Xe-126)\"" << std::endl;
        }else if (peak == 603||peak == 686||peak == 744||peak == 774||peak == 1085||peak == 1186){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.45, 0.44" << std::endl;
            bfilebash << "string char size 1.3" << std::endl;
            bfilebash << "string def \"gated on 1476 keV \\xg\\0 ray (Xe-126)\"" << std::endl;
        }else if (peak == 2678){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.55, 0.43" << std::endl;
            bfilebash << "string char size 1.3" << std::endl;
            bfilebash << "string def \"gated on 331 keV \\xg\\0 ray\"" << std::endl;
        }else if (peak == 306||peak==873||peak==945||peak==1377||peak==1461||peak==1652||peak==1997||peak==1799||peak==371||peak==560||peak==652||peak==1039||peak == 1421||peak==1486||peak == 1811||peak == 1836||peak == 1886||peak==1223||peak == 1736){
            bfilebash << "with string" << std::endl;
            bfilebash << "string on" << std::endl;
            bfilebash << "string loctype view" << std::endl;
            bfilebash << "string 0.55, 0.44" << std::endl;
            bfilebash << "string char size 1.3" << std::endl;
            bfilebash << "string def \"gated on 843 keV \\xg\\0 ray\"" << std::endl;
        }
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
        
        if (exp_max > 75000){
            bfilebash << "yaxis tick major 10000" << std::endl;
            bfilebash << "yaxis tick minor 5000" << std::endl;
            bfilebash << "view 0.17, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max > 50000 && exp_max < 75000){
            bfilebash << "yaxis tick major 10000" << std::endl;
            bfilebash << "yaxis tick minor 5000" << std::endl;
            bfilebash << "view 0.17, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max > 15000 && exp_max < 50000){
            bfilebash << "yaxis tick major 5000" << std::endl;
            bfilebash << "yaxis tick minor 1000" << std::endl;
            bfilebash << "view 0.17, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max > 8000 && exp_max <= 15000){
            bfilebash << "yaxis tick major 2000" << std::endl;
            bfilebash << "yaxis tick minor 500" << std::endl;
            bfilebash << "view 0.17, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max > 4000 && exp_max <= 8000){
            bfilebash << "yaxis tick major 1000" << std::endl;
            bfilebash << "yaxis tick minor 200" << std::endl;
            bfilebash << "view 0.17, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max > 2000 && exp_max <= 4000){
            bfilebash << "yaxis tick major 500" << std::endl;
            bfilebash << "yaxis tick minor 100" << std::endl;
            bfilebash << "view 0.17, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max >= 1000 && exp_max <= 2000){
            bfilebash << "yaxis tick major 200" << std::endl;
            bfilebash << "yaxis tick minor 50" << std::endl;
            bfilebash << "view 0.17, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max > 500 && exp_max < 1000){
            bfilebash << "yaxis tick major 100" << std::endl;
            bfilebash << "yaxis tick minor 50" << std::endl;
            bfilebash << "view 0.14, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max > 250 && exp_max <= 500){
            bfilebash << "yaxis tick major 50" << std::endl;
            bfilebash << "yaxis tick minor 10" << std::endl;
            bfilebash << "view 0.14, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max >= 150 && exp_max <= 250){
            bfilebash << "yaxis tick major 40" << std::endl;
            bfilebash << "yaxis tick minor 20" << std::endl;
            bfilebash << "view 0.14, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max >= 100 && exp_max < 150){
            bfilebash << "yaxis tick major 40" << std::endl;
            bfilebash << "yaxis tick minor 20" << std::endl;
            bfilebash << "view 0.14, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max >= 60 && exp_max < 100){
            bfilebash << "yaxis tick major 10" << std::endl;
            bfilebash << "yaxis tick minor 5" << std::endl;
            bfilebash << "view 0.11, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max >= 40 && exp_max < 60){
            bfilebash << "yaxis tick major 10" << std::endl;
            bfilebash << "yaxis tick minor 2" << std::endl;
            bfilebash << "view 0.11, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max >= 20 && exp_max < 40){
            bfilebash << "yaxis tick major 5" << std::endl;
            bfilebash << "yaxis tick minor 1" << std::endl;
            bfilebash << "view 0.11, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max >= 10 && exp_max < 20){
            bfilebash << "yaxis tick major 2" << std::endl;
            bfilebash << "yaxis tick minor 1" << std::endl;
            bfilebash << "view 0.11, 0.13, 1.26, 0.8" << std::endl;
        }else if (exp_max >= 2 && exp_max < 10){
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
        bfilebash << "saveall \"AC_"<< peak << ".arg\"" << std::endl;
        bfilebash << "print to \"AC_"<< peak << ".eps\"" << std::endl;
        bfilebash << "device \"EPS\" OP \"level2\"" << std::endl;
        bfilebash << "print" << std::endl;
        
        bfilebash.close();
    }
    else
    {
        std::cout << "error: the file is not open";
    }
    
    float yminchi=0;
    float ymaxchi=0;
    if (peak == 818 || peak == 2324 || peak == 2683 || peak == 1811){
        yminchi= 0.098 * chi_min;
        ymaxchi= 4 * chi_max;
    }else if (peak == 884 || peak == 2348 || peak == 2399){
        yminchi= 0.1 * chi_min;
        ymaxchi= 5.0 * chi_max;
    }else if (peak == 1126 || peak == 1421){
        yminchi= 0.3 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 1505){
        yminchi= 0.3 * chi_min;
        ymaxchi= 6.0 * chi_max;
    }else if (peak == 1630){
        yminchi= 0.05 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 1698){
        yminchi= 0.3 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 1822){
        yminchi= 0.3 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 1904){
        yminchi= 0.3 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 2004){
        yminchi= 0.3 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 2130){
        yminchi= 0.3 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 2212){
        yminchi= 0.3 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 2318){
        yminchi= 0.3 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 2421){
        yminchi= 0.4 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 2445){
        yminchi= 0.1 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 2658){
        yminchi= 0.1 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 2747){
        yminchi= 0.1 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 2809){
        yminchi= 0.3 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 2892){
        yminchi= 0.1 * chi_min;
        ymaxchi= 2.5 * chi_max;
    }else if (peak == 603){
        yminchi= 0.2 * chi_min;
        ymaxchi= 2.5 * chi_max;
    }else if (peak == 686){
        yminchi= 0.4 * chi_min;
        ymaxchi= 4.0 * chi_max;
    }else if (peak == 744){
        yminchi= 0.4 * chi_min;
        ymaxchi= 4.5 * chi_max;
    }else if (peak == 774){
        yminchi= 0.5 * chi_min;
        ymaxchi= 4.5 * chi_max;
    }else if (peak == 1085){
        yminchi= 0.5 * chi_min;
        ymaxchi= 4.5 * chi_max;
    }else if (peak == 1186){
        yminchi= 0.6 * chi_min;
        ymaxchi= 1.5 * chi_max;
    }else{
        //yminchi= 0.5 * chi_min;
        yminchi= 0.15 * chi_min;
        ymaxchi= 2.0 * chi_max;
    }
    
    std::stringstream Chi;
    //Chi << Output2 << "Chiscript" << peak << ".bfile";//script file includes xmgrace Chi Square plots
    Chi << datadirfilename << "Chiscript" << peak << ".bfile";//script file includes xmgrace Chi Square plots
    std::string bfileChi = Chi.str();
    std::ofstream bfilebashChi(bfileChi);
    
    if (bfilebashChi.is_open())
    {
        bfilebashChi << "read xy \"chi_square" << peak << "_0-2-0_delta_0_Xe126.dat\"" << std::endl;
        bfilebashChi << "read xy \"chi_square" << peak << "_1-2-0_Xe126.dat\"" << std::endl;
        bfilebashChi << "read xy \"chi_square" << peak << "_2-2-0_Xe126.dat\"" << std::endl;
        bfilebashChi << "read xy \"chi_square" << peak << "_3-2-0_Xe126.dat\"" << std::endl;
        bfilebashChi << "read xy \"chi_square" << peak << "_4-2-0_delta_0_Xe126.dat\"" << std::endl;
        bfilebashChi << "read xy \"confidencelev95_8pi_Xe126.dat\"" << std::endl;
        
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
        
        if (peak == 884 || peak == 1698 ||peak == 1822||peak == 1904||peak == 2130||peak == 2212||peak == 2658||peak == 2809||peak==744||peak == 1085) {
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        }else if (peak == 1186) {
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        }else if (peak == 1630||peak==652||peak == 2004||peak == 2445||peak == 2537||peak == 2747||peak == 2892||peak == 603||peak == 774) {
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        }else if (peak == 1514 || peak == 1933 || peak == 2604||peak == 1836||peak==1716){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
        }else if (peak == 883 || peak == 1126 || peak == 1421||peak == 2318){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
        }else if (peak == 2399){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        }else if (peak==1505){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
        }else if (peak==2421||peak == 686){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
        }else if (peak == 1385||peak==2199){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
        }else if (peak == 1457 || peak == 1380 ){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        }else if (peak == 2324 || peak == 2780 || peak == 2848||peak==2839){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
        }else if (peak == 873 || peak == 945|| peak == 1891 || peak == 1973|| peak==2162 ){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
        }else if (peak == 1073 || peak == 0000){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
        }else if (peak == 818){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
        }else if (peak==2962){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
        }else if (peak == 1799){
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 1))/pow(10, 1) << ")\"" << std::endl;
        }else{
            bfilebashChi << "s0 legend \"0-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min020*pow(10, 2))/pow(10, 2) << ")\"" << std::endl;
            bfilebashChi << "s1 legend \"1-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min120*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s2 legend \"2-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min220*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s3 legend \"3-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min320*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
            bfilebashChi << "s4 legend \"4-2-0 (\\xc\\0\\S2\\N/\\xn=" << round(c_min420*pow(10, 0))/pow(10, 0) << ")\"" << std::endl;
        }
        if (peak == 497 || peak == 884 || peak == 1126  || peak == 1421||peak==1630||peak == 2747||peak == 2809){
            bfilebashChi << "legend 0.25, 0.43" << std::endl;
        }else if (peak==1822||peak==1904||peak == 2004||peak == 2130||peak == 2212||peak == 2318||peak == 2445||peak == 2537||peak == 2658||peak==744||peak == 1085){
            bfilebashChi << "legend 0.20, 0.43" << std::endl;
        }else if (peak==1186){
            bfilebashChi << "legend 0.30, 0.43" << std::endl;
        }else if (peak==1505){
            bfilebashChi << "legend 0.89, 0.43" << std::endl;
        }else if (atan (delta_min) * 180 / PI > 30) {
            bfilebashChi << "legend 0.20, 0.43" << std::endl;
        }else{
            bfilebashChi << "legend 0.90, 0.43" << std::endl;
        }
        bfilebashChi << "legend char size 1.2" << std::endl;
        
        bfilebashChi << "with string" << std::endl;
        bfilebashChi << "string on" << std::endl;
        bfilebashChi << "string loctype view" << std::endl;
        if (peak == 2010 || peak == 2324 || peak == 2399 || peak == 2476 || peak == 2683|| peak == 1514) {
            bfilebashChi << "string 0.65, 0.65" << std::endl;
        }else if (peak == 1505){
            bfilebashChi << "string 0.5, 0.2" << std::endl;
        }else if (peak == 1164|| peak == 1933||peak==1486){
            bfilebashChi << "string 0.65, 0.2" << std::endl;
        }else if (peak == 1822||peak == 2004||peak == 603||peak == 774){
            bfilebashChi << "string 0.75, 0.2" << std::endl;
        }else if (peak == 884||peak == 1126 || peak == 1421||peak==1630||peak == 1698||peak == 1904||peak == 2130||peak == 2212||peak == 2318||peak == 2445||peak == 2537||peak == 2658||peak == 2747||peak == 2809||peak==744||peak == 1085||peak==1186){
            bfilebashChi << "string 0.95, 0.2" << std::endl;
        }else if (peak == 1223){
            bfilebashChi << "string 0.45, 0.2" << std::endl;
        }else if (peak ==2180){
            bfilebashChi << "string 0.3, 0.2" << std::endl;
        }else{
            bfilebashChi << "string 0.35, 0.2" << std::endl;
        }
        bfilebashChi << "string char size 1.3" << std::endl;
        bfilebashChi << "string def \"" << peak << " keV \\xg\\0 ray\"" << std::endl;
        
	bfilebashChi << "with string" << std::endl;
        bfilebashChi << "string on" << std::endl;
        bfilebashChi << "string loctype view" << std::endl;
        bfilebashChi << "string 0.1, 0.23" << std::endl;
        bfilebashChi << "string char size 1.3" << std::endl;
        bfilebashChi << "string def \"95\%\"" << std::endl;
        
        if (peak == 1073 || peak == 884||peak==818||peak==1126||peak == 1421||peak == 1505||peak == 1630||peak == 1698||peak == 1822||peak == 1904||peak == 2004||peak == 2130||peak == 2212||peak == 2318||peak==2421||peak == 2445||peak == 2537||peak == 2658||peak == 2747||peak == 2809||peak == 2892){
            bfilebashChi << "with string" << std::endl;
            bfilebashChi << "string on" << std::endl;
            bfilebashChi << "string loctype view" << std::endl;
            bfilebashChi << "string 0.5, 0.74" << std::endl;
            bfilebashChi << "string char size 1.3" << std::endl;
            bfilebashChi << "string def \"gated on 658 keV \\xg\\0 ray (Xe-126)\"" << std::endl;
        }else if (peak == 603 || peak == 686||peak==744||peak == 774||peak == 1085||peak==1186){
            bfilebashChi << "with string" << std::endl;
            bfilebashChi << "string on" << std::endl;
            bfilebashChi << "string loctype view" << std::endl;
            bfilebashChi << "string 0.50, 0.74" << std::endl;
            bfilebashChi << "string char size 1.3" << std::endl;
            bfilebashChi << "string def \"gated on 1476 keV \\xg\\0 ray (Xe-126)\"" << std::endl;
        }else if (peak == 883){
            bfilebashChi << "with string" << std::endl;
            bfilebashChi << "string on" << std::endl;
            bfilebashChi << "string loctype view" << std::endl;
            bfilebashChi << "string 0.7, 0.2" << std::endl;
            bfilebashChi << "string char size 1.3" << std::endl;
            bfilebashChi << "string def \"gated on 331 keV \\xg\\0 ray\"" << std::endl;
        }else if (peak==2962||peak == 2323||peak==1487){
            bfilebashChi << "with string" << std::endl;
            bfilebashChi << "string on" << std::endl;
            bfilebashChi << "string loctype view" << std::endl;
            bfilebashChi << "string 0.56, 0.68" << std::endl;
            bfilebashChi << "string char size 1.3" << std::endl;
            bfilebashChi << "string def \"gated on 331 keV \\xg\\0 ray\"" << std::endl;
        }
        
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
        bfilebashChi << "saveall \"Chi_"<< peak << ".arg\"" << std::endl;
        bfilebashChi << "print to \"Chi_"<< peak << ".eps\"" << std::endl;
        bfilebashChi << "device \"EPS\" OP \"level2\"" << std::endl;
        bfilebashChi << "print" << std::endl;
        
        bfilebashChi.close();
    }
    else
    {
        std::cout << "error: the file is not open";
    }
    
    std::stringstream chib;
    //chib << Output2<< "AChi" << peak << "bash.txt";
    chib << datadirfilename << "AChi" << peak << "bash.txt";
    std::string bfileb = chib.str();
    std::ofstream bbashAChi(bfileb);
    
    if (bbashAChi.is_open())
    {
        bbashAChi << "xmgrace -batch Chiscript" << peak << ".bfile -nosafe -hardcopy" << std::endl;
        bbashAChi << "xmgrace -batch ACscript" << peak << ".bfile -nosafe -hardcopy" << std::endl;
        bbashAChi << "ps2pdf" << " AC_" << peak << ".eps" << std::endl;
        bbashAChi << "pdfcrop" << " AC_" << peak << ".pdf" << std::endl;
        bbashAChi << "ps2pdf" << " Chi_" << peak << ".eps" << std::endl;
        bbashAChi << "pdfcrop" << " Chi_" << peak << ".pdf" << std::endl;
//        bbashAChi << "cp" << " Chi_" << peak << "-crop.pdf ../Xe126plots2-0/" << std::endl;
//        bbashAChi << "cp" << " AC_" << peak << "-crop.pdf ../Xe126plots2-0/" << std::endl;
//        bbashAChi << "xpdf" << " ../Xe126plots2-0/Chi_" << peak << "Xe126-crop.pdf" << std::endl;
//        bbashAChi << "xpdf" << " ../Xe126plots2-0/AC_" << peak << "Xe126-crop.pdf" << std::endl;
        bbashAChi << "xpdf" << " Chi_" << peak << "-crop.pdf" << std::endl;
        bbashAChi << "xpdf" << " AC_" << peak << "-crop.pdf" << std::endl;
        bbashAChi.close();
    }
    else
    {
        std::cout << "error: the AChi" << peak << "bash.txt file is not open";
    }
    return 0;
    
}
