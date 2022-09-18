//
//  main.cpp
//  ChiXe126_x-2-0
//
//  Created by Sambuu on 2019-07-10.
//  Copyright © 2019 Sambuu. All rights reserved.
//  Modifeid for 126Xe by Farnaz Ghazi Moradi January 2020
//
 
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
