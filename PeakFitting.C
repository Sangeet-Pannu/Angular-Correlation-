void PeakFitting()
{
	auto c1 = new TCanvas("c1", "CANVAS", 500, 500); 
	bool status = true;
	char user;
	double cent;
	double lowProj,highProj;
	while(status)
	{
		asE->Draw();
		cout << "Centriod to be Fit: ";
		cin >> cent;

		cout << "Low limit for fit: ";
		cin >> lowProj;

		cout << "high limit for fit: ";
		cin >> highProj;

		TABPeak *p1 = new TABPeak(cent);
		TABPeak *p2 = new TABPeak(7732);
		TABPeak *p3 = new TABPeak(7743);
		TRWPeak *p4 = new TRWPeak(7280);
		TPeakFitter *pf = new TPeakFitter(lowProj,highProj);

		pf->AddPeak(p1);
		//pf->AddPeak(p2);
	 	//pf->AddPeak(p3);
		//pf->AddPeak(p4);
		//pf->GetBackground()->FixParameter(1,0);
		 //p1->GetFitFunction()->FixParameter(5,0);
		// p2->GetFitFunction()->FixParameter(5,0);
		 //p3->GetFitFunction()->FixParameter(5,0);

		 //p1->GetFitFunction()->FixParameter(1,2460.28);
		// p2->GetFitFunction()->FixParameter(1,7432);
		p1->GetFitFunction()->SetParLimits(1,cent-2,cent+2);
		p2->GetFitFunction()->SetParLimits(1,7732-5,7732+5);
		p3->GetFitFunction()->SetParLimits(1,7743-5,7743+5);
		//p4->GetFitFunction()->SetParLimits(1,7280-4,7280+4);

		//p1->GetFitFunction()->SetParLimits(2,2.1,2.7);
		//p2->GetFitFunction()->SetParLimits(2,2.1,2.7);

		pf->Fit(asE,"EM");
		asE->SetAxisRange(lowProj-5,highProj+5);
		c1->Update();

		//cout << "FWHM: " << p1->GetFitFunction()->GetParameter(2)*2.35 << endl;
		cout << "Another fit? (y/n): ";
		cin >> user;

		if(user == 'y')continue;
		else if(user == 'n')break;

		c1->Close();
	}
}
