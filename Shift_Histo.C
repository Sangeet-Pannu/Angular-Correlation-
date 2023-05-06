#include "TAxis.h"
#include "TH1.h"
#include "TArrayD.h"

Double_t ScaleX(Double_t x)
{
  Double_t v;
  v = x - 511; // "linear scaling" function example
  return v;
}

Double_t ScaleX2(Double_t x)
{
  Double_t v;
  v = x - (511*2); // "linear scaling" function example
  return v;
}


Double_t ScaleY(Double_t y)
{
  Double_t v;
  v = y; // "linear scaling" function example
  return v;
}

Double_t ScaleZ(Double_t z)
{
  Double_t v;
  v = 30 * z + 300; // "linear scaling" function example
  return v;
}

void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t))
{
  if (!a) return; // just a precaution
  if (a->GetXbins()->GetSize())
    {
      // an axis with variable bins
      // note: bins must remain in increasing order, hence the "Scale"
      // function must be strictly (monotonically) increasing
      TArrayD X(*(a->GetXbins()));
      for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i]);
      a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
    }
  else
    {
      // an axis with fix bins
      // note: we modify Xmin and Xmax only, hence the "Scale" function
      // must be linear (and Xmax must remain greater than Xmin)
      a->Set( a->GetNbins(),
              Scale(a->GetXmin()), // new Xmin
              Scale(a->GetXmax()) ); // new Xmax
    }
  return;
}

void ScaleXaxis(TH1 *h, Double_t (*Scale)(Double_t))
{
  if (!h) return; // just a precaution
  ScaleAxis(h->GetXaxis(), Scale);
  return;
}

void ScaleYaxis(TH1 *h, Double_t (*Scale)(Double_t))
{
  if (!h) return; // just a precaution
  ScaleAxis(h->GetYaxis(), Scale);
  return;
}

void ScaleZaxis(TH1 *h, Double_t (*Scale)(Double_t))
{
  if (!h) return; // just a precaution
  ScaleAxis(h->GetZaxis(), Scale);
  return;
}

void ShiftHisto(TH1 *h1){
  //TCanvas *c1;
  //c1->Draw();
  
   auto legend = new TLegend(0.1,0.7,0.28,0.9);



   legend->Draw();
  h1->SetLineColor(kBlue);
  h1->Draw();
 // c1->Update();
  TH1 *histo = (TH1*) h1->Clone();
  TH1 *histo2 = (TH1*) h1->Clone();

  ScaleXaxis(histo, ScaleX);
  ScaleXaxis(histo2,ScaleX2);
  
  histo->SetLineColor(kRed);
  histo->Draw("same");
  histo2->SetLineColor(kMagenta);
  histo2->Draw("same");

  //c1->Update();;

  
   legend->AddEntry(h1, "Original Histogram", "l");
   legend->AddEntry(histo, "Shifted Histogram by 511 keV", "l");
       legend->AddEntry(histo2, "Shifted Histogram by 1022 keV", "l"); 
   legend->Draw();

}
