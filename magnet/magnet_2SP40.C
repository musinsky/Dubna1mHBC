// Author: Jan Musinsky
// 15/05/2019

#include <TGraph.h>
#include <TAxis.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TStyle.h>

void magnet_2SP40()
{
  // 2005-01 measurement with DH10 Hall sensor
  // other Hall sensor than 2004-10 measurement

  const Int_t np = 20;
  Double_t B[np] = { 0.065, 0.09, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
                     0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 };
  Double_t Ux[np] = { 6.28, 8.69, 14.45, 19.27, 24.06, 28.86, 33.65, 38.44,
                      43.22, 48.00, 52.77, 57.55, 62.30, 67.05, 71.80, 76.55,
                      81.29, 86.04, 90.75, 95.49 };

  TGraph *gra = new TGraph(np, Ux, B);
  gra->SetTitle(0);
  gra->SetMarkerStyle(20);
  gra->SetMarkerSize(0.75);
  gra->GetYaxis()->CenterTitle();
  gra->GetYaxis()->SetTitle("B (Tesla)");
  gra->GetXaxis()->SetTitle("U_{x} (mV)");
  gra->Draw("AP");

  gra->Fit("pol1");
  TF1 *pol1 = gra->GetFunction("pol1");
  pol1->SetTitle("Hall sensor DH10");
  pol1->SetLineWidth(1);
  pol1->SetLineColor(kBlack);
  Double_t fpar0 = pol1->GetParameter(0);
  Double_t fpar1 = pol1->GetParameter(1);
  TLegend *leg = new TLegend(0.20, 0.65, 0.55, 0.75);
  //  leg->SetHeader("magnet  2SP-40");
  //  leg->AddEntry(pol1, "B = 0.01048 U_{x} - 0.00213");
  leg->AddEntry(pol1, TString::Format("B = %7.5f U_{x} - %7.5f", fpar1, TMath::Abs(fpar0)));
  leg->SetTextAlign(22);
  leg->SetLineColor(kGray);
  leg->Draw();

  gStyle->SetStripDecimals(kFALSE);
  gStyle->SetGridColor(kGray);
  gPad->SetGrid();
  gPad->GetCanvas()->SetWindowSize(700, 400);

  const Int_t np0 = 11, np1 = 14, np2 = 14;
  Double_t x[3] = { 50, 75, 84 };
  Double_t z0[np0] = { 76, 57, 38, 24, 10, 5, 0, -5, -10, -20, -30 };
  Double_t z1[np1] = { 76, 66, 56, 46, 36, 26, 16, 6, 0, -5, -10, -20, -30, -40 };
  Double_t z2[np2] = { 76, 66, 56, 46, 36, 26, 16, 6, 0, -5, -10, -20, -30, -40 };
  Double_t Ux0[np0] = { 82.09, 81.96, 81.80, 81.60, 81.22, 81.44, 74.15,
                        47.20, 32.35, 20.84, 14.71 };
  Double_t Ux1[np1] = { 81.51, 81.51, 81.47, 81.46, 81.44, 81.40, 81.28,
                        81.58, 71.48, 44.28, 31.98, 20.53, 14.45, 10.19 };
  Double_t Ux2[np2] = { 81.05, 81.06, 81.11, 81.12, 81.12, 81.10, 80.90,
                        80.72, 68.65, 43.61, 31.39, 20.15, 14.21, 9.92 };

  TGraph *g0 = new TGraph();
  g0->SetMarkerStyle(25);
  g0->SetLineStyle(1);
  TGraph *g1 = new TGraph();
  g1->SetMarkerStyle(20);
  g1->SetLineStyle(2);
  TGraph *g2 = new TGraph();
  g2->SetMarkerStyle(24);
  g2->SetLineStyle(3);

  Double_t zc = 75.0, xc = 50.0;
  for (Int_t j = 0; j < np0; j++) g0->SetPoint(j, zc - z0[j], pol1->Eval(Ux0[j]));
  for (Int_t j = 0; j < np1; j++) g1->SetPoint(j, zc - z1[j], pol1->Eval(Ux1[j]));
  for (Int_t j = 0; j < np2; j++) g2->SetPoint(j, zc - z2[j], pol1->Eval(Ux2[j]));

  new TCanvas();
  g0->SetMinimum(0.0);
  g0->SetMaximum(1.0);
  g0->GetXaxis()->SetTitle("z (cm)");
  g0->GetXaxis()->SetNdivisions(310);
  g0->GetYaxis()->SetTitle("B (Tesla)");
  g0->GetYaxis()->CenterTitle();
  g0->GetYaxis()->SetNdivisions(409);
  g0->Draw("APL");
  g1->Draw("PL");
  g2->Draw("PL");

  TLine line;
  line.SetLineStyle(4);
  line.DrawLine(zc, g0->GetMinimum(), zc, g0->GetMaximum());
  leg = new TLegend(0.25, 0.35, 0.55, 0.50);
  leg->SetHeader("B_{max} = 0.85 Tesla, y = 0 cm");
  leg->AddEntry(g0, TString::Format("x = %2g cm", x[0] - xc), "PL");
  leg->AddEntry(g1, TString::Format("x = %2g cm", x[1] - xc), "PL");
  leg->AddEntry(g2, TString::Format("x = %2g cm", x[2] - xc), "PL");
  leg->SetTextAlign(22);
  leg->SetLineColor(kGray);
  leg->Draw();
}
