// Author: Jan Musinsky
// 08/10/2009

{
  gROOT->SetStyle("Plain");
  gStyle->SetStripDecimals(kFALSE);
  const Int_t n_points[3]={11,14,14};
  Double_t x[3]={50,75,84};
  Double_t z0[n_points[0]]={76,57,38,24,10,5,0,-5,-10,-20,-30};
  Double_t z1[n_points[1]]={76,66,56,46,36,26,16,6,0,-5,-10,-20,-30,-40};
  Double_t z2[n_points[2]]={76,66,56,46,36,26,16,6,0,-5,-10,-20,-30,-40};
  Double_t Ux_m0[n_points[0]]={82.09,81.96,81.80,81.60,81.22,81.44,74.15,47.20,32.35,20.84,14.71};
  Double_t Ux_m1[n_points[1]]={81.51,81.51,81.47,81.46,81.44,81.40,81.28,81.58,71.48,44.28,31.98,20.53,14.45,10.19};
  Double_t Ux_m2[n_points[2]]={81.05,81.06,81.11,81.12,81.12,81.10,80.90,80.72,68.65,43.61,31.39,20.15,14.21,9.92};
  Double_t zcenter = 75.0, xcenter = 50.0;
  Int_t count = 0;
  g0 = new TGraph(); g0->SetMarkerStyle(03); g0->SetLineStyle(1);
  g1 = new TGraph(); g1->SetMarkerStyle(20); g1->SetLineStyle(2);
  g2 = new TGraph(); g2->SetMarkerStyle(24); g2->SetLineStyle(3);
  fun = new TF1("calib_DH10", "0.0104793*x-0.00212844");

  for (Int_t j = 0; j < n_points[0]; j++)
    g0->SetPoint(count++, zcenter-z0[j], fun->Eval(Ux_m0[j]));
  x[0] -= xcenter;
  count = 0;
  for (Int_t j = 0; j < n_points[1]; j++)
    g1->SetPoint(count++, zcenter-z1[j], fun->Eval(Ux_m1[j]));
  x[1] -= xcenter;
  count = 0;
  for (Int_t j = 0; j < n_points[2]; j++)
    g2->SetPoint(count++, zcenter-z2[j], fun->Eval(Ux_m2[j]));
  x[2] -= xcenter;

  new TCanvas();
  g0->SetMinimum(0); g0->SetMaximum(1);
  g0->GetXaxis()->SetTitle("z, cm");
  g0->GetXaxis()->SetNdivisions(310);
  g0->GetYaxis()->CenterTitle();
  g0->GetYaxis()->SetTitle("B, Tesla");
  g0->GetYaxis()->SetNdivisions(409);
  g0->Draw("APL"); g1->Draw("PL"); g2->Draw("PL");
  TLine line;
  line.SetLineStyle(3);
  line.DrawLine(zcenter, g0->GetMinimum(), zcenter, g0->GetMaximum());
  leg = new TLegend(0.25,0.35,0.55,0.50);
  leg->SetHeader("B_{max} = 0.85 Tesla, y = 0 cm");
  leg->AddEntry(g0, Form("x = %2g cm", x[0]), "PL");
  leg->AddEntry(g1, Form("x = %2g cm", x[1]), "PL");
  leg->AddEntry(g2, Form("x = %2g cm", x[2]), "PL");
  leg->SetFillColor(0);
  leg->SetTextAlign(22);
  leg->Draw();
  gROOT->LoadMacro("$ROOTSYS/macros/printTeX.C");
  printTeX("magnet_eff");

  const Int_t np = 20;
  Double_t B[np]={0.065,0.09,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
  Double_t Ux_minus[np]={6.28,8.69,14.45,19.27,24.06,28.86,33.65,38.44,43.22,48.00,52.77,57.55,62.30,67.05,71.80,76.55,81.29,86.04,90.75,95.49};
  grad = new TGraph(np, Ux_minus, B);
  grad->SetTitle(0);
  new TCanvas();
  gPad->SetGridx(); gPad->SetGridy();
  grad->SetMarkerStyle(8);
  grad->SetMarkerSize(0.75);
  grad->GetYaxis()->CenterTitle();
  grad->GetYaxis()->SetTitle("B, Tesla");
  grad->GetYaxis()->SetTitleOffset(0.75);
  grad->GetXaxis()->SetTitle("U_{x}, mV");
  grad->Draw("AP");
  grad->Fit("pol1", "0");
  pol1 = grad->GetFunction("pol1");
  pol1->SetLineWidth(1);
  pol1->Draw("same");
  leg = new TLegend(0.20,0.65,0.50,0.75);
  //  (leg->AddEntry("NULL", "magnet  2SP- 40", "h"))->SetTextSize(0.040);
  leg->AddEntry(pol1, "B = 0.01048 U_{x}  -  0.00213");
  leg->SetFillColor(0);
  leg->SetTextAlign(22);
  leg->Draw();
  gPad->SetWindowSize(700,400);
  printTeX("magnet_hall_calib");
}
