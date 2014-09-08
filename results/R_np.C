// Author: Jan Musinsky
// 23/10/2009

{
  gROOT->SetStyle("Plain");
  const Int_t n = 31;
  Double_t x[n] = {
    13.5, 13.9, 30.1, 50.0, 90.0, 95.0, 95.7,135.0,143.9,152.0,200.0,
    270.0,299.7,319.8,339.7,359.6,379.6,380.0,399.7,419.8,440.0,460.1,
    480.4,500.9,521.1,539.4,557.4,647.0,710.0,794.0,800.0};
  for (Int_t i = 0; i < n; i++) x[i]/=1000.; // MeV to GeV

  Double_t sd[n] = {
    0.180,0.185,0.141,0.240,0.397,0.480,0.587,0.652,0.601,0.650,0.553,
    0.710,0.652,0.643,0.637,0.626,0.641,0.200,0.610,0.623,0.630,0.611,
    0.608,0.592,0.604,0.617,0.632,0.600,0.483,0.560,0.660};
  Double_t sde[n] = {
    0.000,0.000,0.035,0.060,0.044,0.030,0.029,0.154,0.057,0.100,0.030,
    0.021,0.033,0.032,0.032,0.031,0.032,0.035,0.031,0.031,0.032,0.031,
    0.030,0.030,0.030,0.031,0.032,0.080,0.080,0.040,0.080};
  Double_t id[n] = {
    2.704,2.604,3.278,1.778,0.679,0.389,0.137,0.023,0.109,0.026,0.206,
    -0.061,0.023,0.037,0.047,0.065,0.040,2.333,0.093,0.070,0.058,0.091,
    0.097,0.126,0.104,0.081,0.055,0.111,0.150,0.191,0.010};
  Double_t ide[n] = {
    0.000,0.000,1.174,0.694,0.186,0.087,0.056,0.241,0.105,0.158,0.065,
    0.029,0.051,0.052,0.053,0.053,0.053,0.583,0.055,0.033,0.053,0.055,
    0.054,0.057,0.055,0.054,0.053,0.148,0.191,0.085,0.122};
  TString r[n] = {
    "b",  "a",  "b",  "b",  "a",  "b",  "b",  "b",  "b",  "a",  "a",
    "a",  "a",  "a",  "a",  "a",  "a",  "a",  "a",  "a",  "a",  "a",
    "a",  "a",  "a",  "a",  "a",  "b",  "a",  "a",  "b"};
  TString lab[n] = {
    "Livermore", "Moscow", "Rutherford", "Rutherford", "UCRL", "Harvard U.",
    "Harwell", "Harwell", "Harwell", "Harvard U.", "JINR DNLP", "UCRL", "PSI",
    "PSI", "PSI", "PSI", "PSI", "INP Dubna", "PSI", "PSI", "PSI", "PSI",
    "PSI", "PSI", "PSI", "PSI", "PSI", "LAMPF", "LRL", "LAMPF", "LAMPF"};

  TGraphErrors *gr = new TGraphErrors(n, x, sd, 0, sde);
  gr->RemovePoint(17);
  gr->GetXaxis()->SetLimits(0.0, 2.1);
  gr->SetTitle(0); gr->SetMinimum(0.0); gr->SetMaximum(1.0);
  gr->GetXaxis()->SetTitle("T_{n }, GeV");
  gr->GetXaxis()->SetNdivisions(409);
  gr->GetYaxis()->SetTitle("R_{np #rightarrow pn}");
  gr->GetYaxis()->CenterTitle(); gr->GetYaxis()->SetTitleOffset(1.25);
  gr->GetYaxis()->SetNdivisions(305);
  gr->SetMarkerStyle(20); gr->SetMarkerSize(1.00);
  gr->SetMarkerColor(kGray+2);
  gr->Draw("AP");

  // Delta-Sigma
  const Int_t n2 = 7;
  Double_t x2[n2] = {0.55, 0.80, 1.00, 1.20, 1.40, 1.80, 2.00};
  Double_t y2[n2] = {0.589, 0.554, 0.553, 0.551, 0.576, 0.568, 0.564};
  Double_t y2e[n2] = {0.046, 0.023, 0.026, 0.022, 0.038, 0.033, 0.045};
  TGraphErrors *gr2 = new TGraphErrors(n2, x2, y2, 0, y2e);
  gr2->SetMarkerStyle(24); gr2->SetMarkerSize(1.00);
  gr2->SetMarkerColor(kGray+2);
  gr2->Draw("P");

  // Strela
  TGraphErrors *gr3 = new TGraphErrors();
  gr3->SetPoint(0, 0.981, 0.55);
  gr3->SetPointError(0, 0, 0.08);
  gr3->SetMarkerStyle(21); gr3->SetMarkerSize(1.15);
  gr3->Draw("P");

  //  gPad->SetGridx();
  //  gPad->SetGridy();
  gPad->SetWindowSize(500,500);
  gROOT->LoadMacro("$ROOTSYS/macros/printTeX.C");
  printTeX("R_np");
}
