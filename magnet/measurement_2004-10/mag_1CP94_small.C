// 2004-10-29 original
// 2019-05-15 update to ROOT6
// Musinsky Jan

void mag_1CP94_small()
{
  // 2004-10 measurement (?! with not correct Hall sensor ?!)
  // other Hall sensor than 2005-01 measurement

  const UShort_t n_point=33;
  Float_t B[n_point];
  for (UShort_t i=0; i < n_point; i++)
    B[i] = 0.15 + i*0.05;

  Float_t Ux_minus[n_point] =
    { 14.24, 18.98, 23.71, 28.44, 33.17, 37.89, 42.61, 47.32, 52.02, 56.72,
      61.42, 66.10, 70.78, 75.46, 80.12, 84.78, 89.46, 94.11, 98.77, 103.42,
      108.07, 112.70, 117.34, 121.97, 126.60, 131.22, 135.85, 140.47, 145.08,
      149.67, 154.28, 158.90, 163.51 };

  TGraph *gra_1CP94 = new TGraph(n_point,Ux_minus,B);
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetGridx();
  c1->SetGridy();

  gra_1CP94->SetMarkerStyle(7);
  gra_1CP94->SetTitle("1CP-94 small magnet (old 2004-10)");
  gra_1CP94->GetYaxis()->SetTitle("B [ Tesla ]");
  gra_1CP94->GetYaxis()->CenterTitle(kTRUE);
  gra_1CP94->GetYaxis()->SetTitleOffset(1.15);
  gra_1CP94->GetXaxis()->SetTitle("Ux- [ mV ]");
  gra_1CP94->GetXaxis()->CenterTitle(kTRUE);
  gra_1CP94->Draw("AP");

  gra_1CP94->Fit("pol1");
  TF1 *pol1 = gra_1CP94->GetFunction("pol1");
  pol1->SetLineColor(kBlack);
  pol1->SetLineWidth(1);

  Double_t q = pol1->GetParameter(0);
  Double_t k = pol1->GetParameter(1);

  TLatex *bla = new TLatex(20, 1.45,
                           TString::Format(" B = %8.6f * Ux   %8.6f", k, q));
  bla->SetTextSize(0.04);
  bla->Draw();

  c1->Update();
  c1->Modified();
}
