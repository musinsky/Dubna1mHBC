// 2004-10-29 original
// 2019-05-15 update to ROOT6
// Musinsky Jan

void mag_2CP40_big()
{
  // 2004-10 measurement (?! with not correct Hall sensor ?!)
  // other Hall sensor than 2005-01 measurement

  const UShort_t n_point=37;
  Float_t B[n_point];
  for (UShort_t i=0; i < n_point; i++)
    B[i] = 0.15 + i*0.05;

  Float_t Ux_minus[n_point] =
    { 14.19, 18.94, 23.68, 28.43, 33.17, 37.89, 42.63, 47.36, 52.07, 56.79,
      61.50, 66.21, 70.90, 75.59, 80.28, 84.97, 89.64, 94.32, 99.00, 103.66,
      108.32, 112.98, 117.63, 122.29, 126.92, 131.58, 136.21, 140.84, 145.48,
      150.11, 154.74, 159.35, 163.98, 168.60, 173.22, 177.83, 182.44 };

  TGraph *gra_2CP40 = new TGraph(n_point,Ux_minus,B);
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetGridx();
  c1->SetGridy();

  gra_2CP40->SetMarkerStyle(7);
  gra_2CP40->SetTitle("2CP-40 big magnet (old 2004-10)");
  gra_2CP40->GetYaxis()->SetTitle("B [ Tesla ]");
  gra_2CP40->GetYaxis()->CenterTitle(kTRUE);
  gra_2CP40->GetYaxis()->SetTitleOffset(1.15);
  gra_2CP40->GetXaxis()->SetTitle("Ux- [ mV ]");
  gra_2CP40->GetXaxis()->CenterTitle(kTRUE);
  gra_2CP40->Draw("AP");

  gra_2CP40->Fit("pol1");
  TF1 *pol1 = gra_2CP40->GetFunction("pol1");
  pol1->SetLineColor(kBlack);
  pol1->SetLineWidth(1);

  Double_t q = pol1->GetParameter(0);
  Double_t k = pol1->GetParameter(1);

  TLatex *bla = new TLatex(20, 1.65,
                           TString::Format(" B = %8.6f * Ux   %8.6f", k, q));
  bla->SetTextSize(0.04);
  bla->Draw();

  c1->Update();
  c1->Modified();
}
