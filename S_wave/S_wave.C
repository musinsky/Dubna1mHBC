// Author: Jan Musinsky
// 11/09/2014

#include <TGraph.h>
#include <TH1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>

void S_wave()
{
  const char *fmt = "%lg %*lg %*lg %*lg %lg";
  TGraph *g1 = new TGraph("wave2d", fmt);
  g1->SetTitle("Paris");
  TGraph *g2 = new TGraph("ba2d", fmt);
  g2->SetTitle("Bonn A");
  TGraph *g3 = new TGraph("bb2d", fmt);
  g3->SetTitle("Bonn B");
  TGraph *g4 = new TGraph("bc2d", fmt);
  g4->SetTitle("Bonn C");

  g1->Draw("AL");
  g1->GetHistogram()->SetBit(TH1::kNoTitle);
  g1->SetMinimum(-0.19);
  g1->SetMaximum(1.19);
  g1->GetXaxis()->SetRangeUser(0.0, 0.5);
  g1->GetXaxis()->SetTitle("p (GeV/c)");
  g1->GetYaxis()->SetTitle("u^{2} / (u^{2} + w^{2})");
  // yaxis in old articles = nucleon polarization
  g1->GetYaxis()->SetTitleOffset(1.25);
  g1->GetYaxis()->CenterTitle();
  g2->SetLineStyle(2);
  g2->Draw("L");
  g3->SetLineStyle(3);
  g3->Draw("L");
  g4->SetLineStyle(4);
  g4->Draw("L");

  TLegend *leg = new TLegend(0.55, 0.70, 0.85, 0.85);
  leg->AddEntry(g1, "", "L");
  leg->AddEntry(g2, "", "L");
  leg->AddEntry(g3, "", "L");
  leg->AddEntry(g4, "", "L");
  leg->SetTextAlign(22);
  leg->SetLineColor(kGray);
  leg->Draw();

  gStyle->SetGridColor(kGray);
  gPad->SetGrid();
  gPad->GetCanvas()->SetWindowSize(500, 500);
}
