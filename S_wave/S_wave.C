// Author: Jan Musinsky
// 01/10/2009

{
  gROOT->SetStyle("Plain");
  g1 = new TGraph("wave2d.csv","%lg, %lg");
  g1->SetTitle("Paris");
  g2 = new TGraph("ba2d.csv","%lg, %lg");
  g2->SetTitle("Bonn A");
  g3 = new TGraph("bb2d.csv","%lg, %lg");
  g3->SetTitle("Bonn B");
  g4 = new TGraph("bc2d.csv","%lg, %lg");
  g4->SetTitle("Bonn C");

  g1->Draw("AL");
  g1->SetTitle(0);
  g1->SetMinimum(-0.19);
  g1->SetMaximum(1.19);
  g1->GetXaxis()->SetRangeUser(0, 0.5);
  g1->GetXaxis()->SetTitle("p, GeV/c");
  g1->GetYaxis()->SetTitle("u^{2} / (u^{2} + w^{2})");
  g1->GetYaxis()->SetTitleOffset(1.2);
  g1->GetYaxis()->CenterTitle();
  g2->Draw("L"); g2->SetLineStyle(2);
  g3->Draw("L"); g3->SetLineStyle(3);
  g4->Draw("L"); g4->SetLineStyle(4);

  leg = new TLegend(0.55,0.7,0.85,0.85);
  leg->AddEntry(g1,"Paris","l");
  leg->AddEntry(g2,g2->GetTitle(),"l");
  leg->AddEntry(g3,g3->GetTitle(),"l");
  leg->AddEntry(g4,g4->GetTitle(),"l");
  leg->SetFillColor(0);
  leg->SetTextAlign(22);
  leg->Draw();

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetWindowSize(500,500);
  gROOT->LoadMacro("$ROOTSYS/macros/printTeX.C");
  printTeX("S_waves");
}
