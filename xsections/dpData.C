// Author: Jan Musinsky
// 28/08/2014

#include <TSystem.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TPad.h>
#include <TStyle.h>
#include <fstream>

Double_t AsymmError(Double_t elow, Double_t ehigh)
{
  return TMath::Sqrt(0.5*(elow*elow + ehigh*ehigh));
}

TGraphErrors *ParsePDGData(const char *fname, const char *refer = "", Bool_t systematic = kFALSE)
{
  std::ifstream file(fname);
  if (!file.good()) {
    Printf("cannot read file %s", fname);
    return 0;
  }

  std::string line;
  Int_t count = 0;
  const char *format = "%*d %lg %lg %lg %lg %lg %lg %lg %lg %s %s"; // ignore first number
  TGraphErrors *gr = new TGraphErrors();

  while(std::getline(file, line, '\n')) {
    Double_t plab, plabMin, plabMax, sig, staErrP, staErrM, sysErrP, sysErrM;
    char refAuthor[64], refYear[8];
    Int_t nwords = sscanf(line.c_str(), format, &plab, &plabMin, &plabMax, &sig,
                          &staErrP, &staErrM, &sysErrP, &sysErrM, refAuthor, refYear);
    TString ref = TString::Format("%s %s", refAuthor, refYear);

    //    Printf("%s", line.c_str());
    //    Printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t'%s'+'%s'\t=> %d", plab, plabMin, plabMax, sig,
    //           staErrP, staErrM, sysErrP, sysErrM, refAuthor, refYear, nwords);

    line.clear();
    if (nwords != 10) continue;         // only real data
    if (!ref.Contains(refer)) continue; // only refer data

    Double_t staErr = AsymmError(staErrP, staErrM);
    Double_t sysErr = AsymmError((sysErrP/100.0)*sig, (sysErrM/100.0)*sig); // in percent
    if (!systematic) sysErr = 0.0;
    Double_t sigErr = TMath::Sqrt(staErr*staErr + sysErr*sysErr);
    gr->SetPoint(count, plab, sig);
    gr->SetPointError(count, (plabMax - plabMin)/2.0, sigErr);
    count++;
  }

  file.close();
  Printf("found %d data points", gr->GetN());
  if (gr->GetN() == 0) {
    delete gr;
    return 0;
  }
  return gr;
}

void dpData()
{
  TString fname = "rpp2014-pdeut_total.dat";
  //  fname = "rpp2014-np_total.dat";
  TString furl = "http://pdg.lbl.gov/2014/hadronic-xsections/" + fname;
  if (gSystem->AccessPathName(fname))
    gSystem->Exec(TString::Format("wget %s", furl.Data()));

  TGraphErrors *gr = ParsePDGData(fname, "", kFALSE);
  if (!gr) return;
  gr->SetTitle("Total cross sections for p-d collisions");
  gr->GetXaxis()->SetTitle("p_{LAB} (GeV/c)");
  gr->GetYaxis()->SetTitle("#sigma (mb)");
  gr->GetYaxis()->CenterTitle();
  gr->SetMarkerStyle(20);
  gr->Draw("AP");

  gStyle->SetGridColor(kGray);
  gPad->SetLogx();
  gPad->SetTicks();
  gPad->SetGrid();
  gr->GetYaxis()->SetRangeUser(44., 96.);

  gr = ParsePDGData(fname, "BUGG", kFALSE);
  if (!gr) return;
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kRed);
  gr->SetLineColor(gr->GetMarkerColor());
  gr->Draw("P");
}
