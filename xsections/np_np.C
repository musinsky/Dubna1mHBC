// Author: Jan Musinsky
// 21/09/2014

#include <TROOT.h>
#include <TSystem.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TVirtualFitter.h>
#include <fstream>

void CheckDirFiles(const char *dir, TList *files)
{
  void *pdir = gSystem->OpenDirectory(dir);
  if (!pdir || !files) return;

  const char *file;
  while ((file = gSystem->GetDirEntry(pdir))) {
    TString dirFile = TString::Format("%s/%s", dir, file);
    if (!dirFile.EndsWith(".np")) continue;   // only data file
    if (files->FindObject(dirFile)) continue; // no duplicate
    files->Add(new TObjString(dirFile));
  }
  files->Sort();
}

void ParseNNOnlineData(const char *fname)
{
  FILE *file = fopen(fname, "r");
  if (!file) return;
  char *line = 0; // must be NULL, otherwise backtrace after call: free(line);
  size_t len;
  const char *format = "%lg %*s %lg %*s %*s %*s %lg %*s %lg";

  while (getline(&line, &len, file) != -1) {
    Double_t tLab, cmAngle, dsg, dsgErr;
    Int_t nwords = sscanf(line, format, &tLab, &cmAngle, &dsg, &dsgErr);

    //    printf("%s", line);
    //    Printf("%f\t%f\t%f\t%f\t=> %d", tLab, cmAngle, dsg, dsgErr, nwords);

    if (nwords != 4) continue; // only real data
    // TODO
  }

  free(line);
  fclose(file);
}

TGraphErrors *CreateGraph(const char *file)
{
  ifstream fin(file);
  if (!fin) {
    Printf("file %s does not exist", file);
    return 0;
  }
  TString word, word_prev;
  Double_t energy_full, energy_kin = -9999, mass_inv, p_lab, p_cms, angle;
  Double_t mass_p = 0.93827, mass_n = 0.93956, sigma, sigma_err;
  Int_t dsg = 0, num = 0;
  TGraphErrors *gr = 0;

  while (fin) {
    fin >> word;
    // kinetic energy
    if (word == "MeV," && word_prev.IsFloat()) {
      if (energy_kin < 0) { // only once (first time)
        energy_kin = word_prev.Atof();
        // transform kinetic energy (LAB) to impuls (CMS) for reaction np->np
        energy_full = energy_kin*0.001 + mass_n; // MeV -> GeV
        mass_inv = mass_n*mass_n + mass_p*mass_p + 2.0*energy_full*mass_p;
        p_lab = TMath::Sqrt(energy_full*energy_full - mass_n*mass_n);
        p_cms = (mass_p*p_lab)/TMath::Sqrt(mass_inv);
      }
      else if (energy_kin != word_prev.Atof()) { // "verify" kinetic energy
        Printf("wrong kinetic energy %g", word_prev.Atof());
        return 0;
      }
    }
    // DSG (number of points)
    if (word == "DSG" && word_prev.IsDigit()) { // only once (first time)
      dsg = word_prev.Atoi();
      if (gr) { // "verify" DSG
        Printf("wrong DSG %d", dsg);
        return 0;
      }
      gr = new TGraphErrors(dsg);
    }
    // data
    if (word == "deg.," && word_prev.IsFloat()) {
      angle = word_prev.Atof();
      fin >> word >> word >> word; // DSG   =   sigma(5.43210)
      if (word.IsFloat())
        sigma = word.Atof();
      else {
        Printf("wrong sigma %s", word.Data());
        return 0;
      }
      fin >> word >> word; //         "+-"   sigma_err(0.12345)
      if (word.IsFloat())
        sigma_err = word.Atof();
      else {
        Printf("wrong sigma_err %s", word.Data());
        return 0;
      }
      // transform dSigma/dAngle => dSigma/dT
      angle = TMath::DegToRad()*(180.0 - angle); // ?! charge exchange ?!
      angle = 2.0*p_cms*p_cms*(1.0 - TMath::Cos(angle)); // => -t
      sigma *= TMath::Pi()/(p_cms*p_cms);
      sigma_err *= TMath::Pi()/(p_cms*p_cms); // impuls_err = 0
      gr->SetPoint(num, angle, sigma);
      gr->SetPointError(num, 0, sigma_err);
      num++;
    }
    word_prev = word;
  }

  fin.close();
  if (dsg != num) Printf("DSG != num => %d != %d", dsg, num);
  gr->SetMarkerStyle(7);
  gr->SetTitle(Form("E_kin =%9.5f, P_lab =%9.5f GeV/c",
                    energy_kin*0.001, p_lab)); // MeV -> GeV
  gr->SetUniqueID(UInt_t(p_lab*100000)); // store in graph p_lab value
  return gr;
}

void FitGraph(const Double_t t_max, TGraphErrors *gr, TGraphErrors *gr0)
{
  TF1 *fite = new TF1("fite","[0]*TMath::Exp([1]*x+[2]*x*x)");
  // funguje velmi dobre na "zmensenom" intrvale (0.00, 0.01) a bez parametrov
  // avsak ma problemy s celym intervalom

  fite->SetLineWidth(1);
  gr->Fit(fite, "Q", "", 0.00, t_max);
  Int_t ip = gr0->GetN();
  // f(0) = p0 => f(0)_err = p0_err
  gr0->SetPoint(ip, gr->GetUniqueID()/100000.0, fite->GetParameter(0));
  gr0->SetPointError(ip, 0.00, fite->GetParError(0));
}

void np_np()
{
  TMultiGraph *bizards  = new TMultiGraph("bizards","Bizard");
  TMultiGraph *bonners  = new TMultiGraph("bonners","Bonner");
  TMultiGraph *shepards = new TMultiGraph("shepards","Shepard");
  TGraphErrors *bizard0  = new TGraphErrors();
  bizard0->SetName("bizard0"); bizard0->SetMarkerStyle(8);
  TGraphErrors *bonner0  = new TGraphErrors();
  bonner0->SetName("bonner0"); bonner0->SetMarkerStyle(24);
  TGraphErrors *shepard0 = new TGraphErrors();
  shepard0->SetName("shepard0"); shepard0->SetMarkerStyle(26);

  TList *files = new TList();
  CheckDirFiles("NN-OnLine", files);
  CheckDirFiles("NN-OnLine/others", files);
  TIter next(files);
  TObjString *obj;
  TGraphErrors *ge;
  while ((obj = (TObjString*)next())) {
    ge = CreateGraph(obj->GetName());
    if (!ge) {
      Printf("wrong data in file %s", obj->GetName());
      continue;
    }
    if ((obj->GetString()).Contains("bizard")) {
      ge->SetTitle(Form("bizard,  %s", ge->GetTitle()));
      FitGraph(0.010, ge, bizard0);
      bizards->Add(ge);
    }
    else if ((obj->GetString()).Contains("bonner")) {
      ge->SetTitle(Form("bonner,  %s", ge->GetTitle()));
      FitGraph(0.015, ge, bonner0);
      bonners->Add(ge);
    }
    else if ((obj->GetString()).Contains("shepard")) {
      ge->SetTitle(Form("shepard, %s", ge->GetTitle()));
      FitGraph(0.030, ge, shepard0);
      shepards->Add(ge);
    }
  }
  files->Delete(); // delete TObjString
  delete files;

  gStyle->SetOptStat(0);
  new TCanvas();
  TLegend *leg = new TLegend(0.60,0.75,0.85,0.85);
  // leg->SetHeader("Bizard (Saturne)");
  // ((TLegendEntry *)leg->GetListOfPrimitives()->Last())->SetTextAlign(22);
  ge = (TGraphErrors *)bizards->GetListOfGraphs()->At(9);
  ge->GetXaxis()->SetRangeUser(0.0, 0.01); ge->SetTitle(0);
  ge->SetMarkerStyle(4); ge->SetMarkerSize(0.70);
  ge->SetMinimum(15.0); ge->SetMaximum(85.0);
  ge->GetXaxis()->SetTitle("#left|t#right|, (GeV/c)^{2}");
  ge->GetYaxis()->SetTitle("#frac{d#sigma}{dt}, #frac{mb}{(GeV/c)^{2}}");
  ge->GetYaxis()->CenterTitle(); ge->GetYaxis()->SetTitleOffset(1.10);
  ge->GetYaxis()->SetNdivisions(905);
  ge->Draw("AP");
  leg->AddEntry(ge,Form("p = %4.2f GeV/c", ge->GetUniqueID()/100000.0),"P");
  ge = (TGraphErrors *)bizards->GetListOfGraphs()->At(18);
  ge->GetXaxis()->SetRangeUser(0.0, 0.01); ge->SetTitle(0);
  ge->SetMarkerStyle(8); ge->SetMarkerSize(0.70);
  ge->Draw("P");
  leg->AddEntry(ge,Form("p = %4.2f GeV/c", ge->GetUniqueID()/100000.0),"P");
  leg->SetFillColor(0);
  leg->Draw();
  //  printTeX("np_two_bizard");

  new TCanvas();
  bizard0->GetXaxis()->SetRangeUser(1.0, 2.0);
  bizard0->SetMinimum(0.0); bizard0->SetMaximum(210.0);
  bizard0->GetXaxis()->SetTitle("p, GeV/c");
  bizard0->GetYaxis()->SetTitle("#frac{d#sigma}{dt} #void8_{ t = 0}, #frac{mb}{(GeV/c)^{2}}");
  bizard0->GetYaxis()->CenterTitle(); bizard0->GetYaxis()->SetTitleOffset(1.10);
  bizard0->GetYaxis()->SetNdivisions(905);
  bizard0->Draw("AP");
  bonner0->Draw("P");
  shepard0->Draw("P");
  leg = new TLegend(0.55,0.685,0.85,0.85);
  leg->AddEntry(bizard0,"Saclay (Bizard)","P");
  leg->AddEntry(bonner0,"Los Alamos (Bonner)","P");
  leg->AddEntry(shepard0,"Princeton (Shepard)","P");
  leg->SetFillColor(0);
  leg->Draw();
  Double_t impuls = 3.35/2.;
  TF1 *fe = new TF1("fe","[0]*TMath::Exp([1]*x)");
  fe->SetParameter(0, 100); fe->SetLineWidth(1);
  bizard0->Fit(fe, "Q");
  Double_t p0 = fe->GetParameter(0), p0_err = fe->GetParError(0);
  Double_t p1_err = fe->GetParError(1);
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();
  Double_t cov = fitter->GetCovarianceMatrixElement(0, 1); // = (1, 0)
  Double_t error = (p0_err*p0_err)/(p0*p0) + (impuls*impuls*p1_err*p1_err);
  error += 2.0*cov*(impuls/p0);
  error = fe->Eval(impuls)*TMath::Sqrt(error);
  Printf("p = %5.3f, sigma(0) = %g +- = %g", impuls, fe->Eval(impuls), error);
  //  printTeX("np_sigma0");

  TObjArray *t0_graphs = new TObjArray();
  t0_graphs->AddVector(bizard0, bonner0, shepard0, 0);
  Double_t t0, t0_err, p, p_min = 1.0;
  const char *latex = "%.3f & %.3f $\\pm$ %.3f & %s \\\\ \\hline";
  TList *entries = new TList();

  TIter nextg(t0_graphs);
  while ((ge = (TGraphErrors *)nextg())) {
    for (Int_t i = 0; i < ge->GetN(); i++) {
      ge->GetPoint(i, p, t0);
      t0_err = ge->GetErrorY(i);
      if (p < p_min) continue;
      entries->Add(new TObjString(Form(latex, p, t0, t0_err, ge->GetName())));
    }
  }

  entries->Sort(); // small -> high energy
  TIter nexte(entries);
  while ((obj = (TObjString *)nexte()))
    Printf("%s", obj->GetName());
}
