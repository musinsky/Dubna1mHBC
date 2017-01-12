// @(#) 18 Jul 2006
// Author: Jan Musinsky

#include "TSystem.h"
#include "TObjString.h"
#include "Riostream.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"

TList *checkfiles(const char *dirname)
{
  // all data from http://nn-online.org
  
  // G. Bizard et al., Nucl. Phys. B85 (1975), 14-30 (SATURNE)
  // http://www.sciencedirect.com/science/article/B6TVC-472T3W3-1TJ/2/370f635a255050718b5f015a32cfa7bb
  
  // B.E. Bonner et al., Phys. Rev. Lett. 41 (1978), 1200-1203 (Los Alamos)
  // http://link.aps.org/abstract/PRL/v41/p1200
  
  // P.F. Shepard et al., Phys. Rev. D 10 (1974), 2735-2762 (Princeton)
  // http://link.aps.org/abstract/PRD/v10/p2735
  
  const char *file;
  void *dir = gSystem->OpenDirectory(dirname);
  TList *list = new TList();
  while ((file = gSystem->GetDirEntry(dir))) {
    TString tmp = Form("%s/%s", dirname, file);
    if (!tmp.EndsWith(".np")) continue;
    list->Add(new TObjString(tmp));
  }
  // sort the files in alphanumeric order (small -> high energy)
  list->Sort();
  return list;
}
TGraphErrors *creategraph(const char *file, const Bool_t noangle = kTRUE)
{
  ifstream fin(file);
  if (!fin) {
    Printf("wrong file %s", file);
    return 0;
  }
  TString tmp, tmpPrev;
  Double_t kenergy = -99, angle, sigma, sigmaErr;
  Int_t dsg = 0, num = 0;
  TGraphErrors *gr = 0;
  Double_t massp = 0.93827, massn = 0.93956, tenergy, invmass, pLab, pCMS;
  
  while (fin) {
    fin >> tmp;
    // kinetic energy
    if (tmp == "MeV," && tmpPrev.IsFloat()) {
      if (kenergy < 0) { // only once (first time)
	kenergy = tmpPrev.Atof();
	// transform kinetic energy (LAB) to impuls (CMS)
	tenergy = kenergy*0.001 + massn;
	invmass = massn*massn + massp*massp + 2.*tenergy*massp;
	pLab = TMath::Sqrt(tenergy*tenergy - massn*massn);
	pCMS = (massp*pLab)/TMath::Sqrt(invmass);
      }
      else if (kenergy != tmpPrev.Atof()) { // "verify" kenergy
	Printf("wrong kenergy %g", tmpPrev.Atof());
	return 0;
      }
    }
    // DSG (number of points)
    if (tmp == "DSG" && tmpPrev.IsDigit()) { // only once (first time)
      dsg = tmpPrev.Atoi();
      if (gr) { // "verify" DSG
	Printf("wrong DSG %d", dsg);
	return 0;
      }
      gr = new TGraphErrors(dsg);
    }
    // data
    if (tmp == "deg.," && tmpPrev.IsFloat()) {
      angle = tmpPrev.Atof();
      fin >> tmp >> tmp >> tmp; // DSG   =   x.xyz
      if (tmp.IsFloat())
	sigma = tmp.Atof();
      else {
	Printf("wrong sigma %s", tmp.Data());
	return 0;
      }
      fin >> tmp >> tmp; // "+-"   x.xyz
      if (tmp.IsFloat())
	sigmaErr = tmp.Atof();
      else {
	Printf("wrong sigmaErr %s", tmp.Data());
	return 0;
      }
      if (noangle) { // transform dSigma/dAngle => dSigma/dT
	angle = 180.0 - angle; // ?! charge exchange ?!
	angle *= TMath::DegToRad();
	angle = 2.*pCMS*pCMS*(1.0 - TMath::Cos(angle)); // => -t
	sigma *= TMath::Pi()/(pCMS*pCMS);
	sigmaErr *= TMath::Pi()/(pCMS*pCMS); // impulsErr = 0
      }
      gr->SetPoint(num, angle, sigma);
      gr->SetPointError(num, 0, sigmaErr);
      num++;
    }
    tmpPrev = tmp;
  }
  fin.close();
  if (dsg != num)
    Printf("DSG != num => %d != %d", dsg, num);
  gr->SetMarkerStyle(7);
  gr->SetTitle(Form("Ekin =%9.5f, Plab =%9.5f GeV/c", kenergy*0.001, pLab));
  return gr;
}
void fitgraph(const Double_t endf, TGraphErrors *gr, TGraphErrors *gr0)
{
  //  TF1 *fite = new TF1("fite","[0]*TMath::Exp([1]*x)+[2]*TMath::Exp([3]*x)");
  //  fite->SetParameter(0,10); // !!! must be set !!! ~ 10-100
  // tato funkcia lepsie opisuje funkciu na celom intervale
  // musi byt aspon prvy parameter nastaveny
  // ak fitujem len interval (0.00, 0.01) tak najprv fitovat cely interval
  // a az potm tento "zmenseny" intrval (t.j. pred-nastavene parametre)
  // znacne nestabilny fit
  
  TF1 *fite = new TF1("fite","[0]*TMath::Exp([1]*x+[2]*x*x)");
  // funguje velmi dobre na "zmensenom" intrvale (0.00, 0.01) a bez parametrov
  // na tomto intervale(rychlo rastucom) ma lepsi chi2/ndf
  // avsak ma problemy s celym intervalom
  
  gr->Fit(fite, "Q", "", 0.00, endf);
  TString tmp(gr->GetTitle());
  tmp.Remove(0,32); // after "Plab ="
  Double_t p0 = fite->GetParameter(0), p0Err = fite->GetParError(0);
  Int_t ip = gr0->GetN();
  // f(0) = p0 => f(0)Err = p0Err
  gr0->SetPoint(ip, tmp.Atof(), p0);
  gr0->SetPointError(ip, 0.00, p0Err);
  
  if (fite->GetNpar() == 4) {
    // f(0) = p0 + p2 => f(0)Err^2 = p0Err^2 + p2Err^2
    Double_t p2 = fite->GetParameter(2), p2Err = fite->GetParError(2);
    gr0->SetPoint(ip, tmp.Atof(), p0 + p2);
    gr0->SetPointError(ip, 0.00, TMath::Sqrt(p0Err*p0Err + p2Err*p2Err));
  }
}
void np_OLD3()
{
  TMultiGraph *bizards  = new TMultiGraph("bizards","Bizard");
  TMultiGraph *bonners  = new TMultiGraph("bonners","Bonner");
  TMultiGraph *shepards = new TMultiGraph("shepards","Shepard");
  TGraphErrors *bizard0  = new TGraphErrors();
  bizard0->SetName("bizard"); bizard0->SetLineColor(kRed);
  TGraphErrors *bonner0  = new TGraphErrors();
  bonner0->SetName("bonner"); bonner0->SetLineColor(kBlue);
  TGraphErrors *shepard0 = new TGraphErrors();
  shepard0->SetName("shepard"); shepard0->SetLineColor(kGreen);
  
  TList *files = checkfiles("NN-OnLine");
  TIter next(files);
  TObjString *obj;
  TGraphErrors *ge;
  while ((obj = (TObjString*)next())) {
    ge = creategraph(obj->GetName(), kTRUE);
    if (!ge) {
      Printf("wrong data in file %s", obj->GetName());
      continue;
    }
    if ((obj->GetString()).Contains("bizard")) {
      ge->SetTitle(Form("bizard,  %s", ge->GetTitle()));
      fitgraph(0.010, ge, bizard0);
      bizards->Add(ge);
    }
    else if ((obj->GetString()).Contains("bonner")) {
      ge->SetTitle(Form("bonner,  %s", ge->GetTitle()));
      fitgraph(0.015, ge, bonner0);
      bonners->Add(ge);
    }
    else if ((obj->GetString()).Contains("shepard")) {
      ge->SetTitle(Form("shepard, %s", ge->GetTitle()));
      fitgraph(0.030, ge, shepard0);
      shepards->Add(ge);
    }
  }
  files->Delete();
  delete files;
  
  //  shepards->Draw("AP"); bizards->Draw("P"); bonners->Draw("P");
  //TGraphErrors *ge = (TGraphErrors *)bizards->GetListOfGraphs()->At(14);

  bizard0->Draw("AP"); bonner0->Draw("P"); shepard0->Draw("P");
  Double_t x, y, e; bizard0->GetPoint(14, x, y); e = bizard0->GetErrorY(14);
  Printf("Plab = %g, dSigma/dT(0) = %g, error = %g", x, y, e);
}
