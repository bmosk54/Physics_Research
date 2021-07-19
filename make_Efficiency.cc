#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TF1.h"
#include <TGraph.h>
#include <TGraphErrors.h>

#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TMinuit.h>
#include <TLatex.h>

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
void SetHistStyle(TH1D* h){
}

TFile *f;


TGraphErrors*  H3Acceptance;
TGraphErrors*  GrThresholdError;
TGraph*  GrThreshold;


std::vector<double> diffRecoParentEnergyParent;
std::vector<double> diffRecoParentEnergyNorm;
std::vector<double> diffRecoParentEnergyMean;
std::vector<double> diffRecoParentEnergyStdDev;

int NEBINS=20;
int NS1BINS=20;
int NRecoEnergyBins=20;

TH1D* tspectrum;
TH1D* tspectrum_all;
TH1D* tspectrum_ss;

//TH1D* energyParentAll;
TH1D* energyParent_All;
TH1D* energyParent_ss;

TH1D* energyRecBeforeFit_keV;
//TH1D* energyParentBeforeFit_keV;

std::vector<TPad*>  make_pads(TCanvas* c, int nx, int ny){

  c->cd();

  std::vector<TPad*>  PadList;

  double dx=1.0/float(nx);
  double dy=1.0/float(ny);

  for(int ix=0;ix<nx;ix++){
    for(int iy=0;iy<ny;iy++){
      double x1=ix*dx;
      double x2=(ix+1)*dx;

      double y1=iy*dy;
      double y2=(iy+1)*dy;

      TPad* pad = new TPad("c","c",x1,y1,x2,y2);
      PadList.push_back(pad);
    }
  }
  return PadList;
}

std::vector<TH1D*> GetTH1DVector(int NBINS,string hname1){

  std::vector<TH1D*> TH1DVector;
  string hname;
  for(int i=0;i<NBINS;i++){
    hname =  hname1+to_string(i);
    TH1DVector.push_back((TH1D*) f->Get(hname.c_str()));
  }
  return TH1DVector;
}

std::vector<TProfile*> GetTProfileVector(int NBINS,string hname1){

  std::vector<TProfile*> TProfileVector;

  string hname;
  for(int i=0;i<NBINS;i++){
    hname =  hname1+to_string(i);
    TProfileVector.push_back((TProfile*) f->Get(hname.c_str()));
  }
  return TProfileVector;
}



void effX() {
   auto effErf = [](double* x, double* p) {
      return (TMath::Erf((x[0] - p[0]) / p[1]) + 1) / 2. * p[2];
   };

   TF1* myErf = new TF1("myErf", effErf, 0., 10.);
   myErf->SetParameter(0, 5.);
   myErf->SetParameter(1, 5.);
   myErf->SetParameter(2, 1.);

   //   eff->Fit(myErf);
}
double FermiFunction(double T){
  double masse = 511;  // mass of electron
  double Zp = 2.0;     // charge of daughter element i.e. He
  double alpha=1.0/137.; // The alpha_em

  double beta=sqrt(2*T/masse);  // velocity

  double xtmp = 2.0*M_PI*alpha*Zp/beta;
  double denom = (1.0-exp(-xtmp));
  double fermi = xtmp/denom;
  return fermi;
}
/*
Minimizer is Linear
Chi2                      =       249.85
NDf                       =          179
p0                        =      3.19291   +/-   0.00741538  
p1                        =    0.0588292   +/-   0.00138176  
p2                        =  -0.00219583   +/-   8.59119e-05 
p3                        =  4.72473e-05   +/-   2.43349e-06 
p4                        = -5.74813e-07   +/-   3.5487e-08  
p5                        =  3.90672e-09   +/-   2.75161e-10 
p6                        = -1.37867e-11   +/-   1.07518e-12 
p7                        =  1.96281e-14   +/-   1.66177e-15 
*/

double GetMedianWimp_old(double x) {
  double p0  =      3.19291;
  double p1  =    0.0588292;
  double p2  =  -0.00219583;
  double p3  =  4.72473e-05;
  double p4  = -5.74813e-07;
  double p5  =  3.90672e-09;
  double p6  = -1.37867e-11;
  double p7  =  1.96281e-14;

  double y=p0+p1*pow(x,1)+p2*pow(x,2)+p3*pow(x,3)+p4*pow(x,4)+p5*pow(x,5)+p6*pow(x,6)+p7*pow(x,7);
  return y;
}

/*
inimizer is Linear
Chi2                      =      1471.05
NDf                       =          351
p0                        =      3.12198   +/-   0.0109153   
p1                        =    0.0818549   +/-   0.00290273  
p2                        =  -0.00454078   +/-   0.000269032 
p3                        =   0.00015603   +/-   1.195e-05   
p4                        = -3.26929e-06   +/-   2.9261e-07  
p5                        =  4.23779e-08   +/-   4.2257e-09  
p6                        = -3.39294e-10   +/-   3.6816e-11  
p7                        =  1.62543e-12   +/-   1.89769e-13 
p8                        = -4.25605e-15   +/-   5.31928e-16 
p9                        =  4.67382e-18   +/-   6.23991e-19 
*/
double GetMedianWimp(double x) {
  double p0                        =      3.12198;
  double p1                        =    0.0818549;
  double p2                        =  -0.00454078;
  double p3                        =   0.00015603;
  double p4                        = -3.26929e-06;
  double p5                        =  4.23779e-08;
  double p6                        = -3.39294e-10;
  double p7                        =  1.62543e-12;
  double p8                        = -4.25605e-15;
  double p9                        =  4.67382e-18;

  double y=p0+p1*pow(x,1)+p2*pow(x,2)+p3*pow(x,3)+p4*pow(x,4)+p5*pow(x,5)+p6*pow(x,6)+p7*pow(x,7)+p8*pow(x,8)+p9*pow(x,9);
  return y;
}

void make_Efficiency(){

  string sampleName="WIMP";
  string pngFile;

  string hname;
  int ix=0;int iy=0; int ipad =0;
  int nx=1;int ny=1;


  TText *t;

  gStyle->SetOptStat(kTRUE);
  gStyle->SetOptStat(11111111);
  gStyle->SetOptTitle(kTRUE);

  f = new TFile("H3_plots.root");f->cd();f->ls();


  TH2D* S1Area_phd_logS2Area_phd_ss = (TH2D*) f->Get("S1Area_phd_logS2Area_phd_ss");
       
  int nxbins =  S1Area_phd_logS2Area_phd_ss->GetNbinsX();
  int nybins =  S1Area_phd_logS2Area_phd_ss->GetNbinsY();

  int sumYAll[nxbins];
  int sumYBelow[nxbins];
  int sumYAbove[nxbins];


  std:: cout <<  " opened file " << nxbins << "  " << nybins << endl;

  nxbins=200;

  double S1bins[nxbins];
  double ex[nxbins];
  double Acceptance[nxbins];
  double AccError[nxbins];

  double Theshold[nxbins];

  for(int ix=1;ix<nxbins;ix++){
    int ix0=ix-1;
    sumYAll[ix0]=0.;
    sumYBelow[ix0]=0.;
    sumYAbove[ix0]=0.;
  
    double x =  S1Area_phd_logS2Area_phd_ss->GetXaxis()->GetBinCenter(ix); //->GetBinCenterX(ix);
    double MeanY = 1.0*GetMedianWimp(x);
  
    S1bins[ix0]=x;
    ex[ix]=0.0;
    Theshold[ix0]=MeanY;

    cout << " ix " << ix << " x " <<  x << " meanY " << MeanY  << endl;

    for(int iy=1;iy<nybins;iy++){
      double y =  double(S1Area_phd_logS2Area_phd_ss->GetYaxis()->GetBinCenter(iy));
      double nevents =  S1Area_phd_logS2Area_phd_ss->GetBinContent(ix,iy);
      //      cout << " ix " << ix  << " iy " << iy << " y " << y << " events " << nevents << endl;

      sumYAll[ix0]+=nevents;
      if(y<=MeanY) {
	sumYBelow[ix0]+=nevents;
      }
      else {
	sumYAbove[ix0]+=nevents;
      }
    }
  }

  for(int ix=0;ix<nxbins-1;ix++){
    cout << " ix " << ix;
    cout << " all " <<  sumYAll[ix];
    cout << " below " <<  sumYBelow[ix];
    cout << " above " <<  sumYAbove[ix];
    if(sumYAll[ix]>0){
      double ratio = double(sumYBelow[ix])/double(sumYAll[ix]);
      double error = sqrt(ratio*(1-ratio)/sumYAll[ix]);
	cout << " ratio " <<  ratio <<  " +/- " << error;
	Acceptance[ix]=ratio;
        AccError[ix]=error;

    }
    cout << endl;
  }

  TCanvas* C1 = new TCanvas("C1","C1",1000,800);
  nx=1;ny=2; std::vector<TPad*> C1Pads =  make_pads(C1,nx,ny);

  gStyle->SetOptStat(kFALSE);

  H3Acceptance = new TGraphErrors(nxbins-1,S1bins,Acceptance,ex,AccError);
  H3Acceptance = new TGraphErrors(nxbins-1,S1bins,Acceptance,ex,AccError);
  GrThreshold  = new TGraph(nxbins-1,S1bins,Theshold);

  H3Acceptance->SetTitle("Fraction of H3 Events below the Mean of WIMP Response");
  H3Acceptance->SetMarkerColor(2);
  H3Acceptance->SetMarkerStyle(20);
  H3Acceptance->SetMarkerSize(1);
  H3Acceptance->SetMaximum(0.2);
  H3Acceptance->SetMinimum(0.0);

  GrThreshold->SetLineColor(2);
  GrThreshold->SetLineWidth(2);


  S1Area_phd_logS2Area_phd_ss->GetXaxis()->SetRangeUser(0,100);
  S1Area_phd_logS2Area_phd_ss->GetYaxis()->SetRangeUser(2.0,6.0);

  ix=0;iy=1; ipad = ny*ix+iy; C1->cd();C1Pads[ipad]->Draw(); C1Pads[ipad]->cd(); S1Area_phd_logS2Area_phd_ss->Draw("colz"); GrThreshold->Draw("L");

  
  TH2D* S1AreaTemp = (TH2D*) S1Area_phd_logS2Area_phd_ss->Clone();
  S1AreaTemp->Reset();
  S1AreaTemp->GetXaxis()->SetRangeUser(0,100);
  S1AreaTemp->GetYaxis()->SetRangeUser(0.0,0.2);

  S1AreaTemp->SetTitle("Fraction of H3 Events below the Mean of WIMP Response");

  ix=0;iy=0; ipad = ny*ix+iy; C1->cd();C1Pads[ipad]->Draw(); C1Pads[ipad]->cd(); S1AreaTemp->Draw(); H3Acceptance->Draw("LP");

  t = new TText(40.0,0.3,"Structure in the distribution is under investigation."); t->SetTextColor(kRed+2); t->Draw();
  

  pngFile = "H3Acceptance.png";  
  C1->SaveAs(pngFile.c_str());

}
