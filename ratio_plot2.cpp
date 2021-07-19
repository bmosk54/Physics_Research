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
TFile *g;


int NS1BINS=20;






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

std::vector<TH1D*> GetTH1DVector(int NBINS,string hname1, TFile *f){

    std::vector<TH1D*> TH1DVector;
    string hname;
    for(int i=0;i<NBINS;i++){
        hname =  hname1+to_string(i);
        TH1DVector.push_back((TH1D*) f->Get(hname.c_str()));
    }
    return TH1DVector;
}

std::vector<TH2D*> GetTH2DVector(int NBINS, string hname1){
    std::vector<TH2D*> TH2DVector;
    string hname;
    for(int i = 0; i < NBINS; i++){
        hname = hname1 + to_string(i);
        TH2DVector.push_back((TH2D*) f->Get(hname.c_str()));
    }

    return TH2DVector;

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

double TritiumDecaySpectrum(double T){

    double Q=18.590;

    if(T>Q) return 0;
    if(T<=0) return 0;

    double m_nue=0.;
    double m_electron=511.0;
    //  double Fermi=1.0;
    double Norm=1.0e-5/1.00121;
    double C=1.0e-3;

    double a1 = sqrt(pow(T,2.)+2*T*m_electron);
    double a2 = T+m_electron;
    double a3 = Q-T;
    double a4 = sqrt(pow(a3,2.)-pow(m_nue,2.));

    double weight = Norm*C*a1*a2*a3*a4*FermiFunction(T);

    return weight;
}

Double_t H3Spectrum(Double_t *x,Double_t *par){
    double T= x[0];
    Double_t fitval = par[0]*TritiumDecaySpectrum(T);
    return fitval;
}



TGraph *getMean(TFile *f, string xVectorHist, string yVectorHist) {
    TGraph *mean;

    //get x and y vectors
    std::vector<TH1D*>  RecoS1AreaVec =  GetTH1DVector(NS1BINS,xVectorHist, f);
    std::vector<TH1D*> RecoS2Area_phdOverS1AreaVec =  GetTH1DVector(16,yVectorHist, f);

    //get x entries
    int npoints = RecoS2Area_phdOverS1AreaVec.size();

    //this is literally the same as the class reference, still dont understand
    const int nq=100;
    double xq[nq];
    for(int i=0;i<nq;i++){ xq[i]= double(i+1)/nq;}
        //S1 mean values over npoints
        double S1Mean[npoints];
        double yq[npoints][nq];


        for(int i=0;i<npoints;i++){
        //set mean value to reconstructed values
        S1Mean[i]= RecoS1AreaVec[i]->GetMean();
        //sets yq to quantile values
        RecoS2Area_phdOverS1AreaVec[i]->GetQuantiles(nq,yq[i],xq);
    }


    double logQ50[npoints];

    for(int i=0;i<npoints;i++){
        logQ50[i] = log10(yq[i][50]);
    }
    cout<<S1Mean[npoints-1]<<endl;
    mean = new TGraph(npoints,S1Mean,logQ50);
    return mean;
}



void ratio_plot2(){

    string hname;
    int ix=0;int iy=0; int ipad =0;
    int nx=1;int ny=1;

    TText *t;

    //get wimpl
    f = new TFile("H3_DataRQ_plotsWIMPCut.root");f->cd();f->ls();
    g = new TFile("H3_DataRQ_plotsNOTWIMPCut.root");g->cd();g->ls();



    
    //get means of graphs and style them
    TGraph *fFileMean = getMean(f, "RecoS1Area_phd_S1Bin", "RecoS2Area_phd_S1Bin");
    TGraph *gFileMean = getMean(g, "RecoS1Area_phd_S1Bin", "RecoS2Area_phd_S1Bin");
    fFileMean->SetLineColor(7);
    gFileMean->SetLineColor(6);
    fFileMean->SetLineWidth(4);
    gFileMean->SetLineWidth(4);
    cout<<gFileMean->GetXaxis()->GetXmax()<<endl;

    //get graph from f graph
    const char* graph= "S1Area_phd_logS2Area_phd_ss";
    TH2D* S1Area_phd_logS2AreaOverS1Area_phd_ss =(TH2D*) f->Get(graph);  
    TH2D* GGS1Area_phd_logS2AreaOverS1Area_phd_ss = (TH2D*) g->Get(graph); 




    TCanvas* C25 = new TCanvas("C25","C25",1000,1000); C25->cd();



    //set axises
    //S1Area_phd_logS2AreaOverS1Area_phd_ss->GetXaxis()->SetRangeUser(0,200);
    //S1Area_phd_logS2AreaOverS1Area_phd_ss->GetYaxis()->SetRangeUser(0,5);
    
    
    //fFileMean->Draw("same");

  


    //TCanvas C_Rat* = new TCanvas("C_Rat","CRat",1000,1000); C_Rat->cd();


    

    /*
    Double_t ratioplotArr[npoints+1];
    for(int i =0;i<npoints+1;i++){ratioplotArr[i]=i;}
    TH1D* ratioPlot = new TH1D("ratioplot", "ratioplot ", npoints, ratioplotArr);
    */

    /*not wimp mean is higher forwhatever that means
    check if graph is higher than that?*/

    Int_t xbins = S1Area_phd_logS2AreaOverS1Area_phd_ss->GetNbinsX();
    xbins = 150; //this is the extent of thje median graph
    Int_t ybins = S1Area_phd_logS2AreaOverS1Area_phd_ss->GetNbinsY();
    int count = 0;
    int total[xbins]; 
    int passed[xbins];
    double ratio[xbins];
    for(int x =0;x<xbins;x++){
        total[x]=0;passed[x]=0;ratio[x]=0.0;
    }


    
    for(Int_t x = 0; x < xbins; x++){
        Double_t delimeter = gFileMean->Eval(x);
        if (x==221){
            cout<<delimeter<<endl;
        }
        for(Int_t y = 0; y<ybins;y++){
            Int_t freq = S1Area_phd_logS2AreaOverS1Area_phd_ss->GetBinContent(x,y);
            int freqAsDouble = (int)freq;
            Int_t bin = S1Area_phd_logS2AreaOverS1Area_phd_ss->GetBin(x,y);
            Double_t someval = S1Area_phd_logS2AreaOverS1Area_phd_ss->GetYaxis()->GetBinUpEdge(y);
            //cout<<someval<< " median: "<< delimeter<<endl;
      
            if (someval>delimeter){
                    
                passed[x]+=freqAsDouble;
                
            }
            total[x]+=freqAsDouble;
        }
    }
   
    for(int x = 0;x<xbins;x++){
        if (total[x]==0){
            ratio[x] = 0;
        }else{
            //cout<<passed[x]<<" "<<total[x]<<endl;
            ratio[x] = (double)(passed[x])/total[x];
            //cout<<"x: "<<x <<": "<<ratio[x]<<endl;
        }
       
        
    }

    S1Area_phd_logS2AreaOverS1Area_phd_ss->Draw("colz");



    //S1Area_phd_logS2AreaOverS1Area_phd_ss->SetBinContent(50,15, 100000);

    //GGS1Area_phd_logS2AreaOverS1Area_phd_ss->Draw("box same");


    gFileMean->Draw("same");

    C25->SaveAs("OverlayWithLines.png");

    

    TCanvas *C_rat = new TCanvas("name", "title", 1000,1000);
    C_rat->SetLogy();
    Double_t ratioplotArr[xbins+1];
    for(int i =0;i<xbins+1;i++){ratioplotArr[i]=i;}
    TH1D* ratioPlot = new TH1D("ratioplot", "ratioplot ", xbins, ratioplotArr);
    for(int i =0;i<xbins;i++){
        //cout<<ratio[i]<<endl;
        ratioPlot->SetBinContent(i, ratio[i]);
        if (ratio[i] != 0)
           ratioPlot->SetBinError(i, sqrt(ratio[i]*(1-ratio[i]))/total[i]); //meaningless error value/10
    }

    ratioPlot->Draw();
    ratioPlot->GetYaxis()->SetRangeUser(0.001,0.1);
    
   
    C_rat->SaveAs("working onit.png");
   
}
