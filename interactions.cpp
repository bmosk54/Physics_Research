#define HiggsAnalysis_cxx
#include "HiggsAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void HiggsAnalysis::Loop()
{
	// In a ROOT session, you can do:
	// Root > .L HiggsAnalysis.C
	// Root > HiggsAnalysis t
	// Root > t.GetEntry(12); // Fill t data members with entry number 12
	// Root > t.Show(); // Show values of entry 12
	// Root > t.Show(16); // Read and show values of entry 16
	// Root > t.Loop(); // Loop on all entries

	// This is the loop skeleton where:
	// jentry is the global entry number in the chain
	// ientry is the entry number in the current Tree
	// Note that the argument to GetEntry must be:
	// jentry for TChain::GetEntry
	// ientry for TTree::GetEntry and TBranch::GetEntry

	// To read only selected branches, Insert statements like:
	// METHOD1:
	// fChain->SetBranchStatus("*",0); // disable all branches
	// fChain->SetBranchStatus("branchname",1); // activate branchname
	// METHOD2: replace line
	// fChain->GetEntry(jentry); //read all branches
	//by b_branchname->GetEntry(ientry); //read only this branch
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	TFile* output = TFile::Open("Dielectron_MC.root", "RECREATE"); // "RECREATE" would produce a new root file with name Dielectron_MC.root every time you run the code

	TTree *tree = new TTree("myTree","A ROOT tree");

	

	TH1F* Z_ee = new TH1F("Z_ee", "Di-electron candidate invariant mass", 200, 0, 200);
	TH1F* H_zz = new TH1F("H_zz", "ZZ candidate invariant mass", 200, 0, 300);
	TH1F* ZH_sum = new TH1F("ZH_sum", "Higgs candidate invariant mass", 200, 0, 300);

	TH1F* el1_pt = new TH1F("el1_pt", "Element 1 pt", 200, 0, 300);
	TH1F* el1_eta = new TH1F("el1_eta", "Element 1 eta", 200, -10, 10);
	TH1F* el1_phi = new TH1F("el1_phi", "Element 1 phi", 200, -10, 10);
	TH1F* el2_pt = new TH1F("el2_pt", "Element 2 pt", 200, 0, 300);
	TH1F* el2_eta = new TH1F("el2_eta", "Element 2 eta", 200, -10, 10);
	TH1F* el2_phi = new TH1F("el2_phi", "Element 2 phi", 200, -10, 10);
	TH1F* el3_pt = new TH1F("el3_pt", "Element 3 pt", 200, 0, 300);
	TH1F* el3_eta = new TH1F("el3_eta", "Element 3 eta", 200, -10, 10);
	TH1F* el3_phi = new TH1F("el3_phi", "Element 3 phi", 200, -10, 10);
	TH1F* el4_pt = new TH1F("el4_pt", "Element 4 pt", 200, 0, 300);
	TH1F* el4_eta = new TH1F("el4_eta", "Element 4 eta", 200, -10, 10);
	TH1F* el4_phi = new TH1F("el4_phi", "Element 4 phi", 200, -10, 10);


	Float_t l1_pt, l1_eta, l1_phi,
				l2_pt, l2_eta, l2_phi,
				l3_pt, l3_eta, l3_phi,
				l4_pt, l4_eta, l4_phi,Z1, Z2, HIGGS;



	tree->Branch("f_lept1_pt", &l1_pt, "l1_pt/F");
	tree->Branch("f_lept1_eta", &l1_eta, "l1_eta/F");
	tree->Branch("f_lept1_phi", &l1_phi, "l1_phi/F");

	tree->Branch("f_lept2_pt", &l2_pt, "l2_pt/F");
	tree->Branch("f_lept2_eta", &l2_eta, "l2_eta/F");
	tree->Branch("f_lept2_phi", &l2_phi, "l2_phi/F");

	tree->Branch("f_lept3_pt", &l3_pt, "l3_pt/F");
	tree->Branch("f_lept3_eta", &l3_eta, "l3_eta/F");
	tree->Branch("f_lept3_phi", &l3_phi, "l3_phi/F");

	tree->Branch("f_lept4_pt", &l4_pt, "l4_pt/F");
	tree->Branch("f_lept4_eta", &l4_eta, "l4_eta/F");
	tree->Branch("f_lept4_phi", &l4_phi, "l4_phi/F");

	tree->Branch("f_Z1mass", &Z1, "f_Z1mass/F");
	tree->Branch("f_Z2mass", &Z2, "f_Z2mass/F");

	tree->Branch("f_mass4l", &HIGGS, "f_mass4l/F");


	
	double el1mt = 0.0;
	double el1pt = 0.0;
	double el1eta = 0.0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry); nbytes += nb;


		 TLorentzVector el1, el2, el3, el4;
		el1.SetPtEtaPhiM(f_lept1_pt, f_lept1_eta, f_lept1_phi, 0.0);
		el2.SetPtEtaPhiM(f_lept2_pt, f_lept2_eta, f_lept2_phi, 0.0);
		TLorentzVector zCandidate = el1 + el2;
		el3.SetPtEtaPhiM(f_lept3_pt, f_lept3_eta, f_lept3_phi, 0.0);
		el4.SetPtEtaPhiM(f_lept4_pt, f_lept4_eta, f_lept4_phi, 0.0);
		TLorentzVector zCandidate2 = el3 + el4;
		 TLorentzVector higgsCandidate = zCandidate + zCandidate2;
		
		 

		HIGGS = higgsCandidate.M();		
		Z1 = zCandidate.M();
		Z2 = zCandidate2.M();

		//el1mt = el1.Mt();
		l1_pt = el1.Mt();
		l1_eta = el1.Eta();
		l1_phi = el1.Phi();

		l2_pt = el2.Mt();
		l2_eta = el2.Eta();
		l2_phi = el2.Phi();

		l3_pt = el3.Mt();
		l3_eta = el3.Eta();
		l3_phi = el3.Phi();

		l4_pt = el4.Mt();
		l4_eta = el4.Eta();
		l4_phi = el4.Phi();

		//cout << el1mt << endl;	

		tree->Fill();
	}

	
	tree->Write();
	
}
