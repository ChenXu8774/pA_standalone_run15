#include "pAuAN.hh"

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <set>
#include <map>
#include <vector>
#include <cmath>

#include "TPDF.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TText.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TAxis.h"
#include "TString.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TPavesText.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TF1.h"

#define  n_phi_bin	16
#define  n_pt_bin 	2 
#define  n_xF_bin	1

#define  Bg1_MinMs  1.5
#define  Bg1_MaxMs  2.4

#define  ncross        120
#define  pi            3.14159265358979323846
#define  twopi         6.28318530717958647692
#define  piby2         1.57079632679489661923
#define  input_spin    "GL2Pdata.root"
#define  input_dimuon  "Run15_pAu.root"
#define  goodrunlist   "run15.list"

using namespace std;

pAuAN::pAuAN()
{
	ReadGoodRun();
	ReadFillTable();
	_file = new TFile(input_dimuon);
	_t = (TTree*)_file->Get("T");
	SetTree();
	cout<<"number of dimuon events:"<<_t->GetEntries()<<endl;
	set<int> fill;

	for(int evt = 0; evt < _t->GetEntries();evt++)
	{
	_t->GetEntry(evt);

	if (evt<20)
	cout<<"test for trigger: "<<(lvl1_trigscaled & 0x00000002)<<"     "<<(lvl1_trigscaled & 0x00100000)<<"		"<<(lvl1_trigscaled & 0x00800000)<<"	"<<(lvl1_trigscaled & 0x00200000)<<endl;	
	if (Cut()) continue;

	fill.insert(_fillTable.find(Run_Number)->second);	
	}
	
	_nFills = fill.size();
	cout << "the number of fills is :"<<_nFills<<endl;
	_Fill = new double[_nFills];

	int count = 0;
	set<int>::iterator itr;
	for(itr = fill.begin(); itr != fill.end(); itr++) _Fill[count++]=*itr;

	GetSpinInfo();
	GetRelLumi();
	GenerateHisto();	
	GetAN();	
}


pAuAN::~pAuAN()
{
	delete [] _Fill;

	map<int, short*>::iterator itr;
	for (itr = _spinB.begin(); itr != _spinB.end(); itr++)
	delete [] itr->second;

	map<int,double*>::iterator itr2;
	for (itr2 = _BBCinB.begin(); itr2 != _BBCinB.end(); itr2++)
	delete [] itr2->second;

	_file->Close();
	delete _file;
}




void pAuAN::ReadGoodRun()
{
	cout<<"Reading the good run"<<endl;
	ifstream fin(goodrunlist);
	int run;
	while (fin>>run) _runList.insert(run);
	_nRuns = _runList.size();

	cout<<"number of good runs: "<<_nRuns<<endl;
}


void pAuAN::ReadFillTable()
{
	//match the run with fill
	int run_number,fill_number;

	cout <<"Read the fill table"<<endl;
	
	ifstream fin("runfilltable.list");
	while (fin>>run_number>>fill_number) 
	{
		_fillTable[run_number]=fill_number;
	//	cout<<_fillTable[run_number]<<"	"<<run_number<<endl;
	}
	fin.close();
	cout<<"fill-run matching finished!"<<endl;
}

	//get spin for each crossing in fill
void pAuAN::GetSpinInfo()
{
	cout<<"Get spin information"<<endl;	
	TFile *fin1 = new TFile(input_spin);
	TTree *t1 = (TTree*)fin1->Get("T");
	
	cout <<"reading...."<< fin1->GetName()<<endl;

	ofstream run_xing_lumi("run_xing_lumi.txt");

	int    runNum, filNum;
	short spinB[ncross];
	float bbc_in[ncross];
	float polB, epolBstat,epolBsyst;
	int XingShift;

	t1->SetBranchAddress("XingShift",  &XingShift);
	t1->SetBranchAddress("RunNumber", &runNum); 
	t1->SetBranchAddress("FillNumber",&filNum);
	t1->SetBranchAddress("SpinB",     spinB);
//	t1->SetBranchAddress("crossing_PHENIX", &crossing);
	t1->SetBranchAddress("BBCin", 	  bbc_in);
	t1->SetBranchAddress("PolB", 	  & polB);
	t1->SetBranchAddress("ePolBstat", &epolBstat); 
	t1->SetBranchAddress("ePolBsyst", &epolBsyst);

	ofstream spin_run("run_spin.txt");
	int xing=0;
	int fill=0, oldfill=0; 
	bool exists = false;
	for (int i =0; i<t1->GetEntries();i++)
	{
		t1->GetEntry(i);

		if(_runList.count(runNum)==0) continue;
	
		fill = filNum;
		exists = false;
		for (int k=0; k<_nFills;k++)
		if(fill==_Fill[k]) {exists=true;break;}		
		
	if (!exists)
	{
	cout<< fill<<"	"<<runNum<<endl;
	continue;
	}

	if (fill !=oldfill)
	{
		oldfill = fill;
		
		_PolB[fill] = polB;
		_ePolB[fill]= float(sqrt(epolBstat*epolBstat + epolBsyst*epolBsyst));
		_spinB[fill] = new short[ncross];
		_BBCinB[fill] = new double[ncross];

		for (int j=0; j<ncross;j++)
		{
			if(j>=XingShift) xing = j-XingShift;			
			else xing = (j+120-XingShift);
			_spinB[fill][xing] = spinB[j];
			spin_run<<runNum<< "	"<<setw(6)<<j<<setw(6)<<spinB[j]<<endl;
			_BBCinB[fill][xing] = 0.;
		}
	}
	
	for (int j=0; j<ncross; j++)
		{
		if(j>=XingShift) xing = j-XingShift;
		else xing = (j+120-XingShift);
		_BBCinB[fill][xing] += double(bbc_in[j]);
		run_xing_lumi<<runNum<<"        "<<setw(5)<<j<<setw(15)<<bbc_in[j]<<endl;
		}
	}
	fin1->Close();
}

void pAuAN::SetTree()
{
  _t->SetBranchAddress("Run_Number", 	    &Run_Number);
  _t->SetBranchAddress("Evt_Number",        &Evt_Number);
  _t->SetBranchAddress("Evt_Nmu",           &Evt_Nmu);
  _t->SetBranchAddress("Evt_bbcZ",          &Evt_bbcZ);
  _t->SetBranchAddress("Evt_vtxchi2",       &Evt_vtxchi2);
  _t->SetBranchAddress("Tr0_DDG0",          &Tr0_DDG0);
  _t->SetBranchAddress("Tr0_DG0",           &Tr0_DG0);
  _t->SetBranchAddress("Tr0_DS3",           &Tr0_DS3);
  _t->SetBranchAddress("Tr0_idchi2",        &Tr0_idchi2);
  _t->SetBranchAddress("Tr0_nidhits",       &Tr0_nidhits);
  _t->SetBranchAddress("Tr0_px",            &Tr0_px);
  _t->SetBranchAddress("Tr0_py",            &Tr0_py);
  _t->SetBranchAddress("Tr0_pz",            &Tr0_pz);
  _t->SetBranchAddress("Tr0_ntrhits",       &Tr0_ntrhits);
  _t->SetBranchAddress("Tr0_lastgap",	    &Tr0_lastgap);
  _t->SetBranchAddress("Tr0_dca_r",         &Tr0_dca_r);
  _t->SetBranchAddress("Tr1_DDG0",          &Tr1_DDG0);
  _t->SetBranchAddress("Tr1_DG0",           &Tr1_DG0);
  _t->SetBranchAddress("Tr1_DS3",           &Tr1_DS3);
  _t->SetBranchAddress("Tr1_idchi2",        &Tr1_idchi2);
  _t->SetBranchAddress("Tr1_nidhits",       &Tr1_nidhits);
  _t->SetBranchAddress("Tr1_px",            &Tr1_px);
  _t->SetBranchAddress("Tr1_py",            &Tr1_py);
  _t->SetBranchAddress("Tr1_pz",            &Tr1_pz);
  _t->SetBranchAddress("Tr1_ntrhits",       &Tr1_ntrhits);
  _t->SetBranchAddress("Tr1_lastgap",       &Tr1_lastgap);
  _t->SetBranchAddress("Tr1_dca_r",	    &Tr1_dca_r);
  _t->SetBranchAddress("charge",            &charge);
  _t->SetBranchAddress("mass",              &mass);
  _t->SetBranchAddress("p",                 &p);
  _t->SetBranchAddress("pT",                &pT);
  _t->SetBranchAddress("px",                  &px);
  _t->SetBranchAddress("py",                &py);
  _t->SetBranchAddress("pz",                &pz);
  _t->SetBranchAddress("rapidity",          &rapidity);
  _t->SetBranchAddress("x1",                &x1);
  _t->SetBranchAddress("x2",                &x2);
  _t->SetBranchAddress("xF",                &xF);
  _t->SetBranchAddress("SpinX_ID",          &SpinX_ID);
  _t->SetBranchAddress("same_event",	    &same_event);
  _t->SetBranchAddress("beamclk0", 	    &beamclk0); 
  _t->SetBranchAddress("lvl1_trigscaled",   &lvl1_trigscaled);

}

void pAuAN::GenerateHisto()
{

	double pT_MinJpsiMs[2][n_pt_bin];
        double pT_MaxJpsiMs[2][n_pt_bin];
	double xF_MinJpsiMs[2];
	double xF_MaxJpsiMs[2];

	double JPSI_SIG_MASS_MAX = -1.;
        double JPSI_SIG_MASS_MIN = -1.;


        pT_MinJpsiMs[North][0] =2.88121;
        pT_MaxJpsiMs[North][0] =3.43513;
        pT_MinJpsiMs[North][1] =2.87472;
        pT_MaxJpsiMs[North][1] =3.45992;

        pT_MinJpsiMs[South][0] =2.84384;
        pT_MaxJpsiMs[South][0] =3.39732;
        pT_MinJpsiMs[South][1] =2.84184;
        pT_MaxJpsiMs[South][1] =3.40799;

	xF_MinJpsiMs[North] = 2.879;
	xF_MaxJpsiMs[North] = 3.423;
	xF_MinJpsiMs[South] = 2.831;
	xF_MaxJpsiMs[South] = 3.405;

	double pt_edge[3] = {0,2,10};
	double xf_edge[2] = {0,10};		

	TCanvas *c[2];
	c[0] = new TCanvas("# of events(inclusive)","# of events(inclusive)",1000,800);
	c[1] = new TCanvas("# of events(background)","# of events(background)",1000,800);

	cout<<"Generating histograms for AN calc"<<endl;

	for (int i =0; i<2; i++)
	{	
		_inc_spinUp[i] = new TH2D(Form("Inc_UP_arm%d",i),Form("Inc_UP_arm%d",i),n_pt_bin,pt_edge,n_phi_bin,-pi,pi); 
		_inc_spinDown[i] = new TH2D(Form("Inc_DN_arm%d",i),Form("Inc_DN_arm%d",i),n_pt_bin,pt_edge,n_phi_bin,-pi,pi);
		_bgr_spinUp[i] = new TH2D(Form("Bgr_UP_arm%d",i),Form("Bgr_UP_arm%d",i),n_pt_bin,pt_edge,n_phi_bin,-pi,pi);
        	_bgr_spinDown[i] = new TH2D(Form("Bgr_DN_arm%d",i),Form("Bgr_DN_arm%d",i),n_pt_bin,pt_edge,n_phi_bin,-pi,pi);

		_xF_inc_spinUp[i] = new TH2D(Form("xF_Inc_UP_arm%d",i),Form("xF_Inc_UP_arm%d",i),n_xF_bin,xf_edge,n_phi_bin,-pi,pi);
                _xF_inc_spinDown[i] = new TH2D(Form("xF_Inc_DN_arm%d",i),Form("xF_Inc_DN_arm%d",i),n_xF_bin,xf_edge,n_phi_bin,-pi,pi);
                _xF_bgr_spinUp[i] = new TH2D(Form("xF_Bgr_UP_arm%d",i),Form("xF_Bgr_UP_arm%d",i),n_xF_bin,xf_edge,n_phi_bin,-pi,pi);
                _xF_bgr_spinDown[i] = new TH2D(Form("xF_Bgr_DN_arm%d",i),Form("xF_Bgr_DN_arm%d",i),n_xF_bin,xf_edge,n_phi_bin,-pi,pi);
	}

	ofstream event("event_output.txt");
	
	TVector3 *vDimu = new TVector3(0.,0.,0.);
	int fill = 0;
	int spin = 0;
	int arm = -9;
	double phi = 0.;

	for (int evt=0; evt<_t->GetEntries();evt++)
	{
		_t->GetEntry(evt);
		if(Cut()) continue;
		vDimu->SetXYZ(px, py, pz);
		phi = double(vDimu->Phi());

		arm = rapidity > 0 ? North : South;	

		fill = _fillTable.find(Run_Number)->second; 
		spin = _spinB[fill][SpinX_ID];

		 event<<
                        setw(15)<<" run num: "<<setw(10)<<Run_Number
                        <<setw(10)<<" beamclk0: "<<setw(15)<<beamclk0
//                        <<setw(10)<<" mass: "<<setw(10)<<mass
//                        <<setw(6)<<" pT: "<<setw(8)<<pT
                        <<setw(8)<<" phi: "<<setw(8)<<phi
                        <<setw(8)<<" Xing: "<<setw(8)<<SpinX_ID
			<<setw(8)<<" spin: "<<setw(8)<<spin 
                        <<endl;

		for(int i=0;i<2;i++)
                        if (pT >pt_edge[i] && pT<pt_edge[i+1])
                        {
                                JPSI_SIG_MASS_MAX = pT_MaxJpsiMs[arm][i];
                                JPSI_SIG_MASS_MIN = pT_MinJpsiMs[arm][i];
                        }
		if (mass > JPSI_SIG_MASS_MIN && mass<JPSI_SIG_MASS_MAX && charge==0)
		{
			if(spin==UP) _inc_spinUp[arm]->Fill(pT,phi);
			if(spin==DN) _inc_spinDown[arm]->Fill(pT,phi);
		}

		if (mass > Bg1_MinMs && mass < Bg1_MaxMs && charge == 0)
		{
			if(spin==UP) _bgr_spinUp[arm]->Fill(pT,phi);
                        if(spin==DN) _bgr_spinDown[arm]->Fill(pT,phi);
		}
		if (mass > xF_MinJpsiMs[arm] && mass < xF_MaxJpsiMs[arm])
		{
			if(spin==UP) _xF_inc_spinUp[arm]->Fill(abs(xF),phi);
			if(spin==DN) _xF_inc_spinDown[arm]->Fill(abs(xF),phi);
		}
		if (mass > Bg1_MinMs && mass < Bg1_MaxMs)
                {
                        if(spin==UP) _xF_bgr_spinUp[arm]->Fill(abs(xF),phi);
			if(spin==DN) _xF_bgr_spinDown[arm]->Fill(abs(xF),phi);
		}
		

	}

	ofstream n_hit_xf("hit_plots_xf_s_sig.txt");
	ofstream n_hit0("hit_plots_s_sig.txt");
        ofstream n_hit1("hit_plots_n_sig.txt");
	for (int k=0; k<n_phi_bin; k++)
	{
		n_hit0<<"event:"<<setw(3)<<k<<" : "
		<<" UP "<<_inc_spinUp[1]->GetBinContent(1,k+1)
		<<" DN "<<_inc_spinDown[1]->GetBinContent(1,k+1)
		<<endl;
		n_hit1<<"event:"<<setw(3)<<k<<" : "
		<<" UP "<<_inc_spinUp[0]->GetBinContent(1,k+1)
                <<" DN "<<_inc_spinDown[0]->GetBinContent(1,k+1)
		<<endl;
		n_hit_xf<<"event:"<<setw(3)<<k<<" : "
		<<" UP "<<_xF_inc_spinUp[1]->GetBinContent(1,k+1)
		<<" DN "<<_xF_inc_spinDown[1]->GetBinContent(1,k+1)
		<<endl;
	}	

	for (int i =0; i <2; i++)
	{

		_inc_spinUp[i]->GetXaxis()->SetTitle("pT");
		_inc_spinUp[i]->GetYaxis()->SetTitle("phi");
		_inc_spinDown[i]->GetXaxis()->SetTitle("pT");
        	_inc_spinDown[i]->GetYaxis()->SetTitle("phi");
		_bgr_spinUp[i]->GetXaxis()->SetTitle("pT");
        	_bgr_spinUp[i]->GetYaxis()->SetTitle("phi");
        	_bgr_spinDown[i]->GetXaxis()->SetTitle("pT");
        	_bgr_spinDown[i]->GetYaxis()->SetTitle("phi");
		c[i]->Divide(2,2);
		c[i]->cd(1);
		_inc_spinUp[i]->Draw("colz");
		c[i]->cd(2);
		_inc_spinDown[i]->Draw("colz");
		c[i]->cd(3);
        	_bgr_spinUp[i]->Draw("colz");
        	c[i]->cd(4);
        	_bgr_spinDown[i]->Draw("colz");
	}
}

void pAuAN::GetRelLumi()
{

	double L_up;
	double L_down;
	double L_up_total=0;
	double L_down_total=0;	
	double L_total=0;
	double bbcin;
	double weight[_nFills];
	short spinB;
	int fill;

	ofstream lumi_fill("lumi_fill.txt");
	
	for(int i=0;i<_nFills;i++)
	{
		weight[i]=0;
		L_up=0;
		L_down=0;
		fill = int(_Fill[i]);
		for(int j=0;j<ncross;j++)
		{
			bbcin = _BBCinB.find(fill)->second[j];
			spinB = _spinB.find(fill)->second[j];
			if(spinB == UP) L_up += bbcin;
			if(spinB == DN) L_down += bbcin;
			weight[i] = weight[i] + bbcin;
		}

		L_up_total    += L_up;
		L_down_total  += L_down;
	}

	_R  = L_up_total/L_down_total;
	_eR = _R*sqrt(1./L_up_total+1./L_down_total);
	cout<<"total Lumi for spin  1 is:"<<L_up_total<<endl;
	cout<<"total Lumi for spin -1 is:"<<L_down_total<<endl;
	cout<<"relative Lumi:"<<_R<<endl;
	cout<<"error for relative Lumi:"<<_eR<<endl;

	_P = 0;
	_eP= 0;

	ofstream polar_fill("polar_fill.txt");

	for(int i = 0; i < _nFills; i++)
	{
		fill = int(_Fill[i]);
		_P = _P+weight[i]* (_PolB.find(fill)->second);
		_eP= _eP+weight[i]*(_ePolB.find(fill)->second);
		L_total += weight[i];
		lumi_fill<<setw(6)<<"fill: "<<fill<<setw(8)<< "weight: "<<weight[i]<<"  polB  "<<_PolB.find(fill)->second<<endl;
	}
	_P = _P/L_total;
	_eP = _eP/L_total;
	cout<<"total Polarization is:"<<_P<<endl;
	cout<<"error for total polarization is:"<<_eP<<endl;	 
}

void pAuAN::GetAN()
{
	gStyle->SetOptFit();	

	double n_up_inc[n_phi_bin], n_down_inc[n_phi_bin];
	double n_up_bgr[n_phi_bin], n_down_bgr[n_phi_bin];

	double xf_n_up_inc, xf_n_down_inc;
        double xf_n_up_bgr, xf_n_down_bgr;

	double A_L_inc[2][n_pt_bin][n_phi_bin], eA_L_inc[2][n_pt_bin][n_phi_bin];
	double A_L_bgr[2][n_pt_bin][n_phi_bin],eA_L_bgr[2][n_pt_bin][n_phi_bin];	

	double xF_A_L_inc[2][n_phi_bin], xF_eA_L_inc[2][n_phi_bin];
        double xF_A_L_bgr[2][n_phi_bin], xF_eA_L_bgr[2][n_phi_bin];

	double phi[n_phi_bin];
	TGraphErrors *Bgr_A_L[2][n_pt_bin];
	TGraphErrors *Inc_A_L[2][n_pt_bin];

	TGraphErrors *xF_Bgr_A_L[2];
        TGraphErrors *xF_Inc_A_L[2];
	ofstream al_check("al_check.txt");

	for (int k = 0; k<2; k++){ 
	for (int i = 0; i<n_pt_bin; i++)
	{
		for (int j = 0; j<n_phi_bin; j++)
		{
			n_up_inc[j] = _inc_spinUp[k]->GetBinContent(i+1,j+1);
			n_down_inc[j] = _inc_spinDown[k]->GetBinContent(i+1,j+1);
			n_up_bgr[j] = _bgr_spinUp[k]->GetBinContent(i+1,j+1);
			n_down_bgr[j] = _bgr_spinDown[k]->GetBinContent(i+1,j+1);	

			if(k ==0&&i ==0)
			phi[j] = _inc_spinUp[k]->GetYaxis()->GetBinCenter(j+1);
	
			A_L_inc[k][i][j] = 1./_P*(n_up_inc[j]-_R*n_down_inc[j])/(n_up_inc[j]+_R*n_down_inc[j]);
			A_L_bgr[k][i][j] =1./_P*(n_up_bgr[j]-_R*n_down_bgr[j])/(n_up_bgr[j]+_R*n_down_bgr[j]);  

			if (i == 0 && k ==0)
			{
				al_check<<"bin "<<setw(3)<<j<<" : "
				<<setw(15)<<A_L_inc[k][i][j]
				<<setw(15)<<A_L_bgr[k][i][j]
				<<endl;
			}
			eA_L_inc[k][i][j] =sqrt(4.*(n_up_inc[j]*_R*_R*n_down_inc[j]*n_down_inc[j]+_eR*_eR*_R*_R*n_up_inc[j]*n_up_inc[j]+n_down_inc[j]*_R*_R*n_up_inc[j]*n_up_inc[j])/(n_up_inc[j]+_R*n_down_inc[j])/(n_up_inc[j]+_R*n_down_inc[j])/(n_up_inc[j]+_R*n_down_inc[j])/(n_up_inc[j]+_R*n_down_inc[j])/_P/_P+1./_P/_P/_P/_P*(n_up_inc[j]-_R*n_down_inc[j])*(n_up_inc[j]-_R*n_down_inc[j])/(n_up_inc[j]+_R*n_down_inc[j])/(n_up_inc[j]+_R*n_down_inc[j])*_eP*_eP);
			eA_L_bgr[k][i][j] =sqrt(4.*(n_up_bgr[j]*_R*_R*n_down_bgr[j]*n_down_bgr[j]+_eR*_eR*_R*_R*n_up_bgr[j]*n_up_bgr[j]+n_down_bgr[j]*_R*_R*n_up_bgr[j]*n_up_bgr[j])/(n_up_bgr[j]+_R*n_down_bgr[j])/(n_up_bgr[j]+_R*n_down_bgr[j])/(n_up_bgr[j]+_R*n_down_bgr[j])/(n_up_bgr[j]+_R*n_down_bgr[j])/_P/_P+1./_P/_P/_P/_P*(n_up_bgr[j]-_R*n_down_bgr[j])*(n_up_bgr[j]-_R*n_down_bgr[j])/(n_up_bgr[j]+_R*n_down_bgr[j])/(n_up_bgr[j]+_R*n_down_bgr[j])*_eP*_eP);
	}
		Inc_A_L[k][i] = new TGraphErrors(n_phi_bin,phi,A_L_inc[k][i],0,eA_L_inc[k][i]);
		Inc_A_L[k][i]->SetTitle(Form("Inc_pT_bin%d_arm%d",i,k));
		Inc_A_L[k][i]->GetXaxis()->SetTitle("phi");
		Inc_A_L[k][i]->GetYaxis()->SetTitle("AL");
		Inc_A_L[k][i]->SetMinimum(-0.4);
		Inc_A_L[k][i]->SetMaximum(0.4);
		Bgr_A_L[k][i] = new TGraphErrors(n_phi_bin,phi,A_L_bgr[k][i],0,eA_L_bgr[k][i]);
		Bgr_A_L[k][i]->SetTitle(Form("Bgr_pT_bin%d_arm%d",i,k));
        	Bgr_A_L[k][i]->GetXaxis()->SetTitle("phi");
        	Bgr_A_L[k][i]->GetYaxis()->SetTitle("AL");
		Bgr_A_L[k][i]->SetMinimum(-0.4);
                Bgr_A_L[k][i]->SetMaximum(0.4);
	}
				}
	ofstream xf_AN("xf_AN_n.txt");
	//xF_AN
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < n_phi_bin; j++)
		{
			xf_n_up_inc = _xF_inc_spinUp[i]->GetBinContent(1,j+1);
			xf_n_down_inc = _xF_inc_spinDown[i]->GetBinContent(1,j+1);
			xf_n_up_bgr = _xF_bgr_spinUp[i]->GetBinContent(1,j+1);
			xf_n_down_bgr = _xF_bgr_spinDown[i]->GetBinContent(1,j+1);
		
			xF_A_L_inc[i][j] = 1./_P*(xf_n_up_inc-_R*xf_n_down_inc)/(xf_n_up_inc+_R*xf_n_down_inc);
			xF_A_L_bgr[i][j] = 1./_P*(xf_n_up_bgr-_R*xf_n_down_bgr)/(xf_n_up_bgr+_R*xf_n_down_bgr);

			xF_eA_L_inc[i][j] =sqrt(4.*(xf_n_up_inc*_R*_R*xf_n_down_inc*xf_n_down_inc+_eR*_eR*_R*_R*xf_n_up_inc*xf_n_up_inc+xf_n_down_inc*_R*_R*xf_n_up_inc*xf_n_up_inc)/(xf_n_up_inc+_R*xf_n_down_inc)/(xf_n_up_inc+_R*xf_n_down_inc)/(xf_n_up_inc+_R*xf_n_down_inc)/(xf_n_up_inc+_R*xf_n_down_inc)/_P/_P+1./_P/_P/_P/_P*(xf_n_up_inc-_R*xf_n_down_inc)*(xf_n_up_inc-_R*xf_n_down_inc)/(xf_n_up_inc+_R*xf_n_down_inc)/(xf_n_up_inc+_R*xf_n_down_inc)*_eP*_eP);
                        xF_eA_L_bgr[i][j] =sqrt(4.*(xf_n_up_bgr*_R*_R*xf_n_down_bgr*xf_n_down_bgr+_eR*_eR*_R*_R*xf_n_up_bgr*xf_n_up_bgr+xf_n_down_bgr*_R*_R*xf_n_up_bgr*xf_n_up_bgr)/(xf_n_up_bgr+_R*xf_n_down_bgr)/(xf_n_up_bgr+_R*xf_n_down_bgr)/(xf_n_up_bgr+_R*xf_n_down_bgr)/(xf_n_up_bgr+_R*xf_n_down_bgr)/_P/_P+1./_P/_P/_P/_P*(xf_n_up_bgr-_R*xf_n_down_bgr)*(xf_n_up_bgr-_R*xf_n_down_bgr)/(xf_n_up_bgr+_R*xf_n_down_bgr)/(xf_n_up_bgr+_R*xf_n_down_bgr)*_eP*_eP);
			if(i == 1)
			{
			xf_AN<<"bin "<<setw(3)<<j<<" : "
			<<setw(15)<<xF_A_L_inc[i][j]
			<<setw(15)<<xF_A_L_bgr[i][j]
			<<endl;
			}
		}	

		xF_Inc_A_L[i] = new TGraphErrors(n_phi_bin,phi,xF_A_L_inc[i],0,xF_eA_L_inc[i]);
		xF_Inc_A_L[i]->GetXaxis()->SetTitle("phi");
		xF_Inc_A_L[i]->GetYaxis()->SetTitle("AL");
		xF_Inc_A_L[i]->SetMinimum(-0.4);
		xF_Inc_A_L[i]->SetMaximum(0.4);
		xF_Inc_A_L[i]->SetTitle(Form("AN vs xF(inclusive), Arm%d",i)); 

		xF_Bgr_A_L[i] = new TGraphErrors(n_phi_bin,phi,xF_A_L_bgr[i],0,xF_eA_L_bgr[i]);
                xF_Bgr_A_L[i]->GetXaxis()->SetTitle("phi");
                xF_Bgr_A_L[i]->GetYaxis()->SetTitle("AL");
                xF_Bgr_A_L[i]->SetMinimum(-0.4);
                xF_Bgr_A_L[i]->SetMaximum(0.4);
                xF_Bgr_A_L[i]->SetTitle(Form("AN vs xF(background), Arm%d",i));	
	}
	TF1 *f_inc[2][n_pt_bin];//[arm] 
	TF1 *f_bgr[2][n_pt_bin];

	TF1 *xF_f_inc[2];
	TF1 *xF_f_bgr[2];	

	for(int i=0; i<2; i++)
		for(int j=0; j<n_pt_bin; j++)		
		{
			f_inc[i][j] = new TF1(Form("f_Inc_pT_bin%d_arm%d",i,j),"[0]*TMath::Cos(x)+[1]*TMath::Sin(x)",-pi,pi);
			f_bgr[i][j] = new TF1(Form("f_Bgr_pT_bin%d_arm%d",i,j),"[0]*TMath::Cos(x)+[1]*TMath::Sin(x)",-pi,pi);
		}
 	
	for(int i=0; i<2; i++)
        {
		xF_f_inc[i] = new TF1(Form("f_Inc_xF_arm%d",i),"[0]*TMath::Cos(x)+[1]*TMath::Sin(x)",-pi,pi);
                xF_f_bgr[i] = new TF1(Form("f_Bgr_xF_arm%d",i),"[0]*TMath::Cos(x)+[1]*TMath::Sin(x)",-pi,pi);
	}

	for(int i =0; i<2; i++)
		for(int j=0; j<n_pt_bin; j++)
		{
			Inc_A_L[i][j]->Fit(f_inc[i][j]);
			Bgr_A_L[i][j]->Fit(f_bgr[i][j]);
		}

	TCanvas *c2[2];
	for (int i =0; i < n_pt_bin; i++)
	{
		c2[i] = new TCanvas(Form("pT_BIN_%d",i),Form("pT_BIN_%d",i),1000,800);
		c2[i]->Divide(2,2);
		c2[i]->cd(1);
		Inc_A_L[0][i]->Draw("AP");
		f_inc[0][i]->Draw("SAME");
		c2[i]->cd(2);
		Bgr_A_L[0][i]->Draw("AP");
		f_bgr[0][i]->Draw("SAME");
 		c2[i]->cd(3);
		Inc_A_L[1][i]->Draw("AP");
                f_inc[1][i]->Draw("SAME");
		c2[i]->cd(4);
		Bgr_A_L[1][i]->Draw("AP");
                f_bgr[1][i]->Draw("SAME");
		c2[i]->Draw();
	}

	for(int i =0; i<2; i++)
        {
                xF_Inc_A_L[i]->Fit(xF_f_inc[i]);
                xF_Bgr_A_L[i]->Fit(xF_f_bgr[i]);
        }
	
	TCanvas *c3;
	c3 = new TCanvas("AN vs xF fit","AN vs xF fit",1000,800);
	c3->Divide(2,2);
	c3->cd(1);
	xF_Inc_A_L[0]->Draw("AP");
        xF_f_inc[0]->Draw("SAME");
	c3->cd(2);
	xF_Inc_A_L[1]->Draw("AP");
        xF_f_inc[1]->Draw("SAME");
	c3->cd(3);
	xF_Bgr_A_L[0]->Draw("AP");
        xF_f_bgr[0]->Draw("SAME");
	c3->cd(4);
	xF_Bgr_A_L[1]->Draw("AP");
        xF_f_bgr[1]->Draw("SAME");
}

bool pAuAN::Cut()
{
	if (_runList.count(Run_Number) == 0) return true;

        if (same_event != true) return true;
        if (mass<=1.0 || mass>=4.0) return true;
        if (Tr0_pz*Tr1_pz < 0.) return true;
        if (charge !=0) return true;
        if (fabs(Evt_bbcZ) == 0.) return true;
        if (fabs(Evt_bbcZ) >= 30.) return true;

        if (Tr0_pz >0 && (fabs(Tr0_DG0) >= 25. || fabs(Tr1_DG0) >= 25.)) return true;
        if (Tr0_pz <0 && (fabs(Tr0_DG0) >= 30. || fabs(Tr1_DG0) >= 30.)) return true;

        if (Tr0_DDG0 >= 10. || Tr1_DDG0 >= 10.) return true;
        if (Tr0_ntrhits<=9 || Tr1_ntrhits<=9) return true;
        if (Tr0_nidhits<=5 || Tr1_nidhits<=5) return true;
        if (Tr0_lastgap<=2 || Tr1_lastgap<=2) return true;
        if (Tr0_dca_r>5. || Tr1_dca_r>5.) return true;

        if ( pT >= 10.) return true;

        if (pz >= 100.) return true;
        if (fabs(rapidity) <= 1.2 || fabs(rapidity) >= 2.2) return true;

/*        if (Tr0_pz >0) 
	{
		if (!((lvl1_trigscaled & 0x00400000)||(lvl1_trigscaled &0x00100000))) 
		return true;
	}

        if (Tr0_pz <0)
	{
		if (!((lvl1_trigscaled & 0x00800000)||(lvl1_trigscaled &0x00200000))) 
		return true;
	}
*/
        if (Evt_vtxchi2 >= 5.) return true;

        return false;	
}
