/** \file main.cpp
 * \brief This is the main routine of the cluster reconstruction algorithm. It reads a macro file, imports its parameters, performs the cluster reconstruction algorithm and outputs the results.
 *
 * \author Sebastian Templ <sebastian.templ@gmail.com>
 * \date 2016
 */


#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <time.h>

#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TError.h>
#include <TF1.h>
#include <TColor.h>

#include "CRprintInfo.h"

#include "CREntry.h"
#include "CREvent.h"

#include "langaufit.h"
#include "langaufun.h"
#include "langaupro.h"


int main(int argc, char** argv) {
  TString argument;
  
  TString macroFileName;
  TString filename, calibrationFileName;

  float SNRThreshold;
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4], fr[2]; // parameters for Landau*Gauss fit

  // default values:
  TString outFileName = "output_";
  size_t nBins = 256;
  unsigned int channel_startPos = 1;
  unsigned int channel_endPos   = 256;
  unsigned int timeCut_lower = 0;
  unsigned int timeCut_upper = 25;
  unsigned int maxClusterSize = 1;
  int scaleCutValue = 10;
  bool isMacroFile = false;
  bool doConvertADC = false;
  bool doSubtractBackground = true;
  bool doCommonModeCorrection = true;
  bool doSetSNRThreshold = false;
  bool doLangausFit = false;
  bool isLangausFitRange_start = false;
  bool isLangausFitRange_end = false;
  bool isLangausXiLower = false;
  bool isLangausXiUpper = false;
  bool isLangausMPVLower = false;
  bool isLangausMPVUpper = false;
  bool isLangausAreaLower = false;
  bool isLangausAreaUpper = false;
  bool isLangausSigmaGLower = false;
  bool isLangausSigmaGUpper = false;
  bool isLangausXiInitial = false;
  bool isLangausMPVInitial = false;
  bool isLangausAreaInitial = false;
  bool isLangausSigmaGInitial = false;
  bool isLangausRanges = false;
  bool isLangausStartValues = false;
  bool isNoiseCutValue = false;
  bool isSetTimeCuts = false;
  bool isSetTimeCut_lower = false;
  bool isSetTimeCut_upper = false;

  if(argc<2) {
    std::cout << "   ERROR: No arguments provided. State the name of a macro file!" << std::endl;
    return 1;
  }else if(argc>2) {
    std::cout << "   ERROR: Too many arguments provided. State the name of a macro file!" << std::endl;
    return 1;
  }else { // i.e., if(argc==2)
    argument = argv[1];
    if(argument.EndsWith(".root")) {
      std::cout << "   ERROR: Invalid file input. Only macro files can be passed." << std::endl;
      return 1;
    }else {
      std::cout << argument << " set as macro file..." << std::endl;
      isMacroFile = true;
      macroFileName = argument;
     }
  }

  if(isMacroFile) {
    std::ifstream macroFile(macroFileName);
    if(!macroFile) {
      std::cout << "   ERROR: Macro file " << macroFileName << " cannot be opened." << std::endl;
      return 2;
    }
    std::cout << "Setting parameters now..." << std::endl << std::endl;
    std::string currentLine_temp;
    while(std::getline(macroFile, currentLine_temp)) {
      TString currentLine = currentLine_temp;
      if(currentLine.BeginsWith("#")) continue;
      if(currentLine_temp.empty()) continue;
      currentLine.ReplaceAll(" ", "");
      TString identifier = currentLine;
      identifier.Remove(identifier.First('='), identifier.Length());
      TString value = currentLine;
      value.Remove(0, value.Last('=')+1);
      value = value.Strip(TString::kBoth, '\"');

      if(identifier == "ROOTFILENAME") {
	filename = value;
	if(!filename.EndsWith(".root")) filename += ".root";
	std::cout << "File " << filename << " set as input ROOT file..." << std::endl;
      }
      if(identifier == "CALIBRATIONFILENAME") {
	calibrationFileName = value;
	if(calibrationFileName != "none") {
	  doConvertADC = true;
	  if(!calibrationFileName.EndsWith(".root")) calibrationFileName += ".root";
	  std::cout << "File " << calibrationFileName << " set as calibration file..." << std::endl;
	}else {
	  doConvertADC = false;
	}
      }
      if(identifier == "CHANNEL_START") {
	if(value.IsDigit()) {
	  channel_startPos = value.Atoi();
	  //std::cout << "First channel number set to " << channel_startPos << "..." << std::endl;
	}else {
	  std::cout << "   ERROR: 'CHANNEL_START = " << value << "' is an invalid input." << std::endl;
	  return 3;
	}
      }
      if(identifier == "CHANNEL_END") {
	if(value.IsDigit()) {
	  channel_endPos = value.Atoi();
	  //std::cout << "Last channel number set to " << channel_endPos << "..." << std::endl;
	}else {
	  std::cout << "   ERROR: 'CHANNEL_END = " << value << "' is an invalid input." << std::endl;
	  return 3;
	}
      }
      if(identifier == "NBINS") {
	if(value.IsDigit()) {
	  nBins = value.Atoi();
	  //std::cout << "Total number of bins set to " << nBins << "..." << std::endl;
	}else {
	  std::cout << "   WARNING: 'NBINS = " << value << "' is an invalid input. The default value will be used." << std::endl;
	}
      }
      if(identifier == "TIMECUT_LOWER") {
	if(value.IsDigit()) {
	  isSetTimeCuts = true;
	  isSetTimeCut_lower = true;
	  timeCut_lower = value.Atoi();
	  //std::cout << "Lower time cut value set to " << timeCut_lower << "..." << std::endl;
	  if(isSetTimeCut_upper && (timeCut_lower >= timeCut_upper)) {
	    std::cout << "   ERROR: Invalid input for the time cut values." << std::endl;
	    return 3;
	  }
	}else {
	  std::cout << "   WARNING: 'TIMECUT_LOWER = " << value << "' is an invalid input. The default value will be used." << std::endl;
	}
      }
      if(identifier == "TIMECUT_UPPER") {
	if(value.IsDigit()) {
	  isSetTimeCuts = true;
	  isSetTimeCut_upper = true;
	  timeCut_upper = value.Atoi();
	  //std::cout << "Upper time cut value set to " << timeCut_upper << "..." << std::endl;
	  if(isSetTimeCut_lower && (timeCut_lower >= timeCut_upper)) {
	    std::cout << "   ERROR: Invalid input for the time cut values." << std::endl;
	    return 3;
	  }
	}else {
	  std::cout << "   WARNING: 'TIMECUT_UPPER = " << value << "' is an invalid input. The default value will be used." << std::endl;
	}
      }
      if(identifier == "MAXCLUSTERSIZE") {
	if(value.IsDigit()) {
	  maxClusterSize = value.Atoi();
	  //std::cout << "Maximum cluster size set to " << maxClusterSize << "..." << std::endl;
	}else {
	  std::cout << "   ERROR: 'MAXCLUSTERSIZE = " << value << "' is an invalid input." << std::endl;
	  return 3;
	}
      }
      if(identifier == "COMMONMODECORRECTION") {
	value.ToLower();
	if(value == "yes") {
	  doCommonModeCorrection = true;
	  //std::cout << "A common-mode correction will be performed..." << std::endl;
	}else if(value == "no") {
	  doCommonModeCorrection = false;
	  //std::cout << "A common-mode correction won't be performed..." << std::endl;
	}else {
	  std::cout << "   WARNING: 'COMMONMODECORRECTION = " << value << "' is an invalid input. The default value will be used." << std::endl;
	}
      }
      if(identifier == "SNRTHRESHOLD") {
	if(value.IsFloat()) {
	  doSetSNRThreshold = true;
	  SNRThreshold = value.Atof();
	  //std::cout << "SNR threshold set to " << SNRThreshold << "..." << std::endl;
	}else {
	  std::cout << "   WARNING: 'SNRTHRESHOLD = " << value << "' is an invalid input. No threshold set." << std::endl;
	}
      }
      if(identifier == "SUBTRACTBACKGROUND") {
	value.ToLower();
	if(value == "yes") {
	  doSubtractBackground = true;
	  //std::cout << "The estimated background will be subtracted..." << std::endl;
	}else if(value == "no") {
	  doSubtractBackground = false;
	  //std::cout << "The estimated background won't be subtracted..." << std::endl;
	}else {
	  std::cout << "   WARNING: 'SUBTRACTBACKGROUND = " << value << "' is an invalid input. The default value will be used." << std::endl;
	}
      }
      if(identifier == "SCALECUTVALUE") {
	if(value.IsDigit()) {
	  scaleCutValue = value.Atoi();
	  //std::cout << "Scale cut value set to " << scaleCutValue << "..." << std::endl;
	}else {
	  std::cout << "   WARNING: 'SCALECUTVALUE = " << value << "' is an invalid input. The default value will be used." << std::endl;
	}
      }
      if(identifier == "LANGAUSFIT") {
	value.ToLower();
	if(value == "yes") {
	  doLangausFit = true;
	  //std::cout << "A fit with a Landau-Gauss convolution function will be performed..." << std::endl;
	}else if(value == "no") {
	  doLangausFit = false;
	  //std::cout << "A fit with a Landau-Gauss convolution function won't be performed..." << std::endl;
	}else {
	  std::cout << "   WARNING: 'LANGAUSFIT = " << value << "' is an invalid input. The default value will be used." << std::endl;
	}
      }
      if(identifier.BeginsWith("LANGAUS")) {
	std::string value_str;
	value_str = value;
	if(identifier == "LANGAUS_FITRANGE_START") {
	  isLangausFitRange_start = true;
	  fr[0] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_FITRANGE_END") {
	  isLangausFitRange_end = true;
	  fr[1] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_XI_LOWER") {
	  isLangausXiLower = true;
	  pllo[0] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_XI_UPPER") {
	  isLangausXiUpper = true;
	  plhi[0] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_MPV_LOWER") {
	  isLangausMPVLower = true;
	  pllo[1] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_MPV_UPPER") {
	  isLangausMPVUpper = true;
	  plhi[1] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_AREA_LOWER") {
	  isLangausAreaLower = true;
	  pllo[2] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_AREA_UPPER") {
	  isLangausAreaUpper = true;
	  plhi[2] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_SIGMAG_LOWER") {
	  isLangausSigmaGLower = true;
	  pllo[3] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_SIGMAG_UPPER") {
	  isLangausSigmaGUpper = true;
	  plhi[3] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_XI_INITIAL") {
	  isLangausXiInitial = true;
	  sv[0] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_MPV_INITIAL") {
	  isLangausMPVInitial = true;
	  sv[1] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_AREA_INITIAL") {
	  isLangausAreaInitial = true;
	  sv[2] = std::atof(value_str.c_str());
	}
	if(identifier == "LANGAUS_SIGMAG_INITIAL") {
	  isLangausSigmaGInitial = true;
	  sv[3] = std::atof(value_str.c_str());
	}
      }
    }
  }
  if(isLangausXiLower && isLangausXiUpper &&
     isLangausMPVLower && isLangausMPVUpper &&
     isLangausAreaLower && isLangausAreaUpper &&
     isLangausSigmaGLower && isLangausSigmaGUpper) {
    isLangausRanges = true;
  }
  if(isLangausXiInitial && isLangausMPVInitial && isLangausAreaInitial && isLangausSigmaGInitial) {
    isLangausStartValues = true;
  }
  if(doLangausFit) {
    if(isLangausRanges && isLangausStartValues) {
      //std::cout << "Parameters for the Landau*Gauss fit set..." << std::endl;
    }else if(!isLangausRanges) {
      std::cout << "   WARNING: Not all ranges for the Landau*Gauss fit were set properly. Default values will be used." << std::endl;
    }else if(!isLangausStartValues) {
      std::cout << "   WARNING: Not all start values for the Landau*Gauss fit were set properly. Default values will be used." << std::endl;
    }
  }
  std::cout << "All parameters set." << std::endl << std::endl;
  std::cout << "Summary of input values:" << std::endl;


  // output status and parameters:
  CRprintInfo(doConvertADC, "ADC->ke conversion using calibration file");
  CRprintInfo(doCommonModeCorrection, "Common-mode correction");
  CRprintInfo(doSubtractBackground, "Background subtraction");
  CRprintInfo(doLangausFit, "Landau*Gauss fit to signal distribution");
  if(isLangausFitRange_start) CRprintInfo((float)fr[0], "Landau*Gauss fit range start");
  if(isLangausFitRange_end) CRprintInfo((float)fr[1], "Landau*Gauss fit range end");
  CRprintInfo((doLangausFit && isLangausRanges && isLangausStartValues), "All parameters for the Landau*Gauss fit set successfully");
  CRprintInfo(doSetSNRThreshold, "SNR threshold setting");
  if(doSetSNRThreshold) {
    CRprintInfo(SNRThreshold, "SNR threshold value");
  }
  CRprintInfo(isSetTimeCuts, "Time cuts setting");
  if(isSetTimeCuts) {
    CRprintInfo(timeCut_lower, "Lower time cut value");
    CRprintInfo(timeCut_upper, "Upper time cut value");
  }
  CRprintInfo(channel_startPos, "First channel number");
  CRprintInfo(channel_endPos, "Last channel number");
  CRprintInfo((int)nBins, "Total number of bins");
  CRprintInfo(maxClusterSize, "Maximally-allowed cluster multiplicity");
  CRprintInfo(scaleCutValue, "Scale cut value");
  

  ///////////////////////////////////////////////////////////////////
  
  const unsigned int nChannels = channel_endPos-channel_startPos+1;
  const int markerVal = -1e4; // for marking unused bins/entries
  const int min       = -30,
            max       = 150,
            binnumber = (abs(min)+abs(max));

  TFile* file = TFile::Open(filename);
  TTree* tree = (TTree*) file->Get("tree");

  TTree* treeCalibr;
  float slope_calibr;
  if(doConvertADC) {
    TFile* fileCalibr = TFile::Open(calibrationFileName);
    treeCalibr = (TTree*) fileCalibr->Get("treeCalibr");
    treeCalibr->SetBranchAddress("slope",&slope_calibr);
  }

  float amplitude[nBins];
  float time;
  tree->SetBranchAddress("amplitude", amplitude);
  tree->SetBranchAddress("time",&time);
  float median;

  std::cout << "Input files read..." << std::endl;
  std::cout << "Making preparations...";

  TH1F* hist_noise = (TH1F*) file->Get("histNoise");
  TH1F* hist_spectrum = new TH1F("hist_spectrum","",binnumber,min,max);
  TH1F* hist_background = new TH1F("hist_background","",binnumber,min,max);

  TH1F *hist_clusterSizeDistr = new TH1F("hist_clusterSizeDistr","",maxClusterSize,1,maxClusterSize+1);
  TH1F *hist_clusterSizeDistr_background = new TH1F("hist_clusterSizeDistr_background","",maxClusterSize,1,maxClusterSize+1);

  std::vector<CREvent> CREventCollection;

  // Import and store ADC->ke conversion values and convert the noise distribution, if the option is enabled:
  float slopePerChannel[nBins];
  if(doConvertADC) {
    for(unsigned int i=0; i<nBins; i++) {
      treeCalibr->GetEvent(i);
      slopePerChannel[i] = slope_calibr;
      hist_noise->SetBinContent(i+1, hist_noise->GetBinContent(i+1)/slopePerChannel[i]);
    }
  }

 ///////////////////////////////////////////////////////////////////

  unsigned int nEv = tree->GetEntries();
  std::cout << "\rMaking preparations... DONE" << std::endl;
  std::cout << "Starting cluster reconstruction algorithm..." << std::endl << std::endl;
  clock_t tStart = clock();
  for(unsigned int ev=0; ev<nEv; ev++) {
    if((ev%100)==0) std::cout << "\rGenerating cluster data... " << ev*100/nEv << "%";
    tree->GetEvent(ev);

    if(doConvertADC) {
      for(unsigned int i=0; i<nBins; i++) {
	amplitude[i] /= slopePerChannel[i];
      }
    }
    
    if(isSetTimeCuts && (time<timeCut_lower||time>timeCut_upper)) continue;

    if(doCommonModeCorrection) {
      float amplitude_part[nChannels];
      unsigned int k=0;
      for(unsigned int j=channel_startPos;j<=channel_endPos;j++) {
	amplitude_part[k] = amplitude[j-1]; //creating an appropriate array for TMath::Median(...)
	k++;
      }
      median = TMath::Median(nChannels,amplitude_part);
      for(unsigned int i=0;i<nBins;i++) amplitude[i] -= median;
    }

    CREvent *currentCREvent = new CREvent(nChannels, amplitude, nBins);

    for(unsigned int channel=channel_startPos; channel<=channel_endPos; channel++) {
      for(unsigned int clusterSize=1; clusterSize<=maxClusterSize; clusterSize++) {
	float amplitudeOfCluster = 0., noiseOfCluster = 0.;
	if(channel+(clusterSize-1) <= channel_endPos) { // prevent array overfloats
	  for(unsigned int nn=0; nn<clusterSize; nn++) {
	    amplitudeOfCluster += amplitude[(channel-1)+nn]; // amplitudeOfCluster = a_1 + a_2 + ...
	    noiseOfCluster += TMath::Power(hist_noise->GetBinContent(channel+nn),2); // noiseOfCluster^2 = delta_a_1^2 + delta_a_2^2 + ...
	  }
	  noiseOfCluster = TMath::Sqrt(noiseOfCluster);
	  CREntry currentCREntry(ev, channel, clusterSize, amplitudeOfCluster, amplitudeOfCluster/noiseOfCluster);
	  currentCREvent->addCREntry(currentCREntry);
	}
      }
    }
    CREventCollection.push_back(*currentCREvent);

    if(doSetSNRThreshold) {
      if(currentCREvent->getMostSignCREntry().getSignificanceOfCluster() >= SNRThreshold) { // set SNR threshold
	hist_clusterSizeDistr->Fill(currentCREvent->getMostSignCREntry().getClusterSize());
	hist_clusterSizeDistr_background->Fill(currentCREvent->getLeastSignCREntry().getClusterSize());
	hist_spectrum->Fill(currentCREvent->getMostSignCREntry().getAmplitudeOfCluster());
	hist_background->Fill(-currentCREvent->getLeastSignCREntry().getAmplitudeOfCluster());
      }
    }else {
      hist_clusterSizeDistr->Fill(currentCREvent->getMostSignCREntry().getClusterSize());
      hist_clusterSizeDistr_background->Fill(currentCREvent->getLeastSignCREntry().getClusterSize());
      hist_spectrum->Fill(currentCREvent->getMostSignCREntry().getAmplitudeOfCluster());
      hist_background->Fill(-currentCREvent->getLeastSignCREntry().getAmplitudeOfCluster());
    }
    delete currentCREvent;
  }
  std::cout << "\rGenerating cluster data... DONE "
	    << "(time taken: " << (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds)" 
	    << std::endl << std::endl;


  ///////////////////////////////////////////////////////////////////
  
  
  unsigned int integral_upperBoundary = scaleCutValue-min;
  float scaleVal = (hist_spectrum->Integral(1,integral_upperBoundary))/(hist_background->Integral(1,integral_upperBoundary));
  if(doSubtractBackground) CRprintInfo(scaleVal, "Calculated scale value for the background distribution");
  hist_background->Scale(scaleVal);
  hist_clusterSizeDistr_background->Scale(scaleVal);
  hist_background->Sumw2(); // necessary after histogram scaling
  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  hist_spectrum->Draw("e1 x0");
  hist_background->SetLineColor(kRed);
  hist_background->Draw("e1 x0 hist same");

  if(doSubtractBackground) {
    hist_spectrum->Add(hist_background,-1);
    hist_clusterSizeDistr->Add(hist_clusterSizeDistr_background,-1);
    for(unsigned int i=1;i<=hist_spectrum->GetNbinsX();i++) {
      hist_spectrum->SetBinError(i,TMath::Sqrt(TMath::Power(hist_spectrum->GetBinError(i),2) + TMath::Power(hist_background->GetBinError(i),2)));
    }
  }

  if(doLangausFit) {
    std::cout << "Performing fit with Landau-Gauss convolution function..." << std::endl;
    tStart = clock();
    TF1 *fitsnr = new TF1();

    // Setting fit range and start values
    if(!isLangausFitRange_start) {
      fr[0] = 0.4*hist_spectrum->GetMean();
    }
    if(!isLangausFitRange_end) {
      fr[1] = 3.0*hist_spectrum->GetMean();
    }

    // "Width","MP","Area","GSigma":
    if(!isLangausRanges) {
      pllo[0]=0.25; pllo[1]=5.0; pllo[2]=1.0; pllo[3]=0.4;
      plhi[0]=80.0; plhi[1]=250.0; plhi[2]=1e6; plhi[3]=50.0;
    }
    if(!isLangausStartValues) {
      sv[0]=2.; sv[1]=20.; sv[2]=5e4; sv[3]=5.;
      // sv[0]=1.8; sv[1]=20.0; sv[2]=50000.0; sv[3]=3.0;
    }
    Double_t chisqr;
    Int_t    ndf;
    fitsnr = langaufit(hist_spectrum,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);

    Double_t SNRPeak, SNRFWHM;
    langaupro(fp,SNRPeak,SNRFWHM);

    std::cout << "   SNRPeak = " << SNRPeak << std::endl;
    std::cout << "   SNRFWHM = " << SNRFWHM << std::endl;
    std::cout << "Performing fit with Landau-Gauss convolution function... DONE " 
	      << "(time taken: " << (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds)" << std::endl;

    // Global style settings
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(111);
    gStyle->SetLabelSize(0.03,"x");
    gStyle->SetLabelSize(0.03,"y");

    if(doConvertADC) {
      hist_spectrum->SetXTitle("Amplitude / (ke)");
    }else {
      hist_spectrum->SetXTitle("Amplitude / (ADC)");
    }
    hist_spectrum->SetYTitle("Counts / (#Delta amplitude)");
    hist_spectrum->GetYaxis()->SetTitleOffset(1.2);
    hist_spectrum->Draw("e1 x0");
    fitsnr->Draw("lsame");
  }

  c->SaveAs("outputfiles/temp_signal.root");
  c->SaveAs("outputfiles/temp_signal.pdf");
  
  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->SetGridy();
  hist_clusterSizeDistr->Scale(100/((float)hist_clusterSizeDistr->GetEntries())); // normalization and representation as percentage
  hist_clusterSizeDistr->SetXTitle("Multiplicity");
  hist_clusterSizeDistr->SetYTitle("Percentage");
  hist_clusterSizeDistr->GetYaxis()->SetRangeUser(0,100);
  hist_clusterSizeDistr->Draw();
  hist_clusterSizeDistr->SaveAs("outputfiles/temp_clustDistr.root");
  c2->SaveAs("outputfiles/temp_clustDistr.pdf");


  // generate further plots related to the charge distribution of two-strip-cluster events:
  std::vector<CREvent> CREventCollection_2SC;
  for(unsigned int ev=0; ev<CREventCollection.size(); ev++) {
    CREventCollection_2SC.push_back(CREventCollection.at(ev).getCREvent_forCS(2));
  }
  int eta_bins = 100, eta_min = 0, eta_max = 1;
  TH1F *hist_eta = new TH1F("hist_eta","",eta_bins,eta_min,eta_max);
  TH1F *hist_eta_background = new TH1F("hist_eta_background","",eta_bins,eta_min,eta_max);
  TH1F *hist_angular = new TH1F("hist_angular","",90,0,TMath::Pi()/2.);
  int crosstalk_bins = 130, 
    crosstalk_min = -10, 
    crosstalk_max = 120;
  if(doConvertADC) {
    crosstalk_bins = 130;
    crosstalk_min = -5;
    crosstalk_max = 55;
  }
  TH2F *hist_crosstalk = new TH2F("hist_crosstalk","",crosstalk_bins,crosstalk_min,crosstalk_max,crosstalk_bins,crosstalk_min,crosstalk_max);
  TH2F *hist_crosstalk_background = new TH2F("hist_crosstalk_background","",crosstalk_bins,crosstalk_min,crosstalk_max,crosstalk_bins,crosstalk_min,crosstalk_max);
  for(unsigned int ev=0; ev<CREventCollection_2SC.size(); ev++) {
    CREntry mostSignCREntry = CREventCollection_2SC.at(ev).getMostSignCREntry();
    CREntry leastSignCREntry = CREventCollection_2SC.at(ev).getLeastSignCREntry();
    float amp_leftChannel = CREventCollection_2SC.at(ev).getAmplitude_at(mostSignCREntry.getStartChannel()-1);
    float amp_rightChannel = CREventCollection_2SC.at(ev).getAmplitude_at(mostSignCREntry.getStartChannel());
    float eta = amp_rightChannel/(amp_leftChannel+amp_rightChannel);
    if(doSetSNRThreshold && CREventCollection_2SC.at(ev).getMostSignCREntry().getSignificanceOfCluster() >= SNRThreshold) { //set SNR threshold
      hist_crosstalk->Fill(amp_rightChannel, amp_leftChannel);
      hist_angular->Fill(TMath::ATan(amp_leftChannel/amp_rightChannel));
      hist_eta->Fill(eta);
    }else if(!doSetSNRThreshold && amp_leftChannel*amp_leftChannel+amp_rightChannel*amp_rightChannel >= 400) { //dirty approach: just cut away background part within a circle of certain radius
      hist_crosstalk->Fill(amp_rightChannel, amp_leftChannel);
      hist_angular->Fill(TMath::ATan(amp_leftChannel/amp_rightChannel));
      hist_eta->Fill(eta);
    }
    if(doSubtractBackground) {
      float amp_leftChannel_background = CREventCollection_2SC.at(ev).getAmplitude_at(leastSignCREntry.getStartChannel()-1);
      float amp_rightChannel_background = CREventCollection_2SC.at(ev).getAmplitude_at(leastSignCREntry.getStartChannel());
      hist_crosstalk_background->Fill(-amp_rightChannel_background, -amp_leftChannel_background);
      float eta_background = amp_rightChannel_background/(amp_leftChannel_background+amp_rightChannel_background);
      hist_eta_background->Fill(eta_background);   
    }
  }

  if(doSubtractBackground) {
    hist_crosstalk_background->Scale(scaleVal);
    hist_crosstalk_background->Sumw2();
    //hist_crosstalk->Add(hist_crosstalk_background,-1);

    hist_eta_background->Scale(scaleVal);
    hist_eta_background->Sumw2();
    //hist_eta->Add(hist_eta_background,-1);
  }
  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->SetGrid();
  hist_eta->SetXTitle("#eta_{ (m=2)}");
  hist_eta->SetYTitle("Counts");
  hist_eta->GetYaxis()->SetTitleOffset(1.2);
  hist_eta->SetStats(kFALSE);
  hist_eta->Draw("e");
  hist_eta->SaveAs("outputfiles/temp_hist_eta.root");
  c3->SaveAs("outputfiles/temp_eta.pdf");
  c3->SaveAs("outputfiles/temp_eta.root");

  TCanvas *c3_2 = new TCanvas("c3_2","",800,600);
  c3_2->SetGrid();
  hist_eta_background->SetXTitle("#eta_{ (m=2)}");
  hist_eta_background->SetYTitle("Counts");
  hist_eta_background->GetYaxis()->SetTitleOffset(1.2);
  hist_eta_background->SetStats(kFALSE);
  hist_eta_background->Draw("e");
  hist_eta_background->SaveAs("outputfiles/temp_hist_eta_background.root");
  c3_2->SaveAs("outputfiles/temp_eta_background.pdf");
  c3_2->SaveAs("outputfiles/temp_eta_background.root");

  TCanvas *c4 = new TCanvas("c4","",800,600);
  c4->SetGrid();
  if(doConvertADC) {
    hist_crosstalk->SetXTitle("a_{R} / (ke)");
    hist_crosstalk->SetYTitle("a_{L} / (ke)");
  }else {
    hist_crosstalk->SetXTitle("a_{R} / (ADC)");
    hist_crosstalk->SetYTitle("a_{L} / (ADC)");
  }
  hist_crosstalk->SetStats(kFALSE);
  hist_crosstalk->Draw("colz");
  c4->SaveAs("outputfiles/temp_crosstalk.pdf");
  c4->SaveAs("outputfiles/temp_crosstalk.root");
  hist_crosstalk->SaveAs("outputfiles/temp_hist_crosstalk.root");

  
  TCanvas *c5 = new TCanvas("c5","",800,600);
  c5->SetGrid();
  if(doConvertADC) {
    hist_crosstalk_background->SetXTitle("a_{R} / (ke)");
    hist_crosstalk_background->SetYTitle("a_{L} / (ke)");
  }else {
    hist_crosstalk_background->SetXTitle("a_{R} / (ADC)");
    hist_crosstalk_background->SetYTitle("a_{L} / (ADC)");
  }
  hist_crosstalk_background->SetStats(kFALSE);
  hist_crosstalk_background->Draw("colz");
  c5->SaveAs("outputfiles/temp_crosstalk_background.pdf");
  c5->SaveAs("outputfiles/temp_crosstalk_background.root");
  hist_crosstalk_background->SaveAs("outputfiles/temp_hist_crosstalk_background.root");
  
  TCanvas *c6 = new TCanvas("c6","",800,600);
  c6->SetGrid();
  hist_angular->SetXTitle("Angle");
  hist_angular->SetYTitle("Counts");
  hist_angular->SetStats(kFALSE);
  TF1 *gausfit_lowertheta = new TF1("gausfit_lowertheta","[0]*exp(-.5*((x-[1])/[2])**2)",.04,.348);
  gausfit_lowertheta->SetParameters(100.,.2,.1);
  hist_angular->Fit("gausfit_lowertheta","","",.04,.348);
  TF1 *gausfit_uppertheta = new TF1("gausfit_uppertheta","[0]*exp(-.5*((x-[1])/[2])**2)",1.224,1.53);
  gausfit_uppertheta->SetParameters(180.,1.4,.1);
  hist_angular->Fit("gausfit_uppertheta","+","",1.224,1.53);
  hist_angular->Draw("e");
  TLatex latex;
  char mu_1[45], mu_2[45], sigma_1[45], sigma_2[45];
  sprintf(mu_1,"#scale[.75]{#color[2]{#mu_{1} = %.3f}}",gausfit_lowertheta->GetParameter(1));
  sprintf(mu_2,"#scale[.75]{#color[2]{#mu_{2} = %.3f}}",gausfit_uppertheta->GetParameter(1));
  sprintf(sigma_1,"#scale[.75]{#color[2]{#sigma_{1} = %.3f}}",gausfit_lowertheta->GetParameter(2));
  sprintf(sigma_2,"#scale[.75]{#color[2]{#sigma_{2} = %.3f}}",gausfit_uppertheta->GetParameter(2));
  latex.DrawLatex(.1,205,mu_1);
  latex.DrawLatex(.1,185,sigma_1);
  latex.DrawLatex(.95,205,mu_2);
  latex.DrawLatex(.95,185,sigma_2);
  c6->SaveAs("outputfiles/temp_angular.pdf");
  c6->SaveAs("outputfiles/temp_angular.root");
  hist_angular->SaveAs("outputfiles/temp_hist_angular.root");

  return 0;
}
