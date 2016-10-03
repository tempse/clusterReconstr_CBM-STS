/** \file main.cpp
 * \brief This is the main routine of the cluster reconstruction algorithm. It reads a macro file, imports its parameters, performs the cluster reconstruction algorithm and outputs the results.
 *
 * \author Sebastian Templ <sebastian.templ@gmail.com>
 * \version 1.0
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
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TError.h>
#include <TF1.h>
#include <TColor.h>

#include "CREntry.h"
#include "CREvent.h"

//#include "read_treeCR.h"

#include "langaufit.h"
#include "langaufun.h"
#include "langaupro.h"


int main(int argc, char** argv) {
  TString argument;
  
  TString macroFileName;
  TString filename, calibrationFileName;

  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4]; // parameters for Landau*Gauss fit

  // default values:
  TString outFileName = "output_";
  size_t nBins = 256;
  unsigned int channel_startPos = 0;
  unsigned int channel_endPos   = 256;
  unsigned int timeCut_lower = 0;
  unsigned int timeCut_upper = 25;
  unsigned int maxClusterSize = 1;
  int scaleCutValue = 10;
  bool isMacroFile = false;
  bool doConvertADC = false;
  bool doSubtractBackground = true;
  bool doCommonModeCorrection = true;
  bool doLangausFit = false;
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
    std::cout << std::endl;
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
	  std::cout << "First channel number set to " << channel_startPos << "..." << std::endl;
	}else {
	  std::cout << "   ERROR: 'CHANNEL_START = " << value << "' is an invalid input." << std::endl;
	  return 3;
	}
      }
      if(identifier == "CHANNEL_END") {
	if(value.IsDigit()) {
	  channel_endPos = value.Atoi();
	  std::cout << "Last channel number set to " << channel_endPos << "..." << std::endl;
	}else {
	  std::cout << "   ERROR: 'CHANNEL_END = " << value << "' is an invalid input." << std::endl;
	  return 3;
	}
      }
      if(identifier == "NBINS") {
	if(value.IsDigit()) {
	  nBins = value.Atoi();
	  std::cout << "Total number of bins set to " << nBins << "..." << std::endl;
	}else {
	  std::cout << "   WARNING: 'NBINS = " << value << "' is an invalid input. The default value will be used." << std::endl;
	}
      }
      if(identifier == "TIMECUT_LOWER") {
	if(value.IsDigit()) {
	  isSetTimeCuts = true;
	  isSetTimeCut_lower = true;
	  timeCut_lower = value.Atoi();
	  std::cout << "Lower time cut value set to " << timeCut_lower << "..." << std::endl;
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
	  std::cout << "Upper time cut value set to " << timeCut_upper << "..." << std::endl;
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
	  std::cout << "Maximum cluster size set to " << maxClusterSize << "..." << std::endl;
	}else {
	  std::cout << "   ERROR: 'MAXCLUSTERSIZE = " << value << "' is an invalid input." << std::endl;
	  return 3;
	}
      }
      if(identifier == "COMMONMODECORRECTION") {
	value.ToLower();
	if(value == "yes") {
	  doCommonModeCorrection = true;
	  std::cout << "A common-mode correction will be performed..." << std::endl;
	}else if(value == "no") {
	  doCommonModeCorrection = false;
	  std::cout << "A common-mode correction won't be performed..." << std::endl;
	}else {
	  std::cout << "   WARNING: 'COMMONMODECORRECTION = " << value << "' is an invalid input. The default value will be used." << std::endl;
	}
      }
      if(identifier == "SUBTRACTBACKGROUND") {
	value.ToLower();
	if(value == "yes") {
	  doSubtractBackground = true;
	  std::cout << "The estimated background will be subtracted..." << std::endl;
	}else if(value == "no") {
	  doSubtractBackground = false;
	  std::cout << "The estimated background won't be subtracted..." << std::endl;
	}else {
	  std::cout << "   WARNING: 'SUBTRACTBACKGROUND = " << value << "' is an invalid input. The default value will be used." << std::endl;
	}
      }
      if(identifier == "SCALECUTVALUE") {
	if(value.IsDigit()) {
	  scaleCutValue = value.Atoi();
	  std::cout << "Scale cut value set to " << scaleCutValue << "..." << std::endl;
	}else {
	  std::cout << "   WARNING: 'SCALECUTVALUE = " << value << "' is an invalid input. The default value will be used." << std::endl;
	}
      }
      if(identifier == "LANGAUSFIT") {
	value.ToLower();
	if(value == "yes") {
	  doLangausFit = true;
	  std::cout << "A fit with a Landau-Gauss convolution function will be performed..." << std::endl;
	}else if(value == "no") {
	  doLangausFit = false;
	  std::cout << "A fit with a Landau-Gauss convolution function won't be performed..." << std::endl;
	}else {
	  std::cout << "   WARNING: 'LANGAUSFIT = " << value << "' is an invalid input. The default value will be used." << std::endl;
	}
      }
      if(identifier.BeginsWith("LANGAUS")) {
	std::string value_str;
	value_str = value;
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
      std::cout << "Parameters for the Landau*Gauss fit set..." << std::endl;
    }else if(!isLangausRanges) {
      std::cout << "   WARNING: Not all ranges for the Landau*Gauss fit were set properly. Default values will be used." << std::endl;
    }else if(!isLangausStartValues) {
      std::cout << "   WARNING: Not all start values for the Landau*Gauss fit were set properly. Default values will be used." << std::endl;
    }
  }
  std::cout << std::endl << "All parameters set..." << std::endl;


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

  // TFile *fileOutTree = new TFile("treeCR.root", "RECREATE");
  // TTree* treeCR = new TTree("treeCR", "Tree containing cluster reconstruction data");
  // int currentEventID[nChannels*maxClusterSize], currentStartChannel[nChannels*maxClusterSize], currentClusterSize[nChannels*maxClusterSize]; // reserve the maximum amount of space
  // float currentAmplitudeOfCluster[nChannels*maxClusterSize], currentSignificanceOfCluster[nChannels*maxClusterSize]; // reserve the maximum amount of space
  // treeCR->Branch("eventID",currentEventID);
  // treeCR->Branch("startChannel",currentStartChannel);
  // treeCR->Branch("clusterSize",currentClusterSize);
  // treeCR->Branch("amplitudeOfCluster",currentAmplitudeOfCluster);
  // treeCR->Branch("significanceOfCluster",currentSignificanceOfCluster);

  TH1F *hist_clusterSizeDistr = new TH1F("hist_clusterSizeDistr","",maxClusterSize,1,maxClusterSize+1);

  std::vector<CREvent> CREventCollection;


 ///////////////////////////////////////////////////////////////////

  unsigned int nEv = tree->GetEntries();
  std::cout << "\rMaking preparations... DONE" << std::endl;
  std::cout << "Starting cluster reconstruction algorithm..." << std::endl << std::endl;
  clock_t tStart = clock();
  for(unsigned int ev=0; ev<nEv; ev++) {
    if((ev%100)==0) std::cout << "\rGenerating cluster data... " << ev*100/nEv << "%";
    tree->GetEvent(ev);

    if(doConvertADC) {
      for(unsigned int channel=channel_startPos; channel<=channel_endPos; channel++) {
	treeCalibr->GetEvent(channel-1);
	amplitude[channel-1] /= slope_calibr;
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

    CREvent *currentCREvent = new CREvent(nChannels, amplitude);

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


    // for(unsigned int clusterSize=1; clusterSize<=maxClusterSize; clusterSize++) {
    //   TString outFileName = "treeCR_";
    //   outFileName += clusterSize;
    //   outFileName += "StripClusters.root";
    //   TFile *outFile = new TFile(outFileName, "UPDATE");

    //   TString treeName = clusterSize;
    //   treeName += "StripClusters";
    //   TTree *treeCR = new TTree(treeName, "");
    //   treeCR->Branch("startChannel",currentStartChannel);
    //   treeCR->Branch("clusterSize",currentClusterSize);
    //   treeCR->Branch("amplitudeOfCluster",currentAmplitudeOfCluster);
    //   treeCR->Branch("significanceOfCluster",currentSignificanceOfCluster);

    //   for(unsigned int i=0; i<currentCREvent->getNumberOfCREntries(); i++) {
    // 	if(currentCREvent->getCREntry(i).getClusterSize() == clusterSize) {
    // 	  currentStartChannel[i] = currentCREvent->getCREntry(i).getStartChannel();
    // 	  currentClusterSize[i] = currentCREvent->getCREntry(i).getClusterSize();
    // 	  currentAmplitudeOfCluster[i] = currentCREvent->getCREntry(i).getAmplitudeOfCluster();
    // 	  currentSignificanceOfCluster[i] = currentCREvent->getCREntry(i).getSignificanceOfCluster();
    // 	}
    //   }
    //   treeCR->Fill();
    //   treeCR->Write();
    //   outFile->Close();
    // }

    /*
      for(unsigned int i=0; i<currentCREvent->getNumberOfCREntries(); i++) {
        currentEventID[i] = currentCREvent->getCREntry(i).getEventID();
	currentStartChannel[i] = currentCREvent->getCREntry(i).getStartChannel();
	currentClusterSize[i] = currentCREvent->getCREntry(i).getClusterSize();
	currentAmplitudeOfCluster[i] = currentCREvent->getCREntry(i).getAmplitudeOfCluster();
	currentSignificanceOfCluster[i] = currentCREvent->getCREntry(i).getSignificanceOfCluster();
	}
      for(unsigned int i=currentCREvent->getNumberOfCREntries(); i<nChannels*maxClusterSize; i++) {
        currentEventID[i] = markerVal;
	currentStartChannel[i] = markerVal;
	currentClusterSize[i] = markerVal;
	currentAmplitudeOfCluster[i] = markerVal;
	currentSignificanceOfCluster[i] = markerVal;
      }
      treeCR->Fill();
    */
    // CREventCollection.push_back(*currentCREvent); //testing

    //if(currentCREvent->getMostSignCREntry().getSignificanceOfCluster() > 3) 
      {
	hist_clusterSizeDistr->Fill(currentCREvent->getMostSignCREntry().getClusterSize());
	hist_spectrum->Fill(currentCREvent->getMostSignCREntry().getAmplitudeOfCluster());
	hist_background->Fill(-currentCREvent->getLeastSignCREntry().getAmplitudeOfCluster());
	// hist_spectrum->Fill(currentCREvent->getMostSignAmp_forCS(3));
	// hist_background->Fill(-currentCREvent->getLeastSignAmp_forCS(3));
      }
    delete currentCREvent;
  }
  std::cout << "\rGenerating cluster data... DONE "
	    << "(time taken: " << (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds)" 
	    << std::endl << std::endl;
  
  // std::cout << "Writing cluster data...";
  // treeCR->Write();
  // fileOutTree->Close();
  // std::cout << " DONE" << std::endl;

  // read_treeCR(nChannels, maxClusterSize, markerVal);

  ///////////////////////////////////////////////////////////////////
  
  
  unsigned int integral_upperBoundary = scaleCutValue-min;
  float scaleVal = (hist_spectrum->Integral(1,integral_upperBoundary))/(hist_background->Integral(1,integral_upperBoundary));
  hist_background->Scale(scaleVal);
  hist_background->Sumw2(); // necessary after histogram scaling
  TCanvas *c = new TCanvas();
  c->SetGrid();
  hist_spectrum->Draw("e1 x0");
  hist_background->SetLineColor(kRed);
  hist_background->Draw("e1 x0 hist same");

  if(doSubtractBackground) {
    hist_spectrum->Add(hist_background,-1);
    for(unsigned int i=1;i<=hist_spectrum->GetNbinsX();i++) {
      hist_spectrum->SetBinError(i,TMath::Sqrt(TMath::Power(hist_spectrum->GetBinError(i),2) + TMath::Power(hist_background->GetBinError(i),2)));
    }
  }

  if(doLangausFit) {
    std::cout << "Performing fit with Landau-Gauss convolution function..." << std::endl;
    tStart = clock();
    TF1 *fitsnr = new TF1();

    // Setting fit range and start values
    Double_t fr[2];
    fr[0]=0.25*hist_spectrum->GetMean();
    fr[1]=3.0*hist_spectrum->GetMean();

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
    hist_spectrum->Draw();
    fitsnr->Draw("lsame");
  }

  c->SaveAs("outputfiles/temp.root");
  c->SaveAs("outputfiles/temp.pdf");
  
  TCanvas *c2 = new TCanvas();
  c2->SetGridy();
  hist_clusterSizeDistr->Scale(100/((float)hist_clusterSizeDistr->GetEntries())); // normalization and representation as percentage
  hist_clusterSizeDistr->SetXTitle("Multiplicity");
  hist_clusterSizeDistr->SetYTitle("Percentage");
  hist_clusterSizeDistr->Draw();
  hist_clusterSizeDistr->SaveAs("outputfiles/temp_clustDistr.root");
  c2->SaveAs("outputfiles/temp_clustDistr.pdf");


  return 0;
}
