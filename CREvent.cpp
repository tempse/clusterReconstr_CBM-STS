#include "CREvent.h"

CREvent::CREvent(unsigned int newNChannels) {
  nChannels = newNChannels;
}

CREvent::CREvent(unsigned int newNChannels, unsigned int newNBins, float *newAmplitudes) {
  nChannels = newNChannels;
  for(unsigned int i=0; i<newNBins; i++) {
    amplitudes.push_back(newAmplitudes[i]);
  }
}

void CREvent::addCREntry(CREntry newCREntry) {
  CREntries.push_back(newCREntry);
}

unsigned int CREvent::getNumberOfCREntries() const {
  return CREntries.size();
}

CREntry CREvent::getMostSignCREntry() const {
  float sign_max = -1e4;
  CREntry mostSignCREntry;
  for(unsigned int i=0; i<CREntries.size(); i++) {
    if(CREntries.at(i).getSignificanceOfCluster() > sign_max) {
      sign_max = CREntries.at(i).getSignificanceOfCluster();
      mostSignCREntry = CREntries.at(i);
    }
  }
  return mostSignCREntry;
}

CREntry CREvent::getLeastSignCREntry() const {
  float sign_min = +1e4;
  CREntry leastSignCREntry;
  for(unsigned int i=0; i<CREntries.size(); i++) {
    if(CREntries.at(i).getSignificanceOfCluster() < sign_min) {
      sign_min = CREntries.at(i).getSignificanceOfCluster();
      leastSignCREntry = CREntries.at(i);
    }
  }
  return leastSignCREntry;
}

float CREvent::getMostSignAmp_forCS(unsigned int clusterSize) const {
  float sign_max = -1e4;
  float mostSignificantAmplitude;
  for(unsigned int i=0; i<CREntries.size(); i++) {
    if((CREntries.at(i).getClusterSize() == clusterSize) && (CREntries.at(i).getSignificanceOfCluster() > sign_max)) {
      sign_max = CREntries.at(i).getSignificanceOfCluster();
      mostSignificantAmplitude = CREntries.at(i).getAmplitudeOfCluster();
    }
  }
  return mostSignificantAmplitude;
}

float CREvent::getLeastSignAmp_forCS(unsigned int clusterSize) const {
  float sign_min = +1e4;
  float leastSignificantAmplitude;
  for(unsigned int i=0; i<CREntries.size(); i++) {
    if((CREntries.at(i).getClusterSize() == clusterSize) && (CREntries.at(i).getSignificanceOfCluster() < sign_min)) {
      sign_min = CREntries.at(i).getSignificanceOfCluster();
      leastSignificantAmplitude = CREntries.at(i).getAmplitudeOfCluster();
    }
  }
  return leastSignificantAmplitude;
}

CREntry CREvent::getCREntry(unsigned int i) const {
  return CREntries.at(i);
}

/*
void CREvent::writeCREvent(TString fileName) {
  TFile *outFile = new TFile(fileName, "RECREATE");
  //unsigned int CREntries_size = (int)CREntries.size();
  unsigned int currentStartChannel[CREntries.size()], currentClusterSize[CREntries.size()];
  float currentAmplitudeOfCluster[CREntries.size()], currentSignificanceOfCluster[CREntries.size()];
  TTree *treeCR = new TTree("treeCR","treeCR");
  treeCR->Branch("startChannel",currentStartChannel);
  treeCR->Branch("clusterSize",currentClusterSize);
  treeCR->Branch("amplitudeOfCluster",currentAmplitudeOfCluster);
  treeCR->Branch("significanceOfCluster",currentSignificanceOfCluster);
  for(unsigned int i=0; i<CREntries.size(); i++) {
    currentStartChannel[i] = CREntries.at(i).getStartChannel();
    currentClusterSize[i] = CREntries.at(i).getClusterSize();
    currentAmplitudeOfCluster[i] = CREntries.at(i).getAmplitudeOfCluster();
    currentSignificanceOfCluster[i] = CREntries.at(i).getSignificanceOfCluster();
  }
  treeCR->Fill();
  // treeCR->Write();
}
*/
void CREvent::printAmplitudes() {
  for(std::vector<float>::iterator it=amplitudes.begin(); it!=amplitudes.end(); it++) {
    std::cout << *it << std::endl;
  }
}
