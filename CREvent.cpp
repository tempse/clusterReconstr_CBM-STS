#include "CREvent.h"

CREvent::CREvent(unsigned int newNChannels, float *newAmplitudes, unsigned int newNBins) {
  nChannels = newNChannels;
  for(unsigned int i=0; i<newNBins; i++) {
    amplitudes.push_back(newAmplitudes[i]);
  }
}

CREvent::CREvent(unsigned int newNChannels, std::vector<float> newAmplitudes) {
  nChannels = newNChannels;
  for(unsigned int i=0; i<newAmplitudes.size(); i++) {
    amplitudes.push_back(newAmplitudes.at(i));
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

CREvent CREvent::getCREvent_forMaxCS(unsigned int maxClusterSize) {
  CREvent *temp_CREvent = new CREvent(nChannels, amplitudes);
  for(unsigned int i=0; i<CREntries.size(); i++) {
    if(CREntries.at(i).getClusterSize() <= maxClusterSize) {
      temp_CREvent->addCREntry(CREntries.at(i));
    }
  }
  return *temp_CREvent;
}

CREvent CREvent::getCREvent_forCS(unsigned int clusterSize) {
  CREvent *temp_CREvent = new CREvent(nChannels, amplitudes);
  for(unsigned int i=0; i<CREntries.size(); i++) {
    if(CREntries.at(i).getClusterSize() == clusterSize) {
      temp_CREvent->addCREntry(CREntries.at(i));
    }
  }
  return *temp_CREvent;
}

unsigned int CREvent::getAmplitudes_size() {
  return (unsigned int)amplitudes.size();
}

float CREvent::getAmplitude_at(unsigned int i) {
  return amplitudes.at(i);
}

void CREvent::printAmplitudes() {
  for(std::vector<float>::iterator it=amplitudes.begin(); it!=amplitudes.end(); it++) {
    std::cout << *it << std::endl;
  }
}
