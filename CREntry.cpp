#include "CREntry.h"


CREntry::CREntry(const CREntry &newCREntry) {
  eventID = newCREntry.getEventID();
  startChannel = newCREntry.getStartChannel();
  clusterSize = newCREntry.getClusterSize();
  amplitudeOfCluster = newCREntry.getAmplitudeOfCluster();
  significanceOfCluster = newCREntry.getSignificanceOfCluster();
}

CREntry::CREntry(unsigned int newEventID, unsigned int newStartChannel, unsigned int newClusterSize, float newAmplitudeOfCluster, float newSignificanceOfCluster) {
  eventID = newEventID;
  startChannel = newStartChannel;
  clusterSize = newClusterSize;
  amplitudeOfCluster = newAmplitudeOfCluster;
  significanceOfCluster = newSignificanceOfCluster;
}

CREntry CREntry::operator=(const CREntry &newCREntry) {
  eventID = newCREntry.getEventID();
  startChannel = newCREntry.getStartChannel();
  clusterSize = newCREntry.getClusterSize();
  amplitudeOfCluster = newCREntry.getAmplitudeOfCluster();
  significanceOfCluster = newCREntry.getSignificanceOfCluster();
}

void CREntry::setEventID(unsigned int newEventID) {
  eventID = newEventID;
}

void CREntry::setStartChannel(unsigned int newStartChannel) {
  startChannel = newStartChannel;
}
void CREntry::setClusterSize(unsigned int newClusterSize) {
  clusterSize = newClusterSize;
}
void CREntry::setAmplitudeOfCluster(float newAmplitudeOfCluster) {
  amplitudeOfCluster = newAmplitudeOfCluster;
}
void CREntry::setSignificanceOfCluster(float newSignificanceOfCluster) {
  significanceOfCluster = newSignificanceOfCluster;
}

unsigned int CREntry::getEventID() const {
  return eventID;
}

unsigned int CREntry::getStartChannel() const {
  return startChannel;
}
unsigned int CREntry::getClusterSize() const {
  return clusterSize;
}
float CREntry::getAmplitudeOfCluster() const {
  return amplitudeOfCluster;
}
float CREntry::getSignificanceOfCluster() const {
  return significanceOfCluster;
}
