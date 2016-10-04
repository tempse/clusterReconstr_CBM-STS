/** \file CREvent.h
 * \brief Objects of this class are supposed to hold all CREntry instances generated during a single event.
 *
 * \author Sebastian Templ <sebastian.templ@gmail.com>
 * \date 2016
 */


#ifndef CREVENT_H
#define CREVENT_H

#include <iostream>
#include <vector>
#include <algorithm>

/* #include <TString.h> */
/* #include <TFile.h> */
/* #include <TTree.h> */

#include "CREntry.h"

//! The class describing how to collect all CREntry objects of a single event and providing a variety of methods to process them.
/*! Objects of this class are supposed to hold all CREntry instances generated during a single event. A variety of methods to process them are provided.
 */
class CREvent : public CREntry {
 public:
  //! Default constructor.
  CREvent() {};
  //! Constructor which initializes the number-of-channels variable.
  CREvent(unsigned int newNChannels /*!< Total number of considered channels. */);
  //! Constructor which initializes the variables nChannels and amplitudes.
  CREvent(unsigned int newNChannels, /*!< Number of considered channels. */
	  unsigned int newNBins, /*!< Total number of bins, i. e., connected read-out channels. */
	  float *newAmplitudes /*!< The original amplitude array with newNBins entries. */
	  );
  //! Default destructor.
  ~CREvent() {};

  //! Add a CREntry object of the given event to the collection.
  /*! Add a CREntry object of the given event to the collection. These are stored in a vector.
   */
  void addCREntry(CREntry newCREntry /*!< Add this CREntry object to the collection of the given event. */);

  //! Function to return the number of the CREntry objects in the collection for a single event (i. e., the size of the vector in which they are stored).
  /*!
   * \return Number of the CREntry objects in the collection for a single event as unsigned int.
   */
  unsigned int getNumberOfCREntries() const;

  //! Function to return the most significant CREntry object.
  /*!
   * This function searches for the CREntry object with the biggest signifcance value and returns it.
   * \return CREntry object having the biggest significance value.
   */
  CREntry getMostSignCREntry() const;

  //! Function to return the least significant CREntry object.
  /*!
   * This function searches for the CREntry object with the smallest signifcance value and returns it.
   * \return CREntry object having the smallest significance value.
   */
  CREntry getLeastSignCREntry() const;

  //! Function to return the amplitude value of the cluster having the biggest significance value and a certain size (i. e., multiplicity).
  /*!
   * This function consideres only clusters of a given size, searches for the most significant cluster of said multiplicity and returns its amplitude.
   * \return Amplitude of the most significant cluster of a given size.
   */
  float getMostSignAmp_forCS(unsigned int clusterSize /*!< Cluster size to be considered. */) const;

//! Function to return the amplitude value of the cluster having the smallest significance value and a certain size (i. e., multiplicity).
  /*!
   * This function consideres only clusters of a given size, searches for the least significant cluster of said multiplicity and returns its amplitude.
   * \return Amplitude of the least significant cluster of a given size.
   */
  float getLeastSignAmp_forCS(unsigned int clusterSize /*!< Cluster size to be considered. */) const;

  //! Function to return a given CREntry object of the collection.
  /*!
   * \return CREntry object of a given position in the collection of CREntry objects of a single event.
   */
  CREntry getCREntry(unsigned int i /*!< Position of the respective CREntry object in the collection (i. e., the entry number of the vector in which they are stored). */) const;

  /* void writeCREvent(TString); */
  void printAmplitudes();

 private:
  unsigned int nChannels; //!< Total number of considered channels.
  std::vector<float> amplitudes; //!< The original amplitude array.
  std::vector<CREntry> CREntries; //!< The collection of CREntry objects of a single event.
};

#endif
