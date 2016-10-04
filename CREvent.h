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

#include "CREntry.h"

//! The class describing how to collect all CREntry objects of a single event and providing a variety of methods to process them.
/*! Objects of this class are supposed to hold all CREntry instances generated during a single event. A variety of methods to process them are provided.
 */
class CREvent : public CREntry {
 public:
  //! Default constructor.
  CREvent() {};
  //! Constructor which initializes the variables nChannels and amplitudes (using an amplitude array).
  CREvent(unsigned int newNChannels, /*!< Number of considered channels. */
	  float *newAmplitudes, /*!< The original amplitude array with newNBins entries. */
	  unsigned int newNBins /*!< Total number of bins (i.e., size of the amplitude array). */
	  );
  //! Constructor which initializes the variables nChannels and amplitudes (using an amplitude vector).
  CREvent(unsigned int newNChannels, /*!< Number of considered channels. */
	  std::vector<float> newAmplitudes /*!< The original amplitude array in form of a vector. */
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

  //! Function to generate a subset of the CREntry collection (i. e., the CREvent object) of a single event with cluster sizes lower than or equal to a given value.
  /*!
   * The function loops over the whole collection of CREntry objects of a given event, selects the ones having cluster sizes lower than or equal to a given value and returns them as a new CREvent object.
   * \return A new CREvent object holding only the CREntry objects with cluster sizes lower than or equal to a given value.
   */
  CREvent getCREvent_forMaxCS(unsigned int maxClusterSize /*!< The maximum cluster size to be considered. All CREntry objects having this cluster size or a smaller one are selected by the function. */);

  //! Function to generate a subset of the CREntry collection (i. e., the CREvent object) of a signle event with given cluster sizes.
  /*!
   * The function loops over the whole collection of CREntry objects of a given event, selects the ones having a given cluster size and returns them as a new CREvent object.
   * \return A new CREvent object holding only the CREntry objects with a given cluster size.
   */
  CREvent getCREvent_forCS(unsigned int clusterSize /*!< The cluster size to be considered. All CREntry objects having this cluster size are selected by the function. */);

  //! Function to return the size of the amplitude vector stored in the current CREvent object.
  /*!
   * \return Size of the amplitude vector as an unsigned int.
   */
  unsigned int getAmplitudes_size();

  //! Function to return the value of a given entry of the original amplitude array
  /*!
   * \return Value of the orignial amplitude array at a given position.
   * \attention The class instance must hold the amplitude array. This can be, e. g., achieved by calling the respective constructor.
   */
  float getAmplitude_at(unsigned int i /*!< Position of the amplitude value to be returned. */);

  void printAmplitudes();

 private:
  unsigned int nChannels; //!< Total number of considered channels.
  std::vector<float> amplitudes; //!< The original amplitude array.
  std::vector<CREntry> CREntries; //!< The collection of CREntry objects of a single event.
};

#endif
