/** \file CREntry.h
 * \brief Objects of this class hold all relevant information of a cluster.
 *
 * \author Sebastian Templ <sebastian.templ@gmail.com>
 * \date 2016
 */


#ifndef CRENTRY_H
#define CRENTRY_H

//! The class describing how to store all relevant information of a cluster.
/*! Objects of this class hold all relevant information of a cluster.
 */
class CREntry {
 public:
  //! Default constructor.
  CREntry() {};
  //! Copy constructor.
  CREntry(const CREntry&);
  //! Constructor which also initializes all relevant variables.
  CREntry(unsigned int newEventID, /*!< Event number. */
	  unsigned int newStartChannel, /*!< Leftmost channel of the cluster. */
	  unsigned int newClusterSize, /*!< Size of the cluster. */
	  float newAmplitudeOfCluster, /*!< Amplitude of the cluster (sum of individual amplitudes) */
	  float newSignificanceOfCluster /*!< Significance of the cluster (cluster amplitude divided by the square root of the sum of the squared noise values in the cluster) */
	  );
  //! Default destructor.
  ~CREntry() {};

  CREntry operator=(const CREntry&); //!< The '=' operator assigns all (private) variables separately.

  //! Function to set the variable eventID.
  void setEventID(unsigned int newEventID /*!< New event number. */);

  //! Function to set the variable startChannel.
  void setStartChannel(unsigned int newStartChannel /*!< New number for the leftmost channel of the cluster. */);

  //! Function to set the variable clusterSize.
  void setClusterSize(unsigned int newClusterSize /*!< New value for the cluster size. */);

  //! Function to set the variable amplitudeOfCluster.
  void setAmplitudeOfCluster(float newAmplitudeOfCluster /*!< New value of the cluster amplitude. */);

  //! Function to set the variable significanceOfCluster.
  void setSignificanceOfCluster(float newSignificanceOfCluster /*!< New value of the cluster significance. */);



  //! Function to return the value of the variable eventID.
  /*!
   * \return Value of the variable eventID as unsigned int.
   */
  unsigned int getEventID() const;

  //! Function to return the value of the variable startChannel.
  /*! 
   * \return Value of the variable startChannel as unsigned int.
   */
  unsigned int getStartChannel() const;

  //! Function to return the value of the variable clusterSize.
  /*!
   * \return Value of the variable clusterSize as unsigned int.
   */
  unsigned int getClusterSize() const;

  //! Function to return the value of the variable amplitudeOfCluster.
  /*!
   * \return Value of the variable amplitudeOfCluster as float.
   */
  float getAmplitudeOfCluster() const;

  //! Function to return the value of the variable significanceOfCluster.
  /*!
   * \return Value of the variable significanceOfCluster as float.
   */
  float getSignificanceOfCluster() const;

 private:
  unsigned int eventID; //!< The event number.
  unsigned int startChannel; //!< The leftmost channel of the cluster.
  unsigned int clusterSize; //!< The size of the cluster.
  float amplitudeOfCluster; //!< The amplitude of the cluster, i. e., the sum of the amplitudes of the individual channels.
  float significanceOfCluster; //!< The significance of the cluster, i. e., the cluster amplitude divided by the square root of the sum of the squared noise values in the cluster.
};

#endif
