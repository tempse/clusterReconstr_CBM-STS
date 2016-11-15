/** \file CRprintInfo.h
 * \brief Function for convenient and formatted terminal output of parameters and variables.
 *
 * \author Sebastian Templ <sebastian.templ@gmail.com>
 * \date 2016
 */


#ifndef CRPRINTINFO_H
#define CRPRINTINFO_H

#include <iostream>
#include <boost/format.hpp>

#include <TString.h>


//! Function to output a string variable with descriptive text.
void CRprintInfo(TString val, TString text);

//! Function to output a boolean variable with descriptive text.
void CRprintInfo(bool val, TString text);

//! Function to output a float/double variable with descriptive text.
void CRprintInfo(float val, TString text);

//! Function to output an integer variable with descriptive text.
void CRprintInfo(int val, TString text);

//! Function to output an unsigned integer variable with descriptive text.
void CRprintInfo(unsigned int val, TString text);

#endif
