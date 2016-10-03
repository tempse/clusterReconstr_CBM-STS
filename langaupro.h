/** \file langaupro.h
 * \brief Function to find the maximum of a Landau-Gauss convolution and its FWHM value.
 */

//-----------------------------------------------------------------------
//
// Convoluted Landau and Gaussian Fitting Function
//         (using ROOT's Landau and Gauss functions)
//
//  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
//  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
//   Markus Friedl (Markus.Friedl@cern.ch)
//
//-----------------------------------------------------------------------

#ifndef LANGAUPRO_H
#define LANGAUPRO_H

#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

#include "langaufun.h"

Int_t langaupro(Double_t *, Double_t &, Double_t &);


#endif
