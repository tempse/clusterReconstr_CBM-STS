/** \file langaufit.h
 * \brief Function to perform a fit with a Landau-Gauss convolution function.
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

#ifndef LANGAUFIT_H
#define LANGAUFIT_H

#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

#include "langaufun.h"

TF1 *langaufit(TH1F *, Double_t *, Double_t *, Double_t *, Double_t *, Double_t *, Double_t *, Double_t *, Int_t *);

#endif
