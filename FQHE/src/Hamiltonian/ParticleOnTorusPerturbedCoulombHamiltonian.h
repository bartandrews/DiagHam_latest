////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                             coulombian interaction                         //
//                                                                            //
//                        last modification : 18/07/2002                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef PARTICLEONTORUSPERTURBEDCOULOMBHAMILTONIAN_H
#define PARTICLEONTORUSPERTURBEDCOULOMBHAMILTONIAN_H

#include "config.h"
#include "HilbertSpace/ParticleOnTorus.h"
#include "Hamiltonian/ParticleOnTorusGenericHamiltonian.h"
#include "Polynomial/Polynomial.h"
#include "Polynomial/SpecialPolynomial.h"

#include <iostream>

#ifdef __GMP__
#include <gmp.h>
#endif

using std::ostream;


class MathematicaOutput;


class ParticleOnTorusPerturbedCoulombHamiltonian : public ParticleOnTorusGenericHamiltonian
{

 protected:

  // Coulomb prefactor
  double CoulombPrefactor;

  // exponential decay parameter of Yukawa potential squared
  double YukawaParameterSqr;

  // landau Level index
  int LandauLevel;

  // form factor of the interaction (a single Laguerre polynomial for the Landau levels of GaAs)
  Polynomial FormFactor;

  // Number of additional Pseudopotentials
  int NbrPseudopotentials;
  // pseudopotential coefficients
  double* Pseudopotentials;
  // Laguerre polynomial for the pseudopotentials
#ifdef __GMP__
  // LaguerrePolynomialRecursion * LaguerrePolynomialsStandard;
  LaguerrePolynomialRecursionAP * LaguerrePolynomials;
#else
  LaguerrePolynomialRecursion * LaguerrePolynomials;
#endif

  // temporary vector for values
  RealVector TmpLaguerre;
  // temporary array for mpf values
  mpf_t *TmpLaguerreArray;
  size_t ArraySize;

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // ratio = ratio between the width in the x direction and the width in the y direction
  // landauLevel = landauLevel to be simulated
  // coulombPrefactor = prefactor multiplying Coulomb terms
  // yukawaParameter = modify Coulomb interaction to a Yukawa form with exponential decay exp(-\lambda r) added
  // nbrPseudopotentials = number of additional pseudopotentials
  // pseudopotentials = additional pseudopotential coefficients
  // perturbationStrength = prefactor multiplying perturbation terms
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  //
  ParticleOnTorusPerturbedCoulombHamiltonian(ParticleOnTorus* particles, int nbrParticles, int maxMomentum, double ratio, int landauLevel, double coulombPrefactor, double yukawaParameter, int nbrPseudopotentials, double* pseudopotentials, double perturbationStrength, AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnTorusPerturbedCoulombHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  //
  AbstractHamiltonian* Clone ();

 private:
 
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  // 
  double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);

  // get Fourier transform of interaction
  //
  // Q2_half = one half of q^2 value
  // return value = Fourier tranform
  //
  double GetVofQ(double Q2_half, double &Precision);

};

#endif
