// -*- C++ -*-
//
// MENewScalars.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2009-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MENewScalars class.
//

#include "MENewScalars.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/Utilities/Interpolator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <fstream>

using namespace Herwig;

DescribeClass<MENewScalars,HwMEBase>
describeHerwigMENewScalars("Herwig::MENewScalars","MENewScalars.so");

MENewScalars::MENewScalars()
  : _ggiota(1.0) {
  massOption({{0,0}});
}

void MENewScalars::doinit() {
  HwMEBase::doinit();
  PDPtr iota0 = getParticleData(99926);
  miota_ = iota0->mass();
  wiota_ = iota0->generateWidth(miota_);
  iotamass_=dynamic_ptr_cast<GenericMassGeneratorPtr>(iota0->massGenerator());

  
}

void MENewScalars::rebind(const TranslationMap & trans) {
  HwMEBase::rebind(trans);
}

IVector MENewScalars::getReferences() {
  IVector ret = HwMEBase::getReferences();
  return ret;
}



IBPtr MENewScalars::clone() const {
    return new_ptr(*this);
  }
  
IBPtr MENewScalars::fullclone() const {
    return new_ptr(*this);
  }

void MENewScalars::persistentOutput(PersistentOStream & os) const {
  os << _ggiota << iotamass_;
}

void MENewScalars::persistentInput(PersistentIStream & is, int) {
  is  >> _ggiota >> iotamass_;
}

Energy2 MENewScalars::scale() const {
  return sHat();
}

void MENewScalars::Init() {

  static ClassDocumentation<MENewScalars> documentation
    ("The MENewScalars class implements a simple gluon-fusion heavy scalar production process"
     " collisions");
 
 static Parameter<MENewScalars, double> interfaceGGiota
    ("GGiota",
     "The gluon-gluon-heavy scalar coupling",
     &MENewScalars::_ggiota, 1.0, -1000000.0, 1000000.0,
     false, false, Interface::limited);


}

Selector<MEBase::DiagramIndex>
MENewScalars::diagrams(const DiagramVector & diags) const {
  // select the diagram, this is easy for us as we have already done it
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    sel.insert(1.0, i);
  }
  return sel;
}

void MENewScalars::getDiagrams() const {
  // get the particle data objects
  PDPtr gluon = getParticleData(ParticleID::g);
  PDPtr iota0 = getParticleData(99926);
  add(new_ptr((Tree2toNDiagram(2), gluon, gluon, 1, iota0, -1)));
}

Selector<const ColourLines *>
MENewScalars::colourGeometries(tcDiagPtr diag) const {
  // gg to iota
  static const ColourLines line1("1 -2,2 -1");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &line1);
  return sel;
}


double MENewScalars::me2() const {
  //enforce onshell production:
  Energy iotamass = meMomenta()[2].m();
  tcPDPtr iota0 = mePartonData()[2];
  Energy mass = iota0->mass();
  Energy deltamass = 1E-4*GeV;
  if (.0*GeV > iotamass) return 0.0;
  using Constants::pi;
  InvEnergy2 bwfact;
  bwfact = mePartonData()[2]->generateWidth(sqrt(sHat()))*sqrt(sHat())/pi/
      (sqr(sHat()-sqr(miota_))+sqr(miota_*wiota_));
  
  if (iotamass > mass + deltamass || iotamass < mass - deltamass) return 0.0;
  return sqr(_ggiota) * double(UnitRemoval::E2 * bwfact);;
   
}

bool MENewScalars::generateKinematics(const double * r) {
  Lorentz5Momentum pout = meMomenta()[0] + meMomenta()[1];
  pout.rescaleMass();

  meMomenta()[2].setMass(pout.mass());
  meMomenta()[2] = LorentzMomentum(pout.x(),pout.y(),pout.z(),pout.t());
  jacobian(1.0);
  // check whether it passes all the cuts: returns true if it does

  vector<LorentzMomentum> out(1,meMomenta()[2]);
  tcPDVector tout(1,mePartonData()[2]);
  
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}
