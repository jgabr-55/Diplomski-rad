#ifndef InterKTPythia8
#define InterKTPythia8
#include "InterKTTraits.h"

namespace InterKT {	

template<>
struct LorentzTraits<Pythia8::Vec4>:
    public LorentzTraitsBase<Pythia8::Vec4, double, double, 
                             Pythia8::RotBstMatrix,
                             LorentzTraits<Pythia8::Vec4> > {

  typedef Pythia8::RotBstMatrix LorentzRotation;
  typedef Pythia8::Vec4 LorentzMomentum;
  typedef double Unit;
  typedef double Unit2;

  /**
   * Return the x component of \a p.
   */
  static Unit x(const LorentzMomentum & p) {
    return p.px();
  }

  /**
   * Return the y component of \a p.
   */
  static Unit y(const LorentzMomentum & p) {
    return p.py();
  }

  /**
   * Return the z component of \a p.
   */
  static Unit z(const LorentzMomentum & p) {
    return p.pz();
  }

  /**
   * Return the (signed) invariant mass of \a p.
   */
  static Unit m(const LorentzMomentum & p) {
    return p.mCalc();
  }

  /**
   * Return the (signed) invariant mass of \a p.
   */
  static Unit m2(const LorentzMomentum & p) {
    return p.m2Calc();
  }

  /**
   * Return the Lorentz momentum \a p transformed according to the
   * Lorentz rotation \a r.
   */
  static LorentzMomentum rotate(const LorentzMomentum & p,
				const LorentzRotation & r) {
    LorentzMomentum ret = p;
    ret.rotbst(r);
    return ret;
  }

  /**
   * Return the Lorentz rotation \a r2 transformed according to the
   * Lorentz rotation \a r1. For a given lorentz vector, p,
   * rotate(roatate(p, r2),r1) is equivalent to rotate(p, rotate(r1, r2)).
   */
  static LorentzRotation rotate(const LorentzRotation & r1,
				const LorentzRotation & r2) {
    LorentzRotation ret = r2;
    ret.rotbst(r1);
    return ret;
  }

  /**
   * Transform a given Lorentz rotation with a rotation in
   * polar \a angle.
   */
  static void rotateX(LorentzRotation & r, double angle) {
    r.rot(angle, 0.0);
  }

  /**
   * Transform a given Lorentz rotation with a rotation in
   * azimuth \a angle.
   */
  static void rotateZ(LorentzRotation & r, double angle) {
    r.rot(0.0, angle);
  }

  /**
   * Return the Lorentz rotation corresponding to the given boost
   * vector.
   */
  static LorentzRotation createBoost(double bx, double by, double bz) {
    LorentzRotation r;
    r.bst(bx, by, bz);
    return r;
  }

  /**
   * Return the inverse of the given LorentzRotation.
   */
  static LorentzRotation inverse(const LorentzRotation & r) {
    LorentzRotation ret(r);
    ret.invert();
    return ret;
  }

};

}


#endif

