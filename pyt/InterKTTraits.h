#ifndef InterKTTraits_H
#define InterKTTraits_H

namespace InterKT {

using namespace std;

/**
 * The LorentzTraitsBase class defines the interface to a system of
 * unit, Lorentz momentum and rotation classes to be used by the
 * InterKT::Clustering class. Specializations should specialize the
 * LorentzTraits and inherit from the corresponding LorentzTraitsBase
 * to use the default implementations of the Traits functions.
 */
template <typename LorentzMomentum, typename UnitType, typename UnitType2,
  typename LorentzRotationType, typename Super>
struct LorentzTraitsBase {

  /**
   * The unit class used by the corresponding LorentzMomentum class.
   */
  typedef UnitType Unit;

  /**
   * The squared unit class used by the corresponding LorentzMomentum
   * class.
   */
  typedef UnitType2 Unit2;

  /**
   * The Lorentz rotation class used by the corresponding
   * LorentzMomentum class.
   */
  typedef LorentzRotationType LorentzRotation;

  /**
   * Return variable in units of GeV.
   */
  static double GeV(Unit x) {
    return x;
  }

  /**
   * Return the x component of \a p.
   */
  static Unit x(const LorentzMomentum & p) {
    return p.x();
  }

  /**
   * Return the y component of \a p.
   */
  static Unit y(const LorentzMomentum & p) {
    return p.y();
  }

  /**
   * Return the z component of \a p.
   */
  static Unit z(const LorentzMomentum & p) {
    return p.z();
  }

  /**
   * Return the energy component of \a p.
   */
  static Unit e(const LorentzMomentum & p) {
    return p.e();
  }

  /**
   * Return the positive light-cone component of \a p.
   */
  static Unit poslc(const LorentzMomentum & p) {
    return Super::e(p) + Super::z(p);
  }

  /**
   * Return the negative light-cone component of \a p.
   */
  static Unit neglc(const LorentzMomentum & p) {
    return Super::e(p) - Super::z(p);
  }

  /**
   * Return the rapidity of \a p.
   */
  static double rap(const LorentzMomentum & p) {
    double rap = -10000.0;
    if ( Super::poslc(p) > Super::neglc(p) ) rap = 10000.0;
    if ( Super::poslc(p)*Super::neglc(p) > Unit2() )
      rap = 0.5*log(Super::poslc(p)/Super::neglc(p));
    return rap;
  }

  /**
   * Return the rapidity of \a p.
   */
  static double phi(const LorentzMomentum & p) {
    return p.phi();
  }

  /**
   * Return the transverse component of \a p.
   */
  static Unit pt(const LorentzMomentum & p) {
    return sqrt(Super::pt2(p));
  }

  /**
   * Return the squared transverse component of \a p.
   */
  static Unit2 pt2(const LorentzMomentum & p) {
    return Super::x(p)*Super::x(p) + Super::y(p)*Super::y(p);
  }

  /**
   * Return the (signed) invariant mass of \a p.
   */
  static Unit m(const LorentzMomentum & p) {
    return p.m();
  }

  /**
   * Return the squared invariant mass of \a p.
   */
  static Unit2 m2(const LorentzMomentum & p) {
    return p.m2();
  }

  /**
   * Return a Lorentz momentum with the given components.
   */
  static LorentzMomentum create(Unit px, Unit py, Unit pz, Unit pe) {
    return LorentzMomentum(px, py, pz, pe);
  }

  /**
   * Return the Lorentz momentum \a p transformed according to the
   * Lorentz rotation \a r.
   */
  static LorentzMomentum rotate(const LorentzMomentum & p,
			     const LorentzRotation & r) {
    return r*p;
  }

  /**
   * Return the Lorentz rotation \a r2 transformed according to the
   * Lorentz rotation \a r1. For a given lorentz vector, p,
   * rotate(roatate(p, r2),r1) is equivalent to rotate(p, rotate(r1, r2)).
   */
  static LorentzRotation rotate(const LorentzRotation & r1,
				const LorentzRotation & r2) {
    return r1*r2;
  }

  /**
   * Return the Lorentz rotation corresponding to the given boost
   * vector.
   */
  static LorentzRotation createBoost(double bx, double by, double bz) {
    return LorentzRotation(bx, by, bz);
  }

  /**
   * Transform a given Lorentz rotation with a rotation in
   * polar \a angle.
   */
  static void rotateX(LorentzRotation & r, double angle) {
    r.rotateX(angle);
  }

  /**
   * Transform a given Lorentz rotation with a rotation in
   * azimuth \a angle.
   */
  static void rotateZ(LorentzRotation & r, double angle) {
    r.rotateZ(angle);
  }

  /**
   * Return the sum of two momenta.
   */
  static LorentzMomentum add(const LorentzMomentum & p1,
			    const LorentzMomentum & p2) {
    return p1 + p2;
  }

  /**
   * Return the difference between two momenta.
   */
  static LorentzMomentum subtract(const LorentzMomentum & p1,
			    const LorentzMomentum & p2) {
    return p1 - p2;
  }

  /**
   * Return a Lorentz rotation corresponding to boosting \a p1, \a p2
   * and \a p3 to their rest frame following by a rotation putting \a
   * p3 along the negative z-axis.
   */
  static LorentzRotation getCMBoost(const LorentzMomentum & p1,
				    const LorentzMomentum & p2,
				    const LorentzMomentum & p3) {
    LorentzMomentum ptot = Super::add(Super::add(p1, p2), p3);
    LorentzRotation r =
      Super::createBoost(-Super::x(ptot)/Super::e(ptot),
			 -Super::y(ptot)/Super::e(ptot),
			 -Super::z(ptot)/Super::e(ptot));

    LorentzMomentum pp = Super::rotate(p3, r);
    double phi = atan2(Super::y(pp)/Super::e(ptot),
		       Super::x(pp)/Super::e(ptot));
    Super::rotateZ(r, -phi);
    double the = atan2(Super::pt(pp)/Super::e(ptot),
		       Super::z(pp)/Super::e(ptot));
    Super::rotateX(r, M_PI - the);
    
    return r;
  }
				    
  /**
   * Return a Lorentz rotation corresponding to boosting \a p1, \a p2
   * and \a p3 to ther rest frame following by a rotation putting \a
   * p2 along the negative z-axis and \a p1 along the positive one.
   * on the xz-plane.
   */
  static LorentzRotation getCMBoost(const LorentzMomentum & p1,
				    const LorentzMomentum & p2) {
    LorentzMomentum ptot = Super::add(p1, p2);
    LorentzRotation r =
      Super::createBoost(-Super::x(ptot)/Super::e(ptot),
			 -Super::y(ptot)/Super::e(ptot),
			 -Super::z(ptot)/Super::e(ptot));


    LorentzMomentum pp = Super::rotate(p1, r);
    double phi = atan2(Super::y(pp)/Super::e(ptot),
		       Super::x(pp)/Super::e(ptot));
    Super::rotateZ(r, -phi);
    double the = atan2(Super::pt(pp)/Super::e(ptot),
		       Super::z(pp)/Super::e(ptot));
    Super::rotateX(r, -the);
    Super::rotateZ(r, phi);

    return r;
  }
				    
  /**
   * Return the inverse of the given LorentzRotation.
   */
  static LorentzRotation inverse(const LorentzRotation & r) {
    return r.inverse();
  }

};

//   typedef typename LorentzMomentum::LorentzRotation LorentzRotation;

//   typedef typename LorentzMomentum::Unit Unit;


template <typename LorentzMomentum>
struct LorentzTraits:
    public LorentzTraitsBase<LorentzMomentum,
                             typename LorentzMomentum::Unit,
                             typename LorentzMomentum::Unit2,
                             typename LorentzMomentum::LorentzRotation,
                             LorentzTraits<LorentzMomentum> > {};


}

#endif
