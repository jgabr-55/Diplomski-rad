#ifndef InterKT_H
#define InterKT_H

#include "InterKTTraits.h"

#include <vector>
#include <set>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

namespace InterKT {

template <typename LorentzMomentum>

class Clustering: public LorentzTraits<LorentzMomentum> {


  /// Alias for the traits class
  typedef LorentzTraits<LorentzMomentum> Traits;

  /// Alias for the Lorentz rotation class
  typedef typename Traits::LorentzRotation LorentzRotation;

  /// Alias for the momentum unit used
  typedef typename Traits::Unit Unit;

  /// Alias for the squared mmentum unit.
  typedef typename Traits::Unit2 Unit2;

  /// For some rason the functions from the base class are not
  /// immediately accesible so we need these 'using' directives.
  using Traits::pt;
  using Traits::pt2;
  using Traits::poslc;
  using Traits::neglc;
  using Traits::x;
  using Traits::y;
  using Traits::z;
  using Traits::e;
  using Traits::m;
  using Traits::create;
  using Traits::createBoost;
  using Traits::add;
  using Traits::subtract;
  using Traits::GeV;
  using Traits::phi;
  using Traits::rap;
  using Traits::getCMBoost;
  using Traits::inverse;
  using Traits::rotate;


  /// Define a vector for sorting
  typedef vector< pair<Unit2,int> > SortVec;

  /// Used for comparison in sorting SortVec
  struct CMP {
    bool operator()(const pair<Unit2, int> & i1,
		    const pair<Unit2, int> & i2) const {
      return i1.first > i2.first;
    }
    bool operator()(Unit2 val, const pair<Unit2, int> & i2) const {
      return val > i2.first;
    }
    bool operator()(const pair<Unit2, int> & i1, Unit2 val) const {
      return i1.first > val;
    }
  };


public:

  /// The default constructor
  Clustering()
    : maxkt2(Unit2()), nmin(2), nmax(-1), maxMIcos(-0.98), maxMIpt2(Unit2()),
      minkt2(Unit2()), cutkt(Unit()), nclus(0), N(0),
      debug(false), speedy(false), sequence(false), sorted(false), aript(false) {}

  /// The destructor
  ~Clustering() {}

private:

  /// Class representing a cluster (pseudo jet)
  struct Cluster {

    /// Default constructor
    Cluster()
      : minkt2(Unit2()), emitter(0), recoiler(0),
	isrdir(0), isleft(true), istouched(true), incluster(-10), track(-1),
	injet(-1) {}

    /// Construct giving a momentum, the parent cluster or parent track
    Cluster(const LorentzMomentum & pin, int in, int tr)
      : p(pin), minkt2(Unit2()), emitter(0), recoiler(0),
	isrdir(0), isleft(true), istouched(true), incluster(in), track(tr),
	injet(-1) {}

    /// The momentum.
    LorentzMomentum p;
    /// The smallest distance found for this cluster.
    Unit2 minkt2;

    /// The emitter corresponding to the smallest distance
    Cluster * emitter;

    /// The recoiler corresponding to the smallest distance
    Cluster * recoiler;

    /// If if not zero, indicate that recoiler is incoming particle
    /// from side a or b (1 and -1 respectively).
    int isrdir;

    /// Is this cluster still around?
    bool isleft;

    /// Has something happened to this cluster since the last clustering.
    bool istouched;

    /// The parent cluster ie. the original particle to which this
    /// cluster belongs, or, if it has been clustered, the original
    /// particle to which it was clustered.
    int incluster;

    /// The original track
    int track;

    /// The final jet to which this cluster belong.
    int injet;

  };


  /// A class corresponding to a final jet
  struct Jet {

    /// The default constructor
    Jet(): cluster(-10) {}

    /// Construct a jet corresponding to a given (index of a) cluster.
    Jet(int ic): cluster(ic) {}

    /// The cluster to which this jet corresponds
    int cluster;

    /// A set of tracks belonging to this jet.
    set<int> tracks;

    /// The momentum of this jet
    LorentzMomentum p;
  };

public:


  /// The maximum distance squared
  Unit2 maxkt2;

  /// The minimum number of jets to be constructed.
  int nmin;

  /// If larger than zero, the maximum number of jets to be constructed.
  int nmax;

  /// The maximum cosine(theta*) of a possible MI clustering to be
  /// considered.
  double maxMIcos;

  /// The maximum summed pt of a pair of clusters to be considered for
  /// MI-clustering.
  Unit2 maxMIpt2;

  /// The smallest distance (transverse momentum) found among the
  /// remaining clusters.
  Unit2 minkt2;

  /// The cutoff in distance. Only tracks above this value will be
  /// considered in the clustering.
  Unit cutkt;

  /// The current number of clusters
  int nclus;

  /// The number of tracks used in the analysis
  int N;

  /// If true, write out information about each clustering
  bool debug;

  ///  If true, write out information about cpu timings.
  bool speedy;

  ///  If true, write out the sequence of clusterings to the file sequence.dat
  bool sequence;

  ///  If true, use the sorted list of indices for each cluster to
  ///  possibly increase speed.
  bool sorted;

  /// Use Ariadne-like transverse momentum as distance for pure
  /// final-state clusterings.
  bool aript;
  

private:

  /// All clusters - alse the ones which have been clustered to others
  vector<Cluster> clusters;

  /// The initial list of tracks
  vector<LorentzMomentum> tracks;

  /// The other final-state particles not included in the jet clustering
  vector<LorentzMomentum> other;

  /// The missing momentum in the event
  LorentzMomentum missing;

  /// The momentum of the incoming parton from one side
  LorentzMomentum Pa;

  /// The momentum of the incoming parton from the other side
  LorentzMomentum Pb;

  /// Total momentum of the hard subsystem
  LorentzMomentum Ptot;

  /// The squared invariant mass of the hard subsystem.
  Unit2 shat;

  /// Flag to indicate that Ptot has changed
  bool hardtouch;

  /// The final jets
  vector<Jet> jets;

  /// Cache calculated invariant masses
  vector< vector<Unit2> > scache;

  ///  A vector of vector of sorted indices for every pair of clusters
  vector<SortVec> ssort;

  /// A vector of the largest possible invariant mass squared for each cluster.
  vector<Unit2> smaxcache;

  /// The emitter in the last clustering
  int lastemitter;

  /// The recoiler in the last clustering
  int lastrecoiler;

  /// The the cluster removed in the last clustering
  int lastrem1;

  /// If two clusters were removed in the last clustering, the second is given here.
  int lastrem2;

public:

  /// Return the elapsed cpu time
  double elapsed();


  /// The main function to be called from the outside. Gives the
  /// momenta of the tracks to be clusterd, and optionally the momenta
  /// of other particles in the event
  int cluster(const vector<LorentzMomentum> & tracksin,
	      ofstream & test,
              const vector<LorentzMomentum> & otherin =
	                             vector<LorentzMomentum>());

  /// Return the final jets 
  vector<LorentzMomentum> getJets();

private:

  /// Initialize the cache of invariant masses squared between all pairs of
  /// clusters.
  void initscache();

  /// Update the cache of invariant masses squared between all pairs
  /// of clusters.
  void updatescache();

  /// Update the cache of invatiant masses w.r.t. beam a.
  void updatescachea();

  ///  Update the cache of invatiant masses w.r.t. beam b.
  void updatescacheb();

  /// For each cluster, sort the list of indices to other cluster
  /// ordered in invariant mass.
  void sorts();

  /// Return the cached value of invariant mass squared of cluster i and j
  Unit2 sij(int i, int j) const {
    return scache[i][j];
  }


  /// Return the cached value of invariant mass squared of cluster i and beam a
  Unit2 saj(int j) const {
    return scache[clusters.size()][j];
  }

  /// Return the cached value of invariant mass squared of cluster i and beam b
  Unit2 sbj(int j) const {
    return scache[clusters.size() + 1][j];
  }

  /// Construct the final jets
  void filljets();

  /// Add this cluster to the jet in which it ended up.
  int addtojets(int ic);

  /// Perform the clustering according to the latest call to update
  /// for this cluster.
  int cluster(Cluster & c, ofstream & test);

  /// Return true if p1 and p2 are approximately the same
  bool approx(LorentzMomentum p1, const LorentzMomentum & p2,
	      Unit scale, double eps = 1.0e-10);

  /// Update the selected clustering for cluster number ic, going
  /// through all clusters that has changed since last
  /// clustering. Only clusterings in which the number of clusters is
  /// decreased less that maxrem are considered
  void update(int ic, int maxrem);

  /// Print a Lorentz momentum for debugging purposes
  void print(const LorentzMomentum & p) const;

  /// Print distance between Lorentz momenta for debugging purposes
  void printdr(const LorentzMomentum & p, const LorentzMomentum & q) const;

  /// Print distance between Lorentz momenta for debugging purposes
  double deltar(const LorentzMomentum & p, const LorentzMomentum & q) const;

  /// Used to avoid very small energies in printouts.
  static double roundGeV(Unit x);

  /// Write out information about a possible clustering before it is performed.
  void debugbefore(Cluster & c, ofstream & test) const;

  /// Write out information about a clustering after it was performed.
  void debugafter(Cluster & c) const;

  /// The transverse momentum, according to Pythia, for a final-state
  /// splitting where the total invariant mass squared is given by s,
  /// the invariant mass squared between the emitted and the emitter
  /// is sec, and the invariant mass squared of the emitter (emitted)
  /// and the recoiler is src (ser).
  static Unit2 invpt2(Unit2 s, Unit2 sce, Unit2 scr, Unit2 ser) {
    return sce*(s - scr)*(s - ser)/(2.0*s - scr - ser)/(2.0*s - scr - ser);
  }

  /// The transverse momentum, according to Pythia, for a final-state
  /// splitting where the invariant mass squared between the emitted
  /// and the emitter is sec, and the invariant mass squared of the
  /// emitter (emitted) and the recoiler is src (ser). This is
  /// equivalent to the previous version with four argument assuming
  /// that all particles are massless, in which case
  /// s = sce + scr + ser.
  static Unit2 invpt2(Unit2 sce, Unit2 scr, Unit2 ser) {
    return sce*(sce + ser)*(sce + scr)/
      ((2.0*sce + scr + ser)*(2.0*sce + scr + ser));
  }

  /// The invariant pt2 according to a modified ariadne measure.
  inline static Unit2 invpt2A(Unit2 sce, Unit2 scr, Unit2 ser) {
    return sce*min(scr, ser)/(sce + scr + ser);
  }

  ///  The invariant pt2 according to the ariadne measure.
  inline static Unit2 invpt2a(Unit2 sce, Unit2 scr, Unit2 ser) {
    return sce*scr/(sce + scr + ser);
  }

  /// Ratio of light-cone components for a momenta
  static double exp2y(LorentzMomentum p) {
    return poslc(p)/neglc(p);
  }


  /// Check if the given clusters should be selected if considered as
  /// a multiple-interaction.
  void checkMI(Cluster & c, Cluster & em);

  /// Check if the cluster c should be considered emitted from the
  /// cluster em as a final-state splitting with an initial state
  /// recoiler.
  void checkIF(Cluster & c, Cluster & em, int ic, int ie);

  /// Check if the cluster c should be considered emitted from the
  /// cluster em as a final-state splitting with an final state
  /// recoiler.
  void checkFS(Cluster & c, Cluster & em, int ic, int ie, bool checkall); 

  /// Alternate version of checkFS
  void checkFS2(Cluster & c, Cluster & em, int ic, int ie, bool checkall); 

  /// Alternate version of checkFS
  void checkFS3(Cluster & c, Cluster & em, int ic, int ie, bool checkall); 

  /// Alternate version of checkFS with ariande measure
  void checkAFS(Cluster & c, Cluster & em, int ic, int ie, bool checkall); 

  /// Alternate version of checkFS with ariande measure
  void checkAFS2(Cluster & c, Cluster & em, int ic, int ie, bool checkall); 

  /// Helper function for checkFS
  void checklastFS(Cluster & c, Cluster & em, int ic, int ie, int last); 

  /// Helper function for checkAFS
  void checklastAFS(Cluster & c, Cluster & em, int ic, int ie, int last); 

};


/// Very crude pre-clustering putting all tracks in a calorimeter with
/// fixed eta-phi resolution and useing the calorimeter towers as
/// tracks.
template <typename LorentzMomentum>
class Calorimeter: public LorentzTraits<LorentzMomentum> {

  typedef LorentzTraits<LorentzMomentum> Traits;

  typedef typename Traits::LorentzRotation LorentzRotation;
  typedef typename Traits::Unit Unit;

public:

  Calorimeter(int phibins, int etabins, double etamin, double etamax):
    cells(vector< vector<Unit> >(etabins, vector<Unit>(phibins, Unit()))),
    neta(etabins), nphi(phibins), mineta(etamin), maxeta(etamax) {}


public:

  void clear() {
    cells = vector< vector<Unit> >(neta, vector<Unit>(nphi, Unit()));
  }

  void fill(const vector<LorentzMomentum> & tracks) {
    for ( int i = 0, N = tracks.size(); i < N; ++i ) {
      Unit ptot = sqrt(pt2(tracks[i]) + z(tracks[i])*z(tracks[i]));
      double eta = log(pt(tracks[i])/abs(z(tracks[i])));
      if ( z(tracks[i]) < Unit() ) eta = -eta;
      if ( eta < mineta ) continue;
      int ieta = int(double(neta)*(eta - mineta)/(maxeta - mineta));
      if ( ieta >= neta ) continue;
      double phi = atan2(x(tracks[i])/e(tracks[i]), y(tracks[i])/e(tracks[i]));
      if ( phi < 0.0 ) phi += 2.0*M_PI;
      if ( phi >= 2.0*M_PI ) phi -= 2.0*M_PI;
      int iphi = int(double(nphi)*phi/(2.0*M_PI));
      cells[ieta][iphi] += e(tracks[i])*pt(tracks[i])/ptot;
    }
  }
   

  vector<LorentzMomentum> clusters() {
    vector<LorentzMomentum> ret;
    for ( int ieta = 0; ieta < neta; ++ieta )
      for ( int iphi = 0; iphi < nphi; ++ iphi ) {
	if ( cells[ieta][iphi] <= Unit() ) continue;
	double eta = mineta +
	  (double(ieta) + 0.5)*(maxeta - mineta)/double(neta);
	double phi = (double(iphi) + 0.5)*2.0*M_PI/double(nphi);
	ret.push_back(create(cells[ieta][iphi]*sin(phi),
			     cells[ieta][iphi]*cos(phi),
			     cells[ieta][iphi]*sinh(eta),
			     cells[ieta][iphi]*cosh(eta)));
      }
    return ret;
  }

private:

  vector< vector<Unit> > cells;

  int neta;
  int nphi;
  double mineta;
  double maxeta;  

};

}

#include "InterKT.cc" /// promjenjeno iz .tcc u .cc

#endif

