
template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::update(int ic, int maxrem) {
  // Take the Cluster with index ic and check if the selected
  // clustering should be updated with a new one. The selected
  // clustering must at most remove maxrem clusters.
  Cluster & c = clusters[ic];

  // First check if the previous selected clustering has been affected
  // and we therefore need to check all possible clustering. Otherwise
  // we only need to check clusterings that may have been affected
  // after the last step.
  bool checkall = c.istouched;
  if ( c.emitter && ( c.emitter->istouched || !c.emitter->isleft ) )
    checkall = true;
  if ( c.recoiler && ( c.recoiler->istouched || !c.recoiler->isleft ) )
    checkall = true;
  if ( ( !c.emitter || !c.recoiler ) && hardtouch ) checkall = true;
  if ( maxrem < 2 && c.emitter && !c.recoiler && !c.isrdir ) checkall = true;

  if ( checkall ) {
    c.emitter = 0;
    c.recoiler = 0;
    c.isrdir = 0;
    c.minkt2 = shat;
  }

  // If the last step involved a global change of the outgoing
  // clusters we need to check possible initial-state clusterings.
  if ( checkall || hardtouch ) {
    // Start out checking ISR from a-side.
    double zfac = 1.0 - m2(subtract(Ptot, c.p))/shat;
    Unit2 kt2i = -zfac*m2(subtract(Pa, c.p));
    if ( kt2i < Unit2() ) kt2i = sqrt(kt2i)*Unit(); // Signal error where it occurs.
    if ( kt2i < c.minkt2 ) {
      c.minkt2 = kt2i;
      c.isrdir = 1;
      c.emitter = 0;
      c.recoiler = 0;
    }

    // Now check ISR from b-side.
    kt2i = -zfac*m2(subtract(Pb, c.p));
    if ( kt2i < Unit2() ) kt2i = sqrt(kt2i)*Unit(); // Signal error where it occurs.
    if ( kt2i < c.minkt2 ) {
      c.minkt2 = kt2i;
      c.isrdir = -1;
      c.emitter = 0;
      c.recoiler = 0;
    }
  }

  if ( checkall || hardtouch || lastemitter >= 0 || lastrecoiler >= 0 ) {
    if ( sorted ) {
      //    for ( int ie = ic + 1; ie < N; ++ie ) {
      for ( int je = ssort[ic].size() -1; je >=0; --je ) {
	int ie = ssort[ic][je].second;
	if ( ie <= ic ) continue;
	Cluster & em = clusters[ie];
	
	if ( checkall || hardtouch ||
	     ie == lastemitter || ie == lastrecoiler )  {
	  // Check MI
	  if ( maxrem > 1 && maxMIcos > -1.0 && maxMIpt2 > Unit2() ) 
	    checkMI(c, em);
	  checkIF(c, em, ic, ie);
	}
	if ( checkall || lastemitter >= 0 || lastrecoiler >= 0 ) {
	  if ( aript )
	    checkAFS2(c, em, ic, ie, checkall);
	  else
	    checkFS3(c, em, ic, ie, checkall);
	}
      }
    } else {
      for ( int ie = ic + 1; ie < N; ++ie ) {
	Cluster & em = clusters[ie];
	if ( !em.isleft ) continue;

	if ( checkall || hardtouch ||
	     ie == lastemitter || ie == lastrecoiler )  {
	  // Check MI
	  if ( maxrem > 1 && maxMIcos > -1.0 && maxMIpt2 > Unit2() ) 
	    checkMI(c, em);
	  checkIF(c, em, ic, ie);
	}
	if ( checkall || lastemitter >= 0 || lastrecoiler >= 0 ) {
	  if ( aript ) 
	    checkAFS(c, em, ic, ie, checkall);
	  else
	    checkFS(c, em, ic, ie, checkall);
	}
      }
    }
  }
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::
checkFS(Cluster & c, Cluster & em, int ic, int ie, bool checkall) {
  Unit2 sce = sij(ic, ie);
  Unit2 smax = sce + smaxcache[ic] + smaxcache[ie];
  Unit2 pt2min = 2.0*sce*sce*(smax - sce)/(smax + sce)/(smax + sce);
  if ( pt2min > c.minkt2 ) return;
  if ( checkall || ie == lastemitter || ie == lastrecoiler ) {
    for ( int ir = 0; ir < N; ++ir ) {
      Cluster & re = clusters[ir];
      const Unit2 scr = sij(ic, ir);
      const Unit2 ser = sij(ie, ir);
      if ( ir != ie && ir != ic && re.isleft && sce <= scr && sce <= ser ) {
	Unit2 pt2i = invpt2(sce, scr, ser);
	if ( pt2i < c.minkt2 ) {
	  c.minkt2 = pt2i;
	  c.isrdir = 0;
	  c.emitter = &em;
	  c.recoiler = &re;
	  if ( pt2min > pt2i ) return;
	}
      }
    }
    return;
  }
  checklastFS(c, em, ic, ie, lastemitter);
  checklastFS(c, em, ic, ie, lastrecoiler);
  
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::
checkAFS(Cluster & c, Cluster & em, int ic, int ie, bool checkall) {
  Unit2 sce = sij(ic, ie);
  Unit2 smax = 2.0*sce + max(smaxcache[ic], smaxcache[ie]);
  Unit2 pt2min = sce*sce/smax;
  if ( pt2min > c.minkt2 ) return;
  if ( checkall || ie == lastemitter || ie == lastrecoiler ) {
    for ( int ir = 0; ir < N; ++ir ) {
      Cluster & re = clusters[ir];
      const Unit2 scr = sij(ic, ir);
      const Unit2 ser = sij(ie, ir);
      if ( ir != ie && ir != ic && re.isleft && sce <= scr && sce <= ser ) {
	Unit2 pt2i = invpt2A(sce, scr, ser);
	if ( pt2i < c.minkt2 ) {
	  c.minkt2 = pt2i;
	  c.isrdir = 0;
	  c.emitter = &em;
	  c.recoiler = &re;
	  if ( pt2min > pt2i ) return;
	}
      }
    }
    return;
  }
  checklastAFS(c, em, ic, ie, lastemitter);
  checklastAFS(c, em, ic, ie, lastrecoiler);
  
}


template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::
checkFS2(Cluster & c, Cluster & em, int ic, int ie, bool checkall) {
  typedef typename SortVec::iterator Sit;

  Unit2 sce = sij(ic, ie);


  if ( checkall || ie == lastemitter || ie == lastrecoiler ) {

    Sit itj = lower_bound(ssort[ic].begin(), ssort[ic].end(), sce, CMP());
    int jcrmin = ssort[ic].end() - itj;
    itj = lower_bound(ssort[ie].begin(), ssort[ie].end(), sce, CMP());
    int jermin = ssort[ie].end() - itj;
    int jcr = ssort[ic].size() - 1;
    int jer = ssort[ie].size() - 1;

    while ( jcr >= jcrmin && jer >= jermin ) {
      int icr = ssort[ic][jcr].second;
      Unit2 scr = ssort[ic][jcr].first;
      if ( scr < sce ) return;
      if ( icr == ie ) {
	--jcr;
	continue;
      }
      int ier = ssort[ie][jer].second;
      Unit2 ser = ssort[ie][jer].first;
      if ( ser < sce ) return ;
      if ( ier == ic || ier == icr ) {
	--jer;
	continue;
      }

      Unit2 smax = sce + scr + ser;
      if ( 2.0*sce*sce*(smax - sce)/(smax + sce)/(smax + sce) > c.minkt2 )
	return;
      if ( scr > ser ) {
	Unit2 pt2i = invpt2(sce, scr, sij(ie, icr));
	if ( pt2i < c.minkt2 ) {
	  c.minkt2 = pt2i;
	  c.isrdir = 0;
	  c.emitter = &em;
	  c.recoiler = &clusters[icr];
	}
	--jcr;
      } else {
	Unit2 pt2i = invpt2(sce, ser, sij(ic, ier));
	if ( pt2i < c.minkt2 ) {
	  c.minkt2 = pt2i;
	  c.isrdir = 0;
	  c.emitter = &em;
	  c.recoiler = &clusters[ier];
	}
	--jer;
      }
    }
    return;
  }
  checklastFS(c, em, ic, ie, lastemitter);
  checklastFS(c, em, ic, ie, lastrecoiler);  
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::
checkFS3(Cluster & c, Cluster & em, int ic, int ie, bool checkall) {
  typedef typename SortVec::iterator Sit;

  Unit2 sce = sij(ic, ie);


  if ( checkall || ie == lastemitter || ie == lastrecoiler ) {

    const SortVec & sv = ssort[ic];

    Unit2 smax = sce + sv[0].first + ssort[ie][0].first;
    Unit2 pt2min = 2.0*sce*sce*(smax - sce)/(smax + sce)/(smax + sce);
    if ( pt2min > c.minkt2 ) return;

    for ( int jcr = 0, Nj = sv.size(); jcr < Nj; ++jcr ) {
      int icr = sv[jcr].second;
      Unit2 scr = sv[jcr].first;
      if ( scr < sce ) return;
      if ( icr == ie || icr == ic ) continue;
      Unit2 ser = sij(ie, icr);
      if ( ser >= sce ) {
	Unit2 pt2i = invpt2(sce, scr, ser);
	if ( pt2i < c.minkt2 ) {
	  c.minkt2 = pt2i;
	  c.isrdir = 0;
	  c.emitter = &em;
	  c.recoiler = &clusters[icr];
	  if ( pt2i < pt2min ) return;
	}
      }
    }
    return;
  }
  checklastFS(c, em, ic, ie, lastemitter);
  checklastFS(c, em, ic, ie, lastrecoiler);  
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::
checkAFS2(Cluster & c, Cluster & em, int ic, int ie, bool checkall) {
  typedef typename SortVec::iterator Sit;

  Unit2 sce = sij(ic, ie);


  if ( checkall || ie == lastemitter || ie == lastrecoiler ) {

    const SortVec & sv = ssort[ic];

    Unit2 smax = 2.0*sce + max(ssort[ie][0].first, sv[0].first);
    Unit2 pt2min = sce*sce/smax;
    if ( pt2min > c.minkt2 ) return;

    for (int jcr = 0, Nj = sv.size(); jcr < Nj; ++jcr ) {
      Unit2 scr = sv[jcr].first;
      if ( scr < sce ) return;
      int icr = sv[jcr].second;
      Unit2 ser = sij(ie, icr);
      if ( icr == ie || ser < sce ) continue;
      Unit2 pt2i = invpt2A(sce, scr, ser);
      if ( pt2i < c.minkt2 ) {
	c.minkt2 = pt2i;
	c.isrdir = 0;
	c.emitter = &em;
	c.recoiler = &clusters[icr];
      }
    }
    return;
  }
  checklastAFS(c, em, ic, ie, lastemitter);
  checklastAFS(c, em, ic, ie, lastrecoiler);    
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::
checklastAFS(Cluster & c, Cluster & em, int ic, int ie, int last) {
 if ( last >= 0 && last != ic && last != ie ) {
    Cluster & re = clusters[last];
    Unit2 sce = sij(ic, ie);
    Unit2 scr = sij(ic, last);
    Unit2 ser = sij(ie, last);
    if ( sce <= scr && sce <= ser ) {
      Unit2 pt2i = invpt2A(sce, scr, ser);
      if ( pt2i < c.minkt2 ) {
	c.minkt2 = pt2i;
	c.isrdir = 0;
	c.emitter = &em;
	c.recoiler = &re;
      }
    }
  }
}
  
template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::
checklastFS(Cluster & c, Cluster & em, int ic, int ie, int last) {
 if ( last >= 0 && last != ic && last != ie ) {
    Cluster & re = clusters[last];
    Unit2 sce = sij(ic, ie);
    Unit2 scr = sij(ic, last);
    Unit2 ser = sij(ie, last);
    if ( sce <= scr && sce <= ser ) {
      Unit2 pt2i = invpt2(sce, scr, ser);
      if ( pt2i < c.minkt2 ) {
	c.minkt2 = pt2i;
	c.isrdir = 0;
	c.emitter = &em;
	c.recoiler = &re;
      }
    }
  }
}
  
template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::checkMI(Cluster & c, Cluster & em) {
  // Check MI
  Unit2 pt2sum = pt2(add(c.p, em.p));
  if  ( pt2sum < maxMIpt2 ) {
    Unit2 pt2a = pt2(em.p);
    Unit2 pt2i = pt2(c.p);
    double cthe = (pt2sum - pt2a - pt2i)/(2.0*sqrt(pt2a*pt2i));
    if ( cthe < maxMIcos && min(pt2a, pt2i) < c.minkt2 ) {
      c.minkt2 = min(pt2a, pt2i);
      c.isrdir = 0;
      c.emitter = &em;
      c.recoiler = 0;
    }
  }
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::
checkIF(Cluster & c, Cluster & em, int ic, int ie) {
  
  // Check initial-final FS dipole (a-side)
  Unit2 scr = saj(ic);
  Unit2 ser = saj(ie);
  Unit2 sce = sij(ic, ie);
  Unit2 pt2tot=pt2(add(c.p,em.p));
  Unit2 pt2i = invpt2(sce, scr, ser);
  if ( pt2i < c.minkt2  && sqrt(pt2i*pt2tot) > sce) {
    c.minkt2 = pt2i;
    c.isrdir = 1;
    c.emitter = &em;
    c.recoiler = 0;
  }
  // Check initial-final FS dipole (b-side)
  scr = sbj(ic);
  ser = sbj(ie);
  pt2i = invpt2(sce, scr, ser);
  if ( pt2i < c.minkt2 && sqrt(pt2i*pt2tot) > sce) {
    c.minkt2 = pt2i;
    c.isrdir = -1;
    c.emitter = &em;
    c.recoiler = 0;
  }
}




template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::initscache() {
  scache = vector< vector<Unit2> >(N + 2, vector<Unit2>(N + 2));
  smaxcache = vector<Unit2>(N);
  ssort = vector<SortVec>(N, SortVec(N));
  for ( int i = 0; i < N; ++i ) {
    for ( int j = i + 1; j < N; ++j ) {
      Unit2 s = m2(add(clusters[i].p, clusters[j].p));
      scache[i][j] = scache[j][i] = s;
      if ( sorted ) ssort[i][j] = make_pair(s, j);
      if ( sorted ) ssort[j][i] = make_pair(s, i);
      smaxcache[i] = max(smaxcache[i], scache[i][j]);
      smaxcache[j] = max(smaxcache[j], scache[i][j]);
    }
    scache[i][i] = m2(clusters[i].p);
    scache[i][N] = scache[N][i] = m2(add(clusters[i].p, Pa));
    scache[i][N+1] = scache[N+1][i] = m2(add(clusters[i].p, Pb));
  }
  scache[N][N+1] = scache[N+1][N] = shat;
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::sorts() {
  typedef typename SortVec::iterator Sit;
  if ( !sorted ) return;
  for ( int i = 0; i < N; ++i )
    if ( clusters[i].isleft ) {
      sort(ssort[i].begin(), ssort[i].end(), CMP());
      while ( ssort[i].back().first <= Unit2() ) ssort[i].pop_back();
    } 
}


template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::updatescache() {

  typedef typename SortVec::iterator Sit;

  if ( lastemitter < 0 && lastrecoiler < 0 && lastrem1 < 0 && lastrem2 < 0 )
    return;
  int i1 = lastemitter;
  int i2 = lastrecoiler;
  int r1 = lastrem1;
  int r2 = lastrem2;
  if ( sorted && i1 >= 0 ) ssort[i1].clear();
  if ( sorted && i2 >= 0 ) ssort[i2].clear();
  
  for ( int j = 0; j < N; ++j ) {
    Sit sitj1 = ssort[j].end();
    Sit sitj2 = ssort[j].end();
    Sit sitrj1 = ssort[j].end();
    Sit sitrj2 = ssort[j].end();
    Unit2 snew1 = Unit2();
    Unit2 snew2 = Unit2();
    if ( i1 >= 0 && i1 != j) {
      Unit2 sold = scache[i1][j];
      snew1 = clusters[j].isleft?
	m2(add(clusters[i1].p, clusters[j].p)): Unit2();
      scache[i1][j] = scache[j][i1] = snew1;
      if ( sorted && snew1 > Unit2() ) ssort[i1].push_back(make_pair(snew1, j));
      if ( sorted && j != i2 && j != r1 && j != r2 ) {
	pair<Sit,Sit>
	  lim = equal_range(ssort[j].begin(), ssort[j].end(), sold, CMP());
	while ( lim.first != lim.second && lim.first->second != i1 )
	  ++lim.first;
	if ( lim.first != lim.second ) sitj1 = lim.first;
      }
      smaxcache[i1] = max(smaxcache[i1], scache[i1][j]);
      smaxcache[j] = max(smaxcache[j], scache[i1][j]);
    }
    if ( i2 >= 0 && i2 != j) {
      Unit2 sold = scache[i2][j];
      snew2 = clusters[j].isleft?
	m2(add(clusters[i2].p, clusters[j].p)): Unit2();
      scache[i2][j] = scache[j][i2] = snew2;
      if ( sorted && snew2 > Unit2() ) ssort[i2].push_back(make_pair(snew2, j));
      if ( sorted && j != i1 && j != r1 && j != r2 ) {
	pair<Sit,Sit>
	  lim = equal_range(ssort[j].begin(), ssort[j].end(), sold, CMP());
	while ( lim.first != lim.second && lim.first->second != i2 ) ++lim.first;
	if ( lim.first != lim.second ) sitj2 = lim.first;
      }
      smaxcache[i2] = max(smaxcache[i2], scache[i2][j]);
      smaxcache[j] = max(smaxcache[j], scache[i2][j]);
    }
    if ( sorted && r1 >= 0 && r1 != j && j != i1 && j != i2 ) {
      Unit2 sold = scache[r1][j];
      pair<Sit,Sit>
	lim = equal_range(ssort[j].begin(), ssort[j].end(), sold, CMP());
      while ( lim.first != lim.second && lim.first->second != r1 ) ++lim.first;
      if ( lim.first != lim.second ) sitrj1 = lim.first;
    }
    if ( sorted && r2 >= 0 && r2 != j && j != i1 && j != i2 ) {
      Unit2 sold = scache[r2][j];
      pair<Sit,Sit>
	lim = equal_range(ssort[j].begin(), ssort[j].end(), sold, CMP());
      while ( lim.first != lim.second && lim.first->second != r2 ) ++lim.first;
      if ( lim.first != lim.second ) sitrj2 = lim.first;
    }
    if ( sorted ) {
      if ( sitj1 != ssort[j].end() ) sitj1->first = snew1;
      if ( sitj2 != ssort[j].end() ) sitj2->first = snew2;
      if ( sitrj1 != ssort[j].end() ) sitrj1->first = Unit2();
      if ( sitrj2 != ssort[j].end() ) sitrj2->first = Unit2();
    }
  }
  if ( i1 < 0 ) return;
  scache[i1][i1] = m2(clusters[i1].p);
  scache[i1][N] = scache[N][i1] = m2(add(clusters[i1].p, Pa));
  scache[i1][N+1] = scache[N+1][i1] = m2(add(clusters[i1].p, Pb));
  if ( i2 >= 0 ) {
    scache[i2][i2] = m2(clusters[i2].p);
    scache[i2][N] = scache[N][i2] = m2(add(clusters[i2].p, Pa));
    scache[i2][N+1] = scache[N+1][i2] = m2(add(clusters[i2].p, Pb));
  }
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::updatescachea() {
  for ( int j = 0; j < N; ++j )
    if ( clusters[j].isleft )
      scache[N][j] = scache[j][N] = m2(add(clusters[j].p, Pa));
  scache[N][N+1] = scache[N+1][N] = shat;
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::updatescacheb() {
  for ( int j = 0; j < N; ++j )
    if ( clusters[j].isleft )
      scache[N+1][j] = scache[j][N+1] = m2(add(clusters[j].p, Pb));
  scache[N][N+1] = scache[N+1][N] = shat;
}

template <typename LorentzMomentum>
int InterKT::Clustering<LorentzMomentum>::cluster(Cluster & c, ofstream & test) {
  debugbefore(c, test);
  if ( c.emitter && c.recoiler ) {
    LorentzRotation r = inverse(getCMBoost(c.p, c.emitter->p, c.recoiler->p));
    LorentzMomentum ptot = add(add(c.p, c.emitter->p), c.recoiler->p);
    Unit ee = 0.5*m(ptot);
    c.emitter->p = rotate(create(Unit(), Unit(), ee, ee), r);
    c.recoiler->p = rotate(create(Unit(), Unit(), -ee, ee), r);
    c.isleft = false;
    lastrem1 = c.incluster;
    c.incluster = c.emitter->incluster;
    c.emitter->istouched = true;
    lastemitter = c.emitter->incluster;
    c.recoiler->istouched = true;
    lastrecoiler = c.recoiler->incluster;
    return 1;
  }
  LorentzMomentum Paold = Pa;
  LorentzMomentum Pbold = Pb;
  LorentzMomentum Ptotnew = Ptot;
  int nrem = 1;
  if ( c.emitter && c.isrdir ) {
    LorentzMomentum & pr = c.isrdir > 0? Pa: Pb;
    LorentzRotation r = inverse(getCMBoost(c.p, c.emitter->p, pr));
    LorentzMomentum psum = add(add(c.p, c.emitter->p), pr);
    Unit ee = 0.5*m(psum);
    c.emitter->p = rotate(create(Unit(), Unit(), ee, ee), r);
    LorentzMomentum prp = rotate(create(Unit(), Unit(), -ee, ee), r);
    pr = subtract(add(pr, pr), prp);//????????
    Ptotnew = add(Pa, Pb);
    c.isleft = false;
    lastrem1 = c.incluster;
    c.incluster = c.emitter->incluster;
    c.emitter->istouched = true;
    lastemitter = c.emitter->incluster;
    hardtouch = true;
  } else if ( c.emitter ) {
    LorentzMomentum pp = add(c.p, c.emitter->p);
    LorentzMomentum ph = subtract(Ptot, pp);
    Ptotnew = ph;
    LorentzMomentum & pe = (exp2y(pp) > exp2y(ph))? Pa: Pb;
    pe = subtract(pe, pp); //????????????????
    c.isleft = false;
    lastrem1 = c.incluster;
    c.incluster = ( exp2y(pp) > exp2y(ph) )? -1: -2;
    c.emitter->isleft = false;
    lastrem2 = c.emitter->incluster;
    c.emitter->incluster = c.incluster;
    hardtouch = true;
    nrem = 2;
  } else {
    LorentzMomentum & pe = c.isrdir > 0? Pa: Pb;
    pe = subtract(pe, c.p);
    c.isleft = false;
    lastrem1 = c.incluster;
    c.incluster = c.isrdir > 0? -1: -2;
    hardtouch = true;
  }
  LorentzRotation R2 = getCMBoost(Pa, Pb);
  shat = m2(add(Pa, Pb));
  Unit ee = 0.5*m(add(Pa, Pb));
  Pa = create(Unit(), Unit(),  ee, ee);
  Pb = create(Unit(), Unit(), -ee, ee);
  LorentzRotation R3;
  if ( c.isrdir > 0 ) {
    R3 = createBoost(0.0, 0.0, -(neglc(Pbold)*neglc(Pbold) - shat)/
		     (neglc(Pbold)*neglc(Pbold) + shat));
  } else if ( c.isrdir < 0 ) {
      R3 = createBoost(0.0, 0.0,  (poslc(Paold)*poslc(Paold) - shat)/
		       (poslc(Paold)*poslc(Paold) + shat));
  } else {
    R3 = createBoost(0.0, 0.0, (exp2y(Ptotnew) - 1.0)/(exp2y(Ptotnew) + 1.0));
  }
  Paold = rotate(Pa, R3);
  Pbold = rotate(Pb, R3);
  LorentzRotation R = rotate(R3, R2);
  Ptot = (missing = rotate(missing, R));
  for ( int i = 0; i < N; ++i )
    if ( clusters[i].isleft )
      Ptot = add(Ptot, clusters[i].p = rotate(clusters[i].p, R));
  for ( int i = 0, No = other.size(); i < No; ++i )
    Ptot = add(Ptot, other[i] = rotate(other[i], R));
  shat = m2(Ptot);
  Pa = create(Unit(), Unit(), 0.5*poslc(Ptot), 0.5*poslc(Ptot));
  Pb = create(Unit(), Unit(), -0.5*neglc(Ptot), 0.5*neglc(Ptot));
  if ( !approx(Pa, Paold, sqrt(shat)) || !approx(Pb, Pbold, sqrt(shat)) ) {
    cerr << "**** Momentum non-conservation:\n     Pa (summed): ";
    print(Pa);
    cerr << "     Pa (boost) : ";
    print(Paold);
    cerr << "     Pb (summed): ";
    print(Pb);
    cerr << "     Pb (boost) : ";
    print(Pbold);
    cerr << "++++" << endl;
  }
  if ( c.isrdir >= 0 ) updatescachea();
  if ( c.isrdir <= 0 ) updatescacheb();
  return nrem;
}


template <typename LorentzMomentum>
bool InterKT::Clustering<LorentzMomentum>::
approx(LorentzMomentum p1, const LorentzMomentum & p2,
       Unit scale, double eps) {
  p1 = subtract(p1, p2);
  if ( abs(x(p1))/scale > eps ) return false;
  if ( abs(y(p1))/scale > eps ) return false;
  if ( abs(z(p1))/scale > eps ) return false;
  if ( abs(e(p1))/scale > eps ) return false;
  return true;
}


template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::
print(const LorentzMomentum & p) const {
  cerr << "("
       << roundGeV(x(p)) << ", "
       << roundGeV(y(p)) << ", "
       << roundGeV(z(p)) << "; "
       << roundGeV(e(p)) << "; "
       << roundGeV(m(p)) << ")\n                            ("
       << rap(p) << ", "
       << phi(p) << ", "
       << roundGeV(pt(p)) << ")" << endl;
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::
printdr(const LorentzMomentum & p, const LorentzMomentum & q) const {
  double deta = abs(rap(p) - rap(q));
  double dphi = abs(phi(p) - phi(q));
  if ( dphi > M_PI ) dphi = 2.0*M_PI - dphi;
  double dR = sqrt(deta*deta + dphi*dphi);
  cerr << dR << " (" << deta << ", " << dphi << ")" << endl;
}

template <typename LorentzMomentum>
double InterKT::Clustering<LorentzMomentum>::
deltar(const LorentzMomentum & p, const LorentzMomentum & q) const {
  double deta = abs(rap(p) - rap(q));
  double dphi = abs(phi(p) - phi(q));
  if ( dphi > M_PI ) dphi = 2.0*M_PI - dphi;
  return sqrt(deta*deta + dphi*dphi);
}

template <typename LorentzMomentum>
double InterKT::Clustering<LorentzMomentum>::roundGeV(Unit xx) {
  double x = GeV(xx);
  if ( abs(x) < 1.0e-6 ) return 0;
  return x;
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::debugbefore(Cluster & c, ofstream & test) const {
  if ( c.emitter && c.recoiler ){
  double zfac = 1.0 - m2(subtract(Ptot, c.p))/shat;
  Unit2 kt2i = -zfac*m2(subtract(Pa, c.p));
  //test << c.minkt2 << "\t";
  //test << kt2i/c.minkt2 << endl;
  kt2i = -zfac*m2(subtract(Pb, c.p));
  //test << kt2i << endl;
  //test << pt(c.p)/c.minkt2 << endl;
  Unit2 pt2i = Unit2();
  pt2i = invpt2(m2(add(c.p,c.emitter->p)),m2(add(c.p,Pa)),m2(add(c.emitter->p,Pa)));
  //test << pt2i/c.minkt2 << endl; 
  test << sqrt(pt2i*pt2(add(c.p,c.emitter->p)))/m2(add(c.p,c.emitter->p)) << endl;
  }

  if ( !debug ) return;
  cerr << ">>>> Start clustering" << endl;
  cerr << "     Ptot before:     ";
  print(Ptot);
  cerr << "     Pa before:       ";
  print(Pa);
  cerr << "     Pb before:       ";
  print(Pb);
  cerr << "     Cluster before:  ";
  print(c.p);
  if ( c.emitter ) {
    cerr << "     Emitter before:  ";
    print(c.emitter->p);
    cerr << "            Delta12:  ";
    printdr(c.emitter->p, c.p);
  }
  if ( c.recoiler ) {
    cerr << "     Recoiler before: ";
    print(c.recoiler->p);
    cerr << "            Delta13:  ";
    printdr(c.recoiler->p, c.p);
    cerr << "            Delta23:  ";
    printdr(c.emitter->p, c.recoiler->p);
  }
  cerr << "     Scale:           " << GeV(sqrt(c.minkt2)) << endl;

  if ( c.emitter && c.recoiler )
    cerr << "---- Final-state clustering" << endl;
  else if ( c.emitter && c.isrdir )
    cerr << "---- Final-initial clustering (" << c.isrdir << ")" << endl;
  else if ( c.emitter )
    cerr << "---- MI clustering" << endl;
  else
    cerr << "---- Initial-state clustering (" << c.isrdir << ")" << endl;
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::debugafter(Cluster & c) const {
  if ( !debug ) return;
  LorentzMomentum ptotref = missing;
  for ( int i = 0; i < N; ++i )
    if ( clusters[i].isleft ) ptotref = add(ptotref, clusters[i].p);
  for ( int i = 0, No = other.size(); i < No; ++i )
    ptotref = add(ptotref, other[i]);
  cerr << "     Ptot after:      ";
  print(Ptot);
  cerr << "     Ptot(ref) after: ";
  print(Ptot);
  cerr << "     Pa after:        ";
  print(Pa);
  cerr << "     Pb after:        ";
  print(Pb);
  if ( c.emitter && c.emitter->isleft ) {
    cerr << "     Emitter after:   ";
    print(c.emitter->p);
  }
  if ( c.recoiler && c.recoiler->isleft ) {
    cerr << "     Recoiler after:  ";
    print(c.recoiler->p);
  }
  cerr << "     " << nclus << " out of " << N
       << " clusters left." << endl;
  cerr << "<<<< Clustering finished" << endl;
}


template <typename LorentzMomentum>
std::vector<LorentzMomentum> InterKT::Clustering<LorentzMomentum>::getJets() { //////// std dodano
  vector<LorentzMomentum> ret;
  if ( jets.empty() ) filljets();
  for ( int i = 2, Ni = jets.size(); i < Ni; ++i )
    ret.push_back(jets[i].p);
  return ret;
}

template <typename LorentzMomentum>
void InterKT::Clustering<LorentzMomentum>::filljets() {
  jets.clear();
  jets.push_back(Jet(-1));
  jets.back().p = Pa;
  jets.push_back(Jet(-2));
  jets.back().p = Pb;
  for (int ic = 0; ic < N; ++ic ) addtojets(ic);
}

template <typename LorentzMomentum>
int InterKT::Clustering<LorentzMomentum>::addtojets(int ic) {
  if ( clusters[ic].injet >= 0 ) return clusters[ic].injet;
  if ( clusters[ic].incluster == ic ) {
    clusters[ic].injet = jets.size() - 1;
    jets.push_back(Jet());
    jets.back().p = clusters[ic].p;
    jets.back().cluster = ic;
    jets.back().tracks.insert(ic);
    return clusters[ic].injet;
  }
  if ( clusters[ic].incluster >= 0 ) {
    int j = addtojets(clusters[ic].incluster);
    jets[j].tracks.insert(ic);
    return j;
  } else {
    int j = -1 - clusters[ic].incluster;
    jets[j].tracks.insert(ic);
    return j;
  }
}

template <typename LorentzMomentum>
double InterKT::Clustering<LorentzMomentum>::elapsed() {
  static long last = 0;
  long el = clock() - last;
  last += el;
  return double(el)*0.000001;
}

template <typename LorentzMomentum>
int InterKT::Clustering<LorentzMomentum>::
cluster(const vector<LorentzMomentum> & tracksin,
	ofstream & test,
	const vector<LorentzMomentum> & otherin) {

  ofstream os;
  if ( sequence ) os.open("sequence.dat");
  static ofstream oss("/tmp/speed.dat");
  //test << "AAA" << endl;
  jets.clear();
  elapsed();
  tracks = tracksin;
  clusters.clear();
  N = tracks.size();
  hardtouch = true;
  
  Unit sumplus = Unit();
  Unit summinus = Unit();
  Unit sumx = Unit();
  Unit sumy = Unit();
  nclus = 0;
  for ( int i = 0; i < N; ++i ) {
    if ( pt(tracks[i]) >= cutkt ) {
      clusters.push_back(Cluster(tracks[i], nclus++, i));
      sumplus += poslc(tracks[i]);
      summinus += neglc(tracks[i]);
      sumx += x(tracks[i]);
      sumy += y(tracks[i]);
    }
  }
  N = nclus;
  other = otherin;
  for ( int i = 0, No = other.size(); i < No; ++i ) {
    sumplus += poslc(other[i]);
    summinus += neglc(other[i]);
    sumx += x(other[i]);
    sumy += y(other[i]);
  }
  Unit missingpt = sqrt(sumx*sumx + sumy*sumy);
  sumplus += missingpt;
  summinus += missingpt;
  missing = create(-sumx, -sumy, Unit(), missingpt);
  Pa = create(Unit(), Unit(), 0.5*sumplus, 0.5*sumplus);
  Pb = create(Unit(), Unit(), -0.5*summinus, 0.5*summinus);
  Ptot = add(Pa, Pb);
  shat = m2(Ptot);

  initscache();
  lastemitter = -1;
  lastrecoiler = -1;
  lastrem1 = -1;
  lastrem2 = -1;

  while ( nclus > nmin ) {
    updatescache();
    sorts();
    for ( int i = 0; i < N; ++i )
      if ( clusters[i].isleft ) update(i, nclus - nmin);

    if ( nclus == N && speedy ) {
      double el = elapsed();
      oss << N << "  " << el << " ";
      cerr << N << "  " << el << " ";
    }

    minkt2 = shat;
    int isel = -1;

    hardtouch = false;
    lastemitter = -1;
    lastrecoiler = -1;
    lastrem1 = -1;
    lastrem2 = -1;
    for ( int i = 0; i < N; ++i ) {
      if ( !clusters[i].isleft ) continue;
      clusters[i].istouched = false;

      if ( isel < 0 || clusters[i].minkt2 < minkt2 ) {//determining min_d_ij
	isel = i;
	minkt2 = clusters[i].minkt2;
      }

    }
     if ( sequence )
{ 
	os<<nclus<<endl;
	os << setprecision(4) << ( clusters[isel].emitter? clusters[isel].emitter->incluster: -1 )
    	 << " "
    	 << ( clusters[isel].recoiler? clusters[isel].recoiler->incluster: -1 )
    	 << " ";
    	// << clusters[isel].isrdir << endl;
          for ( int i = 0; i < N; ++i ) {	
	
	if(i==isel && clusters[i].emitter && clusters[i].recoiler)
	os << pt(clusters[i].p) <<"\t"<< rap(clusters[i].p) <<"\t"<< phi(clusters[i].p) <<"\t"<<"1"<<"\t";
	else if(i==isel && clusters[i].emitter && clusters[i].isrdir)
	os << pt(clusters[i].p) <<"\t"<< rap(clusters[i].p) <<"\t"<< phi(clusters[i].p) <<"\t"<<"2"<<"\t";
	else if(i==isel && clusters[i].emitter)
	os << pt(clusters[i].p) <<"\t"<< rap(clusters[i].p) <<"\t"<< phi(clusters[i].p) <<"\t"<<"3"<<"\t";
	else if(i==isel)
	os << pt(clusters[i].p) <<"\t"<< rap(clusters[i].p) <<"\t"<< phi(clusters[i].p) <<"\t"<<"4"<<"\t";
        else if(clusters[i].isleft)
	os << pt(clusters[i].p) <<"\t"<< rap(clusters[i].p) <<"\t"<< phi(clusters[i].p) <<"\t"<<"0"<<"\t";
        else
        os << "0.00000" <<"\t"<< "0.00000" <<"\t"<< "0.00000" <<"\t"<<"0"<<"\t";
	
    }
    os << endl;
}
   
    if ( minkt2 > maxkt2 && maxkt2 > Unit2() ) break;//checking cut off parameters

    nclus -= cluster(clusters[isel], test);//doing clustering step
    debugafter(clusters[isel]);
  }

  if ( speedy ) {
    double el = elapsed();
    oss << el << endl;
    cerr << el << endl;
  }

  if ( nclus < nmin || ( nmax >= 0 && nclus > nmax ) )
    return -1;//checking if we have reached wanted number of clusters
  
  if ( sequence )
{ 
	os<<nclus<<endl;
	os<<"-1"<<" "<<"-1"<<" ";
          for ( int i = 0; i < N; ++i ) {
        if(clusters[i].isleft)        
	os << pt(clusters[i].p) <<"\t"<< rap(clusters[i].p) <<"\t"<< phi(clusters[i].p) <<"\t"<<"0"<<"\t";
	else
	os << "0.00000" <<"\t"<< "0.00000" <<"\t"<< "0.00000" <<"\t"<<"0"<<"\t";
    }
    os << endl;
}
  
  return nclus;

}

