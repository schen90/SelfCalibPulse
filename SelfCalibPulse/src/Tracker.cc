#ifndef TRACKER_CC
#define TRACKER_CC

#include "Tracker.hh"

Tracker::Tracker(vector<Hit*>* hits){
  Tracker(hits,-1);
}

Tracker::Tracker(vector<Hit*>* hits, double E){
  Tracker(hits, E, TVector3(0,0,0));
}

Tracker::Tracker(vector<Hit*>* hits, double E, TVector3 sourcepos){

  fHits = hits;
  nhits = fHits->size();

  EGamma = E; //keV

  sPos.SetXYZ(sourcepos.X()/10., sourcepos.Y()/10., sourcepos.Z()/10.); // change to cm
  
  double tmppos[3];
  for(int i=0; i<nhits; i++){
    e[i] = fHits->at(i)->GetE() / 1000.; // change to MeV 
    Pos[i] = fHits->at(i)->GetPosition();
    tmppos[0] = Pos[i].X()/10. - sPos.X();
    tmppos[1] = Pos[i].Y()/10. - sPos.Y();
    tmppos[2] = Pos[i].Z()/10. - sPos.Z();
    Pos[i].SetXYZ(tmppos[0], tmppos[1], tmppos[2]); // cm and relative to source
  }

  // tracker parameters
  radius = 23.5008; //center to front surface in cm
  radius2 = 32.9202; //center to back surface in cm
  eres = 3e-3; // energy resolution MeV
  alfared = 1;
  deltaalfa = 0.05;
  kmax = 10;
  alfamin = 0.15;
  nopair = 0;
  minprobtrack = 0.05;
  sigma_thet = 0.8;
}

Tracker::~Tracker(){
}


bool Tracker::CheckOrder(){
  bool isCorrect = true;
  if(track.size()<2) return isCorrect;
  for(int i=0; i<track.size()-1; i++){
    if( fHits->at(track[i])->GetInterid() > fHits->at(track[i+1])->GetInterid() )
      isCorrect = false;
  }
  return isCorrect;
}


// OFT
double sig_pair(double E){ // pair cross section MeV and cm
  /* fitted in 1e-24cm2 units for Ge from atomic data tables 1970 7 page 590 */
  double temp;
  temp = 0.792189*log(E+0.948261-1.1332*E+0.15567*SQ(E));

  if(E<1.022)
    temp = 0;

  if((E < 1.15) & (E >= 1.022))
    temp=(1-((1.15-E)/0.129))*7.55e-28;

  return temp*1e-24;

}

double sig_compt(double E){ // MeV and cm
  /* sigma = 2 pi r_0^2 *[((1+gamma)/gamma^2 *(2(1+gamma)/1+2gamma - ln(1+2gamma)/gamma)) + ln(1+2gamma)/2gamma - 1+3gamma/(1+2gamma)^2] */
  // fits geant data very well
  
  double temp, temp0, temp1, temp2, temp3;
  double gamma;

  temp0 = 1e4*2*PI*SQ(r_0)*Z_ge; /* in cm2/atom */
  gamma = E/mec2;

  temp1 =1+gamma;
  temp2= 1 +(2*gamma);
  temp3= 1+(3*gamma);
  temp = (temp1/SQ(gamma)) * (((2*temp1)/temp2) - (log(temp2)/gamma));
  temp=temp +(log(temp2)/(2*gamma) - (temp3/SQ(temp2)));
  temp=temp*temp0;

  return temp ;  

}

double sig_abs(double E){ // MeV and cm

  /* sigma = 4 Alpha^4*sqrt(2)*Z^5 phi_0*(E/mec2)^(-7/2),phi_0=8/3 pi r_0^2 */
  /* sigma abs K+L = 9/8* sigma_Kshell */

  double temp;
  double gamma;
  double hnu_k;

  hnu_k = SQ(Z_ge - 0.03)*mec2*SQ(Alpha)/2;
  gamma = CB(E/mec2)*CB(E/mec2)*E/mec2;
  gamma = sqrt(gamma);

  temp = 4*CB(Alpha)*Alpha*1.4142*6.651e-25*CB(Z_ge)*SQ(Z_ge);

  temp = sqrt(E/mec2) * temp/gamma; /* en cm2/atom */

  // not well suited for energies below 20 keV
  // removed the 1.125 factor and added sqrt(E/mec2) to fit data

  if(E < 0.025) {

    temp = 2.2*pow((hnu_k/E),2.6666)*6.3e-18/SQ(Z_ge);
    if(E<0.0111)
      temp = temp/8.5;

  }
  return temp;

}

double range_process(double sig){ // cm
  /* SIGMA MACRO = 1/lambda = 6.022e23 * rho_ge/A_ge * sigma (cm2/atom) */
  double temp;

  temp = (sig * N_av * rho_ge) / A_ge;
  temp = 1/(temp); /* in cm */
  
  return temp ;
}


double proba(double range, double distance){
  double temp;
  double nlambda;

  nlambda = distance/range;
  temp = exp(-nlambda);

  return temp;
}


double err_cos(TVector3 &posa, TVector3 &posb, TVector3 &posc, double rab, double rbc, double err){
  /* error on AB.BC/AB x BC */
  double prod;
  double dcosdxa,dcosdxb,dcosdxc;
  double dcosdya,dcosdyb,dcosdyc;
  double dcosdza,dcosdzb,dcosdzc;

  double temp;

  prod = (posa-posb).Dot(posb-posc);

  double xa = posa.x();  double ya = posa.y();  double za = posa.z();
  double xb = posb.x();  double yb = posb.y();  double zb = posb.z();
  double xc = posc.x();  double yc = posc.y();  double zc = posc.z();

  dcosdxa = -(xc-xb)/(rab*rbc) -0.5*(2*xa-2*xb)*prod/(CB(rab)*rbc);
  dcosdxb = (xa-2*xb+xc)/(rab*rbc) -0.5*(2*xb-2*xa)*prod/(CB(rab)*rbc) +0.5*(2*xc-2*xb)*prod/(CB(rbc)*rab);
  dcosdxc = (xb-xa)/(rab*rbc) - 0.5*(2*xc-2*xb)*prod/(CB(rbc)*rab);

  dcosdya = -(yc-yb)/(rab*rbc) -0.5*(2*ya-2*yb)*prod/(CB(rab)*rbc);
  dcosdyb = (ya-2*yb+yc)/(rab*rbc) -0.5*(2*yb-2*ya)*prod/(CB(rab)*rbc) +0.5*(2*yc-2*yb)*prod/(CB(rbc)*rab);
  dcosdyc = (yb-ya)/(rab*rbc) - 0.5*(2*yc-2*yb)*prod/(CB(rbc)*rab);

  dcosdza = -(zc-zb)/(rab*rbc) -0.5*(2*za-2*zb)*prod/(CB(rab)*rbc);
  dcosdzb = (za-2*zb+zc)/(rab*rbc) -0.5*(2*zb-2*za)*prod/(CB(rab)*rbc) +0.5*(2*zc-2*zb)*prod/(CB(rbc)*rab);
  dcosdzc = (zb-za)/(rab*rbc) - 0.5*(2*zc-2*zb)*prod/(CB(rbc)*rab);

  temp = SQ(dcosdxa)+SQ(dcosdxb)+SQ(dcosdxc);

  temp=temp+SQ(dcosdya)+SQ(dcosdyb)+SQ(dcosdyc);

  temp=temp+SQ(dcosdza)+SQ(dcosdzb)+SQ(dcosdzc);

  temp=sqrt(temp);

  temp=temp*err;

  return temp;
}


// calculate Ftot with iteration for OFT tracking
double Tracker::Ftot(vector<TVector3> &pos, vector<int> &intid, vector<double> &energy,
		     double etotale, vector<int> &order){
  if(order.size()==0){
    order.push_back(0); // order[0]=0: source
    return Ftot(pos,intid,energy,etotale,order);
  }

  vector<int> bestorder;  bestorder.resize(order.size());
  vector<int> testorder;  testorder.resize(order.size());
  for(int i=0; i<order.size(); i++) bestorder[i] = order[i];
  double mincosdif = 0;

  if(order.size()==1){ // order[0]=0: source
    for(int m=1; m<energy.size(); m++){
      // next interaction m
      testorder.resize(order.size());
      for(int i=0; i<order.size(); i++) testorder[i] = order[i];
      testorder.push_back(m);

      double probtest = Ftot(pos,intid,energy,etotale,testorder);

      if(probtest>mincosdif){
        mincosdif = probtest;
        testorder.swap(bestorder);
      }
    }

  }else{ // order.size()>=2
    int j=order[order.size()-2], l=order[order.size()-1];
    double r01 = (pos[l]-pos[j]).Mag();
    double a = order.size()==2 ? 2 : 1;

    for(int m=1; m<energy.size(); m++){
      int skip = 0;
      for(int i=1; i<order.size(); i++) if(m==order[i]) skip=1;
      if(skip) continue;

      // next interaction m
      testorder.resize(order.size());
      for(int i=0; i<order.size(); i++) testorder[i] = order[i];
      testorder.push_back(m);

      double r12 = (pos[m]-pos[l]).Mag();
      double costheta = (pos[l]-pos[j]).Dot(pos[m]-pos[l])/r01/r12;
      double ercos = err_cos(pos[j],pos[l],pos[m],r01,r12,sigma_thet);

      double escattern = etotale/(1+(etotale/mec2)*(1-costheta));
      double escatter = etotale - energy[l];

      double deltaescn = SQ(escattern)*ercos/mec2;
      double deltaesc = sqrt((energy.size()+order.size()-2)*SQ(eres)); // ???
      
      double probtest = exp(-a*SQ(escattern-escatter)/(SQ(deltaescn)+SQ(deltaesc)));

      // testorder.size>=3
      if(testorder.size()< maxinteraction-1 && probtest < 0.005) continue; // 2~4 interactions
      if(testorder.size()==maxinteraction-1 && probtest < 0.0005) continue; // the 5th interaction

      if(order.size()==2){ // first interaction
        double cross1  = sig_compt(etotale);
        double cross   = cross1 + sig_pair(etotale) + sig_abs(etotale);
        double coef1   = cross1/cross;
        double lambda1 = range_process(cross1); // ???
        probtest *= SQ(coef1*proba(lambda1,r_ge[intid[l]][intid[l]]));
      }

      double cross2  = sig_compt(escatter);
      double cross   = cross2 + sig_pair(escatter) + sig_abs(escatter);
      if(testorder.size()==energy.size() || testorder.size()==maxinteraction+1) // last interaction
        cross2 = sig_abs(escatter);
      double coef2   = cross2/cross;
      double lambda2 = range_process(cross2); // ???
      probtest *= coef2*proba(lambda2,r_ge[intid[l]][intid[m]]);

      if(testorder.size()<energy.size() && testorder.size()<maxinteraction+1) // next interaction
        probtest *= Ftot(pos,intid,energy,escatter,testorder);

      if(probtest>mincosdif){
        mincosdif = probtest;
        testorder.swap(bestorder);
      }
    }
  }

  order.swap(bestorder);
  bestorder.clear();
  testorder.clear();
  return mincosdif;
}


/***********************************************************/
/**** OFT tracking *****************************************/
/**** energies in MeV and positions in cm ******************/
/***********************************************************/
void Tracker::OFTtracking(){
    
  //******** computing all distances/angle between interactions ********
  for(int i=0; i<nhits; i++){
    TVector3 Ai = Pos[i] + sPos;
    r[i][i] = Pos[i].Mag();

    /* distance from center to line joining points : rcenter = sintheta * distance center-i */
    double r_center = Ai.Cross(Pos[i]).Mag()/r[i][i];
    r_ge[i][i] = r[i][i];
    thetatest[i][i] = PI;

    /* determine effective length in Ge */
    if(sPos.Mag()<radius){ // source inside

      if(r_center <= radius) {
	// from origin to inner surface of shell
	double r_vacuum=sqrt(SQ(radius)-SQ(r_center));
	// distance from perpendicular intersction to i-source
	double r_vacuumi=sqrt(SQ(Ai.Mag())-SQ(r_center));
	//distance inside Ge from origin to i
	r_ge[i][i]=r_vacuumi-r_vacuum;
      }

    }else if(sPos.Mag()>radius2){ // source outside

      if(r_center <= radius2) {
	// from origin to inner surface of shell
	double r_vacuum=sqrt(SQ(radius2)-SQ(r_center));
	// distance from perpendicular intersction to i-source
	double r_vacuumi=sqrt(SQ(Ai.Mag())-SQ(r_center));
	//distance inside Ge from origin to i
	r_ge[i][i]=r_vacuum-r_vacuumi;
      }

      if((r_center < radius)  && (Ai.Dot(Ai-sPos)>0) && (sPos.Dot(sPos-Ai)>0)) {
	/* if < radius, the photon does not go through */
	/* Ge all the way from j to i */
	/* right triangle:  r_vacuum^2 + rcenter^2 = radius^2 */
	/* 2*r_vacuum = distance to take away from rmin */
	double r_vacuum = sqrt(SQ(radius)-SQ(r_center));
	r_ge[i][i] = r_ge[i][i] - 2*r_vacuum;
      }
    }

    if(r_ge[i][i] < 0) r_ge[i][i]=0.03;

    for(int j=i+1; j<nhits; j++){	
      TVector3 Aj = Pos[j] + sPos;
      r[i][j] = (Ai-Aj).Mag();
      r[j][i] = r[i][j];

      /* norm of vector product = sintheta * distance center-i * distance i-j */
      /* distance from center to line joining points : rcenter = sintheta * distance center-i */
      r_center = Ai.Cross(Ai-Aj).Mag()/r[i][j];
      r_ge[i][j] = r[i][j];
      thetatest[i][j] = Pos[i].Angle(Pos[j]);

      /* determine effective length in Ge */
      if((r_center < radius)  && (Ai.Dot(Ai-Aj)>0) && (Aj.Dot(Aj-Ai)>0)) {
	/* if < radius, the photon does not go through */
	/* Ge all the way from j to i */
	/* right triangle:  r_vacuum^2 + rcenter^2 = radius^2 */
	/* 2*r_vacuum = distance to take away from rmin */
	double r_vacuum = sqrt(SQ(radius)-SQ(r_center));
	r_ge[i][j] = r[i][j] - 2*r_vacuum;
      }

      if(r_ge[i][j] < 0) r_ge[i][j]=0.02;

      r_ge[j][i] = r_ge[j][j];
      thetatest[j][i] = thetatest[i][j];
    }
  }
    
  //******** find cluster ********
  vector<double> et;
  vector<vector<int>> sn(MaxNDets);
  int flagu[MaxNDets];

#ifdef ONECLUST // put all hits in one cluster
  et.push_back(0);
  sn[0].clear();
  for(int i=0; i<nhits; i++){
    sn[0].push_back(i);
    et[0] += e[i];
  }
  int n = et.size();

#else // make cluster by angle
  double power = pow((nhits+2)/3., 0.9);
  double alfamax = acos(1-2./power)/alfared;

  int n = et.size(); // counter of cluster
  for(double alfa=alfamin; alfa<alfamax; alfa+=deltaalfa){
    for(int i=0; i<nhits; i++) flagu[i]=0;

    for(int i=0; i<nhits; i++){
      if(flagu[i]!=0) continue;

      sn[n].clear();
      sn[n].push_back(i);
      et.push_back(e[i]); // total energy of cluster
      flagu[i]=1;

      for(int j=0; j<nhits; j++){
	if(j==i || flagu[j]!=0) continue;

	if(thetatest[i][j]<=alfa){
	  sn[n].push_back(j);
	  et[n] += e[j];
	  flagu[j]=1;
	}

	if(sn[n].size() == kmax) break;
      }

      // accept new cluster only if unique
      for(int nn=0; nn<n; nn++){
	if(sn[nn].size()==sn[n].size() && fabs(et[nn]-et[n])<1e-6){
	  et.pop_back();
	  break;
	}
      }

      n=et.size();
    }
  }// end of loop alfa
#endif // end of ONECLUST
  
  //******** compute figure of merit of clusters ********
  int flagpair[MaxNDets];
  vector<vector<int>> interaction(MaxNDets);
  vector<double> probtot;
  for(int i=0; i<n; i++){ //loop clusters
    flagu[i]=0;
    flagpair[i]=0;
    interaction[i].clear();

    if(sn[i].size()==1){
      probtot.push_back(minprobtrack);
      interaction[i].push_back(sn[i][0]);

    }else{ //sn[i].size()>1
      probtot.push_back(0);

      // pair production -----------------
      if(nopair==0){
	for(int j=0; j<sn[i].size(); j++){

	  if( (et[i]-1.022) >= -(sqrt(sn[i].size())*eres) &&
	      fabs(et[i]-e[sn[i][j]]-1.022) <= (sqrt((sn[i].size()+1))*eres)){

	    flagpair[i]=1;
	    interactionpair=sn[i][j];
	    double cross= sig_pair(et[i]);
	    double cross1 = sig_compt(et[i]);
	    double cross2 = sig_abs(et[i]);
	    double coef1 = cross/(cross+cross1+cross2);
	    double lambda=range_process(cross);
	    probtestpair=coef1*proba(lambda,r_ge[sn[i][j]][sn[i][j]]);
	  }
	}
      }

      // calculate Ftot with iteration
      //vector <TVector3> pos;  pos.push_back(sPos);
      vector <TVector3> pos;  pos.push_back(TVector3(0,0,0));
      vector <int> intid;     intid.push_back(-1);
      vector <double> energy; energy.push_back(0);
      vector <int> order;     order.push_back(0);
      for(int j=0; j<sn[i].size(); j++){
	pos.push_back(Pos[sn[i][j]]); // Pos relative to source
	intid.push_back(sn[i][j]);
	energy.push_back(e[sn[i][j]]);
      }

      double mincosdif = Ftot(pos, intid, energy, et[i], order);
      probtot[i] = pow(mincosdif,1./(2*order.size()-3));
      if(order.size()==sn[i].size()+1 || order.size()==maxinteraction+1)
	for(int j=1; j<order.size(); j++) interaction[i].push_back(intid[order[j]]);
    }

    // pair production event and not good compton
    if(flagpair[i]==1 && probtestpair>probtot[i]){
      if(interaction[i].size()<1) interaction[i].resize(1);
      interaction[i][0]=interactionpair;
      probtot[i]=probtestpair;
    }
    
  }// end of loop clusters

  
  //******** sort clusters according to figure of merit ********

  // single interactions are awarded for the time being the minimum figure of merit
  for(int i=0; i<n; i++)
    for(int j=i+1; j<n; j++)
      if(probtot[i] < probtot[j]){
	swap(probtot, i, j);
	swap(et, i, j);
	swap(flagu, i, j);
	swap(flagpair, i, j);
	interaction[i].swap(interaction[j]);
	sn[i].swap(sn[j]);
      }

  for(int i=0; i<n; i++)
    if(flagu[i]==0)
      for(int k=0; k<sn[i].size(); k++)
	for(int l=i+1; l<n; l++)
	  if(flagu[l]==0)
	    for(int m=0; m<sn[l].size(); m++)
	      if(sn[i][k] == sn[l][m]){
		flagu[l]=1; // remove cluster with used interactions
		break;
	      }

  // save the most likely track
  if(interaction.size()<1) return;
  
  for(int i=0; i<interaction[0].size(); i++) track.push_back(interaction[0][i]);

  return;
}



/***********************************************************/
/**** Simple tracking **************************************/
/**** energies in MeV and positions in cm ******************/
/***********************************************************/

Double_t calcComptonAngle(Double_t Ein, Double_t Edep){
  Double_t cosa = 1. + mec2 / Ein - mec2 / (Ein - Edep);
  if(unlikely( cosa > 1 ))
    cosa = 1;
  else if(unlikely( cosa < -1 ))
    cosa = -1;

  return acos(cosa);
}


// calculate best min chi2 with iteration for simple tracking
void Tracker::CalcMinChi2(vector<TVector3> &pos, vector<int> &intid, vector<double> &energy,
			  double etotale, vector<int> &order, double bestchi2[]){
  if(order.size()==0){
    order.push_back(0); // order[0]=0: source
    return CalcMinChi2(pos,intid,energy,etotale,order, bestchi2);
  }

  vector<int> bestorder;  bestorder.resize(order.size());
  vector<int> testorder;  testorder.resize(order.size());
  for(int i=0; i<order.size(); i++) bestorder[i] = order[i];
  double minchi2[2] = {DBL_MAX,DBL_MAX};

  if(order.size()==1){ // order[0]=0: source
    for(int m=1; m<energy.size(); m++){
      // next interaction m
      testorder.resize(order.size());
      for(int i=0; i<order.size(); i++) testorder[i] = order[i];
      testorder.push_back(m);

      double chi2test[2];
      CalcMinChi2(pos,intid,energy,etotale,testorder, chi2test);

      if(chi2test[0]<minchi2[0]){
	minchi2[1] = minchi2[0];
        minchi2[0] = chi2test[0];
        testorder.swap(bestorder);
      }else if(chi2test[0]<minchi2[1]){
	minchi2[1] = chi2test[0];
      }

      if(chi2test[1]<minchi2[1]){
	minchi2[1] = chi2test[1];
      }
    }

  }else{ // order.size()>=2
    int j=order[order.size()-2], l=order[order.size()-1];
    TVector3 mom01 = pos[l]-pos[j];
    
    for(int m=1; m<energy.size(); m++){
      int skip = 0;
      for(int i=1; i<order.size(); i++) if(m==order[i]) skip=1;
      if(skip) continue;

      // next interaction m
      testorder.resize(order.size());
      for(int i=0; i<order.size(); i++) testorder[i] = order[i];
      testorder.push_back(m);

      TVector3 mom12 = pos[m]-pos[l];
      double angle = mom01.Angle(mom12);

      double cangle = calcComptonAngle(etotale, energy[l]);
      double escatter = etotale - energy[l];

      // doi:10.1016/j.nima.2004.06.154
      double chi2tmp = exp(40. * fabs(cos(cangle) - cos(angle)) );

      double chi2test[2] = {0,0};
      if(testorder.size()<energy.size()) // next interaction
        CalcMinChi2(pos,intid,energy,escatter,testorder, chi2test);

      chi2test[0] += chi2tmp;
      chi2test[1] += chi2tmp;
     
      if(chi2test[0]<minchi2[0]){
	minchi2[1] = minchi2[0];
        minchi2[0] = chi2test[0];
        testorder.swap(bestorder);
      }else if(chi2test[0]<minchi2[1]){
	minchi2[1] = chi2test[0];
      }

      if(chi2test[1]<minchi2[1]){
	minchi2[1] = chi2test[1];
      }      
    }
  }

  order.swap(bestorder);
  bestorder.clear();
  testorder.clear();
}


void Tracker::Simpletracking(){
  if(nhits<2) return;

  //******** computing all angle between interactions ********
  for(int i=0; i<nhits; i++){
    thetatest[i][i] = PI;

    for(int j=i+1; j<nhits; j++){	
      thetatest[i][j] = Pos[i].Angle(Pos[j]);
      thetatest[j][i] = thetatest[i][j];
    }
  }
    
  //******** find cluster ********
  vector<double> et;
  vector<vector<int>> sn(MaxNDets);
  int flagu[MaxNDets];

#ifdef ONECLUST // put all hits in one cluster
  et.push_back(0);
  sn[0].clear();
  for(int i=0; i<nhits; i++){
    sn[0].push_back(i);
    et[0] += e[i];
  }
  int n = et.size();

#else // make cluster by angle
  double power = pow((nhits+2)/3., 0.9);
  double alfamax = acos(1-2./power)/alfared;

  int n = et.size(); // counter of cluster
  for(double alfa=alfamin; alfa<alfamax; alfa+=deltaalfa){
    for(int i=0; i<nhits; i++) flagu[i]=0;

    for(int i=0; i<nhits; i++){
      if(flagu[i]!=0) continue;

      sn[n].clear();
      sn[n].push_back(i);
      et.push_back(e[i]); // total energy of cluster
      flagu[i]=1;

      for(int j=0; j<nhits; j++){
	if(j==i || flagu[j]!=0) continue;

	if(thetatest[i][j]<=alfa){
	  sn[n].push_back(j);
	  et[n] += e[j];
	  flagu[j]=1;
	}

	if(sn[n].size() == kmax) break;
      }

      // accept new cluster only if unique
      for(int nn=0; nn<n; nn++){
	if(sn[nn].size()==sn[n].size() && fabs(et[nn]-et[n])<1e-6){
	  et.pop_back();
	  break;
	}
      }

      n=et.size();
    }
  }// end of loop alfa
#endif // end of ONECLUST
  
  //******** compute min chi2 of clusters ********
  vector<vector<int>> interaction(MaxNDets);
  vector<double> chi2tot;
  vector<double> chi2tot2; // secondbestchi2
  for(int i=0; i<n; i++){ //loop clusters
    flagu[i]=0;
    interaction[i].clear();

    if(sn[i].size()==1){
      chi2tot.push_back(DBL_MAX);
      chi2tot2.push_back(DBL_MAX);
      interaction[i].push_back(sn[i][0]);

    }else{ //sn[i].size()>1
      chi2tot.push_back(DBL_MAX);
      chi2tot2.push_back(DBL_MAX);

      // calculate min chi2 with iteration
      vector <TVector3> pos;  pos.push_back(TVector3(0,0,0));
      vector <int> intid;     intid.push_back(-1);
      vector <double> energy; energy.push_back(0);
      vector <int> order;     order.push_back(0);
      for(int j=0; j<sn[i].size(); j++){
	pos.push_back(Pos[sn[i][j]]); // Pos relative to source
	intid.push_back(sn[i][j]);
	energy.push_back(e[sn[i][j]]);
      }
  
      double minchi2[2];
      CalcMinChi2(pos, intid, energy, et[i], order, minchi2);
      chi2tot[i] = minchi2[0] / (order.size()-2);
      chi2tot2[i] = minchi2[1] / (order.size()-2);

      if(order.size()==sn[i].size()+1)
	for(int j=1; j<order.size(); j++) interaction[i].push_back(intid[order[j]]);
    }
    
  }// end of loop clusters

  
  //******** sort clusters according to chi2 ********

  // single interactions are awarded for the time being the max chi2
  for(int i=0; i<n; i++)
    for(int j=i+1; j<n; j++)
      if(chi2tot[i] > chi2tot[j]){
	swap(chi2tot, i, j);
	swap(chi2tot2, i, j);
	swap(et, i, j);
	swap(flagu, i, j);
	interaction[i].swap(interaction[j]);
	sn[i].swap(sn[j]);
      }

  for(int i=0; i<n; i++)
    if(flagu[i]==0)
      for(int k=0; k<sn[i].size(); k++)
	for(int l=i+1; l<n; l++)
	  if(flagu[l]==0)
	    for(int m=0; m<sn[l].size(); m++)
	      if(sn[i][k] == sn[l][m]){
		flagu[l]=1; // remove cluster with used interactions
		break;
	      }

  // save the most likely track
  if(interaction.size()<1) return;
  
  for(int i=0; i<interaction[0].size(); i++) track.push_back(interaction[0][i]);

  FOM1 = chi2tot[0];
  FOM2 = chi2tot2[0];
  
  return;
}

#endif
