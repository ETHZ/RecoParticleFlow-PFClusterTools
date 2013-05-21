#include "RecoParticleFlow/PFClusterTools/interface/LinkByDetId.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/Math/interface/normalizedPhi.h"
#include "TMath.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

// for debugging (AA)
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

// to enable debugs
//#define PFLOW_DEBUG

#define DEBUG_FOR_AOD

double 
LinkByDetId::testTrackAndClusterByDetId ( const reco::PFRecTrack& track, 
					    const reco::PFCluster&  cluster,
						 const CaloGeometry* caloGeometry,
					    bool isBrem,
					    bool debug)  {
 
// want to be able to debug just this piece of code (AA)
#ifdef PFLOW_DEBUG
debug = true;
#endif

#ifdef PFLOW_DEBUG
  if( debug ) 
    std::cout << "entering  LinkByDetId::testTrackAndClusterByDetId" << std::endl;
#endif
  
  //cluster position
  double clustereta  = cluster.positionREP().Eta();
  double clusterphi  = cluster.positionREP().Phi();
  //double clusterX    = cluster.position().X();
  //double clusterY    = cluster.position().Y();
  double clusterZ    = cluster.position().Z();

  bool barrel = false;
  bool hcal = false;
  // double distance = 999999.9;
  double horesolscale=1.0;

  //track extrapolation
  const reco::PFTrajectoryPoint& atVertex 
    = track.extrapolatedPoint( reco::PFTrajectoryPoint::ClosestApproach );
  const reco::PFTrajectoryPoint& atECAL 
    = track.extrapolatedPoint( reco::PFTrajectoryPoint::ECALShowerMax );

  //track at calo's
  double tracketa = 999999.9;
  double trackphi = 999999.9;
  double track_X  = 999999.9;
  double track_Y  = 999999.9;
  double track_Z  = 999999.9;
  double dHEta = 0.;
  double dHPhi = 0.;

  // Quantities at vertex
  double trackPt = isBrem ? 999. : sqrt(atVertex.momentum().Vect().Perp2());
  // double trackEta = isBrem ? 999. : atVertex.momentum().Vect().Eta();

  switch (cluster.layer()) {
  case PFLayer::ECAL_BARREL: barrel = true;
  case PFLayer::ECAL_ENDCAP:
#ifdef PFLOW_DEBUG
    if( debug )
      std::cout << "Fetching Ecal Resolution Maps"
	   << std::endl;
#endif
    // did not reach ecal, cannot be associated with a cluster.
    if( ! atECAL.isValid() ) return -1.;   
    
    tracketa = atECAL.positionREP().Eta();
    trackphi = atECAL.positionREP().Phi();
    track_X  = atECAL.position().X();
    track_Y  = atECAL.position().Y();
    track_Z  = atECAL.position().Z();

/*    distance 
      = std::sqrt( (track_X-clusterX)*(track_X-clusterX)
		  +(track_Y-clusterY)*(track_Y-clusterY)
		  +(track_Z-clusterZ)*(track_Z-clusterZ)
		   );
*/			           
    break;
   
  case PFLayer::HCAL_BARREL1: barrel = true; 
  case PFLayer::HCAL_ENDCAP:  
#ifdef PFLOW_DEBUG
    if( debug )
      std::cout << "Fetching Hcal Resolution Maps"
	   << std::endl;
#endif
    if( isBrem ) {  
      return  -1.;
    } else { 
      hcal=true;
      const reco::PFTrajectoryPoint& atHCAL 
	= track.extrapolatedPoint( reco::PFTrajectoryPoint::HCALEntrance );
      const reco::PFTrajectoryPoint& atHCALExit 
	= track.extrapolatedPoint( reco::PFTrajectoryPoint::HCALExit );
      // did not reach hcal, cannot be associated with a cluster.
      if( ! atHCAL.isValid() ) return -1.;   
      
      // The link is computed between 0 and ~1 interaction length in HCAL
      dHEta = atHCALExit.positionREP().Eta()-atHCAL.positionREP().Eta();
      dHPhi = atHCALExit.positionREP().Phi()-atHCAL.positionREP().Phi(); 
      if ( dHPhi > M_PI ) dHPhi = dHPhi - 2.*M_PI;
      else if ( dHPhi < -M_PI ) dHPhi = dHPhi + 2.*M_PI; 
      tracketa = atHCAL.positionREP().Eta() + 0.1*dHEta;
      trackphi = atHCAL.positionREP().Phi() + 0.1*dHPhi;
      track_X  = atHCAL.position().X();
      track_Y  = atHCAL.position().Y();
      track_Z  = atHCAL.position().Z();
/*      distance 
	= -std::sqrt( (track_X-clusterX)*(track_X-clusterX)
		     +(track_Y-clusterY)*(track_Y-clusterY)
		     +(track_Z-clusterZ)*(track_Z-clusterZ)
		     );
*/			           
    }
    break;

  case PFLayer::HCAL_BARREL2: barrel = true; 
#ifdef PFLOW_DEBUG
    if( debug )
     std::cout << "Fetching HO Resolution Maps"
	       << std::endl;
#endif
    if( isBrem ) {  
      return  -1.;
    } else { 
      hcal=true;
      horesolscale=1.15;
      const reco::PFTrajectoryPoint& atHO 
	= track.extrapolatedPoint( reco::PFTrajectoryPoint::HOLayer );
      // did not reach ho, cannot be associated with a cluster.
      if( ! atHO.isValid() ) return -1.;   
      
      //
      tracketa = atHO.positionREP().Eta();
      trackphi = atHO.positionREP().Phi();
      track_X  = atHO.position().X();
      track_Y  = atHO.position().Y();
      track_Z  = atHO.position().Z();

      // Is this check really useful ?
      if ( fabs(track_Z) > 700.25 ) return -1.;


/*      distance 
	= std::sqrt( (track_X-clusterX)*(track_X-clusterX)
		      +(track_Y-clusterY)*(track_Y-clusterY)
		      +(track_Z-clusterZ)*(track_Z-clusterZ)
		      );
*/
#ifdef PFLOW_DEBUG
/*      if( debug ) {
	std::cout <<"dist "<<distance<<" "<<cluster.energy()<<" "<<cluster.layer()<<" "
		  <<track_X<<" "<<clusterX<<" "
		  <<track_Y<<" "<<clusterY<<" "
		  <<track_Z<<" "<<clusterZ<<std::endl;
      }
*/
#endif
    }
    break;

  case PFLayer::PS1:
  case PFLayer::PS2:
    //Note Alex: Nothing implemented for the
    //PreShower (No resolution maps yet)
#ifdef PFLOW_DEBUG
    if( debug )
      std::cout << "No link by DetId possible for pre-shower yet!"
	   << std::endl;
#endif
    return -1.;
  default:
    return -1.;
  }


  // Check that, if the cluster is in the endcap, 
  // 0) the track indeed points to the endcap at vertex (DISABLED)
  // 1) the track extrapolation is in the endcap too !
  // 2) the track is in the same end-cap !
  // PJ - 10-May-09
  if ( !barrel ) { 
    // if ( fabs(trackEta) < 1.0 ) return -1; 
    if ( !hcal && fabs(track_Z) < 300. ) return -1.;
    if ( track_Z * clusterZ < 0. ) return -1.;
  }
  // Check that, if the cluster is in the barrel, 
  // 1) the track is in the barrel too !
  if ( barrel ) {
    if ( !hcal && fabs(track_Z) > 300. ) return -1.;
  }

  double dist = LinkByDetId::computeDist( clustereta, clusterphi, 
					 tracketa, trackphi);
  
#ifdef PFLOW_DEBUG
  if(debug) std::cout<<"test link by DetId "<< dist <<" "<<std::endl;
  if(debug){
    std::cout<<" clustereta "  << clustereta 
	<<" clusterphi "  << clusterphi 
	<<" tracketa " << tracketa
	<<" trackphi " << trackphi << std::endl;
  }
#endif
  
  //Testing if Track can be linked by DetId to a cluster.
  //A cluster can be linked to a track if the extrapolated position 
  //of the track to the ECAL ShowerMax/HCAL entrance falls within 
  //the boundaries of any cell that belongs to this cluster.

 // const std::vector< reco::PFRecHitFraction >& 
 //   fracs = cluster.recHitFractions();

  // compare the lengths of the vectors storing the detIds and the pfRecHits (AA)
#ifdef DEBUG_FOR_AOD
  if (cluster.recHitFractions().size() != cluster.hitsAndFractions().size()) 
	  std::cout << "\nSize of RecHits (DetIds) " 
	  << cluster.recHitFractions().size() << "  " 
	  << cluster.hitsAndFractions().size() << std::endl;
#endif

  // for debugging (AA)
  std::vector< reco::PFRecHitFraction >::const_iterator  it_fracs = cluster.recHitFractions().begin();


 // first element of the pair is DetId, the second is the fraction of energy	          
  const std::vector< std::pair<DetId, float> >& hitsAndFracs = cluster.hitsAndFractions();
  bool linkedbyDetId = false;
  //loop over hitsAndFractions stored in the cluster
  std::vector< std::pair<DetId, float> >::const_iterator it_haf = hitsAndFracs.begin();
  for(; it_haf != hitsAndFracs.end(); ++it_haf){

	  // for debuggin (AA)
	     //getting rechit center position
///////////////////////    const reco::PFRecHit& rechit_cluster = *(it_fracs->recHitRef());
//    const math::XYZPoint& posxyzRH 
//      = rechit_cluster.position();
//////////////////////    const reco::PFRecHit::REPPoint& posrepRH  = rechit_cluster.positionREP();
    
    //getting rechit corners
//    const std::vector< math::XYZPoint >& 
//      cornersxyzRH = rechit_cluster.getCornersXYZ();
///////////////////    const std::vector<reco::PFRecHit::REPPoint>& cornersRH = rechit_cluster.getCornersREP();

	 // increment for the next entry (AA)
	 it_fracs++;

	 // -------------------------------------------------|



    if (it_haf->second < 1E-4) continue;
	const DetId cellId = it_haf->first;
    
	// first get the cell center position and corners 
	const CaloCellGeometry* cellGeometry = caloGeometry->getGeometry(cellId);
	
	// get the position of the cell
	const GlobalPoint gpPosition = cellGeometry->getPosition();
	// convert to types used in the PF framework
    const math::XYZPoint posxyz(gpPosition.x(), gpPosition.y(), gpPosition.z()); 
    reco::PFRecHit::REPPoint posrep(gpPosition.perp(), gpPosition.eta(), gpPosition.phi());
    
    //get the cell corners
	CaloCellGeometry::CornersVec cellGeoCorners = cellGeometry->getCorners();
	// convert to types used in the PF framework
	// also note that corner indexing of the cell face in PFRecHit is the reverse of the one returned 
	// by CaloCellGeometry!
	// these are just the front corners
    std::vector< math::XYZPoint > cornersxyz;
	std::vector<reco::PFRecHit::REPPoint> corners;
	for (uint i=0; i<4; ++i) {
		cornersxyz.push_back(math::XYZPoint(cellGeoCorners[3-i].x(),
											cellGeoCorners[3-i].y(),
											cellGeoCorners[3-i].z()
											));
		corners.push_back(reco::PFRecHit::REPPoint(cellGeoCorners[3-i].perp(),
													  cellGeoCorners[3-i].eta(),
													  cellGeoCorners[3-i].phi()
													  ));
	} // filling the corners
    
    if( barrel || hcal ){ // barrel case matching in eta/phi 
                          // (and HCAL endcap too!)
      
      //rechit size determination 
      // blown up by 50% (HCAL) to 100% (ECAL) to include cracks & gaps
      // also blown up to account for multiple scattering at low pt.

		 // here we do the manual corner reassignment for iEta=25,26,29 to match the 
		 // printput from pfrechits (AA)
		 // this is for debugging!!!
		 // need to fix in geometry if this is the sole cause of the discrepancey!!!!

		 
		 //		 if (cellId.det()==DetId::Hcal && HcalDetId(cellId).subdet()!=HcalOuter) {
		 if (cellId.det()!=DetId::Ecal) {
			 if (HcalDetId(cellId.rawId()).subdet()==HcalBarrel || HcalDetId(cellId.rawId()).subdet()==HcalEndcap) {

				 CaloTowerDetId ctDetId(cellId);

				 if (ctDetId.ietaAbs()==25) {

					 if (ctDetId.zside()==1) {
						 //std::cout << "\n\nreassigning eta\n";
						 //std::cout << "old eta " << corners[0].Eta();
						 corners[2].SetEta(2.322);
						 corners[3].SetEta(2.322);
					 }
					 //else {
						// corners[0].SetEta(-2.322);
						// corners[1].SetEta(-2.322);
					 //}
					 posrep.SetEta(double (2.247*ctDetId.zside()) );
					 //std::cout << "  new eta " << corners[0].Eta() << std::endl << std::endl;

				 }
				 else if (ctDetId.ietaAbs()==26) {
					 if (ctDetId.zside()==1) {
					 corners[0].SetEta(2.322);
					 corners[1].SetEta(2.322);
					 }
					 //else {
					 //corners[2].SetEta(-2.322);
					 //corners[3].SetEta(-2.322);
					 //}
					 posrep.SetEta(double (2.411*ctDetId.zside()) );
				 }
				 /*
				 else if (ctDetId.ietaAbs()==28) {
					 if (ctDetId.zside()==1) {
						 corners[2].SetEta(3.0);
						 corners[3].SetEta(3.0);
					 }
					 else {
						 corners[0].SetEta(-3.0);
						 corners[1].SetEta(-3.0);
					 }
					 posrep.SetEta(double (2.825*ctDetId.zside()) );
				 }
				 else if (ctDetId.ietaAbs()==29) {
					 if (ctDetId.zside()==1) {
						 corners[0].SetEta(2.65);
						 corners[1].SetEta(2.65);
					 }
					 else {
						 corners[2].SetEta(-2.65);
						 corners[3].SetEta(-2.65);
					 }
					 posrep.SetEta(double (2.825*ctDetId.zside()) );
				 }
				 */
			 }
		 }
//--------------------------------------------------------------



      double rhsizeEta 
	= fabs(corners[0].Eta() - corners[2].Eta());
      double rhsizePhi 
	= fabs(corners[0].Phi() - corners[2].Phi());
      if ( rhsizePhi > M_PI ) rhsizePhi = 2.*M_PI - rhsizePhi;
      if ( hcal ) { 
	rhsizeEta = rhsizeEta * horesolscale * (1.50 + 0.5/hitsAndFracs.size()) + 0.2*fabs(dHEta);
	rhsizePhi = rhsizePhi * horesolscale * (1.50 + 0.5/hitsAndFracs.size()) + 0.2*fabs(dHPhi); 
	
      } else { 
	rhsizeEta *= 2.00 + 1.0/hitsAndFracs.size()/std::min(1.,trackPt/2.);
	rhsizePhi *= 2.00 + 1.0/hitsAndFracs.size()/std::min(1.,trackPt/2.); 
      }

/*
		// debug (AA) -------------------------------------------------------
		double rhsizeEtaRH 
			= fabs(cornersRH[0].Eta() - cornersRH[2].Eta());
		double rhsizePhiRH 
			= fabs(cornersRH[0].Phi() - cornersRH[2].Phi());
		if ( rhsizePhiRH > M_PI ) rhsizePhiRH = 2.*M_PI - rhsizePhiRH;
		if ( hcal ) { 
			rhsizeEtaRH = rhsizeEtaRH * horesolscale * (1.50 + 0.5/hitsAndFracs.size()) + 0.2*fabs(dHEta);
			rhsizePhiRH = rhsizePhiRH * horesolscale * (1.50 + 0.5/hitsAndFracs.size()) + 0.2*fabs(dHPhi); 

		} else { 
			rhsizeEtaRH *= 2.00 + 1.0/hitsAndFracs.size()/std::min(1.,trackPt/2.);
			rhsizePhiRH *= 2.00 + 1.0/hitsAndFracs.size()/std::min(1.,trackPt/2.); 
		}

		bool foundDifference = false;

		for (int ii=0; ii<4; ++ii) {
			double ddPhi = corners[ii].Phi()-cornersRH[ii].Phi();
			while (ddPhi > TMath::Pi()) ddPhi -= (2.0*TMath::Pi());
			while (ddPhi < -TMath::Pi()) ddPhi += (2.0*TMath::Pi());

			if (fabs(ddPhi)>0.0001 || fabs(corners[ii].Eta()-cornersRH[ii].Eta())>0.0001) {
				foundDifference = true;
				break;
			}
		}

		if (foundDifference) {
			std::cout << std::endl;
			std::cout << "iEta = " << CaloTowerDetId(cellId).ieta() << std::endl;
			std::cout << it_haf - hitsAndFracs.begin()  << " Hcal DetId  eta phi " 
				<< posrep.Eta() << " " 
				<< posrep.Phi() << " " 
				<< std::endl;

			std::cout << it_fracs - cluster.recHitFractions().begin() - 1  << " Hcal RecHit eta phi " 
				<< posrepRH.Eta() << " " DQMOffline       RecoBTau            RecoParticleFlow  TrackingTools

				<< posrepRH.Phi() << " " 
				<< std::endl;

			std::cout << "DetId  SizeEta=" << rhsizeEta
				<< " SizePhi=" << rhsizePhi << std::endl;
			std::cout << "RecHit SizeEta=" << rhsizeEtaRH
				<< " SizePhi=" << rhsizePhiRH << std::endl;

			for ( unsigned jc=0; jc<4; ++jc ) {
				std::cout<<"corners detID "<<jc<<" "<<corners[jc].Eta()
				<<" "<<corners[jc].Phi()<<std::endl;
			std::cout<<"corners recHit"<<jc<<" "<<cornersRH[jc].Eta()
				<<" "<<cornersRH[jc].Phi()<<std::endl;
			}

		}
		//---------------------------------------------
		*/

//	 } // if barrel or HCAL


#ifdef PFLOW_DEBUG
      if( debug ) {
	std::cout << it_haf - hitsAndFracs.begin()  << " Hcal RecHit=" 
	     << posrep.Eta() << " " 
	     << posrep.Phi() << " " 
	     << std::endl; 
	for ( unsigned jc=0; jc<4; ++jc ) 
	  std::cout<<"corners "<<jc<<" "<<corners[jc].Eta()
	      <<" "<<corners[jc].Phi()<<std::endl;
	
	std::cout << "RecHit SizeEta=" << rhsizeEta
	     << " SizePhi=" << rhsizePhi << std::endl;
      }
#endif
      
      //distance track-rechit center
      // const math::XYZPoint& posxyz 
      // = rechit_cluster.position();
      double deta = fabs(posrep.Eta() - tracketa);
      double dphi = fabs(posrep.Phi() - trackphi);
      if ( dphi > M_PI ) dphi = 2.*M_PI - dphi;
      
#ifdef PFLOW_DEBUG
      if( debug ){
	std::cout << "distance=" 
	     << deta << " " 
	     << dphi << " ";
	if(deta < (rhsizeEta/2.) && dphi < (rhsizePhi/2.))
	  std::cout << " link here !" << std::endl;
	else std::cout << std::endl;
      }
#endif
      
      if(deta < (rhsizeEta/2.) && dphi < (rhsizePhi/2.)){ 
	linkedbyDetId = true;
	break;
      }
    }
    else { //ECAL & PS endcap case, matching in X,Y
      
#ifdef PFLOW_DEBUG
      if( debug ){

			// posxyz was already filled?? (AA)
//	const math::XYZPoint& posxyz 
//	  = rechit_cluster.position();
	
	std::cout << "RH " << posxyz.X()
	     << " "   << posxyz.Y()
	     << std::endl;
	
	std::cout << "TRACK " << track_X
	     << " "      << track_Y
	     << std::endl;
      }
#endif
      
      double x[5];
      double y[5];
      
      for ( unsigned jc=0; jc<4; ++jc ) {
	math::XYZPoint cornerposxyz = cornersxyz[jc];
	x[jc] = cornerposxyz.X() + (cornerposxyz.X()-posxyz.X())
	  * (1.00+0.50/hitsAndFracs.size()/std::min(1.,trackPt/2.));
	y[jc] = cornerposxyz.Y() + (cornerposxyz.Y()-posxyz.Y())
	  * (1.00+0.50/hitsAndFracs.size()/std::min(1.,trackPt/2.));
	
#ifdef PFLOW_DEBUG
	if( debug ){
	  std::cout<<"corners "<<jc
	      << " " << cornerposxyz.X()
	      << " " << cornerposxyz.Y()
	      << std::endl;
	}
#endif
      }//loop corners
      
      //need to close the polygon in order to
      //use the TMath::IsInside fonction from root lib
      x[4] = x[0];
      y[4] = y[0];
      
      //Check if the extrapolation point of the track falls 
      //within the rechit boundaries
      bool isinside = TMath::IsInside(track_X,
				      track_Y,
				      5,x,y);
      
      if( isinside ){
	linkedbyDetId = true;
	break;
      }
    }//
    
  }//loop rechits
  
  if( linkedbyDetId ) {
#ifdef PFLOW_DEBUG
    if( debug ) 
      std::cout << "Track and Cluster LINKED BY DetId" << std::endl;
#endif
    /*    
    //if ( distance > 40. || distance < -100. ) 
    double clusterr = std::sqrt(clusterX*clusterX+clusterY*clusterY);
    double trackr = std::sqrt(track_X*track_X+track_Y*track_Y);
    if ( distance > 40. ) 
    std::cout << "Distance = " << distance 
    << ", Barrel/Hcal/Brem ? " << barrel << " " << hcal << " " << isBrem << std::endl
    << " Cluster " << clusterr << " " << clusterZ << " " << clusterphi << " " << clustereta << std::endl
    << " Track   " << trackr << " " << track_Z << " " << trackphi << " " << tracketa << std::endl;
    if ( !barrel && fabs(trackEta) < 1.0 ) { 
      double clusterr = std::sqrt(clusterX*clusterX+clusterY*clusterY);
      double trackr = std::sqrt(track_X*track_X+track_Y*track_Y);
      std::cout << "TrackEta/Pt = " << trackEta << " " << trackPt << ", distance = " << distance << std::endl 
		<< ", Barrel/Hcal/Brem ? " << barrel << " " << hcal << " " << isBrem << std::endl
		<< " Cluster " << clusterr << " " << clusterZ << " " << clusterphi << " " << clustereta << std::endl
		<< " Track   " << trackr << " " << track_Z << " " << trackphi << " " << tracketa << " " << trackEta << " " << trackPt << std::endl;
    } 
    */
    return dist;
  } else {
    return -1.;
  }

} // end of LinkByDetId::testTrackAndClusterByDetId



double
LinkByDetId::testECALAndPSByDetId( const reco::PFCluster& clusterECAL, 
				     const reco::PFCluster& clusterPS,
					  const CaloGeometry* caloGeometry,
				     bool debug)  {

  static double resPSpitch = 0.19/sqrt(12.); // 0.19 -> strip_pitch
  static double resPSlength = 6.1/sqrt(12.); // 6.1  -> strip_length

  // Check that clusterECAL is in ECAL endcap and that clusterPS is a preshower cluster
  if ( clusterECAL.layer() != PFLayer::ECAL_ENDCAP ||
       ( clusterPS.layer() != PFLayer::PS1 && 
	 clusterPS.layer() != PFLayer::PS2 ) ) return -1.;

#ifdef PFLOW_DEBUG
  if( debug ) 
    std::cout<<"entering test link by DetId function for ECAL and PS"<<std::endl;
#endif

  //ECAL cluster position
  double zECAL  = clusterECAL.position().Z();
  double xECAL  = clusterECAL.position().X();
  double yECAL  = clusterECAL.position().Y();

  // PS cluster position, extrapolated to ECAL
  double zPS = clusterPS.position().Z();
  double xPS = clusterPS.position().X(); //* zECAL/zPS;
  double yPS = clusterPS.position().Y(); //* zECAL/zPS;
// MDN jan09 : check that zEcal and zPs have the same sign
	if (zECAL*zPS <0.) return -1.;
  double deltaX = 0.;
  double deltaY = 0.;
  double sqr12 = std::sqrt(12.);
  switch (clusterPS.layer()) {
  case PFLayer::PS1:
    // vertical strips, measure x with pitch precision
    deltaX = resPSpitch * sqr12;
    deltaY = resPSlength * sqr12;
    break;
  case PFLayer::PS2:
    // horizontal strips, measure y with pitch precision
    deltaY = resPSpitch * sqr12;
    deltaX = resPSlength * sqr12;
    break;
  default:
    break;
  }

  // Get the DetIds
  const std::vector< std::pair<DetId, float> >& hitsAndFracs = clusterECAL.hitsAndFractions();
  bool linkedbyDetId = false;
  //loop over hitsAndFractions stored in the cluster
  std::vector< std::pair<DetId, float> >::const_iterator it_haf = hitsAndFracs.begin();
  for(; it_haf != hitsAndFracs.end(); ++it_haf){

    if (it_haf->second < 1E-4) continue;
	const DetId cellId = it_haf->first;
	//-----------
		// first get the cell center position and corners 
	const CaloCellGeometry* cellGeometry = caloGeometry->getGeometry(cellId);
	
	// get the position of the cell
	const GlobalPoint gpPosition = cellGeometry->getPosition();
	// convert to types used in the PF framework
	//
	double zPsEcalScale = zPS/zECAL;
    const math::XYZPoint posxyz(zPsEcalScale*gpPosition.x(), zPsEcalScale*gpPosition.y(), zPsEcalScale*gpPosition.z());

 //   const reco::PFRecHit::REPPoint posrep(gpPosition.perp(), gpPosition.eta(), gpPosition.phi());
 //   posrep = posrep * (zPS/zECAL);
	
    //get the cell corners
	CaloCellGeometry::CornersVec cellGeoCorners = cellGeometry->getCorners();
	// convert to types used in the PF framework
	// also note that corner indexing of the cell face in PFRecHit is the reverse of the one returned 
	// by CaloCellGeometry!
	// these are just the front corners
   std::vector< math::XYZPoint > corners;
	for (uint i=0; i<4; ++i) {
		corners.push_back(math::XYZPoint(cellGeoCorners[3-i].x(),
											cellGeoCorners[3-i].y(),
											cellGeoCorners[3-i].z()
											));
	} // filling the corners
	
	
	
	//------
	
	

#ifdef PFLOW_DEBUG
    if( debug ){
      std::cout << "Ecal rechit " << posxyz.X() << " "   << posxyz.Y() << std::endl;
      std::cout << "PS cluster  " << xPS << " " << yPS << std::endl;
    }
#endif
    
    double x[5];
    double y[5];
    for ( unsigned jc=0; jc<4; ++jc ) {
      // corner position projected onto the preshower
      math::XYZPoint cornerpos = corners[jc] * zPS/zECAL;
      // Inflate the size by the size of the PS strips, and by 5% to include ECAL cracks.
      x[jc] = cornerpos.X() + (cornerpos.X()-posxyz.X()) * (0.05 +1.0/fabs((cornerpos.X()-posxyz.X()))*deltaX/2.);
      y[jc] = cornerpos.Y() + (cornerpos.Y()-posxyz.Y()) * (0.05 +1.0/fabs((cornerpos.Y()-posxyz.Y()))*deltaY/2.);
      
#ifdef PFLOW_DEBUG
      if( debug ){
	std::cout<<"corners "<<jc
	    << " " << cornerpos.X() << " " << x[jc] 
	    << " " << cornerpos.Y() << " " << y[jc]
	    << std::endl;
      }
#endif
    }//loop corners
    
    //need to close the polygon in order to
    //use the TMath::IsInside fonction from root lib
    x[4] = x[0];
    y[4] = y[0];
    
    //Check if the extrapolation point of the track falls 
    //within the rechit boundaries
    bool isinside = TMath::IsInside(xPS,yPS,5,x,y);
      
    if( isinside ){
      linkedbyDetId = true;
      break;
    }

  }//loop rechits
  
  if( linkedbyDetId ) {
    if( debug ) std::cout << "Cluster PS and Cluster ECAL LINKED BY DetId" << std::endl;
    double dist = computeDist( xECAL/1000.,yECAL/1000.,
			       xPS/1000.  ,yPS/1000, 
			       false);    
    return dist;
  } else { 
    return -1.;
  }

}



double 
LinkByDetId::testHFEMAndHFHADByDetId(const reco::PFCluster& clusterHFEM, 
				      const reco::PFCluster& clusterHFHAD,
				      bool debug) {
  
  math::XYZPoint posxyzEM = clusterHFEM.position();
  math::XYZPoint posxyzHAD = clusterHFHAD.position();

  double dX = posxyzEM.X()-posxyzHAD.X();
  double dY = posxyzEM.Y()-posxyzHAD.Y();
  double sameZ = posxyzEM.Z()*posxyzHAD.Z();

  if(sameZ<0) return -1.;

  double dist2 = dX*dX + dY*dY; 

  if( dist2 < 0.1 ) {
    // less than one mm
    double dist = sqrt( dist2 );
    return dist;;
  }
  else 
    return -1.;

}

double
LinkByDetId::computeDist( double eta1, double phi1, 
			   double eta2, double phi2,
			   bool etaPhi )  {
  
  double phicor = etaPhi ? normalizedPhi(phi1 - phi2) : phi1 - phi2;
  
  // double chi2 =  
  //  (eta1 - eta2)*(eta1 - eta2) / ( reta1*reta1+ reta2*reta2 ) +
  //  phicor*phicor / ( rphi1*rphi1+ rphi2*rphi2 );

  double dist = std::sqrt( (eta1 - eta2)*(eta1 - eta2) 
			  + phicor*phicor);

  return dist;

}
