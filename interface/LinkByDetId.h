#ifndef RecoParticleFlow_PFClusterTools_LinkByDetId_h
#define RecoParticleFlow_PFClusterTools_LinkByDetId_h 

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"



class LinkByDetId {
 public:
  LinkByDetId() {};
  ~LinkByDetId() {} ; 


  /// computes a chisquare - same as in LinkByRecHit
  static double computeDist( double eta1, double phi1, 
                             double eta2, double phi2,
                             bool etaPhi = true) ;

  
  ///test the association between a track and a cluster using the list of detId's
  static double testTrackAndClusterByDetId( const reco::PFRecTrack& track, 
                                            const reco::PFCluster& cluster,
														  const CaloGeometry* caloGeometry,
                                            bool isBrem = false,
                                            bool debug = false ) ;  

  ///test the association between ECAL and PS clusters using the lists of DetId's
  static double testECALAndPSByDetId( const reco::PFCluster& clusterECAL, 
                                       const reco::PFCluster& clusterPS,
													const CaloGeometry* caloGeometry,
                                       bool debug=false)  ;
  
  /// test the association between HFEM and HFHAD using the lists of DetId's
  static double testHFEMAndHFHADByDetId( const reco::PFCluster& clusterHFEM, 
                                          const reco::PFCluster& clusterHFHAD,
                                          bool debug=false)  ;
  
  // void setTheGeometry(const EcalBarelGeometry& ebGeometry,
					  // const EcalEndcapGeometry& eeGeometry,					  
					  // const CaloSubdetectorGeometry& psGeometry,
					  // const CaloSubdetectorGeometry& hbGeometry,
					  // const CaloSubdetectorGeometry& heGeometry,					  
					  // const CaloSubdetectorGeometry& hoGeometry
					  // ) {
    // ebGeometry_ = ebGeometry;
	// eeGeometry_ = eeGeometry;
	// psGeometry_ = psGeometry;
	// hbGeometry_ = hbGeometry;
	// heGeometry_ = heGeometry;
	// hoGeometry_ = hoGeometry;
  // }
  
 // void setTheGeometry(const CaloGeometry* geometry) {
 //      caloGeometry_ = geometry;
 //}
  
  
  private:
  // check if the geometry is enough or we should follow the approach in PFClusterMaker (AA)
  ////const CaloGeometry* caloGeometry_;
  // const EcalBarelGeometry* ebGeometry_;
  // const EcalEndcapGeometry* eeGeometry_;
  // const CaloSubdetectorGeometry* psGeometry_;
  // const CaloSubdetectorGeometry* hbGeometry_;
  // const CaloSubdetectorGeometry* heGeometry_;
  // const CaloSubdetectorGeometry* hoGeometry_;

  
  
};

#endif