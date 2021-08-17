#include "TOFEstimators.h"
#include "TOFUtils.h"

#include <cmath>
#include <algorithm>

#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/Operators.h>
#include <UTIL/PIDHandler.h>

#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
using EVENT::LCCollection, EVENT::ReconstructedParticle;

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"

#include "TGraphErrors.h"
#include "TF1.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "HelixClass.h"
#include "marlinutil/CalorimeterHitType.h"
#include <UTIL/PIDHandler.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "marlin/ProcessorEventSeeder.h"
#include <marlin/AIDAProcessor.h>
#include <marlin/Global.h>


#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/ICloud2D.h>
#include <AIDA/IHistogram2D.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace lcio ;
using namespace marlin ;
using namespace TOFUtils ;

using ROOT::Math::XYZVector;

TOFEstimators aTOFEstimators ;


TOFEstimators::TOFEstimators() : Processor("TOFEstimators") {
    _description = "TOFEstimators compute some estimators for the time of flight from calorimeter hits" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "ReconstructedParticleCollection",
                            "Name of the ReconstructedParticle collection",
                            _colNamePFO,
                            std::string("PandoraPFOs") );

    registerProcessorParameter( "MaxLayerNumber",
                            "Use only calorimeter hits up to MaxLayerNumber in TOF estimators",
                            _maxLayerNum,
                            int(100) );

    registerProcessorParameter( "CylRadius",
                            "Cut-off hits further away from the shower core than this raduius for TOF calculation",
                            _cylRadiusCut,
                            double(5.) );


    registerProcessorParameter( "TimeResolution",
                            "Assumed time resolution per hit in ps",
                            _resolution,
                            float(0.) );

    registerProcessorParameter( "ProcessorVersion",
                            "Legacy or new TOFEstimator",
                            _procVersion,
                            std::string("idr") );
}


void TOFEstimators::init() {
    Global::EVENTSEEDER->registerProcessor(this);
    if ( _procVersion == "idr" ){
        // initialize gsl random generator
        _rng = gsl_rng_alloc(gsl_rng_ranlxs2);
        _TOFNames = {"TOFFirstHit", "TOFClosestHits", "TOFClosestHitsError", "TOFFlightLength",  "TOFLastTrkHit" , "TOFLastTrkHitFlightLength"};
    }
    else if (_procVersion == "dev"){
        _smearing.param( std::normal_distribution<double>::param_type(0., _resolution/1000.) );
        _TOFNames = {"TOFClosest", "TOFFastest", "TOFCylFit", "TOFClosestFit", "FlightLength", "MomAtCalo"};

        const auto& detector = dd4hep::Detector::getInstance();
        detector.field().magneticField({0., 0., 0.}, _bField);
    }
    else { 
        throw EVENT::Exception(std::string("Invalid ProcessorVersion parameter passed: \'") + _procVersion + std::string("\'\n Viable options are idr (default), dev")); 
    }
}


void TOFEstimators::processEvent( LCEvent * evt ) {

    if ( _procVersion == "idr" ){
        // use the global Marlin random seed for this processor
        gsl_rng_set( _rng, Global::EVENTSEEDER->getSeed(this) ) ;

        streamlog_out(DEBUG ) << "seed set to " << Global::EVENTSEEDER->getSeed(this) << std::endl;

        streamlog_out(DEBUG2) << "   process event: " << evt->getEventNumber()
        		<< "   in run:  " << evt->getRunNumber() << std::endl ;


        // get the PFO collection from the event if it exists
        LCCollection* colPFO = nullptr ;

        try{
        colPFO = evt->getCollection( _colNamePFO ) ;
        }
        catch(lcio::Exception&){
        streamlog_out( DEBUG6 ) << " collection " << _colNamePFO
        		    << " not found in event - nothing to do ... " << std::endl ;
        }

        if( colPFO->getTypeName() != LCIO::RECONSTRUCTEDPARTICLE ) {

        streamlog_out( ERROR ) << " collection " << _colNamePFO
        		   << " not of type LCIO::RECONSTRUCTEDPARTICLE " << std::endl ;

        colPFO = nullptr ;
        }

        if( colPFO != nullptr ){



        PIDHandler pidh( colPFO );
        int algoID = pidh.addAlgorithm( name()  , _TOFNames);


        int nPFO = colPFO->getNumberOfElements()  ;


        for(int i=0; i< nPFO ; ++i){

          ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(  colPFO->getElementAt( i ) ) ;

          bool isCharged = false ;

          if(  pfo->getClusters().size() != 1 ){

        streamlog_out( DEBUG1 ) << " ignore particle w/ cluster number other than one:  " <<  *pfo << std::endl ;
        continue ;
          }

          if( fabs( pfo->getCharge() ) < 0.1  && pfo->getTracks().size() == 0  ) {;}
          else if ( fabs( pfo->getCharge() ) > 0.1  && pfo->getTracks().size() == 1  ) {
              isCharged = true ;}
          else {
            streamlog_out( DEBUG1 ) << " ignore particle w/ track number other than zero or one:  " <<  *pfo << std::endl ;
            continue ;
          }

          streamlog_out( DEBUG1 ) << " --- compute TOF estimators for particle : " << *pfo << std::endl ;



          // =======  use only Ecal hits  (requires the CalorimeterHitType to be set in the digitizer )
          //          with time information ( > 1 ps) and layer <= max layer

          int maxLayerNum = _maxLayerNum ;
          std::function<bool(CalorimeterHit*)> selectHits =  [maxLayerNum](CalorimeterHit* h){

        return ( isEcal( h )            &&
        	 h->getTime() > 1.e-3   &&
        	 layer( h ) <=  maxLayerNum  ) ;
          } ;

          // -------------------------------------------------------------------------------------------

          Cluster* clu = pfo->getClusters()[0] ;
          const CalorimeterHitVec& cluhv = clu->getCalorimeterHits() ;

          // create vectors of extended handle objects for relevant calorimeter hits
          // one w/ unique_ptr for memory handling

          CaloHitUPtrVec uniqueVec ;
          uniqueVec.reserve( cluhv.size()  ) ;

          CaloHitDataVec caloHitVec ;
          caloHitVec.reserve( cluhv.size()  ) ;
          //-------------------------------------------------------------------------------

          CaloHitLayerMap layerMap ;

          for( auto* clh : cluhv ){

        if( selectHits( clh ) ){

          uniqueVec.push_back(  std::unique_ptr<CaloHitData>( new CaloHitData( clh) )  ) ;

          CaloHitData* ch = uniqueVec.back().get() ;

          caloHitVec.push_back( ch ) ;

          ch->layer = layer( ch->lcioHit ) ;
          ch->timeResolution = _resolution ;

          ch->smearedTime  = ( _resolution > 0. ?
        		       ch->lcioHit->getTime() + gsl_ran_gaussian( _rng, _resolution / 1000. ) : // convert ps to ns
        		       ch->lcioHit->getTime() ) ;

          ch->distanceFromIP = dd4hep::rec::Vector3D( clh->getPosition() ).r() ;


          layerMap[ ch->layer ].push_back( ch ) ;
        }
          }


          if( layerMap.empty() ) {

        streamlog_out( DEBUG1 ) << " --- not suitable Ecal hits found for particle " << std::endl ;
        continue ;
          }

          // streamlog_out( DEBUG ) << " --- map with hits per layer : " << std::endl ;
          // for( auto m : layerMap ){
          // 	streamlog_out( DEBUG ) << "  ----- layer " << m.first << " : " << std::endl ;
          // 	for( auto ch : m.second )
          // 	  streamlog_out( DEBUG ) << "            " << caloTypeStr( ch->lcioHit ) << std::endl ;
          // }


          // --- define reference point: track state at calo for charged - hit closest to IP for neutral
          //     and direction of straight line - either from IP or from track state at calo

          dd4hep::rec::Vector3D refPoint ;
          dd4hep::rec::Vector3D unitDir ;
          float flightLength = 0. ;

          TrackerHit* lastTrackerHit = nullptr ;
          float flightLengthTrkHit = 0. ;

          if( isCharged ){

        Track* trk =  pfo->getTracks()[0] ;
        const TrackState* tscalo = trk->getTrackState( TrackState::AtCalorimeter ) ;

        refPoint = tscalo->getReferencePoint() ;

        float tanL = tscalo->getTanLambda() ;
        float theta = atan( 1. / tanL ) ;

        unitDir = dd4hep::rec::Vector3D( 1. ,  tscalo->getPhi() , theta , dd4hep::rec::Vector3D::spherical ) ;

        flightLength = computeFlightLength( trk ) ;


        // also store time and flight length of last tracker hit
        const TrackState* tsIP = trk->getTrackState( TrackState::AtIP ) ;
        	const TrackState* tslh = trk->getTrackState( TrackState::AtLastHit ) ;

        lastTrackerHit = trk->getTrackerHits().back() ;

        flightLengthTrkHit = computeFlightLength( tsIP , tslh ) ;


        #if 1
        if( lastTrackerHit->getTime() > 1e-3 ) {
          dd4hep::rec::Vector3D rpLH = tslh->getReferencePoint() ;
          dd4hep::rec::Vector3D lhp = lastTrackerHit->getPosition()  ;
          streamlog_out( DEBUG3 ) << " *************** referenece point calo     : " << refPoint << std::endl ;
          streamlog_out( DEBUG3 ) << " *************** referenece point last hit : " << rpLH << std::endl ;
          streamlog_out( DEBUG3 ) << " *************** poistion         last hit : " << lhp << std::endl ;
          streamlog_out( DEBUG3 ) << "   distance hit-trkstate: " << (rpLH - lhp ).r() << " --  distance  calo/last hit ref points : " << (refPoint-rpLH).r() << std::endl ;
          streamlog_out( DEBUG3 ) << "   flight lengths:  " << flightLength << "  - " << flightLengthTrkHit << "  -- diff " << flightLength - flightLengthTrkHit <<
            " time diff: " << (flightLength - flightLengthTrkHit) / CLHEP::c_light << std::endl ;
          streamlog_out( DEBUG3 ) << " track state : " << *tslh << std::endl ;
          streamlog_out( DEBUG3 ) << " last hit : " << *lastTrackerHit << std::endl ;
        }
        #endif


          } else {  // neutral particle

        CaloHitDataVec chv =  layerMap.begin()->second ; // only look in first layer w/ hits

        CaloHitData* closestHit =
          *min_element( chv.begin() , chv.end () ,
        		    [](CaloHitData* c0, CaloHitData* c1 ){ return c0->distanceFromIP < c1->distanceFromIP  ; }
            ) ;

        refPoint = closestHit->lcioHit->getPosition() ;

        flightLength = refPoint.r() ;


        dd4hep::rec::Vector3D cluPos = clu->getPosition() ;

        unitDir = cluPos.unit() ;

          }

          streamlog_out( DEBUG2 ) << " ----- use reference point for TOF : " << refPoint << std::endl ;

          streamlog_out( DEBUG ) << " -----  calorimeter hits considered for estimators : " << std::endl ;


          // ------  loop again over hits and fill missing data

          for( auto ch : caloHitVec ){

        CalorimeterHit* calohit = ch->lcioHit ;

        dd4hep::rec::Vector3D pos = ch->lcioHit->getPosition() ;

        ch->distanceFromReferencePoint = ( pos - refPoint ).r()   ;

        ch->distancefromStraightline = computeDistanceFromLine( calohit, refPoint, unitDir ) ;

        streamlog_out( DEBUG ) <<  "     ----- " << caloTypeStr( calohit )
        		       <<  " --  " <<  ch->toString()   <<   std::endl ;
          }

          // ---- now get hits that are closest to the extrapolated line
          CaloHitDataVec tofHits = findHitsClosestToLine( layerMap ) ;


          streamlog_out( DEBUG )   <<  " ***** hits used for the TOF estimator : " << std::endl ;
          for( auto ch : tofHits ){
        streamlog_out( DEBUG ) <<  "     ----- " << ch->toString() << std::endl ;
          }

          auto t_dt     = computeTOFEstimator( tofHits ) ;

          float tof_fh = tofHits[0]->smearedTime - tofHits[0]->distanceFromReferencePoint / CLHEP::c_light ;

          streamlog_out( DEBUG2 ) << "  #### tof ( first ) : " <<  tof_fh << " +/- " << 0 << std::endl ;
          streamlog_out( DEBUG2 ) << "  #### tof ( straight line ) : " << t_dt.first << " +/- " << t_dt.second << std::endl ;


          float trkHitTime = ( lastTrackerHit ?   lastTrackerHit->getTime()   : 0.  ) ;


          FloatVec TOF_params = { tof_fh,
        		      t_dt.first, t_dt.second,
        		      flightLength , trkHitTime , flightLengthTrkHit } ;

          pidh.setParticleID( pfo , 0, 0 , 0.0 , algoID, TOF_params );


        //=========================================================================================

        }

        }
    }
    else if (_procVersion == "dev"){
        _generator.seed( Global::EVENTSEEDER->getSeed(this) );

        LCCollection* colPFO = evt->getCollection(_colNamePFO);
        PIDHandler handler( colPFO );
        int algoID = handler.addAlgorithm( name() , _TOFNames);

        for (int i=0; i<colPFO->getNumberOfElements(); ++i){
            ReconstructedParticle* pfo = dynamic_cast <ReconstructedParticle*> ( colPFO->getElementAt(i) );
            int nClusters = pfo->getClusters().size();
            int nTracks = pfo->getTracks().size();

            // Only simple cases of PFOs
            if( nClusters != 1 || nTracks != 1) continue;

            const Track* track = pfo->getTracks()[0];
            const Cluster* cluster = pfo->getClusters()[0];

	    float momAtCalo = getMomAtCalo(track).r();
            float flightLength = getFlightLength(track);
            float tofClosest = getTOFClosest(track, cluster);
            float tofFastest = getTOFFastest(track, cluster);
            float tofCylFit = getTOFCylFit(track, cluster);
            float tofClosestFit = getTOFClosestFit(track, cluster);

            const std::vector <float> results{tofClosest, tofFastest, tofCylFit, tofClosestFit, flightLength, momAtCalo};
            handler.setParticleID(pfo , 0, 0, 0., algoID, results);
        }
    }
}


void TOFEstimators::check( LCEvent *evt) {
    if( _procVersion == "idr"){
      streamlog_out( DEBUG ) << " --- check called ! " << std::endl ;


      // create some histograms with beta vs momentum for charged particles

      if( isFirstEvent() ){

        // this creates a directory for this processor ....
        AIDAProcessor::tree( this ) ;

        _h.resize(5) ;
        int nBins = 100 ;
        _h[0] = new TH2F( "hbetaFirstHitsChrg", "beta vs momentum - first hit - charged",    nBins, .1 , 10., nBins, 0.93 , 1.03 ) ;
        _h[1] = new TH2F( "hbetaCloseHitsChrg", "beta vs momentum - closest hits - charged", nBins, .1 , 10., nBins, 0.93 , 1.03 ) ;

        _h[2] = new TH2F( "hbetaFirstHitsNeut", "beta vs momentum - first hit - neutral",    nBins, .1 , 10., nBins, 0.93 , 1.03 ) ;
        _h[3] = new TH2F( "hbetaCloseHitsNeut", "beta vs momentum - closest hits - neutral", nBins, .1 , 10., nBins, 0.93 , 1.03 ) ;

        _h[4] = new TH2F( "hbetaLastTrkHit",    "beta vs momentum - last tracker hit", nBins, .1 , 10., nBins, 0.93 , 1.03 ) ;


      }

      // get the PFO collection from the event if it exists
      LCCollection* colPFO = nullptr ;
      try{ colPFO = evt->getCollection( _colNamePFO ) ; } catch(lcio::Exception&){}

      if( colPFO != nullptr ){

        PIDHandler pidh( colPFO );
        int algoID       = pidh.getAlgorithmID( name() );
        int tof_firsthit = pidh.getParameterIndex(algoID,"TOFFirstHit") ;
        int tof_closest  = pidh.getParameterIndex(algoID,"TOFClosestHits") ;
        int tof_length   = pidh.getParameterIndex(algoID,"TOFFlightLength") ;
        int tof_trkhit   = pidh.getParameterIndex(algoID,"TOFLastTrkHit") ;
        int tof_trk_len  = pidh.getParameterIndex(algoID,"TOFLastTrkHitFlightLength") ;

        int nPFO = colPFO->getNumberOfElements()  ;

        for( int i=0 ; i< nPFO ; ++i){

          ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(  colPFO->getElementAt( i ) ) ;

          const ParticleID& tofPID = pidh.getParticleID( pfo , algoID ) ;

          const FloatVec& tofParams = tofPID.getParameters() ;

          streamlog_out( DEBUG ) << " ****  found TOF parameters for pfo w/ size " << tofParams.size() << std::endl ;

          if( !tofParams.empty() ){

    	const double* mom = pfo->getMomentum() ;
    	double momentum = sqrt( mom[0] * mom[0] +  mom[1] * mom[1] +  mom[2] * mom[2] ) ;

    	double length  =  tofParams[ tof_length  ] ;

    	double beta_fh = ( length / tofParams[ tof_firsthit] ) / CLHEP::c_light ;
    	double beta_ch = ( length / tofParams[ tof_closest ] ) / CLHEP::c_light ;


    	if( abs( pfo->getCharge() )  > 0.5 ) {
    	  _h[0]->Fill( momentum , beta_fh );
    	  _h[1]->Fill( momentum , beta_ch );
    	} else {
    	  _h[2]->Fill( momentum , beta_fh );
    	  _h[3]->Fill( momentum , beta_ch );
    	}

    	if( tofParams[ tof_trk_len  ] > 1.e-3 ){ // if TOF from last tracker hit has been set

    	  double beta    = ( tofParams[ tof_trk_len  ] / tofParams[ tof_trkhit  ] ) / CLHEP::c_light ;

    	  _h[4]->Fill( momentum , beta );
    	}

          }
        }
      }
    }
}


void TOFEstimators::end(){
    if( _procVersion == "idr") gsl_rng_free( _rng );
}


dd4hep::rec::Vector3D TOFEstimators::getMomAtCalo(const Track* track){
  const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
  double phi = ts->getPhi();
  double d0 = ts->getD0();
  double z0 = ts->getZ0();
  double omega = ts->getOmega();
  double tanL = ts->getTanLambda();

  HelixClass helix;
  helix.Initialize_Canonical(phi, d0, z0, omega, tanL, _bField[2]/dd4hep::tesla);
  return helix.getMomentum();
}

/*
//Simplistic FligthLength
double TOFEstimators::getFlightLength(const Track* track){
    const TrackState* ts = track->getTrackState(TrackState::AtIP);
    double phiIP = ts->getPhi();
    ts = track->getTrackState(TrackState::AtCalorimeter);
    double phiCalo = ts->getPhi();
    double omegaCalo = ts->getOmega();
    double tanLCalo = ts->getTanLambda();

    return abs((phiIP - phiCalo)/omegaCalo)*sqrt(1. + tanLCalo*tanLCalo);
}
*/

/*
//First attempt of hit by hit length
//copied from Bohdan
double TOFEstimators::getFlightLength(const Track* track){
  double trackLength = 0.;

  //*** track length from IP state to 1st tracker hit state***
  const TrackState* tsIp = track->getTrackState(TrackState::AtIP);
  double phiIP = tsIp->getPhi();
  double omegaIP = tsIp->getOmega();
  double tanLIP = tsIp->getTanLambda();
  const TrackState* tsFirst = track->getTrackState(TrackState::AtFirstHit);
  double phiFirst = tsFirst->getPhi();
  trackLength += std::abs( (phiIP - phiFirst)/omegaIP )*std::sqrt(1. + tanLIP*tanLIP);

  //*** in the tracker region sum hit-by-hit assuming straight line path between hits***
  std::vector <TrackerHit*> trackHits = track->getTrackerHits();

  auto sortByR = [](TrackerHit* a, TrackerHit* b) {
    XYZVector posA = XYZVector(a->getPosition()[0], a->getPosition()[1], a->getPosition()[2]);
    XYZVector posB = XYZVector(b->getPosition()[0], b->getPosition()[1], b->getPosition()[2]);
    return posA.r() < posB.r();
  };
  sort(trackHits.begin(), trackHits.end(), sortByR);

  for (size_t j=1; j < trackHits.size(); ++j){
    XYZVector pos2 = XYZVector(trackHits[j]->getPosition()[0], trackHits[j]->getPosition()[1], trackHits[j]->getPosition()[2]);
    XYZVector pos1 = XYZVector(trackHits[j-1]->getPosition()[0], trackHits[j-1]->getPosition()[1], trackHits[j-1]->getPosition()[2]);
    trackLength += (pos2-pos1).r();
  }


  // track length from last tracker hit to the ECAL surface
  const TrackState* tsLast = track->getTrackState(TrackState::AtLastHit);
  double phiLast = tsLast->getPhi();
  const TrackState* tsCalo = track->getTrackState(TrackState::AtCalorimeter);
  double phiCalo = tsCalo->getPhi();
  double omegaCalo = tsCalo->getOmega();
  double tanLCalo = tsCalo->getTanLambda();

  trackLength += std::abs( (phiLast - phiCalo)/omegaCalo )*std::sqrt(1. + tanLCalo*tanLCalo);
  return trackLength;
}*/



// version of Mitaroff arxiv:2107.02031v1
// approach 1, using only the track states
// equation 20 with phi*=phi
double TOFEstimators::getFlightLength(const Track* track){
  double trackLength = 0.;
  double trackLength2 = 0.;

  double B = _bField[2]/dd4hep::tesla;
  
  //*** track length from IP state to 1st tracker hit state***
  const TrackState* tsIp = track->getTrackState(TrackState::AtIP);
  double phiIp = tsIp->getPhi();
  double d0Ip = tsIp->getD0();
  double z0Ip = tsIp->getZ0();
  double omegaIp = tsIp->getOmega();
  double tanLIp = tsIp->getTanLambda();
  HelixClass helixIp;
  helixIp.Initialize_Canonical(phiIp, d0Ip, z0Ip, omegaIp, tanLIp, _bField[2]/dd4hep::tesla);
  double qIp = helixIp.getCharge();
  double rhIp = helixIp.getRadius();
  dd4hep::rec::Vector3D pIp= helixIp.getMomentum();
  double khIp = 1./rhIp;
 
  XYZVector posIp = XYZVector(tsIp->getReferencePoint()[0], tsIp->getReferencePoint()[1], tsIp->getReferencePoint()[2]);
  streamlog_out(DEBUG) << "AtIP q="<<qIp<<" Bz="<<B<<" Ku="<<CLHEP::c_light<<" z0="<<z0Ip<<" p="<<pIp.r()<<" rh="<<rhIp<<" posR:"<<posIp.r()<<std::endl;

  std::vector <TrackState*> trackSt = track->getTrackStates();
  double phi[1000], k[1000], z[1000];

  phi[0]=phiIp;
  k[0]=khIp;
  z[0]=tsIp->getReferencePoint()[2];

  auto sortByR = [&](TrackState* ts_a, TrackState* ts_b) {
    XYZVector posA = XYZVector(ts_a->getReferencePoint()[0], ts_a->getReferencePoint()[1], ts_a->getReferencePoint()[2]);
    XYZVector posB = XYZVector(ts_b->getReferencePoint()[0], ts_b->getReferencePoint()[1], ts_b->getReferencePoint()[2]);
    return posA.r() < posB.r();
  };
  sort(trackSt.begin(), trackSt.end(), sortByR);
  
  for (size_t j=1; j < trackSt.size(); ++j){
    HelixClass helix;
    double phi_helix = trackSt.at(j)->getPhi();
    double d0 = trackSt.at(j)->getD0();
    double z0 = trackSt.at(j)->getZ0();
    double omega = trackSt.at(j)->getOmega();
    double tanL = trackSt.at(j)->getTanLambda();
    helix.Initialize_Canonical(phi_helix, d0, z0, omega, tanL, _bField[2]/dd4hep::tesla);

    double q = helix.getCharge();
    double rh = helix.getRadius();
    dd4hep::rec::Vector3D p = helix.getMomentum();
    double kh = 1./rh;
    XYZVector pos = XYZVector(trackSt.at(j)->getReferencePoint()[0], trackSt.at(j)->getReferencePoint()[1], trackSt.at(j)->getReferencePoint()[2]);

    streamlog_out(DEBUG) << "TrackState at:"<<trackSt.at(j)->getLocation()<<" q:"<<q<<" Bz="<<B<<" Ku="<<CLHEP::c_light<<" z0="<<z0<<" p="<<p.r()<<" rh="<<rh<<" posR="<<pos.r()<<std::endl;
    phi[j]=phi_helix;
    k[j]=kh;
    z[j]=trackSt.at(j)->getReferencePoint()[2];

    trackLength += sqrt ( pow( (phi[j]-phi[j-1])/k[j-1],2) + pow(z[j]-z[j-1],2) ) ;
  }

  const TrackState* tsCalo= track->getTrackState(TrackState::AtCalorimeter);
  double phiCalo = tsCalo->getPhi();
  double omegaCalo = tsCalo->getOmega();
  double tanLCalo = tsCalo->getTanLambda();
  
  streamlog_out(DEBUG)<< "Simplest one:"<<((phiIp - phiCalo)/omegaCalo)*sqrt(1. + tanLCalo*tanLCalo) << " Mitaroff: "<<trackLength<<std::endl;

  return trackLength;
  }


// version of Mitaroff arxiv:2107.02031v1
// approach 2, using individual hits and the closest track states to define the helix radius
// equation 20 with phi*=phi
// WORK IN PROGRESS.... NOT SATISFACOTRY RESULTS So far
double TOFEstimators::getFlightLengthHitByHit(const Track* track){
  double trackLength = 0.;
  double trackLength2 = 0.;
  
  double trackLengthM1=getFlightLength(track);
  std::vector <TrackerHit*> trackHits = track->getTrackerHits();
  
  auto sortByR = [&](TrackerHit* a, TrackerHit* b) {
    XYZVector posA = XYZVector(a->getPosition()[0], a->getPosition()[1], a->getPosition()[2]);
    XYZVector posB = XYZVector(b->getPosition()[0], b->getPosition()[1], b->getPosition()[2]);
    return posA.r() < posB.r();
  };
  sort(trackHits.begin(), trackHits.end(), sortByR);

  double phi[1000], k[1000], z[1000];
  double phi2[1000];
  
  z[0]=trackHits[0]->getPosition()[2];
  phi2[0]=TMath::ATan(trackHits[0]->getPosition()[1]/trackHits[0]->getPosition()[0]);
  
  const TrackState* trackSt_0 = track->getTrackState(TrackState::AtIP);
  HelixClass helix_0;
  double phi_0 = trackSt_0->getPhi();
  double d0_0 = trackSt_0->getD0();
  double z0_0 = trackSt_0->getZ0();
  double omega_0 = trackSt_0->getOmega();
  double tanL_0 = trackSt_0->getTanLambda();
  helix_0.Initialize_Canonical(phi_0, d0_0, z0_0, omega_0, tanL_0, _bField[2]/dd4hep::tesla);
  double rh_0 = helix_0.getRadius();
  k[0] = 1./rh_0;
  phi[0]=phi_0;
  
  for (size_t j=1; j < trackHits.size(); ++j){
    z[j]=trackHits[j]->getPosition()[2];
    phi2[j]=TMath::ATan(trackHits[j]->getPosition()[1]/trackHits[j]->getPosition()[0]);
    const TrackState* trackSt = track->getClosestTrackState(trackHits[j]->getPosition()[0], trackHits[j]->getPosition()[1], trackHits[j]->getPosition()[2]);
    HelixClass helix;
    double phi_helix = trackSt->getPhi();
    double d0 = trackSt->getD0();
    double z0 = trackSt->getZ0();
    double omega = trackSt->getOmega();
    double tanL = trackSt->getTanLambda();
    helix.Initialize_Canonical(phi_helix, d0, z0, omega, tanL, _bField[2]/dd4hep::tesla);
    double rh = helix.getRadius();
    k[j] = 1./rh;
    phi[j]=phi_helix;
    trackLength += sqrt( pow( (phi[j]-phi[j-1])/k[j-1],2) + pow(z[j]-z[j-1],2) );
    trackLength2 +=  sqrt( pow( (phi2[j]-phi2[j-1])/k[j-1],2) + pow(z[j]-z[j-1],2) );
    streamlog_out(DEBUG)<<j<<" x:"<<trackHits[j]->getPosition()[0]<<" y:"<<trackHits[j]->getPosition()[1]<<" x:"<<trackHits[j]->getPosition()[2]<<" "<<trackSt->getReferencePoint()[2]<<" "<<trackSt->getLocation()<<std::endl;
  }
  //
  const TrackState* trackSt_1 = track->getTrackState(TrackState::AtCalorimeter);
  HelixClass helix_1;
  double phi_1 = trackSt_1->getPhi();
  double d0_1 = trackSt_1->getD0();
  double z0_1 = trackSt_1->getZ0();
  double omega_1 = trackSt_1->getOmega();
  double tanL_1 = trackSt_1->getTanLambda();
  helix_1.Initialize_Canonical(phi_1, d0_1, z0_1, omega_1, tanL_1, _bField[2]/dd4hep::tesla);
  double rh_1 = helix_1.getRadius();
  k[trackHits.size()] = 1./rh_1;
  phi[trackHits.size()]=phi_1;
  z[trackHits.size()]=trackSt_1->getReferencePoint()[2];
  trackLength += sqrt( pow( (phi[trackHits.size()]-phi[trackHits.size()-1])/k[trackHits.size()-1],2) + pow(z[trackHits.size()]-z[trackHits.size()-1],2) );


  //
  const TrackState* ts = track->getTrackState(TrackState::AtIP);
  double phiIP = ts->getPhi();
  const TrackState* tsCalo= track->getTrackState(TrackState::AtCalorimeter);
  double phiCalo = tsCalo->getPhi();
  double omegaCalo = tsCalo->getOmega();
  double tanLCalo = tsCalo->getTanLambda();
  
  streamlog_out(DEBUG)<< "Simplest one:"<<((phiIP - phiCalo)/omegaCalo)*sqrt(1. + tanLCalo*tanLCalo) <<  "  Mitaroff_v1.1:"<<trackLengthM1<<"  Mitaroff_v2.1:"<<trackLength<<"  Mitaroff_v2.2:"<<trackLength2<<std::endl;

  //  trackLength=sqrt(trackLength);
  return trackLength;
}


double TOFEstimators::getTOFClosest(const Track* track, const Cluster* cluster){
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    dd4hep::rec::Vector3D trackAtCaloPos = ts->getReferencePoint();

    double dToImpactMin = std::numeric_limits<double>::max();
    double hitTime = 0.;
    for (const auto& hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isEcal = (hitType.caloID() == CHT::ecal);
        if (!isEcal) continue;

        dd4hep::rec::Vector3D hitPos = hit->getPosition() ;
        double dToImpact = (hitPos - trackAtCaloPos).r();

        if(dToImpact < dToImpactMin){
            dToImpactMin = dToImpact;
            hitTime = hit->getTime();
        }
    }
    hitTime += _smearing(_generator);
    return hitTime - dToImpactMin / CLHEP::c_light;
}


double TOFEstimators::getTOFFastest(const Track* track, const Cluster* cluster){
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    dd4hep::rec::Vector3D trackAtCaloPos = ts->getReferencePoint();

    double dToImpactMin = 0.;
    double hitTimeMin = std::numeric_limits<double>::max();
    for (const auto& hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isEcal = (hitType.caloID() == CHT::ecal);
        if (!isEcal) continue;

        dd4hep::rec::Vector3D hitPos = hit->getPosition();
        double dToImpact = (hitPos - trackAtCaloPos).r();
        double hitTime = hit->getTime() + _smearing(_generator);
        if(hitTime < hitTimeMin){
            dToImpactMin = dToImpact;
            hitTimeMin = hitTime;
        }
    }
    return hitTimeMin - dToImpactMin / CLHEP::c_light;
}


double TOFEstimators::getTOFCylFit(const Track* track, const Cluster* cluster){
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    dd4hep::rec::Vector3D trackAtCaloPos = ts->getReferencePoint();
    dd4hep::rec::Vector3D trackAtCaloMom = getMomAtCalo(track);

    std::vector <double> x, xErr, y, yErr;

    for (const auto& hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isEcal = (hitType.caloID() == CHT::ecal);
        int layer = hitType.layer();
        if (!isEcal || layer >= _maxLayerNum) continue;

        dd4hep::rec::Vector3D hitPos = hit->getPosition();
        double dToLine = (hitPos - trackAtCaloPos).cross(trackAtCaloMom.unit()).r();
        // take only hits from 5 mm cylinders. 5 mm value was ibtain optimizing on photons
        if (dToLine > _cylRadiusCut) continue;

        x.push_back( (hitPos - trackAtCaloPos).r() );
        y.push_back( hit->getTime() + _smearing(_generator) );
        xErr.push_back(0.);
        yErr.push_back( 0.1 );
    }

    if (x.size() <= 1) return 0.;

    TGraphErrors gr(x.size(), &x[0], &y[0], &xErr[0], &yErr[0]);
    gr.Fit("pol1", "Q");

    return gr.GetFunction("pol1")->GetParameter(0);
}


double TOFEstimators::getTOFClosestFit(const Track* track, const Cluster* cluster){
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    dd4hep::rec::Vector3D trackAtCaloPos = ts->getReferencePoint();
    dd4hep::rec::Vector3D trackAtCaloMom = getMomAtCalo(track);

    std::map<int, double> dToImpact;
    std::map<int, double> hitTime;
    std::map<int, double> dToLineMin;
    for(int i=0; i<_maxLayerNum; ++i) dToImpact[i] = 0.;
    for(int i=0; i<_maxLayerNum; ++i) hitTime[i] = 0.;
    for(int i=0; i<_maxLayerNum; ++i) dToLineMin[i] = std::numeric_limits<double>::max();

    for (const auto& hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isEcal = (hitType.caloID() == CHT::ecal);
        int layer = hitType.layer();
        if (!isEcal || layer >= _maxLayerNum) continue;

        dd4hep::rec::Vector3D hitPos = hit->getPosition();
        double dToLine = (hitPos - trackAtCaloPos).cross(trackAtCaloMom.unit()).r();
        if( dToLine < dToLineMin[layer] ){
            dToLineMin[layer] = dToLine;
            dToImpact[layer] = (hitPos - trackAtCaloPos).r();
            hitTime[layer] = hit->getTime();
        }
    }

    std::vector <double> x, xErr, y, yErr;

    for(int i=0; i<_maxLayerNum; ++i){
        if (hitTime[i] <= 0.) continue;
        x.push_back(dToImpact[i]);
        y.push_back( hitTime[i] + _smearing(_generator) );
        xErr.push_back(0.);
        yErr.push_back(0.1);
    }

    if (x.size() <= 1) return 0.;

    TGraphErrors gr(x.size(), &x[0], &y[0], &xErr[0], &yErr[0]);
    gr.Fit("pol1", "Q");

    return gr.GetFunction("pol1")->GetParameter(0);
}
