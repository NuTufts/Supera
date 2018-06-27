#ifndef __SUPERA_MICHEL_IMAGE_CXX__
#define __SUPERA_MICHEL_IMAGE_CXX__

#include "SuperaMichelImage.h"
#include "Instance2ImageMichel.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/DataFormatTypes.h"
#include "larcv/core/DataFormat/EventBBox.h"
#include "larcv/core/DataFormat/BBox.h"


namespace larcv {

  static SuperaMichelImageProcessFactory __global_SuperaMichelImageProcessFactory__;

  SuperaMichelImage::SuperaMichelImage(const std::string name)
    : SuperaBase(name)
  {}
    
  void SuperaMichelImage::configure(const PSet& cfg)
  {
    // We pass the cfg to our base classes
    // We will access to these parameters through inheritance
    SuperaBase::configure(cfg);             // SuperaBase needs the name of LArSoft data products
    supera::ParamsImage2D::configure(cfg);  // ParamsImage2D stores the output image name for this module
    supera::ImageMetaMaker::configure(cfg); // Stores image parameters

    // Example of grabbing a parameter from cfg
    _origin = cfg.get<int>("Origin");
    _boxlabel = cfg.get<std::string>("BoxLabel");
  }

  void SuperaMichelImage::initialize()
  {}

  bool SuperaMichelImage::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    // PulledPork3DSlicer is a way to crop out an image over all three views in a 3D-consistent manner
    // We can activate it by adding a PulledPork block into SuperaMetaMaker
    // We need to keep this crop consistence across all images made by all modules
    // So we need to use the parameters of the crop from that module
    if(supera::PulledPork3DSlicer::Is(supera::ImageMetaMaker::MetaMakerPtr())) {
      auto ptr = (supera::PulledPork3DSlicer*)(supera::ImageMetaMaker::MetaMakerPtr()); // Get the SuperaMetaMaker instance
      ptr->ClearEventData();
      ptr->AddConstraint(LArData<supera::LArMCTruth_t>());
      ptr->GenerateMeta(LArData<supera::LArSimCh_t>(),TimeOffset());
    }
    
    // Get the ImageMeta (the image parameters)
    auto const& meta_v = Meta();
    
    if(meta_v.empty()) {
      LARCV_CRITICAL() << "Meta not created!" << std::endl;
      throw larbys();
    }
    // Get the container for the output images. We will pass our output image to this unit
    auto ev_image = (EventImage2D*)(mgr.get_data("image2d",OutImageLabel()));
    if(!ev_image) {
      LARCV_CRITICAL() << "Output image could not be created!" << std::endl;
      throw larbys();
    }
    if(!(ev_image->image2d_array().empty())) {
      LARCV_CRITICAL() << "Output image array not empty!" << std::endl;
      throw larbys();
    }

    // Get a container for the output boxes
    auto ev_bbox = (larcv::EventBBox2D*)(mgr.get_data("bbox2d",_boxlabel));

    // Example of how we work with truth data
    // We get the container of MCTrack objects -- what LArSoft uses to store the true
    //   trajectories of particles through the active region of the detector
    //   limited to "track-like" particles: muon, pions, protons, recoiling nuclei (and fragments)
    
    LARCV_DEBUG() << "==============================================================================" << std::endl;
    LARCV_DEBUG() << "SuperaMichelImage:: Store map between trackID and it's ancestor track ID"   << std::endl;

    std::set<int> savedSet;

    std::map<int,int> trackid2ancestorid;

    for(auto const& mctrack : LArData<supera::LArMCTrack_t>()) {
      LARCV_DEBUG() << "mctrack: "
		    << " id=" << mctrack.TrackID() 
		    << " ancestorid=" << mctrack.AncestorTrackID() 
		    << " motherid=" << mctrack.MotherTrackID() 
		    << " pdg=" << mctrack.PdgCode() 
		    << " origin=" << mctrack.Origin()
	            << " process=" << mctrack.Process()
		    << std::endl;

      // we use this to select the type of primary particle we store
      // origin=0: store all
      // origin=1: cosmic rays
      // origin=2: neutrinos
      if(_origin && ((unsigned short)(mctrack.Origin())) != _origin) continue;
      
      trackid2ancestorid[mctrack.TrackID()] = mctrack.AncestorTrackID();
    }
    for(auto const& mcshower : LArData<supera::LArMCShower_t>()) {

      LARCV_DEBUG() << "mcshower: "
		    << " id=" << mcshower.TrackID() 
		    << " ancestorid=" << mcshower.AncestorTrackID() 
		    << " motherid=" << mcshower.MotherTrackID() 
		    << " pdg=" << mcshower.PdgCode() 
		    << " origin=" << mcshower.Origin()
	            << " process=" << mcshower.Process()
		    << std::endl;
      
      if(_origin && ((unsigned short)(mcshower.Origin())) != _origin) continue;
      
      if (mcshower.Process() == "Decay") {
	std::cout << "Found Michel in MCTruth: " 
		  << " id=" << mcshower.TrackID()
		  << " mother=" << mcshower.MotherTrackID() 
		  << " ancestor="<< mcshower.AncestorTrackID() 
		  << std::endl;
	savedSet.insert(mcshower.TrackID());
      }

      trackid2ancestorid[mcshower.TrackID()] = mcshower.AncestorTrackID();
    }
    
    // We get the compression factors
    std::vector<float> row_compression_factor;
    std::vector<float> col_compression_factor;
    for ( auto const& meta : meta_v ) {
      row_compression_factor.push_back( RowCompressionFactor().at(meta.id()) );
      col_compression_factor.push_back( ColCompressionFactor().at(meta.id()) );
    }

    // Get container with data that associates charge seen on the wire with the particle that deposited that charge
    auto simch_data  = LArData<supera::LArSimCh_t>();
    
    // We separate this module with the algorithm that forms the truth image
    // We do this so that we can reuse the algorithm in the larlite framework as well
    // We wrote an algorithm in Instance2ImageMichel.h/.cxx
    std::vector<larcv::BBox2D> michel_boxes;
    auto image_v = supera::Instance2ImageMichel(meta_v, savedSet, simch_data, row_compression_factor, col_compression_factor, michel_boxes, TimeOffset() );
    
    ev_image->emplace( std::move(image_v) );
    for ( auto& bb : michel_boxes )
      ev_bbox->emplace_back(  std::move(bb) );

    // std::cout << "PRINTING SAVEDSET" << std::endl;
    // for (auto& id :savedSet){
    //   std::cout << id << std::endl;
    // }
    
    return true;
  }

  void SuperaMichelImage::finalize()
  {}

}
#endif
