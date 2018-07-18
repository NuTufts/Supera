#ifndef __SUPERAINFILL_CXX__
#define __SUPERAINFILL_CXX__

#include "SuperaInfill.h"
#include "Instance2Image.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"

namespace larcv {

  static SuperaInfillProcessFactory __global_SuperaInfillProcessFactory__;

  SuperaInfill::SuperaInfill(const std::string name)
    //: ProcessBase(name)
    : SuperaBase(name)
  {}
    
  void SuperaInfill::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ImageMetaMaker::configure(cfg);
    //supera::ParamsImage2D::configure(cfg);

    LARCV_DEBUG() << "start\n";
    _wire_v.clear();
    _wire_v = cfg.get<std::vector<size_t> >("WireList",_wire_v);
    _image_producer = cfg.get<std::string>("ImageProducer");
    _output_image_label = cfg.get<std::string>("OutputImageLabel");
    _output_image_label2 = cfg.get<std::string>("OutputImageLabel2");
    _output_image_label3 = cfg.get<std::string>("OutputImageLabel3");
    _chstatus_producer = cfg.get<std::string>("ChStatusProducer","");
    _plane_id = cfg.get<int>("PlaneID",-1);
    _mask_val = cfg.get<float>("MaskValue",0);
    _threshold = (chstatus::ChannelStatus_t)(cfg.get<unsigned short>("ChStatusThreshold",chstatus::kGOOD));
    m_ancestor_label = cfg.get<std::string>("AncestorImageLabel");
    m_instance_label = cfg.get<std::string>("InstanceImageLabel");

    if(_wire_v.empty() && _chstatus_producer.empty()) {
      LARCV_CRITICAL() << "Neither wire list nor ChStatus producer name given. Nothing to mask!" << std::endl;
      throw larbys();
    }
  }

  void SuperaInfill::initialize()
  {}

  bool SuperaInfill::process(IOManager& mgr)
  {

    SuperaBase::process(mgr);


    if(supera::PulledPork3DSlicer::Is(supera::ImageMetaMaker::MetaMakerPtr())) {
      auto ptr = (supera::PulledPork3DSlicer*)(supera::ImageMetaMaker::MetaMakerPtr());
      ptr->ClearEventData();
      ptr->AddConstraint(LArData<supera::LArMCTruth_t>());
      ptr->GenerateMeta(LArData<supera::LArSimCh_t>(),TimeOffset());
    }
    //get Meta
    auto const& filler_meta_v = Meta();
    if(filler_meta_v.empty()) {
       LARCV_CRITICAL() << "Filler meta not created!" << std::endl;
       throw larbys();
       }    

    //Calling on outputs from chstatus and wire
    //SuperaWire output is already masked. Need a copy filled back in
    const EventImage2D* ev_wire = (EventImage2D*)(mgr.get_data("image2d",_image_producer));
    LARCV_DEBUG() << "ev_wire : " << ev_wire << "\n";
    if(!ev_wire) {
      LARCV_CRITICAL() << "Invalid image producer name: " << _image_producer << std::endl;
      throw larbys();
    }

    EventChStatus* ev_chstatus = nullptr;
    if(!_chstatus_producer.empty()) {
      ev_chstatus = (EventChStatus*)(mgr.get_data("chstatus",_chstatus_producer));
      if(!ev_chstatus) {
	LARCV_CRITICAL() << "ChStatus by " << _chstatus_producer << " not found!" << std::endl;
	throw larbys();
      }
    }
    //Import Instance
    const EventImage2D* ev_instance = (EventImage2D*)(mgr.get_data("image2d",m_instance_label));
    if(!ev_instance) {
      LARCV_CRITICAL() << "ev_instance image could not be created!" << std::endl;
      throw larbys();
    }
    const std::vector<larcv::Image2D>& instance_v = ev_instance->image2d_array();

    //Make a copy and save it so I have an image with ADC values
    auto input_ADC_image = (EventImage2D*)(mgr.get_data("image2d",_output_image_label));
    std::vector<Image2D> ADC_image_v = ev_wire->image2d_array();
  
 
    // Make a second copy to store labels
    auto input_label_image = (EventImage2D*)(mgr.get_data("image2d",_output_image_label2));
    std::vector<Image2D> wire_image_v = ev_wire->image2d_array();
    std::vector<larcv::Image2D> image_label_v = {wire_image_v[0].meta() , wire_image_v[1].meta() , wire_image_v[2].meta()};

    //Make a third copy for weighted labels
    /*auto input_weighted_image = (EventImage2D*)(mgr.get_data("image2d",_output_image_label3));
    std::vector<larcv::Image2D> image_weighted_v = {wire_image_v[0].meta() , wire_image_v[1].meta() , wire_image_v[2].meta()};*/

    //In order to label, (and weighted label), first set all pixels in image to zero
    for(int plane =0; plane<3; plane++)
      {
        image_label_v[plane].paint(0);
       // image_weighted_v[plane].paint(0);
      }
    
    // make sure plane id input is valid
    if(_plane_id >=0 && (int)(image_label_v.size()) <= _plane_id) {
      LARCV_CRITICAL() << "Could not find plane: " << _plane_id << ". image_v.size(): " << image_label_v.size() << std::endl;
      throw larbys();
    }  
    //intialize a pixel count for weighting images
    /*int NCharge = 0;
    int NEmpty = 0;
    int TotalPixels = 0;*/
    
    //loop that: sets holes to zero for labels and fills with instance image for ADC
    //loop through planes
    for(int plane=0; plane<(int)(image_label_v.size()); plane++) {

      // construct wire numbers to mask
      std::set<size_t> wire_s;
      for(auto const& w : _wire_v) wire_s.insert(w);
      if(ev_chstatus) {
	auto const& status_v = ev_chstatus->status(plane).as_vector();
	for(size_t w=0; w<status_v.size(); ++w) {
	  if(status_v[w] < _threshold) wire_s.insert(w);
         }
       }
     //create an image for labels
     auto& img_label = image_label_v.at(plane);
     auto const& meta_label = img_label.meta();
     std::vector<float> empty_column(img_label.meta().rows(),1);
      
     // Loop over wires, find target column and erase
     for(auto const& ch : wire_s) {
	if(ch < meta_label.min_x() || ch >= meta_label.max_x()) continue;
        auto col = meta_label.col((double)ch);
        //set columns in label image to one
        img_label.copy(1,col,empty_column); 
        //find number of rows and loop over rows
        size_t rows = ADC_image_v[plane].meta().rows();
        for (unsigned int r=0; r<rows; r++){
           //TotalPixels++;
           //if instance >=0 (background values are negative), set to 1. else set to zero
           int pix_id = static_cast<int>(instance_v[plane].pixel(r,col));
           if (pix_id >= 0){
               ADC_image_v[plane].set_pixel(r,col,1);
               //NCharge++;
           }
           //else NEmpty++;
         }
      }
      //std::cout << "Number of Charged pixels in plane " << plane << " is " << NCharge << std::endl;
      //std::cout << "Number of empty pixels in plane " << plane << " is " << NEmpty << std::endl;
      //std::cout << "Number of total pixels in plane " << plane << " is " << TotalPixels << std::endl;
 
    // Now that I have all of these totals, I need another for loop to set the masked images (in the same plane still)
      /*for(auto const& w : _wire_v) wire_s.insert(w);
      if(ev_chstatus) {
        auto const& status_v = ev_chstatus->status(plane).as_vector();
        for(size_t w=0; w<status_v.size(); ++w) {
          if(status_v[w] < _threshold) wire_s.insert(w);
         }
       }
      for(auto const& ch : wire_s) {
        if(ch < image_weighted_v[plane].meta().min_x() || ch >= image_weighted_v[plane].meta().max_x()) continue;
        auto col = image_weighted_v[plane].meta().col((double)ch);
        size_t rows = image_weighted_v[plane].meta().rows();
        for (unsigned int r=0; r<rows; r++){
           int pix_id = static_cast<int>(instance_v[plane].pixel(r,col));
           if (pix_id >= 0){
               image_weighted_v[plane].set_pixel(r,col,1.0/NCharge);
           }
           else image_weighted_v[plane].set_pixel(r,col, 1.0/NEmpty);
         }
      }
     NCharge = 0;
     NEmpty = 0;
     TotalPixels  = 0;*/
    }   
    
    // put back images
    //input_weighted_image->emplace(std::move(image_weighted_v));
    input_label_image->emplace(std::move(image_label_v));
    input_ADC_image->emplace(std::move(ADC_image_v));
    return true;
  }

  void SuperaInfill::finalize()
  {}

}
#endif
