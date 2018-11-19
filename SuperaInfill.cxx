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

    _output_adc = cfg.get<std::string>("OutputADCLabel");
    _output_labelsbasic = cfg.get<std::string>("OutputLabelsBasicLabel");
    _output_adcmasked = cfg.get<std::string>("OutputADCMaskedLabel");
    _output_labels = cfg.get<std::string>("OutputLabelsLabel");
    _adcthreshold = cfg.get<int>("ADCLabelThreshold");   

    _chstatus_producer = cfg.get<std::string>("ChStatusProducer","");
    _plane_id = cfg.get<int>("PlaneID",-1);
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

    //Make a copy and save it so I have an image with ADC values
    auto input_ADC_image = (EventImage2D*)(mgr.get_data("image2d",_output_adc));
    std::vector<Image2D> ADC_image_v = ev_wire->image2d_array();
 
    // Make a second copy to store labels
    // 1 = dead wires
    // 0 = not dead wires
    auto input_labelbasic_image = (EventImage2D*)(mgr.get_data("image2d",_output_labelsbasic));
    std::vector<Image2D> wire_image_v = ev_wire->image2d_array();
    std::vector<larcv::Image2D> image_labelbasic_v = {wire_image_v[0].meta() , wire_image_v[1].meta() , wire_image_v[2].meta()};
 
    //Make a third copy - ADC vales w/ dead channels included
    auto input_ADCMasked_image = (EventImage2D*)(mgr.get_data("image2d",_output_adcmasked));
    std::vector<larcv::Image2D> image_adcmasked_v = ADC_image_v;

    //make a fourth copy - labeled for weights
    //0 = background
    //1 = adc above threshold 10
    //2 = adc >0 but below threshold
    //3 = adc = 0 but within 2 pixels (horizontally) of a pixel with adc>0
    auto input_label_image = (EventImage2D*)(mgr.get_data("image2d",_output_labels));

    //In order to label, first set all pixels in image to zero
    for(int plane =0; plane<3; plane++)
      {
        image_labelbasic_v[plane].paint(0);
      }
    std::vector<larcv::Image2D> image_label_v = image_labelbasic_v;
    
    // make sure plane id input is valid
    if(_plane_id >=0 && (int)(image_labelbasic_v.size()) <= _plane_id) {
      LARCV_CRITICAL() << "Could not find plane: " << _plane_id << ". image_v.size(): " << image_labelbasic_v.size() << std::endl;
      throw larbys();
    }  
   
    //loop that: sets dead channels to one for labels and zero for the masked adc image
    //loop through planes
    for(int plane=0; plane<(int)(image_labelbasic_v.size()); plane++) {

      // construct wire numbers to mask
      std::set<size_t> wire_s;
      for(auto const& w : _wire_v) wire_s.insert(w);
      if(ev_chstatus) {
	auto const& status_v = ev_chstatus->status(plane).as_vector();
	for(size_t w=0; w<status_v.size(); ++w) {
	  if(status_v[w] < _threshold) wire_s.insert(w);
         }
       }

     //create an image for basic labels
     auto& img_labelbasic = image_labelbasic_v.at(plane);
     auto const& meta_labelbasic = img_labelbasic.meta();
     std::vector<float> empty_column(img_labelbasic.meta().rows(),1);
     //create an image for masked ADC
     auto& img_mask = image_adcmasked_v.at(plane);
     std::vector<float> empty_column_masked(img_labelbasic.meta().rows(),0);

     // Loop over wires, set dead channels in labels to one 
     for(auto const& ch : wire_s) {
	if(ch < meta_labelbasic.min_x() || ch >= meta_labelbasic.max_x()) continue;
        auto col = meta_labelbasic.col((double)ch);
        //set columns in label image to one
        img_labelbasic.copy(1,col,empty_column);
	// set columns in adcmasked image to zero
	img_mask.copy(1, col, empty_column_masked); 
      }
      
      //loop over each pixel for label image
      size_t rows = ADC_image_v[plane].meta().rows();
      size_t cols = ADC_image_v[plane].meta().cols();

      for (unsigned int c = 0; c<cols; c++){
        for (unsigned int r=0; r<rows; r++){
          int pix_id = static_cast<int>(ADC_image_v[plane].pixel(r,c));
          // 1: greater than threshold
          if (pix_id >= _adcthreshold){
             image_label_v[plane].set_pixel(r, c, 1);
          }

          // 2: greater than zero, less than threshold
          else if (pix_id > 0){
             image_label_v[plane].set_pixel(r, c , 2);
          }

          // 3: zero, but within 2 pixels in time from a non-zero pixel
          // lots of conditions to avoid seg faults on borders
         else if (( r >= 2) && (r < (rows - 2))){
             if ( (ADC_image_v[plane].pixel(r-1,c) > 0)
                  || (ADC_image_v[plane].pixel(r-2,c) > 0) 
                  || (ADC_image_v[plane].pixel(r+1,c) > 0) 
                  || (ADC_image_v[plane].pixel(r+2,c) > 0)){
                image_label_v[plane].set_pixel(r, c , 3);
             }
          }
          else if (r == 1){
             if ( (ADC_image_v[plane].pixel(r-1,c) > 0)
                  || (ADC_image_v[plane].pixel(r+1,c) > 0)
                  || (ADC_image_v[plane].pixel(r+2,c) > 0)){
                image_label_v[plane].set_pixel(r, c , 3);
             }
          }
          else if (r == 0){
             if ( (ADC_image_v[plane].pixel(r+1,c) > 0)
                  || (ADC_image_v[plane].pixel(r+2,c) > 0)){
                image_label_v[plane].set_pixel(r, c , 3);
             }
          }
          else if (r == rows -2){
             if ( (ADC_image_v[plane].pixel(r-1,c) > 0)
                  || (ADC_image_v[plane].pixel(r-2,c) > 0)
                  || (ADC_image_v[plane].pixel(r+1,c) > 0)){
                image_label_v[plane].set_pixel(r, c , 3);
             }
          }
          else if (r == rows - 1){
             if ( (ADC_image_v[plane].pixel(r-1,c) > 0)
                  || (ADC_image_v[plane].pixel(r-2,c) > 0)){
                image_label_v[plane].set_pixel(r, c , 3);
             }  
          }
        }
      }
    }   
    
    // put back images
    input_labelbasic_image->emplace(std::move(image_labelbasic_v));
    input_label_image->emplace(std::move(image_label_v));
    input_ADC_image->emplace(std::move(ADC_image_v));
    input_ADCMasked_image->emplace(std::move(image_adcmasked_v));
    return true;
  }

  void SuperaInfill::finalize()
  {}

}
#endif
