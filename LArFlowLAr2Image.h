#ifndef __SUPERA_LARFLOW_LAR2IMAGE_H__
#define __SUPERA_LARFLOW_LAR2IMAGE_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "FMWKInterface.h"
#include "larcv/core/DataFormat/DataFormatTypes.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"

namespace supera {

  //
  // SimChannel => PixelFlowMaps
  // 
  //std::vector<larcv::Image2D>
  void  SimCh2LArFlowImages( const std::vector<larcv::ImageMeta>& meta_v,
			     const std::vector<larcv::Image2D>& adc_v,
			     const std::vector<unsigned int>& track2type_v,			 
			     const std::vector<supera::LArSimCh_t>& sch_v,
			     const larcv::EventChStatus& ev_chstatus,
			     const std::vector<float>& row_compression_factor,
			     const std::vector<float>& col_compression_factor,
			     const int time_offset,
			     std::vector<larcv::Image2D>& img_out_v,
			     std::vector<larcv::Image2D>& img_vis_v);
  
}
#endif
//#endif
//#endif
