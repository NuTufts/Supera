#ifndef __SUPERA_INSTANCE_LAR2IMAGE_H__
#define __SUPERA_INSTANCE_LAR2IMAGE_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "FMWKInterface.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include <map>

namespace supera {

  //
  // SimChannel => PixelFlowMaps
  // 
  std::vector<larcv::Image2D>
    Instance2ImageMichel( const std::vector<larcv::ImageMeta>& meta_v,
			  const std::set<int>& savedSet,
			  const std::vector<supera::LArSimCh_t>& sch_v,
			  const std::vector<float>& row_compression_factor,
			  const std::vector<float>& col_compression_factor,
			  std::vector<larcv::BBox2D>& michel_boxes,
			  const int time_offset );
  
}
#endif
//#endif
//#endif
