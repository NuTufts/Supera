#ifndef __SUPERA_INSTANCE_LAR2IMAGE_CXX__
#define __SUPERA_INSTANCE_LAR2IMAGE_CXX__

#include "Instance2ImageMichel.h"
// larlite
//#include "LArUtil/Geometry.h"
// larcv
#include "larcv/core/Base/larcv_logger.h"
#include "larcv/core/DataFormat/EventBBox.h"
#include "larcv/core/DataFormat/BBox.h"

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
			const int time_offset ) {
    
    LARCV_SINFO() << "Instance ID Image ..." << std::endl;
    
    // we pack truth info about the ancestor particle type
    // we label energy deposition by ancestor track id
    // this groups secondaries with their primary parent
    
    // create images we are going to fill
    std::vector<larcv::Image2D> img_v;     // ADC value per pixel
    std::vector<larcv::Image2D> energy_v;  // largest energy deposition
    for ( auto const& meta : meta_v ) {
      LARCV_SINFO() << meta.dump();      
      larcv::Image2D img1(meta);
      img1.paint(-1.0);
      img_v.emplace_back( std::move(img1) );
      larcv::Image2D Eimg(meta);
      Eimg.paint(-1.0);
      energy_v.emplace_back( std::move(Eimg) );
    }

    // loop over sim channel information
    struct michel_box{
      int min_row[3];
      int max_row[3];
      int min_col[3];
      int max_col[3];
    };
    
    std::map<int,michel_box> id_to_michel_box;
    
    int nfilled = 0;
    int num_images = 0;
    for (auto const& sch : sch_v) {
      auto ch = sch.Channel();
      auto const& wid   = ::supera::ChannelToWireID(ch);
      auto const& plane = wid.Plane;

      
      auto& imgmap1    = img_v.at(plane);
      auto& Eimg       = energy_v.at(plane);
      
      auto const& meta = imgmap1.meta();

      // is the channel inside the meta window
      int wire = (int)wid.Wire;
      if (wire < meta.min_x()) continue;
      if (meta.max_x() <= wire) continue;
      if (plane != meta.id()) continue;

      // remove offset to get position inside image
      //col -= (int)(meta.min_x());
      int col = (int)(meta.col(wire));
      
      // loop over energy deposition
      for (auto const tick_ides : sch.TDCIDEMap()) {
	int tick = supera::TPCTDC2Tick((double)(tick_ides.first)) + time_offset; // true deposition tick
	if (tick <= meta.min_y()) continue;
	if (tick >= meta.max_y()) continue;
	// Where is this tick in the image
	int row   = (int)meta.row(tick); // compressed position
	if ( row<0 || row>=(int)meta.rows() ) continue;
	
	// now we loop over the energy depositions in this tick
	double energy = (double)Eimg.pixel(row,col); // use to keep track the most energetic energy deposition at this point
	//int ancestorid   = -1;
	
	// energy deposition at tick
	for (auto const& edep : tick_ides.second) {
	  
	  // for the Y-plane (ipass==0), we store the deposition with the highest energy
	  //if (edep.energy < energy ) continue; // lower deposition than before
	  if ( savedSet.find( edep.trackID )==savedSet.end() ) {
	    num_images++;
	    continue;
	  }
	  
	  auto it_michel_box = id_to_michel_box.find(edep.trackID);
	  if (it_michel_box==id_to_michel_box.end()){
	    std::cout << "making new michel box for trackID " << edep.trackID << std::endl;
	    michel_box mbox;
	    for (int p=0; p<3; p++){
	      mbox.min_row[p]=(unsigned int)1000000;
	      mbox.max_row[p]=(unsigned int)0;
	      mbox.min_col[p]=(unsigned int)1000000;
	      mbox.max_col[p]=(unsigned int)0;
	    }
	    id_to_michel_box.insert(std::pair<int,michel_box>(edep.trackID,mbox));
	  }

	  auto &mbox=(id_to_michel_box.find(edep.trackID))->second;

	  if (row<mbox.min_row[plane]){
	    mbox.min_row[plane]=row;
	  }
	  if (col<mbox.min_col[plane]){
	    mbox.min_col[plane]=col;
	  }
	  if (row>mbox.max_row[plane]){
	    mbox.max_row[plane]=row;
	  }
	  if (col>mbox.max_col[plane]){
	    mbox.max_col[plane]=col;
	  }


	  //std::cout << "michelpixel (r,c)=(" << row  << "," << col << ")" << std::endl;

	  energy = edep.energy;
	  //auto it_map = savedSet.find( edep.trackID );
	  //ancestorid = it_map->second;
	  
	  // we have non-zero energy deposition, valid edep
	  Eimg.set_pixel( row, col, energy );
	  imgmap1.set_pixel( row, col, edep.trackID );
	  nfilled+=1;

	  
	}
      }
      
    }//end of pass loop
    //std::cout << "Number of filled michel pixels: " << nfilled << std::endl;
    
    // for each michel
    for (auto it_mbox : id_to_michel_box) {
      std::vector<larcv::BBox2D> roi_v;
      michel_box& mbox=it_mbox.second;

      std::cout << "mbox_id: " << it_mbox.first << std::endl;

      // for each plane
      for (int i=0;i<3;i++) {
	const larcv::ImageMeta& meta = img_v[i].meta();
	int r=(mbox.max_row[i]-mbox.min_row[i]);
	int c=mbox.max_col[i]-mbox.min_col[i];
	float h=(r*meta.pixel_height());
	float w=(c*meta.pixel_width());
	float x=meta.pos_x(mbox.min_col[i]);
	float y=meta.pos_y(mbox.min_row[i]);

	// std::cout << "c=" << c << " ";
	// std::cout << "w=" << w << " ";
	// std::cout << "col_compression_factor=" << col_compression_factor.at(meta.id()) << std::endl;

	larcv::BBox2D bbox(x,y,x+w,y+h, meta.id());
	//larcv::ImageMeta bbox(w,h,r,c,x,y,meta.id());
	roi_v.emplace_back( std::move(bbox) );
      }

      for ( auto& bb : roi_v )
	michel_boxes.emplace_back(std::move(bb));
    }
    

    // make output, compressed images
    std::vector<larcv::Image2D> img_out_v;
    for ( auto const& img : img_v ) {
      const larcv::ImageMeta& meta = img.meta();
      // larcv1 meta
      // -----------
      // larcv::ImageMeta meta_out(meta.width(), meta.height(), 
      // 				int( meta.rows()/row_compression_factor.at(meta.plane()) ), int( meta.cols()/col_compression_factor.at(meta.plane()) ),
      // 				meta.min_x(), meta.max_y(), 
      // 				meta.plane() );
      // larcv2 meta
      // -----------
      larcv::ImageMeta meta_out(meta.min_x(), meta.min_y(), meta.max_x(), meta.max_y(),
      				int( meta.rows()/row_compression_factor.at(meta.id()) ), int( meta.cols()/col_compression_factor.at(meta.id()) ),
      				meta.id(), meta.unit() );
      
      larcv::Image2D img_out( meta_out );
      img_out.paint(0.0);
      img_out_v.emplace_back( std::move(img_out) );
    }
    
    int compressed_pixels_filled = 0;
    for (size_t iidx=0; iidx<img_out_v.size(); iidx++) {
      const larcv::Image2D& img       = img_v[iidx];
      larcv::Image2D& imgout          = img_out_v[iidx];
      size_t plane = img.meta().id();
      const larcv::Image2D& energyimg = energy_v.at( plane );
      
      for (int rout=0; rout<(int)imgout.meta().rows(); rout++) {
	for (int clout=0; clout<(int)imgout.meta().cols(); clout++) {
	  // find max pixel to transfer
	  int rmax = 0;
	  int cmax = 0;
	  float enmax = 0;
	  for (int dr=0; dr<(int)row_compression_factor.at(plane); dr++) {
	    for (int dc=0; dc<(int)col_compression_factor.at(plane); dc++) {
	      
	      int r = rout *int(row_compression_factor.at(plane)) + dr;
	      int c = clout*int(col_compression_factor.at(plane)) + dc;
	      
	      if ( img.pixel(r,c)<0 )
		continue;
	      
	      float pixenergy = energyimg.pixel( r, c );
	      if ( pixenergy>enmax ) {
		enmax = pixenergy;
		rmax = r;
		cmax = c;
	      }
	      
	    }
	  }
	  
	  // set the output
	  if (enmax>0 ) {
	    imgout.set_pixel( rout, clout, img.pixel( rmax, cmax ) );
	    compressed_pixels_filled++;
	  }
	}
      }
    }//end of loop over index
    
    return img_out_v;
  }
  
  

}
#endif
