#ifndef __SUPERA_LARFLOW_LAR2IMAGE_CXX__
#define __SUPERA_LARFLOW_LAR2IMAGE_CXX__

#include "LArFlowLAr2Image.h"
#include "larcv/core/Base/larcv_logger.h"
#include "ImageMetaMaker.h"


namespace supera {

  //
  // SimChannel => PixelFlowMaps
  // 
  void  SimCh2LArFlowImages( const std::vector<larcv::ImageMeta>& meta_v,
			     const std::vector<larcv::Image2D>& adc_v,
			     const std::vector<unsigned int>& track2type_v,
			     const std::vector<supera::LArSimCh_t>& sch_v,
			     const larcv::EventChStatus& ev_chstatus,
			     const std::vector<float>& row_compression_factor,
			     const std::vector<float>& col_compression_factor,			 
			     const int time_offset,
                             std::vector<larcv::Image2D>& img_out_v,
                             std::vector<larcv::Image2D>& img_vis_v ) {
    
    LARCV_SINFO() << "Filling Pixel-flow truth image..." << std::endl;

    // flow enum
    enum { kU2V=0, kU2Y, kV2U, kV2Y, kY2U, kY2V };
    const float fThreshold = 7.0;
    const float emptyvalue = -5000.0;
    // we create for each plane:
    //  2 images that list column in other images
    //  2 images that give the probability that pixel is visible

    // to do:
    // we need a charge map to be made as well
    // need above to do compression itself

    std::vector<larcv::Image2D> img_v;     // Flow value
    std::vector<larcv::Image2D> energy_v;  // energy deposition of largest particle
    //std::vector<larcv::Image2D> trackid_v; // track id of particle leaving deposition
    for ( auto const& meta : meta_v ) {
      LARCV_SINFO() << meta.dump();      
      larcv::Image2D img1(meta);
      img1.paint(emptyvalue);
      //img1.paint(-1.0);
      img_v.emplace_back( std::move(img1) );
      larcv::Image2D img2(meta);
      img2.paint(emptyvalue);
      img_v.emplace_back( std::move(img2) );
      
      larcv::Image2D energyimg(meta);
      energyimg.paint(0.0);
      energy_v.emplace_back( std::move(energyimg) );
    }
    LARCV_SDEBUG() << "Pixel Flow vector: " << img_v.size() << std::endl;
    LARCV_SDEBUG() << " Filling Image shape: (row,col)=(" << meta_v.front().rows() << "," << meta_v.front().cols() << ")" << std::endl;
    LARCV_SDEBUG() << " Compression Factors: " << std::endl;
    for (int p=0; p<3; p++) {
      LARCV_SDEBUG() << "  p " << p << ": row=" << row_compression_factor[p] << " col=" << col_compression_factor[p] << std::endl;
    }

    // in order to handle overlapping tracks, we have to set a rule to avoid ambiguity
    // we first fill out the y-plane. this plane is used to resolve duplicates

    int numpixfilled_pass0 = 0;
    int numpixfilled_pass1 = 0;
    for (int ipass=0; ipass<2; ipass++) {
      // first pass: Y-only. This is to provide guidance on how to break ambiguities
      // second pass: U,V planes
    
      for (auto const& sch : sch_v) {
	auto ch = sch.Channel();
	auto const& wid   = ::supera::ChannelToWireID(ch);
	auto const& plane = wid.Plane;

	if ( ipass==0 && plane!=2 )
	  continue; // first pass must be Y-plane
	else if ( ipass==1 && plane==2 )
	  continue;
      
	auto& imgmap1    = img_v.at(2*plane+0);
	auto& imgmap2    = img_v.at(2*plane+1);
	auto& energyimg  = energy_v.at(plane);
	//auto& trackidimg = trackid_v.at(plane);
      
	auto const& meta = imgmap1.meta();

	// is the channel inside the meta window
	float wire = wid.Wire;
	if (wire <= meta.min_x()) continue;
	if (meta.max_x() <= wire) continue;
	if (plane != meta.id()) continue;

	// remove offset to get position inside image
	int col = -1;
	try {
	  col = meta.col(wire);
	}
	catch (...) {
	  LARCV_SCRITICAL() << "Out of bound wire=" << wire << std::endl;
	}

	for (auto const tick_ides : sch.TDCIDEMap()) {
	  int tick = supera::TPCTDC2Tick((double)(tick_ides.first)) + time_offset; // true deposition tick
	  if (tick <= meta.min_y()) continue;
	  if (tick >= meta.max_y()) continue;
	  // Where is this tick in column vector?
	  //size_t index = (size_t)(meta.max_y() - tick);
	  int row   = (int)meta.row(tick); // compressed position
	  if ( row<=0 || row>=(int)meta.rows() ) continue;


	  // now we loop over the energy depositions in this tick
	  double energy = (double)energyimg.pixel(row,col); // use to keep track the most energetic energy deposition at this point

	  std::vector<double> pos3d(3);
	  //::larcv::ROIType_t roi_type =::larcv::kROIUnknown; //comment out
	  //int trackid = -1;
	  std::vector<int> imgcoords(4,-1);	  
	  int nedeps = 0;
	  for (auto const& edep : tick_ides.second) {


	    // for the Y-plane (ipass==0), we store the deposition with the highest energy
	    if (edep.energy < energy ) continue; // lower deposition than before


	    energy = edep.energy;

	    // we skip filling this position if Y-plane (pass=0) already filled this information
	    if ( ipass==1 ) {
	      if ( plane==0 && img_v.at( kU2Y ).pixel( row, col )>emptyvalue ) 
		continue;
	      if ( plane==1 && img_v.at( kV2Y ).pixel( row, col )>emptyvalue )
		continue;
	    }

	    // set img coords
	    // SCE position
	    double x,y,z;
	    x = edep.x;
	    y = edep.y;
	    z = edep.z;
	    supera::ApplySCE(x,y,z);
	    const TVector3 poss3d(x,y,z);
	    for (int p=0; p<3; p++) {
	      imgcoords[p+1] = ::supera::NearestWire( poss3d, p );
	      if ( imgcoords[p+1]<0 )
		imgcoords[p+1] = 0;
	      if ( imgcoords[p+1]>=meta.max_x() )
		imgcoords[p+1] = meta.max_x()-1;
	    }
	    
	    //  PID
	    //roi_type = (::larcv::ROIType_t)temp_roi_type; //comment out
	    nedeps++;
	  }

	  if ( nedeps==0 )
	    continue;
	  
	  if ( /*roi_type!=larcv::kROIUnknown &&*/ energy>0) {
	    // we have non-zero energy deposition, valid edep

	    LARCV_SDEBUG() << "setting labels: (r,c)=(" << row << "," << col << ") "
			   << "projection=(" 
			   << imgcoords[1] << "," 
			   << imgcoords[2] << ","
			   << imgcoords[3] << ")"
			   << std::endl;
	    
	    energyimg.set_pixel( row, col, energy );
	    //trackidimg.set_pixel( row, col, trackid );
	        	    
	    switch( plane ) {
	    case 0:
	      imgmap1.set_pixel( row, col, (float)imgcoords[1+1]-(float)col );
	      imgmap2.set_pixel( row, col, (float)imgcoords[2+1]-(float)col );
	      //imgmap1.set_pixel( row, col, (float)imgcoords[1+1] );
	      //imgmap2.set_pixel( row, col, (float)imgcoords[2+1] );
	      break;
	    case 1:
	      imgmap1.set_pixel( row, col, (float)imgcoords[0+1]-(float)col );
	      imgmap2.set_pixel( row, col, (float)imgcoords[2+1]-(float)col );
	      //imgmap1.set_pixel( row, col, (float)imgcoords[0+1] );
	      //imgmap2.set_pixel( row, col, (float)imgcoords[2+1] );
	      break;
	    case 2:
	      // when filling the Y-plane, we set the other planes as well to keep everything consistent
	      imgmap1.set_pixel( row, col, (float)imgcoords[0+1]-(float)col ); // Y->U
	      imgmap2.set_pixel( row, col, (float)imgcoords[1+1]-(float)col ); // Y->V
	      img_v.at( kU2Y ).set_pixel( row, imgcoords[0+1], (float)col-(float)imgcoords[0+1] );
	      img_v.at( kV2Y ).set_pixel( row, imgcoords[1+1], (float)col-(float)imgcoords[1+1] );
	      img_v.at( kU2V ).set_pixel( row, imgcoords[0+1], (float)imgcoords[1+1]-(float)imgcoords[0+1] );
	      img_v.at( kV2U ).set_pixel( row, imgcoords[1+1], (float)imgcoords[0+1]-(float)imgcoords[1+1] );
	      //imgmap1.set_pixel( row, col, (float)imgcoords[0+1] ); // Y->U
	      //imgmap2.set_pixel( row, col, (float)imgcoords[1+1] ); // Y->V
	      //img_v.at( kU2Y ).set_pixel( row, imgcoords[0+1], (float)col );
	      //img_v.at( kV2Y ).set_pixel( row, imgcoords[1+1], (float)col );
	      //img_v.at( kU2V ).set_pixel( row, imgcoords[0+1], (float)imgcoords[1+1] );
	      //img_v.at( kV2U ).set_pixel( row, imgcoords[1+1], (float)imgcoords[0+1] );
	      energy_v.at(0).set_pixel( row, imgcoords[0+1], energy );
	      energy_v.at(1).set_pixel( row, imgcoords[1+1], energy );
	      break;
	    }
	    
	    if (ipass==0)
	      numpixfilled_pass0++;
	    else
	      numpixfilled_pass1++;
	  }//if parti
	  
	  //column[index] = roi_type;
	}
      }
      //break;
    }//end of pass loop
    
    
    LARCV_SDEBUG() << "LArFlowLAr2Image: num pixels filled pass0=" << numpixfilled_pass0 << " pass1=" << numpixfilled_pass1 << std::endl;

    // max pool compression, we do it ourselves
    //std::cout << "Max pool ourselves: compression factors (row,col)=(" << row_compression_factor.front() << "," << col_compression_factor.front() << ")" << std::endl;
    
    // make output, compressed images
    for ( auto const& img : img_v ) {
      const larcv::ImageMeta& meta = img.meta();
      LARCV_SDEBUG() <<"in meta "<< meta.dump() << std::endl;
      larcv::ImageMeta meta_out(meta.min_x(),meta.min_y(),
      				meta.max_x(), meta.max_y(),
      				int( meta.rows()/row_compression_factor.at(meta.id())),int( meta.cols()/col_compression_factor.at(meta.id()) ), 
      				meta.id(), 
      				(larcv::DistanceUnit_t)(larcv::kUnitWireTime));
      LARCV_SDEBUG() <<"out meta "<< meta_out.dump() << std::endl;
      //larcv::Image2D img_out( meta_out );
      //img_out.paint(emptyvalue); //0
      img_out_v.push_back( img );
      larcv::Image2D vis_out( meta );
      vis_out.paint(-1.0);
      img_vis_v.emplace_back( std::move(vis_out) );
    }

    // for ( auto const& img : img_v ) {
    //   const larcv::ImageMeta& meta = img.meta();
    //   LARCV_SDEBUG() <<"in meta "<< meta.dump() << std::endl;
    //   larcv::ImageMeta meta_out(meta.min_x(),meta.min_y(),
    // 				meta.max_x(), meta.max_y(),
    // 				int( meta.rows()/row_compression_factor.at(meta.id())),int( meta.cols()/col_compression_factor.at(meta.id()) ), 
    // 				meta.id(), 
    // 				(larcv::DistanceUnit_t)(larcv::kUnitWireTime));
    //   LARCV_SDEBUG() <<"out meta "<< meta_out.dump() << std::endl;
    //   larcv::Image2D img_out( meta_out );
    //   img_out.paint(emptyvalue); //0
    //   img_out_v.emplace_back( std::move(img_out) );
    //   larcv::Image2D vis_out( meta_out );
    //   vis_out.paint(-1.0);
    //   img_vis_v.emplace_back( std::move(vis_out) );
    // }
    
    int compressed_pixels_filled = 0;
    // for (size_t iidx=0; iidx<img_out_v.size(); iidx++) {
    //   const larcv::Image2D& img       = img_v[iidx];
    //   larcv::Image2D& imgout          = img_out_v[iidx];
    //   size_t plane = img.meta().id();
    //   const larcv::Image2D& energyimg = energy_v.at( plane );

    //   //std::cout << "outimg (row,col)=" << imgout.meta().rows() << "," << imgout.meta().cols() << ")" << std::endl;

    //   for (int rout=0; rout<(int)imgout.meta().rows(); rout++) {
    // 	for (int clout=0; clout<(int)imgout.meta().cols(); clout++) {
    // 	  // find max pixel to transfer
    // 	  int rmax = 0;
    // 	  int cmax = 0;
    // 	  float enmax = 0;
    // 	  for (int dr=0; dr<(int)row_compression_factor.at(plane); dr++) {
    // 	    for (int dc=0; dc<(int)col_compression_factor.at(plane); dc++) {
	      
    // 	      int r = rout *int(row_compression_factor.at(plane)) + dr;
    // 	      int c = clout*int(col_compression_factor.at(plane)) + dc;

    // 	      if ( img.pixel(r,c)<=emptyvalue ) //0
    // 		continue;

    // 	      float pixenergy = energyimg.pixel( r, c );
    // 	      if ( pixenergy>enmax ) {
    // 		enmax = pixenergy;
    // 		rmax = r;
    // 		cmax = c;
    // 		//std::cout << rmax <<" "<< cmax <<std::endl;
    // 	      }

    // 	    }
    // 	  }

    // 	  // set the output
    // 	  if (enmax>0 ) {
    // 	    imgout.set_pixel( rout, clout, img.pixel( rmax, cmax ) );
    // 	    compressed_pixels_filled++;
    // 	  }
    // 	}
    //   }
    // }//end of loop over index
    
    //fill matchability map
    LARCV_SDEBUG() << "Fill Visibility Map" << std::endl;
    const int st[6] = { 1, 2, 0, 2, 0, 1 };
    for (int idx=0; idx<(int)img_out_v.size(); idx++) {
      
      const larcv::Image2D& flowimg    = img_out_v.at(idx);
      larcv::Image2D& visout           = img_vis_v.at(idx);
      const larcv::ChStatus& status1    = ev_chstatus.status(st[idx]);
      const std::vector<short>& stat1   = status1.as_vector();
      
      for (int rout=0; rout<(int)visout.meta().rows(); rout++){
	for(int cout=0; cout<(int)visout.meta().cols(); cout++){
	  if ( flowimg.pixel(rout,cout)>emptyvalue ) {
	    int outpix = cout + flowimg.pixel(rout,cout);
	    if( stat1[outpix] !=4) 
	      visout.set_pixel( rout, cout, 0);
	    else
	      visout.set_pixel( rout, cout, 1);
	  }
	}
      }      
    }//end of matchability map
    
    
    //assign neighbor label 

    
    LARCV_SDEBUG() << "Assign Neighbor Labels" << std::endl;
    for (size_t iidx=0; iidx<img_out_v.size(); iidx++) {
      larcv::Image2D& imgout          = img_out_v[iidx];
      
      for (int clout=0; clout<(int)imgout.meta().cols(); clout++) {

	// loop over time. check for empty labels, give it value of its neighbor
	for (int rout=0; rout<(int)imgout.meta().rows(); rout++) {
	  int prevrow = rout-1;
	  int nextrow = rout+1;
	  if ( prevrow<0 ) continue; 
	  if ( nextrow<0 || nextrow >= (int)imgout.meta().rows() ) {
	    nextrow = -1;
	    break;
	  }

	  if (imgout.pixel( rout, clout )== emptyvalue && adc_v[imgout.meta().id()].pixel(rout,clout)>fThreshold ) {
	    // if empty label, but has adc value
	    if(imgout.pixel( prevrow, clout ) != emptyvalue) 
	      imgout.set_pixel( rout, clout, imgout.pixel( prevrow, clout ) ); // previous has value
	    else if(nextrow>0 && imgout.pixel( prevrow, clout ) ==emptyvalue && imgout.pixel( nextrow, clout ) !=emptyvalue ) 
	      imgout.set_pixel( rout, clout, imgout.pixel( nextrow, clout ) );
	    else 
	      imgout.set_pixel( rout, clout, 0.0 );
	  }
	}
      }
    }//end of neighbor loop over index

    LARCV_SDEBUG() << "compressed pixels filled: " << compressed_pixels_filled << std::endl;

  }


}
#endif
