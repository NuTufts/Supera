/**
 * \file SuperaMichelImage.h
 *
 * \ingroup Package_Name
 * 
 * \brief Example of Class def header 
 *
 * @author taritree
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERAMichelImage_H__
#define __SUPERAMichelImage_H__

#include <string>

#include "SuperaBase.h"
#include "FMWKInterface.h"
#include "ImageMetaMaker.h"
#include "ParamsImage2D.h"
#include "larcv/core/DataFormat/Image2D.h"

namespace larcv {

  /**
     \class SuperaMichelImage
     User defined class SuperaMichelImage ... these comments are used to generate
     doxygen documentation!

     This is a supera module. It's roll is to take LArSoft information and produce an image and/or meta deta about the image.
     The module is suppose to be one unit in a chain of units 

     We inheret from three base classes
     SuperaBase: Provides interface functions for supera processor chain. Also provides storage and get/set for common larsoft data products
                 SuperaBase also inherits from a larcv::ProcessBase class. This allows it to be used in the larcv image processing chain
     ParamsImage2D: Provides interface to output image configuration (height,weight,rows,columns,etc.)
     ImageMetaMaker: Provides interface to more image construction functions (e.g. output down-sampling factors for example)

     As a result, this class will inherit functions you should use to get information about the output image that is to be made.
     See cxx for example.
  */
  class SuperaMichelImage : public SuperaBase,
                        public supera::ParamsImage2D,
                        public supera::ImageMetaMaker {

  public:
    
    /// Default constructor
    SuperaMichelImage(const std::string name="SuperaMichelImage");
    
    /// Default destructor
    ~SuperaMichelImage(){}

    // Grab parameters from PSet and store the relevant parameters
    void configure(const PSet&);

    // Called before the event loop starts. Take parameters and setup the class
    void initialize();

    // Called for each event (in the event loop). Take LArSoft data and store image data (and meta-data)
    // mgr provides access to output file (and access to images from previous Supera modules)
    bool process(IOManager& mgr);

    // Called after event loop
    void finalize();

  private:
    
    int            _origin;                 // origin flag: example of pset parameter
    std::string    _boxlabel;               // label for container holding our bounding boxes
    //std::string    _chstatus_producer_name; // chstatus image container name: example of string to be filled with pset parameter
    
  };

  /**
     \class larcv::SuperaMichelImageFactory
     \brief A concrete factory class for larcv::SuperaMichelImage
  */
  class SuperaMichelImageProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaMichelImageProcessFactory() { ProcessFactory::get().add_factory("SuperaMichelImage",this); }
    /// dtor
    ~SuperaMichelImageProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaMichelImage(instance_name); }
  };

}
#endif
/** @} */ // end of doxygen group 

