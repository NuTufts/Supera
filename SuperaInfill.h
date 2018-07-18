/**
 * \file SuperaInfill.h
 *
 * \ingroup Package_Name
 * 
 * \brief Class def header for a class SuperaInfill
 *
 * @author kmason
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERAINFILL_H__
#define __SUPERAINFILL_H__

#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"
#include "SuperaBase.h"
#include "ImageMetaMaker.h"
#include "Instance2Image.h"

namespace larcv {

  /**
     \class SuperaInfill
     This module is used to mask column pixels that correspond to a set of specified wires.\n
  */
  class SuperaInfill : //public ProcessBase,
    public SuperaBase,
   // public supera::ParamsImage2D,
    public supera::ImageMetaMaker {

  public:
    
    /// Default constructor
    SuperaInfill (const std::string name="SuperaInfill");
    
    /// Default destructor
    ~SuperaInfill(){}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:
    std::string m_ancestor_label; //instance image
    std::string m_instance_label; //instance image   
 
    std::string _output_image_label; //test
    std::string _output_image_label2; //test
    std::string _output_image_label3; //test
    
    std::string _image_producer;  ///< Image to mask
    int         _plane_id;        ///< Plane ID (i.e. EventImage2D index number) to mask wires for. <0 means ALL planes
    std::vector<size_t> _wire_v;  ///< A list of wire numbers to be masked
    float       _mask_val;        ///< Value to be used for masking (default 0)
    std::string _chstatus_producer; ///< ChStatus producer name (if using ChStatus to mask)
    chstatus::ChannelStatus_t _threshold;     ///< Threshold status for ChStatus
  };

  /**
     \class larcv::SuperaInfillFactory
     \brief A concrete factory class for larcv::WireMask
  */
  class SuperaInfillProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaInfillProcessFactory() { ProcessFactory::get().add_factory("SuperaInfill",this); }
    /// dtor
    ~SuperaInfillProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaInfill(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group 

