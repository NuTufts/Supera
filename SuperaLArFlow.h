/**
 * \file SuperaLArFlow.h
 *
 * \ingroup Package_Name
 * 
 * \brief Class def header for a class SuperaLArFlow
 *
 * @author taritree
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERALARFLOW_H__
#define __SUPERALARFLOW_H__

#include <string>

#include "SuperaBase.h"
#include "FMWKInterface.h"
#include "ImageMetaMaker.h"
#include "ParamsImage2D.h"
#include "larcv/core/DataFormat/Image2D.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaLArFlow ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaLArFlow : public SuperaBase,
                        public supera::ParamsImage2D,
                        public supera::ImageMetaMaker {

  public:
    
    /// Default constructor
    SuperaLArFlow(const std::string name="SuperaLArFlow");
    
    /// Default destructor
    ~SuperaLArFlow(){}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

    unsigned int PdgCode2ROIType(int pdgcode){
      switch(pdgcode) {
      
	// electron
      case  11:
      case -11:
	return 3;
	// muon
      case  13:
      case -13:
	return 6;
	// gamma
      case 22:
	return 4;
	// neutrinos
      case 12:
      case 14:
	return 2;
	// proton
      case  2212:
      case -2212:
	return 9;
	// pi0
      case  111:
	return 5;
        // pi +/-
      case  211:
      case -211:
	return 8;
	// K +/-
      case  321:
      case -321:
	return 7;
      default:
	return 0;
      }
      return 0;
    };

  private:

    unsigned short _origin;
    std::string    _chstatus_producer;
    std::string    _adcimage_producer;
    std::string    m_match_label;
  };

  /**
     \class larcv::SuperaLArFlowFactory
     \brief A concrete factory class for larcv::SuperaLArFlow
  */
  class SuperaLArFlowProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaLArFlowProcessFactory() { ProcessFactory::get().add_factory("SuperaLArFlow",this); }
    /// dtor
    ~SuperaLArFlowProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaLArFlow(instance_name); }
  };

}
#endif
/** @} */ // end of doxygen group 

