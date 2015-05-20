// \file falaise/snemo/reconstruction/alpha_finder_driver.cc

// Ourselves:
#include <falaise/snemo/reconstruction/alpha_finder_driver.h>

// Standard library:
#include <sstream>

// Third party:
// - Bayeux/datatools:
#include <datatools/properties.h>
#include <datatools/units.h>
// - Bayeux/geomtools
#include <geomtools/manager.h>

// This project (Falaise):
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include <falaise/snemo/datamodels/particle_track_data.h>

namespace snemo {

  namespace reconstruction {

    const std::string & alpha_finder_driver::get_id()
    {
      static const std::string s("AFD");
      return s;
    }

    void alpha_finder_driver::set_initialized(const bool initialized_)
    {
      _initialized_ = initialized_;
      return;
    }

    bool alpha_finder_driver::is_initialized() const
    {
      return _initialized_;
    }

    void alpha_finder_driver::set_logging_priority(const datatools::logger::priority priority_)
    {
      _logging_priority_ = priority_;
      return;
    }

    datatools::logger::priority alpha_finder_driver::get_logging_priority() const
    {
      return _logging_priority_;
    }

    bool alpha_finder_driver::has_geometry_manager() const
    {
      return _geometry_manager_ != 0;
    }

    void alpha_finder_driver::set_geometry_manager(const geomtools::manager & gmgr_)
    {
      DT_THROW_IF(is_initialized(), std::logic_error, "Driver is already initialized !");
      _geometry_manager_ = &gmgr_;
      return;
    }

    const geomtools::manager & alpha_finder_driver::get_geometry_manager() const
    {
      DT_THROW_IF(! has_geometry_manager(), std::logic_error, "No geometry manager is setup !");
      return *_geometry_manager_;
    }

    /// Constructor
    alpha_finder_driver::alpha_finder_driver()
    {
      _set_defaults();
      return;
    }

    /// Destructor
    alpha_finder_driver::~alpha_finder_driver()
    {
      if (is_initialized()) {
        reset();
      }
      return;
    }

    /// Initialize the driver through configuration properties
    void alpha_finder_driver::initialize(const datatools::properties & setup_)
    {
      DT_THROW_IF(is_initialized(), std::logic_error, "Driver is already initialized !");

      DT_THROW_IF(! has_geometry_manager(), std::logic_error, "Missing geometry manager !");
      DT_THROW_IF(! get_geometry_manager().is_initialized(), std::logic_error,
                  "Geometry manager is not initialized !");

      // Logging priority
      datatools::logger::priority lp = datatools::logger::extract_logging_configuration(setup_);
      DT_THROW_IF(lp == datatools::logger::PRIO_UNDEFINED,
                  std::logic_error,
                  "Invalid logging priority level for geometry manager !");
      set_logging_priority(lp);

      // // Matching distance tolerance for calorimeter association
      // if (setup_.has_key("matching_tolerance")) {
      //   _matching_tolerance_ = setup_.fetch_real("matching_tolerance");
      //   if (! setup_.has_explicit_unit("matching_tolerance")) {
      //     _matching_tolerance_ *= CLHEP::mm;
      //   }
      // }

      set_initialized(true);
      return;
    }

    /// Reset the driver
    void alpha_finder_driver::reset()
    {
      _set_defaults();
      return;
    }

    void alpha_finder_driver::_set_defaults()
    {
      _initialized_      = false;
      _logging_priority_ = datatools::logger::PRIO_WARNING;

      // _matching_tolerance_   = 50 * CLHEP::mm;
      return;
    }

    void alpha_finder_driver::process(const snemo::datamodel::tracker_trajectory_data & tracker_trajectory_data_,
                                      snemo::datamodel::particle_track_data & particle_track_data_)
    {
      DT_LOG_TRACE(get_logging_priority(), "Entering...");
      DT_THROW_IF(! is_initialized(), std::logic_error, "Driver is not initialized !");

      // this->_measure_matching_calorimeters_(calorimeter_hits_, particle_);

      DT_LOG_TRACE(get_logging_priority(), "Exiting.");
      return;
    }

    // static
    void alpha_finder_driver::init_ocd(datatools::object_configuration_description & ocd_)
    {

      // Prefix "AFD" stands for "Alpha Finder Driver" :
      datatools::logger::declare_ocd_logging_configuration(ocd_, "fatal", "AFD.");

      // {
      //   // Description of the 'VED.matching_tolerance' configuration property :
      //   datatools::configuration_property_description & cpd
      //     = ocd_.add_property_info();
      //   cpd.set_name_pattern("VED.matching_tolerance")
      //     .set_from("snemo::reconstruction::alpha_finder_driver")
      //     .set_terse_description("Set matching length tolerance between vertex track and calo. block")
      //     .set_traits(datatools::TYPE_REAL)
      //     .set_mandatory(false)
      //     .set_explicit_unit(true)
      //     .set_unit_label("length")
      //     .set_unit_symbol("mm")
      //     .set_default_value_real(50.)
      //     .add_example("Set the default value::                           \n"
      //                  "                                                  \n"
      //                  "  CAD.matching_tolerance : real as length = 50 mm \n"
      //                  "                                                  \n"
      //                  )
      //     ;
      // }
    }

  }  // end of namespace reconstruction

}  // end of namespace snemo

/* OCD support */
#include <datatools/object_configuration_description.h>
DOCD_CLASS_IMPLEMENT_LOAD_BEGIN(snemo::reconstruction::alpha_finder_driver,ocd_)
{
  ocd_.set_class_name("snemo::reconstruction::alpha_finder_driver");
  ocd_.set_class_description("A driver class for finding short alpha track");
  ocd_.set_class_library("Falaise_ChargedParticleTracking");
  ocd_.set_class_documentation("This driver extracts the short alpha tracks. \n"
                               );

  // Invoke specific OCD support :
  ::snemo::reconstruction::alpha_finder_driver::init_ocd(ocd_);

  ocd_.set_validation_support(true);
  ocd_.lock();
  return;
}
DOCD_CLASS_IMPLEMENT_LOAD_END() // Closing macro for implementation
DOCD_CLASS_SYSTEM_REGISTRATION(snemo::reconstruction::alpha_finder_driver,
                               "snemo::reconstruction::alpha_finder_driver")

// end of falaise/snemo/reconstruction/alpha_finder_driver.cc
