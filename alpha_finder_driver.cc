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
#include <falaise/snemo/datamodels/calibrated_tracker_hit.h>
#include <falaise/snemo/geometry/locator_plugin.h>
#include <falaise/snemo/geometry/gg_locator.h>

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

     // Get geometry locator plugin
      const geomtools::manager & geo_mgr = get_geometry_manager();
      std::string locator_plugin_name;
      if (setup_.has_key("locator_plugin_name")) {
        locator_plugin_name = setup_.fetch_string("locator_plugin_name");
      } else {
        // If no locator plugin name is set, then search for the first one
        const geomtools::manager::plugins_dict_type & plugins = geo_mgr.get_plugins();
        for (geomtools::manager::plugins_dict_type::const_iterator ip = plugins.begin();
             ip != plugins.end();
             ip++) {
          const std::string & plugin_name = ip->first;
          if (geo_mgr.is_plugin_a<snemo::geometry::locator_plugin>(plugin_name)) {
            DT_LOG_DEBUG(get_logging_priority(), "Find locator plugin with name = " << plugin_name);
            locator_plugin_name = plugin_name;
            break;
          }
        }
      }
      // Access to a given plugin by name and type :
      DT_THROW_IF(! geo_mgr.has_plugin(locator_plugin_name) ||
                  ! geo_mgr.is_plugin_a<snemo::geometry::locator_plugin>(locator_plugin_name),
                  std::logic_error,
                  "Found no locator plugin named '" << locator_plugin_name << "'");
      _locator_plugin_ = &geo_mgr.get_plugin<snemo::geometry::locator_plugin>(locator_plugin_name);

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
      DT_LOG_TRACE(get_logging_priority(), "coucou");
      this-> _find_delayed_unfitted_cluster_(tracker_trajectory_data_, particle_track_data_);

      DT_LOG_TRACE(get_logging_priority(), "Exiting.");
      return;
    }
    
    
    
    
		void alpha_finder_driver::_find_delayed_unfitted_cluster_(const snemo::datamodel::tracker_trajectory_data & tracker_trajectory_data_,
     																													snemo::datamodel::particle_track_data & particle_track_data_)
    {
      DT_LOG_TRACE(get_logging_priority(), "Entering...");
			
			// check if the solution exist
			if(! tracker_trajectory_data_.has_solutions())
			{
				DT_LOG_DEBUG(get_logging_priority(), "No solutions have been found !")
				return;
			}
			
			// Store first/last Geiger hits of prompt clusters
			typedef std::pair<snemo::datamodel::calibrated_tracker_hit::handle_type,
												snemo::datamodel::calibrated_tracker_hit::handle_type> gg_pair_type;
			typedef std::vector<gg_pair_type> gg_pair_collection_type;
			gg_pair_collection_type ggpct;
			
			// the default solution is chosen
			const snemo::datamodel::tracker_trajectory_solution & a_solution = tracker_trajectory_data_.get_default_solution();
			const snemo::datamodel::tracker_trajectory_solution::trajectory_col_type & the_trajectories = a_solution.get_trajectories();
				
			// loop on all the trajectories
			for(snemo::datamodel::tracker_trajectory_solution::trajectory_col_type::const_iterator itraj = the_trajectories.begin(); itraj!= the_trajectories.end(); ++itraj)
			{
				const snemo::datamodel::tracker_trajectory & a_trajectory = itraj->get();
				if(! a_trajectory.has_cluster())
				{
					DT_LOG_DEBUG(get_logging_priority(), "No trajectory have been found !")
					continue;
				}
				const snemo::datamodel::tracker_cluster & a_cluster = a_trajectory.get_cluster();
				const snemo::datamodel::calibrated_tracker_hit::collection_type & the_gg_hits = a_cluster.get_hits();
				
				gg_pair_type gg_pair;
				uint32_t layer_min = +std::numeric_limits<uint32_t>::infinity();						
				uint32_t layer_max = -std::numeric_limits<uint32_t>::infinity();						
				// loop on all the geiger hits of the fitted cluster
				for(snemo::datamodel::calibrated_tracker_hit::collection_type::const_iterator igg = the_gg_hits.begin(); igg!= the_gg_hits.end(); ++igg)
				{
					const snemo::datamodel::calibrated_tracker_hit & a_gg_hit = igg->get();
	        const geomtools::geom_id & a_gid = a_gg_hit.get_geom_id();

					// Extract layer
        	const snemo::geometry::gg_locator & gg_locator = _locator_plugin_->get_gg_locator();
        	const uint32_t layer = gg_locator.extract_layer(a_gid);
 					DT_LOG_TRACE(get_logging_priority(), "Layer of prompt Geiger hit = " << layer);
//					double x_wire = a_gg_hit.get_x();
//					double y_wire = a_gg_hit.get_y();
				}				
			}
				
				const snemo::datamodel::tracker_trajectory_solution::cluster_col_type & the_unfitted_clusters = a_solution.get_unfitted_clusters ();
				// loop on all the unfitted cluster
				for(snemo::datamodel::tracker_trajectory_solution::cluster_col_type::const_iterator iclus = the_unfitted_clusters.begin(); iclus!= the_unfitted_clusters.end();
				++iclus)
				{
					const snemo::datamodel::tracker_cluster & a_delayed_cluster = iclus->get();
					
					// the unfitted cluster has to be delayed
					if(! a_delayed_cluster.is_delayed())
					{
						DT_LOG_DEBUG(get_logging_priority(), "The unfitted cluster is not delayed !")
						return;
					}
					const snemo::datamodel::calibrated_tracker_hit::collection_type & delayed_gg_hits = a_delayed_cluster.get_hits();
				
					// loop on all the geiger hits of the unfitted cluster
					for(snemo::datamodel::calibrated_tracker_hit::collection_type::const_iterator idelayedgg = delayed_gg_hits.begin(); idelayedgg!= delayed_gg_hits.end();
					++idelayedgg)
					{
						const snemo::datamodel::calibrated_tracker_hit & a_delayed_gg_hit = idelayedgg->get();
						double x_delayed_wire = a_delayed_gg_hit.get_x();
						double y_delayed_wire = a_delayed_gg_hit.get_y();
					}			
				}
		DT_LOG_TRACE(get_logging_priority(), "Exiting...");
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
