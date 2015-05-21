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
#include <falaise/snemo/datamodels/tracker_clustering_solution.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/datamodels/line_trajectory_pattern.h>
#include <falaise/snemo/datamodels/helix_trajectory_pattern.h>
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

      // Minimal distance in XY coordinate between Geiger hits
      if (setup_.has_key("minimal_cluster_xy_distance")) {
        _minimal_cluster_xy_search_distance_ = setup_.fetch_real("minimal_cluster_xy_distance");
        if (! setup_.has_explicit_unit("minimal_cluster_xy_distance")) {
          _minimal_cluster_xy_search_distance_ *= CLHEP::cm;
        }
      }

      // Minimal distance in Z coordinate between Geiger hits
      if (setup_.has_key("minimal_cluster_z_distance")) {
        _minimal_cluster_z_search_distance_ = setup_.fetch_real("minimal_cluster_z_distance");
        if (! setup_.has_explicit_unit("minimal_cluster_z_distance")) {
          _minimal_cluster_z_search_distance_ *= CLHEP::cm;
        }
      }

      // Minimal distance between vertex of the prompt track and delayed GG hit
      if (setup_.has_key("minimal_vertex_distance")) {
        _minimal_vertex_distance_ = setup_.fetch_real("minimal_vertex_distance");
        if (! setup_.has_explicit_unit("minimal_vertex_distance")) {
          _minimal_vertex_distance_ *= CLHEP::cm;
        }
      }

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

      _minimal_cluster_xy_search_distance_ = 21 * CLHEP::cm;
      _minimal_cluster_z_search_distance_  = 30 * CLHEP::cm;
      _minimal_vertex_distance_            = 30 * CLHEP::cm;
      return;
    }

    void alpha_finder_driver::process(const snemo::datamodel::tracker_trajectory_data & tracker_trajectory_data_,
                                      snemo::datamodel::particle_track_data & particle_track_data_)
    {
      DT_LOG_TRACE(get_logging_priority(), "Entering...");
      DT_THROW_IF(! is_initialized(), std::logic_error, "Driver is not initialized !");

      this-> _find_delayed_unfitted_cluster_(tracker_trajectory_data_, particle_track_data_);
      this-> _find_delayed_unclustered_hit_(tracker_trajectory_data_, particle_track_data_);

      DT_LOG_TRACE(get_logging_priority(), "Exiting.");
      return;
    }

    void alpha_finder_driver::_find_delayed_unfitted_cluster_(const snemo::datamodel::tracker_trajectory_data & tracker_trajectory_data_,
                                                              snemo::datamodel::particle_track_data & particle_track_data_)
    {
      DT_LOG_TRACE(get_logging_priority(), "Entering...");

      // Check if the solution exist
      if (! tracker_trajectory_data_.has_solutions()) {
        DT_LOG_DEBUG(get_logging_priority(), "No solutions have been found !");
        return;
      }

      // The default solution is chosen
      const snemo::datamodel::tracker_trajectory_solution & a_solution = tracker_trajectory_data_.get_default_solution();
      // Get unfitted cluster first
      if (! a_solution.has_unfitted_clusters()) {
        DT_LOG_DEBUG(get_logging_priority(), "No unfitted clusters has been found !");
        return;
      }
      const snemo::datamodel::tracker_trajectory_solution::cluster_col_type & the_unfitted_clusters
        = a_solution.get_unfitted_clusters();
      // Loop on all the unfitted cluster
      for (snemo::datamodel::tracker_trajectory_solution::cluster_col_type::const_iterator
             iclus = the_unfitted_clusters.begin();
           iclus!= the_unfitted_clusters.end(); ++iclus) {
        const snemo::datamodel::tracker_cluster & a_delayed_cluster = iclus->get();

        // The unfitted cluster has to be delayed
        if (! a_delayed_cluster.is_delayed()) {
          DT_LOG_DEBUG(get_logging_priority(), "The unfitted cluster is not delayed !");
          return;
        }
        const snemo::datamodel::calibrated_tracker_hit::collection_type & delayed_gg_hits
          = a_delayed_cluster.get_hits();

        geomtools::vector_3d associated_vertex;
        geomtools::invalidate(associated_vertex);
        this->_find_short_track_(delayed_gg_hits, a_solution, associated_vertex);
        if (! geomtools::is_valid(associated_vertex)) {
          DT_LOG_DEBUG(get_logging_priority(), "No association with prompt cluster has been found !");
          continue;
        }

        // Add short alpha particle track
        snemo::datamodel::particle_track::handle_type hPT(new snemo::datamodel::particle_track);
        snemo::datamodel::particle_track & a_short_alpha = hPT.grab();
        a_short_alpha.set_track_id(particle_track_data_.get_number_of_particles());
        a_short_alpha.set_charge(snemo::datamodel::particle_track::undefined);
        particle_track_data_.add_particle(hPT);

        // "Fit" short alpha trajectory
        geomtools::vector_3d other_vertex;
        geomtools::invalidate(other_vertex);
        this->_fit_short_track_(delayed_gg_hits, associated_vertex, other_vertex);
        if (geomtools::is_valid(other_vertex)) {
          // Create new 'tracker_trajectory' handle:
          snemo::datamodel::tracker_trajectory::handle_type htrajectory;
          htrajectory.reset(new snemo::datamodel::tracker_trajectory);
          snemo::datamodel::tracker_trajectory & a_trajectory = htrajectory.grab();
          a_trajectory.set_cluster_handle(*iclus);
          snemo::datamodel::line_trajectory_pattern * a_line
            = new snemo::datamodel::line_trajectory_pattern;
          a_trajectory.set_pattern_handle(a_line);
          a_line->grab_segment().set_first(associated_vertex);
          a_line->grab_segment().set_last(other_vertex);
          a_short_alpha.set_trajectory_handle(htrajectory);
        }
      } // end of delayed cluster

      DT_LOG_TRACE(get_logging_priority(), "Exiting...");
      return;
    }

    void alpha_finder_driver::_find_delayed_unclustered_hit_(const snemo::datamodel::tracker_trajectory_data & tracker_trajectory_data_,
                                                             snemo::datamodel::particle_track_data & particle_track_data_)
    {
      DT_LOG_TRACE(get_logging_priority(), "Entering...");

      // Check if the solution exist
      if (! tracker_trajectory_data_.has_solutions()) {
        DT_LOG_DEBUG(get_logging_priority(), "No solutions have been found !");
        return;
      }
      // The default solution is chosen
      const snemo::datamodel::tracker_trajectory_solution & a_solution = tracker_trajectory_data_.get_default_solution();
      // Get associated clustering solution (if any)
      if (! a_solution.has_clustering_solution()) {
        DT_LOG_DEBUG(get_logging_priority(), "No clustering solution has been found !");
        return;
      }
      const snemo::datamodel::tracker_clustering_solution & a_clustering_solution
        = a_solution.get_clustering_solution();
      if (! a_clustering_solution.has_unclustered_hits()) {
        DT_LOG_DEBUG(get_logging_priority(), "No unclustered hits !");
        return;
      }
      const snemo::datamodel::calibrated_tracker_hit::collection_type & unclustered_gg_hits
        = a_clustering_solution.get_unclustered_hits();

      geomtools::vector_3d associated_vertex;
      geomtools::invalidate(associated_vertex);
      this->_find_short_track_(unclustered_gg_hits, a_solution, associated_vertex);
      if (! geomtools::is_valid(associated_vertex)) {
        DT_LOG_DEBUG(get_logging_priority(), "No association with prompt cluster has been found !");
        return;
      }

      // Add short alpha particle track
      snemo::datamodel::particle_track::handle_type hPT(new snemo::datamodel::particle_track);
      snemo::datamodel::particle_track & a_short_alpha = hPT.grab();
      a_short_alpha.set_track_id(particle_track_data_.get_number_of_particles());
      a_short_alpha.set_charge(snemo::datamodel::particle_track::undefined);
      particle_track_data_.add_particle(hPT);

      // "Fit" short alpha trajectory
      geomtools::vector_3d other_vertex;
      geomtools::invalidate(other_vertex);
      this->_fit_short_track_(unclustered_gg_hits, associated_vertex, other_vertex);
      if (geomtools::is_valid(other_vertex)) {
        // Create new 'tracker_trajectory' handle:
        snemo::datamodel::tracker_trajectory::handle_type htrajectory;
        htrajectory.reset(new snemo::datamodel::tracker_trajectory);
        snemo::datamodel::tracker_trajectory & a_trajectory = htrajectory.grab();
        snemo::datamodel::line_trajectory_pattern * a_line
          = new snemo::datamodel::line_trajectory_pattern;
        a_trajectory.set_pattern_handle(a_line);
        a_line->grab_segment().set_first(associated_vertex);
        a_line->grab_segment().set_last(other_vertex);
        a_short_alpha.set_trajectory_handle(htrajectory);
      }

      DT_LOG_TRACE(get_logging_priority(), "Exiting...");
      return;
    }

    void alpha_finder_driver::_find_short_track_(const snemo::datamodel::calibrated_tracker_hit::collection_type & hits_,
                                                 const snemo::datamodel::tracker_trajectory_solution & solution_,
                                                 geomtools::vector_3d & associated_vertex_)
    {
      bool has_associated_alpha = false;
      // Loop on all the geiger hits of the unfitted cluster
      for (snemo::datamodel::calibrated_tracker_hit::collection_type::const_iterator
             ihit = hits_.begin(); ihit != hits_.end(); ++ihit) {
        const snemo::datamodel::calibrated_tracker_hit & a_delayed_gg_hit = ihit->get();
        if (! a_delayed_gg_hit.is_delayed()) continue;

        // Get prompt trajectories
        if (! solution_.has_trajectories()) {
          DT_LOG_DEBUG(get_logging_priority(), "No prompt trajectories have been found !");
          return;
        }
        const snemo::datamodel::tracker_trajectory_solution::trajectory_col_type & the_trajectories
          = solution_.get_trajectories();
        // Loop on all the trajectories
        for (snemo::datamodel::tracker_trajectory_solution::trajectory_col_type::const_iterator
               itraj = the_trajectories.begin();
             itraj != the_trajectories.end(); ++itraj) {
          const snemo::datamodel::tracker_trajectory & a_trajectory = itraj->get();
          // Look into properties to find the default trajectory. Here,
          // default means the one with the best chi2. This flag is set by the
          // 'fitting' module.
          if (! a_trajectory.get_auxiliaries().has_flag("default")) continue;

          if (! a_trajectory.has_cluster()) {
            DT_LOG_DEBUG(get_logging_priority(), "No associated clusters have been found !");
            continue;
          }
          const snemo::datamodel::tracker_cluster & a_prompt_cluster = a_trajectory.get_cluster();
          if (a_prompt_cluster.is_delayed()) {
            DT_LOG_DEBUG(get_logging_priority(), "The cluster is not prompt !");
            continue;
          }
          const snemo::datamodel::calibrated_tracker_hit::collection_type & prompt_gg_hits
            = a_prompt_cluster.get_hits();
          for (snemo::datamodel::calibrated_tracker_hit::collection_type::const_iterator
                 ipromptgg = prompt_gg_hits.begin();
               ipromptgg != prompt_gg_hits.end(); ++ipromptgg) {
            const snemo::datamodel::calibrated_tracker_hit & a_prompt_gg_hit = ipromptgg->get();
            const geomtools::vector_2d a_prompt_position(a_prompt_gg_hit.get_x(),
                                                         a_prompt_gg_hit.get_y());
            const geomtools::vector_2d a_delayed_position(a_delayed_gg_hit.get_x(),
                                                          a_delayed_gg_hit.get_y());
            const double distance_xy = (a_delayed_position - a_prompt_position).mag();
            DT_LOG_DEBUG(get_logging_priority(),
                         "Distance between delayed hit " << a_delayed_gg_hit.get_geom_id()
                         << " and prompt hit " << a_prompt_gg_hit.get_geom_id() << " = "
                         << distance_xy/CLHEP::cm << " cm");
            bool has_minimal_xy_search = false;
            if (distance_xy < _minimal_cluster_xy_search_distance_) {
              has_minimal_xy_search = true;
            }

            const double z_delayed       = a_delayed_gg_hit.get_z();
            const double sigma_z_delayed = a_delayed_gg_hit.get_sigma_z();
            const double z_prompt        = a_prompt_gg_hit.get_z();
            const double sigma_z_prompt  = a_prompt_gg_hit.get_sigma_z();
            const double distance_z      = std::abs(z_delayed - z_prompt);
            const double sigma           = sigma_z_delayed + sigma_z_prompt;
            DT_LOG_DEBUG(get_logging_priority(),
                         "Distance between delayed hit " << a_delayed_gg_hit.get_geom_id()
                         << " and prompt hit " << a_prompt_gg_hit.get_geom_id() << " z = "
                         << distance_z/CLHEP::cm << " cm - sigma = " << sigma/CLHEP::cm << " cm");
            bool has_minimal_z_search = false;
            if ((distance_z - _minimal_cluster_z_search_distance_) < sigma) {
              has_minimal_z_search = true;
            }

            if (has_minimal_xy_search && has_minimal_z_search) {
              DT_LOG_DEBUG(get_logging_priority(), "Find a GG candidate close enough to delayed hit");
              has_associated_alpha = true;
              break;
            }
          } // end of prompt gg hits

          // Look for trajectories extremities
          if (has_associated_alpha) {
            DT_LOG_DEBUG(get_logging_priority(), "Look for trajectory extremities");
            geomtools::vector_3d first;
            geomtools::vector_3d last;
            geomtools::invalidate(first);
            geomtools::invalidate(last);
            const snemo::datamodel::base_trajectory_pattern & a_pattern = a_trajectory.get_pattern();
            const std::string & a_pattern_id = a_pattern.get_pattern_id();
            if (a_pattern_id == snemo::datamodel::line_trajectory_pattern::pattern_id()) {
              const snemo::datamodel::line_trajectory_pattern & ltp
                = dynamic_cast<const snemo::datamodel::line_trajectory_pattern &>(a_pattern);
              first = ltp.get_segment().get_first();
              last  = ltp.get_segment().get_last();
            } else if (a_pattern_id == snemo::datamodel::helix_trajectory_pattern::pattern_id()) {
              const snemo::datamodel::helix_trajectory_pattern & htp
                = dynamic_cast<const snemo::datamodel::helix_trajectory_pattern &>(a_pattern);
              first = htp.get_helix().get_first();
              last  = htp.get_helix().get_last();
            }
            const geomtools::vector_3d a_delayed_position(a_delayed_gg_hit.get_x(),
                                                          a_delayed_gg_hit.get_y(),
                                                          a_delayed_gg_hit.get_z());
            if (geomtools::is_valid(first)) {
              const double distance = (first - a_delayed_position).mag();
              if (distance < _minimal_vertex_distance_) {
                associated_vertex_ = first;
              }
            }
            if (geomtools::is_valid(last)) {
              const double distance = (last - a_delayed_position).mag();
              if (distance < _minimal_vertex_distance_) {
                associated_vertex_ = last;
              }
            }
          }
        } // end of trajectories

        if (has_associated_alpha) {
          break;
        }
      } // end of delayed hits
      return;
    }


    void alpha_finder_driver::_fit_short_track_(const snemo::datamodel::calibrated_tracker_hit::collection_type & hits_,
                                                const geomtools::vector_3d & first_vertex_,
                                                geomtools::vector_3d & last_vertex_)
    {
      double max_distance = 0.0 * CLHEP::cm;
      // Loop on all the delayed geiger hits to compute distance between hit
      // and associated vertex
      for (snemo::datamodel::calibrated_tracker_hit::collection_type::const_iterator
             ihit = hits_.begin(); ihit != hits_.end(); ++ihit) {
        const snemo::datamodel::calibrated_tracker_hit & a_hit = ihit->get();
        const geomtools::vector_3d a_hit_position(a_hit.get_x(),
                                                  a_hit.get_y(),
                                                  a_hit.get_z());
        const double distance = (first_vertex_ - a_hit_position).mag();
        if (distance > max_distance) {
          last_vertex_ = a_hit_position;
          max_distance = distance;
        }
      } // end of delayed geiger hits
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
