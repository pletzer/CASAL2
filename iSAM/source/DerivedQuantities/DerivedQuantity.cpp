/**
 * @file DerivedQuantity.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 4/07/2013
 * @section LICENSE
 *
 * Copyright NIWA Science �2013 - www.niwa.co.nz
 *
 */
#include "DerivedQuantity.h"

#include "Selectivities/Manager.h"

namespace isam {

/**
 * default constructor
 */
DerivedQuantity::DerivedQuantity() {
  parameters_.RegisterAllowed(PARAM_LABEL);
  parameters_.RegisterAllowed(PARAM_TYPE);
  parameters_.RegisterAllowed(PARAM_TIME_STEP);
  parameters_.RegisterAllowed(PARAM_INITIALIZATION_TIME_STEPS);
  parameters_.RegisterAllowed(PARAM_CATEGORIES);
  parameters_.RegisterAllowed(PARAM_SELECTIVITIES);

  model_ = Model::Instance();
}

/**
 * Validate the parameters defined in the configuration
 * file for this derived quantity.
 */
void DerivedQuantity::Validate() {
  CheckForRequiredParameter(PARAM_LABEL);
  CheckForRequiredParameter(PARAM_TYPE);
  CheckForRequiredParameter(PARAM_TIME_STEPS);
  CheckForRequiredParameter(PARAM_SELECTIVITIES);
  CheckForRequiredParameter(PARAM_CATEGORIES);

  time_step_label_                  = parameters_.Get(PARAM_TIME_STEPS).GetValue<string>();
  initialisation_time_step_labels_  = parameters_.Get(PARAM_INITIALIZATION_TIME_STEPS).GetValues<string>();
  selectivity_labels_               = parameters_.Get(PARAM_SELECTIVITIES).GetValues<string>();
  category_labels_                  = parameters_.Get(PARAM_CATEGORIES).GetValues<string>();

  if (category_labels_.size() != selectivity_labels_.size())
    LOG_ERROR(parameters_.location(PARAM_SELECTIVITIES) << " count (" << selectivity_labels_.size() << ") "
        << " is not the same as the categories count (" << category_labels_.size() << ")");
}

/**
 * Build the run time relationships between the derived
 * quantity and other components in the model
 */
void DerivedQuantity::Build() {
  partition_  = accessor::CategoriesPtr(new accessor::Categories(category_labels_));

  selectivities::Manager& selectivity_manager = selectivities::Manager::Instance();
  for (string label : selectivity_labels_) {
    selectivities_.push_back(selectivity_manager.GetSelectivity(label));
  }
}

/**
 * Check if this derived quantity has the initialisation phase label
 * assigned too it. This is used to determine if we need to
 * bind this derived quantity to an initialisation phase
 * for calculation
 *
 * @param label The label of the initialisation time step to check
 * @return True if assigned to parameter initialisation time step, false otherwise
 */
bool DerivedQuantity::IsAssignedToInitialisationPhase(const string& label) {
  return std::find(initialisation_time_step_labels_.begin(), initialisation_time_step_labels_.end(), label)
    != initialisation_time_step_labels_.end();
}

/**
 * Return the calculated value stored in this derived quantity
 * for the parameter year. If the year does not exist as a standard
 * value we'll calculate how many years to go back in to
 * the initialisation phase values.
 *
 * Note: We cannot go back more than 1 phase. If this condition
 * is triggered we return the first value from the phase instead
 * of going back.
 *
 * @param year The year to get the derived quantity value for.
 */
Double DerivedQuantity::GetValue(unsigned year) {
  if (values_.find(year) != values_.end())
    return values_[year];

  if (initialisation_values_.size() == 0)
    LOG_ERROR("Trying to get a value from derived quantity " << label_ << " when no value have been calculated");

  // Calculate how many years to go back. At this point
  // either we're in the init phases or we're going back
  // in to the init phases.
  unsigned years_to_go_back = year - model_->start_year();

  Double result = 0.0;
  if (initialisation_values_.rbegin()->size() > years_to_go_back)
    result = (*initialisation_values_.rbegin()->rbegin() + (years_to_go_back - 1));
  else
    result = (*initialisation_values_.rbegin()->begin()); // first value of last init phase

  return result;
}

} /* namespace isam */