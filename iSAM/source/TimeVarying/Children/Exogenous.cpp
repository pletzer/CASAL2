/**
 * @file Exogenous.cpp
 * @author C.Marsh
 * @github https://github.com/Zaita
 * @date 29/09/2015
 * @section LICENSE
 *
 * Copyright NIWA Science �2014 - www.niwa.co.nz
 *
 */

// headers
#include "Exogenous.h"

#include "Utilities/Map.h"
#include "Model/Objects.h"

// namespaces
namespace niwa {
namespace timevarying {

/**
 * Default constructor
 */
Exogenous::Exogenous() : Exogenous(Model::Instance()) { }
Exogenous::Exogenous(ModelPtr model) : TimeVarying(model) {
  parameters_.Bind<Double>(PARAM_A, &a_, "Shift parameter", "");
  parameters_.Bind<Double>(PARAM_EXOGENOUS_VARIABLE, &exogenous_, "Values of exogeneous variable for each year", "");

  RegisterAsEstimable(PARAM_A, &a_);
}

/**
 *    Validate parameters from config file
 */
void Exogenous::DoValidate() {
  if (years_.size() != exogenous_.size())
    LOG_ERROR_P(PARAM_YEARS) << " provided (" << years_.size() << ") does not match the number of values provided (" << exogenous_.size() << ")";

  // Check that the paramter is of type scalar
  string error = "";
  Estimable::Type estimable_type = model_->objects().GetEstimableType(parameter_, error);
  if (estimable_type != Estimable::kSingle)
    LOG_ERROR_P(PARAM_PARAMETER) << "Parameter must be a scalar, other estimable types not supported yet";
}

/**
 *  Calculate mean of exogenous variable, then using the shift paramter. Store the parameters in the map
 *  values_by_year.
 */
void Exogenous::DoBuild() {
  string error = "";
  values_by_year_ = utilities::Map::create(years_, exogenous_);
  Double* param = model_->objects().GetEstimable(parameter_, error);
  LOG_FINEST() << "Parameter value = " << (*param);
  Double total = 0.0;

  for (Double value : exogenous_)
    total += value;

  Double mean = total / exogenous_.size();

  for (unsigned year : years_)
    parameter_by_year_[year] = (*param) + (a_ * (values_by_year_[year] - mean));
}

/**
 *
 */
void Exogenous::DoReset() {

}

/**
 *
 */
void Exogenous::DoUpdate() {
  LOG_FINE() << "Setting Value to: " << parameter_by_year_[Model::Instance()->current_year()];
  (this->*update_function_)(parameter_by_year_[Model::Instance()->current_year()]);
}

} /* namespace timevarying */
} /* namespace niwa */
