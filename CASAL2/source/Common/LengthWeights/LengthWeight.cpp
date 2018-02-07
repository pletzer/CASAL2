/**
 * @file LengthWeight.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 24/07/2013
 * @section LICENSE
 *
 * Copyright NIWA Science �2013 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * << Add Description >>
 */

// headers
#include "LengthWeight.h"
#include "Common/Model/Model.h"

// namespaces
namespace niwa {

/**
 * default constructor
 */
LengthWeight::LengthWeight(Model* model)
: model_(model) {
  parameters_.Bind<string>(PARAM_LABEL, &label_, "The label of the length-weight relationship", "");
parameters_.Bind<string>(PARAM_TYPE, &type_, "The type of the length-weight relationship", "");
}

/**
 *
 */
void LengthWeight::Validate() {
  parameters_.Populate(model_);
  DoValidate();
}

} /* namespace niwa */
