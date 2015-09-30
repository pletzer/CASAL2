/**
 * @file Factory.cpp
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 28/01/2015
 * @section LICENSE
 *
 * Copyright NIWA Science �2014 - www.niwa.co.nz
 *
 */

// headers
#include "Factory.h"

#include "TimeVarying/Children/AnnualShift.h"
#include "TimeVarying/Children/Constant.h"
#include "TimeVarying/Children/Exogenous.h"
#include "TimeVarying/Manager.h"

// namespaces
namespace niwa {
namespace timevarying {

/**
 * Create the instance of our object as defined by the two parameters
 * object_type and sub_type.
 *
 * @param object_type The type of object to create (e.g age_size, process)
 * @param sub_type The child type of the object to create (e.g ageing, schnute)
 * @return shared_ptr to the object we've created
 */
TimeVaryingPtr Factory::Create(const string& object_type, const string& sub_type) {
  TimeVaryingPtr result;

  if (object_type == PARAM_TIME_VARYING) {
    if (sub_type == PARAM_ANNUAL_SHIFT)
      result = TimeVaryingPtr(new AnnualShift());
    else if (sub_type == PARAM_CONSTANT)
      result = TimeVaryingPtr(new Constant());
    else if (sub_type == PARAM_EXOGENOUS)
      result = TimeVaryingPtr(new Exogenous());

    if (result)
      timevarying::Manager::Instance().AddObject(result);
  }

  return result;
}

} /* namespace projects */
} /* namespace niwa */
