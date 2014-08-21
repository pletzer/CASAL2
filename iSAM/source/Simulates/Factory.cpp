/**
 * @file Factory.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 8/08/2014
 * @section LICENSE
 *
 * Copyright NIWA Science �2014 - www.niwa.co.nz
 *
 */

// headers
#include "Factory.h"

#include "Simulates/Children/Constant.h"
#include "Simulates/Manager.h"

// namespaces
namespace isam {
namespace simulates {

/**
 * Create the instance of our object as defined by the two parameters
 * object_type and sub_type.
 *
 * @param object_type The type of object to create (e.g age_size, process)
 * @param sub_type The child type of the object to create (e.g ageing, schnute)
 * @return shared_ptr to the object we've created
 */
SimulatePtr Factory::Create(const string& object_type, const string& sub_type) {
  SimulatePtr result;

  if (object_type == PARAM_PROJECTS || object_type == PARAM_PROJECT) {
    if (sub_type == PARAM_CONSTANT)
      result = SimulatePtr(new Constant());

    if (result)
      simulates::Manager::Instance().AddObject(result);
  }

  return result;
}

} /* namespace simulates */
} /* namespace isam */