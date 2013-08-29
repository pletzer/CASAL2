/**
 * @file Model.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 23/08/2013
 * @section LICENSE
 *
 * Copyright NIWA Science �2013 - www.niwa.co.nz
 *
 */
#ifndef TESTRESOURCES_CONFIGURATION_MODEL_H_
#define TESTRESOURCES_CONFIGURATION_MODEL_H_

#include <string>

const std::string model_no_init_phases =
"@model\n"
"start_year 1994\n"
"final_year 2008\n"
"min_age 1\n"
"max_age 50\n"
"age_plus t\n"
"time_steps time_step_one\n";

#endif /* MODEL_H_ */