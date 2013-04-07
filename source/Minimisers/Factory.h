/**
 * @file Factory.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 28/02/2013
 * @section LICENSE
 *
 * Copyright NIWA Science �2013 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This is the factory class to create minimisers
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
#ifndef MINIMISERS_FACTORY_H_
#define MINIMISERS_FACTORY_H_

// Headers
#include "Minimisers/Minimiser.h"

// Namespaces
namespace isam {
namespace minimisers {

/**
 * Class Definition
 */
class Factory {
public:
  // Methods
  static MinimiserPtr         Create(const string& block_type, const string& object_type);

private:
  // Methods
  Factory() = delete;
  virtual ~Factory() = delete;
};

} /* namespace minimisers */
} /* namespace isam */
#endif /* FACTORY_H_ */