/**
 * @file Manager.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 21/12/2012
 * @section LICENSE
 *
 * Copyright NIWA Science �2012 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * The time class represents a moment of time.
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */
#ifndef MANAGER_H_
#define MANAGER_H_

// Headers
#include "BaseClasses/Manager.h"
#include "Selectivities/Selectivity.h"

// Namespaces
namespace isam {
namespace selectivities {

/**
 * Class defintiion
 */
class Manager : public isam::base::Manager<isam::selectivities::Manager, isam::Selectivity> {
public:
  Manager();
  virtual ~Manager() noexcept(true) {};
};

} /* namespace selectivities */
} /* namespace isam */
#endif /* MANAGER_H_ */