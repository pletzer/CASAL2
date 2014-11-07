/**
 * @file Executor.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 28/08/2014
 * @section LICENSE
 *
 * Copyright NIWA Science �2014 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * This is an ADT that has only 3 methods. It's used to bind
 * objects to time steps or processes for executing afterwards.
 */
#ifndef BASE_EXECUTOR_H_
#define BASE_EXECUTOR_H_

// headers
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>

#include "BaseClasses/Object.h"

// namespaces
namespace isam {
namespace base {

/**
 * Class definition
 */
class Executor : public isam::base::Object, public boost::enable_shared_from_this<Executor> {
public:
  // methods
  Executor() = default;
  virtual ~Executor() = default;

  virtual void                PreExecute()  = 0;
  virtual void                Execute()     = 0;

  // accessors
  boost::shared_ptr<Executor> shared_ptr() { return shared_from_this(); }
};

} /* namespace base */
} /* namespace isam */

typedef boost::shared_ptr<isam::base::Executor> ExecutorPtr;

#endif /* BASE_EXECUTOR_H_ */