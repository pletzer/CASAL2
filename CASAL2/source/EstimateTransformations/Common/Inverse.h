/**
 * @file Inverse.h
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date Jan 8, 2016
 * @section LICENSE
 *
 * Copyright NIWA Science �2016 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * << Add Description >>
 */
#ifndef SOURCE_ESTIMATETRANSFORMATIONS_CHILDREN_INVERSE_H_
#define SOURCE_ESTIMATETRANSFORMATIONS_CHILDREN_INVERSE_H_

// headers
#include "EstimateTransformations/EstimateTransformation.h"

// namespaces
namespace niwa {
namespace estimatetransformations {

// classes
class Inverse : public EstimateTransformation {
public:
  Inverse() = delete;
  explicit Inverse(Model* model);
  virtual ~Inverse() = default;
  void                        TransformForObjectiveFunction() override final;
  void                        RestoreFromObjectiveFunction() override final;
  std::set<string>            GetTargetEstimates() override final;
  Double                      GetScore() override final;

protected:
  // methods
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoTransform() override final;
  void                        DoRestore() override final;
};

} /* namespace estimatetransformations */
} /* namespace niwa */

#endif /* SOURCE_ESTIMATETRANSFORMATIONS_CHILDREN_INVERSE_H_ */
