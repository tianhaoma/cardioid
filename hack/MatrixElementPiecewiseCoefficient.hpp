#pragma once

#include "mfem.hpp"
#include <unordered_map>
#include <memory>

mfem::DenseMatrix quat2rot(const mfem::Vector& q);

class MatrixElementPiecewiseCoefficient : public mfem::MatrixCoefficient
{
public:
  MatrixElementPiecewiseCoefficient() : mfem::MatrixCoefficient(3) {}
  
  MatrixElementPiecewiseCoefficient(std::shared_ptr<mfem::ParGridFunction> x,
      std::shared_ptr<mfem::ParGridFunction> y,
      std::shared_ptr<mfem::ParGridFunction> z)
  : mfem::MatrixCoefficient(3), p_gf_(x), p_gf_y(y), p_gf_z(z) {}

  virtual void Eval(mfem::DenseMatrix &K,
		    mfem::ElementTransformation& T,
		    const mfem::IntegrationPoint &ip);

  std::shared_ptr<mfem::ParGridFunction> p_gf_;
  std::shared_ptr<mfem::ParGridFunction> p_gf_y;
  std::shared_ptr<mfem::ParGridFunction> p_gf_z;
  std::unordered_map<int,mfem::Vector> heartConductivities_;
  std::unordered_map<int,double> bathConductivities_;
};
