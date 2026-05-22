#ifndef CARDIOID_HACK_POD_COARSE_SPACE_HPP
#define CARDIOID_HACK_POD_COARSE_SPACE_HPP

#include "mfem.hpp"

#include <mpi.h>
#include <vector>

class PODCoarseSpace
{
public:
   PODCoarseSpace(mfem::ParFiniteElementSpace *fes,
                  mfem::PetscParMatrix *A2,
                  const mfem::HypreParMatrix &M_mat,
                  const mfem::Vector &ones_true,
                  double M_total,
                  MPI_Comm comm);

   void BuildBasis(const std::vector<mfem::Vector> &U_snapshots,
                   int k_target,
                   double energy_tol,
                   int k_min = 1);

   void ApplyCoarse(const mfem::Vector &r, mfem::Vector &z) const;

   int Rank() const { return k_; }
   const mfem::Vector &GetMode(int j) const { return Phi_[j]; }

private:
   void GaugeMassWeighted(mfem::Vector &v) const;
   void SymmetricEigenDescending(mfem::DenseMatrix &C,
                                 mfem::Vector &sig2,
                                 mfem::DenseMatrix &W);

   mfem::ParFiniteElementSpace *fes_;
   mfem::PetscParMatrix *A2_;
   const mfem::HypreParMatrix *M_mat_;
   mfem::Vector ones_true_;
   double M_total_;
   MPI_Comm comm_;

   std::vector<mfem::Vector> Phi_;
   int k_;
};

#endif
