#include "pod_coarse_space.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace mfem;

extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
                       double *w, double *work, int *lwork, int *info);

PODCoarseSpace::PODCoarseSpace(ParFiniteElementSpace *fes,
                               PetscParMatrix *A2,
                               const HypreParMatrix &M_mat,
                               const Vector &ones_true,
                               double M_total,
                               MPI_Comm comm)
   : fes_(fes),
     A2_(A2),
     M_mat_(&M_mat),
     ones_true_(ones_true),
     M_total_(M_total),
     comm_(comm),
     k_(0)
{
}

void PODCoarseSpace::GaugeMassWeighted(Vector &v) const
{
   Vector Mv(v.Size());
   M_mat_->Mult(v, Mv);
   const double num = InnerProduct(comm_, ones_true_, Mv);
   v.Add(-num / M_total_, ones_true_);
}

void PODCoarseSpace::SymmetricEigenDescending(DenseMatrix &C,
                                              Vector &sig2,
                                              DenseMatrix &W)
{
   const int M = C.Height();
   sig2.SetSize(M);
   W.SetSize(M, M);
   W = C;

   char jobz = 'V';
   char uplo = 'U';
   int n = M;
   int lda = M;
   int lwork = std::max(1, 4 * M);
   int info = 0;
   Vector work(lwork);
   dsyev_(&jobz, &uplo, &n, W.Data(), &lda, sig2.GetData(),
          work.GetData(), &lwork, &info);
   MFEM_VERIFY(info == 0, "LAPACK dsyev failed while building POD basis.");

   for (int i = 0; i < M / 2; ++i)
   {
      std::swap(sig2(i), sig2(M - 1 - i));
      for (int r = 0; r < M; ++r)
      {
         std::swap(W(r, i), W(r, M - 1 - i));
      }
   }
}

void PODCoarseSpace::BuildBasis(const std::vector<Vector> &U,
                                int k_target,
                                double energy_tol,
                                int k_min)
{
   const int M = static_cast<int>(U.size());
   MFEM_VERIFY(M >= 2, "Need at least two snapshots to build a POD basis.");
   MFEM_VERIFY(k_target >= 1, "POD target rank must be positive.");

   int my_rank = 0;
   MPI_Comm_rank(comm_, &my_rank);

   const int n_loc = U[0].Size();
   std::vector<Vector> AU;
   AU.reserve(M);
   for (int m = 0; m < M; ++m)
   {
      MFEM_VERIFY(U[m].Size() == n_loc, "POD snapshot size mismatch.");
      AU.emplace_back(n_loc);
      A2_->Mult(U[m], AU[m]);
   }

   DenseMatrix C(M, M);
   C = 0.0;
   for (int l = 0; l < M; ++l)
   {
      for (int m = l; m < M; ++m)
      {
         const double cij = InnerProduct(comm_, U[l], AU[m]);
         C(l, m) = cij;
         C(m, l) = cij;
      }
   }

   if (my_rank == 0)
   {
      std::cout << "[POD] Correlation matrix C (" << M << "x" << M
                << ") assembled." << std::endl;
   }

   Vector sig2;
   DenseMatrix W;
   SymmetricEigenDescending(C, sig2, W);

   double total = 0.0;
   for (int j = 0; j < M; ++j)
   {
      total += std::max(sig2(j), 0.0);
   }
   MFEM_VERIFY(total > 0.0, "POD snapshot energy is zero.");

   double accum = 0.0;
   int k = 0;
   const int k_cap = std::min(k_target, M);
   for (; k < k_cap; ++k)
   {
      accum += std::max(sig2(k), 0.0);
      if ((accum / total) >= energy_tol)
      {
         ++k;
         break;
      }
   }
   if (k == 0)
   {
      k = k_cap;
      for (int j = 0; j < k; ++j)
      {
         accum += std::max(sig2(j), 0.0);
      }
   }
   k_ = std::min(std::max(k, std::max(1, k_min)), k_cap);

   if (my_rank == 0)
   {
      std::cout << "[POD] Top eigenvalues:";
      for (int j = 0; j < std::min(10, M); ++j)
      {
         std::cout << " " << std::scientific << std::setprecision(3)
                   << sig2(j);
      }
      std::cout << std::defaultfloat << std::endl;
      std::cout << "[POD] selected k = " << k_
                << " (captured energy = "
                << std::setprecision(8) << (accum / total * 100.0)
                << " %)" << std::endl;
   }

   Phi_.assign(k_, Vector(n_loc));
   for (int j = 0; j < k_; ++j)
   {
      Phi_[j] = 0.0;
      const double sigma = std::sqrt(std::max(sig2(j), 1e-30));
      const double inv_sigma = 1.0 / sigma;
      for (int m = 0; m < M; ++m)
      {
         Phi_[j].Add(inv_sigma * W(m, j), U[m]);
      }
      GaugeMassWeighted(Phi_[j]);
   }

   auto check_basis = [&](double &max_diag_err, double &max_off,
                          double &max_mass_gauge) {
      max_diag_err = 0.0;
      max_off = 0.0;
      max_mass_gauge = 0.0;
      for (int j = 0; j < k_; ++j)
      {
         Vector Aphi(n_loc);
         A2_->Mult(Phi_[j], Aphi);
         for (int l = 0; l <= j; ++l)
         {
            const double a0 = InnerProduct(comm_, Phi_[l], Aphi);
            const double target = (l == j) ? 1.0 : 0.0;
            const double err = std::abs(a0 - target);
            if (l == j)
            {
               max_diag_err = std::max(max_diag_err, err);
            }
            else
            {
               max_off = std::max(max_off, err);
            }
         }

         Vector Mphi(n_loc);
         M_mat_->Mult(Phi_[j], Mphi);
         max_mass_gauge = std::max(max_mass_gauge,
                                   std::abs(InnerProduct(comm_, ones_true_, Mphi)));
      }
   };

   double max_diag_err = 0.0;
   double max_off = 0.0;
   double max_mass_gauge = 0.0;
   check_basis(max_diag_err, max_off, max_mass_gauge);

   if (my_rank == 0)
   {
      std::cout << "[POD] A0 - I: max diag err = " << std::scientific
                << max_diag_err << ", max off-diag = " << max_off
                << ", max |1^T M phi| = " << max_mass_gauge
                << std::defaultfloat << std::endl;
   }

   if (max_off > 1e-8 || max_diag_err > 1e-8)
   {
      if (my_rank == 0)
      {
         std::cout << "[POD] applying A2 modified Gram-Schmidt correction."
                   << std::endl;
      }
      for (int j = 0; j < k_; ++j)
      {
         for (int l = 0; l < j; ++l)
         {
            Vector Aphi_j(n_loc);
            A2_->Mult(Phi_[j], Aphi_j);
            const double coef = InnerProduct(comm_, Phi_[l], Aphi_j);
            Phi_[j].Add(-coef, Phi_[l]);
         }
         Vector Aphi_j(n_loc);
         A2_->Mult(Phi_[j], Aphi_j);
         const double nrm2 = InnerProduct(comm_, Phi_[j], Aphi_j);
         MFEM_VERIFY(nrm2 > 0.0, "A2-MGS produced a zero POD mode.");
         Phi_[j] *= (1.0 / std::sqrt(nrm2));
         GaugeMassWeighted(Phi_[j]);
      }
      check_basis(max_diag_err, max_off, max_mass_gauge);
      if (my_rank == 0)
      {
         std::cout << "[POD] after A2-MGS: max diag err = "
                   << std::scientific << max_diag_err
                   << ", max off-diag = " << max_off
                   << ", max |1^T M phi| = " << max_mass_gauge
                   << std::defaultfloat << std::endl;
      }
   }

   if (my_rank == 0)
   {
      std::ofstream f("pod_diagnostics_rank0.log");
      f << "M_snap = " << M << "\n";
      f << "k = " << k_ << "\n";
      f << "energy_fraction = " << (accum / total) << "\n";
      if (k_ < M)
      {
         f << "tau_k = " << (sig2(k_) / sig2(0)) << "\n";
      }
      f << "max_diag_err = " << max_diag_err << "\n";
      f << "max_off_diag = " << max_off << "\n";
      f << "max_mass_gauge = " << max_mass_gauge << "\n";
      f << "eigenvalues:\n";
      for (int j = 0; j < M; ++j)
      {
         f << j << " " << std::setprecision(16) << sig2(j) << "\n";
      }
   }
}

void PODCoarseSpace::ApplyCoarse(const Vector &r, Vector &z) const
{
   if (z.Size() != r.Size())
   {
      z.SetSize(r.Size());
   }
   z = 0.0;

   for (int j = 0; j < k_; ++j)
   {
      const double coeff = InnerProduct(comm_, Phi_[j], r);
      z.Add(coeff, Phi_[j]);
   }
}
