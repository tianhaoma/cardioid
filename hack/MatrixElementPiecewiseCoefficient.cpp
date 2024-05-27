#include "MatrixElementPiecewiseCoefficient.hpp"
#include <cassert>

mfem::DenseMatrix quat2rot(const mfem::Vector &q)
{
   MFEM_ASSERT(q.Size() == 4, "quat2rot: Dimension of quaternion should be 4");
   // MFEM_ASSERT(vecisnorm(q), "quat2rot: quaternion is not normalized");
   mfem::DenseMatrix Q(3);

   double w = q(0);
   double x = q(1);
   double y = q(2);
   double z = q(3);

   double x2 = x * x;
   double y2 = y * y;
   double z2 = z * z;
   double xy = x * y;
   double xz = x * z;
   double yz = y * z;
   double wx = w * x;
   double wy = w * y;
   double wz = w * z;

   Q(0, 0) = 1 - 2 * y2 - 2 * z2;
   Q(1, 0) = 2 * xy - 2 * wz;
   Q(2, 0) = 2 * xz + 2 * wy;

   Q(0, 1) = 2 * xy + 2 * wz;
   Q(1, 1) = 1 - 2 * x2 - 2 * z2;
   Q(2, 1) = 2 * yz - 2 * wx;

   Q(0, 2) = 2 * xz - 2 * wy;
   Q(1, 2) = 2 * yz + 2 * wx;
   Q(2, 2) = 1 - 2 * x2 - 2 * y2;

   return Q;
}

void MatrixElementPiecewiseCoefficient::Eval(
    mfem::DenseMatrix &K,
    mfem::ElementTransformation &T,
    const mfem::IntegrationPoint &ip)
{
   std::unordered_map<int, mfem::Vector>::iterator iter = heartConductivities_.find(T.Attribute);
   if (iter != heartConductivities_.end())
   {
      mfem::Vector direction(3);
      mfem::Vector sheet_direction(3);
      mfem::Vector transverse_direction(3);
      if (1)
      {
         p_gf_->GetVectorValue(T.ElementNo, ip, direction);
         p_gf_y->GetVectorValue(T.ElementNo, ip, sheet_direction);
         p_gf_z->GetVectorValue(T.ElementNo, ip, transverse_direction);
      }
      else
      {
         direction = 0.0;
      }

      // check normalization here
      double mod = 0.0;
      double mod_y = 0.0;
      double mod_z = 0.0;

      for (int i = 0; i < 3; i++)
      {
         mod += direction(i) * direction(i);
         mod_y += sheet_direction(i) * sheet_direction(i);
         mod_z += transverse_direction(i) * transverse_direction(i);
      }

      double mag = std::sqrt(mod);
      double mag_y = std::sqrt(mod_y);
      double mag_z = std::sqrt(mod_z);

      if ((mag < 1e-15) && (mag_y < 1e-15) && (mag_z < 1e-15))
      {
         // std::cout<<"check mag err!\n";
         throw "WARNING: The direction of this fiber/sheet/transverse is a zero vector!";
         // std::logic_error("The direction is a zero vector");
      }

      if ((abs(mag - 1) > 1e-15) || (abs(mag_y - 1) > 1e-15) || (abs(mag_z - 1) > 1e-15))
      {
         for (int i = 0; i < 3; i++)
         {
            direction(i) = direction(i) / mag;
            sheet_direction(i) = sheet_direction(i) / mag_y;
            transverse_direction(i) = transverse_direction(i) / mag_z;
         }
      }
      // normalize here END

      // quaternion transformation

      /*
      mfem::Vector quat(4);
      double w2 = 1;
      for (int ii=0; ii<3; ii++) {
         quat(ii+1) = direction(ii);
         w2 -= direction(ii)*direction(ii);
      }
      quat(0) = sqrt(w2);

      mfem::DenseMatrix VVV = quat2rot(quat);
      MultADAt(VVV,iter->second,K);
      */

      mfem::DenseMatrix K0(3);
      mfem::DenseMatrix K1(3);
      mfem::DenseMatrix K2(3);

      for (int ii = 0; ii < 3; ii++)
      {
         for (int jj = 0; jj < 3; jj++)
         {
            K0(ii, jj) = direction(ii) * direction(jj);
         }
      }
      for (int ii = 0; ii < 3; ii++)
      {
         for (int jj = 0; jj < 3; jj++)
         {
            K1(ii, jj) = sheet_direction(ii) * sheet_direction(jj);
         }
      }
      for (int ii = 0; ii < 3; ii++)
      {
         for (int jj = 0; jj < 3; jj++)
         {
            K2(ii, jj) = transverse_direction(ii) * transverse_direction(jj);
         }
      }

      mfem::Vector hc(3);
      hc = iter->second;

      for (int ii = 0; ii < 3; ii++)
      {
         for (int jj = 0; jj < 3; jj++)
         {
            K0(ii, jj) = hc(0) * K0(ii, jj);
         }
      }
      for (int ii = 0; ii < 3; ii++)
      {
         for (int jj = 0; jj < 3; jj++)
         {
            K1(ii, jj) = hc(1) * K1(ii, jj);
         }
      }
      for (int ii = 0; ii < 3; ii++)
      {
         for (int jj = 0; jj < 3; jj++)
         {
            K2(ii, jj) = hc(2) * K2(ii, jj);
         }
      }

      for (int ii = 0; ii < 3; ii++)
      {
         for (int jj = 0; jj < 3; jj++)
         {
            K(ii, jj) = K0(ii, jj) + K1(ii, jj) + K2(ii, jj);
         }
      }
   }
   else
   {
      std::unordered_map<int, double>::iterator iter = bathConductivities_.find(T.Attribute);
      assert(iter != bathConductivities_.end());
      K = 0.0;
      for (int ii = 0; ii < 3; ii++)
      {
         K(ii, ii) = iter->second;
      }
   }
}
