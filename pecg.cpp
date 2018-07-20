
#include "mfem.hpp"
#include "object.h"
#include "ddcMalloc.h"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cassert>
#include <memory>
#include <set>

using namespace std;
using namespace mfem;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

DenseMatrix quat2rot(const Vector& q)
{
    MFEM_ASSERT(q.Size()==4, "quat2rot: Dimension of quaternion should be 4");
    //MFEM_ASSERT(vecisnorm(q), "quat2rot: quaternion is not normalized");
    DenseMatrix Q(3);

    double w=q(0);
    double x=q(1);
    double y=q(2);
    double z=q(3);

    double x2=x*x;
    double y2=y*y;
    double z2=z*z;
    double xy=x*y;
    double xz=x*z;
    double yz=y*z;
    double wx=w*x;
    double wy=w*y;
    double wz=w*z;

    Q(0,0)=1-2*y2-2*z2;
    Q(1,0)=2*xy-2*wz;
    Q(2,0)=2*xz+2*wy;

    Q(0,1)=2*xy+2*wz;
    Q(1,1)=1-2*x2-2*z2;
    Q(2,1)=2*yz-2*wx;

    Q(0,2)=2*xz-2*wy;
    Q(1,2)=2*yz+2*wx;
    Q(2,2)=1-2*x2-2*y2;

    return Q;
}

class MatrixElementPiecewiseCoefficient : public MatrixCoefficient
{
 public:
   MatrixElementPiecewiseCoefficient() : MatrixCoefficient(3) {}
   MatrixElementPiecewiseCoefficient(shared_ptr<ParGridFunction> x) : MatrixCoefficient(3) {p_gf_=x;}

   virtual void Eval(DenseMatrix &K, ElementTransformation& T, const IntegrationPoint &ip)
   {
      unordered_map<int,Vector>::iterator iter = heartConductivities_.find(T.Attribute);
      if (iter != heartConductivities_.end())
      {
         Vector direction(3);
         if (0)
         {
            p_gf_->GetVectorValue(T.ElementNo, ip, direction);
         } else {
            direction = 0.0;
         }
         Vector quat(4);
         double w2 = 1;
         for (int ii=0; ii<3; ii++)
         {
            quat(ii+1) = direction(ii);
            w2 -= direction(ii)*direction(ii);
         }
         quat(0) = sqrt(w2);
         
         DenseMatrix VVV = quat2rot(quat);
         MultADAt(VVV,iter->second,K);
      }
      else
      {
         unordered_map<int,double>::iterator iter = bathConductivities_.find(T.Attribute);
         assert(iter != bathConductivities_.end());
         K=0.0;
         for (int ii=0; ii<3; ii++)
         {
            K(ii,ii) = iter->second;
         }
      }
   }

   shared_ptr<ParGridFunction> p_gf_;
   unordered_map<int,Vector> heartConductivities_;
   unordered_map<int,double> bathConductivities_;
};


/** The object_get code doesn't handle empty strings for a default value
 *  very well.  The problem is that object_parse attempts to tokenize
 *  the default value and strtok_r returns a NULL pointer.  This results
 *  in a situation where the pointer passed to object_get (i.e., tmp)
 *  isn't set by object_get.  Fortunately, this can be detected by the
 *  fact that object_get will return 0 instead of 1.  We can catch this
 *  case and do the right thing with the value we return to the caller.
 */
void objectGet(OBJECT* obj, const string& name, string& value, const string& defVal)
{
   char* tmp;
   int nFound = object_get(obj, name.c_str(), &tmp, STRING, 1, defVal.c_str());
   if (nFound != 0)
   {
      value = tmp;
      ddcFree(tmp);
   }
   else
      value = "";
}

void objectGet(OBJECT* obj, const string& name, vector<int>& value)
{
   value.clear();
   int* tmp;
   unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, INT, IGNORE_IF_NOT_FOUND);
   value.reserve(n);
   for (unsigned ii=0; ii<n; ++ii)
      value.push_back(tmp[ii]);
   ddcFree(tmp);
}

void objectGet(OBJECT* obj, const string& name, vector<double>& value)
{
   value.clear();
   double* tmp;
   unsigned n = object_getv(obj, name.c_str(), (void**) &tmp, DOUBLE, IGNORE_IF_NOT_FOUND);
   value.reserve(n);
   for (unsigned ii=0; ii<n; ++ii)
      value.push_back(tmp[ii]);
   ddcFree(tmp);
}

std::set<int> readSet(std::istream& stream) {
   std::set<int> retval;
   while(1)
   {
      int value;
      stream >> value;
      if (!stream) {break;}
      retval.insert(value);
   }
   return retval;
}

int main(int argc, char *argv[])
{
   int num_procs, myid;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   // 1. Parse command-line options.
   int order = 1;

   vector<string> objectFilenames;
   if (argc == 1)
   {
      objectFilenames.push_back("ecg.data");
   }
   for (int iargCursor=1; iargCursor<argc; iargCursor++)
   {
      objectFilenames.push_back(argv[iargCursor]);
   }
   for (int ii=0; ii<objectFilenames.size(); ii++)
   {
      object_compilefile(objectFilenames[ii].c_str());
   }

   OBJECT* obj = object_find("ecg", "ECG");
   assert(obj != NULL);
   string mesh_file;
   objectGet(obj, "mesh", mesh_file, "");

   string groundFilename;
   objectGet(obj,"ground",groundFilename,"");
   std::set<int> ground;
   {
      ifstream groundFile(groundFilename);
      ground = readSet(groundFile);
   }

   ifstream mesh_ifs(mesh_file.c_str());
   ParMesh *mesh = new ParMesh(MPI_COMM_WORLD,mesh_ifs);
   int dim = mesh->Dimension();

   {
      int nbe=mesh->GetNBE();
      for(int i=0; i<nbe; i++){
         Element *ele = mesh->GetBdrElement(i);
         const int *v = ele->GetVertices();
         const int nv = ele->GetNVertices();
         bool isGround = true;
         for( int ivert=0; ivert<nv; ivert++) {
            if (ground.find(v[ivert]) == ground.end()) {
               isGround = false;
               break;
            }
         }
         if (isGround) {
            ele->SetAttribute(2);
         } else {
            ele->SetAttribute(1);
         }
      }
   }
   mesh->SetAttributes();
   
   FiniteElementCollection *fec;
   if (mesh->GetNodes())
   {
      fec = mesh->GetNodes()->OwnFEC();
      cout << "Using isoparametric FEs: " << fec->Name() << endl;
   }
   else
   {
      fec = new H1_FECollection(order, dim);
   }
   ParFiniteElementSpace *fespace = new ParFiniteElementSpace(mesh, fec);
   HYPRE_Int size = fespace->GetTrueVSize();
   cout << "Number of finite element unknowns: "
        << fespace->GetTrueVSize() << endl;
   
   // 5. Determine the list of true (i.e. conforming) essential boundary
   //    dofs.
   //    In this example, the boundary conditions are defined by marking all
   //    the boundary attributes from the mesh as essential (Dirichlet) and
   //    converting them to a list of true dofs.
   Array<int> ess_tdof_list;   // Essential true degrees of freedom
                               // "true" takes into account shared vertices.
   if (mesh->bdr_attributes.Size())
   {
      assert(mesh->bdr_attributes.Max() > 1 && "Can't find a ground boundary!");
      Array<int> ess_bdr(mesh->bdr_attributes.Max());
      ess_bdr = 0;
      ess_bdr[1] = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }
   

   // 7. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   ParGridFunction gf_x(fespace);
   ParGridFunction gf_b(fespace);
   gf_x = 0.0;
   

   //Fill in the MatrixElementPiecewiseCoefficients
   vector<int> bathRegions;
   objectGet(obj,"bath_regions",bathRegions);
   vector<double> sigma_b;
   objectGet(obj,"sigma_b",sigma_b);;
   vector<int> heartRegions;
   objectGet(obj,"heart_regions", heartRegions);
   vector<double> sigma_i;
   vector<double> sigma_e;
   objectGet(obj,"sigma_i",sigma_i);
   objectGet(obj,"sigma_e",sigma_e);

   assert(bathRegions.size() == sigma_b.size());
   assert(heartRegions.size()*3 == sigma_i.size());
   assert(heartRegions.size()*3 == sigma_e.size());

   string fiberFileName;
   objectGet(obj,"fibers",fiberFileName,"");

   shared_ptr<ParGridFunction> fiber_quat;
   {
      ifstream fiberStream(fiberFileName);
      fiber_quat = make_shared<ParGridFunction>(mesh, fiberStream);
   }
   
   MatrixElementPiecewiseCoefficient sigma_i_coeffs(fiber_quat);
   MatrixElementPiecewiseCoefficient sigma_ie_coeffs(fiber_quat);
   for (int ii=0; ii<heartRegions.size(); ii++)
   {
      int heartCursor=3*ii;
      Vector sigma_i_vec(&sigma_i[heartCursor],3);
      Vector sigma_e_vec(&sigma_e[heartCursor],3);
      Vector sigma_ie_vec = sigma_e_vec;
      sigma_ie_vec += sigma_i_vec;
      
      sigma_i_coeffs.heartConductivities_[heartRegions[ii]] = sigma_i_vec;
      sigma_ie_coeffs.heartConductivities_[heartRegions[ii]] = sigma_ie_vec;
   }
   for (int ii=0; ii<bathRegions.size(); ii++)
   {
      sigma_ie_coeffs.bathConductivities_[bathRegions[ii]] = sigma_b[ii];
      sigma_i_coeffs.bathConductivities_[bathRegions[ii]] = 0;
   }

   // 8. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.

   ParBilinearForm *b = new ParBilinearForm(fespace);
   b->AddDomainIntegrator(new DiffusionIntegrator(sigma_i_coeffs));
   b->Assemble();
   // This creates the linear algebra problem.
   HypreParMatrix heart_mat;
   b->FormSystemMatrix(ess_tdof_list, heart_mat);
   heart_mat.Finalize();

   string VmFilename;
   objectGet(obj, "Vm", VmFilename, "");
   shared_ptr<ParGridFunction> gf_Vm;
   {
      ifstream VmStream(VmFilename);
      gf_Vm = make_shared<ParGridFunction>(mesh, VmStream);
   }


   heart_mat.Mult(*gf_Vm, gf_b);
   //cout << heart_mat.Empty() << " "
   cout << heart_mat.Finalized() << " "
        << heart_mat.MaxNorm() << " "
        << heart_mat.Height() << " "
        << heart_mat.Width() << " "
        << heart_mat.ActualWidth() << endl;
   cout << gf_Vm->Max() << " " << gf_Vm->Min() << endl;
   cout << gf_b.Max() << " " << gf_b.Min() << endl;

   ParBilinearForm *a = new ParBilinearForm(fespace);   // defines a.
   // this is the Laplacian: grad u . grad v with linear coefficient.
   // we defined "one" ourselves in step 6.
   a->AddDomainIntegrator(new DiffusionIntegrator(sigma_ie_coeffs));
   a->Assemble();   // This creates the loops.
   HypreParMatrix torso_mat;
   Vector phi_e;
   Vector phi_b;
   a->FormLinearSystem(ess_tdof_list,gf_x,gf_b,torso_mat,phi_e,phi_b);
   
// NOTE THE ifdef
#ifndef MFEM_USE_SUITESPARSE
   // 10. Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //     solve the system A X = B with PCG.
   //   GSSmoother M(torso_mat);
   //   PCG(torso_mat, M, phi_b, phi_e, 1, 200, 1e-12, 0.0);
   HypreSolver *amg = new HypreBoomerAMG(torso_mat);
   HyprePCG *pcg = new HyprePCG(torso_mat);
   pcg->SetTol(1e-12);
   pcg->SetMaxIter(200);
   pcg->SetPrintLevel(2);
   pcg->SetPreconditioner(*amg);
   pcg->Mult(phi_b, phi_e);
#else
   // 10. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
   UMFPackSolver umf_solver;
   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   umf_solver.SetOperator(torso_mat);
   umf_solver.Mult(phi_b, phi_e);
   // See parallel version for HypreSolver - which is an LLNL package.
#endif

   cout << phi_e.Max() << " "
        << phi_e.Min() << " "
        << phi_e.Norml2() << " "
        << phi_e.Size() << endl;
   
   // 11. Recover the solution as a finite element grid function.
   a->RecoverFEMSolution(phi_e, phi_b, gf_x);

   cout << gf_x.Max() << " "
        << gf_x.Min() << " "
        << gf_x.Norml2() << " "
        << gf_x.Size() << endl;

   
   // 12. Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol.gf".
   //ofstream mesh_ofs("refined.mesh");
   //mesh_ofs.precision(8);
   //mesh->Print(mesh_ofs);
   //ofstream sol_ofs("sol.gf");
   //sol_ofs.precision(8);
   //gf_x.SaveAsOne(sol_ofs);

   {
      ostringstream mesh_name, sol_name;
      mesh_name << "mesh." << setfill('0') << setw(6) << myid;
      sol_name << "sol." << setfill('0') << setw(6) << myid;

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(8);
      gf_x.Save(sol_ofs);
   }


   // 14. Free the used memory.
   delete amg;
   delete pcg;
   delete a;
   delete b;
   delete fespace;
   if (order > 0) { delete fec; }
   delete mesh;

   MPI_Finalize();

   return 0;
}
