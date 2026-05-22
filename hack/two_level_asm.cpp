#include "two_level_asm.hpp"

#include <iostream>

using namespace mfem;

namespace
{

struct CoarseCtx
{
   PODCoarseSpace *pod;
   int n_loc;
};

extern "C" PetscErrorCode CoarseApply(PC pc, Vec r, Vec z)
{
   CoarseCtx *ctx = nullptr;
   PetscCall(PCShellGetContext(pc, reinterpret_cast<void **>(&ctx)));

   PetscInt n_loc = 0;
   PetscCall(VecGetLocalSize(r, &n_loc));
   if (ctx && ctx->n_loc != static_cast<int>(n_loc))
   {
      SETERRQ(PetscObjectComm(reinterpret_cast<PetscObject>(pc)),
              PETSC_ERR_ARG_SIZ,
              "POD coarse correction Vec local size does not match MFEM true dofs.");
   }

   const PetscScalar *r_arr = nullptr;
   PetscScalar *z_arr = nullptr;
   PetscCall(VecGetArrayRead(r, &r_arr));
   PetscCall(VecGetArray(z, &z_arr));

   Vector r_mfem(const_cast<PetscScalar *>(r_arr), static_cast<int>(n_loc));
   Vector z_mfem(z_arr, static_cast<int>(n_loc));
   ctx->pod->ApplyCoarse(r_mfem, z_mfem);

   PetscCall(VecRestoreArrayRead(r, &r_arr));
   PetscCall(VecRestoreArray(z, &z_arr));
   return PETSC_SUCCESS;
}

extern "C" PetscErrorCode CoarseDestroy(PC pc)
{
   CoarseCtx *ctx = nullptr;
   PetscCall(PCShellGetContext(pc, reinterpret_cast<void **>(&ctx)));
   delete ctx;
   PetscCall(PCShellSetContext(pc, nullptr));
   return PETSC_SUCCESS;
}

} // namespace

void AttachTwoLevelASM(PetscPCGSolver *ksp_wrapper,
                       PODCoarseSpace *pod,
                       ParFiniteElementSpace *fes,
                       MPI_Comm comm,
                       int asm_overlap,
                       int icc_levels)
{
   int my_rank = 0;
   MPI_Comm_rank(comm, &my_rank);

   KSP ksp = static_cast<KSP>(*ksp_wrapper);
   PC pc_top = nullptr;
   PetscErrorCode ierr = KSPGetPC(ksp, &pc_top);
   MFEM_VERIFY(ierr == 0, "KSPGetPC failed while attaching POD ASM.");

   ierr = KSPSetType(ksp, KSPCG);
   MFEM_VERIFY(ierr == 0, "KSPSetType(KSPCG) failed.");

   ierr = PCSetType(pc_top, PCCOMPOSITE);
   MFEM_VERIFY(ierr == 0, "PCSetType(PCCOMPOSITE) failed.");
   ierr = PCCompositeSetType(pc_top, PC_COMPOSITE_ADDITIVE);
   MFEM_VERIFY(ierr == 0, "PCCompositeSetType(additive) failed.");

   ierr = PCCompositeAddPCType(pc_top, PCSHELL);
   MFEM_VERIFY(ierr == 0, "PCCompositeAddPCType(PCSHELL) failed.");
   PC pc_coarse = nullptr;
   ierr = PCCompositeGetPC(pc_top, 0, &pc_coarse);
   MFEM_VERIFY(ierr == 0, "PCCompositeGetPC(coarse) failed.");

   CoarseCtx *ctx = new CoarseCtx;
   ctx->pod = pod;
   ctx->n_loc = fes->GetTrueVSize();
   ierr = PCShellSetContext(pc_coarse, ctx);
   MFEM_VERIFY(ierr == 0, "PCShellSetContext failed.");
   ierr = PCShellSetApply(pc_coarse, CoarseApply);
   MFEM_VERIFY(ierr == 0, "PCShellSetApply failed.");
   ierr = PCShellSetDestroy(pc_coarse, CoarseDestroy);
   MFEM_VERIFY(ierr == 0, "PCShellSetDestroy failed.");
   ierr = PCShellSetName(pc_coarse, "POD_CoarseCorrection_ZZt");
   MFEM_VERIFY(ierr == 0, "PCShellSetName failed.");

   ierr = PCCompositeAddPCType(pc_top, PCASM);
   MFEM_VERIFY(ierr == 0, "PCCompositeAddPCType(PCASM) failed.");
   PC pc_asm = nullptr;
   ierr = PCCompositeGetPC(pc_top, 1, &pc_asm);
   MFEM_VERIFY(ierr == 0, "PCCompositeGetPC(ASM) failed.");
   ierr = PCASMSetType(pc_asm, PC_ASM_BASIC);
   MFEM_VERIFY(ierr == 0, "PCASMSetType(PC_ASM_BASIC) failed.");
   ierr = PCASMSetOverlap(pc_asm, asm_overlap);
   MFEM_VERIFY(ierr == 0, "PCASMSetOverlap failed.");

   ierr = PCSetUp(pc_top);
   MFEM_VERIFY(ierr == 0, "PCSetUp failed for POD composite PC.");
   ierr = PCSetUp(pc_asm);
   MFEM_VERIFY(ierr == 0, "PCSetUp failed for ASM child PC.");
   ierr = KSPSetUp(ksp);
   MFEM_VERIFY(ierr == 0, "KSPSetUp failed after adding POD two-level ASM.");

   KSP *sub_ksps = nullptr;
   PetscInt n_local = 0;
   PetscInt first = 0;
   ierr = PCASMGetSubKSP(pc_asm, &n_local, &first, &sub_ksps);
   MFEM_VERIFY(ierr == 0, "PCASMGetSubKSP failed.");

   for (PetscInt i = 0; i < n_local; ++i)
   {
      ierr = KSPSetType(sub_ksps[i], KSPPREONLY);
      MFEM_VERIFY(ierr == 0, "KSPSetType(KSPPREONLY) failed for ASM sub-KSP.");
      PC sub_pc = nullptr;
      ierr = KSPGetPC(sub_ksps[i], &sub_pc);
      MFEM_VERIFY(ierr == 0, "KSPGetPC failed for ASM sub-KSP.");
      ierr = PCSetType(sub_pc, PCICC);
      MFEM_VERIFY(ierr == 0, "PCSetType(PCICC) failed for ASM sub-PC.");
      ierr = PCFactorSetLevels(sub_pc, icc_levels);
      MFEM_VERIFY(ierr == 0, "PCFactorSetLevels failed for ASM sub-PC.");
      ierr = KSPSetUp(sub_ksps[i]);
      MFEM_VERIFY(ierr == 0, "KSPSetUp failed for ASM sub-KSP.");
   }

   ierr = PCSetUpOnBlocks(pc_asm);
   MFEM_VERIFY(ierr == 0, "PCSetUpOnBlocks failed for ASM child PC.");

   if (my_rank == 0)
   {
      std::cout << "[POD] two-level ASM attached: PCCOMPOSITE(SHELL+ASM), "
                << "overlap = " << asm_overlap
                << ", ICC levels = " << icc_levels << std::endl;
   }
}
