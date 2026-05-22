#ifndef CARDIOID_HACK_TWO_LEVEL_ASM_HPP
#define CARDIOID_HACK_TWO_LEVEL_ASM_HPP

#include "mfem.hpp"
#include "pod_coarse_space.hpp"

#include <mpi.h>
#include <petscksp.h>

void AttachTwoLevelASM(mfem::PetscPCGSolver *ksp_wrapper,
                       PODCoarseSpace *pod,
                       mfem::ParFiniteElementSpace *fes,
                       MPI_Comm comm,
                       int asm_overlap,
                       int icc_levels);

#endif
