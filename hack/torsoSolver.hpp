#pragma once

#include "mfem.hpp"
#include <memory>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>

using namespace mfem;
using namespace std;

int solveTorsoModel(
    ParMesh* pmesh_torso,
    ParFiniteElementSpace* pfespace_torso,
    ParGridFunction& gf_ue_torso,
    double sigma_T,
    Array<int>& ess_tdof_list_torso,
    int print_level = 2);

int torsoSolver(
        Mesh* heart_mesh,
        Mesh* torso_mesh,
        GridFunction& ue,
        GridFunction& ue_torso,
        FiniteElementSpace* heart_fespace,
        FiniteElementSpace* torso_fespace,
        double sigma_T);