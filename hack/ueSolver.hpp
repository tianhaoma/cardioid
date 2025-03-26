#pragma once

#include "mfem.hpp"
#include <memory>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>

using namespace mfem;

void setHeartTorsoBoundarySafe(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    Array<int>& ess_tdof_list,
    ParFiniteElementSpace* pfespace_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    double tolerance = 1e-6);


Array<int> getHeartTorsoInterfaceDofs(
    ParFiniteElementSpace* pfespace_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    double tolerance = 1e-6);

void transferHeartToTorsoBoundary(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    Array<int>& ess_tdof_list_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    ParFiniteElementSpace* pfespace_torso,
    double tolerance = 1e-6,
    bool debug = true);

void setHeartToTorsoBoundaryValues(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    Array<int>& ess_tdof_list_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    ParFiniteElementSpace* pfespace_torso,
    double tolerance = 1e-6,
    bool debug = true);

void setHeartToTorsoBoundaryCorrectly(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    Array<int>& ess_tdof_list_torso,
    ParFiniteElementSpace* pfespace_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    double tolerance = 1e-6);

void transferHeartToTorsoBoundaryWithDebug(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    Array<int>& ess_tdof_list_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    ParFiniteElementSpace* pfespace_torso,
    double tolerance = 1e-6,
    int heart_bdry_attr = 1,
    const std::string& outputDir = ".",
    int timestep = 0);

void setHeartToTorsoBoundaryValuesWithDebug(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    Array<int>& ess_tdof_list_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    ParFiniteElementSpace* pfespace_torso,
    double tolerance = 1e-6,
    bool debug = true); 

