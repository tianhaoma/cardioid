int solvePseudoBidomainForUe(
    ParMesh* pmesh,
    ParFiniteElementSpace* pfespace,
    ParGridFunction& V_m,//gf_Vm
    ParGridFunction& u_e,//gf_ue
    const std::vector<double>& sigma_i_values,//sigma_i_values
    const std::vector<double>& sigma_e_values,//sigma_e_values
    std::shared_ptr<ParGridFunction>& fiber_quat,//fiber_quat
    std::shared_ptr<ParGridFunction>& sheet_quat,//sheet_quat
    std::shared_ptr<ParGridFunction>& transverse_quat,//transverse_quat
    Array<int>& ess_tdof_list,//ess_tdof_list
    const std::vector<int>& heartRegions,
    bool enforce_zero_mean = true,
    int print_level = 2,
    double* t_ksp2_total = nullptr)

                  solvePseudoBidomainForUe(
                  pmesh, pfespace, gf_Vm, gf_ue, 
                  sigma_i_values, sigma_e_values, fiber_quat, sheet_quat, transverse_quat,
                  ess_tdof_list, heartRegions, true,
                  global_print_level, &t_ksp2_total
              );
{



    
    // 计时开始
    double start_time = MPI_Wtime();
    
    if (my_rank == 0) {
        std::cout << "==============================================" << std::endl;
        std::cout << "开始求解伪双域模型方程..." << std::endl;
        std::cout << "进程数: " << num_procs << std::endl;
        std::cout << "打印级别: " << print_level << std::endl;
        std::cout << "强制均值为零: " << (enforce_zero_mean ? "是" : "否") << std::endl;
        std::cout << "==============================================" << std::endl;
    }
    
    // 1. 参数验证
    if (my_rank == 0) std::cout << "[1/8] 验证输入参数..." << std::endl;
    
    
    // 2. 初始化电导率张量
    if (my_rank == 0) std::cout << "[2/8] 初始化电导率张量..." << std::endl;
    
    double tensor_start_time = MPI_Wtime();
    
    MatrixElementPiecewiseCoefficient sigma_i(fiber_quat, sheet_quat, transverse_quat);
    MatrixElementPiecewiseCoefficient sigma_sum(fiber_quat, sheet_quat, transverse_quat);
    
    for (int ii = 0; ii < heartRegions.size(); ii++) {
        int heartCursor = 3 * ii;
        
        Vector sigma_i_vec(3);
        Vector sigma_e_vec(3);
        Vector sigma_sum_vec(3);
        
        for (int jj = 0; jj < 3; jj++) {
            sigma_i_vec[jj] = sigma_i_values[heartCursor + jj];
            sigma_e_vec[jj] = sigma_e_values[heartCursor + jj];
            sigma_sum_vec[jj] = (sigma_i_vec[jj] + sigma_e_vec[jj]);
        }
        
        sigma_i.heartConductivities_[heartRegions[ii]] = sigma_i_vec;
        sigma_sum.heartConductivities_[heartRegions[ii]] = sigma_sum_vec;
    }
    

    
    double tensor_end_time = MPI_Wtime();
    
    // 3. 建立线性系统
    if (my_rank == 0) std::cout << "[3/8] 构建线性系统..." << std::endl;
    
    double system_start_time = MPI_Wtime();
    
    // 设置左侧矩阵: -∇·((σ_i + σ_e)∇u_e)//a_pblf_recoverue
    ParBilinearForm *a_pblf_recoverue = new ParBilinearForm(pfespace);
    if(use_petsc)
    {
    a_pblf_recoverue->SetOperatorType(Operator::PETSC_MATAIJ);
    }
    a_pblf_recoverue->AddDomainIntegrator(new DiffusionIntegrator(sigma_sum));
    a_pblf_recoverue->Assemble();

    HypreBoomerAMG* precond_recoverue_hypre = nullptr;
    HyprePCG* pcg_recoverue_hypre = nullptr;
    HypreParMatrix A_recoverue_hypre;

    PetscPCGSolver* pcg_recoverue_petsc = nullptr;
    PetscParMatrix A_recoverue_petsc;

    Vector B_recoverue, X_recoverue;

    //HypreParMatrix A_recoverue_hypre;
    a_pblf_recoverue->FormSystemMatrix(ess_tdof_list, A_recoverue_hypre);
    
    // 设置右侧向量: ∇·(σ_i∇V_m)
    
    ParBilinearForm temp_form(pfespace);
    if(use_petsc)
    {
    temp_form.SetOperatorType(Operator::PETSC_MATAIJ);
    }
    temp_form.AddDomainIntegrator(new DiffusionIntegrator(sigma_i));
    temp_form.Assemble();
    
    Vector vm_true(pfespace->GetTrueVSize());
    gf_Vm.GetTrueDofs(vm_true);//time loop after assigning the gf_Vm
    
    HypreParMatrix A_temp_hypre;
    PetscParMatrix A_temp_petsc;
    if(!use_petsc)
    {
    temp_form.FormSystemMatrix(ess_tdof_list, A_temp_hypre);
    }
    else
    {
    temp_form.FormSystemMatrix(ess_tdof_list, A_temp_petsc);
    }
    
    Vector rhs_recoverue(pfespace->GetTrueVSize());
    rhs_recoverue = 0.0;
    if(!use_petsc)
    {
    A_temp_hypre.Mult(-1.0, vm_true, 0.0, rhs_recoverue);//
    }
    else
    {
    A_temp_petsc.Mult(-1.0, vm_true, 0.0, rhs_recoverue);//
    }
    
    
    double system_end_time = MPI_Wtime();
    
    // 4. 检查相容性条件
    //if (my_rank == 0) std::cout << "[4/8] 检查相容性条件..." << std::endl;
    
    double compat_start_time = MPI_Wtime();
    double compat_end_time = MPI_Wtime();
    
    // 5. 输出右侧向量详细信息
    if (my_rank == 0) std::cout << "[5/8] 分析右侧向量..." << std::endl;
    
    double analyze_start_time = MPI_Wtime();
    double analyze_end_time = MPI_Wtime();
    
    // 6. 设置求解器和预处理器
    if (my_rank == 0) std::cout << "[6/8] 配置求解器..." << std::endl;
    
    double solver_setup_time = MPI_Wtime();
    
    // 设置求解器
    //HyprePCG pcg(A);
    if(!use_petsc)
    {
        precond_recoverue_hypre = new HypreBoomerAMG;
        pcg_recoverue_hypre = new HyprePCG(MPI_COMM_WORLD);    
        a_pblf_recoverue->FormSystemMatrix(ess_tdof_list, A_recoverue_hypre);
        precond_recoverue_hypre->SetPrintLevel(0);
        pcg_recoverue_hypre->SetPreconditioner(*precond_recoverue_hypre);
        pcg_recoverue_hypre->SetOperator(A_recoverue_hypre);
        pcg_recoverue_hypre->SetTol(1e-6);
        pcg_recoverue_hypre->SetMaxIter(1000);
        pcg_recoverue_hypre->SetPrintLevel(1);
    }
    else
    {
        pcg_recoverue_petsc = new PetscPCGSolver(MPI_COMM_WORLD, "recoverue_", true);
    }

    
    double solver_setup_end_time = MPI_Wtime();
    
    // 7. 求解线性系统
    if (my_rank == 0) std::cout << "[7/8] 求解线性系统..." << std::endl;
    double solve_start_time = MPI_Wtime();
    
    X_recoverue = 0.0;  // 初始猜测为零
    
    // 执行求解
    double t_ksp2_start = MPI_Wtime();
    if(!use_petsc)
    {
        pcg_recoverue_hypre.Mult(rhs_recoverue, X_recoverue);
    }
    else
    {
        pcg_recoverue_petsc.Mult(rhs_recoverue, X_recoverue);
    }

    double t_ksp2_end = MPI_Wtime();
    if (t_ksp2_total != nullptr) {
    *t_ksp2_total += (t_ksp2_end - t_ksp2_start);
    }
    
    double solve_end_time = MPI_Wtime();
    
    // 8. 后处理解向量
    if (my_rank == 0) std::cout << "[8/8] 后处理解向量..." << std::endl;
    
    double postproc_start_time = MPI_Wtime();
    
    // 将解向量设置到网格函数
    gf_ue.SetFromTrueDofs(X_recoverue);
    
    // 如果需要，强制解的均值为零
    if (enforce_zero_mean && ess_tdof_list.Size() == 0) {
        double local_sum = 0.0;
        for (int i = 0; i < X_recoverue.Size(); i++) {
            local_sum += X_recoverue(i);
        }
        
        double global_sum = 0.0;
        MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        double global_count = 0.0;
        local_sum = X_recoverue.Size();
        MPI_Allreduce(&local_sum, &global_count, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        double mean_value = global_sum / global_count;
        
        if (fabs(mean_value) > 1e-6) {
            if (my_rank == 0 && print_level > 0) {
                std::cout << "调整解向量以确保均值为零，当前均值: " << mean_value << std::endl;
            }
            
            // 从解中减去平均值
            for (int i = 0; i < X_recoverue.Size(); i++) {
                X_recoverue(i) -= mean_value;
            }
            
            // 更新网格函数
            gf_ue.SetFromTrueDofs(X_recoverue);
            
            if (my_rank == 0 && print_level > 0) {
                std::cout << "已从解中减去均值: " << mean_value << std::endl;
            }
        } else {
            if (my_rank == 0 && print_level > 0) {
                std::cout << "解向量均值已经接近零: " << mean_value << "，无需调整。" << std::endl;
            }
        }
    }
    
    
    double postproc_end_time = MPI_Wtime();
    if (my_rank == 0 && print_level > 0) {
        std::cout << "后处理完成，用时: " << (postproc_end_time - postproc_start_time) << " 秒" << std::endl;
    }
    
    // 9. 清理资源
    delete a;
    
    // 10. 总结
    double end_time = MPI_Wtime();
    double total_time = end_time - start_time;
    
    if (my_rank == 0) {
        std::cout << "==============================================" << std::endl;
        std::cout << "伪双域求解总结:" << std::endl;
        std::cout << "  总用时: " << total_time << " 秒" << std::endl;
        //std::cout << "  矩阵大小: " << A.Height() << " x " << A.Width() << std::endl;
        //std::cout << "  迭代次数: " << num_iterations << std::endl;
        //std::cout << "  最终残差: " << final_res_norm << std::endl;
        std::cout << "  解范围: [" << X_recoverue.Min() << ", " << X_recoverue.Max() << "]" << std::endl;
        std::cout << "==============================================" << std::endl;
    }
    
}