#include "torsoSolver.hpp"

// 用于在控制台输出颜色
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"
#define RESET "\033[0m"


int torsoSolver(
    Mesh* heart_mesh,
    Mesh* torso_mesh,
    GridFunction& ue,
    GridFunction& ue_torso,
    FiniteElementSpace* heart_fespace,
    FiniteElementSpace* torso_fespace,
    double sigma_T)
{

    int my_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    int heart_interface_attr = 1;  // 心脏交界面的边界属性编号
    int torso_interface_attr = 1;  // 躯干交界面的边界属性编号
    double match_tolerance = 1e-6; // 坐标匹配的容差（非常量）

    // 注意：这里修正了指针传递
    GridFunction interface_marker(torso_fespace);  // 标记交界面的函数
    interface_marker = 0.0;

    // 6. 创建交界面的标记函数
    GridFunction heart_marker(heart_fespace);
    GridFunction torso_marker(torso_fespace);
    heart_marker = 0.0;
    torso_marker = 0.0;

    // 记录交界面的顶点DOF和坐标
    vector<int> heart_interface_dofs;
    vector<int> torso_interface_dofs;
    
    // 7. 收集心脏交界面顶点
    if (my_id == 0) {
        cout << YELLOW "收集心脏交界面顶点..." RESET << endl;
    }

    for (int i = 0; i < heart_mesh->GetNBE(); i++) {
        int attr = heart_mesh->GetBdrAttribute(i);
        if (attr == heart_interface_attr) {
            Array<int> dofs;
            // 使用箭头操作符而不是点操作符
            heart_fespace->GetBdrElementDofs(i, dofs);
            
            for (int j = 0; j < dofs.Size(); j++) {
                int dof = dofs[j];
                heart_interface_dofs.push_back(dof);
                heart_marker(dof) = 1.0;
            }
        }
    }

    // 去重（简单方法）
    sort(heart_interface_dofs.begin(), heart_interface_dofs.end());
    heart_interface_dofs.erase(
        unique(heart_interface_dofs.begin(), heart_interface_dofs.end()),
        heart_interface_dofs.end()
    );
    
    if (my_id == 0) {
        cout << GREEN "心脏交界面顶点数量: " << heart_interface_dofs.size() << RESET << endl;
    }

    // 8. 收集躯干交界面顶点
    if (my_id == 0) {
        cout << YELLOW "收集躯干交界面顶点..." RESET << endl;
    }
    
    for (int i = 0; i < torso_mesh->GetNBE(); i++) {
        int attr = torso_mesh->GetBdrAttribute(i);
        if (attr == torso_interface_attr) {
            Array<int> dofs;
            // 使用箭头操作符
            torso_fespace->GetBdrElementDofs(i, dofs);
            
            for (int j = 0; j < dofs.Size(); j++) {
                int dof = dofs[j];
                torso_interface_dofs.push_back(dof);
                torso_marker(dof) = 1.0;
            }
        }
    }
    
    // 去重（简单方法）
    sort(torso_interface_dofs.begin(), torso_interface_dofs.end());
    torso_interface_dofs.erase(
        unique(torso_interface_dofs.begin(), torso_interface_dofs.end()),
        torso_interface_dofs.end()
    );

    if (my_id == 0) {
        cout << GREEN "躯干交界面顶点数量: " << torso_interface_dofs.size() << RESET << endl;
    }
    
    // 10. 执行值传输（点匹配算法）
    if (my_id == 0) {
        cout << BLUE "\n开始执行值传输..." RESET << endl;
    }
    
    int match_count = 0;
    double total_dist = 0.0;
    
    // 记录已匹配的躯干节点
    // 使用箭头操作符
    vector<bool> matched_torso_nodes(torso_fespace->GetNDofs(), false);
    
    for (size_t i = 0; i < torso_interface_dofs.size(); i++) {
        int torso_dof = torso_interface_dofs[i];
        
        // 获取躯干点坐标
        double* torso_coord = torso_mesh->GetVertex(torso_dof);
        
        // 寻找最近的心脏交界面点
        double min_dist = 1e10;
        int closest_heart_dof = -1;
        
        for (size_t j = 0; j < heart_interface_dofs.size(); j++) {
            int heart_dof = heart_interface_dofs[j];
            
            // 获取心脏点坐标
            double* heart_coord = heart_mesh->GetVertex(heart_dof);
            
            // 计算距离
            double dx = torso_coord[0] - heart_coord[0];
            double dy = torso_coord[1] - heart_coord[1];
            double dz = 0.0;
            if (heart_mesh->Dimension() == 3 && torso_mesh->Dimension() == 3) {
                dz = torso_coord[2] - heart_coord[2];
            }
            
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            
            if (dist < min_dist) {
                min_dist = dist;
                closest_heart_dof = heart_dof;
            }
        }
        
        // 只在找到足够近的匹配点的情况下赋值
        if (closest_heart_dof >= 0 && min_dist < match_tolerance) {
            ue_torso(torso_dof) = ue(closest_heart_dof);
            interface_marker(torso_dof) = 1.0;  // 标记这是交界面点
            matched_torso_nodes[torso_dof] = true;
            match_count++;
            total_dist += min_dist;
            
            // 输出少量匹配信息
            if (my_id == 0 && i < 5) {
                cout << YELLOW "躯干点 " << torso_dof 
                     << " 匹配到心脏点 " << closest_heart_dof 
                     << "，距离: " << min_dist 
                     << "，电位值: " << ue(closest_heart_dof) << RESET << endl;
            }
        }
        else {
            // 对于未匹配点，输出警告但不赋值（保持原值为0）
            if (my_id == 0 && i < 5) {
                cout << RED "警告：躯干点 " << torso_dof 
                     << " 未找到精确匹配点，最近点 " << closest_heart_dof 
                     << "，距离: " << min_dist 
                     << " > " << match_tolerance << RESET << endl;
            }
        }

    }
    
    // 输出匹配统计
    double avg_dist = match_count > 0 ? total_dist / match_count : 0.0;
    
    if (my_id == 0) {
        cout << GREEN "精确匹配点数量: " << match_count << "/" << torso_interface_dofs.size() 
             << " (" << (100.0 * match_count / torso_interface_dofs.size()) << "%)" RESET << endl;
        cout << GREEN "平均匹配距离: " << avg_dist << RESET << endl;
    }
    
    // 11. 保存边界条件
    {
        ofstream bc_out("torso_boundary_condition.vtk");
        bc_out.precision(8);
        ue_torso.Save(bc_out);
        
        if (my_id == 0) {
            cout << GREEN "已保存边界条件到 torso_boundary_condition.vtk" RESET << endl;
        }
    }
    
    // 12. 求解Laplace方程
    if (my_id == 0) {
        cout << BLUE "\n开始求解Laplace方程..." RESET << endl;
    }
    

    ConstantCoefficient sigma_T_coeff(-sigma_T);
    // 创建一个bilinear form: (∇u, ∇v)
    BilinearForm a(torso_fespace);
    a.AddDomainIntegrator(new DiffusionIntegrator(sigma_T_coeff));
    a.Assemble();
    a.Finalize();
    
    // 创建一个linear form: 右侧为0（Laplace方程）
    LinearForm b(torso_fespace);
    b.Assemble();
    
    // 应用本质边界条件 - 只对匹配的交界面点设置

    Array<int> ess_tdof_list;
    for (size_t i = 0; i < torso_interface_dofs.size(); i++) {
        int dof = torso_interface_dofs[i];
        if (matched_torso_nodes[dof]) {
            ess_tdof_list.Append(dof);
        }
    }

    #if 0
    // 修正后的应用边界条件部分：
    Array<int> ess_tdof_list;
    for (size_t i = 0; i < torso_interface_dofs.size(); i++) {
        int dof = torso_interface_dofs[i];
        // 确保所有交界面节点都被处理
        ess_tdof_list.Append(dof);
        // 如果未匹配，可能需要处理错误
        if (!matched_torso_nodes[dof]) {
            if (my_id == 0) {
                cerr << RED "错误：躯干交界面节点 " << dof << " 未找到匹配！" RESET << endl;
            }
            // 可根据需求终止程序或设置默认值
            // return 1;
        }
    }
    #endif
    
    if (my_id == 0) {
        cout << YELLOW "本质边界点数量: " << ess_tdof_list.Size() << RESET << endl;
    }
    
    // 设置线性系统
    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, ue_torso, b, A, X, B, 1);

    cout << "矩阵A实际大小: " << ((SparseMatrix&)(*A)).Height() << " x " 
     << ((SparseMatrix&)(*A)).Width() << endl;
    
    if (my_id == 0) {
        cout << YELLOW "线性系统大小: " << A->Height() << RESET << endl;
        cout << YELLOW "线性系统右侧范数: " << B.Norml2() << RESET << endl;
        cout << YELLOW "开始求解线性系统..." RESET << endl;
    }
    
    // 创建预处理器
    GSSmoother M((SparseMatrix&)(*A));
    
    // 使用PCG求解系统 - 注意参数顺序！
    int print_level = 3;  // 详细输出
    int max_iter = 1000;
    double rtol = 1e-12;
    double atol = 1e-12;
    
    // PCG参数顺序: (A, Preconditioner, RHS, solution, print_level, max_iter, rtol, atol)
    PCG(*A, M, B, X, print_level, max_iter, rtol, atol);
    
    if (my_id == 0) {
        cout << GREEN "线性系统求解完成!" RESET << endl;
    }
    
    // 恢复解
    a.RecoverFEMSolution(X, b, ue_torso);
    
    // 计算解的统计信息
    double min_val = ue_torso.Min();
    double max_val = ue_torso.Max();
    double l2_norm = ue_torso.Norml2();
    
    if (my_id == 0) {
        cout << YELLOW "解的统计信息:" RESET << endl;
        cout << YELLOW "  最小值: " << min_val << RESET << endl;
        cout << YELLOW "  最大值: " << max_val << RESET << endl;
        cout << YELLOW "  L2范数: " << l2_norm << RESET << endl;
    }
    
    // 13. 输出结果
    {
        ofstream torso_out("torso_potential.vtk");
        torso_out.precision(8);
        ue_torso.Save(torso_out);
        
        if (my_id == 0) {
            cout << GREEN "已保存求解结果到 torso_potential.vtk" RESET << endl;
        }
    }
    
    // 14. 清理内存
    //delete heart_mesh;
    //delete torso_mesh;
    
    if (my_id == 0) {
        cout << BLUE "\n测试已完成！" RESET << endl;
    }


    return 0;
}









int solveTorsoModel(
    ParMesh* pmesh_torso,
    ParFiniteElementSpace* pfespace_torso,
    ParGridFunction& gf_ue_torso,
    double sigma_T,
    Array<int>& ess_tdof_list_torso,
    int print_level)
 {
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0 && print_level > 0) {
        std::cout << "求解Torso模型方程..." << std::endl;
    }
    
    // 检查边界条件是否合理
    bool boundary_valid = true;
    double min_boundary_value = 1e10;
    double max_boundary_value = -1e10;
    
    for (int i = 0; i < ess_tdof_list_torso.Size(); i++) {
        int dof = ess_tdof_list_torso[i];
        if (dof >= 0 && dof < gf_ue_torso.Size()) {
            double value = gf_ue_torso(dof);
            min_boundary_value = std::min(min_boundary_value, value);
            max_boundary_value = std::max(max_boundary_value, value);
            
            if (std::isnan(value) || std::isinf(value)) {
                boundary_valid = false;
                break;
            }
        }
    }
    
    if (my_rank == 0 && print_level > 0) {
        std::cout << "  边界条件值范围: [" << min_boundary_value << ", " << max_boundary_value << "]" << std::endl;
    }
    
    if (!boundary_valid) {
        if (my_rank == 0) {
            std::cout << "错误：边界条件无效，中止求解。" << std::endl;
        }
        return -1;
    }
    
    // 设置常数电导率系数
    ConstantCoefficient sigma_T_coeff(-sigma_T);  // 注意符号：Diffusion算子是-div(sigma*grad)
    
    // 设置双线性形式
    ParBilinearForm *a_torso = new ParBilinearForm(pfespace_torso);
    a_torso->AddDomainIntegrator(new DiffusionIntegrator(sigma_T_coeff));
    a_torso->Assemble();
    
    // 创建零线性形式作为右侧向量（无源项）
    ParLinearForm *f_torso = new ParLinearForm(pfespace_torso);
    f_torso->Assemble();
    
    // 正确设置线性系统
    HypreParMatrix A_torso;
    Vector B_torso, X_torso;
    
    // 使用FormLinearSystem正确地构建考虑边界条件的线性系统
    a_torso->FormLinearSystem(ess_tdof_list_torso, gf_ue_torso, *f_torso, 
                           A_torso, X_torso, B_torso);
    
    if (my_rank == 0 && print_level > 1) {
        std::cout << "  线性系统已准备完成" << std::endl;
        std::cout << "  矩阵大小: " << A_torso.Height() << " x " << A_torso.Width() << std::endl;
        std::cout << "  右侧向量范数: " << B_torso.Norml2() << std::endl;
    }
    
    // 设置求解器
    HyprePCG pcg_torso(A_torso);
    pcg_torso.SetTol(1e-12);
    pcg_torso.SetMaxIter(1000);
    pcg_torso.SetPrintLevel(print_level > 1 ? 2 : 0);
    
    // 设置预处理器
    HypreBoomerAMG amg_torso(A_torso);
    amg_torso.SetPrintLevel(0);
    pcg_torso.SetPreconditioner(amg_torso);
    
    // 求解系统
    pcg_torso.Mult(B_torso, X_torso);
    
    // 获取迭代次数
    int num_iterations = 0;
    pcg_torso.GetNumIterations(num_iterations);
    
    // 检查结果是否合理
    double min_value = X_torso.Min();
    double max_value = X_torso.Max();
    
    if (my_rank == 0 && print_level > 0) {
        std::cout << "  解的范围: [" << min_value << ", " << max_value << "]" << std::endl;
        std::cout << "  PCG迭代次数: " << num_iterations << std::endl;
    }
    
    // 将解恢复到网格函数（自动处理边界条件）
    a_torso->RecoverFEMSolution(X_torso, *f_torso, gf_ue_torso);
    
    // 清理
    delete a_torso;
    delete f_torso;
    
    return num_iterations;
 }
 