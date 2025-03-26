void setHeartTorsoBoundarySafe(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    Array<int>& ess_tdof_list,
    ParFiniteElementSpace* pfespace_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    double tolerance = 1e-6)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0) {
        std::cout << "设置心脏-躯干边界条件(安全版)..." << std::endl;
    }
    
    // 1. 标记躯干网格所有边界为Dirichlet边界
    ParMesh* pmesh_torso = pfespace_torso->GetParMesh();
    int max_attr = pmesh_torso->bdr_attributes.Max();
    Array<int> bdr_marker(max_attr);
    bdr_marker = 1;  // 所有边界属性都标记为Dirichlet边界
    
    // 2. 获取边界DOF列表
    pfespace_torso->GetEssentialTrueDofs(bdr_marker, ess_tdof_list);
    
    if (my_rank == 0) {
        std::cout << "找到 " << ess_tdof_list.Size() << " 个边界DOF" << std::endl;
    }
    
    // 3. 用简单的方法设置边界值 - 不使用复杂的系数类
    gf_ue_torso = 0.0;  // 初始化为0
    
    // 只处理本地拥有的顶点
    for (int i = 0; i < torso_mesh->GetNBE(); i++) {
        Array<int> vertices;
        torso_mesh->GetBdrElementVertices(i, vertices);
        
        for (int j = 0; j < vertices.Size(); j++) {
            int torso_vertex_id = vertices[j];
            
            // 确保这个顶点是本地的
            if (torso_vertex_id >= gf_ue_torso.Size()) continue;
            
            double* torso_coords = torso_mesh->GetVertex(torso_vertex_id);
            
            // 查找最近的心脏点
            double min_distance = 1e10;
            double heart_value = 0.0;
            bool found_match = false;
            
            // 只搜索本地拥有的心脏点
            for (int k = 0; k < heart_mesh->GetNBE(); k++) {
                Array<int> heart_vertices;
                heart_mesh->GetBdrElementVertices(k, heart_vertices);
                
                for (int l = 0; l < heart_vertices.Size(); l++) {
                    int heart_vertex_id = heart_vertices[l];
                    
                    if (heart_vertex_id >= gf_ue_heart.Size()) continue;
                    
                    double* heart_coords = heart_mesh->GetVertex(heart_vertex_id);
                    
                    // 计算距离
                    double dx = heart_coords[0] - torso_coords[0];
                    double dy = heart_coords[1] - torso_coords[1];
                    double dz = heart_coords[2] - torso_coords[2];
                    double dist = sqrt(dx*dx + dy*dy + dz*dz);
                    
                    if (dist < min_distance) {
                        min_distance = dist;
                        heart_value = gf_ue_heart(heart_vertex_id);
                        found_match = true;
                    }
                }
            }
            
            if (found_match && min_distance <= tolerance) {
                gf_ue_torso(torso_vertex_id) = heart_value;
            }
        }
    }
    
    // 同步所有进程的值
    gf_ue_torso.ParallelAssemble();
    
    // 验证边界值
    Vector ue_true;
    gf_ue_torso.GetTrueDofs(ue_true);
    
    int nonzero_count = 0;
    double min_val = 1e10, max_val = -1e10;
    
    for (int i = 0; i < ess_tdof_list.Size(); i++) {
        int tdof = ess_tdof_list[i];
        if (tdof >= 0 && tdof < ue_true.Size()) {
            double value = ue_true(tdof);
            if (std::abs(value) > 1e-10) {
                nonzero_count++;
                min_val = std::min(min_val, value);
                max_val = std::max(max_val, value);
            }
        }
    }
    
    // 收集统计信息(小心使用MPI_Allreduce)
    int global_nonzero = 0;
    MPI_Allreduce(&nonzero_count, &global_nonzero, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    double global_min = 0.0, global_max = 0.0;
    if (nonzero_count > 0) {
        // 只在确实有非零值时进行reduce，避免未初始化值的问题
        global_min = min_val;
        global_max = max_val;
        MPI_Allreduce(MPI_IN_PLACE, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    
    if (my_rank == 0) {
        std::cout << "边界条件非零值数量: " << global_nonzero 
                 << ", 范围: [" << global_min << ", " << global_max << "]" << std::endl;
    }
}

/**
 * 基于MFEM官方示例设置心脏-躯干边界条件
 */
 void setHeartTorsoBoundaryMFEMWay(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    Array<int>& ess_tdof_list,
    ParFiniteElementSpace* pfespace_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    double tolerance = 1e-6)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0) {
        std::cout << "使用MFEM官方方法设置心脏-躯干边界条件..." << std::endl;
    }
    
    // 1. 获取躯干网格边界属性的最大值
    ParMesh* pmesh_torso = pfespace_torso->GetParMesh();
    int max_attr = pmesh_torso->bdr_attributes.Max();
    
    // 2. 创建边界标记数组，标记所有边界属性为Dirichlet边界
    Array<int> bdr_marker(max_attr);
    bdr_marker = 1;  // 标记所有边界为Dirichlet边界
    
    if (my_rank == 0) {
        std::cout << "边界属性最大值: " << max_attr << std::endl;
    }
    
    // 3. 获取关联的本质边界DOF
    pfespace_torso->GetEssentialTrueDofs(bdr_marker, ess_tdof_list);
    
    if (my_rank == 0) {
        std::cout << "找到 " << ess_tdof_list.Size() << " 个本质边界DOF" << std::endl;
    }
    
    // 4. 创建心脏-躯干映射
    // 收集心脏边界点及其值
    std::map<int, double> heart_boundary_values;
    std::map<int, std::array<double, 3>> heart_coords;
    
    for (int i = 0; i < heart_mesh->GetNBE(); i++) {
        Array<int> vertices;
        heart_mesh->GetBdrElementVertices(i, vertices);
        for (int j = 0; j < vertices.Size(); j++) {
            int vertex_id = vertices[j];
            if (vertex_id < gf_ue_heart.Size()) {
                heart_boundary_values[vertex_id] = gf_ue_heart(vertex_id);
                double* coords = heart_mesh->GetVertex(vertex_id);
                heart_coords[vertex_id] = {coords[0], coords[1], coords[2]};
            }
        }
    }
    
    // 5. 创建边界值函数
    // 使用FunctionCoefficient来设置边界值
    class HeartTorsoBoundaryCoefficient : public Coefficient
    {
    private:
        const std::map<int, double>& heart_values_;
        const std::map<int, std::array<double, 3>>& heart_coords_;
        mfem::Mesh* torso_mesh_;
        double tolerance_;
        
    public:
        HeartTorsoBoundaryCoefficient(
            const std::map<int, double>& heart_values,
            const std::map<int, std::array<double, 3>>& heart_coords,
            mfem::Mesh* torso_mesh,
            double tolerance)
            : heart_values_(heart_values), heart_coords_(heart_coords),
              torso_mesh_(torso_mesh), tolerance_(tolerance) {}
        
        virtual double Eval(ElementTransformation& T, const IntegrationPoint& ip)
        {
            // 获取评估点的物理坐标
            Vector x(3);
            T.Transform(ip, x);
            
            // 查找最近的心脏点
            double min_distance = std::numeric_limits<double>::max();
            int closest_heart_vertex = -1;
            
            for (const auto& heart_pair : heart_coords_) {
                int heart_vertex_id = heart_pair.first;
                const auto& heart_pos = heart_pair.second;
                
                // 计算欧几里得距离
                double dx = heart_pos[0] - x(0);
                double dy = heart_pos[1] - x(1);
                double dz = heart_pos[2] - x(2);
                double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (distance < min_distance) {
                    min_distance = distance;
                    closest_heart_vertex = heart_vertex_id;
                }
            }
            
            // 如果找到足够近的点，返回心脏值，否则返回0
            if (min_distance <= tolerance_ && closest_heart_vertex != -1) {
                return heart_values_.at(closest_heart_vertex);
            } else {
                return 0.0;
            }
        }
    };
    
    // 创建边界系数
    HeartTorsoBoundaryCoefficient bdryCoef(
        heart_boundary_values, heart_coords, torso_mesh, tolerance);
    
    // 6. 将边界值投影到网格函数
    // 首先清零解向量
    gf_ue_torso = 0.0;
    
    // 投影边界系数到边界
    gf_ue_torso.ProjectBdrCoefficient(bdryCoef, bdr_marker);
    
    // 7. 验证边界值已设置
    // 获取真实自由度
    Vector ue_true(pfespace_torso->GetTrueVSize());
    gf_ue_torso.GetTrueDofs(ue_true);
    
    // 统计非零边界值
    int nonzero_count = 0;
    double min_val = 1e10, max_val = -1e10;
    
    for (int i = 0; i < ess_tdof_list.Size(); i++) {
        int tdof = ess_tdof_list[i];
        if (tdof >= 0 && tdof < ue_true.Size()) {
            double value = ue_true(tdof);
            if (std::abs(value) > 1e-10) {
                nonzero_count++;
                min_val = std::min(min_val, value);
                max_val = std::max(max_val, value);
            }
        }
    }
    
    // 收集统计信息
    int global_nonzero = 0;
    MPI_Allreduce(&nonzero_count, &global_nonzero, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    double global_min = min_val, global_max = max_val;
    if (nonzero_count == 0) {
        global_min = global_max = 0.0;
    } else {
        MPI_Allreduce(&min_val, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&max_val, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    
    if (my_rank == 0) {
        std::cout << "边界条件非零值数量: " << global_nonzero 
                 << ", 范围: [" << global_min << ", " << global_max << "]" << std::endl;
        
        // 输出一些样例值
        if (global_nonzero > 0) {
            std::cout << "边界值样例:" << std::endl;
            int sample_count = 0;
            for (int i = 0; i < std::min(10, ue_true.Size()); i++) {
                if (std::abs(ue_true(i)) > 1e-10) {
                    std::cout << "  DOF " << i << ": " << ue_true(i) << std::endl;
                    if (++sample_count >= 5) break;
                }
            }
        }
    }
}



void setHeartToTorsoBoundaryCorrectly(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    Array<int>& ess_tdof_list_torso,
    ParFiniteElementSpace* pfespace_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    double tolerance = 1e-6)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0) {
        std::cout << "设置心脏-躯干边界条件..." << std::endl;
    }
    
    // 1. 标记躯干网格中的心脏-躯干交界面
    // 创建与躯干边界元素数量相同的标记数组
    Array<int> is_heart_boundary(torso_mesh->GetNBE());
    is_heart_boundary = 0;  // 初始化为0
    
    // 将所有边界元素设置为1（心脏-躯干交界面）
    for (int i = 0; i < torso_mesh->GetNBE(); i++) {
        is_heart_boundary[i] = 1;
    }
    
    // 2. 设置边界条件标记数组
    // 设置所有边界元素属性为1
    for (int i = 0; i < torso_mesh->GetNBE(); i++) {
        torso_mesh->GetBdrElement(i)->SetAttribute(1);
    }
    torso_mesh->SetAttributes();
    
    // 3. 准备Dirichlet边界条件
    ParMesh* pmesh_torso = pfespace_torso->GetParMesh();
    int max_attr = pmesh_torso->bdr_attributes.Max();
    Array<int> ess_bdr(max_attr);
    ess_bdr = 0;
    ess_bdr[0] = 1;  // 标记属性1为Dirichlet边界
    
    // 4. 更新边界DOF列表
    pfespace_torso->GetEssentialTrueDofs(ess_bdr, ess_tdof_list_torso);
    
    if (my_rank == 0) {
        std::cout << "找到 " << ess_tdof_list_torso.Size() << " 个边界DOF" << std::endl;
    }
    
    // 5. 收集心脏边界点及其值
    std::map<int, double> heart_boundary_values;
    std::map<int, std::array<double, 3>> heart_coords;
    
    for (int i = 0; i < heart_mesh->GetNBE(); i++) {
        Array<int> vertices;
        heart_mesh->GetBdrElementVertices(i, vertices);
        for (int j = 0; j < vertices.Size(); j++) {
            int vertex_id = vertices[j];
            if (vertex_id < gf_ue_heart.Size()) {
                heart_boundary_values[vertex_id] = gf_ue_heart(vertex_id);
                double* coords = heart_mesh->GetVertex(vertex_id);
                heart_coords[vertex_id] = {coords[0], coords[1], coords[2]};
            }
        }
    }
    
    // 6. 对于每个躯干边界点，找到最近的心脏点
    int match_count = 0;
    
    // 创建Dirichlet边界条件函数 - 使用ParGridFunction而不是GridFunction
    ParGridFunction gf_torso_bdry(pfespace_torso);
    gf_torso_bdry = 0.0;  // 初始化为0
    
    // 遍历躯干边界顶点
    for (int i = 0; i < torso_mesh->GetNBE(); i++) {
        Array<int> vertices;
        torso_mesh->GetBdrElementVertices(i, vertices);
        
        for (int j = 0; j < vertices.Size(); j++) {
            int torso_vertex_id = vertices[j];
            if (torso_vertex_id >= gf_torso_bdry.Size()) continue;
            
            double* torso_coords = torso_mesh->GetVertex(torso_vertex_id);
            
            // 查找最近的心脏点
            double min_distance = std::numeric_limits<double>::max();
            int closest_heart_vertex = -1;
            
            for (const auto& heart_pair : heart_coords) {
                int heart_vertex_id = heart_pair.first;
                const auto& heart_pos = heart_pair.second;
                
                // 计算欧几里得距离
                double dx = heart_pos[0] - torso_coords[0];
                double dy = heart_pos[1] - torso_coords[1];
                double dz = heart_pos[2] - torso_coords[2];
                double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (distance < min_distance) {
                    min_distance = distance;
                    closest_heart_vertex = heart_vertex_id;
                }
            }
            
            // 如果找到足够近的点，设置值
            if (min_distance <= tolerance && closest_heart_vertex != -1) {
                double heart_value = heart_boundary_values[closest_heart_vertex];
                gf_torso_bdry(torso_vertex_id) = heart_value;
                match_count++;
            }
        }
    }
    
    // 并行同步边界值
    gf_torso_bdry.ParallelAssemble();
    
    // 7. 统计匹配点数量
    int global_match_count = 0;
    MPI_Allreduce(&match_count, &global_match_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (my_rank == 0) {
        std::cout << "设置了 " << global_match_count << " 个交界面点的值" << std::endl;
    }
    
    // 8. 显式地将边界值从网格函数转移到真实自由度
    Vector x_true(pfespace_torso->GetTrueVSize());
    gf_torso_bdry.GetTrueDofs(x_true);
    
    // 验证边界值
    int nonzero_count = 0;
    double min_val = 1e10, max_val = -1e10;
    
    for (int i = 0; i < ess_tdof_list_torso.Size(); i++) {
        int tdof = ess_tdof_list_torso[i];
        if (tdof >= 0 && tdof < x_true.Size()) {
            double value = x_true(tdof);
            if (std::abs(value) > 1e-10) {
                nonzero_count++;
                min_val = std::min(min_val, value);
                max_val = std::max(max_val, value);
            }
        }
    }
    
    // 收集全局统计
    int global_nonzero = 0;
    MPI_Allreduce(&nonzero_count, &global_nonzero, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    double global_min = min_val, global_max = max_val;
    if (nonzero_count == 0) {
        global_min = global_max = 0.0;
    } else {
        MPI_Allreduce(&min_val, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&max_val, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    
    if (my_rank == 0) {
        std::cout << "边界条件非零值数量: " << global_nonzero 
                 << ", 范围: [" << global_min << ", " << global_max << "]" << std::endl;
    }
    
    // 9. 最关键的一步：将边界条件函数赋值给躯干解
    gf_ue_torso = gf_torso_bdry;
}

/**
 * 设置心脏到躯干的边界条件（改进版）
 * 直接操作真实自由度以确保边界值被正确设置
 */
void setHeartToTorsoBoundaryValues(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    Array<int>& ess_tdof_list_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    ParFiniteElementSpace* pfespace_torso,
    double tolerance = 1e-6,
    bool debug = true)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0 && debug) {
        std::cout << "设置心脏-躯干边界值..." << std::endl;
        std::cout << "心脏解范围: [" << gf_ue_heart.Min() << ", " << gf_ue_heart.Max() << "]" << std::endl;
    }

    // 1. 收集心脏边界点及其值
    std::map<int, double> heart_boundary_values;
    std::map<int, std::array<double, 3>> heart_coords;
    
    for (int i = 0; i < heart_mesh->GetNBE(); i++) {
        Array<int> vertices;
        heart_mesh->GetBdrElementVertices(i, vertices);
        for (int j = 0; j < vertices.Size(); j++) {
            int vertex_id = vertices[j];
            if (vertex_id < gf_ue_heart.Size()) {
                heart_boundary_values[vertex_id] = gf_ue_heart(vertex_id);
                double* coords = heart_mesh->GetVertex(vertex_id);
                heart_coords[vertex_id] = {coords[0], coords[1], coords[2]};
            }
        }
    }
    
    // 2. 对于每个躯干边界点，找到最近的心脏点
    std::map<int, int> torso_to_heart_map;
    
    // 这一步只在主进程进行，然后广播结果
    if (my_rank == 0) {
        // 首先收集所有躯干边界点
        std::map<int, std::array<double, 3>> torso_coords;
        for (int i = 0; i < torso_mesh->GetNBE(); i++) {
            Array<int> vertices;
            torso_mesh->GetBdrElementVertices(i, vertices);
            for (int j = 0; j < vertices.Size(); j++) {
                int vertex_id = vertices[j];
                double* coords = torso_mesh->GetVertex(vertex_id);
                torso_coords[vertex_id] = {coords[0], coords[1], coords[2]};
            }
        }
        
        // 对每个躯干边界点，找最近的心脏点
        for (const auto& torso_pair : torso_coords) {
            int torso_id = torso_pair.first;
            const auto& torso_pos = torso_pair.second;
            
            double min_dist = 1e10;
            int closest_heart_id = -1;
            
            for (const auto& heart_pair : heart_coords) {
                int heart_id = heart_pair.first;
                const auto& heart_pos = heart_pair.second;
                
                double dist = std::sqrt(
                    std::pow(torso_pos[0] - heart_pos[0], 2) +
                    std::pow(torso_pos[1] - heart_pos[1], 2) +
                    std::pow(torso_pos[2] - heart_pos[2], 2)
                );
                
                if (dist < min_dist) {
                    min_dist = dist;
                    closest_heart_id = heart_id;
                }
            }
            
            if (min_dist <= tolerance && closest_heart_id >= 0) {
                torso_to_heart_map[torso_id] = closest_heart_id;
            }
        }
        
        if (debug) {
            std::cout << "找到 " << torso_to_heart_map.size() << " 个匹配的心脏-躯干边界点" << std::endl;
        }
    }
    
    // 3. 广播映射关系到所有进程
    // 首先广播映射大小
    int map_size = torso_to_heart_map.size();
    MPI_Bcast(&map_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // 然后广播映射内容
    std::vector<int> torso_vertices(map_size), heart_vertices(map_size);
    if (my_rank == 0) {
        int i = 0;
        for (const auto& pair : torso_to_heart_map) {
            torso_vertices[i] = pair.first;
            heart_vertices[i] = pair.second;
            i++;
        }
    }
    
    MPI_Bcast(torso_vertices.data(), map_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(heart_vertices.data(), map_size, MPI_INT, 0, MPI_COMM_WORLD);
    
    // 非主进程重建映射
    if (my_rank != 0) {
        for (int i = 0; i < map_size; i++) {
            torso_to_heart_map[torso_vertices[i]] = heart_vertices[i];
        }
    }
    
    // 4. 广播心脏边界值到所有进程
    // 首先广播心脏边界点数量
    int heart_bd_size = heart_boundary_values.size();
    MPI_Bcast(&heart_bd_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // 然后广播心脏边界值
    std::vector<int> heart_ids(heart_bd_size);
    std::vector<double> heart_vals(heart_bd_size);
    
    if (my_rank == 0) {
        int i = 0;
        for (const auto& pair : heart_boundary_values) {
            heart_ids[i] = pair.first;
            heart_vals[i] = pair.second;
            i++;
        }
    }
    
    MPI_Bcast(heart_ids.data(), heart_bd_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(heart_vals.data(), heart_bd_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // 非主进程重建心脏边界值
    if (my_rank != 0) {
        for (int i = 0; i < heart_bd_size; i++) {
            heart_boundary_values[heart_ids[i]] = heart_vals[i];
        }
    }
    
    // 5. 设置边界条件
    // 首先获取当前解向量
    Vector true_dofs(pfespace_torso->GetTrueVSize());
    gf_ue_torso.GetTrueDofs(true_dofs);
    
    // 记录已修改的DOF
    int modified_dofs = 0;
    
    // 获取局部到全局DOF映射
    Array<int> vdofs;
    for (const auto& pair : torso_to_heart_map) {
        int torso_id = pair.first;
        int heart_id = pair.second;
        
        // 确保心脏ID在值字典中
        if (heart_boundary_values.find(heart_id) != heart_boundary_values.end()) {
            double heart_val = heart_boundary_values[heart_id];
            
            // 获取与躯干顶点相关的DOF
            pfespace_torso->GetVertexDofs(torso_id, vdofs);
            
            // 设置所有相关DOF的值
            for (int i = 0; i < vdofs.Size(); i++) {
                int vdof = vdofs[i];
                if (vdof >= 0) {  // 本地DOF
                    int tdof = pfespace_torso->GetLocalTDofNumber(vdof);
                    if (tdof >= 0 && tdof < true_dofs.Size()) {
                        true_dofs(tdof) = heart_val;
                        modified_dofs++;
                    }
                }
            }
        }
    }
    
    // 6. 从真实自由度更新解
    gf_ue_torso.SetFromTrueDofs(true_dofs);
    
    // 7. 统计修改的DOF
    int global_modified = 0;
    MPI_Allreduce(&modified_dofs, &global_modified, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (my_rank == 0 && debug) {
        std::cout << "修改了 " << global_modified << " 个交界面DOF的值" << std::endl;
        
        // 检查设置后的解范围
        Vector updated_dofs(pfespace_torso->GetTrueVSize());
        gf_ue_torso.GetTrueDofs(updated_dofs);
        
        int nonzero_count = 0;
        double min_val = 1e10, max_val = -1e10;
        
        for (int i = 0; i < updated_dofs.Size(); i++) {
            if (fabs(updated_dofs(i)) > 1e-10) {
                nonzero_count++;
                min_val = std::min(min_val, updated_dofs(i));
                max_val = std::max(max_val, updated_dofs(i));
            }
        }
        
        std::cout << "解中非零值数量: " << nonzero_count 
                 << ", 范围: [" << min_val << ", " << max_val << "]" << std::endl;
    }
    
    // 8. 根据修改后的解重新计算ess_tdof_list（如果需要）
    if (global_modified > 0) {
        // 给ess_tdof_list添加所有我们直接设置了值的DOF
        std::set<int> ess_dofs_set;
        
        // 添加现有的DOF
        for (int i = 0; i < ess_tdof_list_torso.Size(); i++) {
            ess_dofs_set.insert(ess_tdof_list_torso[i]);
        }
        
        // 对于每个躯干-心脏匹配点，添加对应的DOF
        for (const auto& pair : torso_to_heart_map) {
            int torso_id = pair.first;
            
            // 获取与躯干顶点相关的DOF
            pfespace_torso->GetVertexDofs(torso_id, vdofs);
            
            // 添加所有相关DOF
            for (int i = 0; i < vdofs.Size(); i++) {
                int vdof = vdofs[i];
                if (vdof >= 0) {  // 本地DOF
                    int tdof = pfespace_torso->GetLocalTDofNumber(vdof);
                    if (tdof >= 0) {
                        ess_dofs_set.insert(tdof);
                    }
                }
            }
        }
        
        // 转换回Array
        ess_tdof_list_torso.SetSize(ess_dofs_set.size());
        int idx = 0;
        for (int dof : ess_dofs_set) {
            ess_tdof_list_torso[idx++] = dof;
        }
        
        if (my_rank == 0 && debug) {
            std::cout << "更新后的边界DOF列表大小: " << ess_tdof_list_torso.Size() << std::endl;
        }
    }
}

/**
 * 将心脏解值传递到躯干边界并输出详细调试信息
 * 
 * @param gf_ue_heart 心脏解
 * @param gf_ue_torso 躯干解
 * @param ess_tdof_list_torso 躯干边界DOF列表
 * @param heart_mesh 心脏网格
 * @param torso_mesh 躯干网格
 * @param pfespace_torso 躯干有限元空间
 * @param tolerance 坐标匹配容差
 * @param heart_bdry_attr 心脏边界属性(用于标识接口)
 * @param outputDir 输出目录
 * @param timestep 当前时间步
 */
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
    int timestep = 0)
{
    int my_rank, num_ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    
    if (my_rank == 0) {
        std::cout << "详细调试: 设置心脏-躯干边界值..." << std::endl;
        std::cout << "心脏网格顶点数: " << heart_mesh->GetNV() << std::endl;
        std::cout << "心脏网格边界元素数: " << heart_mesh->GetNBE() << std::endl;
        std::cout << "躯干网格顶点数: " << torso_mesh->GetNV() << std::endl;
        std::cout << "躯干网格边界元素数: " << torso_mesh->GetNBE() << std::endl;
        std::cout << "躯干有限元空间自由度: " << pfespace_torso->GetTrueVSize() << std::endl;
        std::cout << "边界DOF列表大小: " << ess_tdof_list_torso.Size() << std::endl;
    }
    
    // 1. 检查心脏解的值范围
    double heart_min = gf_ue_heart.Min();
    double heart_max = gf_ue_heart.Max();
    if (my_rank == 0) {
        std::cout << "心脏解范围: [" << heart_min << ", " << heart_max << "]" << std::endl;
    }
    
    // 2. 收集心脏边界上的顶点和电位值
    std::map<int, double> heart_boundary_values;
    std::vector<int> heart_bdry_vertices;
    
    for (int i = 0; i < heart_mesh->GetNBE(); i++) {
        Array<int> vertices;
        heart_mesh->GetBdrElementVertices(i, vertices);
        for (int j = 0; j < vertices.Size(); j++) {
            int vertex_id = vertices[j];
            if (vertex_id < gf_ue_heart.Size()) {
                heart_boundary_values[vertex_id] = gf_ue_heart(vertex_id);
                heart_bdry_vertices.push_back(vertex_id);
            }
        }
    }
    
    // 输出一些心脏边界点的值
    if (my_rank == 0) {
        std::cout << "收集到 " << heart_boundary_values.size() << " 个心脏边界点" << std::endl;
        
        // 统计非零值数量
        int nonzero_heart_vals = 0;
        double min_val = 1e10, max_val = -1e10;
        
        for (const auto& pair : heart_boundary_values) {
            if (fabs(pair.second) > 1e-10) {
                nonzero_heart_vals++;
                min_val = std::min(min_val, pair.second);
                max_val = std::max(max_val, pair.second);
            }
        }
        
        std::cout << "心脏边界点非零值数量: " << nonzero_heart_vals 
                 << ", 范围: [" << min_val << ", " << max_val << "]" << std::endl;
                 
        // 输出一些采样点的值
        std::cout << "心脏边界点值样例:" << std::endl;
        int count = 0;
        for (const auto& pair : heart_boundary_values) {
            if (count++ < 10) {
                std::cout << "  顶点 " << pair.first << ": " << pair.second << std::endl;
            } else {
                break;
            }
        }
    }
    
    // 3. 准备接收心脏值的数据结构
    std::vector<std::tuple<int, double, double, double, double>> heart_points;
    for (const auto& pair : heart_boundary_values) {
        int vertex_id = pair.first;
        double value = pair.second;
        double* coords = heart_mesh->GetVertex(vertex_id);
        heart_points.push_back(std::make_tuple(
            vertex_id, coords[0], coords[1], coords[2], value
        ));
    }
    
    // 4. 找到躯干网格上的接口顶点
    int match_count = 0;
    std::vector<std::tuple<int, double, double, double, double, double>> matched_points;
    std::set<int> matched_torso_vertices;
    
    // 遍历躯干边界元素
    for (int i = 0; i < torso_mesh->GetNBE(); i++) {
        // 检查是否是接口边界
        if (torso_mesh->GetBdrAttribute(i) == heart_bdry_attr) {
            Array<int> vertices;
            torso_mesh->GetBdrElementVertices(i, vertices);
            
            for (int j = 0; j < vertices.Size(); j++) {
                int torso_vertex_id = vertices[j];
                if (torso_vertex_id >= gf_ue_torso.Size()) continue;
                
                matched_torso_vertices.insert(torso_vertex_id);
                double* torso_coords = torso_mesh->GetVertex(torso_vertex_id);
                
                // 查找最近的心脏点
                double min_distance = std::numeric_limits<double>::max();
                int closest_heart_vertex = -1;
                double heart_value = 0.0;
                
                for (const auto& heart_point : heart_points) {
                    double heart_x = std::get<1>(heart_point);
                    double heart_y = std::get<2>(heart_point);
                    double heart_z = std::get<3>(heart_point);
                    
                    double distance = std::sqrt(
                        std::pow(heart_x - torso_coords[0], 2) +
                        std::pow(heart_y - torso_coords[1], 2) +
                        std::pow(heart_z - torso_coords[2], 2)
                    );
                    
                    if (distance < min_distance) {
                        min_distance = distance;
                        closest_heart_vertex = std::get<0>(heart_point);
                        heart_value = std::get<4>(heart_point);
                    }
                }
                
                // 如果找到足够近的点，设置值
                if (min_distance <= tolerance && closest_heart_vertex != -1) {
                    double old_value = gf_ue_torso(torso_vertex_id);
                    gf_ue_torso(torso_vertex_id) = heart_value;
                    match_count++;
                    
                    // 存储匹配信息用于调试
                    matched_points.push_back(std::make_tuple(
                        torso_vertex_id, 
                        torso_coords[0], torso_coords[1], torso_coords[2],
                        old_value, heart_value
                    ));
                }
            }
        }
    }
    
    // 同步所有进程上的值
    gf_ue_torso.ParallelAssemble();
    
    // 5. 统计匹配点数量
    int global_match_count = 0;
    MPI_Allreduce(&match_count, &global_match_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (my_rank == 0) {
        std::cout << "设置了 " << global_match_count << " 个交界面点的值" << std::endl;
        
        // 输出匹配点详情
        if (!matched_points.empty()) {
            std::cout << "匹配点样例:" << std::endl;
            int count = 0;
            for (const auto& point : matched_points) {
                if (count++ < 10) {
                    std::cout << "  躯干顶点 " << std::get<0>(point) 
                             << " (" << std::get<1>(point) << ", " << std::get<2>(point) << ", " << std::get<3>(point) << ")"
                             << ", 原值: " << std::get<4>(point)
                             << ", 新值: " << std::get<5>(point) << std::endl;
                } else {
                    break;
                }
            }
        }
    }
    
    // 6. 检查边界DOF是否正确设置
    std::vector<int> dof_to_vertex(pfespace_torso->GetVSize(), -1);
    for (int i = 0; i < matched_torso_vertices.size(); i++) {
        int vertex_id = *std::next(matched_torso_vertices.begin(), i);
        Array<int> dofs;
        pfespace_torso->GetVertexDofs(vertex_id, dofs);
        for (int j = 0; j < dofs.Size(); j++) {
            if (dofs[j] >= 0 && dofs[j] < dof_to_vertex.size()) {
                dof_to_vertex[dofs[j]] = vertex_id;
            }
        }
    }
    
    // 7. 验证ess_tdof_list中的DOF是否存在于匹配的顶点中
    int valid_dofs = 0;
    std::vector<int> valid_torso_dofs;
    
    for (int i = 0; i < ess_tdof_list_torso.Size(); i++) {
        int tdof = ess_tdof_list_torso[i];
        int ldof = pfespace_torso->GetLocalTDofNumber(tdof);
        if (ldof >= 0 && ldof < dof_to_vertex.size()) {
            int vertex_id = dof_to_vertex[ldof];
            if (vertex_id >= 0 && matched_torso_vertices.count(vertex_id) > 0) {
                valid_dofs++;
                valid_torso_dofs.push_back(tdof);
            }
        }
    }
    
    // 统计有效DOF
    int global_valid_dofs = 0;
    MPI_Allreduce(&valid_dofs, &global_valid_dofs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (my_rank == 0) {
        std::cout << "边界DOF列表中有效交界面DOF数量: " << global_valid_dofs 
                 << " (总计: " << ess_tdof_list_torso.Size() << ")" << std::endl;
    }
    
    // 8. 输出到VTK文件以便可视化
    if (my_rank == 0) {
        // 创建目录
        std::string debug_dir = outputDir + "/debug";
        std::string cmd = "mkdir -p " + debug_dir;
        system(cmd.c_str());
        
        // 导出心脏边界点
        std::string heart_points_file = debug_dir + "/heart_boundary_" + std::to_string(timestep) + ".vtk";
        std::ofstream heart_out(heart_points_file);
        
        heart_out << "# vtk DataFile Version 3.0\n";
        heart_out << "Heart boundary points\n";
        heart_out << "ASCII\n";
        heart_out << "DATASET POLYDATA\n";
        heart_out << "POINTS " << heart_points.size() << " float\n";
        
        for (const auto& point : heart_points) {
            heart_out << std::get<1>(point) << " " << std::get<2>(point) << " " << std::get<3>(point) << "\n";
        }
        
        heart_out << "POINT_DATA " << heart_points.size() << "\n";
        heart_out << "SCALARS value float 1\n";
        heart_out << "LOOKUP_TABLE default\n";
        
        for (const auto& point : heart_points) {
            heart_out << std::get<4>(point) << "\n";
        }
        
        heart_out.close();
        
        // 导出匹配的躯干点
        std::string torso_points_file = debug_dir + "/torso_boundary_" + std::to_string(timestep) + ".vtk";
        std::ofstream torso_out(torso_points_file);
        
        torso_out << "# vtk DataFile Version 3.0\n";
        torso_out << "Torso boundary points\n";
        torso_out << "ASCII\n";
        torso_out << "DATASET POLYDATA\n";
        torso_out << "POINTS " << matched_points.size() << " float\n";
        
        for (const auto& point : matched_points) {
            torso_out << std::get<1>(point) << " " << std::get<2>(point) << " " << std::get<3>(point) << "\n";
        }
        
        torso_out << "POINT_DATA " << matched_points.size() << "\n";
        torso_out << "SCALARS value float 1\n";
        torso_out << "LOOKUP_TABLE default\n";
        
        for (const auto& point : matched_points) {
            torso_out << std::get<5>(point) << "\n";
        }
        
        torso_out.close();
        
        std::cout << "已导出调试VTK文件到: " << debug_dir << std::endl;
    }
    
    // 9. 检查设置后的边界值
if (my_rank == 0) {
    std::cout << "检查设置后的边界值:" << std::endl;
    Vector true_dofs;
    gf_ue_torso.GetTrueDofs(true_dofs);
    
    int nonzero_values = 0;
    double min_val = 1e10, max_val = -1e10;
    
    for (int i = 0; i < ess_tdof_list_torso.Size(); i++) {
        int tdof = ess_tdof_list_torso[i];
        if (tdof >= 0 && tdof < true_dofs.Size()) {
            double value = true_dofs(tdof);
            if (fabs(value) > 1e-10) {
                nonzero_values++;
                min_val = std::min(min_val, value);
                max_val = std::max(max_val, value);
            }
        }
    }
    
    std::cout << "设置后非零边界值数量: " << nonzero_values 
             << ", 范围: [" << min_val << ", " << max_val << "]" << std::endl;
}
    
    // 10. 确认边界条件已经传递到解向量
    Vector u_e_torso_true(pfespace_torso->GetTrueVSize());
    gf_ue_torso.GetTrueDofs(u_e_torso_true);
    
    int nonzero_in_vector = 0;
    double min_val = 1e10, max_val = -1e10;
    
    for (int i = 0; i < u_e_torso_true.Size(); i++) {
        if (fabs(u_e_torso_true(i)) > 1e-10) {
            nonzero_in_vector++;
            min_val = std::min(min_val, u_e_torso_true(i));
            max_val = std::max(max_val, u_e_torso_true(i));
        }
    }
    
    int global_nonzero = 0;
    MPI_Allreduce(&nonzero_in_vector, &global_nonzero, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (my_rank == 0) {
        std::cout << "解向量中非零值数量: " << global_nonzero 
                 << ", 范围: [" << min_val << ", " << max_val << "]" << std::endl;
    }
}



/**
 * 获取心脏-躯干接口的边界DOF列表
 * 安全版本 - 处理并行分区
 * 
 * @param pfespace_torso 躯干有限元空间
 * @param heart_mesh 心脏网格
 * @param torso_mesh 躯干网格
 * @param tolerance 坐标匹配容差
 * @return 接口上的DOF列表
 */
Array<int> getHeartTorsoInterfaceDofs(
    ParFiniteElementSpace* pfespace_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    double tolerance = 1e-6)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0) {
        std::cout << "计算心脏-躯干接口边界DOF..." << std::endl;
    }
    
    // 获取并行网格
    ParMesh* pmesh_torso = pfespace_torso->GetParMesh();
    
    // 创建边界属性标记数组
    int max_attr = pmesh_torso->bdr_attributes.Max();
    Array<int> ess_bdr(max_attr);
    ess_bdr = 0;  // 初始化为零
    
    // 标记心脏-躯干接口的属性
    int heart_bdry_attr = 1;  // 假设接口标记为1
    if (heart_bdry_attr <= max_attr) {
        ess_bdr[heart_bdry_attr - 1] = 1;
    } else {
        if (my_rank == 0) {
            std::cout << "警告: 接口属性值 " << heart_bdry_attr 
                      << " 超出最大边界属性值 " << max_attr << std::endl;
        }
    }
    
    // 使用MFEM内置方法获取边界DOF
    Array<int> ess_tdof_list;
    pfespace_torso->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    
    // 统计并输出边界DOF数量
    int local_dofs = ess_tdof_list.Size();
    int global_dofs = 0;
    MPI_Allreduce(&local_dofs, &global_dofs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (my_rank == 0) {
        std::cout << "找到 " << global_dofs << " 个交界面DOF" << std::endl;
    }
    
    return ess_tdof_list;
}

/**
 * 将心脏解值传递到躯干边界
 * 
 * @param gf_ue_heart 心脏解
 * @param gf_ue_torso 躯干解
 * @param heart_mesh 心脏网格
 * @param torso_mesh 躯干网格
 * @param tolerance 坐标匹配容差
 * @param heart_bdry_attr 心脏边界属性(用于标识接口)
 */
void transferHeartToTorsoBoundary(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    double tolerance = 1e-6,
    int heart_bdry_attr = 1)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    // 1. 收集心脏边界上的顶点和电位值
    std::map<int, double> heart_boundary_values;
    
    for (int i = 0; i < heart_mesh->GetNBE(); i++) {
        Array<int> vertices;
        heart_mesh->GetBdrElementVertices(i, vertices);
        for (int j = 0; j < vertices.Size(); j++) {
            int vertex_id = vertices[j];
            if (vertex_id < gf_ue_heart.Size()) {
                heart_boundary_values[vertex_id] = gf_ue_heart(vertex_id);
            }
        }
    }
    
    // 2. 准备接收心脏值的数据结构
    std::vector<std::tuple<int, double, double, double>> heart_points;
    for (const auto& pair : heart_boundary_values) {
        int vertex_id = pair.first;
        double* coords = heart_mesh->GetVertex(vertex_id);
        heart_points.push_back(std::make_tuple(
            vertex_id, coords[0], coords[1], coords[2]
        ));
    }
    
    // 3. 找到躯干网格上的接口顶点并设置边界值
    int match_count = 0;
    
    // 获取躯干边界元素
    for (int i = 0; i < torso_mesh->GetNBE(); i++) {
        // 检查是否是接口边界
        if (torso_mesh->GetBdrAttribute(i) == heart_bdry_attr) {
            Array<int> vertices;
            torso_mesh->GetBdrElementVertices(i, vertices);
            
            for (int j = 0; j < vertices.Size(); j++) {
                int torso_vertex_id = vertices[j];
                if (torso_vertex_id >= gf_ue_torso.Size()) continue;
                
                double* torso_coords = torso_mesh->GetVertex(torso_vertex_id);
                
                // 查找最近的心脏点
                double min_distance = std::numeric_limits<double>::max();
                int closest_heart_vertex = -1;
                
                for (const auto& heart_point : heart_points) {
                    double heart_x = std::get<1>(heart_point);
                    double heart_y = std::get<2>(heart_point);
                    double heart_z = std::get<3>(heart_point);
                    
                    double distance = std::sqrt(
                        std::pow(heart_x - torso_coords[0], 2) +
                        std::pow(heart_y - torso_coords[1], 2) +
                        std::pow(heart_z - torso_coords[2], 2)
                    );
                    
                    if (distance < min_distance) {
                        min_distance = distance;
                        closest_heart_vertex = std::get<0>(heart_point);
                    }
                }
                
                // 如果找到足够近的点，设置值
                if (min_distance <= tolerance && closest_heart_vertex != -1) {
                    gf_ue_torso(torso_vertex_id) = heart_boundary_values[closest_heart_vertex];
                    match_count++;
                }
            }
        }
    }
    
    // 同步所有进程上的值
    gf_ue_torso.ParallelAssemble();
    
    // 统计匹配点数量
    int global_match_count = 0;
    MPI_Allreduce(&match_count, &global_match_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (my_rank == 0) {
        std::cout << "设置了 " << global_match_count << " 个交界面点的值" << std::endl;
    }
}


/**
 * 建立心脏-躯干网格边界映射关系
 * 
 * @param heart_mesh 心脏网格
 * @param torso_mesh 躯干网格
 * @param tolerance 坐标匹配容差
 * @return 返回从躯干网格顶点ID到心脏网格顶点ID的映射
 */
std::map<int, int> buildHeartTorsoMapping(
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    double tolerance = 1e-6)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0) {
        std::cout << "构建心脏-躯干边界映射关系..." << std::endl;
    }
    
    std::map<int, int> torso_to_heart_map;
    
    // 1. 收集心脏边界顶点和坐标
    std::set<int> heart_boundary_vertices;
    std::map<int, double*> heart_boundary_coords;
    
    for (int i = 0; i < heart_mesh->GetNBE(); i++) {
        Array<int> vertices;
        heart_mesh->GetBdrElementVertices(i, vertices);
        for (int j = 0; j < vertices.Size(); j++) {
            int vertex_id = vertices[j];
            heart_boundary_vertices.insert(vertex_id);
            heart_boundary_coords[vertex_id] = heart_mesh->GetVertex(vertex_id);
        }
    }
    
    // 2. 收集躯干边界顶点
    std::set<int> torso_boundary_vertices;
    for (int i = 0; i < torso_mesh->GetNBE(); i++) {
        Array<int> vertices;
        torso_mesh->GetBdrElementVertices(i, vertices);
        for (int j = 0; j < vertices.Size(); j++) {
            torso_boundary_vertices.insert(vertices[j]);
        }
    }
    
    // 3. 对于每个躯干边界顶点，查找最近的心脏边界顶点
    int match_count = 0;
    for (int torso_vertex_id : torso_boundary_vertices) {
        double* torso_coords = torso_mesh->GetVertex(torso_vertex_id);
        
        int closest_heart_vertex = -1;
        double min_distance = tolerance * 2;  // 初始设为大于容差
        
        for (int heart_vertex_id : heart_boundary_vertices) {
            double* heart_coords = heart_boundary_coords[heart_vertex_id];
            
            // 计算欧几里得距离
            double distance = 0.0;
            for (int k = 0; k < 3; k++) {
                double diff = heart_coords[k] - torso_coords[k];
                distance += diff * diff;
            }
            distance = std::sqrt(distance);
            
            if (distance < min_distance) {
                min_distance = distance;
                closest_heart_vertex = heart_vertex_id;
            }
        }
        
        // 如果找到足够近的匹配点，添加到映射
        if (min_distance <= tolerance && closest_heart_vertex != -1) {
            torso_to_heart_map[torso_vertex_id] = closest_heart_vertex;
            match_count++;
        }
    }
    
    // 收集所有进程的匹配点数量
    int global_match_count = 0;
    MPI_Allreduce(&match_count, &global_match_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (my_rank == 0) {
        std::cout << "找到 " << global_match_count << " 个心脏-躯干交界面匹配点" << std::endl;
    }
    
    return torso_to_heart_map;
}

/**
 * 获取Dirichlet边界条件对应的DOF列表
 * 
 * @param torso_to_heart_map 从躯干顶点到心脏顶点的映射
 * @param pfespace_torso 躯干有限元空间
 * @return 返回Dirichlet边界DOF列表
 */
Array<int> getInterfaceDofs(
    const std::map<int, int>& torso_to_heart_map, 
    ParFiniteElementSpace* pfespace_torso)
{
    std::set<int> dof_set;
    
    // 为每个映射顶点获取对应的DOF
    Array<int> vdofs;
    for (const auto& pair : torso_to_heart_map) {
        int torso_vertex_id = pair.first;
        
        pfespace_torso->GetVertexDofs(torso_vertex_id, vdofs);
        for (int i = 0; i < vdofs.Size(); i++) {
            int vdof = vdofs[i];
            if (vdof >= 0) {  // 本地拥有的DOF
                int tdof = pfespace_torso->GetLocalTDofNumber(vdof);
                if (tdof >= 0) {
                    dof_set.insert(tdof);
                }
            }
        }
    }
    
    // 转换为Array
    Array<int> dofs(dof_set.size());
    int idx = 0;
    for (int dof : dof_set) {
        dofs[idx++] = dof;
    }
    
    return dofs;
}

/**
 * 从心脏解传递值到躯干边界
 * 
 * @param torso_to_heart_map 从躯干顶点到心脏顶点的映射
 * @param gf_ue_heart 心脏解
 * @param gf_ue_torso 躯干解
 */
void transferHeartValuesToTorso(
    const std::map<int, int>& torso_to_heart_map,
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0 && torso_to_heart_map.size() > 0) {
        std::cout << "传递心脏解值到躯干边界..." << std::endl;
    }
    
    // 为每个映射顶点设置值
    for (const auto& pair : torso_to_heart_map) {
        int torso_vertex_id = pair.first;
        int heart_vertex_id = pair.second;
        
        if (torso_vertex_id < gf_ue_torso.Size() && heart_vertex_id < gf_ue_heart.Size()) {
            gf_ue_torso(torso_vertex_id) = gf_ue_heart(heart_vertex_id);
        }
    }
    
    // 同步所有进程上的值
    gf_ue_torso.ParallelAssemble();
}

/**
 * 改进版边界相交识别和条件设置
 * 使用精确的几何条件和距离检查来识别真正的交界面点
 * 
 * @param gf_ue_heart 心脏解
 * @param gf_ue_torso 躯干解
 * @param heart_mesh 心脏网格
 * @param torso_mesh 躯干网格
 * @param distance_threshold 距离阈值，只有小于此距离的点才被认为是匹配的
 * @param debug_output 是否输出调试信息
 * @param max_debug_points 最大调试输出点数
 */
void setIntersectionBoundary(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh, 
    double distance_threshold = 1e-6,
    bool debug_output = true,
    int max_debug_points = 980)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0) {
        std::cout << "使用改进算法识别心脏-躯干交界面..." << std::endl;
    }
    
    // 步骤1: 收集心脏表面点的坐标和对应值
    std::map<int, double> heart_boundary_values;  // 顶点ID到值的映射
    
    // 收集心脏边界顶点
    std::set<int> heart_boundary_vertices;
    for (int i = 0; i < heart_mesh->GetNBE(); i++) {
        Array<int> vertices;
        heart_mesh->GetBdrElementVertices(i, vertices);
        for (int j = 0; j < vertices.Size(); j++) {
            heart_boundary_vertices.insert(vertices[j]);
        }
    }
    
    // 获取这些顶点对应的值
    for (int vertex_id : heart_boundary_vertices) {
        if (vertex_id < gf_ue_heart.Size()) {
            heart_boundary_values[vertex_id] = gf_ue_heart(vertex_id);
        }
    }
    
    if (my_rank == 0) {
        std::cout << "收集到心脏边界顶点数量: " << heart_boundary_vertices.size() << std::endl;
    }
    
    // 准备调试输出文件
    std::ofstream debug_file;
    if (debug_output && my_rank == 0) {
        debug_file.open("boundary_points_matching.csv");
        debug_file << "匹配类型,心脏顶点ID,躯干顶点ID,X坐标,Y坐标,Z坐标,距离,心脏值,躯干值" << std::endl;
    }
    
    // 步骤2: 构建心脏表面点的kd树或其他空间索引
    // 这里使用简化的方法：存储所有心脏边界点的坐标，并在后续计算实际距离
    std::vector<std::tuple<int, double, double, double>> heart_boundary_points;
    
    for (int vertex_id : heart_boundary_vertices) {
        double* coords = heart_mesh->GetVertex(vertex_id);
        heart_boundary_points.push_back(std::make_tuple(
            vertex_id, coords[0], coords[1], coords[2]
        ));
    }
    
    // 步骤3: 识别躯干网格上的可能交界面点
    std::set<int> torso_boundary_vertices;
    for (int i = 0; i < torso_mesh->GetNBE(); i++) {
        Array<int> vertices;
        torso_mesh->GetBdrElementVertices(i, vertices);
        for (int j = 0; j < vertices.Size(); j++) {
            torso_boundary_vertices.insert(vertices[j]);
        }
    }
    
    if (my_rank == 0) {
        std::cout << "收集到躯干边界顶点数量: " << torso_boundary_vertices.size() << std::endl;
    }
    
    // 步骤4: 对于每个躯干边界点，查找距离最近的心脏边界点
    int matched_count = 0;
    int debug_count = 0;
    
    // 调试输出的匹配点信息结构
    struct MatchPoint {
        int heart_id;
        int torso_id;
        double x, y, z;
        double distance;
        double heart_val;
        double torso_val;
    };
    std::vector<MatchPoint> match_points;
    std::vector<std::tuple<int, double, double, double>> unmatched_heart_points;
    
    // 只查找躯干边界点
    for (int torso_vertex_id : torso_boundary_vertices) {
        if (torso_vertex_id >= gf_ue_torso.Size()) continue;
        
        double* torso_coords = torso_mesh->GetVertex(torso_vertex_id);
        
        // 查找最近的心脏边界点
        double min_distance = std::numeric_limits<double>::max();
        int closest_heart_vertex = -1;
        
        for (const auto& heart_point : heart_boundary_points) {
            int heart_vertex_id = std::get<0>(heart_point);
            double heart_x = std::get<1>(heart_point);
            double heart_y = std::get<2>(heart_point);
            double heart_z = std::get<3>(heart_point);
            
            // 计算欧几里得距离
            double distance = std::sqrt(
                std::pow(heart_x - torso_coords[0], 2) +
                std::pow(heart_y - torso_coords[1], 2) +
                std::pow(heart_z - torso_coords[2], 2)
            );
            
            if (distance < min_distance) {
                min_distance = distance;
                closest_heart_vertex = heart_vertex_id;
            }
        }
        
        // 如果找到的最近点足够近，则认为是匹配点
        if (min_distance <= distance_threshold && closest_heart_vertex != -1) {
            double heart_value = heart_boundary_values[closest_heart_vertex];
            double torso_value_before = gf_ue_torso(torso_vertex_id);
            
            // 设置躯干点的值为对应心脏点的值
            gf_ue_torso(torso_vertex_id) = heart_value;
            matched_count++;
            
            // 收集调试信息
            if (debug_output && debug_count < max_debug_points) {
                match_points.push_back({
                    closest_heart_vertex, torso_vertex_id,
                    torso_coords[0], torso_coords[1], torso_coords[2],
                    min_distance,
                    heart_value, torso_value_before
                });
                debug_count++;
            }
        }
    }
    
    // 收集一些未匹配的心脏边界点用于调试
    if (debug_output && my_rank == 0) {
        std::set<int> matched_heart_vertices;
        for (const auto& point : match_points) {
            matched_heart_vertices.insert(point.heart_id);
        }
        
        int count = 0;
        for (const auto& heart_point : heart_boundary_points) {
            int heart_vertex_id = std::get<0>(heart_point);
            if (matched_heart_vertices.find(heart_vertex_id) == matched_heart_vertices.end() && count < 10) {
                unmatched_heart_points.push_back(heart_point);
                count++;
            }
        }
    }
    
    // 同步所有进程的匹配点数量
    int global_matched_count = 0;
    MPI_Allreduce(&matched_count, &global_matched_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    // 调试输出
    if (debug_output && my_rank == 0) {
        std::cout << "检测到的匹配点详情 (最多显示" << max_debug_points << "个点):" << std::endl;
        std::cout << "心脏顶点ID\t躯干顶点ID\tX坐标\tY坐标\tZ坐标\t距离\t心脏值\t躯干初始值" << std::endl;
        
        for (const auto& point : match_points) {
            std::cout << point.heart_id << "\t" << point.torso_id << "\t"
                      << std::fixed << std::setprecision(6)
                      << point.x << "\t" << point.y << "\t" << point.z << "\t"
                      << point.distance << "\t"
                      << point.heart_val << "\t" << point.torso_val << std::endl;
            
            if (debug_file.is_open()) {
                debug_file << "匹配点," << point.heart_id << "," << point.torso_id << ","
                          << point.x << "," << point.y << "," << point.z << ","
                          << point.distance << ","
                          << point.heart_val << "," << point.torso_val << std::endl;
            }
        }
        
        // 输出未匹配的心脏边界点(可选)
        std::cout << "\n未匹配的心脏边界点样本：" << std::endl;
        for (const auto& point : unmatched_heart_points) {
            int heart_id = std::get<0>(point);
            double x = std::get<1>(point);
            double y = std::get<2>(point);
            double z = std::get<3>(point);
            
            std::cout << "未匹配心脏点: ID=" << heart_id 
                     << ", 坐标=(" << x << ", " << y << ", " << z << ")" 
                     << ", 值=" << heart_boundary_values[heart_id] << std::endl;
            
            if (debug_file.is_open()) {
                debug_file << "未匹配心脏点," << heart_id << ",,"
                          << x << "," << y << "," << z << ",,,"
                          << heart_boundary_values[heart_id] << std::endl;
            }
        }
        
        // 关闭调试文件
        if (debug_file.is_open()) {
            debug_file.close();
        }
    }
    
    if (my_rank == 0) {
        std::cout << "找到并设置了 " << global_matched_count << " 个交界面点的值" << std::endl;
        if (debug_output) {
            std::cout << "详细匹配信息已保存到 boundary_points_matching.csv" << std::endl;
        }
    }
    
    // 同步所有进程上的值
    gf_ue_torso.ParallelAssemble();
    
    if (my_rank == 0) {
        std::cout << "交界面边界条件设置完成。" << std::endl;
    }
}

/**
 * 为FormLinearSystem准备交界面DOF列表
 * 
 * @param gf_ue_torso 躯干解（已设置交界面值）
 * @param pfespace_torso 躯干有限元空间
 * @param heart_mesh 心脏网格
 * @param torso_mesh 躯干网格
 * @param tolerance 坐标匹配容差
 * @return 交界面DOF列表
 */
Array<int> getIntersectionDofs(
    ParGridFunction& gf_ue_torso,
    ParFiniteElementSpace* pfespace_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    double tolerance = 1e-6)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0) {
        std::cout << "准备交界面DOF列表..." << std::endl;
    }
    
    // 收集心脏边界顶点的坐标
    std::set<int> heart_boundary_vertices;
    for (int i = 0; i < heart_mesh->GetNBE(); i++) {
        Array<int> vertices;
        heart_mesh->GetBdrElementVertices(i, vertices);
        for (int j = 0; j < vertices.Size(); j++) {
            heart_boundary_vertices.insert(vertices[j]);
        }
    }
    
    // 创建离散化坐标查找结构
    std::set<std::tuple<int,int,int>> heart_coords;
    for (int vertex_id : heart_boundary_vertices) {
        double* coords = heart_mesh->GetVertex(vertex_id);
        
        int scale = int(1.0 / tolerance);
        int x_key = int(coords[0] * scale);
        int y_key = int(coords[1] * scale);
        int z_key = int(coords[2] * scale);
        
        heart_coords.insert(std::make_tuple(x_key, y_key, z_key));
    }
    
    // 找出躯干网格上与心脏表面重合的点，并获取对应的DOF
    std::set<int> intersection_dofs_set;
    
    // 首先找出躯干网格上与心脏表面重合的顶点
    std::set<int> intersection_vertices;
    for (int i = 0; i < torso_mesh->GetNV(); i++) {
        double* coords = torso_mesh->GetVertex(i);
        
        int scale = int(1.0 / tolerance);
        int x_key = int(coords[0] * scale);
        int y_key = int(coords[1] * scale);
        int z_key = int(coords[2] * scale);
        
        std::tuple<int,int,int> key(x_key, y_key, z_key);
        
        if (heart_coords.find(key) != heart_coords.end()) {
            intersection_vertices.insert(i);
        }
    }
    
    // 为这些顶点获取真DOF
    Array<int> vdofs;
    for (int vertex_id : intersection_vertices) {
        pfespace_torso->GetVertexDofs(vertex_id, vdofs);
        for (int i = 0; i < vdofs.Size(); i++) {
            int vdof = vdofs[i];
            if (vdof >= 0) {  // 本地拥有的DOF
                int tdof = pfespace_torso->GetLocalTDofNumber(vdof);
                if (tdof >= 0) {
                    intersection_dofs_set.insert(tdof);
                }
            }
        }
    }
    
    // 转换为Array
    Array<int> intersection_dofs(intersection_dofs_set.size());
    int count = 0;
    for (int dof : intersection_dofs_set) {
        intersection_dofs[count++] = dof;
    }
    
    // 同步所有进程的DOF数量
    int local_count = intersection_dofs.Size();
    int global_count = 0;
    MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (my_rank == 0) {
        std::cout << "找到 " << global_count << " 个交界面DOF" << std::endl;
    }
    
    return intersection_dofs;
}

/**
 * 完整流程：设置交界面边界条件并获取DOF列表
 * 
 * @param gf_ue_heart 心脏解
 * @param gf_ue_torso 躯干解
 * @param heart_mesh 心脏网格
 * @param torso_mesh 躯干网格
 * @param pfespace_torso 躯干有限元空间
 * @param tolerance 坐标匹配容差
 * @return 交界面DOF列表
 */
Array<int> setupIntersectionBoundary(
    const ParGridFunction& gf_ue_heart,
    ParGridFunction& gf_ue_torso,
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    ParFiniteElementSpace* pfespace_torso,
    double tolerance = 1e-6)
{
    // 设置边界条件
    setIntersectionBoundary(gf_ue_heart, gf_ue_torso, heart_mesh, torso_mesh, tolerance, true, 50);
    
    // 获取DOF列表
    return getIntersectionDofs(gf_ue_torso, pfespace_torso, heart_mesh, torso_mesh, tolerance);
}




/**
 * 使用全局网格方法识别心脏-躯干交界面DOF
 * 只在主进程(rank 0)执行几何交界面检测，然后将结果广播给所有进程
 * 
 * @param heart_mesh 心脏网格(串行)
 * @param torso_mesh 躯干网格(串行)
 * @param pfespace_torso 躯干有限元空间(并行)
 * @param heart_boundary_marker 标记心脏-躯干交界面的属性值
 * @param tolerance 几何匹配容差
 * @return 交界面上的DOF列表
 */
Array<int> findHeartTorsoIntersectionDofs(
    mfem::Mesh* heart_mesh,
    mfem::Mesh* torso_mesh,
    ParFiniteElementSpace* pfespace_torso,
    int heart_boundary_marker = 1,
    double tolerance = 1e-6)
{
    int my_rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    Array<int> boundary_elements;  // 将存储交界面边界元素ID
    
    // 只在主进程执行几何识别
    if (my_rank == 0) {
        if (heart_mesh == nullptr || torso_mesh == nullptr) {
            std::cerr << "错误：空的网格指针!" << std::endl;
            // 返回空数组
            int size = 0;
            MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
            return boundary_elements;
        }
        
        std::cout << "在主进程上执行心脏-躯干交界面识别..." << std::endl;
        
        // 收集心脏表面点
        std::map<int, Vector> heart_surface_points;
        
        for (int i = 0; i < heart_mesh->GetNBE(); i++) {
            Array<int> vertices;
            heart_mesh->GetBdrElementVertices(i, vertices);
            
            for (int j = 0; j < vertices.Size(); j++) {
                int vertex_id = vertices[j];
                double* coords = heart_mesh->GetVertex(vertex_id);
                
                Vector point(3);
                point(0) = coords[0];
                point(1) = coords[1];
                point(2) = coords[2];
                
                heart_surface_points[vertex_id] = point;
            }
        }
        
        std::cout << "收集到 " << heart_surface_points.size() << " 个心脏表面点" << std::endl;
        
        // 在躯干网格中查找接近心脏表面的边界元素
        std::set<int> intersection_elements;
        
        for (int i = 0; i < torso_mesh->GetNBE(); i++) {
            Array<int> vertices;
            torso_mesh->GetBdrElementVertices(i, vertices);
            
            bool has_intersection = false;
            
            for (int j = 0; j < vertices.Size(); j++) {
                int vertex_id = vertices[j];
                double* coords = torso_mesh->GetVertex(vertex_id);
                
                Vector torso_point(3);
                torso_point(0) = coords[0];
                torso_point(1) = coords[1];
                torso_point(2) = coords[2];
                
                // 查找最近的心脏表面点
                double min_distance = tolerance * 2; // 初始化为大于tolerance的值
                
                for (const auto& heart_pair : heart_surface_points) {
                    const Vector& heart_point = heart_pair.second;
                    
                    // 计算距离
                    double distance = 0.0;
                    for (int k = 0; k < 3; k++) {
                        double diff = heart_point(k) - torso_point(k);
                        distance += diff * diff;
                    }
                    distance = sqrt(distance);
                    
                    if (distance < min_distance) {
                        min_distance = distance;
                    }
                    
                    // 如果找到距离足够近的点，标记为交界面并停止搜索
                    if (distance < tolerance) {
                        has_intersection = true;
                        break;
                    }
                }
                
                if (has_intersection) {
                    break;
                }
            }
            
            if (has_intersection) {
                intersection_elements.insert(i);
            }
        }
        
        std::cout << "识别出 " << intersection_elements.size() << " 个交界面边界元素" << std::endl;
        
        // 保存交界面边界元素ID到数组，用于广播
        boundary_elements.SetSize(intersection_elements.size());
        int idx = 0;
        for (int element_id : intersection_elements) {
            boundary_elements[idx++] = element_id;
        }
    }
    
    // 广播交界面边界元素数量
    int num_boundary_elements = boundary_elements.Size();
    MPI_Bcast(&num_boundary_elements, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // 在非主进程上分配空间
    if (my_rank != 0) {
        boundary_elements.SetSize(num_boundary_elements);
    }
    
    // 广播交界面边界元素ID
    if (num_boundary_elements > 0) {
        MPI_Bcast(boundary_elements.GetData(), num_boundary_elements, MPI_INT, 0, MPI_COMM_WORLD);
    }
    
    // 在所有进程上，使用这些边界元素ID来标记交界面并获取DOF
    
    // 获取ParMesh以便访问边界属性
    ParMesh* pmesh_torso = pfespace_torso->GetParMesh();
    
    // 保存原始边界属性
    Array<int> old_attrs(pmesh_torso->GetNBE());
    for (int i = 0; i < pmesh_torso->GetNBE(); i++) {
        old_attrs[i] = pmesh_torso->GetBdrAttribute(i);
    }
    
    // 由于没有GetBdrElementGlobalIndex，我们需要另一种方法来匹配边界元素
    // 一种替代方法是使用边界元素的几何特征（如顶点坐标）来匹配
    
    // 首先，在主进程上收集交界面边界元素的顶点坐标
    std::vector<std::vector<Vector>> interface_element_coords;
    
    if (my_rank == 0) {
        interface_element_coords.resize(boundary_elements.Size());
        for (int i = 0; i < boundary_elements.Size(); i++) {
            int element_id = boundary_elements[i];
            Array<int> vertices;
            torso_mesh->GetBdrElementVertices(element_id, vertices);
            
            std::vector<Vector> coords;
            for (int j = 0; j < vertices.Size(); j++) {
                double* vertex_coords = torso_mesh->GetVertex(vertices[j]);
                Vector point(3);
                point(0) = vertex_coords[0];
                point(1) = vertex_coords[1];
                point(2) = vertex_coords[2];
                coords.push_back(point);
            }
            
            interface_element_coords[i] = coords;
        }
    }
    
    // 广播交界面元素的顶点数量
    std::vector<int> num_vertices_per_element;
    if (my_rank == 0) {
        for (const auto& coords : interface_element_coords) {
            num_vertices_per_element.push_back(coords.size());
        }
    } else {
        num_vertices_per_element.resize(num_boundary_elements);
    }
    
    MPI_Bcast(num_vertices_per_element.data(), num_boundary_elements, MPI_INT, 0, MPI_COMM_WORLD);
    
    // 广播顶点坐标
    std::vector<double> flat_coords;
    if (my_rank == 0) {
        for (const auto& element_coords : interface_element_coords) {
            for (const auto& point : element_coords) {
                for (int j = 0; j < 3; j++) {
                    flat_coords.push_back(point(j));
                }
            }
        }
    }
    
    int total_coords = 0;
    for (int n : num_vertices_per_element) {
        total_coords += n * 3;  // 每个顶点3个坐标
    }
    
    if (my_rank != 0) {
        flat_coords.resize(total_coords);
    }
    
    MPI_Bcast(flat_coords.data(), total_coords, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // 在每个进程上重建交界面元素坐标
    if (my_rank != 0) {
        interface_element_coords.resize(num_boundary_elements);
        int cursor = 0;
        for (int i = 0; i < num_boundary_elements; i++) {
            int num_vertices = num_vertices_per_element[i];
            std::vector<Vector> coords;
            
            for (int j = 0; j < num_vertices; j++) {
                Vector point(3);
                for (int k = 0; k < 3; k++) {
                    point(k) = flat_coords[cursor++];
                }
                coords.push_back(point);
            }
            
            interface_element_coords[i] = coords;
        }
    }
    
    // 在所有进程上查找匹配的边界元素
    for (int i = 0; i < pmesh_torso->GetNBE(); i++) {
        Array<int> vertices;
        pmesh_torso->GetBdrElementVertices(i, vertices);
        
        std::vector<Vector> element_coords;
        for (int j = 0; j < vertices.Size(); j++) {
            double* vertex_coords = pmesh_torso->GetVertex(vertices[j]);
            Vector point(3);
            point(0) = vertex_coords[0];
            point(1) = vertex_coords[1];
            point(2) = vertex_coords[2];
            element_coords.push_back(point);
        }
        
        // 检查当前边界元素是否与任何交界面元素匹配
        bool is_on_interface = false;
        
        for (const auto& interface_coords : interface_element_coords) {
            // 简单匹配：检查每个顶点是否有匹配
            if (interface_coords.size() != element_coords.size()) {
                continue;  // 顶点数不同，不可能匹配
            }
            
            // 检查是否所有顶点都匹配（考虑容差）
            int matched_vertices = 0;
            for (const auto& element_point : element_coords) {
                for (const auto& interface_point : interface_coords) {
                    double distance = 0.0;
                    for (int k = 0; k < 3; k++) {
                        double diff = element_point(k) - interface_point(k);
                        distance += diff * diff;
                    }
                    distance = sqrt(distance);
                    
                    if (distance < tolerance) {
                        matched_vertices++;
                        break;
                    }
                }
            }
            
            if (matched_vertices == element_coords.size()) {
                is_on_interface = true;
                break;
            }
        }
        
        // 设置边界属性
        if (is_on_interface) {
            pmesh_torso->SetBdrAttribute(i, heart_boundary_marker);
        } else {
            pmesh_torso->SetBdrAttribute(i, heart_boundary_marker + 1);
        }
    }
    
    // 设置边界标记数组
    int max_attr = pmesh_torso->bdr_attributes.Max();
    Array<int> ess_bdr(max_attr);
    ess_bdr = 0;
    if (heart_boundary_marker <= max_attr) {
        ess_bdr[heart_boundary_marker - 1] = 1;  // 只标记交界面
    } else {
        if (my_rank == 0) {
            std::cout << "警告: 心脏-躯干交界面标记(" << heart_boundary_marker 
                      << ")超出了最大边界属性值(" << max_attr << ")" << std::endl;
        }
    }
    
    // 获取交界面上的DOF
    Array<int> dofs_on_intersection;
    pfespace_torso->GetEssentialTrueDofs(ess_bdr, dofs_on_intersection);
    
    // 恢复原始边界属性
    for (int i = 0; i < pmesh_torso->GetNBE(); i++) {
        pmesh_torso->SetBdrAttribute(i, old_attrs[i]);
    }
    
    // 输出结果统计
    int local_dofs = dofs_on_intersection.Size();
    int global_dofs = 0;
    MPI_Allreduce(&local_dofs, &global_dofs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (my_rank == 0) {
        std::cout << "找到心脏-躯干交界面上的DOF数量: " << global_dofs << std::endl;
    }
    
    return dofs_on_intersection;
}

/**
 * 设置躯干Dirichlet边界条件函数
 * 基于空间位置映射建立心脏-躯干交界面的对应关系
 * 
 * @param gf_ue_heart 心脏网格上的解
 * @param gf_ue_torso 躯干网格上的解
 * @param ess_tdof_list 交界面上的DOF列表
 * @param heart_mesh 心脏网格
 * @param torso_mesh 躯干网格
 * @param distance_threshold 点匹配的距离阈值
 */
void setTorsoIntersectionBoundaryConditions(
   const ParGridFunction& gf_ue_heart,
   ParGridFunction& gf_ue_torso,
   const Array<int>& ess_tdof_list,
   mfem::Mesh* heart_mesh,
   mfem::Mesh* torso_mesh,
   double distance_threshold = 1e-4)  // 增大阈值，默认值改为1e-4
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (my_rank == 0) {
        std::cout << "在心脏-躯干交界面上设置Dirichlet边界条件..." << std::endl;
        std::cout << "交界面DOF数量: " << ess_tdof_list.Size() << std::endl;
        std::cout << "使用距离阈值: " << distance_threshold << std::endl;
    }
    
    // 获取有限元空间
    ParFiniteElementSpace* pfespace_heart = gf_ue_heart.ParFESpace();
    ParFiniteElementSpace* pfespace_torso = gf_ue_torso.ParFESpace();
    
    // 收集心脏网格的所有边界点及其值
    std::vector<std::tuple<Vector, double>> heart_boundary_points;
    
    // 首先收集心脏边界顶点
    for (int i = 0; i < heart_mesh->GetNBE(); i++) {
        Array<int> vertices;
        heart_mesh->GetBdrElementVertices(i, vertices);
        
        for (int j = 0; j < vertices.Size(); j++) {
            int vertex_id = vertices[j];
            if (vertex_id < gf_ue_heart.Size()) {
                double* coords = heart_mesh->GetVertex(vertex_id);
                
                Vector point(3);
                point(0) = coords[0];
                point(1) = coords[1];
                point(2) = coords[2];
                
                // 存储该点的电位值
                heart_boundary_points.push_back(std::make_tuple(point, gf_ue_heart(vertex_id)));
            }
        }
    }
    
    if (my_rank == 0) {
        std::cout << "收集到心脏边界点数量: " << heart_boundary_points.size() << std::endl;
        
        // 输出一些心脏边界点的值范围
        if (!heart_boundary_points.empty()) {
            double min_val = std::get<1>(heart_boundary_points[0]);
            double max_val = min_val;
            
            for (const auto& point_data : heart_boundary_points) {
                double val = std::get<1>(point_data);
                min_val = std::min(min_val, val);
                max_val = std::max(max_val, val);
            }
            
            std::cout << "心脏边界点值范围: [" << min_val << ", " << max_val << "]" << std::endl;
        }
    }
    
    // 创建调试输出文件
    std::ofstream debug_file;
    if (my_rank == 0) {
        debug_file.open("boundary_mapping_debug.csv");
        if (debug_file.is_open()) {
            debug_file << "躯干顶点ID,躯干坐标X,躯干坐标Y,躯干坐标Z,心脏值,距离,原躯干值,是否在交界面列表中" << std::endl;
        } else {
            std::cerr << "错误：无法创建调试文件" << std::endl;
        }
    }
    
    // 处理躯干网格上的交界面点
    int local_set_count = 0;
    double local_min = std::numeric_limits<double>::max();
    double local_max = -std::numeric_limits<double>::max();
    bool has_values = false;
    
    // 创建一个集合来存储交界面DOF
    std::set<int> interface_tdofs;
    for (int i = 0; i < ess_tdof_list.Size(); i++) {
        interface_tdofs.insert(ess_tdof_list[i]);
    }
    
    // 遍历躯干网格的所有顶点
    for (int vertex_id = 0; vertex_id < torso_mesh->GetNV(); vertex_id++) {
        // 检查该顶点是否在边界上
        bool is_boundary_vertex = false;
        for (int i = 0; i < torso_mesh->GetNBE(); i++) {
            Array<int> elem_vertices;
            torso_mesh->GetBdrElementVertices(i, elem_vertices);
            for (int j = 0; j < elem_vertices.Size(); j++) {
                if (elem_vertices[j] == vertex_id) {
                    is_boundary_vertex = true;
                    break;
                }
            }
            if (is_boundary_vertex) break;
        }
        
        if (!is_boundary_vertex || vertex_id >= gf_ue_torso.Size()) {
            continue;
        }
        
        // 检查是否为交界面点（检查所有相关DOF）
        Array<int> vertex_dofs;
        pfespace_torso->GetVertexDofs(vertex_id, vertex_dofs);
        
        bool is_on_interface = false;
        for (int k = 0; k < vertex_dofs.Size(); k++) {
            int vdof = vertex_dofs[k];
            if (vdof >= 0) {
                int tdof = pfespace_torso->GetLocalTDofNumber(vdof);
                if (tdof >= 0 && interface_tdofs.count(tdof) > 0) {
                    is_on_interface = true;
                    break;
                }
            }
        }
        
        // 获取顶点坐标
        double* coords = torso_mesh->GetVertex(vertex_id);
        Vector torso_point(3);
        torso_point(0) = coords[0];
        torso_point(1) = coords[1];
        torso_point(2) = coords[2];
        
        // 无论是否在交界面列表中，先找最近的心脏点
        double min_distance = std::numeric_limits<double>::max();
        double heart_value = 0.0;
        bool found_match = false;
        
        for (const auto& heart_point_data : heart_boundary_points) {
            const Vector& heart_point = std::get<0>(heart_point_data);
            
            // 计算距离
            double distance = 0.0;
            for (int k = 0; k < 3; k++) {
                double diff = heart_point(k) - torso_point(k);
                distance += diff * diff;
            }
            distance = std::sqrt(distance);
            
            if (distance < min_distance) {
                min_distance = distance;
                heart_value = std::get<1>(heart_point_data);
                
                if (distance < distance_threshold) {
                    found_match = true;
                }
            }
        }
        
        // 写入调试文件
        if (my_rank == 0 && debug_file.is_open()) {
            debug_file << vertex_id << ","
                     << coords[0] << "," << coords[1] << "," << coords[2] << ","
                     << heart_value << "," << min_distance << "," << gf_ue_torso(vertex_id) 
                     << "," << (is_on_interface ? "是" : "否") << std::endl;
        }
        
        // 如果找到匹配的心脏点，设置边界值
        if (found_match && is_on_interface) {
            double original_value = gf_ue_torso(vertex_id);
            gf_ue_torso(vertex_id) = heart_value;
            local_set_count++;
            
            // 更新统计数据
            if (std::isfinite(heart_value)) {
                local_min = std::min(local_min, heart_value);
                local_max = std::max(local_max, heart_value);
                has_values = true;
            }
            
            // 输出前10个点的调试信息
            if (my_rank == 0 && local_set_count <= 10) {
                std::cout << "设置交界面点 " << local_set_count 
                         << ", 躯干顶点=" << vertex_id 
                         << ", 心脏值=" << heart_value 
                         << ", 距离=" << min_distance
                         << ", 躯干原值=" << original_value << std::endl;
            }
        }
    }
    
    // 关闭调试文件
    if (my_rank == 0 && debug_file.is_open()) {
        debug_file.close();
    }
    
    // 同步所有进程上的值
    gf_ue_torso.ParallelAssemble();
    
    // 二次检查设置的值是否被保留
    if (my_rank == 0 && local_set_count > 0) {
        std::cout << "检查值是否正确同步..." << std::endl;
        int check_count = 0;
        for (int vertex_id = 0; vertex_id < torso_mesh->GetNV() && check_count < 5; vertex_id++) {
            if (std::fabs(gf_ue_torso(vertex_id)) > 1e-10) {
                std::cout << "顶点 " << vertex_id << " 值: " << gf_ue_torso(vertex_id) << std::endl;
                check_count++;
            }
        }
    }
    
    // 收集统计信息
    int global_set_count = 0;
    MPI_Allreduce(&local_set_count, &global_set_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    // 避免MPI错误，只在需要时收集最小和最大值
    int has_local_values = has_values ? 1 : 0;
    int has_global_values = 0;
    MPI_Allreduce(&has_local_values, &has_global_values, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    
    double global_min = 0.0, global_max = 0.0;
    if (has_global_values) {
        // 如果本地没有值，设置为不影响全局min/max的值
        if (!has_values) {
            local_min = std::numeric_limits<double>::max();
            local_max = -std::numeric_limits<double>::max();
        }
        
        MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    
    if (my_rank == 0) {
        std::cout << "共设置了 " << global_set_count << " 个边界点的值" << std::endl;
        if (has_global_values) {
            std::cout << "边界条件值范围: [" << global_min << ", " << global_max << "]" << std::endl;
        } else {
            std::cout << "警告: 没有设置任何有效的边界点值!" << std::endl;
        }
        
        // 验证边界条件是否与心脏解范围一致
        double heart_min = gf_ue_heart.Min();
        double heart_max = gf_ue_heart.Max();
        std::cout << "心脏解范围: [" << heart_min << ", " << heart_max << "]" << std::endl;
        
        if (global_min < heart_min || global_max > heart_max) {
            std::cout << "警告: 边界条件范围超出心脏解范围!" << std::endl;
        }
        
        std::cout << "边界条件详细信息已保存到 boundary_mapping_debug.csv" << std::endl;
        std::cout << "Dirichlet边界条件设置完成。" << std::endl;
    }
    
    // 额外的检查：确保交界面DOF的值不会在后续步骤中被重置
    if (my_rank == 0) {
        std::cout << "注意: 请确保在求解过程中不要重置边界条件值!" << std::endl;
    }
}






// 用于比较两个点坐标是否相等的辅助结构体
struct Point3D {
   double x, y, z;
   
   Point3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
   
   // 重载小于运算符以便用于map
   bool operator<(const Point3D& other) const {
       if (x != other.x) return x < other.x;
       if (y != other.y) return y < other.y;
       return z < other.z;
   }
   
   // 点坐标相等性检查（带容差）
   bool equals(const Point3D& other, double tolerance = 1e-4) const {
       return std::abs(x - other.x) < tolerance && 
              std::abs(y - other.y) < tolerance && 
              std::abs(z - other.z) < tolerance;
   }
};

/**
 * 读取torso网格并识别与heart网格的重叠边界
 * 添加更多调试信息和容错性
 * 
 * @param torso_mesh_file torso网格文件路径
 * @param heart_mesh 已经读取的heart网格
 * @param heart_boundary_marker 需要在torso网格上标记的边界属性值
 * @param tolerance 坐标匹配的容差，默认值增加到1e-6
 * @return 返回读取的torso网格
 */
 mfem::Mesh* readTorsoMeshAndIdentifyBoundary(
   const std::string& torso_mesh_file, 
   mfem::Mesh* heart_mesh,
   int heart_boundary_marker = 1,
   double tolerance = 1e-6)  // 增加默认容差
{
   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   
   // 读取torso网格
   if (my_rank == 0) {
       std::cout << "正在读取torso网格: " << torso_mesh_file << std::endl;
   }
   
   mfem::Mesh* torso_mesh = new mfem::Mesh(torso_mesh_file.c_str(), 1, 1);
   
   if (my_rank == 0) {
       std::cout << "Torso网格信息：" << std::endl;
       std::cout << "  维度: " << torso_mesh->Dimension() << std::endl;
       std::cout << "  顶点数: " << torso_mesh->GetNV() << std::endl;
       std::cout << "  元素数: " << torso_mesh->GetNE() << std::endl;
       std::cout << "  边界元素数: " << torso_mesh->GetNBE() << std::endl;
       
       // 检查网格坐标范围
       double xmin = 1e10, ymin = 1e10, zmin = 1e10;
       double xmax = -1e10, ymax = -1e10, zmax = -1e10;
       
       for (int i = 0; i < torso_mesh->GetNV(); i++) {
           double* coords = torso_mesh->GetVertex(i);
           xmin = std::min(xmin, coords[0]);
           ymin = std::min(ymin, coords[1]);
           zmin = std::min(zmin, coords[2]);
           xmax = std::max(xmax, coords[0]);
           ymax = std::max(ymax, coords[1]);
           zmax = std::max(zmax, coords[2]);
       }
       
       std::cout << "  坐标范围: X[" << xmin << ", " << xmax << "], "
                 << "Y[" << ymin << ", " << ymax << "], "
                 << "Z[" << zmin << ", " << zmax << "]" << std::endl;
   }
   
   // 同样检查heart网格的坐标范围
   if (my_rank == 0) {
       double xmin = 1e10, ymin = 1e10, zmin = 1e10;
       double xmax = -1e10, ymax = -1e10, zmax = -1e10;
       
       for (int i = 0; i < heart_mesh->GetNV(); i++) {
           double* coords = heart_mesh->GetVertex(i);
           xmin = std::min(xmin, coords[0]);
           ymin = std::min(ymin, coords[1]);
           zmin = std::min(zmin, coords[2]);
           xmax = std::max(xmax, coords[0]);
           ymax = std::max(ymax, coords[1]);
           zmax = std::max(zmax, coords[2]);
       }
       
       std::cout << "Heart网格坐标范围: X[" << xmin << ", " << xmax << "], "
                 << "Y[" << ymin << ", " << ymax << "], "
                 << "Z[" << zmin << ", " << zmax << "]" << std::endl;
       
       // 检查两个网格是否有重叠
       std::cout << "容差设置: " << tolerance << std::endl;
   }
   
   // 收集heart网格的所有表面点
   std::map<Point3D, int> heart_surface_points;
   
   // 首先收集heart网格的边界元素和点
   for (int i = 0; i < heart_mesh->GetNBE(); i++) {
       Array<int> vertices;
       heart_mesh->GetBdrElementVertices(i, vertices);
       
       for (int j = 0; j < vertices.Size(); j++) {
           double* coords = heart_mesh->GetVertex(vertices[j]);
           Point3D p(coords[0], coords[1], coords[2]);
           heart_surface_points[p] = vertices[j];
       }
   }
   
   if (my_rank == 0) {
       std::cout << "收集到Heart网格表面点数量: " << heart_surface_points.size() << std::endl;
       
       // 输出一些心脏表面点坐标，帮助调试
       int count = 0;
       std::cout << "前5个心脏表面点坐标:" << std::endl;
       for (std::map<Point3D, int>::const_iterator it = heart_surface_points.begin(); 
            it != heart_surface_points.begin() && count < 5; ++it, ++count) {
           std::cout << "  (" << it->first.x << ", " << it->first.y << ", " << it->first.z << ")" << std::endl;
       }
   }
   
   // 遍历torso网格的边界元素，检查哪些与heart网格的表面重叠
   std::set<int> boundary_elements;
   std::set<int> heart_torso_interface_vertices;
   
   // 创建一个计数器来跟踪测试过的点对数量
   int total_tests = 0;
   int close_points = 0;
   
   for (int i = 0; i < torso_mesh->GetNBE(); i++) {
       Array<int> vertices;
       torso_mesh->GetBdrElementVertices(i, vertices);
       
       bool is_on_interface = false;
       
       // 检查该边界元素的每个顶点是否与heart网格的表面点重合
       for (int j = 0; j < vertices.Size(); j++) {
           double* coords = torso_mesh->GetVertex(vertices[j]);
           Point3D p(coords[0], coords[1], coords[2]);
           
           // 在heart表面点中查找相同位置的点
           for (std::map<Point3D, int>::const_iterator it = heart_surface_points.begin(); 
                it != heart_surface_points.end(); ++it) {
               total_tests++;
               
               if (p.equals(it->first, tolerance)) {
                   is_on_interface = true;
                   heart_torso_interface_vertices.insert(vertices[j]);
                   close_points++;
                   break;
               }
           }
           
           if (is_on_interface) {
               break;
           }
       }
       
       if (is_on_interface) {
           boundary_elements.insert(i);
           // 将该边界元素标记为心脏-躯干交界面
           torso_mesh->GetBdrElement(i)->SetAttribute(heart_boundary_marker);
       }
   }
   
   if (my_rank == 0) {
       std::cout << "识别出的Torso网格与Heart交界面的边界元素数量: " << boundary_elements.size() << std::endl;
       std::cout << "识别出的Torso网格与Heart交界面的顶点数量: " << heart_torso_interface_vertices.size() << std::endl;
       std::cout << "测试了 " << total_tests << " 对点，找到 " << close_points << " 对匹配点" << std::endl;
       
       // 如果没有找到匹配点，输出一些诊断信息
       if (heart_torso_interface_vertices.empty()) {
           std::cout << "警告：未找到任何匹配的Heart-Torso交界面点！" << std::endl;
           std::cout << "可能的原因：" << std::endl;
           std::cout << "  1. 坐标系统不同或有偏移" << std::endl;
           std::cout << "  2. 容差设置不足" << std::endl;
           std::cout << "  3. 网格文件中的心脏和躯干没有共同边界" << std::endl;
           
           // 输出一些躯干边界点的坐标
           std::cout << "5个随机Torso边界点坐标:" << std::endl;
           int count = 0;
           for (int i = 0; i < torso_mesh->GetNBE() && count < 5; i++) {
               Array<int> vertices;
               torso_mesh->GetBdrElementVertices(i, vertices);
               if (vertices.Size() > 0) {
                   double* coords = torso_mesh->GetVertex(vertices[0]);
                   std::cout << "  (" << coords[0] << ", " << coords[1] << ", " << coords[2] << ")" << std::endl;
                   count++;
               }
           }
           
           // 计算并打印心脏表面点与躯干边界点之间的最小距离
           double min_distance = 1e10;
           Point3D closest_heart_point(0, 0, 0);
           Point3D closest_torso_point(0, 0, 0);
           
           // 随机抽样一些点来检查，避免检查所有组合
           int sample_size = std::min(100, int(heart_surface_points.size()));
           
           std::cout << "计算心脏-躯干点的最小距离（抽样" << sample_size << "个心脏点）..." << std::endl;
           
           int heart_count = 0;
           for (std::map<Point3D, int>::const_iterator it = heart_surface_points.begin(); 
                it != heart_surface_points.end() && heart_count < sample_size; ++it, ++heart_count) {
               
               const Point3D& heart_point = it->first;
               
               for (int i = 0; i < torso_mesh->GetNBE(); i++) {
                   Array<int> vertices;
                   torso_mesh->GetBdrElementVertices(i, vertices);
                   
                   for (int j = 0; j < vertices.Size(); j++) {
                       double* coords = torso_mesh->GetVertex(vertices[j]);
                       Point3D torso_point(coords[0], coords[1], coords[2]);
                       
                       double dx = heart_point.x - torso_point.x;
                       double dy = heart_point.y - torso_point.y;
                       double dz = heart_point.z - torso_point.z;
                       double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                       
                       if (distance < min_distance) {
                           min_distance = distance;
                           closest_heart_point = heart_point;
                           closest_torso_point = torso_point;
                       }
                   }
               }
           }
           
           std::cout << "心脏-躯干点的最小距离: " << min_distance << std::endl;
           std::cout << "最近心脏点: (" << closest_heart_point.x << ", " 
                     << closest_heart_point.y << ", " << closest_heart_point.z << ")" << std::endl;
           std::cout << "最近躯干点: (" << closest_torso_point.x << ", " 
                     << closest_torso_point.y << ", " << closest_torso_point.z << ")" << std::endl;
           
           // 建议的解决方案
           std::cout << "建议：" << std::endl;
           std::cout << "  1. 尝试增加容差值（当前设置：" << tolerance << "）" << std::endl;
           std::cout << "  2. 检查两个网格的坐标系统是否一致" << std::endl;
           std::cout << "  3. 可能需要对一个或两个网格进行变换以使它们对齐" << std::endl;
           
           // 自动尝试更宽松的容差
           if (min_distance < 1.0) {  // 如果最小距离小于1.0，可能只是容差问题
               double new_tolerance = min_distance * 1.1;  // 稍微大于最小距离
               std::cout << "自动调整容差至: " << new_tolerance << " 并重试..." << std::endl;
               
               // 重新检查使用新容差
               heart_torso_interface_vertices.clear();
               boundary_elements.clear();
               
               for (int i = 0; i < torso_mesh->GetNBE(); i++) {
                   Array<int> vertices;
                   torso_mesh->GetBdrElementVertices(i, vertices);
                   
                   bool is_on_interface = false;
                   
                   for (int j = 0; j < vertices.Size(); j++) {
                       double* coords = torso_mesh->GetVertex(vertices[j]);
                       Point3D p(coords[0], coords[1], coords[2]);
                       
                       for (std::map<Point3D, int>::const_iterator it = heart_surface_points.begin(); 
                            it != heart_surface_points.end(); ++it) {
                           const Point3D& heart_point = it->first;
                           
                           double dx = p.x - heart_point.x;
                           double dy = p.y - heart_point.y;
                           double dz = p.z - heart_point.z;
                           double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                           
                           if (distance <= new_tolerance) {
                               is_on_interface = true;
                               heart_torso_interface_vertices.insert(vertices[j]);
                               break;
                           }
                       }
                       
                       if (is_on_interface) {
                           break;
                       }
                   }
                   
                   if (is_on_interface) {
                       boundary_elements.insert(i);
                       torso_mesh->GetBdrElement(i)->SetAttribute(heart_boundary_marker);
                   }
               }
               
               std::cout << "使用新容差后，识别出 " << boundary_elements.size() << " 个交界面边界元素，" 
                        << heart_torso_interface_vertices.size() << " 个交界面顶点。" << std::endl;
           }
       }
   }
   
   // 更新torso网格以使新的边界属性生效
   torso_mesh->SetAttributes();
   
   return torso_mesh;
}