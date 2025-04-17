#include "mfem.hpp"
#include "object.h"
#include "object_cc.hh"
#include "ddcMalloc.h"
#include "pio.h"
#include "pioFixedRecordHelper.h"
#include "units.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <cassert>
#include <memory>
#include <set>
#include <dirent.h>
#include <regex.h>
#include <unistd.h>
#include <sys/stat.h>
#include "util.hpp"
#include "MatrixElementPiecewiseCoefficient.hpp"
#include "cardiac_coefficients.hpp"
#include "torsoSolver.hpp"

#include <map>
#include <unordered_set>
#include <algorithm>
#include <cmath>

#include <mpi.h> // Include MPI header
#include <iomanip> // For formatted output

#define StartTimer(x)
#define EndTimer()

using namespace mfem;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;


const double DEFAULT_HEART_POTENTIAL = -83.0;  // 默认心脏电位值（仅作为备用）


void SaveUsingParallelDataCollection(ParMesh *pmesh, ParGridFunction *pgf, const string &output_dir, const string &collection_name, int cycle, MPI_Comm comm)
{
    int myrank;
    MPI_Comm_rank(comm, &myrank);
    
    // 创建并行数据收集对象
    ParaViewDataCollection *pvdc = new ParaViewDataCollection(collection_name, pmesh);
    
    // 设置输出路径
    pvdc->SetPrefixPath(output_dir);
    
    // 注册网格函数
    pvdc->RegisterField("solution", pgf);
    
    // 设置当前周期和可选的物理时间
    pvdc->SetCycle(cycle);
    pvdc->SetTime(cycle * 0.01); // 示例时间步长
    
    // 保存数据
    pvdc->Save();
    
    if (myrank == 0)
    {
        cout << "使用ParaViewDataCollection保存了解，路径: " << output_dir << "/" << collection_name << endl;
    }
    
    delete pvdc;
}

// 定义边界条件类型的映射
enum BoundaryType {
    NEUMANN_ZERO = 1,   // 零Neumann边界(外部边界)
    DIRICHLET = 100,    // Dirichlet边界(心脏-躯干界面)
    SOURCE = 100        // 源边界(与Dirichlet边界相同)
};

// 边界点结构
struct BoundaryPoint {
    Vector coords;
    double value;
    
    BoundaryPoint(const Vector& c, double v) : coords(c), value(v) {}
};

// 辅助类：基于坐标的网格函数插值器
class CoordinateBasedTransfer {
private:
    ParMesh* source_mesh;
    ParFiniteElementSpace* source_fes;
    ParGridFunction* source_gf;
    
    // 用于空间查找的容忍度
    const double spatial_tol = 1e-8;
    
    // 存储源网格边界点的坐标和值
    std::vector<BoundaryPoint> boundary_points;

public:
    CoordinateBasedTransfer(ParMesh* src_mesh, ParFiniteElementSpace* src_fes, 
                           ParGridFunction* src_gf)
        : source_mesh(src_mesh), source_fes(src_fes), source_gf(src_gf) {
        // 构建源网格(心脏)边界点的映射
        BuildBoundaryPointsMap();
    }

    // 构建边界点坐标-值映射
    void BuildBoundaryPointsMap() {
        // 获取源网格的所有边界顶点
        for (int i = 0; i < source_mesh->GetNBE(); i++) {
            Array<int> vertices;
            source_mesh->GetBdrElementVertices(i, vertices);
            
            // 对于每个顶点
            for (int j = 0; j < vertices.Size(); j++) {
                int vdof = vertices[j];
                
                // 获取顶点坐标
                Vector coords(3);
                source_mesh->GetNode(vdof, coords);
                
                // 获取该顶点处的值
                // 注意：我们使用vdof作为索引获取值（与顶点一一对应）
                double value = (*source_gf)(vdof);
                
                // 存储到向量中
                boundary_points.emplace_back(coords, value);
            }
        }
    }

    // 获取目标点的值
    double GetValueAtPoint(const Vector& point) {
        // 查找最近的点
        double min_dist = std::numeric_limits<double>::max();
        double closest_value = DEFAULT_HEART_POTENTIAL;
        
        for (const auto& bp : boundary_points) {
            double dist = 0.0;
            for (int i = 0; i < 3; i++) {
                dist += (bp.coords(i) - point(i)) * (bp.coords(i) - point(i));
            }
            dist = sqrt(dist);
            
            if (dist < min_dist) {
                min_dist = dist;
                closest_value = bp.value;
                
                // 如果距离小于容差，立即返回
                if (dist < spatial_tol) {
                    return closest_value;
                }
            }
        }
        
        // 返回最近点的值
        return closest_value;
    }
};

// 用于设置边界值的系数函数
class BoundaryValuesCoefficient : public Coefficient {
private:
    CoordinateBasedTransfer* transfer;
    
public:
    BoundaryValuesCoefficient(CoordinateBasedTransfer* t) : transfer(t) {}
    
    virtual double Eval(ElementTransformation& T, const IntegrationPoint& ip) {
        Vector physical_coord(3);
        T.Transform(ip, physical_coord);
        return transfer->GetValueAtPoint(physical_coord);
    }
};


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
 * 通过几何比较找出心脏和躯干相交的边界元素
 * 将这些元素标记为特定属性，然后只在这些元素上设置Dirichlet边界条件
 * 
 * @param heart_mesh 心脏网格
 * @param torso_mesh 躯干网格
 * @param torso_space 躯干有限元空间
 * @param heart_boundary_marker 标记心脏-躯干交界面的属性值
 * @param tolerance 几何比较的容差
 * @return 相交边界上的DOF列表
 */
 Array<int> findHeartTorsoIntersectionDofs(
   mfem::Mesh* heart_mesh,
   mfem::Mesh* torso_mesh,
   ParFiniteElementSpace* torso_space,
   int heart_boundary_marker = 1,
   double tolerance = 1e-6)
{
   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   
   if (my_rank == 0) {
       std::cout << "查找心脏-躯干几何相交区域..." << std::endl;
   }
   
   // 收集心脏边界点
   std::set<std::tuple<double, double, double>> heart_boundary_points;
   for (int i = 0; i < heart_mesh->GetNBE(); i++) {
       Array<int> vertices;
       heart_mesh->GetBdrElementVertices(i, vertices);
       
       for (int j = 0; j < vertices.Size(); j++) {
           double* coords = heart_mesh->GetVertex(vertices[j]);
           heart_boundary_points.insert(std::make_tuple(coords[0], coords[1], coords[2]));
       }
   }
   
   if (my_rank == 0) {
       std::cout << "心脏边界点数量: " << heart_boundary_points.size() << std::endl;
   }
   
   // 收集躯干边界点并检查与心脏的相交
   std::set<int> intersection_boundary_elements;
   
   ParMesh* pmesh_torso = torso_space->GetParMesh();
   
   for (int i = 0; i < pmesh_torso->GetNBE(); i++) {
       Array<int> vertices;
       pmesh_torso->GetBdrElementVertices(i, vertices);
       
       bool has_intersection = false;
       for (int j = 0; j < vertices.Size(); j++) {
           double* coords = pmesh_torso->GetVertex(vertices[j]);
           
           // 检查这个点是否接近任何心脏边界点
           for (const auto& heart_point : heart_boundary_points) {
               double heart_x = std::get<0>(heart_point);
               double heart_y = std::get<1>(heart_point);
               double heart_z = std::get<2>(heart_point);
               
               double dx = coords[0] - heart_x;
               double dy = coords[1] - heart_y;
               double dz = coords[2] - heart_z;
               double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
               
               if (dist < tolerance) {
                   has_intersection = true;
                   break;
               }
           }
           
           if (has_intersection) {
               break;
           }
       }
       
       if (has_intersection) {
           intersection_boundary_elements.insert(i);
       }
   }
   
   // 收集所有进程的结果
   int local_count = intersection_boundary_elements.size();
   int global_count = 0;
   MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   
   if (my_rank == 0) {
       std::cout << "找到心脏-躯干相交的边界元素数量: " << global_count << std::endl;
   }
   
   // 收集相交边界上的DOF
   Array<int> dofs_on_intersection;
   
   // 创建临时的局部属性数组，保存当前的边界属性
   Array<int> old_attrs(pmesh_torso->GetNBE());
   for (int i = 0; i < pmesh_torso->GetNBE(); i++) {
       old_attrs[i] = pmesh_torso->GetBdrAttribute(i);
   }
   
   // 临时修改边界属性 - 逐个设置
   for (int i = 0; i < pmesh_torso->GetNBE(); i++) {
       if (intersection_boundary_elements.find(i) != intersection_boundary_elements.end()) {
           // 设置为心脏-躯干交界面标记
           pmesh_torso->SetBdrAttribute(i, heart_boundary_marker);
       } else {
           // 设置为不同的标记
           pmesh_torso->SetBdrAttribute(i, heart_boundary_marker + 1);
       }
   }
   
   // 创建边界标记数组
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
   torso_space->GetEssentialTrueDofs(ess_bdr, dofs_on_intersection);
   
   // 恢复原始边界属性 - 逐个恢复
   for (int i = 0; i < pmesh_torso->GetNBE(); i++) {
       pmesh_torso->SetBdrAttribute(i, old_attrs[i]);
   }
   
   if (my_rank == 0) {
       std::cout << "找到心脏-躯干交界面上的DOF数量: " << dofs_on_intersection.Size() << std::endl;
   }
   
   return dofs_on_intersection;
}

/**
* 设置躯干Dirichlet边界条件函数
* 只在心脏-躯干相交区域上设置
* 
* @param gf_ue_heart 心脏网格上的解
* @param gf_ue_torso 躯干网格上的解
* @param ess_tdof_list 交界面上的DOF列表
*/
void setTorsoIntersectionBoundaryConditions(
   const ParGridFunction& gf_ue_heart,
   ParGridFunction& gf_ue_torso,
   const Array<int>& ess_tdof_list)
{
   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   
   if (my_rank == 0) {
       std::cout << "在心脏-躯干交界面上设置Dirichlet边界条件..." << std::endl;
       std::cout << "交界面DOF数量: " << ess_tdof_list.Size() << std::endl;
   }
   
   // 计算心脏解的值
   double heart_value = 0.0;
   if (gf_ue_heart.Size() > 0) {
       heart_value = gf_ue_heart(0);  // 使用第一个值
   }
   
   // 同步所有进程的心脏值
   double global_heart_value;
   MPI_Allreduce(&heart_value, &global_heart_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   
   if (my_rank == 0) {
       std::cout << "使用心脏值: " << global_heart_value << std::endl;
   }
   
   // 设置交界面上的DOF值
   ParFiniteElementSpace* pfespace_torso = gf_ue_torso.ParFESpace();
   
   for (int i = 0; i < ess_tdof_list.Size(); i++) {
       int true_dof = ess_tdof_list[i];
       int local_dof = pfespace_torso->GetLocalTDofNumber(true_dof);
       if (local_dof >= 0) {  // 本地进程拥有该DOF
           gf_ue_torso(local_dof) = global_heart_value;
       }
   }
   
   // 同步所有进程上的值
   gf_ue_torso.ParallelAssemble();
   
   if (my_rank == 0) {
       std::cout << "Dirichlet边界条件设置完成。" << std::endl;
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






/**
 * 检查电导率张量是否已正确设置
 * 
 * @param sigma 要检查的电导率张量
 * @throws 如果电导率张量为空则抛出异常
 */
void checkConductivityTensors(MatrixElementPiecewiseCoefficient& sigma) {
    // 检查 heartConductivities_ 中是否有条目
    if (sigma.heartConductivities_.empty()) {
        throw std::runtime_error("错误：电导率张量为空!");
    }
}

int solvePseudoBidomainForUe(
    ParMesh* pmesh,
    ParFiniteElementSpace* pfespace,
    ParGridFunction& V_m,
    ParGridFunction& u_e,
    const std::vector<double>& sigma_i_values,
    const std::vector<double>& sigma_e_values,
    std::shared_ptr<ParGridFunction>& fiber_quat,
    std::shared_ptr<ParGridFunction>& sheet_quat,
    std::shared_ptr<ParGridFunction>& transverse_quat,
    Array<int>& ess_tdof_list,
    const std::vector<int>& heartRegions,
    bool enforce_zero_mean = true,
    int print_level = 2,
    double* t_ksp2_total = nullptr)
{
    int my_rank = 0, num_procs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


    
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
    
    if (sigma_i_values.size() < 3 || sigma_e_values.size() < 3) {
        if (my_rank == 0) {
            std::cerr << "错误：电导率数组大小不足。需要至少3个值。" << std::endl;
        }
        return -1;
    }
    
    if (heartRegions.size() * 3 != sigma_i_values.size() || 
        heartRegions.size() * 3 != sigma_e_values.size()) {
        if (my_rank == 0) {
            std::cerr << "错误：电导率数组大小与心脏区域数量不匹配。" << std::endl;
            std::cerr << "心脏区域数量: " << heartRegions.size() << std::endl;
            std::cerr << "细胞内电导率值数量: " << sigma_i_values.size() << std::endl;
            std::cerr << "细胞外电导率值数量: " << sigma_e_values.size() << std::endl;
        }
        return -1;
    }
    
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
    
    // 检查张量是否正确设置
    bool tensors_valid = !sigma_i.heartConductivities_.empty() && 
                          !sigma_sum.heartConductivities_.empty();
    int global_tensors_valid = tensors_valid ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &global_tensors_valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    
    if (global_tensors_valid == 0) {
        if (my_rank == 0) {
            std::cerr << "错误：电导率张量未正确初始化。" << std::endl;
        }
        return -1;
    }
    
    double tensor_end_time = MPI_Wtime();
    if (my_rank == 0 && print_level > 0) {
        std::cout << "电导率张量设置完成，用时: " << (tensor_end_time - tensor_start_time) << " 秒" << std::endl;
        
        // 打印一些电导率值示例
        if (!sigma_i.heartConductivities_.empty()) {
            int example_region = heartRegions[0];
            std::cout << "心脏区域 " << example_region << " 电导率示例:" << std::endl;
            std::cout << "  细胞内电导率: [" 
                     << sigma_i.heartConductivities_[example_region][0] << ", "
                     << sigma_i.heartConductivities_[example_region][1] << ", "
                     << sigma_i.heartConductivities_[example_region][2] << "]" << std::endl;
            std::cout << "  总电导率: [" 
                     << sigma_sum.heartConductivities_[example_region][0] << ", "
                     << sigma_sum.heartConductivities_[example_region][1] << ", "
                     << sigma_sum.heartConductivities_[example_region][2] << "]" << std::endl;
        }
    }
    
    // 3. 建立线性系统
    if (my_rank == 0) std::cout << "[3/8] 构建线性系统..." << std::endl;
    
    double system_start_time = MPI_Wtime();
    
    // 设置左侧矩阵: -∇·((σ_i + σ_e)∇u_e)
    ParBilinearForm *a = new ParBilinearForm(pfespace);
    a->AddDomainIntegrator(new DiffusionIntegrator(sigma_sum));
    
    if (my_rank == 0 && print_level > 1) std::cout << "  组装系统矩阵..." << std::endl;
    a->Assemble();
    
    HypreParMatrix A;
    a->FormSystemMatrix(ess_tdof_list, A);
    
    // 设置右侧向量: ∇·(σ_i∇V_m)
    if (my_rank == 0 && print_level > 1) std::cout << "  计算右侧向量..." << std::endl;
    
    ParBilinearForm temp_form(pfespace);
    temp_form.AddDomainIntegrator(new DiffusionIntegrator(sigma_i));
    temp_form.Assemble();
    
    Vector vm_true(pfespace->GetTrueVSize());
    V_m.GetTrueDofs(vm_true);
    
    HypreParMatrix A_temp;
    temp_form.FormSystemMatrix(ess_tdof_list, A_temp);
    
    Vector rhs(pfespace->GetTrueVSize());
    rhs = 0.0;
    A_temp.Mult(-1.0, vm_true, 0.0, rhs);
    
    double system_end_time = MPI_Wtime();
    if (my_rank == 0 && print_level > 0) {
        std::cout << "线性系统构建完成，用时: " << (system_end_time - system_start_time) << " 秒" << std::endl;
        std::cout << "矩阵大小: " << A.Height() << " x " << A.Width() << std::endl;
        std::cout << "右侧向量大小: " << rhs.Size() << std::endl;
        std::cout << "右侧向量范数: " << rhs.Norml2() << std::endl;
    }
    
    // 4. 检查相容性条件
    //if (my_rank == 0) std::cout << "[4/8] 检查相容性条件..." << std::endl;
    
    double compat_start_time = MPI_Wtime();

    #if 0
    if (enforce_zero_mean && ess_tdof_list.Size() == 0) {
        if (my_rank == 0 && print_level > 0) {
            std::cout << "纯Neumann问题检测到，执行相容性条件检查..." << std::endl;
        }
        
        // 创建常数函数用于强制平均值为零
        ParGridFunction ones(pfespace);
        ones = 1.0;
        
        Vector ones_true(pfespace->GetTrueVSize());
        ones.GetTrueDofs(ones_true);
        
        // 计算RHS积分值
        double local_rhs_sum = 0.0;
        for (int i = 0; i < rhs.Size(); i++) {
            local_rhs_sum += rhs(i);
        }
        
        double global_rhs_sum = 0.0;
        MPI_Allreduce(&local_rhs_sum, &global_rhs_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // 计算域的大小(简化估计)
        double local_count = rhs.Size();
        double global_count = 0.0;
        MPI_Allreduce(&local_count, &global_count, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        if (my_rank == 0 && print_level > 0) {
            std::cout << "右侧项积分值: " << global_rhs_sum << std::endl;
            std::cout << "自由度总数: " << global_count << std::endl;
        }
        
        //if (fabs(global_rhs_sum) > 1e-10) {
        if (0){
            if (my_rank == 0 && print_level > 0) {
                std::cout << "相容性条件不满足，调整右侧项..." << std::endl;
            }
            
            // 调整RHS以满足相容性条件
            double correction = global_rhs_sum / global_count;
            for (int i = 0; i < rhs.Size(); i++) {
                rhs(i) -= correction;
            }
            
            if (my_rank == 0 && print_level > 0) {
                std::cout << "已调整右侧向量，每个元素减去: " << correction << std::endl;
                
                // 重新计算积分，验证调整后相容性条件是否满足
                double new_local_sum = 0.0;
                for (int i = 0; i < rhs.Size(); i++) {
                    new_local_sum += rhs(i);
                }
                
                double new_global_sum = 0.0;
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Allreduce(&new_local_sum, &new_global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
                std::cout << "调整后右侧项积分值: " << new_global_sum << std::endl;
            }
        } else {
            if (my_rank == 0 && print_level > 0) {
                //std::cout << "相容性条件已满足，无需调整。" << std::endl;
                std::cout << "不调整rhs。" << std::endl;
            }
        }
    } else {
        if (my_rank == 0 && print_level > 0) {
            if (ess_tdof_list.Size() > 0) {
                std::cout << "存在Dirichlet边界条件，跳过相容性检查。" << std::endl;
                std::cout << "Dirichlet边界点数量: " << ess_tdof_list.Size() << std::endl;
            } else {
                std::cout << "未启用均值为零约束，跳过相容性检查。" << std::endl;
            }
        }
    }
    #endif
    if (enforce_zero_mean && ess_tdof_list.Size() == 0) {
        // 只有rank 0输出
        if (my_rank == 0 && print_level > 0) {
            std::cout << "纯Neumann问题检测到，执行相容性条件检查..." << std::endl;
        }
        
        // 所有进程执行的代码
        ParGridFunction ones(pfespace);
        ones = 1.0;
        
        Vector ones_true(pfespace->GetTrueVSize());
        ones.GetTrueDofs(ones_true);
        
        double local_rhs_sum = 0.0;
        for (int i = 0; i < rhs.Size(); i++) {
            local_rhs_sum += rhs(i);
        }
        
        double global_rhs_sum = 0.0;
        MPI_Allreduce(&local_rhs_sum, &global_rhs_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        double local_count = rhs.Size();
        double global_count = 0.0;
        MPI_Allreduce(&local_count, &global_count, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // 只有rank 0输出
        if (my_rank == 0 && print_level > 0) {
            std::cout << "右侧项积分值: " << global_rhs_sum << std::endl;
            std::cout << "自由度总数: " << global_count << std::endl;
        }
        
        //if (fabs(global_rhs_sum) > 1e-10) {
        if (0){
            // 只有rank 0输出
            if (my_rank == 0 && print_level > 0) {
                std::cout << "相容性条件不满足，调整右侧项..." << std::endl;
            }
            
            double correction = global_rhs_sum / global_count;
            for (int i = 0; i < rhs.Size(); i++) {
                rhs(i) -= correction;
            }
            
            // 所有进程重新计算局部和
            double new_local_sum = 0.0;
            for (int i = 0; i < rhs.Size(); i++) {
                new_local_sum += rhs(i);
            }
            
            // 所有进程都参与的集体通信
            double new_global_sum = 0.0;
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Allreduce(&new_local_sum, &new_global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            
            // 只有rank 0输出结果
            if (my_rank == 0 && print_level > 0) {
                std::cout << "已调整右侧向量，每个元素减去: " << correction << std::endl;
                std::cout << "调整后右侧项积分值: " << new_global_sum << std::endl;
            }
        } else {
            // 只有rank 0输出
            if (my_rank == 0 && print_level > 0) {
                //std::cout << "相容性条件已满足，无需调整。" << std::endl;
                std::cout << "不调整rhs。" << std::endl;
            }
        }
    } else {
        // 只有rank 0输出
        if (my_rank == 0 && print_level > 0) {
            if (ess_tdof_list.Size() > 0) {
                std::cout << "存在Dirichlet边界条件，跳过相容性检查。" << std::endl;
                std::cout << "Dirichlet边界点数量: " << ess_tdof_list.Size() << std::endl;
            } else {
                std::cout << "未启用均值为零约束，跳过相容性检查。" << std::endl;
            }
        }
    }




    
    double compat_end_time = MPI_Wtime();
    if (my_rank == 0 && print_level > 0) {
        std::cout << "相容性检查完成，用时: " << (compat_end_time - compat_start_time) << " 秒" << std::endl;
    }
    
    // 5. 输出右侧向量详细信息
    if (my_rank == 0) std::cout << "[5/8] 分析右侧向量..." << std::endl;
    
    double analyze_start_time = MPI_Wtime();
    
    if (my_rank == 0 && print_level > 1) {
        // 计算右侧向量统计信息
        int zero_count = 0;
        double max_val = -1e30, min_val = 1e30, sum = 0.0, sum_sq = 0.0;
        
        for (int i = 0; i < rhs.Size(); i++) {
            double val = rhs(i);
            if (fabs(val) < 1e-10) zero_count++;
            max_val = std::max(max_val, val);
            min_val = std::min(min_val, val);
            sum += val;
            sum_sq += val * val;
        }
        
        double mean = sum / rhs.Size();
        double stdev = sqrt(sum_sq/rhs.Size() - mean*mean);
        
        std::cout << "右侧向量统计:" << std::endl;
        std::cout << "  大小: " << rhs.Size() << std::endl;
        std::cout << "  范围: [" << min_val << ", " << max_val << "]" << std::endl;
        std::cout << "  均值: " << mean << std::endl;
        std::cout << "  标准差: " << stdev << std::endl;
        std::cout << "  接近零的值比例: " << (100.0 * zero_count / rhs.Size()) << "%" << std::endl;
        
        if (print_level > 2 && rhs.Size() > 0) {
            // 打印前几个值作为示例
            std::cout << "右侧向量前10个值: ";
            for (int i = 0; i < std::min(10, rhs.Size()); i++) {
                std::cout << rhs(i) << " ";
            }
            std::cout << std::endl;
        }
    }
    
    double analyze_end_time = MPI_Wtime();
    if (my_rank == 0 && print_level > 0) {
        std::cout << "右侧向量分析完成，用时: " << (analyze_end_time - analyze_start_time) << " 秒" << std::endl;
    }
    
    // 6. 设置求解器和预处理器
    if (my_rank == 0) std::cout << "[6/8] 配置求解器..." << std::endl;
    
    double solver_setup_time = MPI_Wtime();
    
    // 设置求解器
    HyprePCG pcg(A);
    pcg.SetTol(1e-8);       // 容差
    pcg.SetAbsTol(1e-8);    // 绝对容差
    pcg.SetMaxIter(1000);    // 最大迭代次数
    
    // 设置求解器打印级别
    // 0 = 无输出, 1 = 仅最终结果, 2 = 每次迭代, 3 = 详细信息
    int pcg_print_level = 0;
    if (print_level > 1) pcg_print_level = 2;  // 显示迭代过程
    if (print_level > 2) pcg_print_level = 3;  // 显示详细信息
    pcg.SetPrintLevel(3);
    
    // 配置预处理器
    if (my_rank == 0 && print_level > 0) {
        std::cout << "配置BoomerAMG预处理器..." << std::endl;
    }
    
    HypreBoomerAMG amg(A);
    amg.SetPrintLevel(print_level > 2 ? 1 : 0);  // AMG打印级别
    
    // 设置AMG参数 - 使用正确的方法名称
    amg.SetCoarsening(10);      // HMIS粗化
    amg.SetMaxLevels(25);       // 最大级数
    amg.SetRelaxType(3);        // 混合Gauss-Seidel
// SetCycleNumSweeps需要两个参数：pre-relaxation和post-relaxation扫描次数
int prerelax_sweeps = 1;
int postrelax_sweeps = 1;
amg.SetCycleNumSweeps(prerelax_sweeps, postrelax_sweeps);
    
    pcg.SetPreconditioner(amg);
    
    if (my_rank == 0 && print_level > 0) {
        std::cout << "求解器配置:" << std::endl;
        std::cout << "  类型: HyprePCG + BoomerAMG" << std::endl;
        std::cout << "  相对容差: 1e-12" << std::endl;
        std::cout << "  绝对容差: 1e-14" << std::endl;
        std::cout << "  最大迭代次数: 2000" << std::endl;
        std::cout << "  打印级别: " << pcg_print_level << std::endl;
    }
    
    double solver_setup_end_time = MPI_Wtime();
    if (my_rank == 0 && print_level > 0) {
        std::cout << "求解器设置完成，用时: " << (solver_setup_end_time - solver_setup_time) << " 秒" << std::endl;
    }
    
    // 7. 求解线性系统
    if (my_rank == 0) std::cout << "[7/8] 求解线性系统..." << std::endl;
    
    double solve_start_time = MPI_Wtime();
    
    Vector X(rhs.Size());
    X = 0.0;  // 初始猜测为零
    
    // 执行求解
 double t_ksp2_start = MPI_Wtime();
pcg.Mult(rhs, X);
double t_ksp2_end = MPI_Wtime();
if (t_ksp2_total != nullptr) {
    *t_ksp2_total += (t_ksp2_end - t_ksp2_start);
}
    
    double solve_end_time = MPI_Wtime();
    
    // 获取迭代信息
    int num_iterations = 0;
    pcg.GetNumIterations(num_iterations);
    
    // 注意：MFEM没有GetFinalNorm方法，我们可以自己计算残差
    double final_res_norm = 0.0;
    // 计算 ||A*X - rhs||
    Vector res(rhs.Size());
    A.Mult(X, res);  // res = A*X
    res -= rhs;      // res = A*X - rhs
    final_res_norm = res.Norml2();
    
    if (my_rank == 0) {
        std::cout << "求解完成，用时: " << (solve_end_time - solve_start_time) << " 秒" << std::endl;
        std::cout << "迭代次数: " << num_iterations << std::endl;
        std::cout << "最终残差范数: " << final_res_norm << std::endl;
        std::cout << "相对残差范数: " << (final_res_norm / rhs.Norml2()) << std::endl;
    }
    
    // 8. 后处理解向量
    if (my_rank == 0) std::cout << "[8/8] 后处理解向量..." << std::endl;
    
    double postproc_start_time = MPI_Wtime();
    
    // 将解向量设置到网格函数
    u_e.SetFromTrueDofs(X);
    
    // 如果需要，强制解的均值为零
    if (enforce_zero_mean && ess_tdof_list.Size() == 0) {
        double local_sum = 0.0;
        for (int i = 0; i < X.Size(); i++) {
            local_sum += X(i);
        }
        
        double global_sum = 0.0;
        MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        double global_count = 0.0;
        local_sum = X.Size();
        MPI_Allreduce(&local_sum, &global_count, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        double mean_value = global_sum / global_count;
        
        if (fabs(mean_value) > 1e-6) {
            if (my_rank == 0 && print_level > 0) {
                std::cout << "调整解向量以确保均值为零，当前均值: " << mean_value << std::endl;
            }
            
            // 从解中减去平均值
            for (int i = 0; i < X.Size(); i++) {
                X(i) -= mean_value;
            }
            
            // 更新网格函数
            u_e.SetFromTrueDofs(X);
            
            if (my_rank == 0 && print_level > 0) {
                std::cout << "已从解中减去均值: " << mean_value << std::endl;
            }
        } else {
            if (my_rank == 0 && print_level > 0) {
                std::cout << "解向量均值已经接近零: " << mean_value << "，无需调整。" << std::endl;
            }
        }
    }
    
    // 分析解向量
    if (my_rank == 0 && print_level > 0) {
        double min_val = X.Min();
        double max_val = X.Max();
        
        // 计算解向量统计信息
        int zero_count = 0, large_count = 0;
        double sum = 0.0, sum_sq = 0.0;
        double threshold = 1000.0;
        
        for (int i = 0; i < X.Size(); i++) {
            double val = X(i);
            if (fabs(val) < 1e-10) zero_count++;
            if (fabs(val) > threshold) large_count++;
            sum += val;
            sum_sq += val * val;
        }
        
        double mean = sum / X.Size();
        double stdev = sqrt(sum_sq/X.Size() - mean*mean);
        
        std::cout << "解向量统计:" << std::endl;
        std::cout << "  范围: [" << min_val << ", " << max_val << "]" << std::endl;
        std::cout << "  均值: " << mean << std::endl;
        std::cout << "  标准差: " << stdev << std::endl;
        std::cout << "  范数: " << X.Norml2() << std::endl;
        std::cout << "  接近零的值比例: " << (100.0 * zero_count / X.Size()) << "%" << std::endl;
        
        if (large_count > 0) {
            std::cout << "警告: 解向量中有 " << large_count 
                      << " 个绝对值大于 " << threshold << " 的值 ("
                      << (100.0 * large_count / X.Size()) << "%)" << std::endl;
        }
        
        if (print_level > 2 && X.Size() > 0) {
            // 打印前几个值作为示例
            std::cout << "解向量前10个值: ";
            for (int i = 0; i < std::min(10, X.Size()); i++) {
                std::cout << X(i) << " ";
            }
            std::cout << std::endl;
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
        std::cout << "  矩阵大小: " << A.Height() << " x " << A.Width() << std::endl;
        std::cout << "  迭代次数: " << num_iterations << std::endl;
        std::cout << "  最终残差: " << final_res_norm << std::endl;
        std::cout << "  解范围: [" << X.Min() << ", " << X.Max() << "]" << std::endl;
        std::cout << "==============================================" << std::endl;
    }
    
    return num_iterations;
}

//Stolen from SingleCell
class Timeline
{
 public:
   Timeline(double dt, double duration)
   {
      maxTimesteps_ = round(duration/dt);
      dt_ = duration/maxTimesteps_;
   }
   int maxTimesteps() const { return maxTimesteps_; };
   double dt() const { return dt_; }
   double maxTime() const { return dt_*maxTimesteps_; }
   double realTimeFromTimestep(int timestep) const
   {
      return timestep*dt_;
   }
   int timestepFromRealTime(double realTime) const
   {
      return round(realTime/dt_);
   }
   std::string outputIdFromTimestep(const int timestep) const
   {
      double resolution = 1e-3;
      int width = 8;
      while (resolution > dt_) {
         resolution /= 10;
         width++;
      }
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(width)
         << int(round(dt_*timestep/resolution));
      return ss.str();
   }

 private:
   double dt_;
   int maxTimesteps_;
};

class OutputCoordinator
{

 private:
   
};

void recursive_mkdir(const std::string dirname, mode_t mode=S_IRWXU|S_IRWXG)
{
   int startSearch=0;
   do
   {
      int endSearch = dirname.find("/", startSearch);
      //if directory doesn't exist
      if (endSearch < 0) {
         endSearch = dirname.length();
      }
      std::string thisDirname = dirname.substr(0, endSearch);
      DIR* dir = opendir(thisDirname.c_str());
      if (dir)
      {
         closedir(dir);
      }
      else if (ENOENT == errno) {
         //make the directory
         int ret = mkdir(thisDirname.c_str(), mode);
         assert(ret == 0);
      }
      startSearch=endSearch+1;
   } while (startSearch < dirname.length());
}

void save1dNumpyArray(const std::string filename, const std::vector<double>& data)
{
   const int NUMPY_HEADER_SIZE = 128;
   char numpyHeaderFull[NUMPY_HEADER_SIZE];
   //FIXME, report the correct endianness.  I'm too lazy ATM.
   char numpyHeader1[] = 
      "\x93NUMPY\x01\x00\x76\x00{'descr': '<f8', 'fortran_order': False, 'shape': (";
   char numpyHeader2[] = ",)}";
   int cursor=0;
   memcpy(numpyHeaderFull+cursor, numpyHeader1, sizeof(numpyHeader1)-1);
   cursor += sizeof(numpyHeader1)-1;
   {
      std::stringstream ss;
      ss << data.size();
      ss.str();
      memcpy(numpyHeaderFull+cursor,ss.str().c_str(),ss.str().size());
      cursor += ss.str().size();
   }
   memcpy(numpyHeaderFull+cursor,numpyHeader2, sizeof(numpyHeader2)-1);
   cursor += sizeof(numpyHeader2)-1;
   for (; cursor<NUMPY_HEADER_SIZE-1; cursor++) {
      numpyHeaderFull[cursor] = ' ';
   }
   numpyHeaderFull[NUMPY_HEADER_SIZE-1] = '\n';
   FILE* numpyFile = fopen(filename.c_str(), "w");
   assert(numpyFile != NULL);
   std::size_t bytesWritten = fwrite(numpyHeaderFull, sizeof(char), NUMPY_HEADER_SIZE, numpyFile);
   assert(bytesWritten == NUMPY_HEADER_SIZE);
   std::size_t doublesWritten = fwrite(&data[0], sizeof(double), data.size(), numpyFile);
   assert(doublesWritten == data.size());
   fclose(numpyFile);
}

int main(int argc, char *argv[])
{
   MPI_Init(NULL,NULL);
   int num_ranks, my_rank;
   MPI_Comm_size(COMM_LOCAL,&num_ranks);
   MPI_Comm_rank(COMM_LOCAL,&my_rank);

   units_internal(1e-3, 1e-9, 1e-3, 1e-3, 1, 1e-9, 1);
   units_external(1e-3, 1e-9, 1e-3, 1e-3, 1, 1e-9, 1);

       // --- Timer Variable Declarations ---
    double t_start, t_end; // Temporary start/end times
    double t_total_elapsed = 0.0;
    double t_total = 0.0;
    double t_problem1_total = 0.0; // e.g., Heart Solve
    double t_problem2_total = 0.0; // e.g., BC Transfer
    double t_problem3_total = 0.0; // e.g., Torso Solve
    double t_ksp1_total = 0.0;     // e.g., Heart KSP
    double t_ksp2_total = 0.0;     // e.g., Torso KSP
    double t_ksp3_total = 0.0;     // e.g., Other KSP (if applicable)
    double t_ionic_model_total = 0.0;
    double t_ionic_start, t_ksp1_start, t_ksp2_start, t_ksp3_start;
    double t_ionic_end, t_ksp1_end, t_ksp2_end, t_ksp3_end; 




   if (my_rank == 0)
   {
      std::cout << "Initializing with " << num_ranks << " MPI ranks." << std::endl;
   }
   
   int order = 1;

   std::vector<std::string> objectFilenames;
   if (argc == 1)
      objectFilenames.push_back("femheart.data");

   for (int iargCursor=1; iargCursor<argc; iargCursor++)
      objectFilenames.push_back(argv[iargCursor]);

   if (my_rank == 0) {
      for (int ii=0; ii<objectFilenames.size(); ii++)
	 object_compilefile(objectFilenames[ii].c_str());
   }
   object_Bcast(0,MPI_COMM_WORLD);

   OBJECT* obj = object_find("femheart", "HEART");
   assert(obj != NULL);

   StartTimer("Read the mesh");
   // Read shared global mesh
   mfem::Mesh *mesh = ecg_readMeshptr(obj, "mesh");
   EndTimer();
   int dim = mesh->Dimension();

   std::string torso_mesh_file;
   objectGet(obj, "torso_mesh", torso_mesh_file, "");

   //Fill in the MatrixElementPiecewiseCoefficients
   std::vector<int> heartRegions;
   objectGet(obj,"heart_regions", heartRegions);

   std::vector<double> sigma_m;
   objectGet(obj,"sigma_m",sigma_m);
   assert(heartRegions.size()*3 == sigma_m.size());

   // 读取细胞内电导率
   std::vector<double> sigma_i_values;
   objectGet(obj, "sigma_i", sigma_i_values);
   assert(heartRegions.size()*3 == sigma_i_values.size());

   // 读取细胞外电导率
   std::vector<double> sigma_e_values;
   objectGet(obj, "sigma_e", sigma_e_values);
   assert(heartRegions.size()*3 == sigma_e_values.size());

   double sigma_torso;
   objectGet(obj, "sigma_torso", sigma_torso, "0.2");  // 默认值0.2 mS/mm

   // 检查是否应该求解细胞外电位
   bool solveForUe;
   objectGet(obj, "solve_for_ue", solveForUe, "1");  // 默认启用

   bool solve_torso_model;
   objectGet(obj, "solve_torso", solve_torso_model, "1");  // 默认开启

   int heart_boundary_marker;
   objectGet(obj, "heart_bdry_marker", heart_boundary_marker, "1");  // 默认为1

   double tolerance;
   objectGet(obj, "mesh_tolerance", tolerance, "1e-6");  // 读取容差设置

   double dt;
   objectGet(obj,"dt",dt,"0.01 ms");
   double Bm;
   objectGet(obj,"Bm",Bm,"140"); // 1/mm
   double Cm;
   objectGet(obj,"Cm",Cm,"0.01"); // 1 uF/cm^2 = 0.01 uF/mm^2
 
   std::string reactionName;
   objectGet(obj, "reaction", reactionName, "BetterTT06");

   std::string outputDir;
   objectGet(obj, "outdir", outputDir, ".");
   
   double endTime;
   objectGet(obj, "end_time", endTime, "0 ms");

   double outputRate;
   objectGet(obj, "output_rate", outputRate, "1 ms");

   //double checkpointRate;
   //objectGet(obj, "checkpoint_rate", checkpointRate, "100 ms");

   double initVm;
   objectGet(obj, "init_vm", initVm, "-83");

   bool useNodalIion;
   objectGet(obj, "nodal_ion", useNodalIion, "1");

   Mesh* torso_mesh = nullptr;
   ParMesh* pmesh_torso = nullptr;
   ParFiniteElementSpace* pfespace_torso = nullptr;
   ParGridFunction* gf_ue_torso = nullptr;

  


 








   StimulusCollection stims(dt);
   {
      std::vector<std::string> stimulusNames;
      objectGet(obj, "stimulus", stimulusNames);
      for (auto name : stimulusNames)
      {
         OBJECT* stimobj = object_find(name.c_str(), "STIMULUS");
         assert(stimobj != NULL);
         int numTimes;
         objectGet(stimobj, "n", numTimes, "1");
         double bcl;
         objectGet(stimobj, "bcl", bcl, "0 ms");
         assert(numTimes == 1 || bcl != 0);
         double startTime;
         objectGet(stimobj, "start", startTime, "0 ms");
         double duration;
         objectGet(stimobj, "duration", duration, "1 ms");
         double strength;
         objectGet(stimobj, "strength", strength, "0"); //uA/uF
         assert(strength >= 0);
         std::string location;
         objectGet(stimobj, "where", location, "");
         assert(!location.empty());
         OBJECT* locobj = object_find(location.c_str(), "REGION");
         assert(locobj != NULL);
         std::string regionType;
         objectGet(locobj, "type", regionType, "");
         assert(!regionType.empty());
         shared_ptr<StimulusLocation> stimLoc;
         if (regionType == "ball")
         {
            std::vector<double> center;
            objectGet(locobj, "center", center);
            assert(center.size() == 3);
            double radius;
            objectGet(locobj, "radius", radius, "-1");
            assert(radius >= 0);
            stimLoc = std::make_shared<CenterBallStimulus>(center[0],center[1],center[2],radius);
         }
         else if (regionType == "box")
         {
            std::vector<double> lower;
            objectGet(locobj, "lower", lower);
            assert(lower.size() == 3);
            vector<double> upper;
            objectGet(locobj, "upper", upper);
            assert(upper.size() == 3);
            stimLoc = std::make_shared<BoxStimulus>
               (lower[0], upper[0],
                lower[1], upper[1],
                lower[2], upper[2]);
         }
         shared_ptr<StimulusWaveform> stimWave(new SquareWaveform());
         stims.add(Stimulus(numTimes, startTime, duration, bcl, strength, stimLoc, stimWave));
      }
   }
   
   Timeline timeline(dt, endTime);   

   StartTimer("Setting Attributes");
   mesh->SetAttributes();
   EndTimer();

   StartTimer("Partition Mesh");
   // If I read correctly, pmeshpart will now point to an integer array
   //  containing a partition ID (rank!) for every element ID.
   int *pmeshpart = mesh->GeneratePartitioning(num_ranks);
   EndTimer();


   //Go through all the elements and label the partitioning for the vertices
   std::vector<set<int> > pvertset(mesh->GetNV());
   for (int ielem=0; ielem<mesh->GetNE(); ielem++)
   {
      Array<int> verts;
      mesh->GetElementVertices(ielem, verts);
      for (int ivert=0; ivert<verts.Size(); ivert++)
      {
         pvertset[verts[ivert]].insert(pmeshpart[ielem]);
      }
   }

   std::vector<int> local_extents(num_ranks+1);
   {
      std::vector<int> local_counts(num_ranks, 0);
      for(int i=0; i<mesh->GetNV(); i++)
      {
         if ( ! pvertset[i].empty())
         {
            local_counts[*(pvertset[i].begin())]++;
         }
      }

      local_extents[0] = 0;
      for (int irank=0; irank<num_ranks; irank++)
      {
         local_extents[irank+1] = local_extents[irank]+local_counts[irank];
      }
   }

   std::vector<int> globalvert_from_ranklookup(local_extents[num_ranks]);
   std::vector<int> ghostlocalvert_from_ranklookup(local_extents[num_ranks]);
   {
      std::vector<int> cursor_ghostlocal_from_rank(num_ranks, 0);
      std::vector<int> cursor_ranklookup_from_rank = local_extents;
      for(int i=0; i<mesh->GetNV(); i++)
      {
         if ( ! pvertset[i].empty())
         {
            int irank = *(pvertset[i].begin());
            int ranklookup = cursor_ranklookup_from_rank[irank]++;
            int globalvert = i;
            int ghostlocal = cursor_ghostlocal_from_rank[irank];
            globalvert_from_ranklookup[ranklookup] = globalvert;
            ghostlocalvert_from_ranklookup[ranklookup] = ghostlocal;
            for (const int used_by_this_rank : pvertset[i])
            {
               cursor_ghostlocal_from_rank[used_by_this_rank]++;
            }
         }
      }
   }

   //Get the element material types for each index.
   std::vector<int> material_from_ranklookup(local_extents[num_ranks]);
   {
      std::vector<int> element_from_globalvert(mesh->GetNV(), mesh->GetNE());
      for (int ielem=0; ielem<mesh->GetNE(); ielem++)
      {
         Array<int> verts;
         mesh->GetElementVertices(ielem, verts);
         for (int ivert=0; ivert<verts.Size(); ivert++)
         {
            element_from_globalvert[verts[ivert]] = std::min(element_from_globalvert[verts[ivert]], ielem);
         }
      }
      std::vector<int> cursor_ranklookup_from_rank = local_extents;
      for(int i=0; i<mesh->GetNV(); i++)
      {
         if ( ! pvertset[i].empty())
         {
            int irank = *(pvertset[i].begin());
            int ranklookup = cursor_ranklookup_from_rank[irank]++;
            int globalvert = i;

            int ielem = element_from_globalvert[globalvert];
            material_from_ranklookup[ranklookup] = mesh->GetElement(ielem)->GetAttribute();
         }
      }
   }
   
   if (my_rank == 0)
   {
      for(int i=0; i<num_ranks; i++) {
         std::cout << "Rank " << i << " has " << local_extents[i+1]-local_extents[i] << " nodes!" << std::endl;
      }
   }
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, pmeshpart);
   
   // Build a new FEC...
   FiniteElementCollection *fec;
   if (my_rank == 0) { std::cout << "Creating new FEC..." << std::endl; }
   fec = new H1_FECollection(order, dim);
   // ...and corresponding FES
   ParFiniteElementSpace *pfespace = new ParFiniteElementSpace(pmesh, fec);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   std::cout << "[" << my_rank << "] Number of finite element unknowns: "
	     << pfespace->GetTrueVSize() << std::endl;

   // 5. Determine the list of true (i.e. conforming) essential boundary DOFs
   Array<int> ess_tdof_list;   // Essential true degrees of freedom
   // "true" takes into account shared vertices.
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = 0;
      pfespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   




   // 7. Define the solution vector x as a finite element grid function
   //    corresponding to pfespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   ParGridFunction gf_Vm(pfespace);
   ParGridFunction gf_ue(pfespace);  // 用于细胞外电位的网格函数
   ParGridFunction gf_b(pfespace);
   gf_Vm = initVm;
   gf_ue = 0.0;  // 初始化为零
   gf_b = 0.0;



   //if (solve_torso_model && !torso_mesh_file.empty()) 
    if (my_rank == 0) {
        std::cout << "\n===== 设置Torso模型 =====\n" << std::endl;
        std::cout << "Torso网格文件: " << torso_mesh_file << std::endl;
        std::cout << "Torso电导率: " << sigma_torso << " mS/mm" << std::endl;
        std::cout << "点匹配容差: " << tolerance << std::endl;
    }
    
    // 使用改进的函数读取torso网格并识别边界
    torso_mesh = new Mesh(torso_mesh_file, 1, 1);

            // 设置躯干网格的边界属性
            for (int i = 0; i < torso_mesh->GetNBE(); i++) {
              int bid = torso_mesh->GetBdrAttribute(i);
              if (bid == 1) {
                  torso_mesh->SetBdrAttribute(i, BoundaryType::NEUMANN_ZERO);
              } else if (bid == 100) {
                  torso_mesh->SetBdrAttribute(i, BoundaryType::DIRICHLET);
              }
          }
    
    if (my_rank == 0) {
        std::cout << "Torso网格准备完成，进行分区..." << std::endl;
    }
    
    // 分区torso网格
    int* ptorso_meshpart = torso_mesh->GeneratePartitioning(num_ranks);
    pmesh_torso = new ParMesh(MPI_COMM_WORLD, *torso_mesh, ptorso_meshpart);
    delete torso_mesh;
        // 释放临时内存
    delete[] ptorso_meshpart;


  // 3. 创建网格值传输器
  CoordinateBasedTransfer* transfer = 
  new CoordinateBasedTransfer(pmesh, pfespace, &gf_ue);
        
    
    if (my_rank == 0) {
        std::cout << "创建Torso有限元空间..." << std::endl;
    }
    
    // 为torso创建有限元空间，与heart使用相同的有限元类型
    FiniteElementCollection* torso_fec = new H1_FECollection(1, 3);
    pfespace_torso = new ParFiniteElementSpace(pmesh_torso, torso_fec);

    
    if (my_rank == 0) {
        std::cout << "Torso有限元空间创建完成，自由度数量: " << pfespace_torso->GetTrueVSize() << std::endl;
    }
    


    
    // 创建torso解向量
    gf_ue_torso = new ParGridFunction(pfespace_torso);
    *gf_ue_torso = 0.0;
    

    
    if (my_rank == 0) {
        std::cout << "\n===== Torso模型设置完成 =====\n" << std::endl;
    }

















   // Load fiber quaternions from file
   std::shared_ptr<GridFunction> flat_fiber_quat;
   ecg_readGF(obj, "fibers", mesh, flat_fiber_quat);
   std::shared_ptr<ParGridFunction> fiber_quat;
   fiber_quat = std::make_shared<mfem::ParGridFunction>(pmesh, flat_fiber_quat.get(), pmeshpart);

  std::shared_ptr<GridFunction> flat_sheet_quat;
   ecg_readGF(obj, "sheet", mesh, flat_sheet_quat);
   std::shared_ptr<ParGridFunction> sheet_quat;
   sheet_quat = std::make_shared<mfem::ParGridFunction>(pmesh, flat_sheet_quat.get(), pmeshpart);

   std::shared_ptr<GridFunction> flat_transverse_quat;
   ecg_readGF(obj, "transverse", mesh, flat_transverse_quat);
   std::shared_ptr<ParGridFunction> transverse_quat;
   transverse_quat = std::make_shared<mfem::ParGridFunction>(pmesh, flat_transverse_quat.get(), pmeshpart);
   
   // Load conductivity data
   MatrixElementPiecewiseCoefficient sigma_m_pos_coeffs(fiber_quat, sheet_quat, transverse_quat);
   MatrixElementPiecewiseCoefficient sigma_m_neg_coeffs(fiber_quat, sheet_quat, transverse_quat);
   for (int ii=0; ii<heartRegions.size(); ii++) {
      int heartCursor=3*ii;
      Vector sigma_m_vec(&sigma_m[heartCursor],3);
      Vector sigma_m_pos_vec(3);
      Vector sigma_m_neg_vec(3);
      for (int jj=0; jj<3; jj++)
      {
         double value = sigma_m[heartCursor+jj]*dt/2/Bm/Cm;
         sigma_m_pos_vec[jj] = value;
         sigma_m_neg_vec[jj] = -value;
      }
    
      sigma_m_pos_coeffs.heartConductivities_[heartRegions[ii]] = sigma_m_pos_vec;
      sigma_m_neg_coeffs.heartConductivities_[heartRegions[ii]] = sigma_m_neg_vec;
   }

   // 准备求解细胞外电位
   if (solveForUe) {
      if (my_rank == 0) {
         std::cout << "准备求解细胞外电位..." << std::endl;
         
         // 打印电导率信息
         std::cout << "心脏区域数量: " << heartRegions.size() << std::endl;
         std::cout << "细胞内电导率值数量: " << sigma_i_values.size() << std::endl;
         std::cout << "细胞外电导率值数量: " << sigma_e_values.size() << std::endl;
         
         // 打印一些样本值
         if (!heartRegions.empty()) {
            std::cout << "第一个心脏区域: " << heartRegions[0] << std::endl;
         }
         if (sigma_i_values.size() >= 3) {
            std::cout << "第一组细胞内电导率: " 
                      << sigma_i_values[0] << ", " 
                      << sigma_i_values[1] << ", " 
                      << sigma_i_values[2] << std::endl;
         }
      }
   }

   StartTimer("Forming bilinear system (RHS)");


   ConstantCoefficient one(1.0);
   ParBilinearForm *b = new ParBilinearForm(pfespace);
   b->AddDomainIntegrator(new DiffusionIntegrator(sigma_m_neg_coeffs));
   b->AddDomainIntegrator(new MassIntegrator(one));
   b->Assemble();
   // This creates the linear algebra problem.
   HypreParMatrix RHS_mat;
   b->FormSystemMatrix(ess_tdof_list, RHS_mat);
   EndTimer();

   StartTimer("Forming bilinear system (LHS)");
   
   // Brought out of loop to avoid unnecessary duplication
   ParBilinearForm *a = new ParBilinearForm(pfespace);   // defines a.
   a->AddDomainIntegrator(new DiffusionIntegrator(sigma_m_pos_coeffs));
   a->AddDomainIntegrator(new MassIntegrator(one));
   a->Update(pfespace);
   a->Assemble();
   HypreParMatrix LHS_mat;
   a->FormSystemMatrix(ess_tdof_list,LHS_mat);
   EndTimer();

   //Set up the solve
   HyprePCG pcg(LHS_mat);
   pcg.SetTol(1e-12);
   pcg.SetMaxIter(2000);
   pcg.SetPrintLevel(2);
   HypreSolver *M_test = new HypreBoomerAMG(LHS_mat);
   pcg.SetPreconditioner(*M_test);


   //Set up the ionic models
   ParLinearForm *c = new ParLinearForm(pfespace);
   //positive dt here because the reaction models use dVm = -Iion
   c->AddDomainIntegrator(new DomainLFIntegrator(stims));


   
   
   ThreadServer& threadServer = ThreadServer::getInstance();
   ThreadTeam defaultGroup = threadServer.getThreadTeam(vector<unsigned>());
   std::vector<std::string> reactionNames;
   reactionNames.push_back(reactionName);
   std::vector<int> cellTypes;

   //int Iion_order = 2*order+3;
   int Iion_order = 2*order-1;
   QuadratureSpace quadSpace(pmesh, Iion_order);
   if (useNodalIion)
   {
      for (int ranklookup=local_extents[my_rank]; ranklookup<local_extents[my_rank+1]; ranklookup++)
      {
         cellTypes.push_back(material_from_ranklookup[ranklookup]);
      }
   }
   else
   {
      for (int i = 0; i < pfespace->GetNE(); ++i)
      {
         ElementTransformation *T = pfespace->GetElementTransformation(i);
         //This is a hack.  There's no way to get access to the offsets() array
         //in Quadrature Space without declaring ourselves to be a friend class.
         //This is broken and I hope it is fixed in 4.0
         Vector localVm;
         int NNN = quadSpace.GetElementIntRule(i).GetNPoints();
         for ( int j=0; j<NNN; j++)
         {
            cellTypes.push_back(T->Attribute);
         }
      }
   }

   ReactionWrapper reactionWrapper(dt,reactionNames,defaultGroup,cellTypes);
   reactionWrapper.Initialize();
   cellTypes.clear();
   reactionNames.clear();
   
   ParBilinearForm *Iion_blf;
   HypreParMatrix Iion_mat;
   ConstantCoefficient dt_coeff(dt);
   ReactionFunction* rf = NULL;
   if (useNodalIion) {
      Iion_blf = new ParBilinearForm(pfespace);
      Iion_blf->AddDomainIntegrator(new MassIntegrator(dt_coeff));
      Iion_blf->Update(pfespace);
      Iion_blf->Assemble();
      Iion_blf->FormSystemMatrix(ess_tdof_list,Iion_mat);
   } else {
      Iion_blf = NULL;
      
      rf = new ReactionFunction(&quadSpace,pfespace,&reactionWrapper); 
      c->AddDomainIntegrator(new QuadratureIntegrator(rf, dt)); 
   }

   Vector actual_Vm(pfespace->GetTrueVSize()), actual_b(pfespace->GetTrueVSize()), actual_old(pfespace->GetTrueVSize());
   Vector actual_Iion(pfespace->GetTrueVSize());
   bool first=true;

   if (useNodalIion)
   {
      actual_Vm = reactionWrapper.getVmReadonly();
   }





        // 5. 创建双线性型和线性型
        ParBilinearForm* a_torsoSolver = new ParBilinearForm(pfespace_torso);
        ConstantCoefficient sigma_torso_coeff(sigma_torso);
        a_torsoSolver->AddDomainIntegrator(new DiffusionIntegrator(sigma_torso_coeff));
        a_torsoSolver->Assemble();
        a_torsoSolver->Finalize();

        // 创建线性型
        ParLinearForm* b_torsoSolver = new ParLinearForm(pfespace_torso);
                // 为源边界创建边界属性数组
                Array<int> source_bdr_torso(pmesh_torso->bdr_attributes.Max());
                source_bdr_torso = 0;
                source_bdr_torso[BoundaryType::SOURCE-1] = 1;
        // 6. 应用Dirichlet边界条件
        Array<int> ess_bdr_torso(pmesh_torso->bdr_attributes.Max());
        ess_bdr_torso = 0;
        ess_bdr_torso[BoundaryType::DIRICHLET-1] = 1; 
            Array<int> ess_tdof_list_torso;
        pfespace_torso->GetEssentialTrueDofs(ess_bdr_torso, ess_tdof_list_torso);
// 创建包含边界值的系数函数
BoundaryValuesCoefficient* bdr_coef = new BoundaryValuesCoefficient(transfer); 


    MPI_Barrier(MPI_COMM_WORLD); // Sync before starting overall timer
    t_start = MPI_Wtime();
    t_total_elapsed = t_start; // Store start time here initially
double t_total_start;
double t_total_end;
   
   int itime=0;
   while (1)
   {
      if (my_rank == 0)
      {
         std::cout << "time = " << timeline.realTimeFromTimestep(itime) << std::endl;
      }
      //output if appropriate
      if ((itime % timeline.timestepFromRealTime(outputRate)) == 0)
      {
         // 添加同步点，确保所有进程在输出前完成计算
         MPI_Barrier(MPI_COMM_WORLD);
         
         if (my_rank == 0)
         {
            std::string timedir = outputDir + "/tm" + timeline.outputIdFromTimestep(itime);
            recursive_mkdir(timedir); 

            std::vector<double> dataBuffer(local_extents[num_ranks]);
            for (int irank=0; irank<num_ranks; irank++)
            {
               int local_size = local_extents[irank+1] - local_extents[irank];
               std::vector<double> rankBuffer(local_size);
               if (irank==0)
               {
                  if (local_extents[irank+1] > 0)
                  {
                     //Since rank 0 is always the least, ranklookup == localvert
                     //the following assertion makes sure this is always true.
                     assert(ghostlocalvert_from_ranklookup[local_extents[irank+1]-1] == local_extents[irank+1]-1);
                     memcpy(&rankBuffer[0], &gf_Vm[0], sizeof(double)*local_size);
                  }
               }
               else
               {
                  MPI_Status dontcare;
                  MPI_Recv(&rankBuffer[0], local_size,
                           MPI_DOUBLE, irank, 455, MPI_COMM_WORLD, &dontcare);
               }
               for (int ii=0; ii<local_size; ii++)
               {
                  dataBuffer[globalvert_from_ranklookup[ii+local_extents[irank]]] = rankBuffer[ii];
               }
            }
            std::string VmFilename = timedir + "/Vm.npy";
            save1dNumpyArray(timedir + "/Vm.npy", dataBuffer);
            
            // 也输出u_e（如果我们正在求解它）
            if (solveForUe) {
               std::vector<double> ueDataBuffer(local_extents[num_ranks]);
               for (int irank = 0; irank < num_ranks; irank++) {
                  int local_size = local_extents[irank+1] - local_extents[irank];
                  std::vector<double> rankBuffer(local_size);
                  if (irank == 0) {
                     if (local_extents[irank+1] > 0) {
                        memcpy(&rankBuffer[0], &gf_ue[0], sizeof(double)*local_size);
                     }
                  } else {
                     MPI_Status dontcare;
                     MPI_Recv(&rankBuffer[0], local_size,
                              MPI_DOUBLE, irank, 456, MPI_COMM_WORLD, &dontcare);
                  }
                  for (int ii = 0; ii < local_size; ii++) {
                     ueDataBuffer[globalvert_from_ranklookup[ii+local_extents[irank]]] = rankBuffer[ii];
                  }
               }
               save1dNumpyArray(timedir + "/Ue.npy", ueDataBuffer);
               int valid_values = 0;
            for (int i = 0; i < ueDataBuffer.size(); i++) {
                if (fabs(ueDataBuffer[i]) > 1e-10) valid_values++;
                }
            std::cout << "有效ueDataBuffer电势值节点数: " << valid_values << " 总节点数: " << ueDataBuffer.size() << std::endl;
            }


         }
         else
         {
            int local_size = local_extents[my_rank+1]-local_extents[my_rank];
            std::vector<double> dataFromLocalRanklookup(local_size);
            for (int ii=0; ii<local_size; ii++)
            {
               int ranklookup = local_extents[my_rank] + ii;
               dataFromLocalRanklookup[ii] = gf_Vm[ghostlocalvert_from_ranklookup[ranklookup]];
            }
            MPI_Send(&dataFromLocalRanklookup[0], local_size,
                     MPI_DOUBLE, 0, 455, MPI_COMM_WORLD);
                     
            // 发送u_e数据（如果我们正在求解它）
            if (solveForUe) {
               std::vector<double> ueDataFromLocalRanklookup(local_size);
               for (int ii = 0; ii < local_size; ii++) {
                  int ranklookup = local_extents[my_rank] + ii;
                  ueDataFromLocalRanklookup[ii] = gf_ue[ghostlocalvert_from_ranklookup[ranklookup]];
               }
               MPI_Send(&ueDataFromLocalRanklookup[0], local_size,
                       MPI_DOUBLE, 0, 456, MPI_COMM_WORLD);
            }
         }
         
         // 添加同步障碍，确保所有进程都完成了数据传输
         MPI_Barrier(MPI_COMM_WORLD);
      }
      //if end time, then exit
      if (itime == timeline.maxTimesteps()) { break; }


t_total_start = MPI_Wtime();


t_ionic_start = MPI_Wtime();
      //calculate the ionic contribution.
      if (useNodalIion) {
         reactionWrapper.getVmReadwrite() = actual_Vm; //should be a memcpy
         reactionWrapper.Calc();
      } else {
         rf->Calc(gf_Vm);
      }
 t_ionic_end = MPI_Wtime();
    t_ionic_model_total += (t_ionic_end - t_ionic_start);
      
      //add stimulii
      stims.updateTime(timeline.realTimeFromTimestep(itime));
      
      //compute the Iion and stimulus contribution
      c->Update();
      c->Assemble();
      a->FormLinearSystem(ess_tdof_list, gf_Vm, *c, LHS_mat, actual_Vm, actual_b, 1);
      //compute the RHS matrix contribution
      RHS_mat.Mult(actual_Vm, actual_old);
      actual_b += actual_old;

      if (useNodalIion)
      {
         Iion_mat.Mult(reactionWrapper.getIionReadonly(), actual_old);
         actual_b += actual_old;
      }
      //solve the matrix
      t_ksp1_start = MPI_Wtime();
      pcg.Mult(actual_b, actual_Vm);
      t_ksp1_end = MPI_Wtime();
      t_ksp1_total += (t_ksp1_end - t_ksp1_start);

      a->RecoverFEMSolution(actual_Vm, *c, gf_Vm);
      
      // 求解伪双域模型以恢复细胞外电位u_e
      if (solveForUe) {
          if (my_rank == 0 && itime % 10 == 0) {
              std::cout << "求解细胞外电位(u_e)..." << std::endl;
          }
          
          // 设置统一的打印级别
          int local_print_level = (my_rank == 0 && itime % 50 == 0) ? 1 : 0;
          int global_print_level;
          MPI_Allreduce(&local_print_level, &global_print_level, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
          
          try {
              // 使用统一的打印级别求解细胞外电位
              solvePseudoBidomainForUe(
                  pmesh, pfespace, gf_Vm, gf_ue, 
                  sigma_i_values, sigma_e_values, fiber_quat, sheet_quat, transverse_quat,
                  ess_tdof_list, heartRegions, true,
                  global_print_level, &t_ksp2_total
              );
          } catch (const std::exception& e) {
              // 确保所有进程都知道发生了异常
              if (my_rank == 0) {
                  std::cerr << "求解细胞外电位时发生错误: " << e.what() << std::endl;
              }
          }
      }
      
// 在时间迭代循环中，找到求解伪双域模型的部分后面
if (solve_torso_model  && pmesh_torso && pfespace_torso && gf_ue_torso) {
   if (my_rank == 0 && itime % 10 == 0) {
       std::cout << "\n===== 求解Torso模型 =====\n" << std::endl;
   }
   
   try {

    // 释放旧的传输器
    if (transfer) {
        delete transfer;
        transfer = nullptr;
    }
    if (bdr_coef) {
        delete bdr_coef;
        bdr_coef = nullptr;
    }
    
    // 创建新的传输器，使用最新的心脏解
    transfer = new CoordinateBasedTransfer(pmesh, pfespace, &gf_ue);
    bdr_coef = new BoundaryValuesCoefficient(transfer);




    
//setTorsoIntersectionBoundaryConditions(gf_ue, *gf_ue_torso, ess_tdof_list_torso);

      
       
       // 设置打印级别
       //int local_print_level = (my_rank == 0 && itime % 50 == 0) ? 2 : 1;
       //int global_print_level;
       //MPI_Allreduce(&local_print_level, &global_print_level, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
       
       // 求解torso模型
 //solveTorsoModel(pmesh_torso, pfespace_torso, *gf_ue_torso,
                  //sigma_torso, ess_tdof_list_torso, global_print_level);
       
//torsoSolver(mesh, torso_mesh, gf_ue, *gf_ue_torso, pfespace, pfespace_torso, sigma_torso);
// 求解torso模型


   
        
//delete b_torsoSolver;
//b_torsoSolver = new ParLinearForm(pfespace_torso);
                // 在源边界上添加边界积分
    //ConstantCoefficient source_coeff(1.0);
    //b_torsoSolver->AddBoundaryIntegrator(new BoundaryLFIntegrator(source_coeff),
                                    //source_bdr_torso);
    b_torsoSolver->Assemble();



        // 应用边界条件 - 仅在Dirichlet边界上使用心脏网格的值
        gf_ue_torso->ProjectBdrCoefficient(*bdr_coef, ess_bdr_torso);

        

        // 7. 求解线性系统
        OperatorPtr A_torsoSolver;
        Vector B_torsoSolver, X_torsoSolver;
        a_torsoSolver->FormLinearSystem(ess_tdof_list_torso,
             *gf_ue_torso, *b_torsoSolver, 
             A_torsoSolver, X_torsoSolver, B_torsoSolver);
        
        // 设置求解器
        HypreBoomerAMG* precond_torsoSolver = new HypreBoomerAMG;
        HyprePCG* pcg_torsoSolver = new HyprePCG(MPI_COMM_WORLD);
        
        HypreParMatrix* parcsr_A = A_torsoSolver.As<HypreParMatrix>();
        
        precond_torsoSolver->SetPrintLevel(0);
        pcg_torsoSolver->SetPreconditioner(*precond_torsoSolver);
        pcg_torsoSolver->SetOperator(*parcsr_A);
        pcg_torsoSolver->SetTol(1e-12);
        pcg_torsoSolver->SetMaxIter(1000);
        pcg_torsoSolver->SetPrintLevel(1);
        
        // 求解线性系统
        t_ksp3_start = MPI_Wtime();
        pcg_torsoSolver->Mult(B_torsoSolver, X_torsoSolver);
      t_ksp3_end = MPI_Wtime();
      t_ksp3_total += (t_ksp3_end - t_ksp3_start);
        
        // 8. 恢复解
        a_torsoSolver->RecoverFEMSolution(X_torsoSolver, *b_torsoSolver, *gf_ue_torso);

 t_total_end = MPI_Wtime();
    t_total += (t_total_end - t_total_start);

delete pcg_torsoSolver;
delete precond_torsoSolver;





        if (itime % timeline.timestepFromRealTime(outputRate) == 0) {
        std::string timedir = outputDir + "/tm" + timeline.outputIdFromTimestep(itime);
        string collection_name = "/torso_solution";
        SaveUsingParallelDataCollection(pmesh_torso, gf_ue_torso,
            timedir, collection_name, 0, MPI_COMM_WORLD);
            }






#if 0
       // 可选：输出torso解为VTK文件以便可视化
       if (itime % timeline.timestepFromRealTime(outputRate) == 0) {
           std::string timedir = outputDir + "/tm" + timeline.outputIdFromTimestep(itime);
           std::string vtk_filename = timedir + "/torso_solution.vtk";
           
           if (my_rank == 0) {
               // 确保目录存在
               recursive_mkdir(timedir);
               
               std::ofstream vtk_ofs(vtk_filename);
               if (vtk_ofs.is_open()) {
                   // 保存网格和解
                   pmesh_torso->PrintVTK(vtk_ofs, 0);
                   gf_ue_torso->SaveVTK(vtk_ofs, "u_T", 0);
                   vtk_ofs.close();
                   
                   std::cout << "Torso解已保存为VTK文件: " << vtk_filename << std::endl;
               } else {
                   std::cerr << "无法打开文件进行写入: " << vtk_filename << std::endl;
               }
           }
       }
#endif


   } catch (const std::exception& e) {
       if (my_rank == 0) {
           std::cerr << "处理Torso模型时发生异常: " << e.what() << std::endl;
       }
   }
   
   if (my_rank == 0 && itime % 10 == 0) {
       std::cout << "\n===== Torso模型处理完成 =====\n" << std::endl;
   }


}



      itime++;
      first=false;
   }



    MPI_Barrier(MPI_COMM_WORLD); // Sync before stopping overall timer
    t_end = MPI_Wtime();
    t_total_elapsed = t_end - t_total_elapsed; // Calculate total duration



    double t_total_max;//, t_problem1_max, t_problem2_max, t_problem3_max;
    double t_ksp1_max, t_ksp2_max, t_ksp3_max, t_ionic_model_max;

    MPI_Reduce(&t_total_elapsed, &t_total_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //MPI_Reduce(&t_problem1_total, &t_problem1_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //MPI_Reduce(&t_problem2_total, &t_problem2_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //MPI_Reduce(&t_problem3_total, &t_problem3_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_ksp1_total, &t_ksp1_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_ksp2_total, &t_ksp2_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_ksp3_total, &t_ksp3_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); // Remember this is placeholder
    MPI_Reduce(&t_ionic_model_total, &t_ionic_model_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);





    if (my_rank == 0) {
        std::cout << "\n--- Parallel Timing Results (Max Across " << num_ranks << " Ranks) ---" << std::endl;
        std::cout << std::fixed << std::setprecision(6); // Format output
        std::cout << "Total Simulation Time: " << t_total_max << " s" << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
        //std::cout << "Problem 1 (e.g., Heart): " << t_problem1_max << " s" << std::endl;
        //std::cout << "Problem 2 (e.g., BC Tx): " << t_problem2_max << " s" << std::endl;
        //std::cout << "Problem 3 (e.g., Torso): " << t_problem3_max << " s" << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << "KSP 1 (e.g., Heart LinSolv): " << t_ksp1_max << " s" << std::endl;
        std::cout << "KSP 2 (e.g., Torso LinSolv): " << t_ksp2_max << " s" << std::endl;
        std::cout << "KSP 3 (e.g., Other LinSolv): " << t_ksp3_max << " s" << std::endl; // Adjust name
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << "Ionic Model Calculation:     " << t_ionic_model_max << " s" << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;

        // Optional: Calculate percentage of total time
        if (t_total_max > 1e-9) { // Avoid division by zero
           double ksp_total_max = t_ksp1_max + t_ksp2_max + t_ksp3_max;
           //double problem_sum_max = t_problem1_max + t_problem2_max + t_problem3_max;
           std::cout << "\n--- Percentage of Total Time (Max) ---" << std::endl;
           //std::cout << "Problem 1: " << (t_problem1_max / t_total_max) * 100.0 << "%" << std::endl;
           //std::cout << "Problem 2: " << (t_problem2_max / t_total_max) * 100.0 << "%" << std::endl;
           //std::cout << "Problem 3: " << (t_problem3_max / t_total_max) * 100.0 << "%" << std::endl;
           std::cout << "Ionic Model: " << (t_ionic_model_max / t_total_max) * 100.0 << "%" << std::endl;
           std::cout << "Total KSP: " << (ksp_total_max / t_total_max) * 100.0 << "%" << std::endl;
           std::cout << "--------------------------------------------------" << std::endl;
           // Note: Sum of percentages might not be 100% due to overhead not timed
           //       and KSP/Ionic times being *part of* Problem times.
           //std::cout << "Debug: Sum of Problems: " << problem_sum_max << " s" << std::endl;
           //std::cout << "Debug: KSP1 within Problem1: " << (t_ksp1_max / t_problem1_max) * 100.0 << "%" << std::endl;
           //std::cout << "Debug: Ionic within Problem1: " << (t_ionic_model_max / t_problem1_max) * 100.0 << "%" << std::endl;
           //std::cout << "Debug: KSP2 within Problem3: " << (t_ksp2_max / t_problem3_max) * 100.0 << "%" << std::endl;

        }
    }






   // 14. Free the used memory.
   delete M_test;
   delete a;
   delete b;
   delete c;
   if (rf) delete rf;
   if (Iion_blf) delete Iion_blf;
   delete pfespace;
   if (order > 0) { delete fec; }
   delete mesh;
   delete pmesh;
   delete[] pmeshpart;

   if (gf_ue_torso) delete gf_ue_torso;
    if (pfespace_torso) delete pfespace_torso;
    if (pmesh_torso) delete pmesh_torso;
    if (torso_mesh) delete torso_mesh;
    // 主函数末尾资源清理部分应该还需要添加：
if (transfer) delete transfer;
if (bdr_coef) delete bdr_coef;
delete torso_fec; // 此行未出现在清理代码中

if (a_torsoSolver) delete a_torsoSolver;
if (b_torsoSolver) delete b_torsoSolver;
   
   return 0;
}