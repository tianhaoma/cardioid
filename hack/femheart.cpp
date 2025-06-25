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

#include <string> 

#define StartTimer(x)
#define EndTimer()

using namespace mfem;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;


const double DEFAULT_HEART_POTENTIAL = -83.0;  // 默认心脏电位值（仅作为备用）




// 改进的CoordinateBasedTransfer类
class ParallelCoordinateBasedTransfer {
private:
    ParMesh* source_mesh;
    ParFiniteElementSpace* source_fes;
    ParGridFunction* source_gf;
    
    const double spatial_tol = 1e-8;
    
    // 存储全局边界点信息
    struct GlobalBoundaryPoint {
        double coords[3];
        double value;
        int owner_rank;  // 拥有该点的进程
        
        GlobalBoundaryPoint() {}
        GlobalBoundaryPoint(double x, double y, double z, double v, int rank) 
            : value(v), owner_rank(rank) {
            coords[0] = x; coords[1] = y; coords[2] = z;
        }
    };
    
    std::vector<GlobalBoundaryPoint> global_boundary_points;

    // 用来存储每个进程拥有的边界点数量
    std::vector<int> points_counts_per_rank_; 
    // 用来存储Allgatherv需要的位移信息
    std::vector<int> points_disps_per_rank_;  
    
public:
    ParallelCoordinateBasedTransfer(ParMesh* src_mesh, ParFiniteElementSpace* src_fes, 
                                  ParGridFunction* src_gf)
        : source_mesh(src_mesh), source_fes(src_fes), source_gf(src_gf) {
        BuildGlobalBoundaryPointsMap();
    }
    
    // 构建全局边界点映射
    void BuildGlobalBoundaryPointsMap() {
        int my_rank, num_procs;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
        
        // 步骤1：收集本地边界点
        std::vector<GlobalBoundaryPoint> local_boundary_points;
        std::set<int> processed_vertices;  // 避免重复处理顶点
        
        for (int i = 0; i < source_mesh->GetNBE(); i++) {
            Array<int> vertices;
            source_mesh->GetBdrElementVertices(i, vertices);
            
            for (int j = 0; j < vertices.Size(); j++) {
                int vdof = vertices[j];
                
                // 避免重复处理
                if (processed_vertices.find(vdof) != processed_vertices.end()) {
                    continue;
                }
                processed_vertices.insert(vdof);
                
                // 获取顶点坐标
                double coords[3];
                source_mesh->GetNode(vdof, coords);
                
                // 获取该顶点处的值
                double value = (*source_gf)(vdof);
                
                // 创建边界点
                local_boundary_points.emplace_back(
                    coords[0], coords[1], coords[2], value, my_rank
                );
            }
        }
        
        // 步骤2：收集每个进程的边界点数量
        int local_count = local_boundary_points.size();
        std::vector<int> counts(num_procs);
        MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
        
        // 计算位移和总数
        std::vector<int> displacements(num_procs);
        int total_count = 0;
        for (int i = 0; i < num_procs; i++) {
            displacements[i] = total_count;
            total_count += counts[i];
        }
        
        // 步骤3：创建MPI数据类型用于GlobalBoundaryPoint
        MPI_Datatype mpi_boundary_point_type;
        {
            const int nitems = 3;
            int blocklengths[3] = {3, 1, 1};
            MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
            MPI_Aint offsets[3];
            
            GlobalBoundaryPoint dummy;
            MPI_Aint base_address;
            MPI_Get_address(&dummy, &base_address);
            MPI_Get_address(&dummy.coords[0], &offsets[0]);
            MPI_Get_address(&dummy.value, &offsets[1]);
            MPI_Get_address(&dummy.owner_rank, &offsets[2]);
            
            for (int i = 0; i < nitems; i++) {
                offsets[i] -= base_address;
            }
            
            MPI_Type_create_struct(nitems, blocklengths, offsets, types, 
                                 &mpi_boundary_point_type);
            MPI_Type_commit(&mpi_boundary_point_type);
        }
        
        // 步骤4：收集所有边界点到所有进程
        global_boundary_points.resize(total_count);
        
        MPI_Allgatherv(local_boundary_points.data(), local_count, mpi_boundary_point_type,
                       global_boundary_points.data(), counts.data(), displacements.data(),
                       mpi_boundary_point_type, MPI_COMM_WORLD);
        
        // 清理MPI类型
        MPI_Type_free(&mpi_boundary_point_type);

            // 将计算好并使用过的 counts 和 displacements 保存到成员变量中，以备后用
    points_counts_per_rank_ = counts;
    points_disps_per_rank_ = displacements;

        
        if (my_rank == 0) {
            std::cout << "全局边界点收集完成，总数: " << total_count << std::endl;
        }
    }
    
    // 更新边界值（不重建映射）
    void UpdateBoundaryValues() {
            int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // 阶段一：每个进程收集自己“拥有”的点的最新电势值
    std::vector<double> local_updated_values;
    // 使用预先存储的数量来预分配内存，提高效率
    local_updated_values.reserve(points_counts_per_rank_[my_rank]); 

    // 这个循环仍然需要遍历全局点来识别哪些是自己的
    for (const auto& bp : global_boundary_points) {
        if (bp.owner_rank == my_rank) {
            // 注意：这里的本地搜索逻辑与原来相同，它本身也可能成为一个CPU瓶颈。
            // 但我们当前的首要目标是解决通信瓶颈。
            bool found = false;
            for (int i = 0; i < source_mesh->GetNBE(); i++) {
                Array<int> vertices;
                source_mesh->GetBdrElementVertices(i, vertices);
                for (int j = 0; j < vertices.Size(); j++) {
                    int vdof = vertices[j];
                    double node_coords[3];
                    source_mesh->GetNode(vdof, node_coords);
                    
                    double dist_sq = 0.0;
                    for (int k = 0; k < 3; k++) {
                        dist_sq += (node_coords[k] - bp.coords[k]) * (node_coords[k] - bp.coords[k]);
                    }
                    
                    if (sqrt(dist_sq) < spatial_tol) {
                        local_updated_values.push_back((*source_gf)(vdof));
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
        }
    }

    // 阶段二：使用一次Allgatherv来分发所有更新的值
    std::vector<double> global_updated_values(global_boundary_points.size());
    
    MPI_Allgatherv(
        local_updated_values.data(),           // 发送缓冲区 (本地更新的值)
        local_updated_values.size(),           // 发送数量
        MPI_DOUBLE,                            // 发送类型
        global_updated_values.data(),          // 接收缓冲区 (全局所有值)
        points_counts_per_rank_.data(),        // 每个进程接收多少 (已缓存)
        points_disps_per_rank_.data(),         // 接收数据的位移 (已缓存)
        MPI_DOUBLE,                            // 接收类型
        MPI_COMM_WORLD
    );

    // 阶段三：用收到的全局最新值更新本地的 `global_boundary_points` 列表
    // 因为 Allgatherv 收集数据的顺序(按rank 0, 1, 2...)与 global_boundary_points
    // 最初建立时的顺序是一致的，所以可以直接按索引赋值。
    for (size_t i = 0; i < global_boundary_points.size(); ++i) {
        global_boundary_points[i].value = global_updated_values[i];
    }
    }
    
    // 获取目标点的值
    double GetValueAtPoint(const Vector& point) {
        double min_dist = std::numeric_limits<double>::max();
        double closest_value = DEFAULT_HEART_POTENTIAL;
        
        for (const auto& bp : global_boundary_points) {
            double dist = 0.0;
            for (int i = 0; i < 3; i++) {
                dist += (bp.coords[i] - point(i)) * (bp.coords[i] - point(i));
            }
            dist = sqrt(dist);
            
            if (dist < min_dist) {
                min_dist = dist;
                closest_value = bp.value;
                
                if (dist < spatial_tol) {
                    return closest_value;
                }
            }
        }
        
        return closest_value;
    }
};



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

// 优化后的CoordinateBasedTransfer类定义
class CoordinateBasedTransfer {
private:
    ParMesh* source_mesh;
    ParFiniteElementSpace* source_fes;
    ParGridFunction* source_gf;
    
    // 用于空间查找的容忍度
    const double spatial_tol = 1e-8;
    
    // 存储源网格边界点的坐标、索引和值
    struct BoundaryPoint {
        Vector coords;
        int vertex_id;  // 添加顶点ID以便更新
        double value;
        
        BoundaryPoint(const Vector& c, int id, double v) : coords(c), vertex_id(id), value(v) {}
    };
    
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
                double value = (*source_gf)(vdof);
                
                // 存储顶点ID、坐标和值
                boundary_points.emplace_back(coords, vdof, value);
            }
        }
    }
    
    // 新方法：更新边界点的值而不重建映射
    void UpdateBoundaryValues() {
        for (auto& bp : boundary_points) {
            bp.value = (*source_gf)(bp.vertex_id);
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



// 修改 BoundaryValuesCoefficient 类，添加对两种传输类的支持
class BoundaryValuesCoefficient : public Coefficient {
private:
    CoordinateBasedTransfer* transfer;
    ParallelCoordinateBasedTransfer* parallel_transfer;
    bool use_parallel;
    
public:
    // 原有构造函数
    BoundaryValuesCoefficient(CoordinateBasedTransfer* t) 
        : transfer(t), parallel_transfer(nullptr), use_parallel(false) {}
    
    // 新增构造函数
    BoundaryValuesCoefficient(ParallelCoordinateBasedTransfer* t) 
        : transfer(nullptr), parallel_transfer(t), use_parallel(true) {}
    
    virtual double Eval(ElementTransformation& T, const IntegrationPoint& ip) {
        Vector physical_coord(3);
        T.Transform(ip, physical_coord);
        
        if (use_parallel) {
            return parallel_transfer->GetValueAtPoint(physical_coord);
        } else {
            return transfer->GetValueAtPoint(physical_coord);
        }
    }
};






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



int main(int argc, char *argv[])
{
   MPI_Init(NULL,NULL);
   int num_ranks, my_rank;
   MPI_Comm_size(COMM_LOCAL,&num_ranks);
   MPI_Comm_rank(COMM_LOCAL,&my_rank);

   units_internal(1e-3, 1e-9, 1e-3, 1e-3, 1, 1e-9, 1);
   units_external(1e-3, 1e-9, 1e-3, 1e-3, 1, 1e-9, 1);

   bool use_petsc = true;//false;true
   const char *petscrc_file = "rc_fem_heart";
   MFEMInitializePetsc(NULL,NULL,petscrc_file,NULL);


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

    static int total_iterations_monodomain = 0;
    static int solve_count_monodomain = 0;
    static int total_iterations_recoverue = 0;
    static int solve_count_recoverue = 0;
    static int total_iterations_torso = 0;
    static int solve_count_torso = 0;




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

   //const char *coll_name = "par-data-" << std::setfill('0') << std::setw(6) << num_ranks;

   Mesh* torso_mesh = nullptr;
   ParMesh* pmesh_torso = nullptr;
   ParFiniteElementSpace* pfespace_torso = nullptr;
   ParGridFunction* gf_ue_torso = nullptr;


std::string data_path = "parData/";
std::ostringstream oss_coll_name;
oss_coll_name << "par-data-" << std::setfill('0') << std::setw(6) << num_ranks;
std::string coll_name = oss_coll_name.str();
VisItDataCollection visit_dc(MPI_COMM_WORLD, coll_name);
visit_dc.SetPrefixPath(data_path);




   ParMesh *pmesh;
   ParGridFunction *saved_fiber;
   ParGridFunction *saved_sheet;
   ParGridFunction *saved_transverse;

   visit_dc.Load();
    //cout << "visit_dc Loaded;" << endl;
    pmesh = dynamic_cast<ParMesh*>(visit_dc.GetMesh());

    int ne_before_ = pmesh->GetNE();
//pmesh->UniformRefinement();
//pmesh->UniformRefinement();
//pmesh->UniformRefinement();
int ne_after_ = pmesh->GetNE();

if (my_rank == 0) {
    cout << "heart细化前单元数: " << ne_before_ << endl;
    cout << "heart细化后单元数: " << ne_after_ << endl;
    cout << "增长倍数: " << (double)ne_after_/ne_before_ << endl;
}

   //saved_fiber = visit_dc.GetParField("fiber");
   //saved_sheet = visit_dc.GetParField("sheet");
   //saved_transverse = visit_dc.GetParField("trans");
   //cout << "visit_dc Got ParField" << endl;


    std::ostringstream oss_coll_name_torso;
    oss_coll_name_torso << "par-torso-data-" << std::setfill('0') << std::setw(6) << num_ranks;
    std::string coll_name_torso = oss_coll_name_torso.str();
    VisItDataCollection visit_dc_torso(MPI_COMM_WORLD, coll_name_torso);
    visit_dc_torso.SetPrefixPath(data_path);
    visit_dc_torso.Load();
    //cout << "visit_dc_torso Loaded;" << endl;
    pmesh_torso = dynamic_cast<ParMesh*>(visit_dc_torso.GetMesh());
    
    // 验证细化前后的单元数量
int ne_before = pmesh_torso->GetNE();
//pmesh_torso->UniformRefinement();
//pmesh_torso->UniformRefinement();
//pmesh_torso->UniformRefinement();
int ne_after = pmesh_torso->GetNE();

if (my_rank == 0) {
    cout << "细化前单元数: " << ne_before << endl;
    cout << "细化后单元数: " << ne_after << endl;
    cout << "增长倍数: " << (double)ne_after/ne_before << endl;
}


   // Read shared global mesh
   //mfem::Mesh *mesh = nullptr;
   if (my_rank == 0)
   {
    //mesh = ecg_readMeshptr(obj, "mesh");
   }
   //mfem::Mesh *mesh = ecg_readMeshptr(obj, "mesh");
   EndTimer();
   int dim = pmesh->Dimension();

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
   pmesh->SetAttributes();
   EndTimer();

   //StartTimer("Partition Mesh");
   // If I read correctly, pmeshpart will now point to an integer array
   //  containing a partition ID (rank!) for every element ID.
   //int *pmeshpart = mesh->GeneratePartitioning(num_ranks);



   
   if (my_rank == 0)
   {
      for(int i=0; i<num_ranks; i++) {
         //std::cout << "Rank " << i << " has " << local_extents[i+1]-local_extents[i] << " nodes!" << std::endl;
      }
   }
   //ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, pmeshpart);
   
   // Build a new FEC...
   FiniteElementCollection *fec;
   if (my_rank == 0) { std::cout << "Creating new FEC..." << std::endl; }
   fec = new H1_FECollection(order, dim);
   // ...and corresponding FES
   ParFiniteElementSpace *pfespace = new ParFiniteElementSpace(pmesh, fec);
   //FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
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

   

// 1. Create a new VECTOR finite element space for the 3D fiber directions.
//    We use the same FE collection but specify a vector dimension of 3.
const int vdim = 3;
ParFiniteElementSpace *pfespace_vec = new ParFiniteElementSpace(pmesh, fec, vdim);

// 2. Create shared_ptrs for the ParGridFunctions that will hold the direction vectors.
//    These are created on the vector FESpace.
auto fiber_quat = std::make_shared<ParGridFunction>(pfespace_vec);
auto sheet_quat = std::make_shared<ParGridFunction>(pfespace_vec);
auto transverse_quat = std::make_shared<ParGridFunction>(pfespace_vec);

// 3. Define the constant vectors for each direction.
Vector fiber_direction(vdim);
fiber_direction(0) = 1.0; fiber_direction(1) = 0.0; fiber_direction(2) = 0.0;
VectorConstantCoefficient fiber_coeff(fiber_direction);

Vector sheet_direction(vdim);
sheet_direction(0) = 0.0; sheet_direction(1) = 1.0; sheet_direction(2) = 0.0;
VectorConstantCoefficient sheet_coeff(sheet_direction);

Vector trans_direction(vdim);
trans_direction(0) = 0.0; trans_direction(1) = 0.0; trans_direction(2) = 1.0;
VectorConstantCoefficient trans_coeff(trans_direction);

// 4. Project these constant vector coefficients onto the grid functions.
//    This assigns the specified vector to every node in the mesh.
fiber_quat->ProjectCoefficient(fiber_coeff);
sheet_quat->ProjectCoefficient(sheet_coeff);
transverse_quat->ProjectCoefficient(trans_coeff);












   // 7. Define the solution vector x as a finite element grid function
   //    corresponding to pfespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   ParGridFunction gf_Vm(pfespace);
   ParGridFunction gf_ue(pfespace);  // 用于细胞外电位的网格函数
   ParGridFunction gf_b(pfespace);
   gf_Vm = initVm;
   gf_ue = 0.0;  // 初始化为零
   gf_b = 0.0;


   
   

    if (my_rank == 0) {
        std::cout << "\n===== 设置Torso模型 =====\n" << std::endl;
        std::cout << "Torso电导率: " << sigma_torso << " mS/mm" << std::endl;
        std::cout << "点匹配容差: " << tolerance << std::endl;
    }
    
    // 使用改进的函数读取torso网格并识别边界
    //torso_mesh = new Mesh(torso_mesh_file, 1, 1);

            // 设置躯干网格的边界属性
           // for (int i = 0; i < torso_mesh->GetNBE(); i++) {
           //   int bid = torso_mesh->GetBdrAttribute(i);
            //  if (bid == 1) {
             //     torso_mesh->SetBdrAttribute(i, BoundaryType::NEUMANN_ZERO);
             // } else if (bid == 100) {
              //    torso_mesh->SetBdrAttribute(i, BoundaryType::DIRICHLET);
             // }}
             // 设置躯干网格的边界属性
for (int i = 0; i < pmesh_torso->GetNBE(); i++) {
    int bid = pmesh_torso->GetBdrAttribute(i);
    if (bid == 1) {
        pmesh_torso->SetBdrAttribute(i, BoundaryType::NEUMANN_ZERO);
    } else if (bid == 100) {
        pmesh_torso->SetBdrAttribute(i, BoundaryType::DIRICHLET);
    }
}
    
    if (my_rank == 0) {
        std::cout << "Torso网格准备完成，进行分区..." << std::endl;
    }
    



    // 在时间循环之前初始化交界面处理对象 - 添加这段代码
    //CoordinateBasedTransfer* transfer = nullptr;
    //BoundaryValuesCoefficient* bdr_coef = nullptr;
    ParallelCoordinateBasedTransfer* parallel_transfer = nullptr;
    BoundaryValuesCoefficient* bdr_coef = nullptr;
        
    
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



   //std::shared_ptr<ParGridFunction> fiber_quat  = std::make_shared<mfem::ParGridFunction>(*saved_fiber);
   //std::shared_ptr<ParGridFunction> sheet_quat  = std::make_shared<mfem::ParGridFunction>(*saved_sheet);
   //std::shared_ptr<ParGridFunction> transverse_quat  = std::make_shared<mfem::ParGridFunction>(*saved_transverse);

   
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
         double value = sigma_m[jj]*dt/2/Bm/Cm;
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
if(use_petsc)//use petsc
{
   a->SetOperatorType(Operator::PETSC_MATAIJ);
}
   a->AddDomainIntegrator(new DiffusionIntegrator(sigma_m_pos_coeffs));
   a->AddDomainIntegrator(new MassIntegrator(one));
   a->Update(pfespace);
   a->Assemble(use_petsc ? 0 : 1);

   HypreParMatrix LHS_mat;
   HyprePCG* pcg = nullptr;
   HypreSolver *M_test = nullptr;

   PetscPCGSolver* pcg_monodomain_petsc = nullptr;
   PetscParMatrix LHS_monodomain_petsc;
if(!use_petsc)//use petsc
{
   a->FormSystemMatrix(ess_tdof_list,LHS_mat);
   pcg = new HyprePCG(LHS_mat);
   pcg->SetTol(1e-6);
   pcg->SetMaxIter(1000);
   pcg->SetPrintLevel(2);
   M_test = new HypreBoomerAMG(LHS_mat);
   pcg->SetPreconditioner(*M_test);
}
else
{
    a->FormSystemMatrix(ess_tdof_list, LHS_monodomain_petsc);
    pcg_monodomain_petsc = new PetscPCGSolver(MPI_COMM_WORLD);
   pcg_monodomain_petsc->SetOperator(LHS_monodomain_petsc);
   pcg_monodomain_petsc->SetRelTol(1e-6);
   //pcg_monodomain_petsc->SetAbsTol(1e-12);
   pcg_monodomain_petsc->SetMaxIter(1000);
   pcg_monodomain_petsc->SetPrintLevel(2);
    pcg_monodomain_petsc->iterative_mode = true; 
}
   EndTimer();








   //Set up the ionic models
   ParLinearForm *c = new ParLinearForm(pfespace);
   //positive dt here because the reaction models use dVm = -Iion
   c->AddDomainIntegrator(new DomainLFIntegrator(stims));


   
   
   ThreadServer& threadServer = ThreadServer::getInstance();
   ThreadTeam defaultGroup = threadServer.getThreadTeam(vector<unsigned>());
   std::vector<std::string> reactionNames;
   objectGet(obj, "reaction", reactionNames);
   //reactionNames.push_back(reactionName);
   std::vector<int> cellTypes;

   //int Iion_order = 2*order+3;
   int Iion_order = 2*order-1;
   QuadratureSpace quadSpace(pmesh, Iion_order);
   if (useNodalIion)
   {
      //for (int ranklookup=local_extents[my_rank]; ranklookup<local_extents[my_rank+1]; ranklookup++)
      for (int i = 0; i < pfespace->GetNE(); i++)
      {
         //cellTypes.push_back(material_from_ranklookup[ranklookup]);
         cellTypes.push_back(1);
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



 // 创建电导率系数（在时间循环外，只创建一次）
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

// 设置左侧矩阵: -∇·((σ_i + σ_e)∇u_e)
ParBilinearForm *a_pblf_recoverue = new ParBilinearForm(pfespace);
if(use_petsc)
{
    a_pblf_recoverue->SetOperatorType(Operator::PETSC_MATAIJ);
}
a_pblf_recoverue->AddDomainIntegrator(new DiffusionIntegrator(sigma_sum));
a_pblf_recoverue->Assemble(use_petsc ? 0 : 1);  // 注意这里的参数，与torso一致
a_pblf_recoverue->Finalize();

// 设置右侧向量的临时双线性形式: ∇·(σ_i∇V_m)
ParBilinearForm* temp_form = new ParBilinearForm(pfespace);  // 改为指针，与torso风格一致
if(use_petsc)
{
    temp_form->SetOperatorType(Operator::PETSC_MATAIJ);
}
temp_form->AddDomainIntegrator(new DiffusionIntegrator(sigma_i));
temp_form->Assemble(use_petsc ? 0 : 1);
temp_form->Finalize();

// 声明矩阵和求解器（与torso完全一致的风格）
HypreBoomerAMG* precond_recoverue_hypre = nullptr;
HyprePCG* pcg_recoverue_hypre = nullptr;
HypreParMatrix A_recoverue_hypre;  // 栈对象，与torso一致

PetscPCGSolver* pcg_recoverue_petsc = nullptr;
PetscParMatrix A_recoverue_petsc;  // 栈对象，与torso一致

// 临时矩阵也用相同风格
HypreParMatrix A_temp_hypre;
PetscParMatrix A_temp_petsc;

Vector B_recoverue, X_recoverue;
Vector vm_true(pfespace->GetTrueVSize());
Vector rhs_recoverue(pfespace->GetTrueVSize());

// 初始化求解器
if (!use_petsc)
{
    precond_recoverue_hypre = new HypreBoomerAMG;
    pcg_recoverue_hypre = new HyprePCG(MPI_COMM_WORLD);
    
    // 形成系统矩阵
    a_pblf_recoverue->FormSystemMatrix(ess_tdof_list, A_recoverue_hypre);
    temp_form->FormSystemMatrix(ess_tdof_list, A_temp_hypre);
    
    precond_recoverue_hypre->SetPrintLevel(0);
    pcg_recoverue_hypre->SetPreconditioner(*precond_recoverue_hypre);
    pcg_recoverue_hypre->SetOperator(A_recoverue_hypre);
    pcg_recoverue_hypre->SetTol(1e-6);
    pcg_recoverue_hypre->SetMaxIter(1000);
    pcg_recoverue_hypre->SetPrintLevel(1);
}
else
{
    // 形成系统矩阵（PETSc版本）
    a_pblf_recoverue->FormSystemMatrix(ess_tdof_list, A_recoverue_petsc);
    temp_form->FormSystemMatrix(ess_tdof_list, A_temp_petsc);
    
    pcg_recoverue_petsc = new PetscPCGSolver(MPI_COMM_WORLD, "recoverue_", true);
    // 如果PetscPCGSolver需要SetOperator，添加：
    pcg_recoverue_petsc->SetOperator(A_recoverue_petsc);
}

X_recoverue.SetSize(pfespace->GetTrueVSize());
X_recoverue = 0.0;  // 初始猜测为零





        // 5. 创建双线性型和线性型
        ParBilinearForm* a_pblf_torso = new ParBilinearForm(pfespace_torso);
if(use_petsc)
{
   a_pblf_torso->SetOperatorType(Operator::PETSC_MATAIJ);
}

        ConstantCoefficient sigma_torso_coeff(sigma_torso);
        a_pblf_torso->AddDomainIntegrator(new DiffusionIntegrator(sigma_torso_coeff));
        a_pblf_torso->Assemble(use_petsc ? 0 : 1);
        a_pblf_torso->Finalize();

        // 创建线性型
        ParLinearForm* b_plf_torso = new ParLinearForm(pfespace_torso);
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

    if (solve_torso_model && pmesh && pfespace && pfespace_torso) {
        if (my_rank == 0) {
            std::cout << "初始化心脏-躯干交界面传输对象..." << std::endl;
        }
        
        // 创建对象一次，之后只更新值
        //transfer = new CoordinateBasedTransfer(pmesh, pfespace, &gf_ue);
        //bdr_coef = new BoundaryValuesCoefficient(transfer);
            parallel_transfer = new ParallelCoordinateBasedTransfer(pmesh, pfespace, &gf_ue);
            bdr_coef = new BoundaryValuesCoefficient(parallel_transfer);
        
        if (my_rank == 0) {
            std::cout << "心脏-躯干交界面传输对象初始化完成。" << std::endl;
        }
    }

//Set matrix and vector for Torso model
HypreBoomerAMG* precond_hypre = nullptr;
HyprePCG* pcg_hypre = nullptr;
HypreParMatrix A_torso_hypre;

PetscPCGSolver* pcg_petsc = nullptr;
PetscParMatrix A_torso_petsc;

Vector B_torso, X_torso;

   if (!use_petsc)
   {
        precond_hypre = new HypreBoomerAMG;
        pcg_hypre = new HyprePCG(MPI_COMM_WORLD);    
        //parcsr_A_hypre = parcsr_A_hypre.As<HypreParMatrix>();
        a_pblf_torso->FormSystemMatrix(ess_tdof_list_torso, A_torso_hypre);

        precond_hypre->SetPrintLevel(0);
        pcg_hypre->SetPreconditioner(*precond_hypre);
        pcg_hypre->SetOperator(A_torso_hypre);
        pcg_hypre->SetTol(1e-4);
        pcg_hypre->SetMaxIter(1000);
        pcg_hypre->SetPrintLevel(1);
   }
   else
   {
        pcg_petsc = new PetscPCGSolver(MPI_COMM_WORLD, "torso_", true);
   }




    MPI_Barrier(MPI_COMM_WORLD); // Sync before starting overall timer
    t_start = MPI_Wtime();
    t_total_elapsed = t_start; // Store start time here initially
double t_total_start;
double t_total_end;



ParaViewDataCollection paraview_dc("potentials_data", pmesh);
paraview_dc.SetPrefixPath(outputDir);
paraview_dc.RegisterField("Vm", &gf_Vm);
paraview_dc.RegisterField("ue", &gf_ue);
paraview_dc.SetDataFormat(VTKFormat::ASCII);
paraview_dc.SetCycle(0);
paraview_dc.SetTime(0.0);

ParaViewDataCollection paraview_dc_torso("torso_data", pmesh_torso);
paraview_dc_torso.SetPrefixPath(outputDir);
paraview_dc_torso.RegisterField("uT", gf_ue_torso);
paraview_dc_torso.SetDataFormat(VTKFormat::ASCII);
paraview_dc_torso.SetCycle(0);
paraview_dc_torso.SetTime(0.0);

   
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
         std::string timedir = outputDir + "/gf_Data" + "/tm" + timeline.outputIdFromTimestep(itime);
         if (my_rank == 0)
         { recursive_mkdir(timedir); }

      {
         std::string gf_filename = timedir + "/gf_Vm";
         std::string gf_ue_filename = timedir + "/gf_ue";
         gf_Vm.SaveAsOne(gf_filename.c_str());
         gf_ue.SaveAsOne(gf_ue_filename.c_str());

        std::string gf_ue_torso_filename = timedir + "/gf_ue_torso";
        gf_ue_torso->SaveAsOne(gf_ue_torso_filename.c_str());


      }

      paraview_dc.SetCycle(itime);
      paraview_dc.SetTime(timeline.realTimeFromTimestep(itime));
      paraview_dc.Save();

      paraview_dc_torso.SetCycle(itime);
      paraview_dc_torso.SetTime(timeline.realTimeFromTimestep(itime));
      paraview_dc_torso.Save();


         
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

   if (!use_petsc)
   {
      a->FormLinearSystem(ess_tdof_list, gf_Vm, *c, LHS_mat, actual_Vm, actual_b, 1);
    }
    else
    {
        a->FormLinearSystem(ess_tdof_list, gf_Vm, *c, LHS_monodomain_petsc, actual_Vm, actual_b, 1);
    }


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
         if (!use_petsc)
   {
      pcg->Mult(actual_b, actual_Vm);
    }
    else
    {
        
        pcg_monodomain_petsc->Mult(actual_b, actual_Vm);
                int current_iterations = pcg_monodomain_petsc->GetNumIterations();
    total_iterations_monodomain += current_iterations;
    solve_count_monodomain++;
    

    }
      t_ksp1_end = MPI_Wtime();
      t_ksp1_total += (t_ksp1_end - t_ksp1_start);

      a->RecoverFEMSolution(actual_Vm, *c, gf_Vm);
      
      // 求解伪双域模型以恢复细胞外电位u_e
      if (solveForUe) {
          if (my_rank == 0 && itime % 10 == 0) {
              std::cout << "求解细胞外电位(u_e)..." << std::endl;
          }
          
          // 设置统一的打印级别
          int local_print_level = (my_rank == 0 && itime % 50 == 0) ? 3 : 0;
          int global_print_level;
          MPI_Allreduce(&local_print_level, &global_print_level, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
          



    gf_Vm.GetTrueDofs(vm_true);
    
    // 计算右侧向量: -∇·(σ_i∇V_m)
    if (!use_petsc) {
        A_temp_hypre.Mult(-1.0, vm_true, 0.0, rhs_recoverue);
    } else {
        A_temp_petsc.Mult(-1.0, vm_true, 0.0, rhs_recoverue);
    }
    
    // 记录求解时间
    t_ksp2_start = MPI_Wtime();
    
    // 求解线性系统
    if (!use_petsc) {
        pcg_recoverue_hypre->Mult(rhs_recoverue, X_recoverue);
    } else {
        pcg_recoverue_petsc->Mult(rhs_recoverue, X_recoverue);
        int current_iterations = pcg_recoverue_petsc->GetNumIterations();
    total_iterations_recoverue += current_iterations;
    solve_count_recoverue++;
    
    }
    
    t_ksp2_end = MPI_Wtime();
    t_ksp2_total += (t_ksp2_end - t_ksp2_start);
    
    // 更新网格函数
    gf_ue.SetFromTrueDofs(X_recoverue);

    // 如果需要，强制解的均值为零
    if (ess_tdof_list.Size() == 0) {//enforce zero mean
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
            if (my_rank == 0) {
                std::cout << "调整解向量以确保均值为零，当前均值: " << mean_value << std::endl;
            }
            
            // 从解中减去平均值
            for (int i = 0; i < X_recoverue.Size(); i++) {
                X_recoverue(i) -= mean_value;
            }
            
            // 更新网格函数
            gf_ue.SetFromTrueDofs(X_recoverue);
            
            if (my_rank == 0) {
                std::cout << "已从解中减去均值: " << mean_value << std::endl;
            }
        } else {
            if (my_rank == 0) {
                std::cout << "解向量均值已经接近零: " << mean_value << "，无需调整。" << std::endl;
            }
        }
    }


      }
      
// 在时间迭代循环中，找到求解伪双域模型的部分后面
if (solve_torso_model  && pmesh_torso && pfespace_torso && gf_ue_torso) {
   if (my_rank == 0 && itime % 10 == 0) {
       std::cout << "\n===== 求解Torso模型 =====\n" << std::endl;
   }
   
double t_boundary_start_1, t_boundary_end_1, t_boundary_duration_1;
double t_boundary_start_2, t_boundary_end_2, t_boundary_duration_2;
t_boundary_start_1 = MPI_Wtime();

 parallel_transfer->UpdateBoundaryValues();
 t_boundary_end_1 = MPI_Wtime();
t_boundary_duration_1 = t_boundary_end_1 - t_boundary_start_1;

    b_plf_torso->Assemble();

    t_boundary_start_2 = MPI_Wtime();
  // 应用边界条件 - 仅在Dirichlet边界上使用心脏网格的值
    gf_ue_torso->ProjectBdrCoefficient(*bdr_coef, ess_bdr_torso);

 t_boundary_end_2 = MPI_Wtime();
t_boundary_duration_2 = t_boundary_end_2 - t_boundary_start_2;

if (my_rank == 0) {
    std::cout << "UpdateBoundaryValues 运行时间: " << t_boundary_duration_1 << " 秒" << std::endl;
    std::cout << "ProjectBdrCoefficient 运行时间: " << t_boundary_duration_2 << " 秒" << std::endl;
}



   if (!use_petsc)
   {
        a_pblf_torso->FormLinearSystem(ess_tdof_list_torso, *gf_ue_torso, *b_plf_torso, 
             A_torso_hypre, X_torso, B_torso);
        t_ksp3_start = MPI_Wtime();
        pcg_hypre->Mult(B_torso, X_torso);
        t_ksp3_end = MPI_Wtime();
        t_ksp3_total += (t_ksp3_end - t_ksp3_start);
   }
   else
   {
      a_pblf_torso->FormLinearSystem(ess_tdof_list_torso, *gf_ue_torso, *b_plf_torso,
                A_torso_petsc, X_torso, B_torso);
      //pcg_petsc = new PetscPCGSolver(*A_torso_petsc);
        
        pcg_petsc->SetOperator(A_torso_petsc);

              //pcg_petsc->iterative_mode = true; // iterative_mode is true by default with CGSolver
      pcg_petsc->SetRelTol(1e-4);
      pcg_petsc->SetAbsTol(1e-10);
      pcg_petsc->SetMaxIter(3000);
      pcg_petsc->SetPrintLevel(3); // 提高打印级别

      t_ksp3_start = MPI_Wtime();
      pcg_petsc->Mult(B_torso, X_torso);
      t_ksp3_end = MPI_Wtime();
      t_ksp3_total += (t_ksp3_end - t_ksp3_start);
          // 获取PETSc求解器的迭代次数
    int current_iterations = pcg_petsc->GetNumIterations();
    total_iterations_torso += current_iterations;
    solve_count_torso++;
    


   }

        
        // 8. 恢复解
        a_pblf_torso->RecoverFEMSolution(X_torso, *b_plf_torso, *gf_ue_torso);

 t_total_end = MPI_Wtime();
    t_total += (t_total_end - t_total_start);





   
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
        std::cout << "monodomain 平均迭代次数: " << (double)total_iterations_monodomain / solve_count_monodomain << std::endl;
        std::cout << "re ue 平均迭代次数: " << (double)total_iterations_recoverue / solve_count_recoverue << std::endl;
        std::cout << "torso 平均迭代次数: " << (double)total_iterations_torso / solve_count_torso << std::endl;
    }



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
   delete pcg;
   delete a;
   delete b;
   delete c;
   if (rf) delete rf;
   if (Iion_blf) delete Iion_blf;
   delete pfespace;
   //delete fespace;
delete torso_fec;
if (pcg_monodomain_petsc) delete pcg_monodomain_petsc;
   if (order > 0) { delete fec; }
   //delete mesh;
   delete pmesh;
   //delete[] pmeshpart;

   if (gf_ue_torso) delete gf_ue_torso;
    if (pfespace_torso) delete pfespace_torso;
    if (pmesh_torso) delete pmesh_torso;
    //if (torso_mesh) delete torso_mesh;
    // 主函数末尾资源清理部分应该还需要添加：
//if (transfer) delete transfer;
if (parallel_transfer) delete parallel_transfer;
if (bdr_coef) delete bdr_coef;
delete torso_fec; // 此行未出现在清理代码中

if (a_pblf_torso) delete a_pblf_torso;
if (b_plf_torso) delete b_plf_torso;
delete pcg_hypre;
delete precond_hypre;
delete pcg_petsc;

delete pfespace_vec;


MFEMFinalizePetsc();
MPI_Finalize();
   
   return 0;
}