#ifndef PhysicsList_h
#define PhysicsList_h 1

// 引入 Geant4 官方的 Shielding 物理列表作为基类
// Shielding 专门用于地下实验室、屏蔽设计和低本底实验
#include "Shielding.hh"
#include "globals.hh"

// 类定义：我们继承自 "Shielding" 而不是空的 "G4VUserPhysicsList"
namespace MKID_reactor_sim{
class PhysicsList: public Shielding
{
public:
  // 构造函数：我们将在这里进行“器官移植”（替换 EM 物理过程）
  PhysicsList();
  
  // 析构函数
  virtual ~PhysicsList();

  // 设置截断值 (Cuts)：决定粒子何时停止产生次级粒子
  // 对于你的 10eV 需求，这里非常关键
  virtual void SetCuts();
};
}

#endif