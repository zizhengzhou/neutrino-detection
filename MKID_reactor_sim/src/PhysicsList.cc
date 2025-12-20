#include "PhysicsList.hh"
#include "G4EmParameters.hh"

// 引入 Penelope 电磁物理包 (专门用于低能 < 1 keV)
#include "G4EmPenelopePhysics.hh"
// 引入标准电磁物理包 (Shielding 默认用的，我们需要它的头文件来识别并移除它)
#include "G4EmStandardPhysics.hh"

// 引入单位系统
#include "G4SystemOfUnits.hh"

// ---- 构造函数 ----
namespace MKID_reactor_sim
{
PhysicsList::PhysicsList() : Shielding()
{
  // 1. 设置详细级别 (Verbose Level)
  // 设为 1 可以在运行开始时看到 Geant4 打印加载了哪些物理过程，方便检查
  SetVerboseLevel(1);

  // 2. 【核心手术】替换电磁物理列表
  // Shielding 默认加载的是 G4EmStandardPhysics (适合高能物理，不适合 MKID 低能区)
  // 我们使用 ReplacePhysics 方法，强制将其替换为 Penelope
  // G4EmPenelopePhysics 对原子壳层效应、荧光、俄歇电子(Auger)的处理非常精确
  this->ReplacePhysics(new G4EmPenelopePhysics());

  // [关键] 允许 Geant4 追踪极低能粒子
  G4EmParameters* emParams = G4EmParameters::Instance();
  
  // 1. 设置能量下限 (默认为 100 eV 或 1 keV)
  // 为了模拟 10 eV 信号，我们需要把它降得更低
  emParams->SetMinEnergy(10 * eV);
  emParams->SetLowestElectronEnergy(10 * eV);
  
  // 2. 开启荧光和俄歇电子 (MKID 重要信号来源)
  emParams->SetFluo(true);
  emParams->SetAuger(true);
  emParams->SetPixe(true); // 如果关注粒子激发X射线
  // 注意：你不需要手动添加中子或衰变物理，
  // 因为基类 Shielding() 的构造函数里已经自动帮你加载了：
  // - G4HadronPhysicsShielding (含中子 HP 模型)
  // - G4RadioactiveDecayPhysics (含同位素衰变)
}

// ---- 析构函数 ----
PhysicsList::~PhysicsList()
{
}
}