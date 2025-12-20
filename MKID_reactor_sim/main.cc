// =================== 1. 头文件引入 ===================
#include "G4RunManagerFactory.hh"  // 运行管理器工厂（负责创建引擎）
#include "G4UImanager.hh"          // UI管理器（负责处理命令）
#include "G4VisExecutive.hh"       // 可视化主管（负责画图）
#include "G4UIExecutive.hh"        // UI主管（负责弹窗界面）

#include "Randomize.hh"            // 随机数工具

// 引入你自己写的三个核心类
#include "DetectorConstruction.hh" // 几何
#include "PhysicsList.hh"          // 物理
#include "ActionInitialization.hh" // 动作

// =================== 2. 主函数 ===================
int main(int argc, char** argv)
{
  // 3. 初始化随机数引擎 (可选但推荐)
  // 避免每次运行结果完全一样（伪随机）
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // =================== 4. 创建运行管理器 (RunManager) ===================
  // 这是 Geant4 的“心脏”。
  // SerialOnly   = 强制单线程 (调试用)
  // Default      = 自动检测 (如果电脑支持多线程，自动开启 MT)
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // =================== 5. 注册用户类 (User Initialization) ===================
  // 必须告诉心脏：世界长什么样？物理定律是什么？我们要干什么？

  // A. 注册几何 (DetectorConstruction)
  // 你之前写的 DetectorConstruction.cc 在这里被加载
  runManager->SetUserInitialization(new MKID_reactor_sim::DetectorConstruction());

  // B. 注册物理列表 (PhysicsList)
  // 你写的 Shielding + Penelope 混合列表在这里被加载
  runManager->SetUserInitialization(new MKID_reactor_sim::PhysicsList());

  // C. 注册动作初始化 (ActionInitialization)
  // 这里面包含了 PrimaryGenerator, Run, Event, Stepping Actions
  runManager->SetUserInitialization(new MKID_reactor_sim::ActionInitialization());

  // =================== 6. 初始化可视化 (VisManager) ===================
  // 如果没有它，你只能跑数据，看不到图形界面
  G4VisManager* visManager = new G4VisExecutive;
  // Initialize() 会检查显卡驱动、OpenGL 等
  visManager->Initialize();

  // =================== 7. 获取 UI 指针 ===================
  // 用来执行 "/run/beamOn" 等命令
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // =================== 8. 模式选择：界面 vs 批处理 ===================
  // argc 是参数个数。
  // ./mkid_sim          -> argc=1 (交互模式)
  // ./mkid_sim run.mac  -> argc=2 (批处理模式)

  if (argc == 1) {
    // ---------------- [交互模式] ----------------
    // 创建一个 UI 会话 (Qt 界面或终端界面)
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);

    // 执行初始化宏文件 (通常用于设置默认视角、画出几何体)
    // 你的项目里应该有一个 init_vis.mac 文件
    UImanager->ApplyCommand("/control/execute init_vis.mac");

    // 开始事件循环 (程序停在这里，直到你点退出)
    ui->SessionStart();

    delete ui; // 退出后清理内存
  }
  else {
    // ---------------- [批处理模式] ----------------
    // 适合在服务器上跑通宵
    G4String command = "/control/execute ";
    G4String fileName = argv[1]; // 获取文件名 (例如 run.mac)
    
    // 执行这个宏文件
    UImanager->ApplyCommand(command + fileName);
  }

  // =================== 9. 清理现场 ===================
  // 程序结束，释放内存
  delete visManager;
  delete runManager;
  
  return 0;
}