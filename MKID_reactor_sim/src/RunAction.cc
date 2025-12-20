//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file MKID_reactor_sim/src/RunAction.cc
/// \brief Implementation of the MKID_reactor_sim::RunAction class

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4AnalysisManager.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

namespace MKID_reactor_sim
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
// 1. 获取分析管理器单例
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root"); // 设置输出格式为 ROOT
  analysisManager->SetNtupleMerging(true);     // 【关键】开启多线程 Ntuple 合并

  // 2. 定义 Ntuple (数据表结构)
  // 就像在 Excel 里定义表头：
  analysisManager->CreateNtuple("FluxData", "Flux entering Sample Box");
  analysisManager->CreateNtupleIColumn("PDG");    // 列 0: 粒子ID
  analysisManager->CreateNtupleDColumn("Energy"); // 列 1: 能量
  analysisManager->CreateNtupleDColumn("Time");   // 列 2: 时间
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

// 每次 Run 开始前，打开一个新的数据文件
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile("MKID_Background_Sim");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
// 1. 保存并关闭文件
  // Geant4 会自动处理多线程数据的写入和合并
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  
  // 我们不需要打印任何东西，数据都在 ROOT 文件里，回去用 Python/ROOT 画图即可。
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace MKID_reactor_sim
