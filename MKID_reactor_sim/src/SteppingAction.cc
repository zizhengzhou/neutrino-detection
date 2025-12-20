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
/// \file MKID_reactor_sim/src/SteppingAction.cc
/// \brief Implementation of the MKID_reactor_sim::SteppingAction class

#include "SteppingAction.hh"
#include "G4AnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"

namespace MKID_reactor_sim
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction() : G4UserSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // 1. 获取 前一步(Pre) 和 后一步(Post) 的点
  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  G4StepPoint* postStepPoint = step->GetPostStepPoint();

  // 安全检查：如果粒子飞出世界边界，PostStep 可能是空的
  if (!postStepPoint->GetTouchableHandle()->GetVolume()) return;

  // 2. 获取物理体积的名字 (Physical Volume Name)
  // 我们比较名字字符串，虽然比比较指针稍微慢一点点，但逻辑最清楚
  G4String preVolName = preStepPoint->GetTouchableHandle()->GetVolume()->GetName();
  G4String postVolName = postStepPoint->GetTouchableHandle()->GetVolume()->GetName();

  // 3. 定义你的样品盒名字 (必须和 DetectorConstruction 一致)
  G4String targetName = "SampleBoxPV"; 

  // 4. 【核心逻辑：边界跨越】
  // 我们不关心粒子在里面怎么乱窜，
  // 我们只关心：上一刻在外面，下一刻在里面的那一瞬间（入射）！
  if (preVolName != targetName && postVolName == targetName) 
  {
      // 5. 【数据获取】
      // 获取 AnalysisManager 单例 (不需要指针传递)
      auto analysisManager = G4AnalysisManager::Instance();

      // 获取粒子类型 (PDG)
      G4int pdg = step->GetTrack()->GetDefinition()->GetPDGEncoding();
      
      // 获取入射动能 (注意：是 PreStep 的能量，即入射前的能量)
      // 如果你用 B1 的 GetTotalEnergyDeposit，那是能量损失，不是粒子能量！
      G4double kinEnergy = preStepPoint->GetKineticEnergy();
      
      // 获取全局时间
      G4double globalTime = preStepPoint->GetGlobalTime();

      // 6. 【数据写入】直接填入 ROOT 文件的 Ntuple
      analysisManager->FillNtupleIColumn(0, pdg);
      analysisManager->FillNtupleDColumn(1, kinEnergy);
      analysisManager->FillNtupleDColumn(2, globalTime);
      analysisManager->AddNtupleRow();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace MKID_reactor_sim
