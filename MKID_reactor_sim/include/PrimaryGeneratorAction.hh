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
/// \file MKID_reactor_sim/include/PrimaryGeneratorAction.hh
/// \brief Definition of the MKID_reactor_sim::PrimaryGeneratorAction class

#ifndef MKID_reactor_simPrimaryGeneratorAction_h
#define MKID_reactor_simPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"

class G4Event;


namespace MKID_reactor_sim
{

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued
/// in front of the phantom across 80% of the (X,Y) phantom size.

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();    
    virtual ~PrimaryGeneratorAction();

    // 这是 Geant4 核心会在每个 Event 开始时调用的函数
    virtual void GeneratePrimaries(G4Event*);         

    // 提供一个接口，万一其他类想获取 GPS 指针（通常不需要，但保留是个好习惯）
    const G4GeneralParticleSource* GetParticleGun() const { return fParticleGun; }
  
  private:
    // 【核心变量】
    // 注意：变量名我依然保留了 fParticleGun 以符合习惯，
    // 但它的类型已经是 G4GeneralParticleSource 了！
    G4GeneralParticleSource* fParticleGun;
};

}  // namespace MKID_reactor_sim

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
