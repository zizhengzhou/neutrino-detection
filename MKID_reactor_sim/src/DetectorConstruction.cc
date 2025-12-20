#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

// 材料/同位素
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"

// Region 和 Cuts (用于低能物理控制)
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh" // 【新增】必须引入这个头文件

// 可视化
#include "G4VisAttributes.hh"
#include "G4Colour.hh"


namespace MKID_reactor_sim 
{

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(nullptr)
{}

DetectorConstruction::~DetectorConstruction() = default;

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // ------------------------ NIST 管理器 ------------------------
  auto* nist = G4NistManager::Instance();
  G4bool checkOverlaps = true;

  // ================================================================
  // 0. 用同位素自然丰度定义关键材料 (保持你原有的代码不变)
  // ================================================================

  auto* elC = new G4Element("Carbon",   "C",  6., 12.011*g/mole);
  auto* elH = new G4Element("Hydrogen", "H",  1., 1.00794*g/mole);

  auto* isoB10 = new G4Isotope("B10",  5, 10, 10.0129*g/mole);
  auto* isoB11 = new G4Isotope("B11",  5, 11, 11.0093*g/mole);
  auto* elB = new G4Element("Boron", "B", 2);
  elB->AddIsotope(isoB10, 19.9*perCent);
  elB->AddIsotope(isoB11, 80.1*perCent);

  auto* isoPb204 = new G4Isotope("Pb204", 82, 204, 203.973*g/mole);
  auto* isoPb206 = new G4Isotope("Pb206", 82, 206, 205.974*g/mole);
  auto* isoPb207 = new G4Isotope("Pb207", 82, 207, 206.976*g/mole);
  auto* isoPb208 = new G4Isotope("Pb208", 82, 208, 207.977*g/mole);
  auto* elPb = new G4Element("Lead", "Pb", 4);
  elPb->AddIsotope(isoPb204,  1.4*perCent);
  elPb->AddIsotope(isoPb206, 24.1*perCent);
  elPb->AddIsotope(isoPb207, 22.1*perCent);
  elPb->AddIsotope(isoPb208, 52.4*perCent);

  auto* isoCu63 = new G4Isotope("Cu63", 29, 63, 62.9296*g/mole);
  auto* isoCu65 = new G4Isotope("Cu65", 29, 65, 64.9278*g/mole);
  auto* elCu = new G4Element("Copper", "Cu", 2);
  elCu->AddIsotope(isoCu63, 69.15*perCent);
  elCu->AddIsotope(isoCu65, 30.85*perCent);

  auto* elAl = new G4Element("Aluminum", "Al", 13., 26.9815*g/mole);

  auto* matPE = new G4Material("Polyethylene", 0.94*g/cm3, 2);
  matPE->AddElement(elC, 1);
  matPE->AddElement(elH, 2);

  auto* matBoratedPE = new G4Material("BoratedPolyethylene", 1.0*g/cm3, 2);
  matBoratedPE->AddMaterial(matPE, 95.*perCent);
  matBoratedPE->AddElement(elB,     5.*perCent);

  auto* matScint = new G4Material("PlasticScintillator", 1.032*g/cm3, 2);
  matScint->AddElement(elC, 9);
  matScint->AddElement(elH,10);

  auto* matPb = new G4Material("ShieldLead", 11.34*g/cm3, 1);
  matPb->AddElement(elPb, 1.0);

  auto* matCu = new G4Material("ShieldCopper", 8.96*g/cm3, 1);
  matCu->AddElement(elCu, 1.0);

  auto* matAl = new G4Material("SampleAl", 2.70*g/cm3, 1);
  matAl->AddElement(elAl, 1.0);

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* env_mat   = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* concrete_mat = nist->FindOrBuildMaterial("G4_CONCRETE");

  // ================================================================
  // 1. World + Shielding Geometry (保持不变)
  // ================================================================

  G4double roomSizeXY = 3.0 * m;
  G4double roomSizeZ  = 3.0 * m;
  G4double concreteThickness = 30.0 * cm;
  G4double worldMargin = 1.0 * m;
  G4double env_sizeXY = roomSizeXY;
  G4double env_sizeZ  = roomSizeZ;
  G4double world_sizeXY = env_sizeXY + 2.0 * concreteThickness + 2.0 * worldMargin;
  G4double world_sizeZ  = env_sizeZ  + 2.0 * concreteThickness + 2.0 * worldMargin;

  auto* solidWorld = new G4Box("World", 0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);
  auto* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  auto* physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false, 0, checkOverlaps);

  G4double concreteOuterXY = env_sizeXY + 2.0 * concreteThickness;
  G4double concreteOuterZ  = env_sizeZ  + 2.0 * concreteThickness;
  auto* solidConcrete = new G4Box("ConcreteShield", 0.5 * concreteOuterXY, 0.5 * concreteOuterXY, 0.5 * concreteOuterZ);
  auto* logicConcrete = new G4LogicalVolume(solidConcrete, concrete_mat, "ConcreteShield");
  new G4PVPlacement(nullptr, G4ThreeVector(), logicConcrete, "ConcreteShield", logicWorld, false, 0, checkOverlaps);

  auto* solidEnv = new G4Box("Envelope", 0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);
  auto* logicEnv = new G4LogicalVolume(solidEnv, env_mat, "Envelope");
  new G4PVPlacement(nullptr, G4ThreeVector(), logicEnv, "Envelope", logicConcrete, false, 0, checkOverlaps);

  // ================================================================
  // 2. 样品盒 + 多层屏蔽 (保持不变)
  // ================================================================

  G4double sampleSizeX = 20.0 * cm; G4double sampleSizeY = 20.0 * cm; G4double sampleSizeZ = 2.0  * cm;
  G4double hx = 0.5 * sampleSizeX; G4double hy = 0.5 * sampleSizeY; G4double hz = 0.5 * sampleSizeZ;

  G4double tPb = 5.0 * cm; G4double tPE = 5.0 * cm; G4double tSc = 5.0 * cm; G4double tCu = 1.0 * cm;
  G4ThreeVector detCenter(0., 0., 0.);

  // 计算各层尺寸 (代码保持不变，为了阅读简洁省略部分重复计算逻辑，直接用你原本的)
  G4double hx_sample = hx; G4double hy_sample = hy; G4double hz_sample = hz;
  G4double hx_Pb1 = hx_sample + tPb; G4double hy_Pb1 = hy_sample + tPb; G4double hz_Pb1 = hz_sample + tPb;
  G4double hx_PE1 = hx_Pb1 + tPE; G4double hy_PE1 = hy_Pb1 + tPE; G4double hz_PE1 = hz_Pb1 + tPE;
  G4double hx_Pb2 = hx_PE1 + tPb; G4double hy_Pb2 = hy_PE1 + tPb; G4double hz_Pb2 = hz_PE1 + tPb;
  G4double hx_PE2 = hx_Pb2 + tPE; G4double hy_PE2 = hy_Pb2 + tPE; G4double hz_PE2 = hz_Pb2 + tPE;
  G4double hx_Pb3 = hx_PE2 + tPb; G4double hy_Pb3 = hy_PE2 + tPb; G4double hz_Pb3 = hz_PE2 + tPb;
  G4double hx_Pb4 = hx_Pb3 + tPb; G4double hy_Pb4 = hy_Pb3 + tPb; G4double hz_Pb4 = hz_Pb3 + tPb;
  G4double hx_Sc  = hx_Pb4 + tSc; G4double hy_Sc  = hy_Pb4 + tSc; G4double hz_Sc  = hz_Pb4 + tSc;
  G4double hx_Pb5 = hx_Sc + tPb; G4double hy_Pb5 = hy_Sc + tPb; G4double hz_Pb5 = hz_Sc + tPb;
  G4double hx_Cu1 = hx_Pb5 + tCu; G4double hy_Cu1 = hy_Pb5 + tCu; G4double hz_Cu1 = hz_Pb5 + tCu;
  G4double hx_Cu2 = hx_Cu1 + tCu; G4double hy_Cu2 = hy_Cu1 + tCu; G4double hz_Cu2 = hz_Cu1 + tCu;

  G4LogicalVolume* mother = logicEnv;

  // 构建各层 (代码保持不变)
  auto* solidCu2  = new G4Box("Cu2", hx_Cu2, hy_Cu2, hz_Cu2);
  auto* logicCu2  = new G4LogicalVolume(solidCu2, matCu, "Cu2");
  new G4PVPlacement(nullptr, detCenter, logicCu2, "Cu2", mother,false,0,checkOverlaps);
  mother = logicCu2;

  auto* solidCu1  = new G4Box("Cu1", hx_Cu1, hy_Cu1, hz_Cu1);
  auto* logicCu1  = new G4LogicalVolume(solidCu1, matCu, "Cu1");
  new G4PVPlacement(nullptr, detCenter, logicCu1, "Cu1", mother,false,0,checkOverlaps);
  mother = logicCu1;

  auto* solidPb5  = new G4Box("Pb5", hx_Pb5, hy_Pb5, hz_Pb5);
  auto* logicPb5  = new G4LogicalVolume(solidPb5, matPb, "Pb5");
  new G4PVPlacement(nullptr, detCenter, logicPb5, "Pb5", mother,false,0,checkOverlaps);
  mother = logicPb5;

  auto* solidScint = new G4Box("PlasticScintillator", hx_Sc, hy_Sc, hz_Sc);
  auto* logicScint = new G4LogicalVolume(solidScint, matScint, "PlasticScintillator");
  new G4PVPlacement(nullptr, detCenter, logicScint, "PlasticScintillator", mother,false,0,checkOverlaps);
  mother = logicScint;

  auto* solidPb4  = new G4Box("Pb4", hx_Pb4, hy_Pb4, hz_Pb4);
  auto* logicPb4  = new G4LogicalVolume(solidPb4, matPb, "Pb4");
  new G4PVPlacement(nullptr, detCenter, logicPb4, "Pb4", mother,false,0,checkOverlaps);
  mother = logicPb4;

  auto* solidPb3  = new G4Box("Pb3", hx_Pb3, hy_Pb3, hz_Pb3);
  auto* logicPb3  = new G4LogicalVolume(solidPb3, matPb, "Pb3");
  new G4PVPlacement(nullptr, detCenter, logicPb3, "Pb3", mother,false,0,checkOverlaps);
  mother = logicPb3;

  auto* solidPE2  = new G4Box("BoratedPE2", hx_PE2, hy_PE2, hz_PE2);
  auto* logicPE2  = new G4LogicalVolume(solidPE2, matBoratedPE, "BoratedPE2");
  new G4PVPlacement(nullptr, detCenter, logicPE2, "BoratedPE2", mother,false,0,checkOverlaps);
  mother = logicPE2;

  auto* solidPb2  = new G4Box("Pb2", hx_Pb2, hy_Pb2, hz_Pb2);
  auto* logicPb2  = new G4LogicalVolume(solidPb2, matPb, "Pb2");
  new G4PVPlacement(nullptr, detCenter, logicPb2, "Pb2", mother,false,0,checkOverlaps);
  mother = logicPb2;

  auto* solidPE1  = new G4Box("BoratedPE1", hx_PE1, hy_PE1, hz_PE1);
  auto* logicPE1  = new G4LogicalVolume(solidPE1, matBoratedPE, "BoratedPE1");
  new G4PVPlacement(nullptr, detCenter, logicPE1, "BoratedPE1", mother,false,0,checkOverlaps);
  mother = logicPE1;

  auto* solidPb1  = new G4Box("Pb1", hx_Pb1, hy_Pb1, hz_Pb1);
  auto* logicPb1  = new G4LogicalVolume(solidPb1, matPb, "Pb1");
  new G4PVPlacement(nullptr, detCenter, logicPb1, "Pb1", mother,false,0,checkOverlaps);
  mother = logicPb1;

  // ================================================================
  // 【修改点 2】样品盒逻辑：修正 PV 名称并添加 Production Cuts
  // ================================================================

  auto* solidSample = new G4Box("SampleRegion", hx_sample, hy_sample, hz_sample);
  auto* logicSample = new G4LogicalVolume(solidSample, matAl, "SampleRegion");

  new G4PVPlacement(nullptr,
                    detCenter,
                    logicSample,
                    "SampleBoxPV", // <--- 【关键】这里改成了 SampleBoxPV，与 SteppingAction 里的判断一致
                    mother,
                    false,
                    0,
                    checkOverlaps);

  fScoringVolume = logicSample;

  // ----------------------------------------------------------------
  // 【修改点 3】设置 MKID 区域和超低阈值 (Cuts)
  // ----------------------------------------------------------------
  // 我们不需要每次都 new，先检查 RegionStore 里是否已经有了
  G4String regName = "MKID_Region";
  G4Region* sampleRegion = G4RegionStore::GetInstance()->GetRegion(regName, false);
  
  if (!sampleRegion) {
      sampleRegion = new G4Region(regName);
      sampleRegion->AddRootLogicalVolume(logicSample); // 将样品盒逻辑体加入该区域

      // 定义极小的生产阈值 (例如 250 nm)
      // 在 Al 中，250nm 对应的电子能量极低，远小于 1keV，这允许 Penelope 模拟 10-1000 eV 信号
      G4ProductionCuts* cuts = new G4ProductionCuts();
      G4double mkidCutValue = 10.0 * nm;

      cuts->SetProductionCut(mkidCutValue, "gamma");
      cuts->SetProductionCut(mkidCutValue, "e-");
      cuts->SetProductionCut(mkidCutValue, "e+");
      cuts->SetProductionCut(mkidCutValue, "neutron");

      // 将 Cuts 应用到该 Region
      sampleRegion->SetProductionCuts(cuts);
  }

  // ================================================================
  // 3. 可视化属性 (保持不变)
  // ================================================================
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto visConcrete = new G4VisAttributes(G4Colour(0.6, 0.6, 0.6, 0.05));
  visConcrete->SetForceSolid(true);
  logicConcrete->SetVisAttributes(visConcrete);

  auto visEnv = new G4VisAttributes(G4Colour(0.3, 0.3, 1.0, 0.10));
  visEnv->SetForceSolid(true);
  logicEnv->SetVisAttributes(visEnv);

  auto visCu2 = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0, 0.25));
  visCu2->SetForceSolid(true);
  logicCu2->SetVisAttributes(visCu2);
  // ... (其余可视化代码省略，保持原样即可) ...
  auto visSample = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.50));
  visSample->SetForceSolid(true);
  logicSample->SetVisAttributes(visSample);

  return physWorld;
}
}

// } // namespace MKID_reactor_sim  <--- 【注意】这里注释掉了 namespace 结束符