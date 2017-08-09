#include "MedWithFOAM.h"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "MEDFileMesh.hxx"
#include "MEDCoupling.hxx"

MedWithFOAM::MedWithFOAM() {}
MedWithFOAM::MedWithFOAM (
  ParaMEDMEM::MEDCouplingUMesh *FoamMesh,
  const ParaMEDMEM::MEDCouplingUMesh *FemusMesh
) {
  SetMeshes (FoamMesh, FemusMesh);
}
// MedWithFOAM::~MedWithFOAM(){}
MedWithFOAM::~MedWithFOAM() {}
void MedWithFOAM::Clear() {
  // CLEARING MESHES
  _FoamMesh  = NULL;
  _FemusMesh = NULL;
  _FoamFMesh = NULL;
  _OFLocal3D  = NULL;
  _OFLocal2D  = NULL;

  // CLEARING COUPLING MAPS
  _Global2Dto3D = NULL;
  _Global3Dto2D = NULL;
  _ParToMedMap = NULL;
  _LOFtoLMED3D = NULL;
  _LOFtoLMED2D = NULL;
  _LMED3DtoLOF = NULL;
  _LMED2DtoLOF = NULL;
  return;
}
void MedWithFOAM::SetMeshes (
  ParaMEDMEM::MEDCouplingUMesh *FoamMesh,
  const ParaMEDMEM::MEDCouplingUMesh *FemusMesh
) {
  _FoamMesh  = FoamMesh;
  _FemusMesh = FemusMesh;
  _areOFandFemusMeshSet = true;

  const int FemusCells = _FemusMesh->getNumberOfCells();
  const int FemusNodesPerCell = _FemusMesh->getNumberOfNodesInCell (0);
  const int FemusDim = _FemusMesh->getSpaceDimension();
  _Global2Dto3D = ParaMEDMEM::DataArrayDouble::New();
  _Global2Dto3D->alloc (FemusCells,1);
  _Global3Dto2D = ParaMEDMEM::DataArrayDouble::New();
  _Global3Dto2D->alloc (FemusCells,1);

  std::vector<int> FemusCellConn;
  std::vector<double> NodeCoord;
  double *FemusNodesCoordinates = new double[3*FemusCells];
  for (int i_cell=0; i_cell<FemusCells; i_cell++) {
    _FemusMesh->getNodeIdsOfCell (i_cell,FemusCellConn);
    for (int dim=0; dim<3; dim++) FemusNodesCoordinates[i_cell*3 + dim] = 0.;
    for (int i_cnode=0; i_cnode<FemusNodesPerCell; i_cnode++) {
      _FemusMesh->getCoordinatesOfNode (FemusCellConn[i_cnode],NodeCoord);
      FemusNodesCoordinates[i_cell*3 + 0] += NodeCoord[0]/FemusNodesPerCell;
      FemusNodesCoordinates[i_cell*3 + 1] += NodeCoord[1]/FemusNodesPerCell;
      if (FemusDim==3) FemusNodesCoordinates[i_cell*3 + 2] += NodeCoord[2]/FemusNodesPerCell;
      NodeCoord.clear();
    }// Mid point of cell i_cell
    FemusCellConn.clear();
  }// now FemusNodesCoordinates contains all the cell midpoints
  ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayInt> elts;
  ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayInt> eltsIndex;
  _FoamMesh->getCellsContainingPoints (FemusNodesCoordinates,FemusCells,1.e-5,elts,eltsIndex);

  for (int i_cell=0; i_cell<FemusCells; i_cell++) {
    int NumPossibleCells = eltsIndex->getIJ (i_cell+1,0) - eltsIndex->getIJ (i_cell,0);
    int NewCell = elts->getIJ (eltsIndex->getIJ (i_cell,0) + 0,0);
    _Global2Dto3D->setIJ (i_cell,0, NewCell);
    _Global3Dto2D->setIJ (NewCell,0, i_cell);
  }// Coupling map is now completed

//   {//print the map
//    ParaMEDMEM::MEDCouplingFieldDouble* MAP=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS);
//    MAP->setMesh(_FemusMesh);
//    MAP->setArray(_CouplingMap);
//    MAP->setName("Map");
//    MEDLoader::WriteField("RESU_MED/CouplingMap.med",MAP,true);
//    MAP->decrRef();
//   }


  delete [] FemusNodesCoordinates;
  return;
}


void MedWithFOAM::SetParallelMap (
  Foam::fvMesh *FoamFMesh,
  ParaMEDMEM::MEDCouplingUMesh *Foam2D,
  const int proc
) {
  printf ("proc: %d Setting coupling maps\n",proc);

  _FoamFMesh = FoamFMesh;
  const int ParFoamCells = _FoamFMesh->nCells();
  const int FoamNodesPerCell = _FoamFMesh->cellPoints() [0].size();

  _ParToMedMap = ParaMEDMEM::DataArrayDouble::New();
  _LOFtoLMED3D = ParaMEDMEM::DataArrayDouble::New();
  _LOFtoLMED2D = ParaMEDMEM::DataArrayDouble::New();
  _LMED3DtoLOF = ParaMEDMEM::DataArrayDouble::New();
  _LMED2DtoLOF = ParaMEDMEM::DataArrayDouble::New();
  _ParToMedMap->alloc (ParFoamCells,1);
  _LOFtoLMED3D->alloc (ParFoamCells,1);
  _LOFtoLMED2D->alloc (ParFoamCells,1);
  _LMED3DtoLOF->alloc (ParFoamCells,1);
  _LMED2DtoLOF->alloc (ParFoamCells,1);

  std::vector<int> FoamCellConn;
  std::vector<double> NodeCoord;
  double *FoamNodesCoordinates3D = new double[3*ParFoamCells]; // Array with cell midpoints
  double *FoamNodesCoordinates2D = new double[2*ParFoamCells]; // Array with cell midpoints -> cells contracted on z=0 plane

  for (int i=0; i<ParFoamCells; i++) {
    for (int j=0; j<3; j++) {
      double coord = _FoamFMesh->C().internalField() [i][j];
      FoamNodesCoordinates3D[i*3 + j] = coord;
    }
    for (int j=0; j<2; j++) {
      double coord = _FoamFMesh->C().internalField() [i][j];
      FoamNodesCoordinates2D[i*2 + j] = coord;
    }
  }

  printf ("proc: %d Building local to global 3D med map\n",proc);

  ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayInt> elts2;
  ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayInt> eltsIndex2;
  _FoamMesh->getCellsContainingPoints (FoamNodesCoordinates3D,ParFoamCells,1.e-5,elts2,eltsIndex2);
  int *CellArray3D = new int[ParFoamCells];// Array containing the global 3d med mesh cells constituting the local mesh
  for (int i_cell=0; i_cell<ParFoamCells; i_cell++) {   // BUILDING MAP - FROM LOCAL FOAM TO GLOBAL 3D MED
    int NumPossibleCells = eltsIndex2->getIJ (i_cell+1,0) - eltsIndex2->getIJ (i_cell,0);
    int NewCell = elts2->getIJ (eltsIndex2->getIJ (i_cell,0) + 0,0);
    _ParToMedMap->setIJ (i_cell,0, NewCell);
    CellArray3D[i_cell] = NewCell; //
  }// Coupling map is now completed

  // BUILDING THE 3D MED LOCAL MESH - _OFLocal3D ======================================

  printf ("proc: %d Building local 3d mesh \n",proc);
  _OFLocal3D = _FoamMesh->buildPartOfMySelf (CellArray3D, CellArray3D + ParFoamCells);
  _OFLocal3D->setName ("ParMesh_"+std::to_string (proc));
//   MEDLoader::WriteUMesh ("ParMesh3D_"+std::to_string (proc) +".med", _OFLocal3D,false);
//   _OFLocal3D->sortCellsInMEDFileFrmt();

  printf ("proc: %d Building local to global 2D med map\n",proc);
  _FemusMesh->getCellsContainingPoints (FoamNodesCoordinates2D,ParFoamCells,1.e-5,elts2,eltsIndex2);
  int *CellArray2D = new int[ParFoamCells];// Array containing the global 2d med mesh cells constituting the local mesh
  for (int i_cell=0; i_cell<ParFoamCells; i_cell++) {   // BUILDING MAP - FROM LOCAL FOAM TO GLOBAL 2D MED
    int NumPossibleCells = eltsIndex2->getIJ (i_cell+1,0) - eltsIndex2->getIJ (i_cell,0);
    int NewCell = elts2->getIJ (eltsIndex2->getIJ (i_cell,0) + 0,0);
//         _ParToMedMap->setIJ ( i_cell,0, NewCell );
    CellArray2D[i_cell] = NewCell; //
  }// Coupling map is now completed

  printf ("proc: %d Building local 2d mesh \n",proc);
  _OFLocal2D = _FemusMesh->buildPartOfMySelf (CellArray2D, CellArray2D + ParFoamCells);
  _OFLocal2D->setName ("2DInterface_"+std::to_string (proc));
//   MEDLoader::WriteUMesh ("ParMesh2D_"+std::to_string (proc) +".med", _OFLocal2D,false);
  //   _OFLocal2D->sortCellsInMEDFileFrmt();

  printf ("proc: %d Building local to local 3D med map\n",proc);
  _OFLocal3D->getCellsContainingPoints (FoamNodesCoordinates3D,ParFoamCells,1.e-5,elts2,eltsIndex2);
  for (int i_cell=0; i_cell<ParFoamCells; i_cell++) {   // BUILDING MAP - FROM LOCAL FOAM TO LOCAL 3D MED
    int NumPossibleCells = eltsIndex2->getIJ (i_cell+1,0) - eltsIndex2->getIJ (i_cell,0);
    int NewCell = elts2->getIJ (eltsIndex2->getIJ (i_cell,0) + 0,0);
    _LOFtoLMED3D->setIJ (i_cell,0, NewCell);
    _LMED3DtoLOF->setIJ (NewCell,0, i_cell);
  }// Coupling map is now completed

  printf ("proc: %d Building local to local 2D med map\n",proc);
  _OFLocal2D->getCellsContainingPoints (FoamNodesCoordinates2D,ParFoamCells,1.e-5,elts2,eltsIndex2);
  for (int i_cell=0; i_cell<ParFoamCells; i_cell++) {   // BUILDING MAP - FROM LOCAL FOAM TO LOCAL 2D MED
    int NumPossibleCells = eltsIndex2->getIJ (i_cell+1,0) - eltsIndex2->getIJ (i_cell,0);
    int NewCell = elts2->getIJ (eltsIndex2->getIJ (i_cell,0) + 0,0);
    _LOFtoLMED2D->setIJ (i_cell,0, NewCell);
    _LMED2DtoLOF->setIJ (NewCell,0, i_cell);
  }// Coupling map is now completed

  delete [] FoamNodesCoordinates3D;
  delete [] FoamNodesCoordinates2D;
  delete [] CellArray3D;
  delete [] CellArray2D;
  return;
}

// ParaMEDMEM::MEDCouplingUMesh* MedWithFOAM::GetParMesh ( const int proc ) {
//     const std::string FileName = "RESU_MED/Par3DMedMesh"+std::to_string ( proc ) +".med";
//     std::vector<std::string> MeshNames = MEDLoader::GetMeshNames ( FileName );
//     std::vector<std::string> FamiliNames = MEDLoader::GetMeshFamiliesNames ( FileName,MeshNames[0] );
//
//     ParaMEDMEM::MEDFileUMesh * meshMEDFile = ParaMEDMEM::MEDFileUMesh::New ( FileName );
//     std::vector<int> IDS = meshMEDFile->getFamiliesIds ( FamiliNames );
//
//     for ( int i=0; i<IDS.size(); i++ ) printf ( "\n\n Proc: %d  Family Id: %d \n\n",proc, IDS[i] );
//     for ( int i=0; i<FamiliNames.size(); i++ ) printf ( "\n\n Proc: %d  Family name: %s \n\n",proc, FamiliNames[i].c_str() );
//
//     const std::string MeshName = MeshNames[0];
// //   const std::string MeshName = MeshNames[0];
// //   const std::vector<std::string> Families = {"Fam"+std::to_string(proc)};
//     const std::vector<std::string> Families = {"Family_"+std::to_string ( - ( proc+1 ) ) };
//     const std::vector<std::string> Groups = {"group_"+std::to_string ( ( proc+1 ) ) };
//
// //   Families.resize(1);
// //   Families[0] = std::to_string(proc);
// //
// //   ParaMEDMEM::MEDCouplingUMesh * ParMesh = MEDLoader::ReadUMeshFromFamilies(FileName, MeshName, 0, Families);
//     ParaMEDMEM::MEDCouplingUMesh * ParMesh = MEDLoader::ReadUMeshFromGroups ( FileName, MeshName, 0, Groups );
//     MEDLoader::WriteUMesh ( "ParallelMesh"+std::to_string ( proc ) +".med",ParMesh, true );
//     return ParMesh;
//     }

// ParaMEDMEM::MEDCouplingUMesh* MedWithFOAM::Get2DParMesh ( const int proc ) {
//     const std::string FileName = "RESU_MED/Par2DMedMesh"+std::to_string ( proc ) +".med";
//     std::vector<std::string> MeshNames = MEDLoader::GetMeshNames ( FileName );
//     ParaMEDMEM::MEDFileUMesh * meshMEDFile = ParaMEDMEM::MEDFileUMesh::New ( FileName );
//
//     const std::string MeshName = MeshNames[0];
//     const std::vector<std::string> Groups = {"group_"+std::to_string ( ( proc+1 ) ) };
//
//     ParaMEDMEM::MEDCouplingUMesh * ParMesh = MEDLoader::ReadUMeshFromGroups ( FileName, MeshName, 0, Groups );
//     MEDLoader::WriteUMesh ( "ParallelMesh"+std::to_string ( proc ) +".med",ParMesh, true );
// //   MEDLoader::
// //   const std::vector<std::string> Families = {"Family_"+std::to_string(-(proc+1))};
// //     std::vector<std::string> FamiliNames = MEDLoader::GetMeshFamiliesNames(FileName,MeshNames[0]);
// //   std::vector<int> IDS = meshMEDFile->getFamiliesIds(FamiliNames);
// //   for(int i=0; i<IDS.size(); i++) printf("\n\n Proc: %d  Family Id: %d \n\n",proc, IDS[i]);
// //   for(int i=0; i<FamiliNames.size(); i++) printf("\n\n Proc: %d  Family name: %s \n\n",proc, FamiliNames[i].c_str());
//     return ParMesh;
//     }

void MedWithFOAM::PrintMedMesh (Foam::fvMesh *FoamMesh, int proc) {

  std::ostringstream stm ;
  stm << proc;
  std::string Proc = stm.str();

  int OFToMedLocalCon[] = {0,1,3,2,4,5,7,6};
  const int Dim = 3;
  const int NCells = FoamMesh->nCells();    // number of cells
  const int NNodes = FoamMesh->nPoints();   // number of nodes
  int ElemConn[8];                     // CellConn

  ParaMEDMEM::MEDCouplingUMesh *MedMesh=ParaMEDMEM::MEDCouplingUMesh::New ("My2DMesh"+Proc,3);
  MedMesh->allocateCells (NCells);
  ParaMEDMEM::DataArrayDouble *coordsArr=ParaMEDMEM::DataArrayDouble::New();
  coordsArr->alloc (NNodes,3);

  for (int i=0; i<NCells; i++) {   // Loop over cells -> adding cells to med mesh ------------
    for (int j=0; j<8; j++) {   // Loop over cell nodes -> reading coordinates ---------------
      ElemConn[j] = FoamMesh->cellPoints() [i][OFToMedLocalCon[j]];
      for (int dim=0; dim<Dim; dim++)
        coordsArr->setIJ (ElemConn[j],dim,FoamMesh->points() [ElemConn[j]][dim]);   // storing nodes coordinates inside coordsArr med array
    }// end loop over cell nodes --------------------------------------------------------
    MedMesh->insertNextCell (INTERP_KERNEL::NORM_HEXA8,8,ElemConn);
  }// end loop over cells ---------------------------------------------------------------

  MedMesh->finishInsertingCells(); // med mesh is complete
  MedMesh->setCoords (coordsArr);   // adding coordinates to mesh nodes
  MEDLoader::WriteUMesh ("MedMeshOF"+Proc+".med",MedMesh,true);

  coordsArr->decrRef();
  MedMesh->decrRef();

  return;
}

void MedWithFOAM::CreateMedBoundaryMesh (Foam::fvMesh *FoamMesh, std::string patch, int proc) {

  std::ostringstream stm ;
  stm << proc;
  std::string Proc = stm.str();

  Foam::label MovingWallPatchId = FoamMesh->boundaryMesh().findPatchID (patch);   // -> label (ovvero numero della patch)
  const Foam::polyPatch &cPatch = FoamMesh->boundaryMesh() [MovingWallPatchId];          // info della patch: tipo, numero facce, faccia iniziale

  const int NumOfFaces  = cPatch.size();
  if (NumOfFaces==0) {
    std::cout<<"patch: "<<patch<<"\t proc "<<proc<<"\t 0 faces, now exiting \n";
    return;
  }
  const int NumOfPoints = cPatch.nPoints();
  std::map <int, int> GlobalToLocal;
  Foam::labelList GlobalNodesNum = cPatch.meshPoints();
  Foam::labelList LocalNodesNum  = cPatch.boundaryPoints();
  forAll (GlobalNodesNum, j) GlobalToLocal[GlobalNodesNum[j]] = LocalNodesNum[j];

  ParaMEDMEM::MEDCouplingUMesh *BoundMedMesh=ParaMEDMEM::MEDCouplingUMesh::New (patch+"_"+Proc,2);
  BoundMedMesh->allocateCells (NumOfFaces);
  ParaMEDMEM::DataArrayDouble *coordsArr=ParaMEDMEM::DataArrayDouble::New();
  coordsArr->alloc (NumOfPoints,3);

  Foam::faceList FACCE = FoamMesh->faces();         // list with faces ids, nodes per face and node ids
  Foam::labelList PatchFaces  = cPatch.faceCells(); // list of patch faces, local numbering: 0 -> cPatch.size()
  int FaceConn[4];
  forAll (PatchFaces, index) {
    int FaceId = index + cPatch.start();
    forAll (FACCE[FaceId], l) {
      FaceConn[l] = GlobalToLocal[FACCE[FaceId][l]];
      for (int dim=0; dim<3; dim++) coordsArr->setIJ (FaceConn[l],dim,FoamMesh->points() [FACCE[FaceId][l]][dim]);
    }
    BoundMedMesh->insertNextCell (INTERP_KERNEL::NORM_QUAD4,4,FaceConn);
  }

  BoundMedMesh->finishInsertingCells(); // med mesh is complete
  BoundMedMesh->setCoords (coordsArr);   // adding coordinates to mesh nodes
  MEDLoader::WriteUMesh ("Boundary_"+patch+"_"+Proc+".med",BoundMedMesh,true);
  coordsArr->decrRef();
  BoundMedMesh->decrRef();
  return;
}

ParaMEDMEM::MEDCouplingUMesh *MedWithFOAM::GetMedMesh (Foam::fvMesh *FoamMesh, int proc) {

  std::ostringstream stm ;
  stm << proc;
  std::string Proc = stm.str();

  int OFToMedLocalCon[] = {0,1,3,2,4,5,7,6};
//   int OFToMedLocalCon[] ={0,1,3,2,5,4,7,6};
//   int OFToMedLocalCon[] ={0,1,2,3,4,5,6,7};
  const int Dim = 3;
  const int NCells = FoamMesh->nCells();    // number of cells
  const int NNodes = FoamMesh->nPoints();   // number of nodes
  int ElemConn[8];                     // CellConn

  ParaMEDMEM::MEDCouplingUMesh *MedMesh=ParaMEDMEM::MEDCouplingUMesh::New ("My2DMesh"+Proc,3);
  MedMesh->allocateCells (NCells);
  ParaMEDMEM::DataArrayDouble *coordsArr=ParaMEDMEM::DataArrayDouble::New();
  coordsArr->alloc (NNodes,3);

  for (int i=0; i<NCells; i++) {   // Loop over cells -> adding cells to med mesh ------------
    for (int j=0; j<8; j++) {   // Loop over cell nodes -> reading coordinates ---------------
      ElemConn[j] = FoamMesh->cellPoints() [i][OFToMedLocalCon[j]];
      if (i==8) Foam::Info<<FoamMesh->points() [FoamMesh->cellPoints() [i][j]] <<"\n";
      for (int dim=0; dim<Dim; dim++)
        coordsArr->setIJ (ElemConn[j],dim,FoamMesh->points() [ElemConn[j]][dim]);   // storing nodes coordinates inside coordsArr med array
    }// end loop over cell nodes --------------------------------------------------------
    MedMesh->insertNextCell (INTERP_KERNEL::NORM_HEXA8,8,ElemConn);
  }// end loop over cells ---------------------------------------------------------------

  MedMesh->finishInsertingCells(); // med mesh is complete
  MedMesh->setCoords (coordsArr);   // adding coordinates to mesh nodes
  coordsArr->decrRef();

  return MedMesh;
}


ParaMEDMEM::MEDCouplingUMesh *MedWithFOAM::GetMedBoundaryMesh (Foam::fvMesh *FoamMesh, std::string patch, int proc) {

  std::ostringstream stm ;
  stm << proc;
  std::string Proc = stm.str();

  Foam::label MovingWallPatchId = FoamMesh->boundaryMesh().findPatchID (patch);   // -> label (ovvero numero della patch)
  const Foam::polyPatch &cPatch = FoamMesh->boundaryMesh() [MovingWallPatchId];          // info della patch: tipo, numero facce, faccia iniziale

  const int NumOfFaces  = cPatch.size();
  if (NumOfFaces==0) {
    std::cout<<"No faces for patch "<<patch<<" on processor "<<proc<<"; exiting from function GetMedBoundaryMesh \n";
    return NULL;
  }
  const int NumOfPoints = cPatch.nPoints();
  std::map <int, int> GlobalToLocal;
  Foam::labelList GlobalNodesNum = cPatch.meshPoints();
  Foam::labelList LocalNodesNum  = cPatch.boundaryPoints();
  forAll (GlobalNodesNum, j) GlobalToLocal[GlobalNodesNum[j]] = LocalNodesNum[j];

  ParaMEDMEM::MEDCouplingUMesh *BoundMedMesh=ParaMEDMEM::MEDCouplingUMesh::New (patch+"_"+Proc,2);
  BoundMedMesh->allocateCells (NumOfFaces);
  ParaMEDMEM::DataArrayDouble *coordsArr=ParaMEDMEM::DataArrayDouble::New();
  coordsArr->alloc (NumOfPoints,3);

  Foam::faceList FACCE = FoamMesh->faces();         // list with faces ids, nodes per face and node ids
  Foam::labelList PatchFaces  = cPatch.faceCells(); // list of patch faces, local numbering: 0 -> cPatch.size()
  int FaceConn[4];
  forAll (PatchFaces, index) {
    int FaceId = index + cPatch.start();
    forAll (FACCE[FaceId], l) {
      FaceConn[l] = GlobalToLocal[FACCE[FaceId][l]];
      for (int dim=0; dim<3; dim++) coordsArr->setIJ (FaceConn[l],dim,FoamMesh->points() [FACCE[FaceId][l]][dim]);
    }
    BoundMedMesh->insertNextCell (INTERP_KERNEL::NORM_QUAD4,4,FaceConn);
  }

  BoundMedMesh->finishInsertingCells(); // med mesh is complete
  BoundMedMesh->setCoords (coordsArr);   // adding coordinates to mesh nodes

  coordsArr->decrRef();

  return BoundMedMesh;
}

void MedWithFOAM::PrintMedSolution (Foam::volVectorField *Field,
                                    std::string FieldName,
                                    int proc,
                                    ParaMEDMEM::MEDCouplingUMesh *MedMesh,
                                    bool WriteFromScratch
                                   ) {
  std::ostringstream stm ;
  stm << proc;
  std::string Proc = stm.str();

  std::string directions[] = {"x","y","z"};

  // Save last solution in .med format - a file for each processor
  const int NCells = MedMesh->getNumberOfCells();    // number of cells
  ParaMEDMEM::MEDCouplingFieldDouble *f = ParaMEDMEM::MEDCouplingFieldDouble::New (ParaMEDMEM::ON_CELLS);
  f->setMesh (MedMesh);

  ParaMEDMEM::DataArrayDouble *Vel=ParaMEDMEM::DataArrayDouble::New();
  Vel->alloc (NCells,3);
  for (int i=0; i<NCells; i++)
    for (int dim=0; dim<MedMesh->getMeshDimension(); dim ++)
      Vel->setIJ (i,dim,Field->internalField() [i][dim]) ;

  for (int dim=0; dim<MedMesh->getMeshDimension(); dim ++)
    Vel->setInfoOnComponent (dim,directions[dim]);

  f->setArray (Vel);
  f->setName (FieldName);
  MEDLoader::WriteField ("RESU_MED/Solution"+Proc+".med",f,WriteFromScratch);
  f->decrRef();

  return;
}

ParaMEDMEM::MEDCouplingFieldDouble *MedWithFOAM::GetMedSolution (Foam::volVectorField *Field,
    std::string FieldName,
    int proc,
    ParaMEDMEM::MEDCouplingUMesh *MedMesh
                                                                ) {
  std::ostringstream stm ;
  stm << proc;
  std::string Proc = stm.str();

  std::string directions[] = {"x","y","z"};

  // Save last solution in .med format - a file for each processor
  const int NCells = MedMesh->getNumberOfCells();    // number of cells
  ParaMEDMEM::MEDCouplingFieldDouble *f = ParaMEDMEM::MEDCouplingFieldDouble::New (ParaMEDMEM::ON_CELLS);
  f->setMesh (MedMesh);

  ParaMEDMEM::DataArrayDouble *Vel=ParaMEDMEM::DataArrayDouble::New();
  const int Dim = MedMesh->getMeshDimension();
  Vel->alloc (NCells,Dim);
  
  for (int i=0; i<NCells; i++){
   if(Dim==3) for (int dim=0; dim<MedMesh->getMeshDimension(); dim ++)
      Vel->setIJ (i,dim,Field->internalField() [_LMED3DtoLOF->getIJ (i,0)][dim]) ;
   if(Dim==2) for (int dim=0; dim<MedMesh->getMeshDimension(); dim ++)
      Vel->setIJ (i,dim,Field->internalField() [_LMED2DtoLOF->getIJ (i,0)][dim]) ; 
  }

  for (int dim=0; dim<MedMesh->getMeshDimension(); dim ++)
    Vel->setInfoOnComponent (dim,directions[dim]);

  f->setArray (Vel);
  f->setName (FieldName);
//   MEDLoader::WriteField("RESU_MED/Solution"+Proc+".med",f,WriteFromScratch);
  return f;
}


void MedWithFOAM::PrintMedSolution (Foam::volScalarField *Field,
                                    std::string FieldName,
                                    int proc,
                                    ParaMEDMEM::MEDCouplingUMesh *MedMesh,
                                    bool WriteFromScratch) {
  std::ostringstream stm ;
  stm << proc;
  std::string Proc = stm.str();

  // Save last solution in .med format - a file for each processor
  const int NCells = MedMesh->getNumberOfCells();    // number of cells
  ParaMEDMEM::MEDCouplingFieldDouble *f = ParaMEDMEM::MEDCouplingFieldDouble::New (ParaMEDMEM::ON_CELLS);
  f->setMesh (MedMesh);
  const int Dim = MedMesh->getMeshDimension();
  ParaMEDMEM::DataArrayDouble *Vel=ParaMEDMEM::DataArrayDouble::New();
  Vel->alloc (NCells,1);
  for (int i=0; i<NCells; i++)  Vel->setIJ (i,0,Field->internalField() [i]) ;
  f->setArray (Vel);
  f->setName (FieldName);
  MEDLoader::WriteField ("RESU_MED/Solution"+Proc+".med",f,WriteFromScratch);
  f->decrRef();

  return;
}

ParaMEDMEM::MEDCouplingFieldDouble *MedWithFOAM::GetMedSolution (Foam::volScalarField *Field,
    std::string FieldName,
    int proc,
    ParaMEDMEM::MEDCouplingUMesh *MedMesh) {
  std::ostringstream stm ;
  stm << proc;
  std::string Proc = stm.str();

  // Save last solution in .med format - a file for each processor
  const int NCells = MedMesh->getNumberOfCells();    // number of cells
  ParaMEDMEM::MEDCouplingFieldDouble *f = ParaMEDMEM::MEDCouplingFieldDouble::New (ParaMEDMEM::ON_CELLS);
  f->setMesh (MedMesh);

  const int Dim = MedMesh->getMeshDimension();
  
  ParaMEDMEM::DataArrayDouble *Vel=ParaMEDMEM::DataArrayDouble::New();
  Vel->alloc (NCells,1);
  if(Dim==3)  for (int i=0; i<NCells; i++)  Vel->setIJ (i,0,Field->internalField() [_LMED3DtoLOF->getIJ (i,0)]) ;
  if(Dim==2)  for (int i=0; i<NCells; i++)  Vel->setIJ (i,0,Field->internalField() [_LMED2DtoLOF->getIJ (i,0)]) ;  
  f->setArray (Vel);
  f->setName (FieldName);
//   MEDLoader::WriteField("RESU_MED/Solution"+Proc+".med",f,WriteFromScratch);
//   f->decrRef();

  return f;
}
void MedWithFOAM::PrintMedSolution (Foam::volScalarField *Field,
                                    std::string FieldName,
                                    int proc,
                                    Foam::fvMesh *OFMesh,
                                    std::string PatchName,
                                    bool WriteFromScratch
                                   ) {
  std::ostringstream stm ;
  stm << proc;
  std::string Proc = stm.str();
  std::string SolName = "Solution_"+PatchName+"_"+Proc+".med";
  // Save last solution in .med format - a file for each processor

  ParaMEDMEM::MEDCouplingUMesh *MedBoundMesh = GetMedBoundaryMesh (OFMesh, PatchName, proc);
  if (MedBoundMesh == NULL) {
    std::cout<<"No faces for patch "<<PatchName<<" on processor "<<proc<<"; exiting from function PrintMedSolution \n";
    return;
  }

  const int NCells = MedBoundMesh->getNumberOfCells();    // number of cells

  ParaMEDMEM::MEDCouplingFieldDouble *f = ParaMEDMEM::MEDCouplingFieldDouble::New (ParaMEDMEM::ON_CELLS);
  f->setMesh (MedBoundMesh);

  Foam::label PatchId = OFMesh->boundaryMesh().findPatchID (PatchName);   // -> label (ovvero numero della patch)
  const Foam::polyPatch &cPatch = OFMesh->boundaryMesh() [PatchId];          // info della patch: tipo, numero facce, faccia iniziale

  ParaMEDMEM::DataArrayDouble *Vel=ParaMEDMEM::DataArrayDouble::New();
  Vel->alloc (NCells,1);
  forAll (cPatch ,faceI)  Vel->setIJ (faceI,0,Field->boundaryField() [PatchId][faceI]) ;
  f->setArray (Vel);
  f->setName (FieldName);
  MEDLoader::WriteField (SolName,f,WriteFromScratch);
  f->decrRef();
  MedBoundMesh->decrRef();
  return;
}

ParaMEDMEM::MEDCouplingFieldDouble *MedWithFOAM::GetMedSolution (Foam::volScalarField *Field,
    std::string FieldName,
    int proc,
    Foam::fvMesh *OFMesh,
    std::string PatchName
                                                                ) {
  std::ostringstream stm ;
  stm << proc;
  std::string Proc = stm.str();
  std::string SolName = "Solution_"+PatchName+"_"+Proc+".med";
  // Save last solution in .med format - a file for each processor

  ParaMEDMEM::MEDCouplingUMesh *MedBoundMesh = GetMedBoundaryMesh (OFMesh, PatchName, proc);
  if (MedBoundMesh == NULL) {
    std::cout<<"No faces for patch "<<PatchName<<" on processor "<<proc<<"; exiting from function PrintMedSolution \n";
  }

  const int NCells = MedBoundMesh->getNumberOfCells();    // number of cells

  ParaMEDMEM::MEDCouplingFieldDouble *f = ParaMEDMEM::MEDCouplingFieldDouble::New (ParaMEDMEM::ON_CELLS);
  f->setMesh (MedBoundMesh);

  Foam::label PatchId = OFMesh->boundaryMesh().findPatchID (PatchName);   // -> label (ovvero numero della patch)
  const Foam::polyPatch &cPatch = OFMesh->boundaryMesh() [PatchId];          // info della patch: tipo, numero facce, faccia iniziale

  ParaMEDMEM::DataArrayDouble *Vel=ParaMEDMEM::DataArrayDouble::New();
  Vel->alloc (NCells,1);
  forAll (cPatch ,faceI)  Vel->setIJ (faceI,0,Field->boundaryField() [PatchId][faceI]) ;
  f->setArray (Vel);
  f->setName (FieldName);
//   MEDLoader::WriteField(SolName,f,WriteFromScratch);
//   f->decrRef();
//   MedBoundMesh->decrRef();
  return f;
}

void MedWithFOAM::PrintMedSolution (Foam::volVectorField *Field,
                                    std::string FieldName,
                                    int proc,
                                    Foam::fvMesh *OFMesh,
                                    std::string PatchName,
                                    bool WriteFromScratch
                                   ) {
  std::ostringstream stm ;
  stm << proc;
  std::string Proc = stm.str();
  std::string SolName = "Solution_"+PatchName+"_"+Proc+".med";
  // Save last solution in .med format - a file for each processor
  std::string directions[] = {"x","y","z"};

  ParaMEDMEM::MEDCouplingUMesh *MedBoundMesh = GetMedBoundaryMesh (OFMesh, PatchName, proc);
  if (MedBoundMesh == NULL) {
    std::cout<<"No faces for patch "<<PatchName<<" on processor "<<proc<<"; exiting from function PrintMedSolution \n";
    return;
  }

  const int NCells = MedBoundMesh->getNumberOfCells();    // number of cells

  ParaMEDMEM::MEDCouplingFieldDouble *f = ParaMEDMEM::MEDCouplingFieldDouble::New (ParaMEDMEM::ON_CELLS);
  f->setMesh (MedBoundMesh);

  Foam::label PatchId = OFMesh->boundaryMesh().findPatchID (PatchName);   // -> label (ovvero numero della patch)
  const Foam::polyPatch &cPatch = OFMesh->boundaryMesh() [PatchId];          // info della patch: tipo, numero facce, faccia iniziale

  ParaMEDMEM::DataArrayDouble *Vel=ParaMEDMEM::DataArrayDouble::New();
  Vel->alloc (NCells,3);
  forAll (cPatch ,faceI)
  for (int dim=0; dim<3; dim++)
    Vel->setIJ (faceI,dim,Field->boundaryField() [PatchId][faceI][dim]) ;

  for (int dim=0; dim<3; dim++)  Vel->setInfoOnComponent (dim,directions[dim]);
  f->setArray (Vel);
  f->setName (FieldName);
  MEDLoader::WriteField (SolName,f,WriteFromScratch);
  f->decrRef();
  MedBoundMesh->decrRef();
  return;
}

void MedWithFOAM::PrintMedTimeSolution (ParaMEDMEM::MEDCouplingFieldDouble *Field,
                                        std::string FileName,
                                        double time,
                                        int iteration,
                                        int order
                                       ) {
  Field->setTime (time,iteration,order);
  MEDLoader::WriteFieldUsingAlreadyWrittenMesh (FileName,Field);
  return;
}

ParaMEDMEM::MEDCouplingFieldDouble *MedWithFOAM::GetGlobalSolution (Foam::List< Foam::List< Foam::scalar > > ScalarField,
    Foam::List< Foam::List< Foam::scalar > > GlobalMap,
    ParaMEDMEM::MEDCouplingUMesh *GlobalMesh,
    std::string FieldName) {
  int NTotCells = GlobalMesh->getNumberOfCells();
  ParaMEDMEM::DataArrayDouble *Tglobal = ParaMEDMEM::DataArrayDouble::New();
  Tglobal->alloc (NTotCells,1);
  for (int pp=0; pp<ScalarField.size(); pp++) {
    for (int jj=0; jj<ScalarField[pp].size(); jj++) {
      Tglobal->setIJ (GlobalMap[pp][jj], 0, ScalarField[pp][jj]);
    }
  }

  ParaMEDMEM::MEDCouplingFieldDouble *TFieldGlobal = ParaMEDMEM::MEDCouplingFieldDouble::New (ParaMEDMEM::ON_CELLS);
  TFieldGlobal->setMesh (GlobalMesh);
  TFieldGlobal->setArray (Tglobal);
  TFieldGlobal->setName (FieldName);

  Tglobal->decrRef();
  return TFieldGlobal;
}

ParaMEDMEM::MEDCouplingFieldDouble *MedWithFOAM::GetGlobalSolution (Foam::List< Foam::List< Foam::vector > > VectorField,
    Foam::List< Foam::List< Foam::scalar > > GlobalMap,
    ParaMEDMEM::MEDCouplingUMesh *GlobalMesh,
    std::string FieldName) {
  int NTotCells = GlobalMesh->getNumberOfCells();
  ParaMEDMEM::DataArrayDouble *Tglobal = ParaMEDMEM::DataArrayDouble::New();
  Tglobal->alloc (NTotCells,3);
  for (int pp=0; pp<VectorField.size(); pp++) {
    for (int jj=0; jj<VectorField[pp].size(); jj++) {
      Tglobal->setIJ (GlobalMap[pp][jj], 0, VectorField[pp][jj][0]);
      Tglobal->setIJ (GlobalMap[pp][jj], 1, VectorField[pp][jj][1]);
      Tglobal->setIJ (GlobalMap[pp][jj], 2, VectorField[pp][jj][2]);
    }
  }
  Tglobal->setInfoOnComponent (0,"x");
  Tglobal->setInfoOnComponent (1,"y");
  Tglobal->setInfoOnComponent (2,"z");

  ParaMEDMEM::MEDCouplingFieldDouble *TFieldGlobal = ParaMEDMEM::MEDCouplingFieldDouble::New (ParaMEDMEM::ON_CELLS);
  TFieldGlobal->setMesh (GlobalMesh);
  TFieldGlobal->setArray (Tglobal);
  TFieldGlobal->setName (FieldName);

  Tglobal->decrRef();
  return TFieldGlobal;
}

Foam::List< Foam::List< Foam::scalar > > MedWithFOAM::GatherCouplingMap (int proc, int procs) {

  int NCells = _OFLocal3D->getNumberOfCells();
  Foam::List<Foam::scalar> LocalToGlobalMap (NCells);
  for (unsigned int I=0; I<NCells; ++I) LocalToGlobalMap[I] =  _Global3Dto2D->getIJ (_ParToMedMap->getIJ (I,0),0);
  Foam::List< Foam::List<Foam::scalar> > gatheredMap (procs);
  gatheredMap[proc] = LocalToGlobalMap;
  Foam::Pstream::gatherList (gatheredMap);

  _IsMapGathered = true;
  return gatheredMap;
}

Foam::List< Foam::List< Foam::scalar > > MedWithFOAM::GatherScalarData (Foam::volScalarField *V,
    int proc, int procs) {

  int NCells = _OFLocal3D->getNumberOfCells();
  Foam::List<Foam::scalar> GlobalMap (NCells);
  for (unsigned int I=0; I<NCells; ++I) GlobalMap[I] =  V->internalField() [I];
  Foam::List< Foam::List<Foam::scalar> > gatheredMap (procs);
  gatheredMap[proc] = GlobalMap;
  Foam::Pstream::gatherList (gatheredMap);
  return gatheredMap;
}

Foam::List< Foam::List< Foam::vector > > MedWithFOAM::GatherVectorData (Foam::volVectorField *V,
    int proc, int procs) {

  int NCells = _OFLocal3D->getNumberOfCells();
  Foam::List<Foam::vector> GlobalMap (NCells);
  for (unsigned int I=0; I<NCells; ++I) GlobalMap[I] =  V->internalField() [I];
  Foam::List< Foam::List<Foam::vector> > gatheredMap (procs);
  gatheredMap[proc] = GlobalMap;
  Foam::Pstream::gatherList (gatheredMap);
  return gatheredMap;
}

ParaMEDMEM::MEDCouplingFieldDouble *MedWithFOAM::GetGlobalSolution (Foam::volScalarField *V,
    int proc, // Actual processor number
    int procs, // Number of all processors
    ParaMEDMEM::MEDCouplingUMesh *GlobalMesh, // 2D med mesh
    std::string FieldName) {
  if (!_IsMapGathered) _GatheredCouplingMap = GatherCouplingMap (proc, procs);
  Foam::List< Foam::List<Foam::scalar> > gatheredSolution = GatherScalarData (V, proc, procs);
  ParaMEDMEM::MEDCouplingFieldDouble *GatheredField = NULL;
  if (proc==0) GatheredField = GetGlobalSolution (gatheredSolution, _GatheredCouplingMap, GlobalMesh, FieldName);
  return GatheredField;
}
ParaMEDMEM::MEDCouplingFieldDouble *MedWithFOAM::GetGlobalSolution (Foam::volVectorField *V,
    int proc, // Actual processor number
    int procs, // Number of all processors
    ParaMEDMEM::MEDCouplingUMesh *GlobalMesh, // 2D med mesh
    std::string FieldName) {
  if (!_IsMapGathered) _GatheredCouplingMap = GatherCouplingMap (proc, procs);
  Foam::List< Foam::List<Foam::vector> > gatheredSolution = GatherVectorData (V, proc, procs);
  ParaMEDMEM::MEDCouplingFieldDouble *GatheredField = NULL;
  if (proc==0) GatheredField = GetGlobalSolution (gatheredSolution, _GatheredCouplingMap, GlobalMesh, FieldName);
  return GatheredField;
}



// ************************************************************************* //
