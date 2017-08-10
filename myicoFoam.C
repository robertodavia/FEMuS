/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

//  MED
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

// FEMUS
#include   "Printinfo_conf.h"

// solver library -------------------------------------
#include  "Solverlib_conf.h"  // Solver library options 
// Petsc
#ifdef HAVE_PETSCM
#include "petsc.h" // for Petsc solver
#endif
// Mpi
#ifdef HAVE_MPI
#include <mpi.h>   //For MPI_COMM_WORLD
#endif

// include local class
#include "MGFemusInit.h"
#include "MGUtils.h"
#include "MGSystem.h"
#include "MGGeomEl.h"
#include "MGMesh.h"
#include "MGFEMap.h"
#include "MGFE.h"
#include "MGEquationsSystem.h"
#include "MGTimeLoop.h"
#include "FEMUS.h"
#include "MEDLoader.hxx"
#include "MEDCouplingRemapper.hxx"
// #include "MMed.h"
// #include "DEC.hxx"
#include "InterfaceProjection.h"
// MED print class
#include "MedWithFOAM.h"

// PARALLEL MED
#include "ParaFIELD.hxx"
#include "ParaMESH.hxx"
#include "CommInterface.hxx"
#include "InterpKernelDEC.hxx"
#include "MPIProcessorGroup.hxx"

#include "MyParMedCom.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void PrintSystem(int color, std::string message);

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

  std::vector<MGUtils*> mgutils;
  std::string mesh_nameP[NUM_MESH];
  std::ostringstream filenameP[2];  std::ostringstream osfilenameP[2];

  for(int i_mesh=0; i_mesh< NUM_MESH; i_mesh++) {
    // MGUtils constructor ----------------------------------------------------
    mgutils.push_back(new MGUtils(i_mesh+1));
    // mesh names -------------------------------------------------------------
    mesh_nameP[i_mesh]= mgutils[i_mesh]->get_file("F_MESH_READ"); // name mesh
    int posP = mesh_nameP[i_mesh].find(".");  // position of "live" in str
    filenameP[i_mesh] <<   mesh_nameP[i_mesh].substr(0,posP)  << "_MedToMg.med" ;
    osfilenameP[i_mesh]<< mgutils[i_mesh]->_mesh_dir <<filenameP[i_mesh].str();
    std::cout<<" \n P mesh file "<< i_mesh+1 << "= "<< osfilenameP[i_mesh].str().c_str() <<"\n "; 
  }
  
  MGFEMap *mgfemap; mgfemap=new MGFEMap();
  MGFE *dfe_q;    dfe_q=new MGFE(2,ELTYPE); dfe_q->init_qua();
  mgfemap->set_FE(dfe_q);                                        // initialize quadratic fem
  MGFE *dfe_l;  dfe_l=new MGFE(1,ELTYPE); dfe_l->init_lin();
  mgfemap->set_FE(dfe_l);                                        // initialize linear fem
  MGFE *dfe_k; dfe_k=new MGFE(0,ELTYPE);  dfe_k->init_pie();
  mgfemap->set_FE(dfe_k);      
  
  MGGeomEl *mggeomel;  mggeomel=new  MGGeomEl();
  std::vector<FIELDS> myproblemP; myproblemP.resize(2);
  // Problem to solve for each mesh
   myproblemP[0]=NS_F; 
   myproblemP[1]=T_F;
  // MGFemusInit --------------------------------------------------------------
  FEMUS P;                                        // constructor
  P.init_param(*mgutils[0],0);          // init parameter
  P.init_fem(*mggeomel,*mgfemap);                 // init fem      
  P.setMesh();                                    // set MGmesh   
  P.setSystem(myproblemP);                        // set system
  P.setMeshName(mesh_nameP[0]);
//   
  P.init_interface(11,1,2,filenameP[0].str().c_str());
  MedWithFOAM extractor;
  
//   InterpKernelDEC dec(groupA, groupB);

  
  std::map<std::string, bool> YesNo;
  YesNo["yes"] = true;
  YesNo["no"]  = false;  
  
  mgutils[0]->read("/DATA/par2.in");
  std::string Foam2DMesh      = mgutils[0]->get_file("FoamMesh2D");
  std::string Foam3DMesh      = mgutils[0]->get_file("FoamMesh3D");
  std::string SolveFem        = mgutils[0]->get_file("SolveFemus");
  std::string UseFem          = mgutils[0]->get_file("UseFemus");
  std::string dist            = mgutils[0]->get_file("DistributedPrinting");
  std::string P2interpolation = mgutils[0]->get_file("P2Interp");

  const bool SolveFemus = YesNo[SolveFem];
  const bool UseFemus = YesNo[UseFem];
  const bool P2Interp = YesNo[P2interpolation];
  const bool Distributed = YesNo[dist];
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  const int proc = Pstream::myProcNo(); 
  std::ostringstream stm ;
  stm << proc;
  std::string Proc = stm.str();
  std::string TimeName = "RESU_MED/Time_p"+Proc+".med";
  
  
  const ParaMEDMEM::MEDCouplingUMesh* FemusMesh = P.getUMesh(11);                                // MESH OF FEMUS INTERFACE
  ParaMEDMEM::MEDCouplingUMesh * Foam2D = MEDLoader::ReadUMeshFromFile("MESH/"+Foam2DMesh, 0);   // 2D MED MESH FOR OPENFOAM COUPLING - 3D MESH CONTRACTION ON Z=0 PLANE
  ParaMEDMEM::MEDCouplingUMesh * Foam3D = MEDLoader::ReadUMeshFromFile("MESH/"+Foam3DMesh, 0);   // 3D MED MESH OF OPENFOAM DOMAIN
  
  ParaMEDMEM::MEDCouplingUMesh* FemusMesh2 = MEDLoader::ReadUMeshFromFile("MESH/femus_20x20_MedToMg.med", 0);         // MESH OF FEMUS INTERFACE
  ParaMEDMEM::MEDCouplingUMesh * FemusPart = MEDLoader::ReadUMeshFromFile("MESH/femus_20x20_info.med", 0);
  ParaMEDMEM::MEDCouplingFieldDouble *ProcField = MEDLoader::ReadField(ParaMEDMEM::ON_CELLS, "MESH/femus_20x20_info.med", "Mesh_1", 0, "Proc", -1,-1);
  
  std::vector<int> CellIdProc;
  for(int i=0; i<ProcField->getNumberOfTuples(); i++){
    if((int) ProcField->getIJ(i,0) == proc) CellIdProc.push_back(i);
  }
  int *CellId = new int[CellIdProc.size()];
  for(int i=0; i<CellIdProc.size(); i++) CellId[i] = CellIdProc[i];
  ParaMEDMEM::MEDCouplingUMesh *FemusPar = FemusPart->buildPartOfMySelf (CellId, CellId + CellIdProc.size());
  FemusPar->zipCoords(); // It changes nodes numbering -> new numbering is from 0 to FemusPar->getNumberOfNodes
  FemusPar->setName("FemusProcMesh");

  // CREATING A MAP FROM LOCAL NODES NUMBERING TO GLOBAL NODES NUMBERING - RELATIVE TO MED FILE 
  int *NodesId = new int[FemusPar->getNumberOfNodes()];
  
  std::vector<int> FemusCellConn;
  std::vector<int> FemusCellConn2;
  int FemusCells = FemusPar->getNumberOfCells();
  int FemusNodesPerCell = FemusPar->getNumberOfNodesInCell(0);
  for (int i_cell=0; i_cell<FemusCells; i_cell++) {
    FemusPar->getNodeIdsOfCell (i_cell,FemusCellConn);
    FemusMesh->getNodeIdsOfCell (CellId[i_cell],FemusCellConn2);
    for (int i_cnode=0; i_cnode<FemusNodesPerCell; i_cnode++)  {
      NodesId[FemusCellConn[i_cnode]] = FemusCellConn2[i_cnode];
    }
    FemusCellConn.clear();
    FemusCellConn2.clear();
  }// END MAP ---------------------------------------------------------------------------------

  
  MedWithFOAM Printer(Foam3D, Foam2D);  
  Printer.SetParallelMap(&mesh, Foam2D, proc);
  
  ParaMEDMEM::MEDCouplingUMesh *Par3DFoamMesh = Printer.GetParMesh(proc);   // LOCAL 3D MED MESH OF OPENFOAM DOMAIN
  ParaMEDMEM::MEDCouplingUMesh *Par2DFoamMesh = Printer.Get2DParMesh(proc); // LOCAL 2D MED MESH OF OPENFOAM COUPLING MESH
  
  const int NCells = mesh.nCells();
  
//   Info<< mesh.cellPoints()[0] << endl;   
  MMed * FEMuSMed = new BoundInterp();
  FEMuSMed->InitFe(2);
  FEMuSMed->FillParameters(FemusMesh, Par2DFoamMesh, Volume);  // INTERPOLATION WILL BE EXECUTED FROM FEMUS WHOLE MESH TO OPENFOAM LOCAL MESH
  MMed OFMed;
  OFMed.InitFe(3);
  
  
  int    n_steps = mgutils[0]->get_par("nsteps");
  double      dt = mgutils[0]->get_par("dt");
  int print_step = mgutils[0]->get_par("printstep");
  int    itime_0 = mgutils[0]->get_par("itime");
  double time    = 0.;
  int itime      = 0;
  int iter       = 0;

  
  if(Distributed) MEDLoader::WriteUMesh(TimeName,Par3DFoamMesh,true);
  else if(proc==0) MEDLoader::WriteUMesh(TimeName,Foam2D,true);
  
  P.solve_setup(itime_0,time);                    // initial time loop (t=0)
  
  Info<< "\nStarting time loop\n" << endl;    
    
  ParaMEDMEM::MEDCouplingFieldDouble *FieldFromFemus=NULL, *FemusAvCellField=NULL, *TFieldOnFoamMesh=NULL, *P2OnFoamMesh=NULL;
  ParaMEDMEM::MEDCouplingRemapper remapper;
   
  remapper.setPrecision(1e-11);
  remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
  remapper.prepare(FemusMesh,Par2DFoamMesh,"P0P0");  
  
  
  FieldFromFemus =P.getValuesOnBoundary(11,"T",1,0);
  P2OnFoamMesh = FEMuSMed->InterpolatedField(FieldFromFemus);   
  FemusAvCellField = FEMuSMed->GetCellField(P2OnFoamMesh); 
  FemusAvCellField->setName("CellT");
  FemusAvCellField->setNature(ParaMEDMEM::ConservativeVolumic);
  
  
  if(UseFemus){              
    if(!SolveFemus){            
      if(P2Interp){	
	for(int i=0; i<NCells; i++) { Buoyant[i][1] = 9.81*10.*(FemusAvCellField->getIJ(Printer._LOFtoLMED2D->getIJ(i,0),0) -573.); }
	FEMuSMed->PrintMed(P2OnFoamMesh, "Estratto"+std::to_string(proc), false);
      }
      else{	
        TFieldOnFoamMesh=remapper.transferField(FemusAvCellField,1000.);
        TFieldOnFoamMesh->setName("CellT");
        TFieldOnFoamMesh->setNature(ParaMEDMEM::ConservativeVolumic);      
        for(int i=0; i<NCells; i++)  { Buoyant[i][1] = 9.81*10.*(TFieldOnFoamMesh->getIJ(Printer._LOFtoLMED2D->getIJ(i,0),0) -573.); }
      }     
    }
  }
  
  
  while (runTime.loop()){
    Info<< "Time = " << runTime.timeName() << nl << endl;

    Info<< "\033[038;5;202m \n======================================================================\n"<<
       "\t\tOpenFOAM Iteration"<<
       "\n======================================================================\n\033[0m" << endl;
  
#include "readPISOControls.H"
#include "CourantNo.H"

    if(SolveFemus && UseFemus){// READ FEMUS FIELD AT EVERY TIME STEP   
      FieldFromFemus = P.getValuesOnBoundary(11,"T",1,0);
      
      if(P2Interp){// P2P2 INTERPOLATION ON FOAM SUPPORT AND THEN CALCULATION OF CELL MEAN INTEGRAL VALUES
	P2OnFoamMesh = FEMuSMed->InterpolatedField(FieldFromFemus);      
        FemusAvCellField = FEMuSMed->GetCellField(P2OnFoamMesh); 
	FemusAvCellField->setName("FemusT");
        FemusAvCellField->setNature(ParaMEDMEM::ConservativeVolumic);
	for(int i=0; i<NCells; i++) { Buoyant[i][1] = 9.81*10.*(FemusAvCellField->getIJ(Printer._LOFtoLMED2D->getIJ(i,0),0) -573.); }
      }
      else{// CELL MEAN INTEGRAL VALUES CALCULATION AND THEN P0P0 INTERPOLATION ON FOAM SUPPORT
	FemusAvCellField = FEMuSMed->GetCellField(FieldFromFemus); 
	FemusAvCellField->setName("FemusT");
        FemusAvCellField->setNature(ParaMEDMEM::ConservativeVolumic);
	
        TFieldOnFoamMesh=remapper.transferField(FemusAvCellField,1000.);
        TFieldOnFoamMesh->setName("cellT+");
        TFieldOnFoamMesh->setNature(ParaMEDMEM::ConservativeVolumic);      
        for(int i=0; i<NCells; i++) { Buoyant[i][1] = 9.81*10.*(TFieldOnFoamMesh->getIJ(Printer._LOFtoLMED2D->getIJ(i,0),0) -573.); }
      }   
    } 
    
    if(!UseFemus) {for(int i=0; i<NCells; i++) Buoyant[i][1] = 9.81*10.*(T[i] -573.);}
           
    fvVectorMatrix UEqn(
        fvm::ddt(U)
      + fvm::div(phi, U)
      - fvm::laplacian(nu, U)
      ==
        Buoyant
    );
    PrintSystem(46, "Velocity"); 
    solve(UEqn == -fvc::grad(p) );
    
    // --- PISO loop

    for (int corr = 0; corr < nCorr; corr++){
      volScalarField rUA = 1.0/UEqn.A();

      U = rUA*UEqn.H();
      phi = (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, U, phi);

      adjustPhi(phi, U, p);

      for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++){
        fvScalarMatrix pEqn(
            fvm::laplacian(rUA, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
	PrintSystem(117, "Pressure");
        pEqn.solve();

        if (nonOrth == nNonOrthCorr){
          phi -= pEqn.flux();
        }
      }

#include "continuityErrs.H"

      U -= rUA*fvc::grad(p);
      U.correctBoundaryConditions();
    }
        
    fvScalarMatrix TEqn(
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(DT, T)
    );

    PrintSystem(209, "Temperature");
    TEqn.solve();
	
    iter ++;
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
	    
	    
    Info<< "\033[038;5;39m \n======================================================================\n"<<
           "\t\tFEMuS Iteration"<<
	   "\n======================================================================\n\033[0m" << endl;	
	
    if(SolveFemus) P.solve_onestep(itime_0,itime,print_step,time,dt);    // solving P   
    itime ++;    
	   
    int wt = (int) runTime.writeInt();
    runTime.write();
    
    if(itime % wt == 0){//MED TIME WRITING
      double itTime = stod(runTime.timeName());
      if(Distributed){// EVERY PROC PRINTS ITS OWN PART OF SOLUTION ON ITS LOCAL DOMAIN
        ParaMEDMEM::MEDCouplingFieldDouble * Field = Printer.GetMedSolution(&U,"Velocity",proc,Par3DFoamMesh); 	 
        Printer.PrintMedTimeSolution(Field,TimeName, itTime, iter, 0); 
        Field = Printer.GetMedSolution(&T,"Temperature",proc,Par3DFoamMesh); 
        Printer.PrintMedTimeSolution(Field,TimeName, itTime, iter, 0);  
        Field = Printer.GetMedSolution(&Buoyant,"Buoyant",proc,Par3DFoamMesh); 
        Printer.PrintMedTimeSolution(Field,TimeName, itTime, iter, 0);  
	Field->decrRef();
      }else{// PROC 0 PRINTS THE WHOLE SOLUTION ON THE GLOBAL DOMAIN AFTER SOLUTION RECOLLECTION
	ParaMEDMEM::MEDCouplingFieldDouble *Field1 = Printer.GetGlobalSolution(&T, proc, Pstream::nProcs(), Foam2D, "T");
	ParaMEDMEM::MEDCouplingFieldDouble *Field2 = Printer.GetGlobalSolution(&U, proc, Pstream::nProcs(), Foam2D, "V");
	ParaMEDMEM::MEDCouplingFieldDouble *Field3 = Printer.GetGlobalSolution(&Buoyant, proc, Pstream::nProcs(), Foam2D, "Buoyant");
	if(proc==0) Printer.PrintMedTimeSolution(Field1,TimeName, itTime, iter, 0); 
	if(proc==0) Printer.PrintMedTimeSolution(Field2,TimeName, itTime, iter, 0); 
	if(proc==0) Printer.PrintMedTimeSolution(Field3,TimeName, itTime, iter, 0); 
	Field1=NULL;
	Field2=NULL;
	Field3=NULL;
      }
    }   
  }
  
  MyParMedCom *EX = new MyParMedCom(proc, Pstream::nProcs(), CellId, NodesId, FemusPar);
  if(proc==0){
    ParaMEDMEM::MEDCouplingFieldDouble *TF = P.getValuesOnBoundary(11,"T",1,0);
    MEDLoader::WriteField("Tbefore.med",TF,true);
    EX->LoadMedFieldToScatter(TF);
  }
  ParaMEDMEM::MEDCouplingFieldDouble *ScatterProc = EX->ScatterMedFieldFromProc0("TS");
  printf("%d \n",ScatterProc->getNumberOfTuples());
  
  MEDLoader::WriteField("T_"+std::to_string(proc)+".med",ScatterProc,true);
  P.setExtField("NS0X",ScatterProc);
  P.setExtField("NS0Y",ScatterProc);
  P.solve_onestep(itime_0,itime,print_step,time,dt);    // solving P  
  if(proc==0){
    ParaMEDMEM::MEDCouplingFieldDouble *TF = P.getValuesOnBoundary(11,"NS0X",1,0);
    MEDLoader::WriteField("Tfinal.med",TF,true);
  }

//   const ParaMEDMEM::CommInterface comm_interface;
//   std::set<int> proc_receiving;
//   std::set<int> proc_sending;
// 
//   proc_receiving.insert(proc);
//   for(int i=0; i<Pstream::nProcs(); i++) proc_sending.insert(i);
//   
//   
//   int PCells = FemusPar->getNumberOfCells();
// 
//   double *LMED2DtoLOF = Printer._LMED2DtoLOF->getPointer();
//   double *ParToMedMap = Printer._ParToMedMap->getPointer();
//   double *Global3Dto2D = Printer._Global3Dto2D->getPointer();   
//    
//   int *foamint = (int *)malloc(sizeof(int) * Par2DFoamMesh->getNumberOfCells());
//   for(int k=0; k<Par2DFoamMesh->getNumberOfCells(); k++) foamint[k] = (int) Global3Dto2D[(int) ParToMedMap[(int) LMED2DtoLOF[k]]];
//    
//   MyParMedCom *EX = new MyParMedCom(proc, Pstream::nProcs(), foamint, Par2DFoamMesh);
//   ParaMEDMEM::MEDCouplingFieldDouble *Field = Printer.GetMedSolution(&T,"Temperature",proc,Par2DFoamMesh); 
//   MEDLoader::WriteField("Par_prel"+std::to_string(proc)+".med",Field, true);
//   EX->GatherMedFieldOnProc0(Field, Foam2D);
//   ParaMEDMEM::MEDCouplingFieldDouble *GatheredField;
//    
//   if(proc==0) {
//     GatheredField = EX->GetGatheredFieldOnProc0();    
//     MEDLoader::WriteField("Gathered.med",GatheredField, true);
//     EX->LoadMedFieldToScatter(GatheredField);
//   }
// 
//   ParaMEDMEM::MEDCouplingFieldDouble *ff = EX->ScatterMedFieldFromProc0(Field->getName());
//   MEDLoader::WriteField("Par"+std::to_string(proc)+".med",ff, true);
//   
//   
//   remapper.setPrecision(1e-11);
//   remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
//   remapper.prepare(Foam2D,FemusMesh,"P0P0");  
//   ParaMEDMEM::MEDCouplingFieldDouble *Field1 = Printer.GetGlobalSolution(&T, proc, Pstream::nProcs(), Foam2D, "T");
//   ParaMEDMEM::MEDCouplingFieldDouble *FieldForFemus;
  
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
  
//   
//   if(P2Interp) {   
//     const int FoamCells = Foam2D->getNumberOfCells();
//     const int FemusCells = FemusMesh->getNumberOfCells();
//     
//     ParaMEDMEM::MEDCouplingFieldDouble *TFoamField2D=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS);
//     TFoamField2D->setMesh(Foam2D);
//     
//     std::string Name = "OF" + to_string(FoamCells) + "_FE" + to_string(FemusCells) + "_" + UseFem + "prova";
//     
//     ParaMEDMEM::DataArrayDouble * CellArray = ParaMEDMEM::DataArrayDouble::New();
//     CellArray->alloc(FoamCells,1);
//     for(int i=0; i<FoamCells; i++) { CellArray->setIJ(i,0,TField->getIJ(Printer._CouplingMap->getIJ(i,0),0)); }
//         
//     TFoamField2D->setArray(CellArray);
//     TFoamField2D->setName("FoamT");
//     TFoamField2D->setNature(ParaMEDMEM::ConservativeVolumic);
//     
//     T_Norm_Fe_i = FEMuSMed->Integrate(P2OnFoamMesh,2,1,0, NormL2);
//     
// /*    ParaMEDMEM::MEDCouplingFieldDouble *TDifference = (*FemusAvCellField) - (*TFoamField2D);
//     TDifference->setName("TDiff");
//     norm_diff = FEMuSMed->Integrate(TDifference,0,1,0, NormL2);   */ 
//     
//     FEMuSMed->PrintMed(TFoamField2D,Name,false);
//     FEMuSMed->PrintMed(FemusAvCellField,Name,true);    
// //     FEMuSMed->PrintMed(TDifference,Name,true);   
//     
// //     TDifference->decrRef();  
//     
//     
//     ParaMEDMEM::MEDCouplingFieldDouble * VelField = FEMuSMed->GetVelocityField(&P, 11);
//     ParaMEDMEM::MEDCouplingFieldDouble * VelFieldInt = FEMuSMed->InterpolatedField(VelField);     
//     ParaMEDMEM::MEDCouplingFieldDouble * VelFieldCell = FEMuSMed->GetCellField(VelFieldInt); 
//     VelFieldCell->setName("FE_VEL");
//     VelFieldCell->getArray()->setInfoOnComponent(0,"x");
//     VelFieldCell->getArray()->setInfoOnComponent(1,"y");
//     VelFieldCell->setNature(ParaMEDMEM::ConservativeVolumic);
//     
//     ParaMEDMEM::MEDCouplingFieldDouble * VField = Printer.GetMedSolution(&U,"OF_VEL",proc,Foam3D);
//     ParaMEDMEM::DataArrayDouble * VArray = ParaMEDMEM::DataArrayDouble::New();
//     VArray->alloc(FoamCells,2);
//     for(int j=0; j<2; j++)
//       for(int i=0; i<FoamCells; i++){
// 	VArray->setIJ(i,j,VField->getIJ(Printer._CouplingMap->getIJ(i,0),j)); 
//       }
//     VArray->setInfoOnComponent(0,"x");
//     VArray->setInfoOnComponent(1,"y");
//     ParaMEDMEM::MEDCouplingFieldDouble *VFoamField2D=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS);
//     VFoamField2D->setMesh(Foam2D);
//     VFoamField2D->setArray(VArray);
//     VFoamField2D->setName("OF_VEL");
//     VFoamField2D->setNature(ParaMEDMEM::ConservativeVolumic);
//     
//     FEMuSMed->PrintMed(VelFieldCell,Name,true); 
//     FEMuSMed->PrintMed(VFoamField2D,Name,true); 
//     
//     
//     ParaMEDMEM::DataArrayInt * NodConn = Foam2D->getNodalConnectivity();
//     ParaMEDMEM::DataArrayInt * NodConnInd = Foam2D->getNodalConnectivityIndex();
//     ParaMEDMEM::DataArrayDouble * CONN_TABLE= ParaMEDMEM::DataArrayDouble::New();
//     CONN_TABLE->alloc(Foam2D->getNumberOfCells(), Foam2D->getNumberOfNodesInCell(0));
//     for(int j=0; j<Foam2D->getNumberOfCells(); j++){
//       const int in = NodConnInd->getIJ(j,0);
//       const int fin = NodConnInd->getIJ(j+1,0);
//       for(int l=in+1; l<fin; l++){
// 	CONN_TABLE->setIJ(j,l-in-1, NodConn->getIJ(l,0));   
//       }
//     }
//     for(int i=0; i<Foam2D->getNumberOfNodesInCell(0); i++) CONN_TABLE->setInfoOnComponent(i,to_string(i));
//     ParaMEDMEM::MEDCouplingFieldDouble *CONN=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS);
//     CONN->setMesh(Foam2D);
//     CONN->setArray(CONN_TABLE);
//     CONN->setName("Connectivity");
//     FEMuSMed->PrintMed(CONN,Name,true); 
//     
// //     ParaMEDMEM::MEDCouplingField * NODE ;//= ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES);
// //     NODE->setArray(NodConn);
//     
// //     FEMuSMed->PrintMed(NODE,Name,true); 
//     
//     TFoamField2D->decrRef();  CellArray->decrRef();
//     VelField->decrRef();  VelFieldInt->decrRef();
//     VelFieldCell->decrRef();
//   }
//   
//   std::cout<<"Norma l2 della differenza fra i campi: "<<norm_diff<<std::endl;
//   if(UseFemus){
//     std::cout<<"\033[038;5;117;1m Temperature field - L2 norm from FEMuS: "<<T_Norm_Fe
//              <<"\n Norm of interpolated field on FOAM support "<<T_Norm_Fe_i
//              <<"\n and from OpenFOAM, using FEMuS source term "<<T_Norm_OF/0.01              
//              <<"\033[0m \n";
//   }
//   if(!UseFemus){
//     std::cout<<"\033[038;5;117;1m Temperature field - L2 norm from FEMuS: "<<T_Norm_Fe
//              <<" and from OpenFOAM, using OpenFOAm source term "<<T_Norm_OF/0.01 
//              <<"\033[0m \n";
//   }
  
  Info<<"\n";
  Info<< "End MyIcoFoam \n" << endl;

 
  
  P.terminate(); 
  mgutils.clear();  
  
// Cleaning med structures
//   remapper.~MEDCouplingRemapper();  
//   TFieldOnFoamMesh->decrRef();
  
  Foam3D->decrRef();
  Foam2D->decrRef();
  FemusMesh->decrRef();
  P2OnFoamMesh->decrRef();
//   TField->decrRef();
  
  FieldFromFemus->decrRef();
  FemusAvCellField->decrRef();
  
  if(TFieldOnFoamMesh != NULL)TFieldOnFoamMesh->decrRef();
  
  delete dfe_q;    
  delete dfe_l; 
  delete dfe_k; 
  
  delete mggeomel; 
  delete mgfemap; 
//     
  Printer.Clear();  
  return 0;
}


void PrintSystem(int color, std::string message)
{
  std::cout  << "\033[038;5;"<<color<<";1m \n"
	<< "---------------------------------------------------\n\t"
        << message
        << "\n---------------------------------------------------\n  \033[0m";	
  
return;
}

