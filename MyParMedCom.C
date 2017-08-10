#include "MyParMedCom.h"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "MEDFileMesh.hxx"
#include "MEDCoupling.hxx"

MyParMedCom::MyParMedCom():__proc(0), __nprocs(1){}

MyParMedCom::MyParMedCom(int proc, 
			 int procs, 
			 int LocalToGlobalCellId[], 
			 int LocalToGlobalNodesId[], 
			 ParaMEDMEM::MEDCouplingUMesh * ProcMesh):
			 __proc(proc), __nprocs(procs)
			 {
  
  __ProcMesh = ProcMesh;
  int NumOfCells = __ProcMesh->getNumberOfCells();
  int NumOfNodes = __ProcMesh->getNumberOfNodes();
  _NumOfCells = NumOfCells;
  _ProcCellsLocToGlob = (int *) malloc(sizeof(int) * NumOfCells);
  _ProcNodesLocToGlob = (int *) malloc(sizeof(int) * NumOfNodes);
  for(int k=0; k<NumOfCells; k++)  _ProcCellsLocToGlob[k] = LocalToGlobalCellId[k];
  for(int k=0; k<NumOfNodes; k++)  _ProcNodesLocToGlob[k] = LocalToGlobalNodesId[k];
  
  if(__proc==0) {
    _ProcCells = (int *) malloc(sizeof(int) * __nprocs);
    _ProcNodes = (int *) malloc(sizeof(int) * __nprocs);
  }
  
  MPI_Comm_dup(MPI_COMM_WORLD, &_ClassComm);
  
  // Gathering the number of proc mesh cells and nodes
  MPI_Gather(&NumOfCells, 1, MPI_INT, _ProcCells, 1, MPI_INT, 0, _ClassComm); 
  MPI_Gather(&NumOfNodes, 1, MPI_INT, _ProcNodes, 1, MPI_INT, 0, _ClassComm); 
  

  if(__proc==0) {
     _Celldispls = (int *) malloc(sizeof(int) * __nprocs);
     _Nodesdispls = (int *) malloc(sizeof(int) * __nprocs);
     _Celldispls[0] = 0;
     _Nodesdispls[0] = 0;
     _globcells = _ProcCells[0];
     _globnodes = _ProcNodes[0];
     for(int i = 1; i<__nprocs; i++) {
       _globcells += _ProcCells[i];
       _globnodes += _ProcNodes[i];
       _Celldispls[i] = _Celldispls[i-1] + _ProcCells[i-1];
       _Nodesdispls[i] = _Nodesdispls[i-1] + _ProcNodes[i-1];
    }     
    _AllProcCellsLocToGlob = (int*) malloc(sizeof(int)*_globcells);
    _AllProcNodesLocToGlob = (int*) malloc(sizeof(int)*_globnodes);
  }

  // Gathering all proc local to global cell numbering -> map used to fill field on the global mesh
  MPI_Gatherv(LocalToGlobalCellId,  NumOfCells, MPI_INT, _AllProcCellsLocToGlob, _ProcCells, _Celldispls, MPI_INT, 0, _ClassComm);   
  MPI_Gatherv(LocalToGlobalNodesId,  NumOfNodes, MPI_INT, _AllProcNodesLocToGlob, _ProcNodes, _Nodesdispls, MPI_INT, 0, _ClassComm);   
}

MyParMedCom::~MyParMedCom(){
  free(_ProcCellsLocToGlob);
  if(__proc==0){
    free(_AllProcCellsLocToGlob);
    free(_ProcCells);
    free(_ProcNodes);
    free(_Celldispls);
    free(_Nodesdispls);
    free(__GlobFieldToScatter);
  }
  __ProcMesh = NULL;  
}

void MyParMedCom::GatherMedFieldOnProc0 (ParaMEDMEM::MEDCouplingFieldDouble *ProcField, ParaMEDMEM::MEDCouplingUMesh *GlobMesh)
{  
  ParaMEDMEM::TypeOfField SourceType = ProcField->getTypeOfField();   
  double *ArrayFromMedField = NULL;
  double *GatheredField = NULL;
  ParaMEDMEM::DataArrayDouble *GlobArray = ParaMEDMEM::DataArrayDouble::New(); 
  
  // SIZE OF MED FIELD AND GATHERING ARRAY ================================
  int NbOfTuples, NbOfGatheredValues;
  if(SourceType == ParaMEDMEM::ON_CELLS) {
    NbOfTuples = _globcells;
    NbOfGatheredValues = _globcells;
  }  
  if(SourceType == ParaMEDMEM::ON_NODES) {
    NbOfTuples = GlobMesh->getNumberOfNodes();
    NbOfGatheredValues = _globnodes;
  }// END SIZE ============================================================
  
  if(__proc==0){    
     GlobArray->alloc(NbOfTuples,ProcField->getNumberOfComponents());
     ArrayFromMedField = const_cast<double*>(GlobArray->getPointer()); 
     GatheredField = (double *) malloc(sizeof(double) * NbOfGatheredValues);
  } 
   
  // GATHERING PROCESS: FROM ALL PROCS TO PROC 0 ========================== 
  if(SourceType == ParaMEDMEM::ON_CELLS){  
    MPI_Gatherv(ProcField->getArray()->getPointer(),  
		ProcField->getMesh()->getNumberOfCells(), 
		MPI_DOUBLE, 
		GatheredField, 
		_ProcCells, 
		_Celldispls, 
		MPI_DOUBLE, 
		0, _ClassComm);  
  }
  if(SourceType == ParaMEDMEM::ON_NODES){  
    MPI_Gatherv(ProcField->getArray()->getPointer(),  
		ProcField->getMesh()->getNumberOfNodes(), 
		MPI_DOUBLE, 
		GatheredField, 
		_ProcNodes, 
		_Nodesdispls, 
		MPI_DOUBLE, 
		0, _ClassComm);  
  }  
  // END GATHERING PROCESS ================================================
  
  // REORDERING GATHERED FIELD INSIDE MED FIELD ===========================
  if(__proc==0){
   if(SourceType == ParaMEDMEM::ON_CELLS) 
      for(int j=0; j<GlobMesh->getNumberOfCells(); j++) 
	ArrayFromMedField[_AllProcCellsLocToGlob[j]] = GatheredField[j];
      
   if(SourceType == ParaMEDMEM::ON_NODES)
     for(int j=0; j<GlobMesh->getNumberOfNodes(); j++) 
       ArrayFromMedField[_AllProcNodesLocToGlob[j]] = GatheredField[j];
     
      if(__GlobFieldToGatherFilled) __GlobFieldToGather = NULL; // If __GlobFieldToGather is already filled we erase it   
      __GlobFieldToGather = ParaMEDMEM::MEDCouplingFieldDouble::New(SourceType);
      __GlobFieldToGather->setArray(GlobArray);
      __GlobFieldToGather->setMesh(GlobMesh);
      __GlobFieldToGather->setName(ProcField->getName());
      
      ArrayFromMedField = NULL;
      GlobArray->decrRef();
      free(GatheredField);      
      __GlobFieldToGatherFilled = true;
  }// END REORDERING =====================================================    
}

ParaMEDMEM::MEDCouplingFieldDouble *MyParMedCom::GetGatheredFieldOnProc0()
{
  if(__proc==0 && __GlobFieldToGatherFilled)  return __GlobFieldToGather;
  
  else if(__proc != 0) {
    printf("\033[1;31m WARNING: PROC %d DOESN'T HAVE THE GATHERED FIELD \n FUNCTION MUST BE CALLED BY PROC 0\n\033[0m",__proc);
    return NULL;
  } 
  else {
    printf("\033[1;31m WARNING: GLOB FIELD IS NOT FILLED: CALL FIRST GatherMedFieldOnProc0 FUNCTION \033[0m");
    return NULL;
  }
}

void MyParMedCom::LoadMedFieldToScatter (ParaMEDMEM::MEDCouplingFieldDouble *GlobField)
{
  if(__GlobFieldToScatterFilled) __GlobFieldToScatter = NULL;
  if(__proc==0) {
    double* TmpArray = const_cast<double*>(GlobField->getArray()->getPointer());  
    _FieldToScatterType = (GlobField->getTypeOfField() == ParaMEDMEM::ON_CELLS)? 0:1; 
    _CompToScatter = GlobField->getNumberOfComponents();
    // SIZE OF MED FIELD AND GATHERING ARRAY ===============================
    int NbOfTuples;
    if(!_FieldToScatterType) NbOfTuples = _globcells;  
    if(_FieldToScatterType) NbOfTuples = _globnodes;
    // END SIZE ============================================================  
    
    /* REORDERING FIELD VALUES
       __GlobFieldToScatter is then organized as
         proc0    proc1        procn
       [0...Nb0, 0...Nb1,     0...Nbn]
     */
    __GlobFieldToScatter = (double*) malloc(sizeof(double) * NbOfTuples);
    if(!_FieldToScatterType)
      for(int k=0; k<NbOfTuples; k++) 
	__GlobFieldToScatter[k] = TmpArray[_AllProcCellsLocToGlob[k]];
      
    if(_FieldToScatterType)
      for(int k=0; k<NbOfTuples; k++) 
	__GlobFieldToScatter[k] = TmpArray[_AllProcNodesLocToGlob[k]];  
    
    __GlobFieldToScatterFilled = true;
    TmpArray = NULL;
  }
  else printf("\[033[1;31m LoadMedFieldToScatter FUNCTION CANNOT BE CALLED BY PROC %d; \n IT MUST BE CALLED BY PROC 0 \n",__proc);
  return;
}

ParaMEDMEM::MEDCouplingFieldDouble *MyParMedCom::ScatterMedFieldFromProc0 (std::string FieldName)
{ 
  MPI_Bcast(&_FieldToScatterType, 1, MPI_INT, 0, _ClassComm);
  MPI_Bcast(&_CompToScatter, 1, MPI_INT, 0, _ClassComm);

  ParaMEDMEM::MEDCouplingFieldDouble *ScatteredField;
  if(!_FieldToScatterType) ScatteredField = ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS);
  else ScatteredField = ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES);

  // SIZE OF MED FIELD AND GATHERING ARRAY ================================
  int NbOfTuples;
  if(!_FieldToScatterType) NbOfTuples = __ProcMesh->getNumberOfCells();  
  if(_FieldToScatterType) NbOfTuples = __ProcMesh->getNumberOfNodes();
  // END SIZE ============================================================ 
  
  // DataArray where we store the scattered global field values
  ParaMEDMEM::DataArrayDouble *ProcArray =  ParaMEDMEM::DataArrayDouble::New();
  ProcArray->alloc(NbOfTuples, _CompToScatter); 
  double* ReceivingArray=const_cast<double*>(ProcArray->getPointer());    
  
  if(!_FieldToScatterType){
    MPI_Scatterv(__GlobFieldToScatter, 
		 _ProcCells, 
		 _Celldispls, 
		 MPI_DOUBLE, 
		 ReceivingArray, 
		 NbOfTuples, 
		 MPI_DOUBLE, 
		 0, 
		 _ClassComm);   
  }
  if(_FieldToScatterType){
    MPI_Scatterv(__GlobFieldToScatter, 
		 _ProcNodes, 
		 _Nodesdispls, 
		 MPI_DOUBLE, 
		 ReceivingArray, 
		 NbOfTuples, 
		 MPI_DOUBLE, 
		 0, 
		 _ClassComm);   
  }  

  ScatteredField->setMesh(__ProcMesh);   
  ScatteredField->setArray(ProcArray);  
  ScatteredField->setName(FieldName);
  ProcArray->decrRef();
  ReceivingArray = NULL;
  return ScatteredField;
}
