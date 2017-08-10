#ifndef __MyParMedCom__
#define __MyParMedCom__

// #ifdef HAVE_MPI
#include <mpi.h>   //For MPI_COMM_WORLD
// #endif

#ifdef HAVE_PETSCM
#include "petsc.h" // for Petsc solver
#endif

#include <sstream>
#include <vector>
#include <map>

namespace ParaMEDMEM
{
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class DataArrayInt;
  class DataArrayDouble;
}

class MyParMedCom
{
private:
  const int __proc;
  const int __nprocs;
  int _globcells=0;
  // total of nodes belonging to each proc mesh -> greater than number of nodes of global mesh 
  int _globnodes=0;
  // Specific proc mesh
  ParaMEDMEM::MEDCouplingUMesh * __ProcMesh = NULL;
  // Med Field we use to store the gathered field
  ParaMEDMEM::MEDCouplingFieldDouble *__GlobFieldToGather = NULL;
  
  int _CompToScatter;
  int _FieldToScatterType=0; // 0: On Cells, 1: On Nodes
  // Array we use to store the field to scatter - after values reordering
  double *__GlobFieldToScatter = NULL;
  // Flag we use to check whether __GlobFieldToGather is filled or not (every proc sees it)
  bool __GlobFieldToGatherFilled = false;
  // Flag we use to check whether __GlobFieldToScatter is filled or not (every proc sees it)
  bool __GlobFieldToScatterFilled = false;  
  // MPI communicator which is a duplicate of MPI_COMM_WORLD  
  MPI_Comm _ClassComm;
  
public:
  int _NumOfCells = 0;
  int * _ProcCellsLocToGlob = NULL;
  int * _ProcNodesLocToGlob = NULL;
  
  // structures filled by proc 0
  int * _ProcCells = NULL;
  int * _ProcNodes = NULL;
  int * _Celldispls = NULL;
  int * _Nodesdispls = NULL;  
  int * _AllProcCellsLocToGlob = NULL;
  int * _AllProcNodesLocToGlob = NULL;
  
  MyParMedCom();
  MyParMedCom(int proc, int procs, int LocalToGlobalCellId[], ParaMEDMEM::MEDCouplingUMesh * ProcMesh);
  MyParMedCom(int proc, int procs, int LocalToGlobalCellId[], int LocalToGlobalNodeId[], ParaMEDMEM::MEDCouplingUMesh * ProcMesh);
  ~MyParMedCom();
  
 /*! \brief The following functions allow to gather a med field on proc 0 and to get it.  
     First \c GatherMedFieldOnProc0 should be called by all procs, 
     then \c GetGatheredFieldOnProc0 should be called only by proc 0 in order to get 
     the gathered field
  */
  void GatherMedFieldOnProc0(ParaMEDMEM::MEDCouplingFieldDouble *ProcField, ParaMEDMEM::MEDCouplingUMesh *GlobMesh);
  ParaMEDMEM::MEDCouplingFieldDouble *GetGatheredFieldOnProc0();
  
 /*! \brief The following functions allow to scatter a med field from proc 0 on all procs.  
     First \c LoadMedFieldToScatter should be called by proc 0, giving the global field as input parameter, 
     then \c ScatterMedFieldFromProc0 should be called only by all procs in order to scatter the field
  */  
  void LoadMedFieldToScatter(ParaMEDMEM::MEDCouplingFieldDouble *GlobField);
  ParaMEDMEM::MEDCouplingFieldDouble *ScatterMedFieldFromProc0(std::string FieldName);
};

#endif