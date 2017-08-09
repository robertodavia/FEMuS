#ifndef __MedWithFOAM__
#define __MedWithFOAM__
#include <sstream>
#include <vector>
#include <map>
// #include "../../../../plat_base_packs/Salome-V7_8_0-public/tools/Medcoupling-V7_8_0/include/MEDCouplingAutoRefCountObjectPtr.hxx"

// MED
// #include "MEDLoader.hxx"
// #include "MEDCouplingUMesh.hxx"
// #include "MEDCouplingFieldDouble.hxx"

// OpenFOAM
#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "volFields.H"


// class fvMesh;
// class volVectorField;
// class volScalarField;

namespace ParaMEDMEM
{
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class DataArrayInt;
  class DataArrayDouble;
}

class MedWithFOAM
{
private:
  ParaMEDMEM::MEDCouplingUMesh * _FoamMesh = NULL;
  const ParaMEDMEM::MEDCouplingUMesh * _FemusMesh = NULL;
    Foam::fvMesh * _FoamFMesh = NULL;

  bool _areOFandFemusMeshSet = false;
public:
    ParaMEDMEM::DataArrayDouble * _Global2Dto3D = NULL;	// From global 2d med mesh to global 3d med mesh
    ParaMEDMEM::DataArrayDouble * _Global3Dto2D = NULL;	// From global 3d med mesh to global 2d med mesh
    ParaMEDMEM::DataArrayDouble * _ParToMedMap = NULL;	// From local OF mesh to global 3d med mesh (cell numbering)
    ParaMEDMEM::DataArrayDouble * _LOFtoLMED3D = NULL;	// From local OF mesh to local 3d med mesh
    ParaMEDMEM::DataArrayDouble * _LOFtoLMED2D = NULL;	// From local OF mesh to local 2d med mesh  
    ParaMEDMEM::DataArrayDouble * _LMED3DtoLOF = NULL;	// From local 3d med mesh to local OF mesh
    ParaMEDMEM::DataArrayDouble * _LMED2DtoLOF = NULL;	// From local 2d med mesh to local OF mesh

    // Map from local 2d med to global 2d med: (int) _Global3Dto2D[(int) _ParToMedMap[(int) _LMED2DtoLOF[k]]];
    
    
    ParaMEDMEM::MEDCouplingUMesh * _OFLocal3D = NULL;   // 3D MED mesh of local domain 
    ParaMEDMEM::MEDCouplingUMesh * _OFLocal2D = NULL;   // 2D MED mesh of local coupling interface mesh

    Foam::List < Foam::List < Foam::scalar >> _GatheredCouplingMap;
  bool _IsMapGathered = false;
    MedWithFOAM ();
    MedWithFOAM (ParaMEDMEM::MEDCouplingUMesh * FoamMesh,
		 const ParaMEDMEM::MEDCouplingUMesh * FemusMesh);
   ~MedWithFOAM ();
  void Clear ();

  // This function calculates all the needed map between local and global
  // OpenFOAM mesh and 3D and 2D MED meshes
  void SetParallelMap (Foam::fvMesh * FoamMesh,
		       ParaMEDMEM::MEDCouplingUMesh * Foam2D, const int proc);

  // This function returns the local 3D MED mesh of OpenFOAM parallel domain
  inline ParaMEDMEM::MEDCouplingUMesh * GetParMesh (const int proc)
  {
    return _OFLocal3D;
  }

  // This function returns the local 2D MED mesh of OpenFOAM parallel coupling interface
  inline ParaMEDMEM::MEDCouplingUMesh * Get2DParMesh (const int proc)
  {
    return _OFLocal2D;
  }

  void PrintMedMesh (Foam::fvMesh * mesh, int proc);
  ParaMEDMEM::MEDCouplingUMesh * GetMedMesh (Foam::fvMesh * mesh, int proc);
  void PrintMedSolution (Foam::volVectorField * Vel,	// The field
			 std::string FieldName,	// Field name
			 int proc,	// Processor number
			 ParaMEDMEM::MEDCouplingUMesh * MedMesh,	// Med mesh
			 bool WriteFromScratch = false	// If true it first creates the solution file
			 //                    int NComp,                             // Number of components
			 //                    std::vector<string>                    // Vector field component labels
    );
  void PrintMedSolution (Foam::volScalarField * Temp,	// The field
			 std::string FieldName,	// Field name
			 int proc,	// Processor number
			 ParaMEDMEM::MEDCouplingUMesh * MedMesh,	// Med mesh
			 bool WriteFromScratch = false	// If true it first creates the solution file
    );
  void PrintMedSolution (Foam::volScalarField * Field,
			 std::string FieldName,
			 int proc,
			 Foam::fvMesh * OFMesh,
			 std::string PatchName, bool WriteFromScratch);
  
  void PrintMedSolution (Foam::volVectorField * Field,
			 std::string FieldName,
			 int proc,
			 Foam::fvMesh * OFMesh,
			 std::string PatchName, bool WriteFromScratch);

  // This function returns a MED field for a given local Foam vector field
  ParaMEDMEM::MEDCouplingFieldDouble * GetMedSolution (Foam::volVectorField *Field,
						       std::string FieldName,
						       int proc,
						       ParaMEDMEM::
						       MEDCouplingUMesh *
						       MedMesh);

  // This function returns a MED field for a given local Foam scalar field
  ParaMEDMEM::MEDCouplingFieldDouble * GetMedSolution (Foam::volScalarField *Field,
						       std::string FieldName,
						       int proc,
						       ParaMEDMEM::
						       MEDCouplingUMesh *
						       MedMesh);

  // This function returns a MED field for a given local Foam scalar field on a given OpenFOAM boundary
  ParaMEDMEM::MEDCouplingFieldDouble * GetMedSolution (Foam::volScalarField *
						       Field,
						       std::string FieldName,
						       int proc,
						       Foam::fvMesh * OFMesh,
						       std::string PatchName);

  // This function returns a MED field of a global ScalarField - after solution recollection has been performed
  ParaMEDMEM::MEDCouplingFieldDouble * GetGlobalSolution (Foam::List <Foam::List <Foam::scalar >> ScalarField,
							  Foam::List <Foam::List <Foam::scalar >> GlobalMap,
							  ParaMEDMEM::MEDCouplingUMesh *GlobalMesh,
							  std::string FieldName);

  // This function returns a MED field of a global VectorField - after solution recollection has been performed
  ParaMEDMEM::MEDCouplingFieldDouble * GetGlobalSolution (Foam::List <Foam::List <Foam::vector >>VectorField,
							  Foam::List <Foam::List <Foam::scalar >> GlobalMap,
							  ParaMEDMEM::MEDCouplingUMesh *GobalMesh,
							  std::string FieldName);

  // This function returns a MED field of a global ScalarField - it recollects the global solution
  ParaMEDMEM::MEDCouplingFieldDouble * GetGlobalSolution (Foam::volScalarField * V,	// Field to print 
							  int proc,	// Actual processor id
							  int procs,	// Total number of processors
							  ParaMEDMEM::MEDCouplingUMesh * GlobalMesh,	//* 2D med mesh
							  std::string FieldName);

  // This function returns a MED field of a global VectorField - it recollects the global solution
  ParaMEDMEM::MEDCouplingFieldDouble * GetGlobalSolution (Foam::volVectorField * V, 
							  int proc, 
							  int procs,
		                                          ParaMEDMEM::MEDCouplingUMesh * GlobalMesh,
		                                          std::string FieldName);

  // This function gathers the global solution of a given vector field
  Foam::List<Foam::List<Foam::vector>> GatherVectorData (Foam::volVectorField * V, 
							  int proc,
							  int procs);
  
  // This function gathers the global solution of a given scalar field 
  Foam::List<Foam::List<Foam::scalar>> GatherScalarData (Foam::volScalarField * V, 
							 int proc,
							 int procs);
  
  // This function gathers the global coupling map
  Foam::List<Foam::List<Foam::scalar>> GatherCouplingMap (int proc,
							  int procs);

  void PrintMedTimeSolution (ParaMEDMEM::MEDCouplingFieldDouble * Field,
			     std::string FileName,
			     double time, 
			     int iteration, 
			     int order);
  
  void CreateMedBoundaryMesh (Foam::fvMesh * FoamMesh,
			      std::string patch, 
			      int proc);

  void SetMeshes (ParaMEDMEM::MEDCouplingUMesh * FoamMesh,
		  const ParaMEDMEM::MEDCouplingUMesh * FemusMesh);

  ParaMEDMEM::MEDCouplingUMesh * GetMedBoundaryMesh (Foam::fvMesh * FoamMesh,
						     std::string patch,
						     int proc);

};

#endif
