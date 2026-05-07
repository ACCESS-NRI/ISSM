/*!\file: issm_nuopc_capi.cpp
 * \brief: Minimal handle-based ISSM coupling API for an in-process NUOPC cap.
 */

#include "./issm_nuopc_capi.h"

#include <string>

#include <ESMC.h>

#include "../modules/GetVectorFromInputsx/GetVectorFromInputsx.h"
#include "../modules/InputUpdateFromVectorx/InputUpdateFromVectorx.h"
#include "../shared/Enum/Enum.h"

struct ISSM_NUOPC_Model {
	FemModel* femmodel;
};

static std::string EnsureTrailingSlash(const char* path){/*{{{*/
	std::string root(path ? path : ".");
	if(root.empty()) root = ".";
	if(root[root.size()-1] != '/') root += "/";
	return root;
}/*}}}*/

static ISSM_MPI_Comm CommFromFortran(int mpi_comm_f){/*{{{*/
#if defined(_HAVE_MPI_)
	return MPI_Comm_f2c(mpi_comm_f);
#else
	return mpi_comm_f;
#endif
}/*}}}*/

static void InitializePetscForCoupling(ISSM_MPI_Comm comm){/*{{{*/
#ifdef _HAVE_PETSC_
	PetscBool initialized = PETSC_FALSE;
	PetscErrorCode ierr;
	PetscInitialized(&initialized);
	if(!initialized){
		PETSC_COMM_WORLD = comm;
		ierr = PetscInitializeNoArguments();
		if(ierr) _error_("Could not initialize Petsc for ISSM_NUOPC");
	}
#else
	(void)comm;
#endif
}/*}}}*/

static void BuildCaseFilePaths(const char* case_dir, const char* model_name, std::string& rootpath, std::string& binfilename, std::string& outbinfilename, std::string& toolkitsfilename, std::string& lockfilename, std::string& restartfilename){/*{{{*/
	rootpath         = EnsureTrailingSlash(case_dir);
	binfilename      = rootpath + model_name + ".bin";
	outbinfilename   = rootpath + model_name + ".outbin";
	toolkitsfilename = rootpath + model_name + ".toolkits";
	lockfilename     = rootpath + model_name + ".lock";
	restartfilename  = rootpath + model_name + "_rank" + std::to_string(IssmComm::GetRank()) + ".rst";
}/*}}}*/

static ISSM_NUOPC_Model* GetModel(void* model_handle){/*{{{*/
	if(!model_handle) _error_("ISSM_NUOPC received a null model handle");
	ISSM_NUOPC_Model* model = reinterpret_cast<ISSM_NUOPC_Model*>(model_handle);
	if(!model->femmodel) _error_("ISSM_NUOPC model handle is missing a FemModel instance");
	return model;
}/*}}}*/

static void ExportVertexField(FemModel* femmodel, int input_enum, double* values, int size){/*{{{*/
	IssmDouble* field_values = NULL;
	GetVectorFromInputsx(&field_values, femmodel, input_enum, VertexSIdEnum);
	for(int i=0;i<size;i++) values[i] = static_cast<double>(field_values[i]);
	xDelete<IssmDouble>(field_values);
}/*}}}*/

void* ISSM_NUOPC_CreateFromCase(const char* case_dir, const char* model_name, const char* solution_name, int mpi_comm_f){/*{{{*/
	ISSM_MPI_Comm comm = CommFromFortran(mpi_comm_f);
	InitializePetscForCoupling(comm);
	IssmComm::SetComm(comm);

	std::string rootpath;
	std::string binfilename;
	std::string outbinfilename;
	std::string toolkitsfilename;
	std::string lockfilename;
	std::string restartfilename;
	BuildCaseFilePaths(case_dir, model_name, rootpath, binfilename, outbinfilename, toolkitsfilename, lockfilename, restartfilename);

	int solution_type = StringToEnumx(solution_name);
	if(solution_type != TransientSolutionEnum){
		_error_("ISSM_NUOPC first version only supports TransientSolution");
	}

	ISSM_NUOPC_Model* model = new ISSM_NUOPC_Model();
	model->femmodel = new FemModel(const_cast<char*>(rootpath.c_str()), const_cast<char*>(binfilename.c_str()), const_cast<char*>(outbinfilename.c_str()), const_cast<char*>(toolkitsfilename.c_str()), const_cast<char*>(lockfilename.c_str()), const_cast<char*>(restartfilename.c_str()), const_cast<char*>(model_name), comm, solution_type, NULL);

	int legacy_ocean_coupling = 0;
	model->femmodel->parameters->FindParam(&legacy_ocean_coupling, TransientIsoceancouplingEnum);
	if(legacy_ocean_coupling != 0){
		_error_("ISSM_NUOPC first version expects md.transient.isoceancoupling = 0 so that the cap owns the coupling exchange");
	}

	return reinterpret_cast<void*>(model);
}/*}}}*/

void ISSM_NUOPC_Destroy(void* model_handle){/*{{{*/
	ISSM_NUOPC_Model* model = GetModel(model_handle);
	OutputResultsx(model->femmodel);
	model->femmodel->CleanUp();
	delete model->femmodel;
	model->femmodel = NULL;
	delete model;
}/*}}}*/

void ISSM_NUOPC_WriteRestart(void* model_handle){/*{{{*/
	ISSM_NUOPC_Model* model = GetModel(model_handle);
	model->femmodel->CheckPoint();
}/*}}}*/

void ISSM_NUOPC_GetMeshCounts(void* model_handle, int* num_nodes, int* num_elements){/*{{{*/
	ISSM_NUOPC_Model* model = GetModel(model_handle);
	*num_nodes    = model->femmodel->vertices->Size();
	*num_elements = model->femmodel->elements->Size();
}/*}}}*/

void ISSM_NUOPC_GetMeshNodes(void* model_handle, int* node_ids, int* node_owners, double* node_coords){/*{{{*/
	ISSM_NUOPC_Model* model = GetModel(model_handle);
	int rank = IssmComm::GetRank();

	for(int i=0;i<model->femmodel->vertices->Size();i++){
		Vertex* vertex = xDynamicCast<Vertex*>(model->femmodel->vertices->GetObjectByOffset(i));
		node_ids[i]           = vertex->Sid() + 1;
		node_owners[i]        = rank;
		node_coords[2*i]      = static_cast<double>(vertex->x);
		node_coords[2*i + 1]  = static_cast<double>(vertex->y);
	}
}/*}}}*/

void ISSM_NUOPC_GetMeshElements(void* model_handle, int* elem_ids, int* elem_types, int* elem_conn){/*{{{*/
	ISSM_NUOPC_Model* model = GetModel(model_handle);

	for(int i=0;i<model->femmodel->elements->Size();i++){
		Element* element = xDynamicCast<Element*>(model->femmodel->elements->GetObjectByOffset(i));
		elem_ids[i]      = element->Sid() + 1;
		elem_types[i]    = ESMC_MESHELEMTYPE_TRI;
		elem_conn[3*i]     = element->vertices[0]->Lid() + 1;
		elem_conn[3*i + 1] = element->vertices[1]->Lid() + 1;
		elem_conn[3*i + 2] = element->vertices[2]->Lid() + 1;
	}
}/*}}}*/

void ISSM_NUOPC_ImportFloatingMelt(void* model_handle, const double* melt_rate, int size){/*{{{*/
	ISSM_NUOPC_Model* model = GetModel(model_handle);
	if(size != model->femmodel->vertices->Size()){
		_error_("ISSM_NUOPC_ImportFloatingMelt received an array with the wrong size");
	}

	IssmDouble* local_values = xNew<IssmDouble>(size);
	for(int i=0;i<size;i++) local_values[i] = static_cast<IssmDouble>(melt_rate[i]);
	InputUpdateFromVectorx(model->femmodel, local_values, BasalforcingsFloatingiceMeltingRateEnum, VertexSIdEnum);
	xDelete<IssmDouble>(local_values);
}/*}}}*/

void ISSM_NUOPC_ExportThickness(void* model_handle, double* thickness, int size){/*{{{*/
	ISSM_NUOPC_Model* model = GetModel(model_handle);
	if(size != model->femmodel->vertices->Size()){
		_error_("ISSM_NUOPC_ExportThickness received an array with the wrong size");
	}
	ExportVertexField(model->femmodel, ThicknessEnum, thickness, size);
}/*}}}*/

void ISSM_NUOPC_ExportSurface(void* model_handle, double* surface, int size){/*{{{*/
	ISSM_NUOPC_Model* model = GetModel(model_handle);
	if(size != model->femmodel->vertices->Size()){
		_error_("ISSM_NUOPC_ExportSurface received an array with the wrong size");
	}
	ExportVertexField(model->femmodel, SurfaceEnum, surface, size);
}/*}}}*/

void ISSM_NUOPC_ExportMask(void* model_handle, double* mask, int size){/*{{{*/
	ISSM_NUOPC_Model* model = GetModel(model_handle);
	if(size != model->femmodel->vertices->Size()){
		_error_("ISSM_NUOPC_ExportMask received an array with the wrong size");
	}
	IssmDouble* mask_levelset = NULL;
	GetVectorFromInputsx(&mask_levelset, model->femmodel, MaskIceLevelsetEnum, VertexSIdEnum);
	for(int i=0;i<size;i++) mask[i] = (mask_levelset[i] <= 0.0 ? 1.0 : 0.0);
	xDelete<IssmDouble>(mask_levelset);
}/*}}}*/

void ISSM_NUOPC_Advance(void* model_handle, double dt_seconds){/*{{{*/
	ISSM_NUOPC_Model* model = GetModel(model_handle);
	IssmDouble start_time;
	IssmDouble final_time;

	model->femmodel->parameters->FindParam(&start_time, TimeEnum);
	final_time = start_time + static_cast<IssmDouble>(dt_seconds);
	model->femmodel->parameters->SetParam(final_time, TimesteppingFinalTimeEnum);
	model->femmodel->Solve();
	model->femmodel->parameters->SetParam(final_time, TimesteppingStartTimeEnum);
}/*}}}*/
