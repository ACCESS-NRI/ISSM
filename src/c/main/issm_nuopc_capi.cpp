/*!\file: issm_nuopc_capi.cpp
 * \brief: Minimal handle-based ISSM coupling API for an in-process NUOPC cap.
 */

#define _DO_NOT_LOAD_GLOBALS_
#include "./issm_nuopc_capi.h"

#include <string>
#include <vector>

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
	std::string execution_dir(case_dir ? case_dir : ".");

	int solution_type = StringToEnumx(solution_name);
	if(solution_type != TransientSolutionEnum){
		_error_("ISSM_NUOPC first version only supports TransientSolution");
	}

	std::vector<std::string> arg_storage;
	arg_storage.emplace_back("issm_nuopc");
	arg_storage.emplace_back(solution_name ? solution_name : "TransientSolution");
	arg_storage.emplace_back(execution_dir);
	arg_storage.emplace_back(model_name ? model_name : "");

	std::vector<char*> argv;
	argv.reserve(arg_storage.size());
	for(size_t i = 0; i < arg_storage.size(); ++i){
		argv.push_back(const_cast<char*>(arg_storage[i].c_str()));
	}

	ISSM_NUOPC_Model* model = new ISSM_NUOPC_Model();
	model->femmodel = new FemModel(static_cast<int>(argv.size()), argv.data(), comm);
	model->femmodel->parameters->AddObject(new IntParam(IsSlcCouplingEnum,0));

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
	InputUpdateFromVectorx(model->femmodel, local_values, BasalforcingsFloatingiceMeltingRateEnum, VertexLIdEnum);
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
