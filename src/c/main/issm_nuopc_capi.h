#ifndef _ISSM_NUOPC_CAPI_H_
#define _ISSM_NUOPC_CAPI_H_

#include "./issm.h"

#ifdef __cplusplus
extern "C" {
#endif

void* ISSM_NUOPC_CreateFromCase(const char* case_dir, const char* model_name, const char* solution_name, int mpi_comm_f);
void  ISSM_NUOPC_Destroy(void* model_handle);
void  ISSM_NUOPC_WriteRestart(void* model_handle);

void  ISSM_NUOPC_GetMeshCounts(void* model_handle, int* num_nodes, int* num_elements);
void  ISSM_NUOPC_GetMeshNodes(void* model_handle, int* node_ids, int* node_owners, double* node_coords);
void  ISSM_NUOPC_GetMeshElements(void* model_handle, int* elem_ids, int* elem_types, int* elem_conn);

void  ISSM_NUOPC_ImportFloatingMelt(void* model_handle, const double* melt_rate, int size);
void  ISSM_NUOPC_ExportThickness(void* model_handle, double* thickness, int size);
void  ISSM_NUOPC_ExportSurface(void* model_handle, double* surface, int size);
void  ISSM_NUOPC_ExportMask(void* model_handle, double* mask, int size);
void  ISSM_NUOPC_Advance(void* model_handle, double dt_seconds);

#ifdef __cplusplus
}
#endif

#endif
