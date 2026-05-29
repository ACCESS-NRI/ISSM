/*!\file VecToMPISerial.cpp
 * \brief gather a Petsc Vector spread across the cluster, onto node 0, and then broadcast to all nodes. 
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../petscincludes.h"
#include "../../../shared/shared.h"

//template<typename doubletype, typename vectype>
//int VecToMPISerialNew(doubletype** pserial_vector, vectype vector,ISSM_MPI_Comm comm,bool broadcast){
//
//	/*Output*/
//	const doubletype *vec_array     = NULL;
//	doubletype       *serial_vector = NULL;
//
//	/*Sequential Vector*/
//	int        n;
//	vectype    vector_seq = NULL;
//	VecScatter ctx        = NULL;
//
//	if(broadcast){
//		VecScatterCreateToAll(vector, &ctx, &vector_seq);
//	}
//	else{
//		VecScatterCreateToZero(vector, &ctx, &vector_seq);
//	}
//
//  /*scatter as many times as you need*/
//  VecScatterBegin(ctx, vector, vector_seq, INSERT_VALUES, SCATTER_FORWARD);
//  VecScatterEnd(  ctx, vector, vector_seq, INSERT_VALUES, SCATTER_FORWARD);
//
//  /*Get pointer to array and copy*/
//  VecGetArrayRead(vector_seq, &vec_array);
//
//  /* Use memcpy to copy data*/
//  VecGetSize(vector_seq, &n);
//  memcpy(serial_vector, vec_array, n*sizeof(doubletype));
//
//  /* Restore and destroy the PETSc Vec array*/
//  VecRestoreArrayRead(vector_seq, &vec_array);
//
//  /* destroy scatter context and local vector when no longer needed*/
//  VecScatterDestroy(&ctx);
//  VecDestroy(&vector_seq);
//
//  /*Assign output pointer*/
//  *pserial_vector = serial_vector;
//}

template<typename doubletype, typename vectype>
int VecToMPISerial(doubletype** pgathered_vector, vectype vector,ISSM_MPI_Comm comm,bool broadcast){

	/*Output*/
	doubletype*        gathered_vector = NULL;
	const PetscScalar* vec_array       = NULL;

	/*Sequential vector and scatter context*/
	vectype    vector_seq = NULL;
	VecScatter ctx        = NULL;

	/*Check for empty vector*/
	int vector_size;
	VecGetSize(vector,&vector_size);
	if(vector_size==0){
		*pgathered_vector=NULL;
		return 1;
	}

	/*Create scatter: ToAll replicates on every rank; ToZero gathers only on rank 0.
	 * VecScatterBegin/End are CUDA-aware: PETSc handles GPU->CPU transfer internally,
	 * avoiding the deprecated VecGetValues path which stalls the GPU pipeline.*/
	if(broadcast){
		VecScatterCreateToAll(vector,&ctx,&vector_seq);
	}
	else{
		VecScatterCreateToZero(vector,&ctx,&vector_seq);
	}

	VecScatterBegin(ctx,vector,vector_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(  ctx,vector,vector_seq,INSERT_VALUES,SCATTER_FORWARD);

	/*Copy data from sequential vector (always CPU-resident after scatter)*/
	int n;
	VecGetSize(vector_seq,&n);
	if(n>0){
		gathered_vector=xNew<doubletype>(n);
		VecGetArrayRead(vector_seq,&vec_array);
		for(int i=0;i<n;i++) gathered_vector[i]=(doubletype)vec_array[i];
		VecRestoreArrayRead(vector_seq,&vec_array);
	}

	/*Destroy scatter context and sequential vector*/
	VecScatterDestroy(&ctx);
	VecDestroy(&vector_seq);

	/*Assign output pointer*/
	*pgathered_vector=gathered_vector;

	return 1;
}

//template int VecToMPISerialNew(IssmDouble** pserial_vector, PVec vector,ISSM_MPI_Comm comm,bool broadcast);
template int VecToMPISerial(IssmDouble** pgathered_vector, PVec vector,ISSM_MPI_Comm comm,bool broadcast);
#if _HAVE_CODIPACK_
//template int VecToMPISerialNew(IssmPDouble** pserial_vector, Vec vector,ISSM_MPI_Comm comm,bool broadcast);
template int VecToMPISerial(IssmPDouble** pgathered_vector, Vec vector,ISSM_MPI_Comm comm,bool broadcast);
#endif
