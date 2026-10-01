module ISSM_NUOPC_CapMod

use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, c_int, c_null_char, c_null_ptr, c_ptr
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

use ESMF, only: ESMF_Clock, ESMF_ClockGet, ESMF_ClockSet, ESMF_COORDSYS_SPH_DEG
use ESMF, only: ESMF_DistGrid, ESMF_DistGridCreate
use ESMF, only: ESMF_Field, ESMF_FieldCreate, ESMF_FieldGet
use ESMF, only: ESMF_Grid, ESMF_GridCreate
use ESMF, only: ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSetEntryPoint
use ESMF, only: ESMF_LOGMSG_ERROR, ESMF_LOGMSG_INFO, ESMF_LogWrite
use ESMF, only: ESMF_Mesh, ESMF_MeshCreate
use ESMF, only: ESMF_MESHLOC_ELEMENT, ESMF_METHOD_INITIALIZE
use ESMF, only: ESMF_State, ESMF_StateGet
use ESMF, only: ESMF_SUCCESS, ESMF_TimeInterval, ESMF_TimeIntervalGet, ESMF_TimeIntervalSet
use ESMF, only: ESMF_KIND_R8, ESMF_TYPEKIND_R8, ESMF_VM, ESMF_VMGet, ESMF_VMGetCurrent

use NUOPC, only: NUOPC_AddNestedState, NUOPC_Advertise, NUOPC_CompAttributeSet, NUOPC_CompDerive
use NUOPC, only: NUOPC_CompFilterPhaseMap, NUOPC_CompAttributeGet
use NUOPC, only: NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
use NUOPC, only: NUOPC_IsConnected, NUOPC_Realize, NUOPC_SetAttribute
use NUOPC_Model, only: NUOPC_ModelGet
use NUOPC_Model, only: model_label_Advance        => label_Advance
use NUOPC_Model, only: model_label_CheckImport    => label_CheckImport
use NUOPC_Model, only: model_label_DataInitialize => label_DataInitialize
use NUOPC_Model, only: model_label_Finalize       => label_Finalize
use NUOPC_Model, only: model_routine_SS           => SetServices

implicit none
private

public :: SetServices

integer, parameter :: coord_dim = 2
integer, parameter :: nodes_per_element = 3
character(len=*), parameter :: glc_cplset_name = 'GLC1'
character(len=*), parameter :: default_scalar_field_name = 'cpl_scalars'
character(len=*), parameter :: import_melt_name = 'IceSheetBasalMeltRate'
character(len=*), parameter :: import_melt_stdname = 'IceSheetBasalMeltRate'
real(c_double), parameter :: max_valid_melt_rate = 1.0e-3_c_double
character(len=*), parameter :: export_thickness_name = 'iceThickness'
character(len=*), parameter :: export_thickness_stdname = 'IceSheetThickness'
character(len=*), parameter :: export_surface_name = 'iceSurface'
character(len=*), parameter :: export_surface_stdname = 'IceSheetSurfaceElevation'
character(len=*), parameter :: export_mask_name = 'iceMask'
character(len=*), parameter :: export_mask_stdname = 'IceSheetMask'
character(len=*), parameter :: export_sg_area_name = 'Sg_area'
character(len=*), parameter :: export_sg_area_stdname = 'Sg_area'
character(len=*), parameter :: export_sg_icemask_name = 'Sg_icemask'
character(len=*), parameter :: export_sg_icemask_stdname = 'Sg_icemask'
character(len=*), parameter :: export_sg_icemask_fluxes_name = 'Sg_icemask_coupled_fluxes'
character(len=*), parameter :: export_sg_icemask_fluxes_stdname = 'Sg_icemask_coupled_fluxes'
character(len=*), parameter :: export_sg_covered_name = 'Sg_ice_covered'
character(len=*), parameter :: export_sg_covered_stdname = 'Sg_ice_covered'
character(len=*), parameter :: export_sg_topo_name = 'Sg_topo'
character(len=*), parameter :: export_sg_topo_stdname = 'Sg_topo'
character(len=*), parameter :: export_flgg_hflx_name = 'Flgg_hflx'
character(len=*), parameter :: export_flgg_hflx_stdname = 'Flgg_hflx'
character(len=*), parameter :: export_fgrg_rofl_name = 'Fgrg_rofl'
character(len=*), parameter :: export_fgrg_rofl_stdname = 'Fgrg_rofl'
character(len=*), parameter :: export_fgrg_rofi_name = 'Fgrg_rofi'
character(len=*), parameter :: export_fgrg_rofi_stdname = 'Fgrg_rofi'

type issm_cap_state_type
  type(c_ptr) :: handle = c_null_ptr
  type(ESMF_Mesh) :: mesh
  type(ESMF_State) :: glc_import_state
  type(ESMF_State) :: glc_export_state
  integer(c_int) :: mpi_comm = 0_c_int
  integer(c_int) :: node_count = 0_c_int
  integer(c_int) :: element_count = 0_c_int
  integer :: scalar_field_count = 0
  integer :: scalar_idx_grid_nx = 0
  integer :: scalar_idx_grid_ny = 0
  integer :: scalar_idx_next_sw_cday = 0
  integer :: scalar_idx_precip_factor = 0
  logical :: mesh_created = .false.
  logical :: write_restart = .false.
  logical :: require_melt_import = .true.
  logical :: melt_diagnostics = .false.
  integer :: advance_count = 0
  integer :: melt_diagnostics_interval = 1
  integer :: issm_advance_seconds = 0
  character(len=64) :: scalar_field_name = default_scalar_field_name
  integer(c_int), allocatable :: node_ids(:)
  integer(c_int), allocatable :: node_owners(:)
  integer(c_int), allocatable :: element_ids(:)
  integer(c_int), allocatable :: element_types(:)
  integer(c_int), allocatable :: element_conn(:)
  real(c_double), allocatable :: node_coords(:)
  real(c_double), allocatable :: element_coords(:)
  real(c_double), allocatable :: import_melt(:)
  real(c_double), allocatable :: import_melt_elements(:)
  real(c_double), allocatable :: retained_melt(:)
  real(c_double), allocatable :: export_thickness(:)
  real(c_double), allocatable :: export_surface(:)
  real(c_double), allocatable :: export_mask(:)
  real(c_double), allocatable :: export_thickness_elements(:)
  real(c_double), allocatable :: export_surface_elements(:)
  real(c_double), allocatable :: export_mask_elements(:)
end type issm_cap_state_type

type(issm_cap_state_type), save :: cap_state

interface
  function ISSM_NUOPC_CreateFromCase(case_dir, model_name, solution_name, mpi_comm_f) &
    bind(C, name='ISSM_NUOPC_CreateFromCase') result(handle)
    import :: c_char, c_int, c_ptr
    character(kind=c_char), intent(in) :: case_dir(*)
    character(kind=c_char), intent(in) :: model_name(*)
    character(kind=c_char), intent(in) :: solution_name(*)
    integer(c_int), value, intent(in) :: mpi_comm_f
    type(c_ptr) :: handle
  end function ISSM_NUOPC_CreateFromCase

  subroutine ISSM_NUOPC_Destroy(handle) bind(C, name='ISSM_NUOPC_Destroy')
    import :: c_ptr
    type(c_ptr), value, intent(in) :: handle
  end subroutine ISSM_NUOPC_Destroy

  subroutine ISSM_NUOPC_WriteRestart(handle) bind(C, name='ISSM_NUOPC_WriteRestart')
    import :: c_ptr
    type(c_ptr), value, intent(in) :: handle
  end subroutine ISSM_NUOPC_WriteRestart

  subroutine ISSM_NUOPC_GetMeshCounts(handle, num_nodes, num_elements) bind(C, name='ISSM_NUOPC_GetMeshCounts')
    import :: c_ptr, c_int
    type(c_ptr), value, intent(in) :: handle
    integer(c_int), intent(out) :: num_nodes
    integer(c_int), intent(out) :: num_elements
  end subroutine ISSM_NUOPC_GetMeshCounts

  subroutine ISSM_NUOPC_GetMeshNodes(handle, node_ids, node_owners, node_coords) bind(C, name='ISSM_NUOPC_GetMeshNodes')
    import :: c_ptr, c_int, c_double
    type(c_ptr), value, intent(in) :: handle
    integer(c_int), intent(out) :: node_ids(*)
    integer(c_int), intent(out) :: node_owners(*)
    real(c_double), intent(out) :: node_coords(*)
  end subroutine ISSM_NUOPC_GetMeshNodes

  subroutine ISSM_NUOPC_GetMeshElements(handle, elem_ids, elem_types, elem_conn) &
    bind(C, name='ISSM_NUOPC_GetMeshElements')
    import :: c_ptr, c_int
    type(c_ptr), value, intent(in) :: handle
    integer(c_int), intent(out) :: elem_ids(*)
    integer(c_int), intent(out) :: elem_types(*)
    integer(c_int), intent(out) :: elem_conn(*)
  end subroutine ISSM_NUOPC_GetMeshElements

  subroutine ISSM_NUOPC_ImportFloatingMelt(handle, melt_rate, size) bind(C, name='ISSM_NUOPC_ImportFloatingMelt')
    import :: c_ptr, c_double, c_int
    type(c_ptr), value, intent(in) :: handle
    real(c_double), intent(in) :: melt_rate(*)
    integer(c_int), value, intent(in) :: size
  end subroutine ISSM_NUOPC_ImportFloatingMelt

  subroutine ISSM_NUOPC_RequireExternalFloatingMelt(handle) bind(C, name='ISSM_NUOPC_RequireExternalFloatingMelt')
    import :: c_ptr
    type(c_ptr), value, intent(in) :: handle
  end subroutine ISSM_NUOPC_RequireExternalFloatingMelt

  subroutine ISSM_NUOPC_ExportFloatingMelt(handle, melt_rate, size) bind(C, name='ISSM_NUOPC_ExportFloatingMelt')
    import :: c_ptr, c_double, c_int
    type(c_ptr), value, intent(in) :: handle
    real(c_double), intent(out) :: melt_rate(*)
    integer(c_int), value, intent(in) :: size
  end subroutine ISSM_NUOPC_ExportFloatingMelt

  subroutine ISSM_NUOPC_ExportThickness(handle, thickness, size) bind(C, name='ISSM_NUOPC_ExportThickness')
    import :: c_ptr, c_double, c_int
    type(c_ptr), value, intent(in) :: handle
    real(c_double), intent(out) :: thickness(*)
    integer(c_int), value, intent(in) :: size
  end subroutine ISSM_NUOPC_ExportThickness

  subroutine ISSM_NUOPC_ExportSurface(handle, surface, size) bind(C, name='ISSM_NUOPC_ExportSurface')
    import :: c_ptr, c_double, c_int
    type(c_ptr), value, intent(in) :: handle
    real(c_double), intent(out) :: surface(*)
    integer(c_int), value, intent(in) :: size
  end subroutine ISSM_NUOPC_ExportSurface

  subroutine ISSM_NUOPC_ExportMask(handle, mask, size) bind(C, name='ISSM_NUOPC_ExportMask')
    import :: c_ptr, c_double, c_int
    type(c_ptr), value, intent(in) :: handle
    real(c_double), intent(out) :: mask(*)
    integer(c_int), value, intent(in) :: size
  end subroutine ISSM_NUOPC_ExportMask

  subroutine ISSM_NUOPC_Advance(handle, dt_seconds) bind(C, name='ISSM_NUOPC_Advance')
    import :: c_ptr, c_double
    type(c_ptr), value, intent(in) :: handle
    real(c_double), value, intent(in) :: dt_seconds
  end subroutine ISSM_NUOPC_Advance
end interface

contains

subroutine SetServices(gcomp, rc)
  type(ESMF_GridComp) :: gcomp
  integer, intent(out) :: rc

  rc = ESMF_SUCCESS

  call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
  if (rc /= ESMF_SUCCESS) return

  call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, userRoutine=InitializeP0, phase=0, rc=rc)
  if (rc /= ESMF_SUCCESS) return

  call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/'IPDv03p1'/), userRoutine=InitializeAdvertise, rc=rc)
  if (rc /= ESMF_SUCCESS) return

  call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/'IPDv03p3'/), userRoutine=InitializeRealize, rc=rc)
  if (rc /= ESMF_SUCCESS) return

  call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, specRoutine=DataInitialize, rc=rc)
  if (rc /= ESMF_SUCCESS) return

  call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, specRoutine=CheckImportNoOp, rc=rc)
  if (rc /= ESMF_SUCCESS) return

  call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, specRoutine=ModelAdvance, rc=rc)
  if (rc /= ESMF_SUCCESS) return

  call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, specRoutine=FinalizeModel, rc=rc)
end subroutine SetServices

subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
  type(ESMF_GridComp) :: gcomp
  type(ESMF_State) :: importState
  type(ESMF_State) :: exportState
  type(ESMF_Clock) :: clock
  integer, intent(out) :: rc

  rc = ESMF_SUCCESS
  call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/'IPDv03p'/), rc=rc)
end subroutine InitializeP0

subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
  type(ESMF_GridComp) :: gcomp
  type(ESMF_State) :: importState
  type(ESMF_State) :: exportState
  type(ESMF_Clock) :: clock
  integer, intent(out) :: rc

  type(ESMF_VM) :: vm
  logical :: is_present
  logical :: is_set
  integer :: mpi_comm_local
  character(len=512) :: message
  character(len=256) :: case_dir
  character(len=256) :: model_name
  character(len=256) :: solution_name
  character(len=32) :: value

  rc = ESMF_SUCCESS
  case_dir = ''
  model_name = ''
  solution_name = 'TransientSolution'
  mpi_comm_local = 0

  call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call ESMF_VMGet(vm, mpiCommunicator=mpi_comm_local, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  cap_state%mpi_comm = int(mpi_comm_local, c_int)

  call RequireAttribute(gcomp, 'case_dir', case_dir, rc)
  if (rc /= ESMF_SUCCESS) return
  call RequireAttribute(gcomp, 'model_name', model_name, rc)
  if (rc /= ESMF_SUCCESS) return

  call NUOPC_CompAttributeGet(gcomp, name='solution_name', value=solution_name, &
    isPresent=is_present, isSet=is_set, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  if (.not. is_present .or. .not. is_set) solution_name = 'TransientSolution'

  call NUOPC_CompAttributeGet(gcomp, name='write_restart', value=value, isPresent=is_present, isSet=is_set, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  cap_state%write_restart = is_present .and. is_set .and. IsTrue(value)

  call NUOPC_CompAttributeGet(gcomp, name='require_melt_import', value=value, isPresent=is_present, isSet=is_set, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  if (is_present .and. is_set) cap_state%require_melt_import = IsTrue(value)

  call NUOPC_CompAttributeGet(gcomp, name='melt_diagnostics', value=value, isPresent=is_present, isSet=is_set, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  if (is_present .and. is_set) cap_state%melt_diagnostics = IsTrue(value)
  call ReadOptionalIntegerAttribute(gcomp, 'melt_diagnostics_interval', &
    cap_state%melt_diagnostics_interval, rc)
  if (rc /= ESMF_SUCCESS) return
  if (cap_state%melt_diagnostics_interval <= 0) then
    rc = 1
    call LogError('ISSM_NUOPC: melt_diagnostics_interval must be positive')
    return
  end if
  call ReadOptionalIntegerAttribute(gcomp, 'issm_advance_seconds', cap_state%issm_advance_seconds, rc)
  if (rc /= ESMF_SUCCESS) return
  if (cap_state%issm_advance_seconds < 0) then
    rc = 1
    call LogError('ISSM_NUOPC: issm_advance_seconds must be positive when set')
    return
  end if

  call NUOPC_CompAttributeGet(gcomp, name='ScalarFieldName', value=cap_state%scalar_field_name, &
    isPresent=is_present, isSet=is_set, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  if (.not. is_present .or. .not. is_set) cap_state%scalar_field_name = default_scalar_field_name
  call ReadOptionalIntegerAttribute(gcomp, 'ScalarFieldCount', cap_state%scalar_field_count, rc)
  if (rc /= ESMF_SUCCESS) return
  call ReadOptionalIntegerAttribute(gcomp, 'ScalarFieldIdxGridNX', cap_state%scalar_idx_grid_nx, rc)
  if (rc /= ESMF_SUCCESS) return
  call ReadOptionalIntegerAttribute(gcomp, 'ScalarFieldIdxGridNY', cap_state%scalar_idx_grid_ny, rc)
  if (rc /= ESMF_SUCCESS) return
  call ReadOptionalIntegerAttribute(gcomp, 'ScalarFieldIdxNextSwCday', cap_state%scalar_idx_next_sw_cday, rc)
  if (rc /= ESMF_SUCCESS) return
  call ReadOptionalIntegerAttribute(gcomp, 'ScalarFieldIdxPrecipFactor', cap_state%scalar_idx_precip_factor, rc)
  if (rc /= ESMF_SUCCESS) return

  cap_state%handle = ISSM_NUOPC_CreateFromCase( &
    trim(case_dir)//c_null_char, &
    trim(model_name)//c_null_char, &
    trim(solution_name)//c_null_char, &
    cap_state%mpi_comm)
  if (.not. c_associated(cap_state%handle)) then
    call ESMF_LogWrite('ISSM_NUOPC: failed to create ISSM model handle', ESMF_LOGMSG_ERROR, rc=rc)
    return
  end if
  if (cap_state%require_melt_import) then
    call ISSM_NUOPC_RequireExternalFloatingMelt(cap_state%handle)
    call LogInfo('ISSM_NUOPC: validated generic prescribed floating melt forcing')
  end if

  call ISSM_NUOPC_GetMeshCounts(cap_state%handle, cap_state%node_count, cap_state%element_count)
  call AllocateCapState(rc)
  if (rc /= ESMF_SUCCESS) return

  write(message,'(a,a,a,a,a,i0,a,i0,a,l1)') &
    'ISSM_NUOPC: initialized model ', trim(model_name), &
    ' solution ', trim(solution_name), &
    ' nodes=', int(cap_state%node_count), &
    ' elements=', int(cap_state%element_count), &
    ' require_melt_import=', cap_state%require_melt_import
  call LogInfo(trim(message))

  if (cap_state%scalar_field_count > 0) then
    call NUOPC_AddNestedState(importState, CplSet=glc_cplset_name, &
      nestedState=cap_state%glc_import_state, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    call NUOPC_AddNestedState(exportState, CplSet=glc_cplset_name, &
      nestedState=cap_state%glc_export_state, rc=rc)
    if (rc /= ESMF_SUCCESS) return
  else
    cap_state%glc_import_state = importState
    cap_state%glc_export_state = exportState
  end if

  if (cap_state%scalar_field_count > 0) then
    call NUOPC_Advertise(cap_state%glc_import_state, standardName=trim(cap_state%scalar_field_name), &
      name=trim(cap_state%scalar_field_name), rc=rc)
    if (rc /= ESMF_SUCCESS) return
    call NUOPC_Advertise(cap_state%glc_export_state, standardName=trim(cap_state%scalar_field_name), &
      name=trim(cap_state%scalar_field_name), rc=rc)
    if (rc /= ESMF_SUCCESS) return
  end if

  if (cap_state%require_melt_import) then
    call NUOPC_Advertise( &
      cap_state%glc_import_state, &
      standardName=import_melt_stdname, &
      name=import_melt_name, &
      TransferOfferGeomObject='will provide', &
      SharePolicyField='share', &
      SharePolicyGeomObject='share', &
      rc=rc &
    )
    if (rc /= ESMF_SUCCESS) return
  end if
  call NUOPC_Advertise(cap_state%glc_export_state, standardName=export_thickness_stdname, &
    name=export_thickness_name, TransferOfferGeomObject='will provide', &
    SharePolicyField='share', SharePolicyGeomObject='share', rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call NUOPC_Advertise(cap_state%glc_export_state, standardName=export_surface_stdname, &
    name=export_surface_name, TransferOfferGeomObject='will provide', &
    SharePolicyField='share', SharePolicyGeomObject='share', rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call NUOPC_Advertise(cap_state%glc_export_state, standardName=export_mask_stdname, &
    name=export_mask_name, TransferOfferGeomObject='will provide', &
    SharePolicyField='share', SharePolicyGeomObject='share', rc=rc)
  if (rc /= ESMF_SUCCESS) return
  if (cap_state%scalar_field_count > 0) then
    call NUOPC_Advertise(cap_state%glc_export_state, standardName=export_sg_area_stdname, &
      name=export_sg_area_name, TransferOfferGeomObject='will provide', &
      SharePolicyField='share', SharePolicyGeomObject='share', rc=rc)
    if (rc /= ESMF_SUCCESS) return
    call NUOPC_Advertise(cap_state%glc_export_state, standardName=export_sg_icemask_stdname, &
      name=export_sg_icemask_name, TransferOfferGeomObject='will provide', &
      SharePolicyField='share', SharePolicyGeomObject='share', rc=rc)
    if (rc /= ESMF_SUCCESS) return
    call NUOPC_Advertise(cap_state%glc_export_state, standardName=export_sg_icemask_fluxes_stdname, &
      name=export_sg_icemask_fluxes_name, TransferOfferGeomObject='will provide', &
      SharePolicyField='share', SharePolicyGeomObject='share', rc=rc)
    if (rc /= ESMF_SUCCESS) return
    call NUOPC_Advertise(cap_state%glc_export_state, standardName=export_sg_covered_stdname, &
      name=export_sg_covered_name, TransferOfferGeomObject='will provide', &
      SharePolicyField='share', SharePolicyGeomObject='share', rc=rc)
    if (rc /= ESMF_SUCCESS) return
    call NUOPC_Advertise(cap_state%glc_export_state, standardName=export_sg_topo_stdname, &
      name=export_sg_topo_name, TransferOfferGeomObject='will provide', &
      SharePolicyField='share', SharePolicyGeomObject='share', rc=rc)
    if (rc /= ESMF_SUCCESS) return
    call NUOPC_Advertise(cap_state%glc_export_state, standardName=export_flgg_hflx_stdname, &
      name=export_flgg_hflx_name, TransferOfferGeomObject='will provide', &
      SharePolicyField='share', SharePolicyGeomObject='share', rc=rc)
    if (rc /= ESMF_SUCCESS) return
    call NUOPC_Advertise(cap_state%glc_export_state, standardName=export_fgrg_rofl_stdname, &
      name=export_fgrg_rofl_name, TransferOfferGeomObject='will provide', &
      SharePolicyField='share', SharePolicyGeomObject='share', rc=rc)
    if (rc /= ESMF_SUCCESS) return
    call NUOPC_Advertise(cap_state%glc_export_state, standardName=export_fgrg_rofi_stdname, &
      name=export_fgrg_rofi_name, TransferOfferGeomObject='will provide', &
      SharePolicyField='share', SharePolicyGeomObject='share', rc=rc)
    if (rc /= ESMF_SUCCESS) return
  end if
end subroutine InitializeAdvertise

subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
  type(ESMF_GridComp) :: gcomp
  type(ESMF_State) :: importState
  type(ESMF_State) :: exportState
  type(ESMF_Clock) :: clock
  integer, intent(out) :: rc

  rc = ESMF_SUCCESS

  call ISSM_NUOPC_GetMeshNodes(cap_state%handle, cap_state%node_ids, cap_state%node_owners, &
    cap_state%node_coords)
  call ISSM_NUOPC_GetMeshElements(cap_state%handle, cap_state%element_ids, &
    cap_state%element_types, cap_state%element_conn)
  ! The MOM mesh uses MISMIP km coordinates with degree metadata; expose matching
  ! coordinates to ESMF while keeping ISSM's internal model state unchanged.
  cap_state%node_coords(:) = cap_state%node_coords(:) / 1000.0_c_double
  call BuildElementCoordinates(rc)
  if (rc /= ESMF_SUCCESS) return

  cap_state%mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, coordSys=ESMF_COORDSYS_SPH_DEG, &
       nodeIds=cap_state%node_ids, nodeCoords=cap_state%node_coords, &
       nodeOwners=cap_state%node_owners, elementIds=cap_state%element_ids, &
       elementTypes=cap_state%element_types, elementConn=cap_state%element_conn, &
       elementCoords=cap_state%element_coords, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  cap_state%mesh_created = .true.

  if (cap_state%scalar_field_count > 0) then
    call RealizeScalarField(cap_state%glc_import_state, rc)
    if (rc /= ESMF_SUCCESS) return
    call RealizeScalarField(cap_state%glc_export_state, rc)
    if (rc /= ESMF_SUCCESS) return
  end if

  if (cap_state%require_melt_import) then
    call RealizeField(cap_state%glc_import_state, import_melt_name, cap_state%mesh, rc)
    if (rc /= ESMF_SUCCESS) return
  end if
  call RealizeField(cap_state%glc_export_state, export_thickness_name, cap_state%mesh, rc)
  if (rc /= ESMF_SUCCESS) return
  call RealizeField(cap_state%glc_export_state, export_surface_name, cap_state%mesh, rc)
  if (rc /= ESMF_SUCCESS) return
  call RealizeField(cap_state%glc_export_state, export_mask_name, cap_state%mesh, rc)
  if (rc /= ESMF_SUCCESS) return
  if (cap_state%scalar_field_count > 0) then
    call RealizeField(cap_state%glc_export_state, export_sg_area_name, cap_state%mesh, rc)
    if (rc /= ESMF_SUCCESS) return
    call RealizeField(cap_state%glc_export_state, export_sg_icemask_name, cap_state%mesh, rc)
    if (rc /= ESMF_SUCCESS) return
    call RealizeField(cap_state%glc_export_state, export_sg_icemask_fluxes_name, cap_state%mesh, rc)
    if (rc /= ESMF_SUCCESS) return
    call RealizeField(cap_state%glc_export_state, export_sg_covered_name, cap_state%mesh, rc)
    if (rc /= ESMF_SUCCESS) return
    call RealizeField(cap_state%glc_export_state, export_sg_topo_name, cap_state%mesh, rc)
    if (rc /= ESMF_SUCCESS) return
    call RealizeField(cap_state%glc_export_state, export_flgg_hflx_name, cap_state%mesh, rc)
    if (rc /= ESMF_SUCCESS) return
    call RealizeField(cap_state%glc_export_state, export_fgrg_rofl_name, cap_state%mesh, rc)
    if (rc /= ESMF_SUCCESS) return
    call RealizeField(cap_state%glc_export_state, export_fgrg_rofi_name, cap_state%mesh, rc)
    if (rc /= ESMF_SUCCESS) return
  end if

  if (cap_state%scalar_field_count > 0) then
    call SetExportScalars(cap_state%glc_export_state, rc)
    if (rc /= ESMF_SUCCESS) return
  end if

  call RefreshExports(cap_state%glc_export_state, rc)
end subroutine InitializeRealize

subroutine DataInitialize(gcomp, rc)
  type(ESMF_GridComp) :: gcomp
  integer, intent(out) :: rc

  type(ESMF_Clock) :: clock
  type(ESMF_State) :: exportState

  rc = ESMF_SUCCESS
  call ESMF_GridCompGet(gcomp, clock=clock, exportState=exportState, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call ApplyAdvanceSecondsClock(gcomp, clock, rc)
  if (rc /= ESMF_SUCCESS) return
  call RefreshExports(cap_state%glc_export_state, rc)
  if (rc /= ESMF_SUCCESS) return
  call LogExportSummary('after initialize', rc)
  if (rc /= ESMF_SUCCESS) return
  call NUOPC_CompAttributeSet(gcomp, name='InitializeDataComplete', value='true', rc=rc)
end subroutine DataInitialize

subroutine ApplyAdvanceSecondsClock(gcomp, clock, rc)
  type(ESMF_GridComp) :: gcomp
  type(ESMF_Clock), intent(inout) :: clock
  integer, intent(out) :: rc

  type(ESMF_TimeInterval) :: time_step
  character(len=64) :: advance_seconds_text
  integer :: advance_seconds
  integer :: log_rc
  integer :: read_status
  logical :: is_present
  logical :: is_set

  rc = ESMF_SUCCESS

  call NUOPC_CompAttributeGet(gcomp, name='advance_seconds', value=advance_seconds_text, &
    isPresent=is_present, isSet=is_set, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  if (.not. is_present .or. .not. is_set) return

  read(advance_seconds_text, *, iostat=read_status) advance_seconds
  if (read_status /= 0 .or. advance_seconds <= 0) then
    rc = 1
    call ESMF_LogWrite('ISSM_NUOPC: invalid advance_seconds attribute', ESMF_LOGMSG_ERROR, rc=log_rc)
    return
  end if

  call ESMF_TimeIntervalSet(time_step, s=advance_seconds, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call ESMF_ClockSet(clock, timeStep=time_step, rc=rc)
end subroutine ApplyAdvanceSecondsClock

subroutine CheckImportNoOp(gcomp, rc)
  type(ESMF_GridComp) :: gcomp
  integer, intent(out) :: rc

  rc = ESMF_SUCCESS
end subroutine CheckImportNoOp

subroutine ModelAdvance(gcomp, rc)
  type(ESMF_GridComp) :: gcomp
  integer, intent(out) :: rc

  type(ESMF_Clock) :: clock
  type(ESMF_State) :: importState
  type(ESMF_State) :: exportState
  type(ESMF_TimeInterval) :: timeStep
  type(ESMF_Field) :: field
  real(ESMF_KIND_R8), pointer :: field_ptr(:)
  character(len=128) :: context
  character(len=256) :: message
  character(len=64) :: advance_seconds_text
  integer :: clock_dt_seconds
  integer :: dt_seconds
  integer :: attr_dt_seconds
  integer :: log_rc
  integer :: read_status
  logical :: is_present
  logical :: is_set

  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, exportState=exportState, rc=rc)
  if (rc /= ESMF_SUCCESS) return

  call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call ESMF_TimeIntervalGet(timeStep, s=clock_dt_seconds, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  dt_seconds = clock_dt_seconds

  call NUOPC_CompAttributeGet(gcomp, name='advance_seconds', value=advance_seconds_text, &
    isPresent=is_present, isSet=is_set, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  if (is_present .and. is_set) then
    read(advance_seconds_text, *, iostat=read_status) attr_dt_seconds
    if (read_status /= 0 .or. attr_dt_seconds <= 0) then
      rc = 1
      call ESMF_LogWrite('ISSM_NUOPC: invalid advance_seconds attribute', ESMF_LOGMSG_ERROR, rc=log_rc)
      return
    end if
    dt_seconds = attr_dt_seconds
  end if
  if (cap_state%issm_advance_seconds > 0) dt_seconds = cap_state%issm_advance_seconds

  if (cap_state%require_melt_import) then
    call ESMF_StateGet(cap_state%glc_import_state, itemName=import_melt_name, field=field, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    call ESMF_FieldGet(field, farrayPtr=field_ptr, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    if (size(field_ptr) /= int(cap_state%element_count)) then
      rc = 1
      call LogError('ISSM_NUOPC: IceSheetBasalMeltRate import has unexpected element count')
      return
    end if
    cap_state%import_melt_elements(:) = field_ptr(:)
    call ValidateMeltValues('IceSheetBasalMeltRate element import', cap_state%import_melt_elements, rc)
    if (rc /= ESMF_SUCCESS) return
    call AverageElementToNode(cap_state%import_melt_elements, cap_state%import_melt, rc)
    if (rc /= ESMF_SUCCESS) return
    call ValidateMeltValues('IceSheetBasalMeltRate nodal import', cap_state%import_melt, rc)
    if (rc /= ESMF_SUCCESS) return
  else
    cap_state%import_melt(:) = 0.0_c_double
  end if
  call ISSM_NUOPC_ImportFloatingMelt(cap_state%handle, cap_state%import_melt, cap_state%node_count)

  cap_state%advance_count = cap_state%advance_count + 1
  if (ShouldLogMeltDiagnostics()) then
    if (cap_state%require_melt_import) then
      call LogMeltSummary('received element melt', cap_state%import_melt_elements, rc)
      if (rc /= ESMF_SUCCESS) return
    end if
    call LogMeltSummary('applied nodal melt before advance', cap_state%import_melt, rc)
    if (rc /= ESMF_SUCCESS) return
  end if

  write(message,'(a,i0,a,i0,a,i0,a,l1)') &
    'ISSM_NUOPC: advance ', cap_state%advance_count, &
    ' clock_dt_seconds=', clock_dt_seconds, &
    ' issm_dt_seconds=', dt_seconds, &
    ' require_melt_import=', cap_state%require_melt_import
  call LogInfo(trim(message))

  call ISSM_NUOPC_Advance(cap_state%handle, real(dt_seconds, c_double))
  if (cap_state%require_melt_import) then
    call ISSM_NUOPC_ExportFloatingMelt(cap_state%handle, cap_state%retained_melt, cap_state%node_count)
    call VerifyRetainedMelt(rc)
    if (rc /= ESMF_SUCCESS) return
    if (ShouldLogMeltDiagnostics()) then
      call LogMeltSummary('retained nodal melt after advance', cap_state%retained_melt, rc)
      if (rc /= ESMF_SUCCESS) return
    end if
  end if
  call RefreshExports(cap_state%glc_export_state, rc)
  if (rc /= ESMF_SUCCESS) return
  write(context,'(a,i0)') 'after advance ', cap_state%advance_count
  call LogExportSummary(trim(context), rc)
end subroutine ModelAdvance

subroutine FinalizeModel(gcomp, rc)
  type(ESMF_GridComp) :: gcomp
  integer, intent(out) :: rc

  character(len=128) :: message

  rc = ESMF_SUCCESS

  write(message,'(a,i0)') 'ISSM_NUOPC: finalize after advances=', cap_state%advance_count
  call LogInfo(trim(message))

  if (c_associated(cap_state%handle)) then
    if (cap_state%write_restart) call ISSM_NUOPC_WriteRestart(cap_state%handle)
    call ISSM_NUOPC_Destroy(cap_state%handle)
  end if

  call ResetCapState()
end subroutine FinalizeModel

subroutine RequireAttribute(gcomp, name, value, rc)
  type(ESMF_GridComp), intent(in) :: gcomp
  character(len=*), intent(in) :: name
  character(len=*), intent(out) :: value
  integer, intent(out) :: rc

  logical :: is_present
  logical :: is_set

  value = ''
  rc = ESMF_SUCCESS
  call NUOPC_CompAttributeGet(gcomp, name=trim(name), value=value, isPresent=is_present, isSet=is_set, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  if (.not. is_present .or. .not. is_set) then
    rc = 1
    call ESMF_LogWrite('ISSM_NUOPC: missing required component attribute '//trim(name), ESMF_LOGMSG_ERROR, rc=rc)
  end if
end subroutine RequireAttribute

subroutine ReadOptionalIntegerAttribute(gcomp, name, target, rc)
  type(ESMF_GridComp), intent(in) :: gcomp
  character(len=*), intent(in) :: name
  integer, intent(inout) :: target
  integer, intent(out) :: rc

  character(len=64) :: value
  integer :: log_rc
  integer :: read_status
  logical :: is_present
  logical :: is_set

  rc = ESMF_SUCCESS
  value = ''
  call NUOPC_CompAttributeGet(gcomp, name=trim(name), value=value, isPresent=is_present, isSet=is_set, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  if (is_present .and. is_set) then
    read(value, *, iostat=read_status) target
    if (read_status /= 0) then
      rc = 1
      call ESMF_LogWrite('ISSM_NUOPC: invalid integer component attribute '//trim(name), &
        ESMF_LOGMSG_ERROR, rc=log_rc)
    end if
  end if
end subroutine ReadOptionalIntegerAttribute

subroutine AllocateCapState(rc)
  integer, intent(out) :: rc

  rc = ESMF_SUCCESS

  if (allocated(cap_state%node_ids)) call ResetCapState()

  allocate(cap_state%node_ids(cap_state%node_count))
  allocate(cap_state%node_owners(cap_state%node_count))
  allocate(cap_state%node_coords(coord_dim * cap_state%node_count))
  allocate(cap_state%element_ids(cap_state%element_count))
  allocate(cap_state%element_types(cap_state%element_count))
  allocate(cap_state%element_conn(nodes_per_element * cap_state%element_count))
  allocate(cap_state%element_coords(coord_dim * cap_state%element_count))
  allocate(cap_state%import_melt(cap_state%node_count))
  allocate(cap_state%import_melt_elements(cap_state%element_count))
  allocate(cap_state%retained_melt(cap_state%node_count))
  allocate(cap_state%export_thickness(cap_state%node_count))
  allocate(cap_state%export_surface(cap_state%node_count))
  allocate(cap_state%export_mask(cap_state%node_count))
  allocate(cap_state%export_thickness_elements(cap_state%element_count))
  allocate(cap_state%export_surface_elements(cap_state%element_count))
  allocate(cap_state%export_mask_elements(cap_state%element_count))

  cap_state%import_melt = 0.0_c_double
  cap_state%element_coords = 0.0_c_double
  cap_state%import_melt_elements = 0.0_c_double
  cap_state%retained_melt = 0.0_c_double
  cap_state%export_thickness = 0.0_c_double
  cap_state%export_surface = 0.0_c_double
  cap_state%export_mask = 0.0_c_double
  cap_state%export_thickness_elements = 0.0_c_double
  cap_state%export_surface_elements = 0.0_c_double
  cap_state%export_mask_elements = 0.0_c_double
end subroutine AllocateCapState

subroutine ResetCapState()
  if (allocated(cap_state%node_ids)) deallocate(cap_state%node_ids)
  if (allocated(cap_state%node_owners)) deallocate(cap_state%node_owners)
  if (allocated(cap_state%node_coords)) deallocate(cap_state%node_coords)
  if (allocated(cap_state%element_ids)) deallocate(cap_state%element_ids)
  if (allocated(cap_state%element_types)) deallocate(cap_state%element_types)
  if (allocated(cap_state%element_conn)) deallocate(cap_state%element_conn)
  if (allocated(cap_state%element_coords)) deallocate(cap_state%element_coords)
  if (allocated(cap_state%import_melt)) deallocate(cap_state%import_melt)
  if (allocated(cap_state%import_melt_elements)) deallocate(cap_state%import_melt_elements)
  if (allocated(cap_state%retained_melt)) deallocate(cap_state%retained_melt)
  if (allocated(cap_state%export_thickness)) deallocate(cap_state%export_thickness)
  if (allocated(cap_state%export_surface)) deallocate(cap_state%export_surface)
  if (allocated(cap_state%export_mask)) deallocate(cap_state%export_mask)
  if (allocated(cap_state%export_thickness_elements)) deallocate(cap_state%export_thickness_elements)
  if (allocated(cap_state%export_surface_elements)) deallocate(cap_state%export_surface_elements)
  if (allocated(cap_state%export_mask_elements)) deallocate(cap_state%export_mask_elements)

  cap_state%handle = c_null_ptr
  cap_state%mpi_comm = 0_c_int
  cap_state%node_count = 0_c_int
  cap_state%element_count = 0_c_int
  cap_state%scalar_field_count = 0
  cap_state%scalar_idx_grid_nx = 0
  cap_state%scalar_idx_grid_ny = 0
  cap_state%scalar_idx_next_sw_cday = 0
  cap_state%scalar_idx_precip_factor = 0
  cap_state%mesh_created = .false.
  cap_state%write_restart = .false.
  cap_state%require_melt_import = .true.
  cap_state%melt_diagnostics = .false.
  cap_state%advance_count = 0
  cap_state%melt_diagnostics_interval = 1
  cap_state%issm_advance_seconds = 0
  cap_state%scalar_field_name = default_scalar_field_name
end subroutine ResetCapState

subroutine LogInfo(message)
  character(len=*), intent(in) :: message

  integer :: log_rc

  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO, rc=log_rc)
  write(*,'(a)') trim(message)
end subroutine LogInfo

subroutine LogError(message)
  character(len=*), intent(in) :: message

  integer :: log_rc

  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR, rc=log_rc)
  write(*,'(a)') trim(message)
end subroutine LogError

subroutine ValidateExports(rc)
  integer, intent(out) :: rc

  rc = ESMF_SUCCESS

  if (.not. allocated(cap_state%export_thickness) .or. &
      .not. allocated(cap_state%export_surface) .or. &
      .not. allocated(cap_state%export_mask)) then
    rc = 1
    call LogError('ISSM_NUOPC: export arrays are not allocated')
    return
  end if

  if (.not. all(ieee_is_finite(cap_state%export_thickness))) then
    rc = 1
    call LogError('ISSM_NUOPC: non-finite IceSheetThickness export')
    return
  end if
  if (.not. all(ieee_is_finite(cap_state%export_surface))) then
    rc = 1
    call LogError('ISSM_NUOPC: non-finite IceSheetSurfaceElevation export')
    return
  end if
  if (.not. all(ieee_is_finite(cap_state%export_mask))) then
    rc = 1
    call LogError('ISSM_NUOPC: non-finite IceSheetMask export')
    return
  end if
end subroutine ValidateExports

subroutine LogExportSummary(context, rc)
  character(len=*), intent(in) :: context
  integer, intent(out) :: rc

  character(len=512) :: message

  call ValidateExports(rc)
  if (rc /= ESMF_SUCCESS) return

  write(message,'(a,a,a,es14.6,a,es14.6,a,es14.6,a,es14.6,a,es14.6,a,es14.6,a)') &
    'ISSM_NUOPC: exports ', trim(context), &
    ' thickness=[', minval(cap_state%export_thickness), ',', maxval(cap_state%export_thickness), &
    '] surface=[', minval(cap_state%export_surface), ',', maxval(cap_state%export_surface), &
    '] mask=[', minval(cap_state%export_mask), ',', maxval(cap_state%export_mask), ']'
  call LogInfo(trim(message))
end subroutine LogExportSummary

logical function ShouldLogMeltDiagnostics()

  ShouldLogMeltDiagnostics = cap_state%melt_diagnostics .and. &
    mod(cap_state%advance_count, cap_state%melt_diagnostics_interval) == 0
end function ShouldLogMeltDiagnostics

subroutine LogMeltSummary(label, values, rc)
  character(len=*), intent(in) :: label
  real(c_double), intent(in) :: values(:)
  integer, intent(out) :: rc

  character(len=512) :: message

  rc = ESMF_SUCCESS
  if (.not. all(ieee_is_finite(values))) then
    rc = 1
    call LogError('ISSM_NUOPC: non-finite '//trim(label))
    return
  end if

  write(message,'(a,a,a,es14.6,a,es14.6,a,i0)') &
    'ISSM_NUOPC: ', trim(label), ' [m s-1] min=', minval(values), &
    ' max=', maxval(values), ' nonzero=', count(values /= 0.0_c_double)
  call LogInfo(trim(message))
end subroutine LogMeltSummary

subroutine VerifyRetainedMelt(rc)
  integer, intent(out) :: rc

  character(len=512) :: message
  real(c_double) :: max_difference
  real(c_double) :: scale
  real(c_double) :: tolerance

  rc = ESMF_SUCCESS
  if (.not. all(ieee_is_finite(cap_state%retained_melt))) then
    rc = 1
    call LogError('ISSM_NUOPC: non-finite retained floating melt after advance')
    return
  end if
  call ValidateMeltValues('retained floating melt after advance', cap_state%retained_melt, rc)
  if (rc /= ESMF_SUCCESS) return

  max_difference = maxval(abs(cap_state%retained_melt - cap_state%import_melt))
  scale = max(1.0e-12_c_double, maxval(abs(cap_state%import_melt)))
  tolerance = 1.0e-14_c_double + 1.0e-10_c_double * scale
  if (max_difference > tolerance) then
    rc = 1
    write(message,'(a,es14.6,a,es14.6)') &
      'ISSM_NUOPC: imported floating melt was overwritten; max_difference=', &
      max_difference, ' tolerance=', tolerance
    call LogError(trim(message))
  end if
end subroutine VerifyRetainedMelt

subroutine ValidateMeltValues(label, values, rc)
  character(len=*), intent(in) :: label
  real(c_double), intent(in) :: values(:)
  integer, intent(out) :: rc

  character(len=512) :: message

  rc = ESMF_SUCCESS
  if (.not. all(ieee_is_finite(values))) then
    rc = 1
    call LogError('ISSM_NUOPC: non-finite '//trim(label))
    return
  end if
  if (any(values < 0.0_c_double) .or. any(values > max_valid_melt_rate)) then
    rc = 1
    write(message,'(a,a,a,es14.6,a,es14.6)') &
      'ISSM_NUOPC: implausible ', trim(label), ' [m s-1] min=', minval(values), &
      ' max=', maxval(values)
    call LogError(trim(message))
  end if
end subroutine ValidateMeltValues

subroutine RealizeField(state, field_name, mesh, rc)
  type(ESMF_State), intent(inout) :: state
  character(len=*), intent(in) :: field_name
  type(ESMF_Mesh), intent(in) :: mesh
  integer, intent(out) :: rc

  type(ESMF_Field) :: field
  real(ESMF_KIND_R8), pointer :: field_ptr(:)

  rc = ESMF_SUCCESS

  field = ESMF_FieldCreate(mesh=mesh, typekind=ESMF_TYPEKIND_R8, &
    meshloc=ESMF_MESHLOC_ELEMENT, name=trim(field_name), rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call ESMF_FieldGet(field, farrayPtr=field_ptr, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  field_ptr(:) = 0.0
  call NUOPC_Realize(state, field=field, rc=rc)
end subroutine RealizeField

subroutine RealizeScalarField(state, rc)
  type(ESMF_State), intent(inout) :: state
  integer, intent(out) :: rc

  type(ESMF_DistGrid) :: distgrid
  type(ESMF_Field) :: field
  type(ESMF_Grid) :: grid

  rc = ESMF_SUCCESS
  if (cap_state%scalar_field_count <= 0) then
    rc = 1
    call LogError('ISSM_NUOPC: ScalarFieldCount must be positive')
    return
  end if

  distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
  if (rc /= ESMF_SUCCESS) return
  grid = ESMF_GridCreate(distgrid, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  field = ESMF_FieldCreate(name=trim(cap_state%scalar_field_name), grid=grid, &
    typekind=ESMF_TYPEKIND_R8, ungriddedLBound=(/1/), &
    ungriddedUBound=(/cap_state%scalar_field_count/), gridToFieldMap=(/2/), rc=rc)
  if (rc /= ESMF_SUCCESS) return

  call NUOPC_Realize(state, field=field, rc=rc)
end subroutine RealizeScalarField

subroutine RefreshExports(exportState, rc)
  type(ESMF_State), intent(inout) :: exportState
  integer, intent(out) :: rc

  rc = ESMF_SUCCESS
  if (.not. c_associated(cap_state%handle)) return

  call ISSM_NUOPC_ExportThickness(cap_state%handle, cap_state%export_thickness, cap_state%node_count)
  call ISSM_NUOPC_ExportSurface(cap_state%handle, cap_state%export_surface, cap_state%node_count)
  call ISSM_NUOPC_ExportMask(cap_state%handle, cap_state%export_mask, cap_state%node_count)
  call ValidateExports(rc)
  if (rc /= ESMF_SUCCESS) return

  call AverageNodeToElement(cap_state%export_thickness, cap_state%export_thickness_elements, rc)
  if (rc /= ESMF_SUCCESS) return
  call AverageNodeToElement(cap_state%export_surface, cap_state%export_surface_elements, rc)
  if (rc /= ESMF_SUCCESS) return
  call AverageNodeToElement(cap_state%export_mask, cap_state%export_mask_elements, rc)
  if (rc /= ESMF_SUCCESS) return

  if (cap_state%scalar_field_count > 0) then
    call SetExportScalars(exportState, rc)
    if (rc /= ESMF_SUCCESS) return
  end if

  call SetFieldData(exportState, export_thickness_name, cap_state%export_thickness_elements, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetFieldUpdated(exportState, export_thickness_name, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetFieldData(exportState, export_surface_name, cap_state%export_surface_elements, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetFieldUpdated(exportState, export_surface_name, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetFieldData(exportState, export_mask_name, cap_state%export_mask_elements, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetFieldUpdated(exportState, export_mask_name, rc)
  if (rc /= ESMF_SUCCESS) return
  if (cap_state%scalar_field_count > 0) then
    call SetFieldConstant(exportState, export_sg_area_name, 1.0_ESMF_KIND_R8, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldUpdated(exportState, export_sg_area_name, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldData(exportState, export_sg_icemask_name, cap_state%export_mask_elements, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldUpdated(exportState, export_sg_icemask_name, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldData(exportState, export_sg_icemask_fluxes_name, cap_state%export_mask_elements, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldUpdated(exportState, export_sg_icemask_fluxes_name, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldData(exportState, export_sg_covered_name, cap_state%export_mask_elements, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldUpdated(exportState, export_sg_covered_name, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldData(exportState, export_sg_topo_name, cap_state%export_surface_elements, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldUpdated(exportState, export_sg_topo_name, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldConstant(exportState, export_flgg_hflx_name, 0.0_ESMF_KIND_R8, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldUpdated(exportState, export_flgg_hflx_name, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldConstant(exportState, export_fgrg_rofl_name, 0.0_ESMF_KIND_R8, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldUpdated(exportState, export_fgrg_rofl_name, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldConstant(exportState, export_fgrg_rofi_name, 0.0_ESMF_KIND_R8, rc)
    if (rc /= ESMF_SUCCESS) return
    call SetFieldUpdated(exportState, export_fgrg_rofi_name, rc)
  end if
end subroutine RefreshExports

subroutine SetExportScalars(state, rc)
  type(ESMF_State), intent(inout) :: state
  integer, intent(out) :: rc

  rc = ESMF_SUCCESS

  call SetScalarValue(state, cap_state%scalar_idx_grid_nx, real(cap_state%element_count, ESMF_KIND_R8), rc)
  if (rc /= ESMF_SUCCESS) return
  call SetScalarValue(state, cap_state%scalar_idx_grid_ny, 1.0_ESMF_KIND_R8, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetScalarValue(state, cap_state%scalar_idx_next_sw_cday, 0.0_ESMF_KIND_R8, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetScalarValue(state, cap_state%scalar_idx_precip_factor, 0.0_ESMF_KIND_R8, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetFieldUpdated(state, cap_state%scalar_field_name, rc)
end subroutine SetExportScalars

subroutine SetScalarValue(state, scalar_id, value, rc)
  type(ESMF_State), intent(inout) :: state
  integer, intent(in) :: scalar_id
  real(ESMF_KIND_R8), intent(in) :: value
  integer, intent(out) :: rc

  type(ESMF_Field) :: field
  type(ESMF_VM) :: vm
  integer :: local_pet
  real(ESMF_KIND_R8), pointer :: scalar_ptr(:,:)

  rc = ESMF_SUCCESS
  if (scalar_id <= 0) return
  if (scalar_id > cap_state%scalar_field_count) then
    rc = 1
    call LogError('ISSM_NUOPC: scalar index exceeds ScalarFieldCount')
    return
  end if

  call ESMF_VMGetCurrent(vm, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call ESMF_VMGet(vm, localPet=local_pet, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call ESMF_StateGet(state, itemName=trim(cap_state%scalar_field_name), field=field, rc=rc)
  if (rc /= ESMF_SUCCESS) return

  if (local_pet == 0) then
    call ESMF_FieldGet(field, farrayPtr=scalar_ptr, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    scalar_ptr(scalar_id, 1) = value
  end if
end subroutine SetScalarValue

subroutine BuildElementCoordinates(rc)
  integer, intent(out) :: rc

  integer :: coord
  integer :: elem
  integer :: node
  integer :: vertex

  rc = ESMF_SUCCESS
  if (.not. allocated(cap_state%element_coords)) then
    rc = 1
    call LogError('ISSM_NUOPC: element coordinate array is not allocated')
    return
  end if

  cap_state%element_coords(:) = 0.0_c_double
  do elem = 1, int(cap_state%element_count)
    do vertex = 1, nodes_per_element
      node = cap_state%element_conn(nodes_per_element * (elem - 1) + vertex)
      if (node < 1 .or. node > int(cap_state%node_count)) then
        rc = 1
        call LogError('ISSM_NUOPC: element connectivity references an invalid node')
        return
      end if
      do coord = 1, coord_dim
        cap_state%element_coords(coord_dim * (elem - 1) + coord) = &
          cap_state%element_coords(coord_dim * (elem - 1) + coord) + &
          cap_state%node_coords(coord_dim * (node - 1) + coord)
      end do
    end do
    do coord = 1, coord_dim
      cap_state%element_coords(coord_dim * (elem - 1) + coord) = &
        cap_state%element_coords(coord_dim * (elem - 1) + coord) / real(nodes_per_element, c_double)
    end do
  end do
end subroutine BuildElementCoordinates

subroutine AverageNodeToElement(node_values, element_values, rc)
  real(c_double), intent(in) :: node_values(:)
  real(c_double), intent(out) :: element_values(:)
  integer, intent(out) :: rc

  integer :: elem
  integer :: node
  integer :: vertex

  rc = ESMF_SUCCESS
  if (size(node_values) /= int(cap_state%node_count) .or. &
      size(element_values) /= int(cap_state%element_count)) then
    rc = 1
    call LogError('ISSM_NUOPC: node-to-element averaging received unexpected array size')
    return
  end if

  do elem = 1, int(cap_state%element_count)
    element_values(elem) = 0.0_c_double
    do vertex = 1, nodes_per_element
      node = cap_state%element_conn(nodes_per_element * (elem - 1) + vertex)
      if (node < 1 .or. node > int(cap_state%node_count)) then
        rc = 1
        call LogError('ISSM_NUOPC: element connectivity references an invalid node')
        return
      end if
      element_values(elem) = element_values(elem) + node_values(node)
    end do
    element_values(elem) = element_values(elem) / real(nodes_per_element, c_double)
  end do
end subroutine AverageNodeToElement

subroutine AverageElementToNode(element_values, node_values, rc)
  real(c_double), intent(in) :: element_values(:)
  real(c_double), intent(out) :: node_values(:)
  integer, intent(out) :: rc

  integer, allocatable :: node_counts(:)
  integer :: elem
  integer :: node
  integer :: vertex

  rc = ESMF_SUCCESS
  if (size(element_values) /= int(cap_state%element_count) .or. &
      size(node_values) /= int(cap_state%node_count)) then
    rc = 1
    call LogError('ISSM_NUOPC: element-to-node averaging received unexpected array size')
    return
  end if

  allocate(node_counts(cap_state%node_count))
  node_values(:) = 0.0_c_double
  node_counts(:) = 0

  do elem = 1, int(cap_state%element_count)
    do vertex = 1, nodes_per_element
      node = cap_state%element_conn(nodes_per_element * (elem - 1) + vertex)
      if (node < 1 .or. node > int(cap_state%node_count)) then
        rc = 1
        call LogError('ISSM_NUOPC: element connectivity references an invalid node')
        deallocate(node_counts)
        return
      end if
      node_values(node) = node_values(node) + element_values(elem)
      node_counts(node) = node_counts(node) + 1
    end do
  end do

  do node = 1, int(cap_state%node_count)
    if (node_counts(node) > 0) then
      node_values(node) = node_values(node) / real(node_counts(node), c_double)
    end if
  end do
  deallocate(node_counts)

  if (.not. all(ieee_is_finite(node_values))) then
    rc = 1
    call LogError('ISSM_NUOPC: non-finite IceSheetBasalMeltRate import')
  end if
end subroutine AverageElementToNode

subroutine SetFieldData(state, field_name, values, rc)
  type(ESMF_State), intent(inout) :: state
  character(len=*), intent(in) :: field_name
  real(c_double), intent(in) :: values(:)
  integer, intent(out) :: rc

  type(ESMF_Field) :: field
  real(ESMF_KIND_R8), pointer :: field_ptr(:)

  rc = ESMF_SUCCESS

  call ESMF_StateGet(state, itemName=trim(field_name), field=field, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call ESMF_FieldGet(field, farrayPtr=field_ptr, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  field_ptr(:) = values(:)
end subroutine SetFieldData

subroutine SetFieldConstant(state, field_name, value, rc)
  type(ESMF_State), intent(inout) :: state
  character(len=*), intent(in) :: field_name
  real(ESMF_KIND_R8), intent(in) :: value
  integer, intent(out) :: rc

  type(ESMF_Field) :: field
  real(ESMF_KIND_R8), pointer :: field_ptr(:)

  rc = ESMF_SUCCESS
  call ESMF_StateGet(state, itemName=trim(field_name), field=field, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call ESMF_FieldGet(field, farrayPtr=field_ptr, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  field_ptr(:) = value
end subroutine SetFieldConstant

subroutine SetFieldUpdated(state, field_name, rc)
  type(ESMF_State), intent(inout) :: state
  character(len=*), intent(in) :: field_name
  integer, intent(out) :: rc

  type(ESMF_Field) :: field

  rc = ESMF_SUCCESS
  call ESMF_StateGet(state, itemName=trim(field_name), field=field, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call NUOPC_SetAttribute(field, name='Updated', value='true', rc=rc)
end subroutine SetFieldUpdated

logical function IsTrue(value)
  character(len=*), intent(in) :: value

  character(len=len(value)) :: lowered
  integer :: i

  lowered = adjustl(value)
  do i = 1, len_trim(lowered)
    select case (lowered(i:i))
    case ('A':'Z')
      lowered(i:i) = achar(iachar(lowered(i:i)) + 32)
    end select
  end do

  IsTrue = trim(lowered) == 'true' .or. trim(lowered) == '1' .or. trim(lowered) == 'yes'
end function IsTrue

end module ISSM_NUOPC_CapMod
