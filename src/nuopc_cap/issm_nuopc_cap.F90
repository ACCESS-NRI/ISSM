module ISSM_NUOPC_CapMod

use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, c_int, c_null_char, c_null_ptr, c_ptr

use ESMF, only: ESMF_Clock, ESMF_ClockGet, ESMF_COORDSYS_CART
use ESMF, only: ESMF_Field, ESMF_FieldCreate, ESMF_FieldGet
use ESMF, only: ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSetEntryPoint
use ESMF, only: ESMF_LOGMSG_ERROR, ESMF_LOGMSG_INFO, ESMF_LogWrite
use ESMF, only: ESMF_Mesh, ESMF_MeshCreate
use ESMF, only: ESMF_MESHLOC_NODE, ESMF_METHOD_INITIALIZE
use ESMF, only: ESMF_State, ESMF_StateGet
use ESMF, only: ESMF_SUCCESS, ESMF_TimeInterval, ESMF_TimeIntervalGet
use ESMF, only: ESMF_KIND_R8, ESMF_TYPEKIND_R8, ESMF_VM, ESMF_VMGet

use NUOPC, only: NUOPC_Advertise, NUOPC_CompAttributeSet, NUOPC_CompDerive
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
character(len=*), parameter :: import_melt_name = 'floatingIceMeltRate'
character(len=*), parameter :: import_melt_stdname = 'IceSheetBasalMeltRate'
character(len=*), parameter :: export_thickness_name = 'iceThickness'
character(len=*), parameter :: export_thickness_stdname = 'IceSheetThickness'
character(len=*), parameter :: export_surface_name = 'iceSurface'
character(len=*), parameter :: export_surface_stdname = 'IceSheetSurfaceElevation'
character(len=*), parameter :: export_mask_name = 'iceMask'
character(len=*), parameter :: export_mask_stdname = 'IceSheetMask'

type issm_cap_state_type
  type(c_ptr) :: handle = c_null_ptr
  type(ESMF_Mesh) :: mesh
  integer(c_int) :: mpi_comm = 0_c_int
  integer(c_int) :: node_count = 0_c_int
  integer(c_int) :: element_count = 0_c_int
  logical :: mesh_created = .false.
  logical :: write_restart = .false.
  integer(c_int), allocatable :: node_ids(:)
  integer(c_int), allocatable :: node_owners(:)
  integer(c_int), allocatable :: element_ids(:)
  integer(c_int), allocatable :: element_types(:)
  integer(c_int), allocatable :: element_conn(:)
  real(c_double), allocatable :: node_coords(:)
  real(c_double), allocatable :: import_melt(:)
  real(c_double), allocatable :: export_thickness(:)
  real(c_double), allocatable :: export_surface(:)
  real(c_double), allocatable :: export_mask(:)
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

  cap_state%handle = ISSM_NUOPC_CreateFromCase( &
    trim(case_dir)//c_null_char, &
    trim(model_name)//c_null_char, &
    trim(solution_name)//c_null_char, &
    cap_state%mpi_comm)
  if (.not. c_associated(cap_state%handle)) then
    call ESMF_LogWrite('ISSM_NUOPC: failed to create ISSM model handle', ESMF_LOGMSG_ERROR, rc=rc)
    return
  end if

  call ISSM_NUOPC_GetMeshCounts(cap_state%handle, cap_state%node_count, cap_state%element_count)
  call AllocateCapState(rc)
  if (rc /= ESMF_SUCCESS) return

  call NUOPC_Advertise( &
    importState, &
    standardName=import_melt_stdname, &
    name=import_melt_name, &
    TransferOfferGeomObject='will provide', &
    SharePolicyField='share', &
    SharePolicyGeomObject='share', &
    rc=rc &
  )
  if (rc /= ESMF_SUCCESS) return
  call NUOPC_Advertise(exportState, standardName=export_thickness_stdname, name=export_thickness_name, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call NUOPC_Advertise(exportState, standardName=export_surface_stdname, name=export_surface_name, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call NUOPC_Advertise(exportState, standardName=export_mask_stdname, name=export_mask_name, rc=rc)
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

  cap_state%mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, coordSys=ESMF_COORDSYS_CART, &
       nodeIds=cap_state%node_ids, nodeCoords=cap_state%node_coords, &
       nodeOwners=cap_state%node_owners, elementIds=cap_state%element_ids, &
       elementTypes=cap_state%element_types, elementConn=cap_state%element_conn, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  cap_state%mesh_created = .true.

  call RealizeField(importState, import_melt_name, cap_state%mesh, rc)
  if (rc /= ESMF_SUCCESS) return
  call RealizeField(exportState, export_thickness_name, cap_state%mesh, rc)
  if (rc /= ESMF_SUCCESS) return
  call RealizeField(exportState, export_surface_name, cap_state%mesh, rc)
  if (rc /= ESMF_SUCCESS) return
  call RealizeField(exportState, export_mask_name, cap_state%mesh, rc)
  if (rc /= ESMF_SUCCESS) return

  call RefreshExports(exportState, rc)
end subroutine InitializeRealize

subroutine DataInitialize(gcomp, rc)
  type(ESMF_GridComp) :: gcomp
  integer, intent(out) :: rc

  type(ESMF_State) :: exportState

  rc = ESMF_SUCCESS
  call ESMF_GridCompGet(gcomp, exportState=exportState, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call RefreshExports(exportState, rc)
  if (rc /= ESMF_SUCCESS) return
  call NUOPC_CompAttributeSet(gcomp, name='InitializeDataComplete', value='true', rc=rc)
end subroutine DataInitialize

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
  integer :: dt_seconds

  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, exportState=exportState, rc=rc)
  if (rc /= ESMF_SUCCESS) return

  call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call ESMF_TimeIntervalGet(timeStep, s=dt_seconds, rc=rc)
  if (rc /= ESMF_SUCCESS) return

  call ESMF_StateGet(importState, itemName=import_melt_name, field=field, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call ESMF_FieldGet(field, farrayPtr=field_ptr, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  cap_state%import_melt(:) = field_ptr(:)
  call ISSM_NUOPC_ImportFloatingMelt(cap_state%handle, cap_state%import_melt, cap_state%node_count)

  call ISSM_NUOPC_Advance(cap_state%handle, real(dt_seconds, c_double))
  call RefreshExports(exportState, rc)
end subroutine ModelAdvance

subroutine FinalizeModel(gcomp, rc)
  type(ESMF_GridComp) :: gcomp
  integer, intent(out) :: rc

  rc = ESMF_SUCCESS

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
  allocate(cap_state%import_melt(cap_state%node_count))
  allocate(cap_state%export_thickness(cap_state%node_count))
  allocate(cap_state%export_surface(cap_state%node_count))
  allocate(cap_state%export_mask(cap_state%node_count))

  cap_state%import_melt = 0.0_c_double
  cap_state%export_thickness = 0.0_c_double
  cap_state%export_surface = 0.0_c_double
  cap_state%export_mask = 0.0_c_double
end subroutine AllocateCapState

subroutine ResetCapState()
  if (allocated(cap_state%node_ids)) deallocate(cap_state%node_ids)
  if (allocated(cap_state%node_owners)) deallocate(cap_state%node_owners)
  if (allocated(cap_state%node_coords)) deallocate(cap_state%node_coords)
  if (allocated(cap_state%element_ids)) deallocate(cap_state%element_ids)
  if (allocated(cap_state%element_types)) deallocate(cap_state%element_types)
  if (allocated(cap_state%element_conn)) deallocate(cap_state%element_conn)
  if (allocated(cap_state%import_melt)) deallocate(cap_state%import_melt)
  if (allocated(cap_state%export_thickness)) deallocate(cap_state%export_thickness)
  if (allocated(cap_state%export_surface)) deallocate(cap_state%export_surface)
  if (allocated(cap_state%export_mask)) deallocate(cap_state%export_mask)

  cap_state%handle = c_null_ptr
  cap_state%mpi_comm = 0_c_int
  cap_state%node_count = 0_c_int
  cap_state%element_count = 0_c_int
  cap_state%mesh_created = .false.
  cap_state%write_restart = .false.
end subroutine ResetCapState

subroutine RealizeField(state, field_name, mesh, rc)
  type(ESMF_State), intent(inout) :: state
  character(len=*), intent(in) :: field_name
  type(ESMF_Mesh), intent(in) :: mesh
  integer, intent(out) :: rc

  type(ESMF_Field) :: field
  real(ESMF_KIND_R8), pointer :: field_ptr(:)

  rc = ESMF_SUCCESS

  field = ESMF_FieldCreate(mesh=mesh, typekind=ESMF_TYPEKIND_R8, &
    meshloc=ESMF_MESHLOC_NODE, name=trim(field_name), rc=rc)
  if (rc /= ESMF_SUCCESS) return
  call ESMF_FieldGet(field, farrayPtr=field_ptr, rc=rc)
  if (rc /= ESMF_SUCCESS) return
  field_ptr(:) = 0.0
  call NUOPC_Realize(state, field=field, rc=rc)
end subroutine RealizeField

subroutine RefreshExports(exportState, rc)
  type(ESMF_State), intent(inout) :: exportState
  integer, intent(out) :: rc

  rc = ESMF_SUCCESS
  if (.not. c_associated(cap_state%handle)) return

  call ISSM_NUOPC_ExportThickness(cap_state%handle, cap_state%export_thickness, cap_state%node_count)
  call ISSM_NUOPC_ExportSurface(cap_state%handle, cap_state%export_surface, cap_state%node_count)
  call ISSM_NUOPC_ExportMask(cap_state%handle, cap_state%export_mask, cap_state%node_count)

  call SetFieldData(exportState, export_thickness_name, cap_state%export_thickness, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetFieldUpdated(exportState, export_thickness_name, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetFieldData(exportState, export_surface_name, cap_state%export_surface, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetFieldUpdated(exportState, export_surface_name, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetFieldData(exportState, export_mask_name, cap_state%export_mask, rc)
  if (rc /= ESMF_SUCCESS) return
  call SetFieldUpdated(exportState, export_mask_name, rc)
end subroutine RefreshExports

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
