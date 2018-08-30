!> This defines some constants and procedures that are useful to the IO server and clients that call it. By using the contents
!! then the client can guarantee consistency against what the server expects
module io_server_client_mod
  use datadefn_mod, only : DEFAULT_PRECISION, SINGLE_PRECISION, DOUBLE_PRECISION, STRING_LENGTH
  use collections_mod, only : hashmap_type, iterator_type, mapentry_type, c_get_iterator, c_has_next, &
       c_next_mapentry, c_get_generic
  use conversions_mod, only : conv_to_string
  use mpi, only : MPI_CHARACTER, MPI_INT, MPI_LOGICAL, MPI_REAL, MPI_DOUBLE_PRECISION, MPI_ADDRESS_KIND, &
    MPI_BYTE, MPI_COMM_WORLD
  use edat
  use  ISO_C_BINDING, only : C_NULL_CHAR
  implicit none

#ifndef TEST_MODE
  private
#endif

  !< Data structure used to hold a sizing description of a field
  type data_sizing_description_type
     character(len=STRING_LENGTH) :: field_name !< Name of the field that this describes
     integer :: dimensions, dim_sizes(4) !< The number of dimensions and size in each dimension
  end type data_sizing_description_type

  type field_description_type
     character(len=STRING_LENGTH) :: definition_name, field_name
     integer :: field_type, data_type
     logical :: optional
  end type field_description_type

  type definition_description_type
     character(len=STRING_LENGTH) :: definition_name
     logical :: send_on_terminate
     integer :: number_fields, frequency
  end type definition_description_type

  ! Constants used in sending and receiving IO data
  integer, parameter :: COMMAND_TAG=9, DATA_TAG=10, REGISTER_COMMAND=1, DEREGISTER_COMMAND=3, DATA_COMMAND_START=4, &
       INTER_IO_COMMUNICATION=2

  !< Field type identifiers
  integer, parameter :: SCALAR_FIELD_TYPE = 1, ARRAY_FIELD_TYPE=2, MAP_FIELD_TYPE=3
  !< Field data type identifiers
  integer, parameter :: INTEGER_DATA_TYPE = 1, BOOLEAN_DATA_TYPE=2, STRING_DATA_TYPE=3, FLOAT_DATA_TYPE=4, &
       DOUBLE_DATA_TYPE=5

  character(len=STRING_LENGTH), parameter :: LOCAL_SIZES_KEY="local_sizes", LOCAL_START_POINTS_KEY="local_start_points", &
       LOCAL_END_POINTS_KEY="local_end_points", NUMBER_Q_INDICIES_KEY="num_q_indicies"

  public COMMAND_TAG, DATA_TAG, REGISTER_COMMAND, DEREGISTER_COMMAND, DATA_COMMAND_START, INTER_IO_COMMUNICATION, &
       SCALAR_FIELD_TYPE, ARRAY_FIELD_TYPE, MAP_FIELD_TYPE, INTEGER_DATA_TYPE, BOOLEAN_DATA_TYPE, &
       STRING_DATA_TYPE, FLOAT_DATA_TYPE, LOCAL_SIZES_KEY, LOCAL_START_POINTS_KEY, LOCAL_END_POINTS_KEY, NUMBER_Q_INDICIES_KEY, &
       DOUBLE_DATA_TYPE, data_sizing_description_type, populate_type_extents, definition_description_type, &
       field_description_type, build_mpi_type_field_description, build_mpi_type_definition_description, &
       pack_scalar_field, pack_array_field, pack_map_field, get_data_description_from_name, send_data_as_edat_event

contains

  subroutine send_data_as_edat_event(data_type, num_elements, target_rank, event_id, int_array_to_send, &
      char_array_to_send, int_to_send, data_sizing_to_send)
    integer, dimension(:), intent(in), optional :: int_array_to_send
    integer, intent(in), optional :: int_to_send
    character, dimension(:), intent(in), optional :: char_array_to_send
    type(data_sizing_description_type), intent(in), optional :: data_sizing_to_send(:)
    integer, intent(in) :: data_type, num_elements, target_rank
    character(len=*), intent(in) :: event_id

    character, dimension(:), allocatable :: buffer
    integer :: element_length, current_point, eid_len, i, ierr, my_rank

    if (data_type == EDAT_NOTYPE .or. data_type == EDAT_NONE) then
      element_length=0
    else if (data_type == EDAT_INT .or. data_type == EDAT_FLOAT) then
      element_length=4
    else if (data_type == EDAT_DOUBLE) then
      element_length=8
    else if (data_type == EDAT_BYTE) then
      element_length=1
    end if

    eid_len=len(trim(event_id))

    allocate(buffer(eid_len + 2 + (element_length * num_elements) + 12))

    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

    current_point=pack_scalar_field(buffer, 1, int_value=data_type)
    current_point=pack_scalar_field(buffer, current_point, int_value=my_rank)
    current_point=pack_scalar_field(buffer, current_point, int_value=eid_len)
    current_point=pack_scalar_field(buffer, current_point, char_value=char(0))
    do i=1, eid_len
      current_point=pack_scalar_field(buffer, current_point, char_value=event_id(i:i))
    end do
    current_point=pack_scalar_field(buffer, current_point, char_value=C_NULL_CHAR)
    if (present(int_array_to_send)) then
      current_point=pack_array_field(buffer, current_point, int_array=int_array_to_send)
    else if (present(char_array_to_send)) then
      current_point=pack_array_field(buffer, current_point, char_array=char_array_to_send)
    else if (present(int_to_send)) then
      current_point=pack_scalar_field(buffer, current_point, int_value=int_to_send)
    else if (present(data_sizing_to_send)) then
      buffer(current_point:(current_point+(element_length * num_elements))-1)=&
        transfer(data_sizing_to_send, buffer(current_point:(current_point+(element_length * num_elements))-1))
    end if

    call MPI_Send(buffer, size(buffer), MPI_BYTE, target_rank, 16384, MPI_COMM_WORLD, ierr)
  end subroutine send_data_as_edat_event

  !> Provides the type extents of the types that we are using in construction of the MPI data type. This is done
  !! once for fast look up in the actual construction phase.
  !! @returns An array of type extents keyed by the DATA_TYPE constants of the configuration parser module)
  function populate_type_extents()
    integer :: populate_type_extents(5)

    integer :: int_data
    real :: float_data
    real*8 :: double_data
    logical :: bool_data
    character :: char_data
    integer(kind=8) :: large_number_extents(5)

    large_number_extents(INTEGER_DATA_TYPE)=kind(int_data)
    large_number_extents(BOOLEAN_DATA_TYPE)=kind(bool_data)
    large_number_extents(STRING_DATA_TYPE)=kind(char_data)
    large_number_extents(FLOAT_DATA_TYPE)=kind(float_data)
    large_number_extents(DOUBLE_DATA_TYPE)=kind(double_data)

    populate_type_extents=int(large_number_extents)
  end function populate_type_extents

  !> Builds the MPI data type for sending data descriptions to registree MONCs
  integer function build_mpi_type_definition_description()
    integer :: new_type, ierr, block_counts(4), old_types(4), offsets(4)
    integer(MPI_ADDRESS_KIND) :: num_addr, base_addr

    type(definition_description_type) :: basic_type

    call mpi_get_address(basic_type, base_addr, ierr)
    old_types(1) = MPI_CHARACTER
    block_counts(1) = STRING_LENGTH
    offsets(1)=0

    call mpi_get_address(basic_type%send_on_terminate, num_addr, ierr)
    old_types(2) = MPI_LOGICAL
    block_counts(2) = 1
    offsets(2)=int(num_addr-base_addr)

    call mpi_get_address(basic_type%number_fields, num_addr, ierr)
    old_types(3) = MPI_INT
    block_counts(3) = 1
    offsets(3)=int(num_addr-base_addr)

    call mpi_get_address(basic_type%frequency, num_addr, ierr)
    old_types(4) = MPI_INT
    block_counts(4) = 1
    offsets(4)=int(num_addr-base_addr)

    call mpi_type_struct(4, block_counts, offsets, old_types, new_type, ierr)
    call mpi_type_commit(new_type, ierr)
    build_mpi_type_definition_description=new_type
  end function build_mpi_type_definition_description

  !> Builds the MPI data type for sending field descriptions to registree MONCs
  integer function build_mpi_type_field_description()
    integer :: new_type, ierr, old_types(5), block_counts(5), offsets(5)
    integer(MPI_ADDRESS_KIND) :: num_addr, base_addr

    type(field_description_type) :: basic_type

    call mpi_get_address(basic_type, base_addr, ierr)
    old_types(1) = MPI_CHARACTER
    block_counts(1) = STRING_LENGTH
    offsets(1)=0

    call mpi_get_address(basic_type%field_name, num_addr, ierr)
    old_types(2) = MPI_CHARACTER
    block_counts(2) = STRING_LENGTH
    offsets(2)=int(num_addr-base_addr)

    call mpi_get_address(basic_type%field_type, num_addr, ierr)
    old_types(3) = MPI_INT
    block_counts(3) = 1
    offsets(3)=int(num_addr-base_addr)

    call mpi_get_address(basic_type%data_type, num_addr, ierr)
    old_types(4) = MPI_INT
    block_counts(4) = 1
    offsets(4)=int(num_addr-base_addr)

    call mpi_get_address(basic_type%optional, num_addr, ierr)
    old_types(5) = MPI_LOGICAL
    block_counts(5) = 1
    offsets(5)=int(num_addr-base_addr)

    call mpi_type_struct(5, block_counts, offsets, old_types, new_type, ierr)
    call mpi_type_commit(new_type, ierr)
    build_mpi_type_field_description=new_type
  end function build_mpi_type_field_description

    !> Packs a map into the send buffer
  !! @param buffer The buffer to pack the field into
  !! @param start_offset The starting offset to write into the buffer
  !! @param map_to_pack The map to pack into the buffer
  !! @returns The next location in the buffer to write to (next start offset)
  integer function pack_map_field(buffer, start_offset, map_to_pack)
    character, dimension(:), intent(inout) :: buffer
    integer, intent(in) :: start_offset
    type(hashmap_type) :: map_to_pack

    integer :: i, target_end, current_offset
    character(len=STRING_LENGTH) :: temp_string
    character(len=STRING_LENGTH), pointer :: sized_raw_character
    class(*), pointer :: raw_data, raw_to_string
    type(iterator_type) :: map_iterator
    type(mapentry_type) :: specific_mapentry

    current_offset=start_offset
    map_iterator=c_get_iterator(map_to_pack)
    do while (c_has_next(map_iterator))
      specific_mapentry=c_next_mapentry(map_iterator)
      temp_string=specific_mapentry%key
      target_end=current_offset+STRING_LENGTH-1
      buffer(current_offset:target_end)=transfer(temp_string, buffer(current_offset:target_end))
      current_offset=target_end+1

      raw_data=>c_get_generic(specific_mapentry)
      raw_to_string=>raw_data
      select type (raw_data)
      type is(integer)
        temp_string=conv_to_string(raw_data)
      type is(real(4))
        temp_string=conv_to_string(raw_data)
      type is (real(8))
        temp_string=conv_to_string(raw_data)
      type is(logical)
        temp_string=conv_to_string(raw_data)
      type is(character(len=*))
        sized_raw_character=>conv_to_string(raw_to_string, .false., STRING_LENGTH)
        temp_string=sized_raw_character
      end select
      target_end=current_offset+STRING_LENGTH-1
      buffer(current_offset:target_end)=transfer(temp_string, buffer(current_offset:target_end))
      current_offset=target_end+1
    end do
    pack_map_field=current_offset
  end function pack_map_field

  !> Packs an array field into the sending buffer
  !! @param buffer The buffer to pack the field into
  !! @param start_offset The starting offset to write into the buffer
  !! @param int_value (Optional) integer array values to pack
  !! @param real_value (Optional) default precision real array values to pack
  !! @returns The next location in the buffer to write to (next start offset)
  integer function pack_array_field(buffer, start_offset, int_array, char_array, real_array_1d, real_array_2d, &
      real_array_3d, real_array_4d)
    character, dimension(:), intent(inout) :: buffer
    integer, intent(in) :: start_offset
    integer, dimension(:), intent(in), optional :: int_array
    character, dimension(:), intent(in), optional :: char_array
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in), optional :: real_array_1d
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in), optional :: real_array_2d
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in), optional :: real_array_3d
    real(kind=DEFAULT_PRECISION), dimension(:,:,:,:), intent(in), optional :: real_array_4d

    integer :: target_end

    if (present(int_array)) then
      target_end=start_offset+kind(int_array)*size(int_array)-1
      buffer(start_offset:target_end) = transfer(int_array, buffer(start_offset:target_end))
    else if (present(char_array)) then
      target_end=start_offset+kind(char_array)*size(char_array)-1
      buffer(start_offset:target_end) = transfer(char_array, buffer(start_offset:target_end))
    else if (present(real_array_1d)) then
      target_end=start_offset+kind(real_array_1d)*size(real_array_1d)-1
      buffer(start_offset:target_end) = transfer(real_array_1d, buffer(start_offset:target_end))
    else if (present(real_array_2d)) then
      target_end=start_offset+kind(real_array_2d)*size(real_array_2d)-1
      buffer(start_offset:target_end) = transfer(real_array_2d, buffer(start_offset:target_end))
    else if (present(real_array_3d)) then
      target_end=start_offset+kind(real_array_3d)*size(real_array_3d)-1
      buffer(start_offset:target_end) = transfer(real_array_3d, buffer(start_offset:target_end))
    else if (present(real_array_4d)) then
      target_end=start_offset+kind(real_array_4d)*size(real_array_4d)-1
      buffer(start_offset:target_end) = transfer(real_array_4d, buffer(start_offset:target_end))
    end if
    pack_array_field=target_end+1
  end function pack_array_field

  !> Packs the data of a scalar field into a buffer
  !! @param buffer The buffer to pack the field into
  !! @param start_offset The starting offset to write into the buffer
  !! @param int_value (Optional) integer scalar value to pack
  !! @param real_value (Optional) default precision real scalar value to pack
  !! @param single_real_value (Optional) single precision real scalar value to pack
  !! @param double_real_value (Optional) double precision real scalar value to pack
  !! @returns The next location in the buffer to write to (next start offset)
  integer function pack_scalar_field(buffer, start_offset, int_value, real_value, single_real_value, double_real_value, &
       string_value, logical_value, char_value)
    character, dimension(:), intent(inout) :: buffer
    integer, intent(in) :: start_offset
    integer, intent(in), optional :: int_value
    real(kind=DEFAULT_PRECISION), intent(in), optional :: real_value
    real(kind=SINGLE_PRECISION), intent(in), optional :: single_real_value
    real(kind=DOUBLE_PRECISION), intent(in), optional :: double_real_value
    character(len=*), intent(in), optional :: string_value
    logical, intent(in), optional :: logical_value
    character, intent(in), optional :: char_value

    integer :: target_end
    character(len=STRING_LENGTH) :: string_to_insert

    if (present(int_value)) then
      target_end=start_offset+kind(int_value)-1
      buffer(start_offset:target_end) = transfer(int_value, buffer(start_offset:target_end))
    else if (present(real_value)) then
      target_end=start_offset+kind(real_value)-1
      buffer(start_offset:target_end) = transfer(real_value, buffer(start_offset:target_end))
    else if (present(single_real_value)) then
      target_end=start_offset+kind(single_real_value)-1
      buffer(start_offset:target_end) = transfer(single_real_value, buffer(start_offset:target_end))
    else if (present(double_real_value)) then
      target_end=start_offset+kind(double_real_value)-1
      buffer(start_offset:target_end) = transfer(double_real_value, buffer(start_offset:target_end))
    else if (present(string_value)) then
      target_end=start_offset+STRING_LENGTH-1
      string_to_insert=string_value
      buffer(start_offset:target_end) = transfer(string_to_insert, buffer(start_offset:target_end))
    else if (present(logical_value)) then
      target_end=start_offset+kind(logical_value)-1
      buffer(start_offset:target_end) = transfer(logical_value, buffer(start_offset:target_end))
    else if (present(char_value)) then
      target_end=start_offset+kind(char_value)-1
      buffer(start_offset:target_end) = transfer(char_value, buffer(start_offset:target_end))
    else
      target_end=start_offset-1
    end if
    pack_scalar_field=target_end+1
  end function pack_scalar_field

  !> Look up the data description that corresponds to a specific field keyed by its name
  !! @param descriptions The data descriptions
  !! @param name The field name to look up
  !! @param field_description The resulting field description if the field is found
  !! @returns Whether or not the field, based upon its name lookup, was found
  logical function get_data_description_from_name(descriptions, name, field_description)
    type(data_sizing_description_type), dimension(:), intent(in) :: descriptions
    type(data_sizing_description_type), intent(out), optional :: field_description
    character(len=*), intent(in) :: name

    integer :: i
    do i=1,size(descriptions)
      if (descriptions(i)%field_name == name) then
        if (present(field_description)) field_description=descriptions(i)
        get_data_description_from_name=.true.
        return
      end if
    end do
    get_data_description_from_name=.false.
  end function get_data_description_from_name
end module io_server_client_mod
