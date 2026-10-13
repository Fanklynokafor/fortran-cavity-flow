module cavity_model
  use, intrinsic :: iso_fortran_env, only: real64, error_unit
  implicit none
  private

  public :: simulation_settings, flow_field, read_settings
  public :: initialize_flow, advance_flow, write_results

  type :: simulation_settings
    integer :: nx = 41
    integer :: ny = 41
    integer :: steps = 4000
    integer :: poisson_iterations = 240
    real(real64) :: length_x = 1.0_real64
    real(real64) :: length_y = 1.0_real64
    real(real64) :: reynolds_number = 100.0_real64
    real(real64) :: lid_speed = 1.0_real64
    real(real64) :: time_step = 0.0025_real64
    real(real64) :: poisson_tolerance = 1.0e-7_real64
    real(real64) :: sor_factor = 1.65_real64
  end type simulation_settings

  type :: flow_field
    ! Boundary nodes are included in these arrays; there are no ghost cells.
    real(real64), allocatable :: streamfunction(:, :)
    real(real64), allocatable :: vorticity(:, :)
    real(real64), allocatable :: next_vorticity(:, :)
    real(real64), allocatable :: u(:, :)
    real(real64), allocatable :: v(:, :)
    integer :: poisson_nonconverged = 0
    real(real64) :: poisson_worst_change = 0.0_real64
  end type flow_field

contains

  subroutine read_settings(filename, settings)
    character(len=*), intent(in) :: filename
    type(simulation_settings), intent(inout) :: settings
    integer :: input_unit, io_status
    integer :: nx, ny, steps, poisson_iterations
    real(real64) :: length_x, length_y, reynolds_number, lid_speed
    real(real64) :: time_step, poisson_tolerance, sor_factor
    real(real64) :: dx, dy, viscosity, advective_limit, diffusive_limit
    namelist /simulation/ nx, ny, steps, poisson_iterations, length_x, length_y, &
      reynolds_number, lid_speed, time_step, poisson_tolerance, sor_factor

    ! Seed NAMELIST variables from the defaults, then override only supplied values.
    nx = settings%nx
    ny = settings%ny
    steps = settings%steps
    poisson_iterations = settings%poisson_iterations
    length_x = settings%length_x
    length_y = settings%length_y
    reynolds_number = settings%reynolds_number
    lid_speed = settings%lid_speed
    time_step = settings%time_step
    poisson_tolerance = settings%poisson_tolerance
    sor_factor = settings%sor_factor

    open(newunit=input_unit, file=filename, status='old', action='read', iostat=io_status)
    if (io_status /= 0) call fail('Cannot open simulation.nml.')
    read(input_unit, nml=simulation, iostat=io_status)
    close(input_unit)
    if (io_status /= 0) call fail('Could not read the &simulation namelist.')

    settings%nx = nx
    settings%ny = ny
    settings%steps = steps
    settings%poisson_iterations = poisson_iterations
    settings%length_x = length_x
    settings%length_y = length_y
    settings%reynolds_number = reynolds_number
    settings%lid_speed = lid_speed
    settings%time_step = time_step
    settings%poisson_tolerance = poisson_tolerance
    settings%sor_factor = sor_factor

    if (settings%nx < 5 .or. settings%ny < 5) call fail('Use at least 5 nodes in each direction.')
    if (settings%steps < 1 .or. settings%poisson_iterations < 1) call fail('Iteration counts must be positive.')
    if (settings%length_x <= 0.0_real64 .or. settings%length_y <= 0.0_real64) &
      call fail('Domain dimensions must be positive.')
    if (settings%reynolds_number <= 0.0_real64 .or. settings%lid_speed <= 0.0_real64) &
      call fail('Reynolds number and lid speed must be positive.')
    if (settings%time_step <= 0.0_real64) call fail('The time step must be positive.')
    if (settings%sor_factor <= 0.0_real64 .or. settings%sor_factor >= 2.0_real64) &
      call fail('SOR factor must be greater than zero and less than two.')
    if (settings%poisson_tolerance <= 0.0_real64) call fail('Poisson tolerance must be positive.')

    dx = settings%length_x / real(settings%nx - 1, real64)
    dy = settings%length_y / real(settings%ny - 1, real64)
    viscosity = settings%lid_speed * settings%length_y / settings%reynolds_number
    advective_limit = 0.35_real64 / &
      (settings%lid_speed * (1.0_real64 / dx + 1.0_real64 / dy))
    diffusive_limit = 0.20_real64 / &
      (2.0_real64 * viscosity * (1.0_real64 / dx**2 + 1.0_real64 / dy**2))
    if (settings%time_step > min(advective_limit, diffusive_limit)) then
      call fail('Time step exceeds the configured advection/diffusion stability limit.')
    end if
  end subroutine read_settings

  subroutine initialize_flow(settings, flow)
    type(simulation_settings), intent(in) :: settings
    type(flow_field), intent(out) :: flow

    allocate(flow%streamfunction(settings%nx, settings%ny))
    allocate(flow%vorticity(settings%nx, settings%ny))
    allocate(flow%next_vorticity(settings%nx, settings%ny))
    allocate(flow%u(settings%nx, settings%ny))
    allocate(flow%v(settings%nx, settings%ny))
    flow%streamfunction = 0.0_real64
    flow%vorticity = 0.0_real64
    flow%next_vorticity = 0.0_real64
    flow%u = 0.0_real64
    flow%v = 0.0_real64
    flow%poisson_nonconverged = 0
    flow%poisson_worst_change = 0.0_real64
  end subroutine initialize_flow

  subroutine advance_flow(settings, flow)
    type(simulation_settings), intent(in) :: settings
    type(flow_field), intent(inout) :: flow
    real(real64) :: dx, dy, viscosity, dwdx, dwdy, laplacian
    integer :: i, j

    dx = settings%length_x / real(settings%nx - 1, real64)
    dy = settings%length_y / real(settings%ny - 1, real64)
    viscosity = settings%lid_speed * settings%length_y / settings%reynolds_number

    ! With u = d(psi)/dy and v = -d(psi)/dx, continuity is satisfied
    ! automatically. This convention gives laplacian(psi) = omega here,
    ! where omega is du/dy - dv/dx (the negative of the usual z-vorticity).
    call solve_streamfunction(settings, flow)
    call apply_wall_vorticity(settings, flow)
    call calculate_velocity(settings, flow)

    ! Keep old and new vorticity separate during the update: every node should
    ! read the same time level, rather than a mix of old and already-updated values.
    flow%next_vorticity = flow%vorticity

    ! The transport equation is d(omega)/dt + u*d(omega)/dx + v*d(omega)/dy
    ! = nu*laplacian(omega). Donor-cell advection is more diffusive than a central stencil,
    ! but it behaves better near the sharp moving-lid corners in this teaching model.
    do j = 2, settings%ny - 1
      do i = 2, settings%nx - 1
        if (flow%u(i, j) >= 0.0_real64) then
          dwdx = (flow%vorticity(i, j) - flow%vorticity(i - 1, j)) / dx
        else
          dwdx = (flow%vorticity(i + 1, j) - flow%vorticity(i, j)) / dx
        end if
        if (flow%v(i, j) >= 0.0_real64) then
          dwdy = (flow%vorticity(i, j) - flow%vorticity(i, j - 1)) / dy
        else
          dwdy = (flow%vorticity(i, j + 1) - flow%vorticity(i, j)) / dy
        end if
        laplacian = (flow%vorticity(i + 1, j) - 2.0_real64 * flow%vorticity(i, j) &
          + flow%vorticity(i - 1, j)) / dx**2 &
          + (flow%vorticity(i, j + 1) - 2.0_real64 * flow%vorticity(i, j) &
          + flow%vorticity(i, j - 1)) / dy**2
        flow%next_vorticity(i, j) = flow%vorticity(i, j) + settings%time_step &
          * (-flow%u(i, j) * dwdx - flow%v(i, j) * dwdy + viscosity * laplacian)
      end do
    end do

    flow%vorticity = flow%next_vorticity
  end subroutine advance_flow

  subroutine solve_streamfunction(settings, flow)
    type(simulation_settings), intent(in) :: settings
    type(flow_field), intent(inout) :: flow
    real(real64) :: dx, dy, denominator, candidate, largest_change
    integer :: i, j, iteration

    dx = settings%length_x / real(settings%nx - 1, real64)
    dy = settings%length_y / real(settings%ny - 1, real64)
    denominator = 2.0_real64 / dx**2 + 2.0_real64 / dy**2

    ! Rearranging laplacian(psi) = omega gives the weighted neighbour update.
    ! SOR reuses the previous time step's solution as its initial guess.
    ! That warm start is why a modest iteration cap is usually enough here.
    do iteration = 1, settings%poisson_iterations
      largest_change = 0.0_real64
      do j = 2, settings%ny - 1
        do i = 2, settings%nx - 1
          candidate = ((flow%streamfunction(i + 1, j) + flow%streamfunction(i - 1, j)) / dx**2 &
            + (flow%streamfunction(i, j + 1) + flow%streamfunction(i, j - 1)) / dy**2 &
            - flow%vorticity(i, j)) / denominator
          candidate = flow%streamfunction(i, j) + settings%sor_factor &
            * (candidate - flow%streamfunction(i, j))
          largest_change = max(largest_change, abs(candidate - flow%streamfunction(i, j)))
          flow%streamfunction(i, j) = candidate
        end do
      end do
      if (largest_change < settings%poisson_tolerance) exit
    end do
    ! TODO: This SOR sweep gets the job done, but it's not exactly speedy on
    ! big grids. Try red-black SOR or multigrid before cranking up the size.
    if (largest_change >= settings%poisson_tolerance) then
      flow%poisson_nonconverged = flow%poisson_nonconverged + 1
      flow%poisson_worst_change = max(flow%poisson_worst_change, largest_change)
    end if
    ! TODO: Update size is a handy quick check, but let's add a real Poisson
    ! residual too before trusting this on a much bigger grid. No guesswork.
  end subroutine solve_streamfunction

  subroutine apply_wall_vorticity(settings, flow)
    type(simulation_settings), intent(in) :: settings
    type(flow_field), intent(inout) :: flow
    real(real64) :: dx, dy
    integer :: i, j

    dx = settings%length_x / real(settings%nx - 1, real64)
    dy = settings%length_y / real(settings%ny - 1, real64)

    ! Thom's boundary relation estimates wall vorticity from nearby psi values
    ! and the wall's tangential speed. The top wall is the only moving boundary.
    do i = 2, settings%nx - 1
      flow%vorticity(i, 1) = 2.0_real64 * flow%streamfunction(i, 2) / dy**2
      flow%vorticity(i, settings%ny) = 2.0_real64 * flow%streamfunction(i, settings%ny - 1) / dy**2 &
        + 2.0_real64 * settings%lid_speed / dy
    end do
    do j = 2, settings%ny - 1
      flow%vorticity(1, j) = 2.0_real64 * flow%streamfunction(2, j) / dx**2
      flow%vorticity(settings%nx, j) = 2.0_real64 &
        * flow%streamfunction(settings%nx - 1, j) / dx**2
    end do

    ! The moving-lid velocity jumps at the two top corners; set corner
    ! vorticity to the adjacent wall average instead of inventing a corner flow.
    flow%vorticity(1, 1) = 0.5_real64 * (flow%vorticity(2, 1) + flow%vorticity(1, 2))
    flow%vorticity(settings%nx, 1) = 0.5_real64 &
      * (flow%vorticity(settings%nx - 1, 1) + flow%vorticity(settings%nx, 2))
    flow%vorticity(1, settings%ny) = 0.5_real64 &
      * (flow%vorticity(2, settings%ny) + flow%vorticity(1, settings%ny - 1))
    flow%vorticity(settings%nx, settings%ny) = 0.5_real64 &
      * (flow%vorticity(settings%nx - 1, settings%ny) &
      + flow%vorticity(settings%nx, settings%ny - 1))
  end subroutine apply_wall_vorticity

  subroutine calculate_velocity(settings, flow)
    type(simulation_settings), intent(in) :: settings
    type(flow_field), intent(inout) :: flow
    real(real64) :: dx, dy
    integer :: i, j

    dx = settings%length_x / real(settings%nx - 1, real64)
    dy = settings%length_y / real(settings%ny - 1, real64)
    flow%u = 0.0_real64
    flow%v = 0.0_real64

    do j = 2, settings%ny - 1
      do i = 2, settings%nx - 1
        flow%u(i, j) = (flow%streamfunction(i, j + 1) &
          - flow%streamfunction(i, j - 1)) / (2.0_real64 * dy)
        flow%v(i, j) = -(flow%streamfunction(i + 1, j) &
          - flow%streamfunction(i - 1, j)) / (2.0_real64 * dx)
      end do
    end do

    flow%u(2:settings%nx - 1, settings%ny) = settings%lid_speed
    ! TODO: The CSV is useful, but a streamline plot would make this flow way
    ! easier to eyeball. Compare it with a published Re=100 benchmark too.
  end subroutine calculate_velocity

  subroutine write_results(settings, flow)
    type(simulation_settings), intent(in) :: settings
    type(flow_field), intent(in) :: flow
    integer :: field_unit, profile_unit, io_status, i, j, center_i
    real(real64) :: dx, dy, x, y

    dx = settings%length_x / real(settings%nx - 1, real64)
    dy = settings%length_y / real(settings%ny - 1, real64)
    center_i = (settings%nx + 1) / 2

    open(newunit=field_unit, file='cavity_field.csv', status='replace', &
      action='write', iostat=io_status)
    if (io_status /= 0) call fail('Could not create cavity_field.csv.')
    write(field_unit, '(A)') 'x_m,y_m,u_m_s,v_m_s,streamfunction_m2_s,vorticity_1_s'
    do j = 1, settings%ny
      y = real(j - 1, real64) * dy
      do i = 1, settings%nx
        x = real(i - 1, real64) * dx
        write(field_unit, '(6(ES16.8,:,","))') x, y, flow%u(i, j), flow%v(i, j), &
          flow%streamfunction(i, j), flow%vorticity(i, j)
      end do
    end do
    close(field_unit)

    open(newunit=profile_unit, file='centerline_profile.csv', status='replace', &
      action='write', iostat=io_status)
    if (io_status /= 0) call fail('Could not create centerline_profile.csv.')
    write(profile_unit, '(A)') 'y_m,u_at_center_x_m_s,v_at_center_x_m_s'
    do j = 1, settings%ny
      y = real(j - 1, real64) * dy
      write(profile_unit, '(3(ES16.8,:,","))') y, flow%u(center_i, j), flow%v(center_i, j)
    end do
    close(profile_unit)
  end subroutine write_results

  subroutine fail(message)
    character(len=*), intent(in) :: message
    write(error_unit, '(A)') trim(message)
    write(error_unit, '(A)') 'Check simulation.nml and try again.'
    stop 1
  end subroutine fail

end module cavity_model
