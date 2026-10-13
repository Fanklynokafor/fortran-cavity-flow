program cavity_flow
  use, intrinsic :: iso_fortran_env, only: real64, output_unit, error_unit
  use cavity_model, only: simulation_settings, flow_field, read_settings, &
                          initialize_flow, advance_flow, write_results
  implicit none

  type(simulation_settings) :: settings
  type(flow_field) :: flow
  integer :: step, report_interval
  real(real64) :: time_now

  ! Keep the run sequence visible here; the numerical details live in the module.
  call read_settings('simulation.nml', settings)
  call initialize_flow(settings, flow)
  report_interval = max(1, settings%steps / 10)

  write(output_unit, '(A)') '2D lid-driven cavity simulation'
  write(output_unit, '(A,I0,A,I0)') 'Grid: ', settings%nx, ' x ', settings%ny
  write(output_unit, '(A,ES10.3)') 'Reynolds number: ', settings%reynolds_number
  write(output_unit, '(A,I0)') 'Time steps: ', settings%steps

  do step = 1, settings%steps
    time_now = real(step, real64) * settings%time_step
    call advance_flow(settings, flow)

    if (mod(step, report_interval) == 0 .or. step == settings%steps) then
      write(output_unit, '(A,I0,A,I0,A,F8.3,A,ES10.3)') &
        'Step ', step, '/', settings%steps, '  time=', time_now, &
        '  max |u|=', maxval(abs(flow%u))
    end if
  end do

  call write_results(settings, flow)
  if (flow%poisson_nonconverged > 0) then
    write(error_unit, '(A,I0,A,ES10.3)') &
      'Warning: Poisson iteration cap reached ', flow%poisson_nonconverged, &
      ' times; largest capped update was ', flow%poisson_worst_change
  end if
  write(output_unit, '(A)') 'Finished. See cavity_field.csv and centerline_profile.csv.'
end program cavity_flow
