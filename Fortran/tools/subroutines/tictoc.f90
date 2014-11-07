subroutine tic()
    call system_clock(count_rate=clock_rate)
    call system_clock(start_clock)
    call cpu_time(cpu_time_start)
end subroutine tic

subroutine toc(times)
    real(kind=8), intent(inout) :: times(2)
    call cpu_time(cpu_time_end)
    call system_clock(end_clock)
    times(1) = cpu_time_end - cpu_time_start
    times(2) = (end_clock - start_clock)/real(clock_rate,8)
end subroutine toc