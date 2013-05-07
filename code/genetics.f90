! Name:      Michiel Bakker and Chiel Donkers
! Course:    International Course on Computational Physics
! Project:   Genetics algorithm
! 
!Program Summary
!
!  Input:    
!
!  Process:  
!
!  Output:   
!



program genetics

  use model
  use plot

  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer,parameter :: Nc = 25, Npop = 100, time = 10000
  integer :: count, i
  integer :: population(Npop, Nc), new_pop(Npop, Nc)
  real(8) :: city(Nc, 2), min_distance(time), avg_distance(time)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call plot_init()
  call system_clock(count)
  call srand(count)
        
  call initiate(Nc, Npop, city, population)
  do i = 1, time
    call crossover(Nc, Npop, city, population, new_pop) 
    call mutation(Nc, Npop, new_pop)
    call selection(Nc, Npop, city, population, new_pop)
    avg_distance(i) = sum(calcdistance(Nc, Npop, population, city)) / Npop
    min_distance(i) = minval(calcdistance(Nc, Npop, population, city))
  end do

  call plot_path(Nc, get_path(Nc, Npop, population, city), city)
  call plot_distance(avg_distance, min_distance, time)
!  call plot_close()
  call plend()
end program genetics
