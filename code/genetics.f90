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

  integer,parameter :: Nc = 50, Npop = 5000, time = 1000
  integer :: count, i, j, t
  integer :: population(Npop, Nc), new_pop(Npop, Nc)
  real(8) :: city(Nc, 2), min_distance(time), avg_distance(time)
  real(8) :: distance(Nc/10,10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call plot_init()
  call system_clock(count)
  call srand(count)
        
  do i = 10, Nc, 10
  do j = 1, 10
  call initiate(i, Npop, city, population)
  do t = 1, time
    call crossover(i, Npop, city, population, new_pop) 
    call mutation(i, Npop, new_pop)
    call selection(i, Npop, city, population, new_pop)
    avg_distance(t) = sum(calcdistance(i, Npop, population, city)) / Npop
    min_distance(t) = minval(calcdistance(i, Npop, population, city))
  end do
  distance(i, j) = min_distance(t)
  end do
  print*,i
  end do

  call plot_path(Nc, get_path(Nc, Npop, population, city), city)
  call plot_distance(avg_distance, min_distance, time)
  call plot_dist(distance, Nc/10)
!  call plot_close()
  call plend()
end program genetics
