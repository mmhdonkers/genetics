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

  integer,parameter :: Nc = 50, Npop = 1000, time = 200
  integer :: count, t, N, j
  integer :: population(Npop, Nc), new_pop(Npop, Nc)
  real(8) :: city(Nc, 2), min_distance(time), avg_distance(time)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call system_clock(count)
  call srand(count)
        
N = 10
do j = 1, N
  call initiate(Nc, Npop, city, population, count)
  !call readdata(city, Nc)
  do t = 1, time
    print*, t
    call crossover(Nc, Npop, city, population, new_pop) 
    call mutation(Nc, Npop, new_pop)
    call selection(Nc, Npop, city, population, new_pop)
    avg_distance(t) = sum(calcdistance(Nc, Npop, population, city)) / Npop
    min_distance(t) = minval(calcdistance(Nc, Npop, population, city))
  end do
  call writedata(Nc, Npop, population, city, avg_distance, min_distance, time, j)
end do
end program genetics
