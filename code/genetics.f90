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

  integer,parameter :: Nc = 50, N = 10
  integer :: count, t, i, j, k, m, time, Npop
  integer,allocatable :: population(:, :), new_pop(:, :)
  real(8) :: city(Nc, 2)
  real(8),allocatable ::  min_distance(:), avg_distance(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call system_clock(count)
  call srand(count)
        
do i = 1, N
do j = 1, 40
Npop = 250 * j
allocate(population(Npop, Nc))
allocate(new_pop(Npop, Nc))
do k = 1, 40
time = 250 * k
allocate(min_distance(time))
allocate(avg_distance(time))
do m = 0, 50, 5
  call initiate(Nc, Npop, city, population, count)
  !call readdata(city, Nc)
  do t = 1, time
    call crossover(Nc, Npop, city, population, new_pop) 
    call mutation(Nc, Npop, new_pop, m)
    call selection(Nc, Npop, city, population, new_pop)
    avg_distance(t) = sum(calcdistance(Nc, Npop, population, city)) / Npop
    min_distance(t) = minval(calcdistance(Nc, Npop, population, city))
  end do
  print*, i, j, k, m
  call writedata(Nc, Npop, population, city, avg_distance, min_distance, time, i, m)
end do
end do
end do
end do
end program genetics
