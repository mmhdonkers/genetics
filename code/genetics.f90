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

  integer,parameter :: Nc = 194, Npop = 2000, time = 200000
  integer :: count, t, population(Npop, Nc), new_pop(Npop, Nc)
  real(8) :: city(Nc, 2), min_distance(time), avg_distance(time)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call system_clock(count)
  call srand(count)
        
  call initiate(Nc, Npop, city, population, count)
  call readdata(city, Nc)
  do t = 1, time
    call crossover(Nc, Npop, city, population, new_pop) 
    call mutation(Nc, Npop, new_pop, 50)
    call selection(Nc, Npop, city, population, new_pop)
    call mutation(Nc, Npop, new_pop, 30)
    call selection(Nc, Npop, city, population, new_pop)
    call mutation(Nc, Npop, new_pop, 10)
    call selection(Nc, Npop, city, population, new_pop)
    avg_distance(t) = sum(calcdistance(Nc, Npop, population, city)) / Npop
    min_distance(t) = minval(calcdistance(Nc, Npop, population, city))
    if (mod(time, 50) == 0) call plot_path(Nc, population(maxloc(calcfitness(Nc, Npop, population, city), 1), :), city)
    print*, t
  end do
  
  call writedata(Nc, Npop, population, city, avg_distance, min_distance, time)
end program genetics
