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

  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer,parameter :: Nc = 4, Npop = 10
  integer :: count
  integer :: i, j, k, population(Npop, Nc), temp_pop(Nc), temp, parents(2,
Nc),child(Nc)
  real(8) :: city(Nc, 2), distance(Npop), fitness(Npop), cum_fitness(Npop)
  real(8) :: r1, r2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  call system_clock(count)
  call srand(count)
 
  ! Generate city locations
  do i = 1, Nc
    city(i,:) = [rand(), rand()]
  end do

  ! Generate populations
  do i = 1, Npop
    temp_pop = [(j, j = 1, Nc)]
    do j = Nc, 1, -1
      k = floor(rand() * j) + 1
      temp = temp_pop(k)
      temp_pop(k) = temp_pop(j)
      temp_pop(j) = temp 
    end do
    population(i,:) = temp_pop
  end do
  
  ! Calculate cumulative fitness
  distance = 0
  do i = 1, Npop
    do j = 2, Nc
      distance(i) = distance(i) + sqrt((city(population(i, j - 1), 1) - city(population(i, j), 1))**2 &
                    + (city(population(i, j - 1), 2) - city(population(i, j), 2))**2)
    end do
  end do
  fitness = 1 / distance
  cum_fitness(1) = fitness(1) / sum(fitness)
  do i = 2, Npop
    cum_fitness(i) = cum_fitness(i - 1) + fitness(i) / sum(fitness)
  end do

  ! Select random parents
  r1 = rand()
  i = 1
  do while (r1 .gt. cum_fitness(i))
    i = i + 1
  end do
  parents(1, :) = population(i, :)
  r1 = rand()
  i = 1
  do while (r1 .gt. cum_fitness(i))
    i = i + 1
  end do
  parents(2, :) = population(i, :)

  r1 = floor(rand() * Nc) + 1
  r2 = floor(rand() * Nc) + 1
  Child=0
  if (r1 .lt. r2) then
    child(r1:r2)=parents(1,r1:r2)
    do i=1,r1-1
      do 
        if (Child .eq. 0 .and. parent(2,j) )
          Child(i)=parent(2,j)
        end if
      end do
    end do
  else
    
  end if


 print*, parents
end program genetics
