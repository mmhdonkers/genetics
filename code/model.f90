! Name:      Michiel Bakker and Chiel Donkers
! Course:    International Course on Computational Physics
! Project:   Genetics algorithm

module model

  implicit none

  private calcfitness
  public initiate, crossover

contains

  subroutine initiate(Nc, Npop, city, population)
    integer,intent(in) :: Nc, Npop
    integer,intent(out) :: population(Npop, Nc)
    real(8),intent(out) :: city(Nc, 2)

    integer :: i, j, k, temp_pop(Nc), temp    

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
  end subroutine
  
  function calcfitness(Nc, Npop, population, city) result(cum_fitness)
    integer,intent(in) :: Nc, Npop, population(Npop, Nc)
    real(8),intent(in) :: city(Nc, 2)
  
    integer :: i, j
    real(8) :: distance(Npop), fitness(Npop), cum_fitness(Npop)

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
  end function

  subroutine crossover(Nc, Npop, city, population, new_pop)
  integer,intent(in) :: Nc, Npop, population(Npop, Nc)
  real(8),intent(in) :: city(Nc, 2)
  integer,intent(out) :: new_pop(Npop, Nc)

  integer :: i, j, k, r1, r2, parents(2, Nc)
  real(8) :: r, cum_fitness(Npop)
  
  new_pop = -1
  cum_fitness = calcfitness(Nc, Npop, population, city)

  do i = 1, Npop
    ! Select random parents
    r = rand()
    j = 1
    do while (r .gt. cum_fitness(j))
      j = j + 1
    end do
    parents(1, :) = population(j, :)
    r = rand()
    j = 1
    do while (r .gt. cum_fitness(j))
      j = j + 1
    end do
    parents(2, :) = population(j, :)

    ! Create children
    r1 = floor(rand() * Nc) + 1
    r2 = floor(rand() * Nc) + 1
    if (r1 .lt. r2) then
      new_pop(i,r1:r2) = parents(1, r1:r2)
    else
      new_pop(i,r2:r1) = parents(1, r2:r1) 
    end if

    do j = 1, Nc
      if (.not.any(parents(2,j) .eq. new_pop(i,:))) then
        k = 1
        do while (new_pop(i,k) .ne. -1)
          k = k + 1
        end do
        new_pop(i,k) = parents(2,j)
      end if
    end do
  end do
  end subroutine
end module
