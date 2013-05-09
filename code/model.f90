! Name:      Michiel Bakker and Chiel Donkers
! Course:    International Course on Computational Physics
! Project:   Genetics algorithm

module model

  implicit none

  private calccumfitness, calcfitness
  public initiate, crossover, mutation, get_path, calcdistance, selection, &
         readdata, writedata

contains

  subroutine initiate(Nc, Npop, city, population, cout)
    integer,intent(in) :: Nc, Npop, cout
    integer,intent(out) :: population(Npop, Nc)
    real(8),intent(out) :: city(Nc, 2)

    integer :: i, j, k, temp_pop(Nc), temp    

    call srand(3)
    ! Generate city locations
    do i = 1, Nc
      city(i,:) = [rand(), rand()]
    end do

    call srand(cout)
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
  
  function calcdistance(Nc, Npop, population, city) result(distance)
    integer,intent(in) :: Nc, Npop, population(Npop, Nc)
    real(8),intent(in) :: city(Nc, 2)
  
    integer :: i, j
    real(8) :: distance(Npop)

    ! Calculate cumulative fitness
    distance = 0
    do i = 1, Npop
      distance(i) = sqrt((city(population(i, 1), 1) - city(population(i, Nc), 1))**2 &
                    + (city(population(i, 1), 2) - city(population(i, Nc), 2))**2)
      do j = 2, Nc
        distance(i) = distance(i) + sqrt((city(population(i, j - 1), 1) - city(population(i, j), 1))**2 &
                    + (city(population(i, j - 1), 2) - city(population(i, j), 2))**2)
      end do
    end do
  end function

  function calcfitness(Nc, Npop, population, city) result(fitness)
    integer,intent(in) :: Nc, Npop, population(Npop, Nc)
    real(8),intent(in) :: city(Nc, 2)
  
    integer :: i, j
    real(8) :: distance(Npop), fitness(Npop)

    ! Calculate cumulative fitness
    distance = 0
    do i = 1, Npop
      distance(i) = sqrt((city(population(i, 1), 1) - city(population(i, Nc), 1))**2 &
                    + (city(population(i, 1), 2) - city(population(i, Nc), 2))**2)
      do j = 2, Nc
        distance(i) = distance(i) + sqrt((city(population(i, j - 1), 1) - city(population(i, j), 1))**2 &
                    + (city(population(i, j - 1), 2) - city(population(i, j), 2))**2)
      end do
    end do
    fitness = 1 / distance
  end function

  function calccumfitness(Nc, Npop, population, city) result(cum_fitness)
    integer,intent(in) :: Nc, Npop, population(Npop, Nc)
    real(8),intent(in) :: city(Nc, 2)
  
    integer :: i, j
    real(8) :: distance(Npop), fitness(Npop), cum_fitness(Npop)

    ! Calculate cumulative fitness
    distance = 0
    do i = 1, Npop
      distance(i) = sqrt((city(population(i, 1), 1) - city(population(i, Nc), 1))**2 &
                    + (city(population(i, 1), 2) - city(population(i, Nc), 2))**2)
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
    cum_fitness = calccumfitness(Nc, Npop, population, city)

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
        new_pop(i, r1:r2) = parents(1, r1:r2)
      else
        new_pop(i, r2:r1) = parents(1, r2:r1) 
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

  subroutine selection(Nc, Npop, city, population, new_pop)
    integer,intent(in) :: Nc, Npop, new_pop(Npop, Nc)
    real(8),intent(in) :: city(Nc, 2)
    integer,intent(inout) :: population(Npop, Nc)

    integer :: i, j, temp_pop(2 * Npop, Nc)
    real(8) :: r, cum_fitness(2 * Npop), fitness(2 * Npop)

    cum_fitness(1:Npop) = calccumfitness(Nc, Npop, population, city)
    cum_fitness(Npop + 1:2 * Npop) = 1 + calccumfitness(Nc, Npop, new_pop, city)
    fitness(1:Npop) = calcfitness(Nc, Npop, population, city)
    fitness(Npop + 1:2 * Npop) = calcfitness(Nc, Npop, population, city)
    temp_pop(1:Npop, :) = population
    temp_pop(Npop + 1:2 * Npop, :) = new_pop

    do i = 1, Npop / 20 
      population(i, :) = temp_pop(maxloc(fitness, 1), :)
    end do

    do i = Npop / 20 + 1, Npop
      r = rand() * 2
      j = 1
      do while (r .gt. cum_fitness(j))
        j = j + 1
      end do
      population(i, :) = temp_pop(j, :)
    end do
  end subroutine

  subroutine mutation(Nc, Npop, population, m)
    integer,intent(in) :: Nc, Npop, m
    integer,intent(inout) :: population(Npop, Nc)

    integer :: i, j, k, temp, r1, r2, r3, new_pop(Nc)
    real(8) :: r

    do i = 1, Npop
      r = rand()
      if (r .lt. m / 100d0) then
        r1 = floor(rand() * Nc) + 1
        r2 = floor(rand() * Nc) + 1
        temp = population(i, r1)
        population(i, r1) = population(i, r2)
        population(i, r2) = temp
      end if
      r = rand()
      if (r .lt. m / 100d0) then
        r1 = floor(rand() * Nc) + 1
        r2 = floor(rand() * Nc) + 1
        new_pop = population(i, :)
        do j = 0, abs(r1 - r2)
          population(i, min(r1, r2) + j) = new_pop(max(r1, r2) - j)
        end do
      end if
      r = rand()
      if (r .lt. m / 100d0) then
        new_pop = -1
        r1 = floor(rand() * Nc) + 1
        r2 = floor(rand() * Nc) + 1
        r3 = floor(rand() * (Nc - abs(r1 - r2))) + 1
        do j = 0, abs(r1 - r2)
          new_pop(r3 + j) = population(i, min(r1, r2) + j)
        end do

        do j = 1, Nc
          if (.not.any(population(i, j) .eq. new_pop)) then
            k = 1
            do while (new_pop(k) .ne. -1)
              k = k + 1
            end do
            new_pop(k) = population(i, j)
          end if
        end do
        population(i, :) = new_pop
      end if
    end do
  end subroutine

  function get_path(Nc, Npop, population, city) result(path)
    integer,intent(in) :: Nc, Npop, population(Npop, Nc)
    real(8),intent(in) :: city(Nc)

    integer :: path(Nc)
    real(8) :: fitness(Npop)

    fitness = calcfitness(Nc, Npop, population, city)

    path = population(maxloc(fitness, 1), :)
  end function
  
  subroutine readdata(city, Nc)
    integer,intent(in) :: Nc
    real(8),intent(out) :: city(Nc,2)

    integer :: i
    real(8) ::  dummycity(Nc, 3)
      
    open(unit=14, file='Qatar.dat', status='old', action='read')
      
    do i=1,Nc
      read(14,*),dummycity(i, :)
      city(i, 1) = dummycity(i, 3)
      city(i, 2) = dummycity(i, 2)
    end do
      
    close(14)  
  end subroutine

  subroutine writedata(Nc, Npop, population, city, avg_distance, min_distance, time, N, l)
    integer,intent(in) :: Nc, Npop, population(Npop, Nc), time, N, l
    real(8),intent(in) :: city(Nc, 2), avg_distance(time), min_distance(time)

    integer :: i
    real(8) :: fitness(Npop)

    fitness = calcfitness(Nc, Npop, population, city)

    open(unit=15, file='path_' // trim(numtostr(N)) // '_' // trim(numtostr(Npop)) // '_' // &
              trim(numtostr(time)) // '_' // trim(numtostr(l)) // '.dat', status='replace')
    open(unit=16, file='dist_' // trim(numtostr(N)) // '_' // trim(numtostr(Npop)) // '_' // &
              trim(numtostr(time)) // '_' // trim(numtostr(l)) // '.dat', status='replace')
    open(unit=17, file='city.dat', status='replace')

    do i = 1, Nc
      write(15,*) population(maxloc(fitness, 1), i)
    end do

    do i = 1, time
      write(16,*) min_distance(i), avg_distance(i)
    end do

    do i = 1, Nc
      write(17,*) city(i, 1), city(i, 2)
    end do

    close(15)
    close(16)
    close(17)
  end subroutine
  
  character(len=25) function numtostr(num) result(str)
    integer,intent(in) :: num

    write(str, '(I5.2)') num
  end function
end module
