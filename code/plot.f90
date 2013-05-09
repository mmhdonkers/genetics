module plot
  use plplot

  implicit none

  private numtostr
  public plot_init, plot_close, plot_spin, plot_path, plot_distance

contains

  subroutine plot_init() 
    call plsdev("xcairo")

    call plscol0(0, 255, 255, 255)  ! white
    call plscol0(1, 255, 0, 0)      ! red
    call plscol0(2, 0, 255, 0)      ! green
    call plscol0(3, 0, 0, 255)      ! blue
    call plscol0(4, 255, 0, 255)    ! magenta
    call plscol0(5, 0, 255, 255)    ! cyan
    call plscol0(6, 255, 255, 0)    ! yellow
    call plscol0(7, 0, 0, 0)        ! black
    call plscol0(8, 255, 77, 0)     ! orange
    call plscol0(9, 128, 128, 128)  ! gray

    call plinit()
  end subroutine

  subroutine plot_close()
    call plspause(.false.)
    call plend()
  end subroutine

  subroutine plot_spin(spin, SIZE, temp)
    integer,intent(in) :: SIZE
    integer,intent(in) :: spin(0:SIZE-1,0:SIZE-1)
    real(8),intent(in) :: temp

    integer :: i, j

    call plcol0(7)
    call plenv(0d0, SIZE*1d0, 0d0, SIZE*1d0, 0, 0)
    call pllab("x", "y", "coupling constant: " // numtostr(1d0/temp))

    do i = 0, SIZE - 1
      do j = 0, SIZE - 1
        if (spin(i,j) .eq. -1) then
          call plcol0(1)
          call plpoin([i + 0.5d0], [j + 0.5d0], 31)
        else if (spin(i,j) .eq. 1) then
          call plcol0(2)
          call plpoin([i + 0.5d0], [j + 0.5d0], 30)
        end if
      end do
    end do

    call plspause(.false.)
  end subroutine

  subroutine plot_profile(Lx, Ly, vel)
    integer,intent(in) :: Lx, Ly
    real(8),intent(in) :: vel(Lx, Ly)

    integer :: i
    real(8) :: y(Ly)

    y = [(i, i = 1, Ly)]

    call plcol0(7)
    call plenv(0d0, maxval(vel), 1d0, Ly * 1d0, 0, 0)
    call pllab("v", "y", "velocity profile")
    
    call plcol0(1)
    call plline(vel(Lx,:),y)

!    call plspause(.false.)
  end subroutine

  character(len=25) function numtostr(num) result(str)
    real(8),intent(in) :: num

    write(str, '(g12.5)') num
  end function

  subroutine plot_path(Nc, path, city)
    integer,intent(in) :: Nc, path(Nc)
    real(8),intent(in) :: city(Nc, 2)

    integer :: i
    real(8) :: x(2), y(2)

    call plcol0(7)
    call plenv(0d0, 1d0, 0d0, 1d0, 0, 0)
    call pllab("x", "y", "Path")

    call plcol0(1)
    x = [city(path(1), 1), city(path(Nc), 1)]
    y = [city(path(1), 2), city(path(Nc), 2)]
    call plline(x, y)
    do i = 2, Nc
      x = [city(path(i - 1), 1), city(path(i), 1)]
      y = [city(path(i - 1), 2), city(path(i), 2)]
      call plline(x, y)
    end do
    
    call plcol0(2)
    call plpoin(city(:, 1), city(:, 2), 1)
    
    call plspause(.false.)
  end subroutine

  subroutine plot_distance(avg_distance, min_distance, time)
    integer,intent(in) :: time
    real(8),intent(in) :: avg_distance(time), min_distance(time)

    integer :: i
    real(8) :: plottime(time)

    plottime = [(i * 1d0, i = 1, time)]

    call plcol0(7)
    call plenv(1d0, time * 1d0, 0d0, maxval(avg_distance), 0, 0)
    call pllab("t", "distance", "average and minimum distance")

    call plcol0(1)
    call plline(plottime, min_distance)

    call plcol0(2)
    call plline(plottime, avg_distance)

    call plspause(.true.)

    print*,min_distance(time)
  end subroutine
    
  subroutine plot_dist(distance, N)
    integer,intent(in) :: N
    real(8),intent(in) :: distance(N, 10)

    integer :: i
    real(8) :: plotN(N), plotD(N)

    plotN = [(i * 1d0, i = 1, N)]
    do i = 1, N
      plotD(i) = sum(distance(i, :)) / (N * 1d0)
    end do 

    call plcol0(7)
    call plenv(1d0, N * 1d0, 0d0, maxval(distance), 0, 0)
    call pllab("# cities", "average distance", "average distance")

    call plcol0(1)
    call plline(plotN, plotD)
  end subroutine
end module
