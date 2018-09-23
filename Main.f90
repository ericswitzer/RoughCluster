PROGRAM main
    USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
    USE rough_cluster_module, ONLY : read_config_atoms, allocate_arrays, deallocate_arrays, shape_rough_cluster, ran3, gaussian
    
    IMPLICIT NONE
  
    INTEGER*8                              :: i, n, nsteps, n_mantle
    REAL*8                                 :: rcutoff
    REAL*8,    DIMENSION(3)                :: r_com, v_com                 ! Dummy center of mass
    INTEGER*8, DIMENSION(:), ALLOCATABLE   :: mantle_array
    REAL*8,    DIMENSION(:,:), ALLOCATABLE :: r, v, f                      ! Postion, Velocity, Force
     
    
    ! Read number of particles
    CALL read_config_atoms('.\Input\cluster01',n)
    
    ! Allocate position, velocity, and force arrays
    CALL allocate_arrays(n,r,v,f)
    
    ! Read position and velocity arrays
    CALL read_config_atoms('.\Input\cluster01',n,r,v)
    
    ! Ensure the cluster's c.o.m. position is zero
    r_com(:) = SUM(r(:,:),dim=2)/DBLE(n) ! Calculate cluster c.o.m. position
    r(:,:) = r(:,:) - SPREAD(r_com(:),dim=2,ncopies=n) ! Move cluster to (0,0,0)

    ! Ensure the cluster's c.o.m. velocity is zero
    v_com(:) = SUM(v(:,:),dim=2)/DBLE(n) ! Calculate cluster's c.o.m. velocity
    v(:,:) = v(:,:) - SPREAD(v_com(:),dim=2,ncopies=n) ! Move cluster in velocity space to (0,0,0)
    
    ! Set initial forces equal to zero
    f = 0.0d0
    
    ! Shape the cluster to the sigma percent of mantle
    CALL shape_rough_cluster(n,0.1d0,r,v,n_mantle,mantle_array)  
    
    ! Write VMD snapshot for visual check
    OPEN(unit=2,file='./Output/cluster01.xyz',status='replace')
    WRITE(unit=2,fmt=*)n
    WRITE(unit=2,fmt=*)
    DO i=1,n
        WRITE(unit=2,fmt=*),'Ar ',r(1,i),r(2,i),r(3,i)
    END DO
    CLOSE(2)
    
    ! Write new position values
    OPEN(unit=2,file='./Output/cluster01_pos.dat',status='replace')
    WRITE(unit=2,fmt=*)n
    DO i=1,n
        WRITE(unit=2,fmt=*),r(1,i),r(2,i),r(3,i)
    END DO
    CLOSE(2)
    
    ! Write new velocity values
    OPEN(unit=2,file='./Output/cluster01_vel.dat',status='replace')
    WRITE(2,*)n
    DO i=1,n
        WRITE(unit=2,fmt=*),v(1,i),v(2,i),v(3,i)
    END DO
    CLOSE(2)

    ! Write mantle positions
    OPEN(unit=2,file='./Output/cluster01_mantle.dat',status='replace')
    WRITE(unit=2,fmt=*)n_mantle
    WRITE(unit=2,fmt=*)
    DO i=1,n_mantle
        WRITE(unit=2,fmt=*),mantle_array(i)
    END DO
    CLOSE(2)
    
    WRITE(unit=output_unit,fmt='(a)') 'Rough surface generation complete!'
    WRITE(unit=output_unit,fmt='(a)') 'Press any key to exit'
    READ(unit=input_unit,fmt=*)
    
    CALL deallocate_arrays(r,v,f)
    
    STOP
    
END PROGRAM main
