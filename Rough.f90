MODULE rough_cluster_module
    
    USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
    
    IMPLICIT NONE
    
    PRIVATE
    
    PUBLIC :: read_config_atoms, allocate_arrays, deallocate_arrays, shape_rough_cluster, ran3, gaussian
    
    CONTAINS
    
        SUBROUTINE read_config_atoms(filename,n,r,v)
            IMPLICIT NONE
            
            CHARACTER(len=*), INTENT(in) :: filename
            REAL*8, DIMENSION(:,:), OPTIONAL, INTENT(out) :: r, v            
            INTEGER*8, INTENT(inout) :: n
            
            INTEGER :: cnf_unit, ioerr, i
            
            ! Read r and v if they exist
            IF (PRESENT(r) .and. PRESENT (v)) THEN
                
                ! Open the position file
                OPEN(newunit=cnf_unit,file=filename//'_pos.dat',status='old',action='read',iostat=ioerr)
                
                IF (ioerr /= 0) THEN
                    WRITE(unit=error_unit,fmt='(a,a,i15)') 'Error opening ', filename//'_pos.dat', ioerr
                    READ(unit=input_unit,fmt=*)
                    STOP 'Error in read_config_atoms'
                END IF
                
                ! Read the position file
                READ(unit=cnf_unit,fmt=*,iostat=ioerr)
                
                DO i = 1, n
                    
                    READ(unit=cnf_unit,fmt=*,iostat=ioerr) r(:,i)
                    
                    IF (ioerr /= 0) THEN
                        
                        WRITE(unit=error_unit,fmt='(a,a,i15)') 'Error reading r from ', filename//'_pos.dat', ioerr
                        READ(unit=input_unit,fmt=*)
                        
                        IF (ioerr == iostat_eor) THEN
                            WRITE(unit=error_unit, fmt='(a)') 'End of record'
                            READ(unit=input_unit,fmt=*)
                        END IF
                        
                        IF (ioerr == iostat_end) THEN
                            WRITE(unit=error_unit, fmt='(a)') 'End of file'
                            READ(unit=input_unit,fmt=*)
                        END IF
                        
                        STOP 'Error in read_config_atoms'
                        
                    END IF
                END DO
                
                ! Close the position file
                CLOSE (unit=cnf_unit)
                
                ! Open the velocity file
                OPEN(newunit=cnf_unit,file=filename//'_vel.dat',status='old',action='read',iostat=ioerr)
                
                IF (ioerr /= 0) THEN
                    WRITE(unit=error_unit,fmt='(a,a,i15)') 'Error opening ', filename//'_vel.dat', ioerr
                    READ(unit=input_unit,fmt=*)
                    STOP 'Error in read_config_atoms'
                END IF
                
                ! Read the velocity file
                READ(unit=cnf_unit,fmt=*,iostat=ioerr)
                
                DO i = 1, n
                    
                    READ (unit=cnf_unit,fmt=*,iostat=ioerr) v(:,i)
                    
                    IF (ioerr /= 0) THEN
                        
                        WRITE(unit=error_unit,fmt='(a,a,i15)') 'Error reading v from ', filename//'_vel.dat', ioerr
                        
                        IF (ioerr == iostat_eor) THEN
                            WRITE(unit=error_unit,fmt='(a)') 'End of record'
                            READ(unit=input_unit,fmt=*)
                        END IF
                        
                        IF (ioerr == iostat_end) THEN
                            WRITE(unit=error_unit,fmt='(a)') 'End of file'
                            READ(unit=input_unit,fmt=*)
                        END IF
                        
                        STOP 'Error in read_config_atoms'
                    
                    END IF
                
                END DO
                
                ! Close the velocity file
                CLOSE (unit=cnf_unit)
                
            ELSE
                
                ! Read n
                OPEN(newunit=cnf_unit,file=filename//'_pos.dat',status='old',action='read',iostat=ioerr)
                
                READ(unit=cnf_unit,fmt=*,iostat=ioerr) n
               
                IF (ioerr /= 0) THEN
                    
                    WRITE(unit=error_unit,fmt='(a,a,i15)') 'Error reading n from ', filename//'_pos.dat', ioerr
                    
                    READ(unit=input_unit,fmt=*)
                    
                    IF (ioerr == iostat_eor) THEN
                        WRITE(unit=error_unit,fmt='(a)') 'End of record'
                        READ(unit=input_unit,fmt=*)
                    END IF
                   
                    IF (ioerr == iostat_end) THEN
                        WRITE(unit=error_unit,fmt='(a)') 'End of file'
                        READ(unit=input_unit,fmt=*)
                    END IF
                    
                    STOP 'Error in read_config_atoms'
               
                END IF
            
            END IF
            
        END SUBROUTINE read_config_atoms
    
        SUBROUTINE allocate_arrays(n,r,v,f)
        
            IMPLICIT NONE
            
            INTEGER*8, INTENT(in) :: n ! Number of particles
            REAL*8, DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: r, v, f ! Position, velocity, and force arrays
            
            ALLOCATE(r(3,n),v(3,n),f(3,n))

        END SUBROUTINE allocate_arrays
        
        SUBROUTINE deallocate_arrays(r,v,f)
        
            IMPLICIT NONE
            
            REAL*8, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: r, v, f
        
            DEALLOCATE(r,v,f)
            
        END SUBROUTINE deallocate_arrays        
        
        SUBROUTINE shape_rough_cluster(n,desired_sigma_percent,r,v,n_mantle,mantle_array)
        
            IMPLICIT NONE
            
            INTEGER*8 :: i, j, num_surface_target, random_number, n_original
            INTEGER*8, INTENT(inout) :: n
            INTEGER*8, DIMENSION(n) :: track
            REAL*8, INTENT(in) :: desired_sigma_percent
            REAL*8 :: mantle_percent, r_max, r_core_max
            REAL*8 :: r_gaussian, pi, theta_random, phi_random
            REAL*8 :: sigma_target, num_sample, sum_y, sum_y_sq, avg_y, avg_y_sq, sigma
            REAL*8, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: r, v
            REAL*8, DIMENSION(n) :: r_modulus, height
            REAL*8, DIMENSION(:), ALLOCATABLE :: r_random_array
            REAL*8, DIMENSION(:,:), ALLOCATABLE :: temp
            INTEGER*8, INTENT(out) :: n_mantle
            REAL*8, DIMENSION(:), ALLOCATABLE, INTENT(out) :: mantle_array
            REAL*8 :: pi_agnost
            
            REAL*8 :: radius_mantle_avg, radius_mantle_avg_sq, rho_mantle_avg, radius_new_max, radius_new_max_sq
                       
            n_original = n
            
            pi_agnost = 4.0d0*ATAN(1.0d0)
            
            random_number = 10250501230
            
            mantle_percent = desired_sigma_percent
            
            ! Find the maximum r value
            DO i=1,n
                r_modulus(i) = (r(1,i)*r(1,i)) + (r(2,i)*r(2,i)) + (r(3,i)*r(3,i))
            END DO
            r_modulus = DSQRT(r_modulus)
            r_max = MAXVAL(r_modulus)
            
            ! Set the core's maximum r value at a % of the cluster's max r value
            r_core_max = (1.0d0 - mantle_percent) * r_max
            
            ! Set the target sigma (and thus the targeted 68% of the total thickness
            ! for the normalized surface)
            sigma_target = mantle_percent * r_max
                       
            ! Identify mantle particles and number of particles
            track = 0
            num_sample = 0
            DO i=1,n
                IF(r_modulus(i) > r_core_max) THEN
                    track(num_sample+1) = i       
                    num_sample = num_sample + 1
                END IF
            END DO
            
            ! Determine the average radius of the mantle, and then the reduced surface particle density for that radius
            radius_mantle_avg = ((r_max - r_core_max)/2.0d0)+r_core_max
            radius_mantle_avg_sq = radius_mantle_avg**2
            rho_mantle_avg = num_sample/radius_mantle_avg_sq
            
            ! Determine the new max radius of the mantle, and then set the target number of particles that matches
            ! the prior calculated reduced surface particle density
            
            radius_new_max = (sigma_target/0.68d0)+r_max
            radius_new_max_sq = radius_new_max**2
            
            num_surface_target = DINT(rho_mantle_avg*radius_new_max_sq)

            ! Allocate a new random radial array for the new particles
            ALLOCATE(r_random_array(num_surface_target))
            
            ! Generate a normally distributed random radii, based on num_surface_target
            ! Scale this to the desired thickness, and add it on top of r_max
            CALL gaussian(random_number,num_surface_target,r_random_array)
            r_random_array = r_random_array * (sigma_target/0.68d0)
            r_random_array = r_random_array + r_max
            OPEN(unit=2,file='./Output/stat.dat',status='replace')
            DO i=1,num_surface_target
                WRITE(unit=2,fmt=*)r_random_array(i)
            END DO
            CLOSE(unit=2)
            
            ! Grow the original position array so that we can add these new particles
            ALLOCATE(temp(3,n+num_surface_target))
            DO i=1,n
                temp(:,i) = r(:,i)
            END DO
            temp(:,n+1:num_surface_target) = 0.0d0 ! Dummy value; this will be changed later in the routine
            DEALLOCATE(r)
            ALLOCATE(r(3,n+num_surface_target))
            r = temp
            DEALLOCATE(temp)
            
            ! Grow the velocity array to accomidate
            ALLOCATE(temp(3,n+num_surface_target))
            DO i=1,n
                temp(:,i) = v(:,i)
            END DO
            temp(:,n+1:n+num_surface_target) = 0.0d0 ! Dummy value; this will be changed later in the routine
            DEALLOCATE(v)
            ALLOCATE(v(3,n+num_surface_target))
            v = temp
            DEALLOCATE(temp)
        
            ! Copy over the proper position array          
            DO i=1,num_surface_target
                theta_random = ran3(random_number-123256897+i)*pi_agnost
                phi_random = ran3(random_number+i)*2.0d0*pi_agnost
                r(1,n+i) = r_random_array(i) * DSIN(theta_random) * DCOS(phi_random)
                r(2,n+i) = r_random_array(i) * DSIN(theta_random) * DSIN(phi_random)
                r(3,n+i) = r_random_array(i) * DCOS(theta_random)
            END DO
            
            ! Calculate sigma for the new cluster level
            height = 0
            num_sample = 0.0d0
            sum_y = 0.0d0
            sum_y_sq = 0.d0
            
            DO i=1,num_surface_target
                height(i) = (r(1,n+i)**2) + (r(2,n+i)**2) + (r(3,n+i)**2)
                height(i) = height(i) - r_max
                sum_y = sum_y + height(i)
                sum_y_sq = sum_y_sq + height(i)**2
                num_sample = num_sample + 1
            END DO
            
            avg_y = sum_y/num_sample
            avg_y_sq = sum_y_sq/num_sample
            sigma = DSQRT(avg_y_sq - avg_y**2)
            
            OPEN(unit=2,file='./Output/stat.dat',status='replace')
            WRITE(unit=2,fmt=*)'Number of samples: ',num_sample
            WRITE(unit=2,fmt=*)'Average delta y: ',avg_y
            WRITE(unit=2,fmt=*)'Average delta y sq: ',avg_y_sq
            WRITE(unit=2,fmt=*)'Sigma: ',sigma
            WRITE(unit=2,fmt=*)'Mantle delta y values below'
            DO i=1,num_surface_target
                WRITE(unit=2,fmt=*) height(i)
            END DO
            CLOSE(unit=2)
            
            ! Change to new n
            n = n + num_surface_target
            
            ! Report back number of mantle particles and new mantle array
            n_mantle = num_surface_target + num_sample
            ALLOCATE(mantle_array(n))
            mantle_array = 0
            j=1
            ! Add original mantle atoms
            DO i=1,num_sample
                IF(track(i) /= 0) THEN
                    index = track(i)
                    mantle_array(j) = index
                    j = j + 1
                END IF
            END DO
            ! Add new mantle atoms
            DO i=1,num_surface_target
                mantle_array(j) = n_original + i
            END DO
                
                
                

        END SUBROUTINE shape_rough_cluster
        
        SUBROUTINE gaussian(iseed,npart,a)
        ! 'a' is the array we want Gaussian numbers for
        ! 
        ! Note - this is scaled to 1. You'll need to rescale to your desired length scale
            IMPLICIT NONE
            INTEGER*8 :: i, iseed, idum, npart
            REAL*8 :: u1, u2, v1, v2, s, r
            REAL*8, DIMENSION(npart) :: a
            
            idum = -2*iseed-1
      
            ! Box-Muller transform for Gaussian numbers
            DO i=1,npart,2
1               u1 = ran3(idum)
                u2 = ran3(idum)
                v1 = 2.0d0*u1-1.0d0
                v2 = 2.0d0*u2-1.0d0
                s = v1*v1 + v2*v2
                IF (s .gt. 1.0) GO TO 1
                r = -2.0d00*DLOG(s)/s
                a(i) = v1*DSQRT(r)
                IF(i .eq. npart) EXIT
            ! Note: if npart is odd, loop goes to
            ! i=npart, so vel(n+1) is outside of the array
            ! and the program breaks. I created an exit
            ! loop condition to avoid this.
                a(i+1) = v2*DSQRT(r)
            END DO
        
        END SUBROUTINE Gaussian
        
        FUNCTION ran3(idum)
        ! ======================================================================
        ! This function generates a random number based on the idum seed.
        ! ======================================================================
              implicit none
              INTEGER*8 idum
              double precision MBIG,MSEED,MZ
              double precision ran3,FAC
              PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
              INTEGER*8 i,iff,ii,inext,inextp,k
              double precision mj,mk,ma(55)
              SAVE iff,inext,inextp,ma
              DATA iff /0/
      
              if(idum.lt.0.or.iff.eq.0)then
                iff=1
                mj=MSEED-iabs(idum)
                mj=mod(mj,MBIG)
                ma(55)=mj
                mk=1
                do 11 i=1,54
                  ii=mod(21*i,55)
                  ma(ii)=mk
                  mk=mj-mk
                  if(mk.lt.MZ)mk=mk+MBIG
                  mj=ma(ii)
11      continue
                do 13 k=1,4
                  do 12 i=1,55
                    ma(i)=ma(i)-ma(1+mod(i+30,55))
                    if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
                inext=0
                inextp=31
                idum=1
              endif
              inext=inext+1
              if(inext.eq.56)inext=1
              inextp=inextp+1
              if(inextp.eq.56)inextp=1
              mj=ma(inext)-ma(inextp)
              if(mj.lt.MZ)mj=mj+MBIG
              ma(inext)=mj
              ran3=mj*FAC
      
              return
        END FUNCTION ran3        
        
END MODULE rough_cluster_module