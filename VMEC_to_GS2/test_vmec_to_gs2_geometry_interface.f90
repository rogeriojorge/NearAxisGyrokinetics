program test_vmec_to_gs2_geometry_interface

  use vmec_to_gs2_geometry_interface_mod

  implicit none

  !*********************************************************************
  ! Input parameters
  !*********************************************************************

  character(len=2000) :: vmec_filename, gs2_output, geometry_output,number_of_field_periods_to_includestr, desired_normalized_toroidal_fluxstr
  integer, parameter :: nalpha = 1
  integer, parameter :: nzgrid = 851
  integer, parameter :: nlambda = 50
  real :: zeta_center = 0.0
  real :: number_of_field_periods_to_include, desired_normalized_toroidal_flux
  integer :: vmec_surface_option = 0
  logical :: verbose = .true.

  !*********************************************************************
  ! Output arrays
  !*********************************************************************
  
  real :: normalized_toroidal_flux_used, safety_factor_q, shat, L_reference, B_reference, bMax, bMin
  real, dimension(nalpha) :: alpha
  real, dimension(-nzgrid:nzgrid) :: zeta
  real, dimension(nalpha, -nzgrid:nzgrid) :: bmag, gradpar, gds2, gds21, gds22, gbdrift, gbdrift0, cvdrift, cvdrift0, grho
  ! This code uses normalizations in which kxfac is always 1, so kxfac is not presently returned.

  !*********************************************************************
  ! Variables used internally by this program
  !*********************************************************************

  integer :: j, iunit, i
  real, dimension(nalpha, nlambda) :: lambda

  !*********************************************************************
  ! Beginning of executable statements
  !*********************************************************************

  call get_command_argument(1,vmec_filename)
  call get_command_argument(2,gs2_output)
  call get_command_argument(3,geometry_output)
  call get_command_argument(4,number_of_field_periods_to_includestr)
  read(number_of_field_periods_to_includestr , *) number_of_field_periods_to_include
  call get_command_argument(5,desired_normalized_toroidal_fluxstr)
  read(desired_normalized_toroidal_fluxstr , *) desired_normalized_toroidal_flux

  call vmec_to_gs2_geometry_interface(vmec_filename, nalpha, nzgrid, zeta_center, number_of_field_periods_to_include, &
       desired_normalized_toroidal_flux, vmec_surface_option, verbose, &
       normalized_toroidal_flux_used, safety_factor_q, shat, L_reference, B_reference, &
       alpha, zeta, bmag, gradpar, gds2, gds21, gds22, gbdrift, gbdrift0, cvdrift, cvdrift0)

  print *,"-------------- Input parameters ------------------"
  print *,"vmec_filename: ",trim(vmec_filename)
  print *,"nalpha:",nalpha
  print *,"nzgrid:",nzgrid
  print *,"zeta_center:",zeta_center
  print *,"number_of_field_periods_to_include:",number_of_field_periods_to_include
  print *,"desired_normalized_toroidal_flux:",desired_normalized_toroidal_flux
  print *,"vmec_surface_option:",vmec_surface_option

  print *,"-------------- Output parameters -----------------"
  print *,"normalized_toroidal_flux_used:",normalized_toroidal_flux_used
  print *,"safety_factor_q:",safety_factor_q
  print *,"shat:",shat
  print *,"L_reference:",L_reference
  print *,"B_reference:",B_reference
  print *,"alpha:"
  print *,alpha
  !print *,"zeta:"
  !print *,zeta

  ! print *,"bmag:"
  ! do j=1,nalpha
  !    print *,bmag(j,:)
  ! end do

  ! print *,"gradpar:"
  ! do j=1,nalpha
  !    print *,gradpar(j,:)
  ! end do

  ! print *,"gds2:"
  ! do j=1,nalpha
  !    print *,gds2(j,:)
  ! end do

  ! print *,"gds21:"
  ! do j=1,nalpha
  !    print *,gds21(j,:)
  ! end do

  ! print *,"gds22:"
  ! do j=1,nalpha
  !    print *,gds22(j,:)
  ! end do

  ! print *,"gbdrift:"
  ! do j=1,nalpha
  !    print *,gbdrift(j,:)
  ! end do

  ! print *,"gbdrift0:"
  ! do j=1,nalpha
  !    print *,gbdrift0(j,:)
  ! end do

  ! print *,"cvdrfit:"
  ! do j=1,nalpha
  !    print *,cvdrift(j,:)
  ! end do

  ! print *,"cvdrift0:"
  ! do j=1,nalpha
  !    print *,cvdrift0(j,:)
  ! end do

  ! print *,"lambda:"
  do j=1,nalpha
   bMax=maxval(bmag(j,:))
   bMin=minval(bmag(j,:))
   do i=1,nlambda
      lambda(j,i)=(bMax - bMax*i + bMin*i - bMin*nlambda)/(bMax*bMin - bMax*bMin*nlambda)
   end do
  ! print *,lambda(j,:)
  end do

  iunit = 6
  open(file=gs2_output,unit=iunit)
  !write (iunit,*) 'nalpha nzgrid'
  !write (iunit,*) nalpha, nzgrid
  !write (iunit,*) 'alpha'
  !write (iunit,*) alpha
  !write (iunit,*) 'zeta'
  !write (iunit,*) zeta

  write (iunit,*) 'nlambda'
  write (iunit,*) nlambda

  write (iunit,*) 'lambda'
  do j=1,nalpha
     write (iunit,'(1(F12.8,:,1X))') lambda(j,:)
  end do

  write (iunit,*) 'ntgrid nperiod ntheta drhodpsi rmaj shat kxfac q'
  write (iunit,*) nzgrid, 1, 2*nzgrid+1, 1., 1., shat, 1, safety_factor_q

  write (iunit,*) 'gbdrift gradpar grho tgrid'
  do j=1,nalpha
      do i=1,2*nzgrid+1
         write (iunit,*) gbdrift(j,i-nzgrid-1), gradpar(j,i-nzgrid-1), 1., zeta(i-nzgrid-1)
      end do
  end do

  write (iunit,*) 'cvdrift gds2 bmag tgrid'
  do j=1,nalpha
      do i=1,2*nzgrid+1
         write (iunit,*) cvdrift(j,i-nzgrid-1), gds2(j,i-nzgrid-1), bmag(j,i-nzgrid-1), zeta(i-nzgrid-1)
      end do
  end do

  write (iunit,*) 'gds21 gds22 tgrid'
  do j=1,nalpha
      do i=1,2*nzgrid+1
         write (iunit,*) gds21(j,i-nzgrid-1), gds22(j,i-nzgrid-1), zeta(i-nzgrid-1)
      end do
  end do

  write (iunit,*) 'cvdrift0 gbdrift0 tgrid'
  do j=1,nalpha
      do i=1,2*nzgrid+1
         write (iunit,*) cvdrift0(j,i-nzgrid-1), gbdrift0(j,i-nzgrid-1), zeta(i-nzgrid-1)
      end do
  end do

  write (iunit,*) 'Rplot Rprime tgrid'
  do j=1,nalpha
      do i=1,2*nzgrid+1
         write (iunit,*) 0., 0., zeta(i-nzgrid-1)
      end do
  end do

  write (iunit,*) 'Zplot Zprime tgrid'
  do j=1,nalpha
      do i=1,2*nzgrid+1
         write (iunit,*) 0., 0., zeta(i-nzgrid-1)
      end do
  end do

  write (iunit,*) 'aplot aprime tgrid'
  do j=1,nalpha
      do i=1,2*nzgrid+1
         write (iunit,*) 0., 0., zeta(i-nzgrid-1)
      end do
  end do

  close(iunit)

  iunit = 6
  open(file=geometry_output,unit=iunit)

  write (iunit,*) 'tgrid(',j,')=',zeta
  do j=1,nalpha
   write (iunit,*) 'gbdrift(',j,')=',gbdrift
   write (iunit,*) 'gradpar(',j,')=',gradpar
   write (iunit,*) 'cvdrift(',j,')=',cvdrift
   write (iunit,*) 'gds2(',j,')=',gds2
   write (iunit,*) 'bmag(',j,')=',bmag
   write (iunit,*) 'gds21(',j,')=',gds21
   write (iunit,*) 'gds22(',j,')=',gds22
   write (iunit,*) 'cvdrift0(',j,')=',cvdrift0
   write (iunit,*) 'gbdrift0(',j,')=',gbdrift0
   write (iunit,*) 'lambda(',j,')=',lambda
  end do

  write (iunit,*) 'shat=',shat
  write (iunit,*) 'alpha=',alpha
  write (iunit,*) 'nlambda=',nlambda
  write (iunit,*) 'normalized_flux=',normalized_toroidal_flux_used

  close(iunit)

end program test_vmec_to_gs2_geometry_interface
