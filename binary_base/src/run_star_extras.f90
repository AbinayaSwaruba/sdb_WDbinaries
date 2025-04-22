! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras 

      use star_lib
      use star_def
      use auto_diff
      use const_def
      use math_lib
      use binary_def
      use utils_lib, only: mesa_error
      
      implicit none
      
    contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% other_energy => Heat_WD
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         
      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup


      integer function extras_start_step(id)
         use chem_def
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         integer :: k, nz, h1,he4
         real(dp) :: eps_nuc,h1_xa
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
         
      end function extras_start_step
      
      subroutine Heat_WD(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            integer :: k, kmid, nz, mid_k
            real(dp) :: L_extra, dM, Mr,star_mass, total_mass, cell_mass,mass_mid 
            real(dp), parameter :: tolerance = 1.0e-2_dp  
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            
            nz = s% nz
            star_mass = 0
            mass_mid = s%star_mass/2
            print *, "mass_mid", mass_mid
            do k = 1, nz
               
               cell_mass = s% dm(k) / msun
               star_mass = star_mass + cell_mass
               if (star_mass - mass_mid > 0) then
                  kmid = k
                  exit 
               end if
            end do

            print *, "total_mass, star_mass", star_mass, s%star_mass/2
            print *, "mass_co-ordinate, kmid", kmid

            L_extra = 10**(-1.9)*Lsun
            
            do k = kmid , s% nz
                  s% extra_heat(k) = L_extra / (mass_mid*Msun)
            end do

            ! output rate at which energy is added
            write(*,*) "Added"
            !*) "Added ", dot_product(s% extra_heat(1:s%nz), s% dm(1:s%nz))/Lsun, &
                  ! x_ctrl(2)

      end subroutine Heat_WD


      integer function extras_check_model(id)
         integer, intent(in) :: id
         extras_check_model = keep_going
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: dt
         ierr = 0
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0
      end subroutine data_for_extra_profile_columns
      

      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      

      end module run_star_extras
      
