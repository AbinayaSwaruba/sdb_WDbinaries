! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
      use const_def
      use math_lib
      
      implicit none
      include "test_suite_extras_def.inc"

      
      contains

      include "test_suite_extras.inc"
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      integer function extras_check_model(id)
        use chem_def
         integer, intent(in) :: id
         real(dp) :: h_env_mass, h_upper_limit, h_lower_limit
         integer :: ierr
         integer :: k, nz, h1
         real(dp) :: total_mass_h, cell_mass,Nowind_time
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if(s% x_integer_ctrl(1) == 1) then

            nz = s% nz
            h1 = s% net_iso(ih1)
            total_mass_h=0
            do k=1,nz
               cell_mass = s% dm(k)/msun
               total_mass_h = total_mass_h + cell_mass*s% xa(h1, k)
            end do
            print *, "Star1base total_mass_h = ", total_mass_h
            print *, "lg_dt", s% star_age

            h_upper_limit = s% x_ctrl(1)
            h_lower_limit = (1d0 - s% x_ctrl(2))*h_upper_limit
            
            h_env_mass = s% star_mass - s% he_core_mass
            print *, "h_upper_limit", h_upper_limit
            print *, "Reimers scaling factor", s% Reimers_scaling_factor
            

            if(total_mass_h > h_upper_limit) then
               extras_check_model = keep_going
            else if(total_mass_h < h_lower_limit .and. s% Reimers_scaling_factor > 0d0) then
               write(*,*) "M_env/target =", total_mass_h/h_lower_limit
               write(*,*) "Envelope mass outside tolerance, retrying."
               extras_check_model = retry
            else ! Turn off winds and finish the run
               s% Reimers_scaling_factor = 0d0
               s% Blocker_scaling_factor = 0d0
               Nowind_time = s%star_age
               s%x_ctrl(3) = s%star_age +1e6
               print *, "N  check_modelNowind_limit", s%x_ctrl(3)
               print *, "check_modelNowind_limit", Nowind_time
            end if
            
         end if

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model

      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
         

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
       ! returns either keep_going, retry, or terminate.
      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         use chem_def
         integer, intent(in) :: id
         real(dp) :: h_env_mass, h_upper_limit, h_lower_limit
         integer :: ierr
         integer :: k, nz, h1
         real(dp) :: total_mass_h, cell_mass,Nowind_time,time_limit
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
            extras_finish_step = keep_going
         if(s% x_integer_ctrl(1) == 1) then
            print *, "N finishing check_modelNowind_limit", s%x_ctrl(3)
            if (s%star_age > s%x_ctrl(3)) then
               print *, "finishing with ctrl,age", s%x_ctrl(3),s%star_age
               extras_finish_step = terminate
               s% termination_code = t_xtra1
               termination_code_str(t_xtra1) = 'hydrogen envelope removed'
            end if
         end if
         
      end function extras_finish_step

      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve

      end module run_star_extras
      
