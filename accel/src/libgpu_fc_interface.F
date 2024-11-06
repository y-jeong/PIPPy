      module libgpu_fc_interface
      use iso_c_binding
      implicit none
      
      interface
         type (c_ptr) function libgpu_fc_create_device()
     *         bind(C, name="libgpu_create_device")
         use iso_c_binding
         implicit none
         end function libgpu_fc_create_device
      
         subroutine libgpu_fc_destroy_device(ptr)
     *     bind(C,name="libgpu_destroy_device")
         use iso_c_binding
         implicit none
         type (c_ptr), value :: ptr
         end subroutine libgpu_fc_destroy_device

         integer (c_int) function libgpu_fc_get_num_devices(ptr)
     *     bind(C,name="libgpu_get_num_devices")
         use iso_c_binding
         implicit none
         type (c_ptr), value :: ptr
         end function libgpu_fc_get_num_devices

         subroutine libgpu_fc_set_device(ptr, id)
     *     bind(C,name="libgpu_set_device")
         use iso_c_binding
         implicit none
         type (c_ptr), value :: ptr
         integer (c_int), value :: id
         end subroutine libgpu_fc_set_device

         subroutine libgpu_fc_setup(ptr, ind_1d, coef, ibasis,
     *     nterms, ncoefs, npairs, nsurfs) bind(C,name="libgpu_setup")
         use iso_c_binding
         implicit none
         
         type (c_ptr), value :: ptr
         integer (c_int), value :: nterms, ncoefs, npairs, nsurfs
         integer (c_int), dimension(*) :: ind_1d
         integer (c_int), dimension(*) :: ibasis
         real (c_double), dimension(*) :: coef
         end subroutine libgpu_fc_setup
     
         subroutine libgpu_fc_compute(ptr, r, v, dvdr, isurf)
     *     bind(C, name="libgpu_compute")
         use iso_c_binding
         implicit none
         
         type (c_ptr), value :: ptr
         real (c_double), dimension(*) :: r
         type (c_ptr), value :: v
         real (c_double), dimension(*) :: dvdr
         integer (c_int), value :: isurf
         end subroutine libgpu_fc_compute

         subroutine libgpu_fc_compute_v(ptr, r, v, isurf)
     *     bind(C, name="libgpu_compute_v")
         use iso_c_binding
         implicit none
     
         type (c_ptr), value :: ptr
         real (c_double), dimension(*) :: r
         type (c_ptr), value :: v
         integer (c_int), value :: isurf
         end subroutine libgpu_fc_compute_v

      end interface
      end module libgpu_fc_interface
