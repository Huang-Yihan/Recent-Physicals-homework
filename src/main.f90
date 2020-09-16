!Author Huang Yihan
!This is the main program of this homework
!1.read inputs
!2.read resonance table
!3.calculate
Program main
    use input_mod
    use table_mod
    use calculate_mod
    implicit none
    call read_inputs
    call res_table_init
    call U8_xs_init
    call H1_xs_init
    call calculate_sigmae(radius,sigmae)
    call calculate_xs
    call deallocate_array
end program main