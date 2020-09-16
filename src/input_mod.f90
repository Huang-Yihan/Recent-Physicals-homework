!Author: Huang Yihan
!This module is used to read some input parameters
!include: Radius, length of side, dancoff(This could be calculated)
!         Ressonance table file name
!The 69 group of energy structure of WIMS is given in this module
module input_mod
    implicit none
    real*8,parameter::WIMS_boundary(70)=(/&
    1.00E+07,6.0655e6,3.679e6,2.231e6,1.353e6,&
    8.21E+05,5.e5,3.025e5,1.83e5,1.11e5,6.734e4, &
    4.09E+04,2.478e4,1.503e4,9118.e0,5530.e0,&
    3.52E+03,2239.45e0,1425.1e0,906.898e0,367.262e0,&
    1.49E+02,75.501e0,48.052e0,27.7e0,15.968e0,9.88E+00,&
    4.e0,3.3e0,2.6e0,2.1e0,1.5e0,1.3e0,1.15E+00,1.123e0,&
    1.097e0,1.071e0,1.045e0,1.02e0,9.96E-01,.972e0,.95e0,&
    .91e0,.85e0,.78e0,6.25E-01,.5e0,.4e0,.35e0,.32e0,.3e0,&
    .28e0,2.50E-01,.22e0,.18e0,.14e0,.1e0,.08e0,6.70E-02,&
    .058e0,.05e0,.042e0,.035e0,.03e0,2.50E-02,.02e0,.015e0,&
     .01e0,.005e0,1.e-5/)    !the energy structure of WIMS 69 groups
    character*5,parameter    ::   input='input'!the filename of inputcard
    character*30             ::   libname      !the lib name
    real*8                   ::   radius       !the radius of pin
    real*8                   ::   length       !the length of side
    real*8                   ::   dancoff(2)   !Read dancoff
    real*8                   ::   density_H1   !the density of H1
    real*8                   ::   density_U8   !the density of U8
    integer*4,parameter      ::   N1=14        !the number of fast groups
    integer*4,parameter      ::   N2=13        !the number of res groups
    integer*4,parameter      ::   first_tag=15 !the first resonance group tag
    contains
        !===========================================
        !<Read inputs
        subroutine read_inputs
            call read_lib(libname)
            call read_radius(radius)
            call read_dancoff(dancoff(1),dancoff(2))
            call read_densityU8(density_U8)
            call read_densityU8(density_H1)
        end subroutine read_inputs
        !===========================================
        !<Read Filename subroutine
        subroutine read_lib(libname)
            character*30,intent(inout)   :: libname
            character*80                 :: words
            integer*4                    :: k
            open(21,file=input)
            do
                read(21,'(A80)') words
                k=index(words,'Lib:')
                if (k .gt. 0) then
                    read(words(k+4:80),*) libname
                    exit
                end if
            end do
            close(21)
        end subroutine read_lib
        !>===========================================
        !<Read Radius subroutine 
        subroutine read_radius(radius)
            real*8,intent(inout)         :: radius
            character*80                 :: words
            integer*4                    :: k
            open(21,file=input)
            do
                read(21,'(A80)') words
                k=index(words,'Radius:')
                if (k .gt. 0) then
                    read(words(k+7:80),*) radius
                    exit
                end if  
            end do      
            close(21)
        end subroutine read_radius
        !>============================================
        !<Read Length of side subroutine
        subroutine read_length(length)
            real*8,intent(inout)         :: length
            character*80                 :: words
            integer*4                    :: k
            open(21,file=input)
            do
                read(21,'(A80)') words
                k=index(words,'Length:')
                if (k .gt. 0) then
                    read(words(k+7:80),*) length
                    exit
                end if   
            end do     
            close(21)
        end subroutine read_length
        !>============================================
        !<Read Temperature subroutine
        subroutine read_temperature(temperature)
            real*8,intent(inout)         :: temperature
            character*80                 :: words
            integer*4                    :: k
            open(21,file=input)
            do
                read(21,'(A80)') words
                k=index(words,'Temperature:')
                if (k .gt. 0) then
                    read(words(k+13:80),*) temperature
                    exit
                end if   
            end do     
            close(21)
        end subroutine read_temperature
        !>============================================
        !<Read Dancoff subroutine
        subroutine read_Dancoff(dancoff1,dancoff2)
            real*8,intent(inout)         :: dancoff1
            real*8,intent(inout)         :: dancoff2
            character*80                 :: words
            integer*4                    :: k
            open(21,file=input)
            do
                read(21,'(A80)') words
                k=index(words,'Dancoff:')
                if (k .gt. 0) then
                    read(words(k+8:80),*) dancoff1,dancoff2
                    exit
                end if   
            end do     
            close(21)
        end subroutine read_Dancoff
        !>============================================
        !<Read density_U8 subroutine
        subroutine read_densityU8(density_U8)
            real*8,intent(inout)         :: density_U8
            character*80                 :: words
            integer*4                    :: k
            open(21,file=input)
            do
                read(21,'(A80)') words
                k=index(words,'Density_U8:')
                if (k .gt. 0) then
                    read(words(k+11:80),*) density_U8
                    exit
                end if   
            end do     
            close(21)
        end subroutine read_densityU8
        !>============================================
        !<Read density_H1 subroutine
        subroutine read_densityH1(density_H1)
            real*8,intent(inout)         :: density_H1
            character*80                 :: words
            integer*4                    :: k
            open(21,file=input)
            do
                read(21,'(A80)') words
                k=index(words,'Density_U8:')
                if (k .gt. 0) then
                    read(words(k+11:80),*) density_H1
                    exit
                end if   
            end do     
            close(21)
        end subroutine read_densityH1
end module input_mod