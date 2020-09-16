!Author:Huang Yihan
!This module is used to read Resonance score table with WIMS-69 structure
!You should give the filename and the first group
!If the first group didn't given,the first group is 4.00eV
module table_mod
    use input_mod
    implicit none
    type res_table_para_type
        integer*4               :: energy_tag
        real*8                  :: low_boundary
        real*8                  :: high_boundary
        real*8                  :: RIN          !RIN
        integer*4               :: M1           !M1
        integer*4               :: M2           !M2
        real*8,allocatable      :: T(:)         !T
        real*8,allocatable      :: SIGP(:)      !SIGP
        real*8,allocatable      :: RSIG(:,:)    !RSIG
        real*8,allocatable      :: U8XSP
        real*8,allocatable      :: H1XSP
        real*8,allocatable      :: U8XSA
        real*8,allocatable      :: U8lamda
        real*8,allocatable      :: H1lamda
    end type
    type(res_table_para_type)   :: group(N2)
    contains
    !==============================================================
    !table_init
    !this subroutine is used to read resonance table and init it 
    !==============================================================
        subroutine res_table_init
            integer*4           :: i
            integer*4           :: j
            integer*4           :: r
            integer*4           :: c
            real*8,allocatable  :: real_array(:) !the arrary to read data from lib file
            open(31,file=libname,action='read')
            do i=1,12
                read(31,*)
            end do
            do i= 1 , N2
                !>init tag and boundary
                group(i)%energy_tag=N1+i
                group(i)%high_boundary=WIMS_boundary(group(i)%energy_tag)
                group(i)%low_boundary=WIMS_boundary(group(i)%energy_tag+1)
                !> read M1,M2
                read(31,*) group(i)%RIN,group(i)%M1,group(i)%M2
                !>allocate T SIGP and RSIG
                allocate(group(i)%T(group(i)%M1))
                allocate(group(i)%SIGP(group(i)%M2))
                allocate(group(i)%RSIG(group(i)%M1,group(i)%M2))
                group(i)%T=0
                group(i)%SIGP=0
                group(i)%RSIG=0
                allocate(real_array(group(i)%M1+group(i)%M2+group(i)%M1*group(i)%M2))
                real_array=0
                call read_real_array(31,group(i)%M1+group(i)%M2+group(i)%M1*group(i)%M2,&
                    real_array)
                do j = 1, group(i)%M1+group(i)%M2+group(i)%M1*group(i)%M2
                    if(j .le. group(i)%M1 ) then
                        group(i)%T(j)=real_array(j)
                    else if (j .gt. group(i)%M1 .and. j .le. group(i)%M1+group(i)%M2) then
                        group(i)%SIGP(j-group(i)%M1)=real_array(j)
                    else if (j .gt. group(i)%M1+group(i)%M2 .and. j .le. &
                    group(i)%M1+group(i)%M2+group(i)%M1*group(i)%M2) then
                        r=(j-group(i)%M1-group(i)%M2)/group(i)%M2+1
                        c=j-group(i)%M1-group(i)%M2-(r-1)*group(i)%M2
                        if (c .eq. 0 ) then
                            c=c+group(i)%M2
                            r=r-1
                        end if
                        group(i)%RSIG(r,c)=real_array(j)
                    end if
                end do
                deallocate(real_array)
            end do
            close(31)
        end subroutine res_table_init
        !==================================================
        !read_real_array
        !this subroutine is used to read a long array in 
        !resonance table
        !==================================================

        subroutine read_real_array(nin,number,a)
            integer*4,intent(in)                       :: nin     !the libfile tag
            integer*4,intent(in)                       :: number  !the number of a
            real*8,intent(inout),allocatable           :: a(:)    !the array
            integer*4                                  :: i
            if (allocated(a) .eqv. .false.) then 
                allocate(a(number))
                a=0 
            end if
            do i=1,number
                if ( mod(i,5) .ne. 0 .and. i .ne. number) then
                    read(nin,'(x,Es14.8)',advance='no') a(i)
                else
                    read(nin,'(x,Es14.8)',advance='yes') a(i)
                end if
            end do
        end subroutine read_real_array
        !==================================================
        !read U8 sigmap and init
        !==================================================
        subroutine U8_xs_init
            character*30          :: words
            integer*4             :: i
            real*8,allocatable    :: real_array(:) !the arrary to read data from lib file
            integer*4             :: J
            real*8                :: AW
            integer*4             :: IAN
            integer*4             :: NF
            integer*4             :: NT
            integer*4             :: NZZ
            open(31,file=libname,action='read')
            rewind(31)
            do i=1,12
                read(31,*)
            end do
            do
                read(31,*) words
                if (words == '999999999') then
                    read(31,*) J,AW,IAN,NF,NT,NZZ
                    if (J .eq. 8238 .and. IAN .eq. 92) exit
                end if
            end do
            if (J .eq. 8238 .and. IAN .eq. 92) then
                allocate(real_array(2*N1+6*N2))
                call read_real_array(31,2*N1+6*N2,real_array)
                do i=1,N2
                    group(i)%U8XSP=real_array(i)
                end do
                do i=1,N2
                    group(i)%U8lamda=real_array(i+2*N1+5*N2)
                    if (group(i)%U8lamda .lt. 1.000000) then
                        group(i)%U8XSP=group(i)%U8XSP/group(i)%U8lamda
                    end if
                end do
                do i=1,N2
                    group(i)%U8XSA=real_array(i+2*N1+3*N2)
                end do
                deallocate(real_array)
            end if
            close(31)
        end subroutine U8_xs_init
        !==================================================
        !read H1 sigmap and init
        !==================================================
        subroutine H1_xs_init
            character*30          :: words
            integer*4             :: i
            real*8,allocatable    :: real_array(:) !the arrary to read data from lib file
            integer*4             :: J
            real*8                :: AW
            integer*4             :: IAN
            integer*4             :: NF
            integer*4             :: NT
            integer*4             :: NZZ
            open(31,file=libname,action='read')
            rewind(31)
            do i=1,12
                read(31,*)
            end do
            do
                read(31,*) words
                if (words == '999999999') then
                    read(31,*) J,AW,IAN,NF,NT,NZZ
                    if (J .eq. 3001 .and. IAN .eq. 1) exit
                end if
            end do
            if (J .eq. 3001 .and. IAN .eq. 1) then
                allocate(real_array(2*N1+6*N2))
                call read_real_array(31,2*N1+6*N2,real_array)
                do i=1,N2
                    group(i)%H1XSP=real_array(i)
                end do
                do i=1,N2
                    group(i)%H1lamda=real_array(i+2*N1+5*N2)
                    if (group(i)%H1lamda .lt. 1.000000) then
                        group(i)%H1XSP=group(i)%H1XSP/group(i)%H1lamda
                    end if
                end do
                deallocate(real_array)
            end if
            close(31)
        end subroutine H1_xs_init
        !==================================================
        !interpolation
        !this subroutine is a interpolation function to get 
        !the cross sections
        !==================================================

        subroutine interpole(energy_tag,energy_T,sigmap,interpolation)
            integer*4,intent(in) :: energy_tag    ! the energy group tag
            real*4,intent(in)    :: energy_T      ! T
            real*8,intent(in)    :: sigmap        ! sigp
            real*8,intent(inout) :: interpolation ! interpolation
            integer*4            :: tag
            integer*4            :: T
            integer*4            :: i
            !> find energy tag
            do i=1,N2
                if(group(i)%energy_tag .eq. energy_tag) tag=i
            end do
            do i=1,group(tag)%M1
                if(group(tag)%T(i)==energy_T) T=i
            end do
            do i=1,group(tag)%M2-1
                if(group(tag)%SIGP(i) .le. sigmap .and. &
                 group(tag)%SIGP(i+1) .ge. sigmap) then
                    interpolation=(sigmap-group(tag)%SIGP(i))/(group(tag)%SIGP(i+1)- &
                    group(tag)%SIGP(i))
                    interpolation=interpolation*(group(tag)%RSIG(T,i+1)-group(tag)%RSIG(T,i))+group(tag)%RSIG(T,i)
                end if
            end do
        end subroutine interpole
        !==================================================
        !deallocate_array
        !This is used to deallocate
        !==================================================
        subroutine deallocate_array
            integer*4               ::i
            do i=1,N2
                    deallocate(group(i)%T)         
                    deallocate(group(i)%SIGP)     
                    deallocate(group(i)%RSIG)    
                    deallocate(group(i)%U8XSP)
                    deallocate(group(i)%H1XSP)
                    deallocate(group(i)%U8XSA)
                    deallocate(group(i)%U8lamda)
                    deallocate(group(i)%H1lamda)
            end do
        end subroutine deallocate_array
end module table_mod
