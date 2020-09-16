!Author:Huang Yihan
!This module is used to calculate this problem
module calculate_mod
    use input_mod
    use table_mod
    implicit none
    real*8                  :: alpha(2)
    real*8                  :: beta(2)
    real*8                  :: sigmae
    integer*4               :: ii
    real*8                  :: sigma0(2)
    real*8                  :: ieff(2)
    real*8                  :: xseff(N2)
    integer*4               :: kk
    integer*4               :: jj
    contains
    !==========================================
    !<This subroutine is used to calculate sigmae
        subroutine calculate_sigmae(radius,sigmae)
            real*8,intent(in)     ::radius
            real*8,intent(inout)  ::sigmae
            sigmae=1/(2*radius)
        end subroutine calculate_sigmae
    !==========================================
    !<This subroutine is used to calculate alpha and beta,
    !<You should get a dancoff to subroutine
        subroutine calculate_ab(dancoff,alpha1,alpha2,beta1,beta2)
            real*8,intent(in)     :: dancoff
            real*8,intent(inout)  :: alpha1
            real*8,intent(inout)  :: alpha2
            real*8,intent(inout)  :: beta1
            real*8,intent(inout)  :: beta2
            real*8                :: A
            A=(1-dancoff)/dancoff
            alpha1=((5*A+6)-sqrt(A*A+36*A+36))/(2*(A+1))
            alpha2=((5*A+6)+sqrt(A*A+36*A+36))/(2*(A+1))
            beta1=((4*A+6)/(A+1)-alpha1)/(alpha2-alpha1)
            beta2=1-beta1
        end subroutine calculate_ab
    !===========================================
    !This subroutine is used to calculate sigma0
        subroutine calculate_sigma0(n,sigmae,density,sigma0,k)
            integer*4,intent(in)                           :: n
            real*8,intent(in)                              :: sigmae
            real*8,intent(in)                              :: density
            real*8,intent(inout)                           :: sigma0
            integer*4,intent(in)                           :: k
            sigma0=(group(n)%H1XSP*density_H1+alpha(k)*sigmae)/density
        end subroutine calculate_sigma0
    !===========================================
    !This subroutine is used to calculate Ieff
        subroutine calculate_i(i1,i2,i)
            real*8,intent(in)                              :: i1
            real*8,intent(in)                              :: i2
            real*8,intent(inout)                           :: i            
            i=beta(1)*i1+beta(2)*i2
        end subroutine calculate_i
    !===========================================
    !This subroutine is used to calculate Xseff
        subroutine calculate_xseff(i1,i2,sigma01,sigma02,xseff)
            real*8,intent(in)                              :: i1
            real*8,intent(in)                              :: i2
            real*8,intent(in)                              :: sigma01
            real*8,intent(in)                              :: sigma02
            real*8,intent(inout)                           :: xseff
            xseff=(beta(1)*i1+beta(2)*i2)/(1-(beta(1)*i1/(group(ii)%U8XSP+sigma01))&
            -(beta(2)*i2/(group(ii)%U8XSP+sigma02)))
        end subroutine calculate_xseff
    !===========================================
    !This subroutine is the main subroutine to calculate Xseff
        subroutine calculate_xs
            do jj=1,2
                call calculate_ab(dancoff(jj),alpha(1),alpha(2),beta(1),beta(2))
                do ii=1,N2
                    do kk=1,2
                        call calculate_sigma0(ii,sigmae,density_U8,sigma0(kk),kk)
                        call interpole(ii+N1,293.0,sigma0(kk),ieff(kk))
                    end do
                    call calculate_xseff(ieff(1),ieff(2),sigma0(1),sigma0(2),xseff(ii))
                end do
                call write_xseff
            end do
        end subroutine calculate_xs
    !============================================
    !This subroutine is used to write results
        subroutine write_xseff
            write(*,*) 'pin=',jj
            write(*,*) 'Group  High Boundary Low Boundary  Xs_absorb'
            do ii=1,N2
                write(*,'(2x,I2,4x,Es10.4,4x,Es10.4,4x,Es12.6)') ii+N1,group(ii)%high_boundary,group(ii)%low_boundary,xseff(ii)
            end do
        end subroutine write_xseff
end module calculate_mod