program cubicSpline

        real (kind=8), dimension(1524) :: xpt, ypt
        real (kind=8) :: s
        integer :: cntr
       
        ! Get our data from a file
        call getData(xpt, ypt, cntr)
       
        ! Open Output file
        open(110, file = 'NCSData.txt')

        ! Interpolate for points at 1-110 in increments of 1 and output x and y
        ! values int NCSData.txt
        do i= 10, 110, 1
               
                call NCS(cntr, xpt, ypt, real(i, 8), s)
                write(110, *) i, s
        enddo
  
      
end program cubicSpline

subroutine getData(xpt, ypt, cntr)

        implicit none
        
        real (kind=8), dimension(1524) :: xpt, ypt
        real (kind=8) :: x, y
        integer :: errval, cntr
        
        open(unit=15, file="data.txt", status="old")

        errval = 0
        cntr = 0

        do while (errval /= -1)

                read(15, *, iostat=errval) x, y
                
                if(errval == 0) then
                        cntr = cntr + 1
                        xpt(cntr) = x
                        ypt(cntr) = y
                else
                        close(15)
                endif
        enddo

end subroutine getData


subroutine NCS(cntr, x, y, xNew, s)
        
        ! Parameter for LAPACK
        parameter (MAXDIM=20)
        ! Variables for LAPACK
        real (kind=8), dimension(MAXDIM) :: A, ASUB, ASUP, ASUP2, B        
        integer, dimension(MAXDIM) :: IPIV
        integer :: INFO, N

        ! Variables for algorithm
        integer :: cntr, hNum, bNum
        real (kind=8), dimension(1524) :: x, y, h, mB, bV, dV
        real (kind=8) :: xNew, s      

        N = MAXDIM        
        
        ! Get h values to make the b matrix values
        do i=1, cntr-1, 1
                h(i) = x(i+1) - x(i) 
        enddo
        
 
        !  Get b matrix values
        do i=1, cntr-1, 1
                mB(i) = ((3/h(i+1))*(y(i+2)-y(i+1)))-((3/h(i))*(y(i+1)-y(i)))
        enddo


        ! Fill up LAPACK Matrix and call Tridiagonal solver
        
        ! Diagonal of A matrix
        hNum = 1
        A(1) = 1
        do j=2, cntr-1, 1
                A(j) = 2*(h(hNum)+h(hNum+1))
                hNum = hNum + 1
        enddo
        A(cntr) = 1
                        
        ! Sup-Diagonal of A matrix
        hNum = 2
        ASUP(1) = 0
        do j=2, cntr-1, 1
                ASUP(j) = h(hNum)
                hNum = hNum + 1
        enddo
        ASUP(cntr) = 0        
        
        ! Sub-Diagonal of A matrix
        hNum = 1
        do j=1, cntr-2, 1
                ASUB(j) = h(hNum)
                hNum = hNum + 1
        enddo
        ASUB(cntr-1) = 0
        ASUB(cntr) = 0

        ! B matrix
        bNum = 1
        B(1) = 0 
        do j=2, cntr-1, 1
                B(j) = mB(bNum)
                bNum = bNum + 1
        enddo
        B(cntr) = 0
        
        ! Call tridiagonal solver functions
        CALL DGTTRF( N, ASUB, A, ASUP, ASUP2, IPIV, INFO)
        CALL DGTTRS( 'N', N, 1, ASUB, A, ASUP, ASUP2, IPIV, B, N, INFO)
         
        ! Build bV and dV from B() we got from LAPACK as cV and the y is aV
        
        do i=1, cntr-1

                bV(i) = (1/h(i))*(y(i+1)-y(i))-(h(i)/3)*(2*B(i)+B(i+1))
                dV(i) = (1/(3*h(i)))*(B(i+1)-B(i))
        enddo
                  
        do i=1, cntr-1
        
                print *, ' '
                print *, y(i), bV(i), B(i), dV(i)
        enddo

        ! Find the right piecewise s function and evaluate for our new x                        
        do i=1, cntr-1
                
                
                if( x(i) < xNew .AND. xNew < x(i+1) ) then
                                       
                        s = y(i) + bV(i)*(xNew-x(i)) + B(i)*(xNew-x(i))*(xNew-x(i)) + dV(i)*&
                            (xNew-x(i))*(xNew-x(i))*(xNew-x(i))                                
                        
                endif
        enddo
                                   
end subroutine NCS
