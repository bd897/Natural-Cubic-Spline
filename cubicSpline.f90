program cubicSpline

        real (kind=8), dimension(1524) :: xpt, ypt
        real (kind=8) :: s
        integer :: cntr, nNum
        
        
        call getData(xpt, ypt, cntr)
        
        nNum = NINT(xpt(cntr))-NINT(xpt(1))

        open(nNum, file = 'NCSData.txt')

        do i=NINT(xpt(1)), NINT(xpt(cntr)), 1

                call NCS(cntr, xpt, ypt, real(i, 8), s)
                write(nNum, *) i, s
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
        hNum = 1
        bNum = 1

        ! Get h values to make the b matrix values
        do i=0, cntr-1, 1
                h(i) = x(i+1) - x(i) 
        enddo

        !  Get b matrix values
        do i=1, cntr-1, 1
                mB(i) = ((3/h(i))*(y(i+1)-y(i)))-((3/h(i-1))*(y(i)-a(i-1)))
        enddo


        ! Fill up LAPACK Matrix and call Tridiagonal solver
        
        ! Diagonal of A matrix
        A(1) = 1
        do j=2, cntr-1, 1
                A(j) = 2*(h(hNum)+h(hNum+1))
                hNum = hNum + 1
        enddo
        A(cntr) = 1
                        
        ! Sup-Diagonal of A matrix
        hNum = 1

        ASUP(1) = 0
        do j=2, cntr-1, 1
                ASUP(j) = h(hNUM)
                hNum = hNum + 1
        enddo
        
        
        ! Sub-Diagonal of A matrix
        hNum = 1

        do j=1, cntr-2, 1
                ASUB = h(hNum)
                hNum = hNum + 1
        enddo
        ASUB(cntr-1) = 0

        ! B matrix
        B(1) = 0 
        do j=2, cntr-1, 1
                B(j) = mB(bNum)
                bNum = bNum + 1
        enddo
        B(cntr) = 0

        ! Call functions
        CALL DGTTRF( N, ASUB, A, ASUP, ASUP2, IPIV, INFO)
        CALL DGTTRS( 'N', N, 1, ASUB, A, ASUP, ASUP2, IPIV, B, N, INFO)
         
        ! Build bV and dV from B() we got from LAPACK as cV and the y is aV
        
                do i=1, cntr-1

                        bV(i) = (1/h(i))*(y(i+1)-y(i))-(h(i)/3)*(2*B(i)+B(i+1))
                        dV(i) = (1/(3*h(i)))*(B(i+1)-B(i))
                enddo
                
                do i=1, cntr-1
                        
                        if ( x(i) < xNew .AND. xNew < x(i+1) ) then
                                s = y(i) + bV(i)*xNew + B(i)*xNew*xNew + dV(i)* &
                                xNew*xNew*xNew
                                exit
                        else if( xNew < x(1) ) then
                                s = y(1) + bV(1)*xNew + B(1)*xNew*xNew + dV(1)* &
                                xNew*xNew*xNew
                                exit
                        else if( x(cntr-1) < xNew ) then
                                s = y(cntr-1) + bV(cntr-1)*xNew + B(cntr-1)*xNew*&
                                xNew + dV(cntr-1)*xNew*xNew*xNew
                                exit
                        endif
                enddo
                
        print *, s         
        ! Read out to file

end subroutine NCS












