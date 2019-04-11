      program dvr_2d
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      parameter (npt = ns*na)
      common/ke/pi,us,ua
      dimension ham(npt,npt),eng(npt),tmp(npt)
      dimension dip(npt)
      dimension ipos(2,npt)
      external get_pot
      open (unit=20,file='pot_cut.dat',status='unknown')
      open (unit=21,file='wf1.dat',status='unknown')
      open (unit=22,file='wf2.dat',status='unknown')
      open (unit=23,file='wf3.dat',status='unknown')
      call setup_pot
      
c		pulled from setup_pot.f, we have
c
c
c		convm = 6.0221367e23*9.1093897e-28
c														! the second part is the electron mass in g; the first part is NA; so this takes a molar mass to me
c		uh = 1.007825/convm				! this then converts 1.007825 g/mol to electron mass
c 		us = .5*(.5/uh)							! then for the symmetric stretch each of those .5 will be converted to a 2 when we do 1/us
c		ua = .5*(.5/uh+2./uh)				! for the asymmetric stretch we use the h+h2 unit
c		
       
      s = 2.81										! initial s value?
      smin = .4*sqrt(2.)/.529177		! min s grid dim
      smax = 3.*sqrt(2.)/.529177		! max s grid dim
      sinc = (smax-smin)/float(ns-1)	! mesh spacing in s
      scoef = us/sinc**2						! m*dx^2 in s

      amin = -1.75*sqrt(2.)/.529177	! min a grid dim
      amax = 1.75*sqrt(2.)/.529177	! max a grid dim
      ainc = (amax-amin)/float(na-1)	! mesh spacing in a
      acoef = ua/ainc**2					! m*dx^2 in a

      print*,'tcoef',scoef,pi
      ip = 0.
      a = amin-ainc							! initialize a to one knotch below its min value
      do i = 1,na									! loop over all the a
         a = a + ainc								! increment a
         s = smin-sinc							! initialize s to one knotch below its min value
         do j = 1, ns								! loop over the s
            s = s + sinc							! increment s
            v = get_pot(s, a, d)
            write(20,*)s,a,v*219474.6
            ip = ip + 1							! increment gridpoint counter
            dip(ip) = d							! no idea what d is...
            											! diagonal element
            ham(ip,ip) = v + (scoef+acoef)*pi*pi/3.
            ipos(1,ip) = i						! add in a coordinate for matrix holding indices for the grid points
            ipos(2,ip) = j						! add in s coordinate for matrix holding indices for the grid points
            do k = 1,ip-1						! loop over previous grid points
               k1 = ipos(1,k)					! get a index value for grid point k
               k2 = ipos(2,k)					! get s index value for grid point k
               if (k1.eq.i) then					! diagonal part in a
               											! off-diagonal parts of classic CM ham in s
                  ham(k,ip) = scoef*2./float(k2-j)**2*(-1)**(k2-j)
               elseif(k2.eq.j) then			! diagonal part in s
               											! off-diagonal parts of classic CM ham in a
                  ham(k,ip) = acoef*2./float(k1-i)**2*(-1)**(k1-i)
               endif
               ham(ip,k) = ham(k, ip)		! enforce symmetry
            enddo
         enddo
      enddo
      													! diagonalize hamiltonian?
      call house(ham,npt,npt,eng,tmp)
      e0 = eng(1)*219474.6				! ground-state energy
      do i = 1,30									! pull out 30 wavefunctions?
         d = 0.										
         do j = 1,npt
            d = d + dip(j)*ham(j,i)*ham(j,1)
         enddo
         print*,i,eng(i)*219474.6,eng(i)*219474.6-e0,d**2
      enddo
      do i = 1,npt
         a = amin+(ipos(1,i)-1)*ainc
         s = smin+(ipos(2,i)-1)*sinc
         write(21,90)a,s,(ham(i,j),j=1,10)
         write(22,90)a,s,(ham(i,j),j=11,20)
         write(23,90)a,s,(ham(i,j),j=21,30)
      enddo
  90  format(2f10.6,10f12.8)
      stop
      end
         
