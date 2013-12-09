program POSCAR
implicit none
!Open Units

integer,PARAMETER		:: wp = selected_real_kind(12,70) !(12,70)
integer,PARAMETER               :: dwp = selected_real_kind(12,70)

integer				:: i,j
real(kind=wp)                   :: theta, psi
real(kind=wp)	                :: half, lattice_lenght !correct length
real(kind=wp), dimension(3,3)  :: rotmat
	type atomo
	real(kind=wp)   		:: x, y, z
	end type atomo

type(atomo) :: atomPb, atomC,  atomN,  atomH1,  atomH2,  atomH3,  atomH4,  atomH5,  atomH6,  atomBr1,  atomBr2,  atomBr3
type(atomo) :: atomPbf,atomCf, atomNf, atomH1f, atomH2f, atomH3f, atomH4f, atomH5f, atomH6f, atomBr1f, atomBr2f, atomBr3f
type(atomo) :: middlep, centercell

open(unit = 10, file = "input.dat", status = "old"    )
open(unit = 20, file = "POSCARc",    status = "unknown")
open(unit = 30, file = "POSCAR",   status = "unknown")

!START

!read(10,"(f4.3)")half
read(10,"(f5.4)")lattice_lenght
read(10,"(f5.2,f5.2)")theta,psi
write(*,*)theta,psi,lattice_lenght
!lattice_lenght=half+half
!write(*,*)half,lattice_lenght
!lattice_lenght=half*2
centercell = atomo((0.5d0*lattice_lenght),(0.5d0*lattice_lenght),(0.5d0*lattice_lenght))
!centercell = atomo(half,half,half)
!it creates the rotation matrix
!call matrice (theta,psi,rotmat)

atomPbf  = atomo( 0.000d0,	 0.00d0,	 0.00d0	)
atomBr1f = atomo( 0.500d0,	 0.00d0,	 0.00d0	)
atomBr2f = atomo( 0.000d0,	 0.50d0,	 0.00d0	)
atomBr3f = atomo( 0.000d0,	 0.00d0,	 0.50d0	)
!write(*,*)atomBr1f
! Transform Pb Br in cartesian coordinates

atomPb   = atomo( 0.000d0*lattice_lenght,	 0.00d0*lattice_lenght,	 0.00d0*lattice_lenght	)

atomBr1  = atomo( 0.50000000d0*lattice_lenght,	 0.00d0*lattice_lenght,	 0.00d0*lattice_lenght	)
atomBr2  = atomo( 0.000d0*lattice_lenght,	 0.50d0*lattice_lenght,	 0.00d0*lattice_lenght	)
atomBr3  = atomo( 0.000d0*lattice_lenght,	 0.00d0*lattice_lenght,	 0.50d0*lattice_lenght	)

!atomBr1  = atomo( 2d0*half,	0.00d0,		0.00d0 )
!atomBr2  = atomo( 0.000d0,	2d0*half,	0.00d0  )
!atomBr3  = atomo( 0.000d0,	0.00d0,		2d0*half  )


!write(*,*)atomBr1
! CH3NH3 coordinate in cartesian
atomC	 = atomo(-0.72790,	 0.00000,	 0.00000)
atomN	 = atomo( 0.72790,	 0.00000,	 0.00000)
atomH1	 = atomo(-1.11210,	 0.00000,	-1.04309)
atomH2	 = atomo(-1.10130,	 0.90298,	 0.52744)
atomH3	 = atomo(-1.10148,	-0.90283,	 0.52757)
atomH4	 = atomo( 1.06158,	 0.84491,	-0.51776)
atomH5	 = atomo( 1.06178,	-0.84969,	-0.50994)
atomH6	 = atomo( 1.06178,	 0.00052,	 0.99096)


!theta	= 0.0
!psi	= 0.0

call matrice (theta,psi,rotmat)

call rotation (atomC%x,atomC%y,atomC%z,rotmat)
call rotation (atomN%x,atomN%y,atomN%z,rotmat)
call rotation (atomH1%x,atomH1%y,atomH1%z,rotmat)
call rotation (atomH2%x,atomH2%y,atomH2%z,rotmat)
call rotation (atomH3%x,atomH3%y,atomH3%z,rotmat)
call rotation (atomH4%x,atomH4%y,atomH4%z,rotmat)
call rotation (atomH5%x,atomH5%y,atomH5%z,rotmat)
call rotation (atomH6%x,atomH6%y,atomH6%z,rotmat)

!Middle point of the molecule, than of the cell

middlep = atomo((atomC%x+atomN%x)/2,(atomC%x+atomN%x)/2,(atomC%x+atomN%x)/2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!STARTING OF TRANSLATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!atomC%x = atomC%x + centercell%x
!atomC%y = atomC%y + centercell%y
!atomC%z = atomC%z + centercell%z
atomC  = atomo( atomC%x  + centercell%x, atomC%y  + centercell%y, atomC%z  + centercell%z)
atomN  = atomo( atomN%x  + centercell%x, atomN%y  + centercell%y, atomN%z  + centercell%z)
atomH1 = atomo( atomH1%x + centercell%x, atomH1%y  + centercell%y, atomH1%z  + centercell%z)
atomH2 = atomo( atomH2%x + centercell%x, atomH2%y  + centercell%y, atomH2%z  + centercell%z)
atomH3 = atomo( atomH3%x + centercell%x, atomH3%y  + centercell%y, atomH3%z  + centercell%z)
atomH4 = atomo( atomH4%x + centercell%x, atomH4%y  + centercell%y, atomH4%z  + centercell%z)
atomH5 = atomo( atomH5%x + centercell%x, atomH5%y  + centercell%y, atomH5%z  + centercell%z)
atomH6 = atomo( atomH6%x + centercell%x, atomH6%y  + centercell%y, atomH6%z  + centercell%z)
!atomN%x = atomN%x + centercell%x
!atomN%y = atomN%y + centercell%y
!atomN%z = atomN%z + centercell%z

!atomH1%x = atomH1%x + centercell%x
!atomH1%y = atomH1%y + centercell%y
!atomH1%z = atomH1%z + centercell%z

!atomH2%x = atomH2%x + centercell%x
!atomH2%y = atomH2%y + centercell%y
!atomH2%z = atomH2%z + centercell%z

!atomH3%x = atomH3%x + centercell%x
!atomH3%y = atomH3%y + centercell%y
!atomH3%z = atomH3%z + centercell%z

!atomH4%x = atomH4%x + centercell%x
!atomH4%y = atomH4%y + centercell%y
!atomH4%z = atomH4%z + centercell%z

!atomH5%x = atomH5%x + centercell%x
!atomH5%y = atomH5%y + centercell%y
!atomH5%z = atomH5%z + centercell%z

!atomH6%x = atomH6%x + centercell%x
!atomH6%y = atomH6%y + centercell%y
!atomH6%z = atomH6%z + centercell%z



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!END OF TRANSLATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Transform CH3NH3 in fractional coordinates

atomCf	 = atomo( atomC%x / lattice_lenght,  atomC%y / lattice_lenght, atomC%z / lattice_lenght)
atomNf	 = atomo( atomN%x / lattice_lenght,  atomN%y / lattice_lenght, atomN%z / lattice_lenght)
atomH1f	 = atomo(atomH1%x / lattice_lenght, atomH1%y / lattice_lenght, atomH1%z / lattice_lenght)
atomH2f	 = atomo(atomH2%x / lattice_lenght, atomH2%y / lattice_lenght, atomH2%z / lattice_lenght)
atomH3f	 = atomo(atomH3%x / lattice_lenght, atomH3%y / lattice_lenght, atomH3%z / lattice_lenght)
atomH4f	 = atomo(atomH4%x / lattice_lenght, atomH4%y / lattice_lenght, atomH4%z / lattice_lenght)
atomH5f	 = atomo(atomH5%x / lattice_lenght, atomH5%y / lattice_lenght, atomH5%z / lattice_lenght)
atomH6f	 = atomo(atomH6%x / lattice_lenght, atomH6%y / lattice_lenght, atomH6%z / lattice_lenght)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WRITE OUTPUT POSCAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CARTESIAN POSCAR

!write(20,*)"#Generated with fortran program by F.Brivio"		!Comment
!write(20,*)"1.0" 							!Scaling factor
!write(20,"(f6.2,f6.2,f6.2)")lattice_lenght,0.00000000,0.00000000
!write(20,"(f6.2,f6.2,f6.2)")0.00000000,lattice_lenght,0.00000000
!write(20,"(f6.2,f6.2,f6.2)")0.00000000,0.00000000,lattice_lenght
!write(20,*)"   I    C   H   N   Mm"					!Atom List
!write(20,*)"   3   1   6   1   1"					!Number of atom per species
!write(20,*)"   Cartesian"						!Type of coordinates			
!write(20,"(f12.7,f12.7,f12.7)")atomBr1
!write(20,"(f12.7,f12.7,f12.7)")atomBr2
!write(20,"(f12.7,f12.7,f12.7)")atomBr3
!write(20,"(f12.7,f12.7,f12.7)")atomC
!write(20,"(f12.7,f12.7,f12.7)")atomH1
!write(20,"(f12.7,f12.7,f12.7)")atomH2
!write(20,"(f12.7,f12.7,f12.7)")atomH3
!write(20,"(f12.7,f12.7,f12.7)")atomH4
!write(20,"(f12.7,f12.7,f12.7)")atomH5
!write(20,"(f12.7,f12.7,f12.7)")atomH6
!write(20,"(f12.7,f12.7,f12.7)")atomN
!write(20,"(f12.7,f12.7,f12.7)")atomPb

!FRACTIONAL POSCAR

write(30,*)"#Generated with fortran program by F.Brivio"		!Comment
write(30,*)"1.0" 							!Scaling factor
write(30,"(f6.3,f6.3,f6.3)")lattice_lenght,0.00000000,0.00000000
write(30,"(f6.3,f6.3,f6.3)")0.00000000,lattice_lenght,0.00000000
write(30,"(f6.3,f6.3,f6.3)")0.00000000,0.00000000,lattice_lenght
write(30,*)"   I    C   H   N   Mm"					!Atom List
write(30,*)"   3   1   6   1   1"					!Number of atom per species
write(30,*)"   Direct"							!Type of coordinates			
write(30,"(f12.7,f12.7,f12.7)")atomBr1f
write(30,"(f12.7,f12.7,f12.7)")atomBr2f
write(30,"(f12.7,f12.7,f12.7)")atomBr3f
write(30,"(f12.7,f12.7,f12.7)")atomCf
write(30,"(f12.7,f12.7,f12.7)")atomH1f
write(30,"(f12.7,f12.7,f12.7)")atomH2f
write(30,"(f12.7,f12.7,f12.7)")atomH3f
write(30,"(f12.7,f12.7,f12.7)")atomH4f
write(30,"(f12.7,f12.7,f12.7)")atomH5f
write(30,"(f12.7,f12.7,f12.7)")atomH6f
write(30,"(f12.7,f12.7,f12.7)")atomNf
write(30,"(f12.7,f12.7,f12.7)")atomPbf

close(10)
close(20)
close(30)

end program POSCAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine matrice (theta,psi,rotmat)
implicit none
integer,PARAMETER           			:: wp = selected_real_kind(12,70)!(6,35) !(12,70)
integer						:: i, j
real(kind=wp),			intent(InOut)	:: theta, psi
real(kind=wp), dimension(3,3),	intent(Out  )	:: rotmat
real(kind=wp), parameter			:: phi = 0.0
real(kind=wp)					:: pi,sintheta,costheta,sinpsi,cospsi,sinphi,cosphi
do i=1,3
	do j=1,3
	rotmat(i,j)=0.0
	end do
end do
!I convert the degrees value in radians
pi	= 4.0_wp*atan(1.0_wp)
theta	= (theta*pi/180.0_wp)
!psi	= (phi*pi/180.0_wp)
psi     = (psi*pi/180.0_wp)

!Instead of inserting the whole expresion i calculate first sin and cosine of angles, than I used them to build the matrix for the rot
sintheta=sin(theta)
costheta=cos(theta)
sinpsi=sin(psi)
cospsi=cos(psi)
sinphi=sin(phi)
cosphi=cos(phi)
!FIRST ROW (xx),(xy),(xz)
rotmat(1,1)=costheta*cospsi
rotmat(1,2)=(sinphi*sintheta*cospsi)-(cosphi*sinpsi)
rotmat(1,3)=(sinphi*sinpsi)+(cosphi*sintheta*cospsi)

!SECOND ROW
rotmat(2,1)=costheta*sinpsi
rotmat(2,2)=(cosphi*cospsi)+(sinphi*sintheta*sinpsi)
rotmat(2,3)=(cosphi*sintheta*sinpsi)-(sinphi*cospsi)

!THIRD ROW
rotmat(3,1)=-sintheta
rotmat(3,2)=sinphi*costheta
rotmat(3,3)=cosphi*costheta

 !  pi=4.0_wp*atan(1.0_wp)
 !  rads=(degrees*pi/180.0_wp)

return
end subroutine matrice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rotation (xx,yy,zz,rotmat)
implicit none
integer,PARAMETER            		   	:: wp = selected_real_kind(12,70)!(6,35) !(12,70)
real(kind=wp),                  intent(InOut)   ::xx,yy,zz
real(kind=wp), dimension(3,3),  intent(In  )    ::rotmat
integer						:: i
real(kind=wp)                  			::xxs,yys,zzs

!I insert the coordinates in a dummy vector, then i rotate
!do i=1,3
!write(*,*)rotmat(i,1),rotmat(i,2),rotmat(i,3)
!enddo
!write(*,*)"ora fa i conti"
!write(*,*)"vettore",xx,yy,zz
xxs=xx
yys=yy
zzs=zz
xx=(xxs*rotmat(1,1))+(yys*rotmat(1,2))+(zzs*rotmat(1,3))
yy=(xxs*rotmat(2,1))+(yys*rotmat(2,2))+(zzs*rotmat(2,3))
zz=(xxs*rotmat(3,1))+(yys*rotmat(3,2))+(zzs*rotmat(3,3))

!is possible also to use matmull
return
end subroutine rotation





























