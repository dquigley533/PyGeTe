! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                               E N E R G Y                                   !
!=============================================================================!
! Contains C-compatible Fortran 2008 routines which calculate the energy of a !
! system of particles interacting via the Zipoli Tersoff-based GeTe model.    !
!                                                                             !
! Key variables                                                               !
!    species(natoms)  : Array of species identifiers (1=Ge, 2=Te)             !
!    pos(3,natoms)    : Array of 3d particle positions (absolute coordinates) !
!    hmatrix(3,3)     : 3 vectors which define the supercell perioic lattice  !
!                                                                             !
! These will correspond to numpy arrays when invoked via a Python interface.  !
!=============================================================================!
module energy

  use iso_c_binding
  implicit None                                 ! Impose strong typing

  integer,parameter :: dp = c_double            ! C compatible double precision
  integer,parameter :: it = c_int               ! C compatible int

  Private                                       ! Everything is private
  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
                                                      !... unless exposed here.
  public :: compute_model_energy
  public :: compute_local_real_energy
  public :: compute_neighbour_list
  public :: fs,fc,tijk_fn,tijk_fast

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  public :: model_energy
  public :: natoms
  
  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  ! Set when number of atoms has been imposed
  logical,save :: initialised = .false.              

  ! Number of atoms in the lattice and number of lattices (for consistency checking)
  integer(kind=it),bind(c),save :: natoms = -1
  
  ! Current energy due to the model Hamiltonian
  real(kind=dp),bind(c),save :: model_energy = 0.0_dp

  ! Lattice translation vectors used to compute interactions with images
  integer :: nivect
  real(kind=dp),allocatable,dimension(:,:),save :: ivect

  !=============================================================!
  ! Constants defining the Zipoli Tersoff-based GeTe Potential  !
  ! Units of energy are eV, units of length are Angstrom        !
  !=============================================================!
  ! A & B energies indexed by species type. 
  real(kind=dp),dimension(0:1,0:1),save :: zipA = reshape((/6216.387_dp,5368.015_dp,5368.015_dp,4237.566_dp/), (/2,2/))
  real(kind=dp),dimension(0:1,0:1),save :: zipB = reshape((/417.086_dp,900.597_dp,900.597_dp,1001.053_dp/), (/2,2/))

  ! R & S matrices - inner and outer cut-off distances in the tapering function fc
  real(kind=dp),dimension(0:1,0:1),save :: zipR = reshape((/2.83633_dp,2.74930_dp,2.74930_dp,2.66494_dp/), (/2,2/))
  real(kind=dp),dimension(0:1,0:1),save :: zipS = reshape((/3.78306_dp,4.44433_dp,4.44433_dp,4.08355_dp/), (/2,2/))


  ! lambda and mu matices controlling the rate of decay of attractive and repsulsive terms with distance
  real(kind=dp),dimension(0:1,0:1),save :: zipLambda = reshape((/3.37634_dp,2.84763_dp,2.84763_dp,2.58006_dp/), (/2,2/))
  real(kind=dp),dimension(0:1,0:1),save :: zipMu     = reshape((/2.11859_dp,2.03009_dp,2.03009_dp,1.91781_dp/), (/2,2/))

  ! theta_0 - idea angle when each atom type is at the vertex
  real(kind=dp),dimension(0:1),save :: zipTheta0 = (/-0.9_dp,-0.9_dp/)

  ! Beta and Chi - appears in definition of coordination number
  real(kind=dp),dimension(0:1),save   :: zipBeta = (/6.4361E-7_dp,0.34202_dp/)
  !real(kind=dp),dimension(0:1),save   :: zipBeta = 0.0_dp
  real(kind=dp),dimension(0:1,0:1),save :: zipChi  = reshape((/1.0_dp,1.0_dp,1.0_dp,1.0_dp/), (/2,2/))

  ! d & h parameter in tijk (equation 6 of Zipoli et al). N.B. Table 1 defines these for
  ! X-Ge-Y and X-Te-Y which would mean X is atom i and Y atom k.
  real(kind=dp),dimension(0:1,0:1,0:1),save :: zipd = -1.0_dp, ziph = -1.0_dp, invzipd = -1.0_dp

  ! c & n parameters in tijk (equation 6 of Zipoli et al)
  real(kind=dp),dimension(0:1),save :: zipC = (/1.0832E6_dp,1.2028E3_dp/)
  real(kind=dp),dimension(0:1),save :: zipN = (/1.20657_dp,0.91154_dp/)
  !real(kind=dp),dimension(0:1),save :: zipN = 1.0_dp

  ! Parameters defining energy penalty for non-ideal coordination
  real(kind=dp),dimension(0:1),save :: zipC1 = (/0.01877_dp,0.10449_dp/)
  real(kind=dp),dimension(0:1),save :: zipC2 = (/0.09932_dp,0.22749_dp/)
  !real(kind=dp),dimension(0:1),save :: zipC1 = 0.0_dp
  !real(kind=dp),dimension(0:1),save :: zipC2 = 0.0_dp
  
  ! Defined in table 1 of Zipoli et al - relating to calculation of non-ideal coordination numbers
  real(kind=dp),parameter :: zT = 0.5_dp
  real(kind=dp),parameter :: zB = 0.5_dp
  real(kind=dp),parameter :: Zcut = 3.0_dp
  real(kind=dp),dimension(0:1),parameter :: z0 = (/4.41569_dp,3.45086_dp/)  ! species dependent

  ! Constant core energy per atom to match DFT data
  real(kind=dp),dimension(0:1),save :: E0 = (/-105.56237,-221.62927/)

  ! Number of atoms of each species type
  integer,dimension(0:1) :: nsp = (/-1,-1/)

  ! Adjacency matrix used to compute coordination number
  real(kind=dp),allocatable,dimension(:,:) :: coord
  
  ! Which imaging method to use (image vectors vs minimum image)
  real(kind=dp) :: r_cut     ! All atom-atom interactions zero beyond this range
  integer :: image_flag
  integer,parameter :: flag_iv = 1 , flag_mi = 2
  character(13),dimension(2) :: image_method = (/"Image vectors","Minimum image"/)
  
  ! Neighbour list. Maxneigh controls max no. neighbours per atom
  real(kind=dp),bind(c) :: skin = 1.0_dp ! Angstoms
  integer,save  :: maxneigh = 50
  integer,allocatable,dimension(:),save :: nn
  integer,allocatable,dimension(:,:),save ::jn,vn

contains

  subroutine energy_init(n,pos,hmatrix,species)
    !==============================================================================!
    ! Initialises all arrays and shared variables used in computing the energy of  !
    ! the system (or of a single particle). Called on first invocation of either   !
    ! compute_model_energy or compute_local_real_energy.                           !
    !------------------------------------------------------------------------------!
    ! D.Quigley February 2018                                                      !
    !------------------------------------------------------------------------------!
    use util,       only : util_determinant,util_images,util_recipmatrix,pi
    implicit none

    ! input variables
    integer,intent(in) :: n
    real(kind=dp),dimension(3,N),intent(in) :: pos
    real(kind=dp),dimension(3,3),intent(in) :: hmatrix
    integer,dimension(N),intent(in)         :: species
    
    integer :: im,jm,km,iat    ! loop counters
    integer :: ierr ! error flag

    real(kind=dp) :: density
    
    if (initialised) stop 'Error in energy_init: already initialised!'

    ! Populate 3-body interaction parameters (to move to reshaped parameter array)
!!$    zipd(0,0,0) = 25.61910_dp   !  Ge - Ge - Ge
!!$    zipd(0,0,1) = 28.18170_dp   !  Ge - Ge - Te
!!$    zipd(1,0,0) = 28.18170_dp   !  Te - Ge - Ge
!!$    zipd(1,0,1) = 28.07053_dp   !  Te - Ge - Te
!!$
!!$    zipd(0,1,0) = 27.86151_dp   !  Ge - Te - Ge
!!$    zipd(0,1,1) = 27.88612_dp   !  Ge - Te - Te
!!$    zipd(1,1,0) = 27.88612_dp   !  Te - Te - Ge
!!$    zipd(1,1,1) = 35.75317_dp   !  Te - Te - Te
!!$
!!$    ziph(0,0,0) = -0.44982_dp   !  Ge - Ge - Ge
!!$    ziph(0,0,1) = -0.35175_dp   !  Ge - Ge - Te
!!$    ziph(1,0,0) = -0.35175_dp   !  Te - Ge - Ge
!!$    ziph(1,0,1) = -0.34448_dp   !  Te - Ge - Te
!!$       
!!$    ziph(0,1,0) = -0.33105_dp   !  Ge - Te - Ge
!!$    ziph(0,1,1) = -0.37131_dp   !  Ge - Te - Te
!!$    ziph(1,1,0) = -0.37131_dp   !  Te - Te - Ge
!!$    ziph(1,1,1) = -0.46333_dp   !  Te - Te - Te


    ! Corrected ordering Sept 2018
    zipd(0,0,0) = 25.61910_dp   !  Ge - Ge - Ge
    zipd(0,0,1) = 28.18170_dp   !  Ge - Ge - Te
    zipd(0,1,0) = 28.18170_dp   !  Ge - Te - Ge
    zipd(0,1,1) = 28.07053_dp   !  Ge - Te - Te

    zipd(1,0,0) = 27.86151_dp   !  Te - Ge - Ge
    zipd(1,0,1) = 27.88612_dp   !  Te - Ge - Te
    zipd(1,1,0) = 27.88612_dp   !  Te - Te - Ge
    zipd(1,1,1) = 35.75317_dp   !  Te - Te - Te

    ziph(0,0,0) = -0.44982_dp   !  Ge - Ge - Ge
    ziph(0,0,1) = -0.35175_dp   !  Ge - Ge - Te
    ziph(0,1,0) = -0.35175_dp   !  Ge - Te - Ge
    ziph(0,1,1) = -0.34448_dp   !  Ge - Te - Te
       
    ziph(1,0,0) = -0.33105_dp   !  Te - Ge - Ge
    ziph(1,0,1) = -0.37131_dp   !  Te - Ge - Te
    ziph(1,1,0) = -0.37131_dp   !  Te - Te - Ge
    ziph(1,1,1) = -0.46333_dp   !  Te - Te - Te


    invzipd(:,:,:) = 1.0_dp/zipd(:,:,:)

    ! Set number of atoms for sanity checks on subsequent calls
    natoms = n

    ! Count the number of each species
    nsp = (/0,0/)
    do iat = 1,n
       nsp(species(iat)) = nsp(species(iat)) + 1
    end do

    ! Determine the image method and neighbour list info from the largest cut off
    r_cut = maxval(zipS)

    ! Check which imaging method to use
    image_flag = flag_mi ! minimum image
!    write(0,'("Forcing image vectors")')

    if ( 2.0*(r_cut + skin) >= sqrt(dot_product(hmatrix(:,1),hmatrix(:,1)))) image_flag = flag_iv
    if ( 2.0*(r_cut + skin) >= sqrt(dot_product(hmatrix(:,2),hmatrix(:,2)))) image_flag = flag_iv
    if ( 2.0*(r_cut + skin) >= sqrt(dot_product(hmatrix(:,3),hmatrix(:,3)))) image_flag = flag_iv  

    write(0,'("Using image method : ",A13)')image_method(image_flag)
    if ( image_flag == flag_iv ) then
       write(0,'("Warning - system too small for minimum image convention.")')
       write(0,'("Using image method only suitable for solids.")')
    end if

    if ( image_flag == flag_iv ) then
    
       ! Find out how many lattice vectors we need to consider all images within cutoff
       im = floor((r_cut+skin)/sqrt(dot_product(hmatrix(:,1),hmatrix(:,1))))+1
       jm = floor((r_cut+skin)/sqrt(dot_product(hmatrix(:,2),hmatrix(:,2))))+1
       km = floor((r_cut+skin)/sqrt(dot_product(hmatrix(:,3),hmatrix(:,3))))+1
       
       nivect = (2*im+1)*(2*jm+1)*(2*km+1)
       !write(0,'("DEBUG - number of ivects ",I4)')nivect

       if (allocated(ivect)) deallocate(ivect)
       allocate(ivect(1:3,1:nivect),stat=ierr)
       if (ierr/=0) stop 'Error in energy_init: Could not allocate ivect'
       
       ! compute current image translation vectors
       call compute_ivects(hmatrix)

    end if
       
    ! Neighbours

    ! Compute maxneigh based on density and a safety factor of 50%
    density = natoms/abs(util_determinant(hmatrix))
    maxneigh = ceiling((density*4.0*pi*(r_cut+skin)**3)/3.0_dp * 1.5_dp)

    if (allocated(nn)) deallocate(nn)
    allocate(nn(1:natoms),stat=ierr)
    if (ierr/=0) stop 'Error in energy_init: Could not allocate nn array'
    if (allocated(jn)) deallocate(jn)
    allocate(jn(1:maxneigh,1:natoms),stat=ierr)
    if (ierr/=0) stop 'Error in energy_init: Could not allocat jn array'

    if ( image_flag == flag_iv ) then
       if (allocated(vn)) deallocate(vn)
       allocate(vn(1:maxneigh,1:natoms),stat=ierr)
       if (ierr/=0) stop 'Error in energy_init: Could not allocate vn array'
    end if


    ! Adjacency matrix
    if (allocated(coord)) deallocate(coord)
    allocate(coord(1:natoms,1:natoms),stat=ierr)
    if (ierr/=0) stop 'Error allocating adjacency matrix when initialising module'
       
    ! Set initialised flag
    initialised = .true.

    ! Make sure first neighbour list is calculated 
    call compute_neighbour_list(natoms,3,pos,3,3,hmatrix,natoms,species)

    ! Make sure energy if up-to-date
    model_energy = compute_model_energy(natoms,3,pos,3,3,hmatrix,natoms,species)
    
    return

  end subroutine energy_init

  
  subroutine compute_ivects(hmatrix)
    !------------------------------------------------------------------------------!
    ! Computes the translation vectors needed to include all images such that each !
    ! atom sees all the images (including those of itself) within the cut off.     !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!
    implicit none
    real(kind=dp),dimension(3,3),intent(in) :: hmatrix
 
    !loop counters
    integer :: icell,jcell,kcell,im,jm,km,k
    real(kind=dp),dimension(3) :: sx,sy,sz

    im = floor((r_cut+skin)/sqrt(dot_product(hmatrix(:,1),hmatrix(:,1))))+1
    jm = floor((r_cut+skin)/sqrt(dot_product(hmatrix(:,2),hmatrix(:,2))))+1
    km = floor((r_cut+skin)/sqrt(dot_product(hmatrix(:,3),hmatrix(:,3))))+1
    
    nivect = (2*im+1)*(2*jm+1)*(2*km+1)

    ! we'd like the central cell to be entry 0
    ! we can flag it as non-self interacting
    ivect(:,1) = 0.0_dp
    
    k = 2
    do icell = -im,im
       sx = real(icell,kind=dp)*hmatrix(:,1)
       do jcell = -jm,jm
          sy = real(jcell,kind=dp)*hmatrix(:,2) 
          do kcell = -km,km
             sz = real(kcell,kind=dp)*hmatrix(:,3)
             
             if ( abs(icell)+abs(jcell)+abs(kcell) == 0 ) cycle
             ivect(:,k)  = sx + sy + sz
             k = k + 1

          end do
       end do
    end do

    return

  end subroutine compute_ivects

  real(kind=dp) function compute_local_real_energy(dmol,n,d,pos,dh2,dh1,hmatrix,n2,species) bind(c)
    !------------------------------------------------------------------------------!
    ! Call either the minimum image or image vector routine to compute the energy  !
    ! of a single particle interacting with its neighbours.                        !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!    
    implicit none
    integer(kind=it),value,intent(in)           :: dmol,d,n,dh1,dh2,n2
    real(kind=dp),intent(in),dimension(d,n)     :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix
    integer(kind=it),intent(in),dimension(n2)   :: species
    
    integer :: imol
    real(kind=dp) :: loc_energy

    ! Sanity checks
    if (d/=3) stop 'Error in compute_local_real_energy: position array must be 3xn'
    if ((dh1/=3).or.(dh2/=3)) stop 'Error in compute_local_real_energy: cell array must be 3x3'

    if (initialised) then
       if (n/=natoms) then
          write(0,'("WARNING in compute_local_real_energy: Number of entries in position array has changed!")')
          write(0,'("Re-initialising module")')
          initialised = .false.
          call energy_init(n,pos,hmatrix,species)
       end if
    else
       call energy_init(n,pos,hmatrix,species)
    end if

    ! We assume this routine will only ever be called from C/Python (never Fortran)
    ! and so we correct imol into 1-based indexing
    imol = dmol + 1
    if ( (imol<1).or.(imol>natoms) ) stop 'Error in compute_local_real_energy : atom index out of range!'
    
    select case (image_flag)
    case (flag_iv)
       call compute_ivects(hmatrix)
       ! Cells need to be so small here that we might as well just recompute everything
       loc_energy = compute_model_energy_iv(n,d,pos,dh2,dh1,hmatrix,n,species)
       !loc_energy = compute_local_real_energy_iv(imol,n,d,pos,dh2,dh1,hmatrix,n2,species)
    case (flag_mi)
       loc_energy = compute_local_real_energy_mi(imol,n,d,pos,dh2,dh1,hmatrix,n2,species)
    case default
       stop 'Error in compute_local_real_energy - unknown image flag'       
    end select

    compute_local_real_energy = loc_energy

  end function compute_local_real_energy

  real(kind=dp) function compute_model_energy(n,d,pos,dh2,dh1,hmatrix,n2,species) bind(c)
    !------------------------------------------------------------------------------!
    ! Call either the minimum image or image vector routine to compute the energy  !
    ! of the entire system.                                                        !
    !------------------------------------------------------------------------------!
    ! D.Quigley September 2017                                                     !
    !------------------------------------------------------------------------------!    
    implicit none
    integer(kind=it),value,intent(in)           :: d,n,dh1,dh2,n2
    real(kind=dp),intent(in),dimension(d,n)     :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix
    integer(kind=it),intent(in),dimension(n2  ) :: species

    real(kind=dp) :: loc_energy

    ! Sanity checks
    if (d/=3) stop 'Error in compute_local_real_energy: position array must be 3xn'
    if ((dh1/=3).or.(dh2/=3)) stop 'Error in compute_local_real_energy: cell array must be 3x3'

    if (initialised) then
       if (n/=natoms) then
          write(0,'("WARNING in compute_local_real_energy: Number of entries in position array has changed!")')
          write(0,'("Re-initialising module")')
          initialised=.false.
          call energy_init(n,pos,hmatrix,species)
       end if
    else
       call energy_init(n,pos,hmatrix,species)
    end if

    select case (image_flag)
    case (flag_iv)
       loc_energy = compute_model_energy_iv(n,d,pos,dh2,dh1,hmatrix,n,species)
    case (flag_mi)
       loc_energy = compute_model_energy_mi(n,d,pos,dh2,dh1,hmatrix,n,species)
    case default
       stop 'Error in compute_model_energy - unknown image flag'       
    end select

    compute_model_energy = loc_energy

  end function compute_model_energy


  real(kind=dp) function compute_model_energy_mi(n,d,pos,dh2,dh1,hmatrix,n2,species) 
    !------------------------------------------------------------------------------!
    ! Computes the energy of the entire simulation cell, using a pair potential.   !
    ! To be used when sampling the energy of the entire system or                  !
    ! when computing the energy change from a volume/cell trial move.              !
    ! This version uses the minimum image convention.                              !
    !------------------------------------------------------------------------------!
    ! D.Quigley February 2018                                                      !
    !------------------------------------------------------------------------------!
    use util, only : Util_determinant,Pi,util_recipmatrix,util_images
    implicit none
    integer,value,intent(in) :: d,n,dh1,dh2,n2
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix
    integer,intent(in),dimension(n2) :: species

    real(kind=dp),dimension(3)   :: tmpvect,tmpvect2
    real(kind=dp)                :: Eij,Ezip,eijk,izc
    real(kind=dp)                :: tijk,r2_ij,rcsq
    real(kind=dp)                :: tbpre,zeta,bij
    real(kind=dp)                :: brak,ctheta,deltaZ
    real(kind=dp)                :: r_ij,r2_ik,r_ik
    real(kind=dp)                :: inv_rij,argf,fcij
    
    real(kind=dp),dimension(3) :: ilj,jlj,klj

    integer :: imol,jmol,kmol ! loop counters
    integer :: ln,ln2,isp,jsp,ksp

    real(kind=dp),dimension(3,3) :: recip_matrix
    
    call util_recipmatrix(hmatrix,recip_matrix)
    
    Ezip  = 0.0_dp

    ! Zero coordination matrix
    !coord(:,:) = 0.0_dp

    !-----------------------------------------!
    !         Zipoli model                    !
    !-----------------------------------------!
    do imol = 1,natoms  ! loop over central atom imol

       isp = species(imol)    ! species
       ilj(:) = pos(:,imol)   ! position

       izc = 0.0_dp

       argf = Pi/(Pi-zipTheta0(isp))
       
       do ln = 1,nn(imol) ! loop over other atom jmol

          jmol = jn(ln,imol)   ! atom
          jsp = species(jmol)  ! species
          jlj(:) = pos(:,jmol) ! position

          rcsq = zipS(jsp,isp)*zipS(jsp,isp) ! squared cut-off for this species pair
          
          ! compute separation vector
          tmpvect = jlj(:) - ilj(:)

          call util_images(tmpvect,hmatrix,recip_matrix)
          
          r2_ij = tmpvect(1)**2 + tmpvect(2)**2 + tmpvect(3)**2
             
          ! compute interactions if in range 
          if ( r2_ij < rcsq ) then

             r_ij = sqrt(r2_ij)
             inv_rij = 1.0_dp/r_ij
             fcij = fc(jsp,isp,r_ij)
             
             ! Two body repulsive part
             Eij = zipA(jsp,isp)*exp(-zipLambda(jsp,isp)*r_ij)  
             
             ! Prefactor to attractive part
             tbpre = -zipB(jsp,isp)*exp(-zipMu(jsp,isp)*r_ij)
             
             zeta = 0.0_dp

             !-------------------------------------------------!
             ! Tersoff three body interaction jmol--imol--kmol !
             !-------------------------------------------------!
             do ln2 = nn(imol),1,-1
             
                kmol = jn(ln2,imol)
                if (kmol==jmol) cycle
                klj(:) = pos(:,kmol)    
                ksp = species(kmol)

                rcsq = zipS(ksp,isp)*zipS(ksp,isp)
                
                tmpvect2(:) = klj(:) - ilj(:)

                call util_images(tmpvect2,hmatrix,recip_matrix)
                
                r2_ik      = tmpvect2(1)**2 + tmpvect2(2)**2 + tmpvect2(3)**2

                ! compute three body term if in range
                if ( r2_ik < rcsq ) then
                   
                   r_ik = sqrt(r2_ik)

                   ! Note exp(mu*r_ij) also used above. Potential optimisation
                   eijk = zipMu(isp,jsp)*r_ij - zipMu(isp,ksp)*r_ik

                   ! LAMMPS has an overflow/underflow check here. Returns zero if the argument
                   ! of the exponential is less than -69.0776 and 1E30 if greater than 69.0776 
                   if (eijk > 69.0776_dp ) then
                      eijk = 1E30_dp
                   else if (eijk < -69.0776_dp) then
                      eijk = 0.0_dp
                   else
                      eijk = exp(eijk)
                   end if
                   
                   ! Angle between tmpvect and tmpvect2
                   ctheta = dot_product(tmpvect,tmpvect2)*inv_rij/r_ik
                   ctheta = ctheta*(1.0_dp-1.0E-12_dp)
                   
!!$                   theta  = acos(ctheta)                   
!!$                   arg    = argf*(theta - zipTheta0(isp)) ! constant denominator wrapped into argf
!!$                                                         
!!$                   ! Modified form of angle function - eq (11) of Zipoli paper
!!$                   tijk = 1.0_dp + zipC(isp)*zipC(isp)*invzipd(isp,jsp,ksp)*invzipd(isp,jsp,ksp)
!!$                   !tijk = 1.0_dp + zipC(isp)*zipC(isp)/(zipd(isp,jsp,ksp)*zipd(isp,jsp,ksp))
!!$                   denom = zipd(isp,jsp,ksp)*zipd(isp,jsp,ksp) + (ziph(isp,jsp,ksp)-cos(arg))**2
!!$                   tijk = tijk - zipC(isp)*zipC(isp)/denom
      
                   tijk = tijk_fn(ctheta,isp,jsp,ksp)
                   zeta = zeta + fc(isp,ksp,r_ik)*eijk*tijk
                   
                end if

             end do ! third body kmol

             brak = (zipBeta(isp)*zeta)**zipN(isp)

             ! Comparing with Tersoff implementation in LAMMPS, the Zipoli et al
             ! paper is missing a minus sign in the exponent here.
             bij = zipChi(jsp,isp)*(1.0_dp+brak)**(-0.5_dp/zipN(isp))
             Eij = 0.5_dp*(Eij + tbpre*bij)
             !-------------------------------------------------!
             ! End Tersoff three body interaction              !
             !-------------------------------------------------!

             ! Using fc here, suspecting typo in Zipoli paper which suggests f_s instead.
             !coord(jmol,imol) =  bij*fc(isp,jsp,r_ij)
             izc = izc + bij*fcij!(isp,jsp,r_ij)
             
             ! Add the energy of this pair into the total
             Ezip = Ezip + Eij*fcij!(isp,jsp,r_ij)
             
            
          end if ! jmol within range of jmol                  
          
       end do  ! end loop over neighbours of imol

       ! Total coordination number for atom imol
       !izc = sum(coord(:,imol))
       deltaZ = fs(isp,izc)

       ! Penalise deviations from idea coordination
       Ezip = Ezip + zipc1(isp)*deltaZ + zipc2(isp)*min(deltaZ*deltaZ,zcut*zcut)
       !Ezip = Ezip + zipc1(isp)*deltaZ + zipc2(isp)*deltaZ*deltaZ
    end do ! end loop over imol

    ! Add in constant term which only depends on number of atoms
    model_energy = Ezip + dot_product(nsp,E0)
    
    compute_model_energy_mi = model_energy
    
    return 

  end function compute_model_energy_mi

  real(kind=dp) function compute_local_real_energy_mi(imol,n,d,pos,dh2,dh1,hmatrix,n2,species)
    !------------------------------------------------------------------------------!
    ! Calculates the real-space contribution to the energy due particle dmol.      !
    ! To be used when computing the changes in energy due to a trial move.         !
    ! This version uses the minimum image convention
    !------------------------------------------------------------------------------!
    ! D.Quigley February 2018                                                      !
    !------------------------------------------------------------------------------!    
    use util,       only : util_determinant,util_images,util_recipmatrix,pi
    implicit none
    integer,value,intent(in) :: imol,d,n,dh1,dh2,n2
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix
    integer,intent(in),dimension(n2) :: species

    ! local variables
    real(kind=dp),dimension(3)   :: tmpvect,tmpvect2,tmpvect3
    real(kind=dp)                :: r2_ij,bijx,rcsq
    real(kind=dp)                :: Eij,Ezip,eijk,izc,Eji
    real(kind=dp)                :: tijk,theta,jzc,jzcx,denom
    real(kind=dp)                :: tbpre,zeta,arg,bij
    real(kind=dp)                :: brak,ctheta,deltaZ
    real(kind=dp)                :: r_ij,r2_ik,r_ik,r2_jl
    real(kind=dp)                :: r2_jk,r_jk,r_jl
    real(kind=dp)                :: tbpre2,zeta2,zeta3
   
    real(kind=dp),dimension(3) :: ilj,jlj,klj,llj

    real(kind=dp),dimension(3,3) :: recip_matrix

    integer :: jmol,kmol,lmol ! loop counters
    integer :: ln,ln2,ln3,isp,jsp,ksp,lsp,ti=0
    
    call util_recipmatrix(hmatrix,recip_matrix)

    ti = 0
    
    Ezip  = 0.0_dp

    !coord(:,imol) = 0.0_dp

    izc = 0.0_dp
    
    !-----------------------------------------!
    !         Zipoli model                    !
    !-----------------------------------------!
    ilj(:) = pos(:,imol)    ! position
    isp    = species(imol)  ! species

    do ln = 1,nn(imol) ! loop over other atom jmol

       jmol = jn(ln,imol)   ! atom
       jsp = species(jmol)  ! species
       jlj(:) = pos(:,jmol) ! position

       rcsq = zipS(jsp,isp)*zipS(jsp,isp) ! squared cut-off for this species pair
          
       ! compute separation vector
       tmpvect = jlj(:) - ilj(:)

       call util_images(tmpvect,hmatrix,recip_matrix)
                 
       r2_ij = tmpvect(1)**2 + tmpvect(2)**2 + tmpvect(3)**2
             
       ! compute interactions if in range 
       if ( r2_ij < rcsq ) then


          !coord(:,jmol) = 0.0_dp
          jzc  = 0.0_dp
          jzcx = 0.0_dp
          
          r_ij = sqrt(r2_ij)
            
          ! Two body repulsive part
          Eij = zipA(jsp,isp)*exp(-zipLambda(jsp,isp)*r_ij)  
             
          ! Prefactor to attractive part
          tbpre = -zipB(jsp,isp)*exp(-zipMu(jsp,isp)*r_ij)  

          zeta = 0.0_dp

          !-------------------------------------------------!
          ! Tersoff three body interaction jmol--imol--kmol !
          !-------------------------------------------------!
          do ln2 = 1,nn(imol)
          
             kmol = jn(ln2,imol)
             if (kmol==jmol) cycle
             
             klj(:) = pos(:,kmol)    
             ksp = species(kmol)
             
             rcsq = zipS(ksp,isp)*zipS(ksp,isp)
             
             tmpvect2(:) = klj(:) - ilj(:)

             call util_images(tmpvect2,hmatrix,recip_matrix)
             
             r2_ik      = tmpvect2(1)**2 + tmpvect2(2)**2 + tmpvect2(3)**2
             
             ! compute three body term if in range
             if ( r2_ik < rcsq ) then
                
                r_ik = sqrt(r2_ik)
                
                ! Note exp(mu*r_ij) also used above. Potential optimisation
                eijk = exp(zipMu(isp,jsp)*r_ij - zipMu(isp,ksp)*r_ik)
                
                ! Angle between tmpvect and tmpvect2
                ctheta = dot_product(tmpvect,tmpvect2)/(r_ik*r_ij)
                ctheta = ctheta*(1.0_dp-1.0E-12_dp)
                theta  = acos(ctheta)                   
                arg    = Pi*(theta - zipTheta0(isp))/(Pi-zipTheta0(isp)) ! constant denominator here
                
                ! Modified form of angle function - eq (11) of Zipoli paper
                tijk = 1.0_dp + zipC(isp)*zipC(isp)/(zipd(isp,jsp,ksp)*zipd(isp,jsp,ksp))
                denom = zipd(isp,jsp,ksp)*zipd(isp,jsp,ksp) + (ziph(isp,jsp,ksp)-cos(arg))*(ziph(isp,jsp,ksp)-cos(arg))
                tijk = tijk - zipC(isp)*zipC(isp)/denom
                                      
                zeta = zeta + fc(isp,ksp,r_ik)*eijk*tijk                 
                
             end if
             
          end do ! third body kmol
          
          brak = (zipBeta(isp)*zeta)**zipN(isp)
          bij = zipChi(jsp,isp)*(1.0_dp+brak)**(-0.5_dp/zipN(isp))
          Eij = 0.5_dp*(Eij + tbpre*bij) 
          !-------------------------------------------------!
          ! End Tersoff three body interaction              !
          !-------------------------------------------------!
          !coord(jmol,imol) =  bij*fc(isp,jsp,r_ij)
          izc = izc + bij*fc(isp,jsp,r_ij)

          tmpvect = -tmpvect ! now vector from j to i

          ! Two body repulsive part
          Eji = zipA(isp,jsp)*exp(-zipLambda(isp,jsp)*r_ij) 
             
          ! Prefactor to 3-body part
          tbpre = -zipB(isp,jsp)*exp(-zipMu(isp,jsp)*r_ij)

          zeta = 0.0_dp

          !-------------------------------------------------!
          ! Tersoff three body interaction imol--jmol--kmol !
          !-------------------------------------------------!
          !write(0,'(I5," has ",I5," neighbours")')jmol,nn(jmol)
          do ln2 = 1,nn(jmol)
          
             kmol = jn(ln2,jmol)
             if (kmol==imol) cycle
             klj(:) = pos(:,kmol)  
             ksp = species(kmol)
             
             rcsq = zipS(ksp,jsp)*zipS(ksp,jsp)

             tmpvect2(:) = klj(:) - jlj(:)

             call util_images(tmpvect2,hmatrix,recip_matrix)
             
             r2_jk      = tmpvect2(1)**2 + tmpvect2(2)**2 + tmpvect2(3)**2
             
             ! compute three body term if in range
             if ( r2_jk <= rcsq ) then

                r_jk = sqrt(r2_jk)
                
                !write(0,'(I5," is in range of ",I5," at distance ",F15.6)')kmol,jmol,r_jk
                


                ! Note exp(mu*r_ij) also used above. Potential optimisation
                eijk = exp(zipMu(jsp,isp)*r_ij - zipMu(jsp,ksp)*r_jk)
                
                ! Angle between tmpvect and tmpvect2
                ctheta = dot_product(tmpvect,tmpvect2)/(r_jk*r_ij)
                ctheta = ctheta * (1.0_dp-1.0E-12_dp)
                theta  = acos(ctheta)                
                arg    = Pi*(theta - zipTheta0(jsp))/(Pi-zipTheta0(jsp)) ! constant denominator here
                
                ! Modified form of angle function - eq (11) of Zipoli paper
                tijk = 1.0_dp + zipC(jsp)*zipC(jsp)/(zipd(jsp,isp,ksp)*zipd(jsp,isp,ksp))
                denom = zipd(jsp,isp,ksp)*zipd(jsp,isp,ksp) + (ziph(jsp,isp,ksp)-cos(arg))*(ziph(jsp,isp,ksp)-cos(arg))
                tijk = tijk - zipC(jsp)*zipC(jsp)/denom
                
                zeta = zeta + fc(jsp,ksp,r_jk)*eijk*tijk

                tbpre2 = -zipB(jsp,ksp)*exp(-zipMu(jsp,ksp)*r_jk)
                zeta2  = 0.0_dp
                zeta3  = 0.0_dp

                !-------------------------------------------------!
                ! Tersoff three body interaction kmol--jmol--lmol !
                !-------------------------------------------------!
                ! Loop over other neighbours of jmol
                ! (third body lmol might include imol)
                do ln3 = 1,nn(jmol)

                   lmol = jn(ln3,jmol)
                   if (lmol==kmol) cycle
                   llj(:) = pos(:,lmol)    
                   lsp = species(lmol)

                   rcsq = zipS(lsp,jsp)*zipS(lsp,jsp)

                   tmpvect3(:) = llj(:) - jlj(:)

                   call util_images(tmpvect3,hmatrix,recip_matrix)
                                
                   r2_jl      = tmpvect3(1)**2 + tmpvect3(2)**2 + tmpvect3(3)**2

                   ! compute three body term if in range
                   if ( r2_jl < rcsq ) then
                      
                      r_jl = sqrt(r2_jl)

                      eijk = exp(zipMu(ksp,jsp)*r_jk - zipMu(jsp,lsp)*r_jl)
                
                      ! Angle between tmpvect and tmpvect2
                      ctheta = dot_product(tmpvect2,tmpvect3)/(r_jk*r_jl)
                      ctheta = ctheta * (1.0_dp-1.0E-12_dp)
                      theta  = acos(ctheta)                
                      arg    = Pi*(theta - zipTheta0(jsp))/(Pi-zipTheta0(jsp)) ! constant denominator here
                      
                      ! Modified form of angle function - eq (11) of Zipoli paper
                      tijk = 1.0_dp + zipC(jsp)*zipC(jsp)/(zipd(jsp,ksp,lsp)*zipd(jsp,ksp,lsp))
                      denom = zipd(jsp,ksp,lsp)*zipd(jsp,ksp,lsp) + (ziph(jsp,ksp,lsp)-cos(arg))*(ziph(jsp,ksp,lsp)-cos(arg))
                      tijk = tijk - zipC(jsp)*zipC(jsp)/denom
                      
                      zeta2 = zeta2 + fc(jsp,lsp,r_jl) *eijk*tijk

                      ! Could just compute imol case seperately outside loop and add in later?
                      ! Need to accumulate zeta without imol contribution also to get difference in bij due to i
                      if (lmol/=imol) zeta3 = zeta3 + fc(jsp,lsp,r_jl)*eijk*tijk
                      
                   end if

                end do

                brak = (zipBeta(jsp)*zeta2)**zipN(jsp)
                bij = zipChi(jsp,ksp)*(1.0_dp+brak)**(-0.5_dp/zipN(jsp))

                jzc = jzc + bij*fc(jsp,ksp,r_jk)
                
                brak = (zipBeta(jsp)*zeta3)**zipN(jsp)
                bijx = zipChi(jsp,ksp)*(1.0_dp+brak)**(-0.5_dp/zipN(jsp))

                jzcx = jzcx + bijx*fc(jsp,ksp,r_jk)
                
                Ezip = Ezip + 0.5_dp*(Tbpre2*(bij-bijx))*fc(jsp,ksp,r_jk)
                !write(0,*)ti,0.5_dp*(Tbpre2*bij)*fc(jsp,ksp,r_jk)
                ti = ti + 1
  
             end if
                                
          end do ! third body kmol
          
          brak = (zipBeta(jsp)*zeta)**zipN(jsp)
          bij = zipChi(isp,jsp)*(1.0_dp+brak)**(-0.5_dp/zipN(jsp))

          Eji = 0.5_dp*(Eji + Tbpre*bij)

          !-------------------------------------------------!
          ! End Tersoff three body interaction              !
          !-------------------------------------------------!

!          coord(imol,jmol) =  bij*fc(jsp,isp,r_ij)
          jzc = jzc + bij*fc(jsp,isp,r_ij)
          
          ! Add in the two energy contributions computed above
          Ezip = Ezip + (Eij+Eji)*fc(isp,jsp,r_ij)

          !jzc = sum(coord(:,jmol))

          
          ! Penalise deviations from idea coordination
          deltaZ = fs(jsp,jzc)
          Ezip = Ezip + zipc1(jsp)*deltaZ + zipc2(jsp)*min(deltaZ*deltaZ,zcut*zcut)

          deltaZ = fs(jsp,jzcx)
          Ezip = Ezip - zipc1(jsp)*deltaZ - zipc2(jsp)*min(deltaZ*deltaZ,zcut*zcut)
          
          
       end if ! jmol within range of imol          

    end do  ! end loop over neighbours of imol

    ! Coordination number for imol
    !izc = sum(coord(:,imol))
    deltaZ = fs(isp,izc)

    ! Penalise deviations from idea coordination
    Ezip = Ezip + zipc1(isp)*deltaZ + zipc2(isp)*min(deltaZ*deltaZ,zcut*zcut)

    ! set return value - all energy involving position of atom imol
    compute_local_real_energy_mi = Ezip
    
    return 

  end function compute_local_real_energy_mi


  real(kind=dp) function compute_model_energy_iv(n,d,pos,dh2,dh1,hmatrix,n2,species)
    !------------------------------------------------------------------------------!
    ! Computes the energy of the entire simulation cell, using a pair potential.   !
    ! To be used when sampling the energy of the entire system or                  !
    ! when computing the energy change from a volume/cell trial move.              !
    ! This version uses stored vectors to consider multiple images.                !
    !------------------------------------------------------------------------------!
    ! D.Quigley February 2018                                                      !
    !------------------------------------------------------------------------------!
    use util, only : Util_determinant, Pi
    implicit none
    integer,value,intent(in) :: d,n,dh1,dh2,n2
    real(kind=dp),intent(in),dimension(d,n) :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix
    integer,intent(in),dimension(n2) :: species
    

    real(kind=dp),dimension(3)   :: tmpvect,tmpvect2
    real(kind=dp)                :: Eij,Ezip,eijk,izc
    real(kind=dp)                :: tijk,theta,denom
    real(kind=dp)                :: tbpre,zeta,arg,bij
    real(kind=dp)                :: brak,ctheta,deltaZ
    real(kind=dp)                :: rcsq,r2_ij
    real(kind=dp)                :: r_ij,r2_ik,r_ik
    
    real(kind=dp),dimension(3) :: ilj,jlj,klj

    integer :: imol,jmol,kmol ! loop counters
    integer :: ji,ki,ln,ln2,isp,jsp,ksp
    
    call compute_ivects(hmatrix)

    Ezip  = 0.0_dp

    ! Zero coordination matrix
    !coord(:,:) = 0.0_dp
    
    !-----------------------------------------!
    !         Zipoli model                    !
    !-----------------------------------------!
    do imol = 1,natoms  ! loop over central atom imol

       isp = species(imol)    ! species
       ilj(:) = pos(:,imol)   ! position

       izc = 0.0_dp

       do ln = 1,nn(imol) ! loop over other atom jmol

          jmol = jn(ln,imol)   ! atom
          ji   = vn(ln,imol)   ! image
          jsp = species(jmol)  ! species

          ! position
          jlj(:) = pos(:,jmol) + ivect(:,ji)

          rcsq = zipS(jsp,isp)*zipS(jsp,isp) ! squared cut-off for this species pair

          ! compute separation vector
          tmpvect = jlj(:) - ilj(:)
          r2_ij   = tmpvect(1)**2 + tmpvect(2)**2 + tmpvect(3)**2
             
          ! compute interactions if in range with the current image
          if ( r2_ij < rcsq ) then

             r_ij = sqrt(r2_ij)
            
             ! Two body repulsive part
             Eij = zipA(jsp,isp)*exp(-zipLambda(jsp,isp)*r_ij)  
             
             ! Prefactor to 3-body part
             tbpre = -zipB(jsp,isp)*exp(-zipMu(jsp,isp)*r_ij)

             zeta = 0.0_dp

             !-------------------------------------------------!
             ! Tersoff three body interaction jmol--imol--kmol !
             !-------------------------------------------------!
             !do ln2 = ln+1,nn(imol) ! third body kmol
             do ln2 = 1,nn(imol)                
             
                kmol = jn(ln2,imol)  ! atom
                ki   = vn(ln2,imol)  ! image

                if ( (kmol==jmol) .and. (ki==ji) ) cycle
                
                ksp  = species(kmol) ! species

                rcsq = zipS(ksp,isp)*zipS(ksp,isp)
                
                ! position
                klj(:) = pos(:,kmol) + ivect(:,ki)   
                      
                tmpvect2(:) = klj(:) - ilj(:)
                r2_ik      = tmpvect2(1)**2 + tmpvect2(2)**2 + tmpvect2(3)**2
                
                ! compute three body term if in range
                if ( r2_ik < rcsq ) then
                   
                   r_ik = sqrt(r2_ik)

                   ! Note exp(mu*r_ij) also used above. Potential optimisation
                   eijk = exp(zipMu(isp,jsp)*r_ij - zipMu(isp,ksp)*r_ik)

                   ! Angle between tmpvect and tmpvect2
                   ctheta = dot_product(tmpvect,tmpvect2)/(r_ik*r_ij)
                   ctheta = ctheta*(1.0_dp-1.0E-12_dp)
                   theta  = acos(ctheta)                   
                   arg    = Pi*(theta - zipTheta0(isp))/(Pi-zipTheta0(isp)) ! constant denominator here
                                                         
                   ! Modified form of angle function - eq (11) of Zipoli paper
                   tijk = 1.0_dp + zipC(isp)*zipC(isp)/(zipd(isp,jsp,ksp)*zipd(isp,jsp,ksp))
                   denom = zipd(isp,jsp,ksp)*zipd(isp,jsp,ksp) + (ziph(isp,jsp,ksp)-cos(arg))*(ziph(isp,jsp,ksp)-cos(arg))
                   tijk = tijk - zipC(isp)*zipC(isp)/denom
                                      
                   zeta = zeta + fc(isp,ksp,r_ik)*eijk*tijk               
                   
                end if

             end do ! third body kmol

             brak = (zipBeta(isp)*zeta)**zipN(isp)

             ! Comparing with Tersoff implementation in LAMMPS, the Zipoli et al
             ! paper is missing a minus sign in the exponent here.
             bij = zipChi(jsp,isp)*(1.0_dp+brak)**(-0.5_dp/zipN(isp))
             Eij = 0.5_dp*(Eij + tbpre*bij)
             !-------------------------------------------------!
             ! End Tersoff three body interaction              !
             !-------------------------------------------------!

             ! Need to increment rather than set this as imol could
             ! interact more than once with jmol
             !coord(jmol,imol) = coord(jmol,imol) + bij*fc(isp,jsp,r_ij)
             izc = izc + bij*fc(isp,jsp,r_ij)
             Ezip = Ezip + Eij*fc(isp,jsp,r_ij)
             
          end if ! jmol within range of jmol inside image ji

       end do  ! end loop over neighbours of imol

       ! Here we vist every atom as imol, but for single particle updates we'll need
       ! to compute this term for every jmol as well as imol.
       !izc = sum(coord(:,imol))
       deltaZ = fs(isp,izc)

       !print*,"Coordination number of atom ",imol," = ",izc
       
       ! Penalise deviations from idea coordination
       Ezip = Ezip + zipc1(isp)*deltaZ + zipc2(isp)*min(deltaZ*deltaZ,zcut*zcut)

       
    end do ! end loop over imol

    ! Add in constant term which only depends on number of atoms
    model_energy = Ezip + dot_product(nsp,E0)
    
    compute_model_energy_iv = model_energy
    
    return 

  end function compute_model_energy_iv


  subroutine  compute_neighbour_list(n,d,pos,dh2,dh1,hmatrix,n2,species) bind(c)
    !------------------------------------------------------------------------------!
    ! Call either the minimum image or image vector routine to compute the energy  !
    ! of the entire system.                                                        !
    !------------------------------------------------------------------------------!
    ! D.Quigley February 2018                                                     !
    !------------------------------------------------------------------------------!    
    implicit none
    integer(kind=it),value,intent(in)           :: d,n,dh1,dh2,n2
    real(kind=dp),intent(in),dimension(d,n)     :: pos
    real(kind=dp),intent(in),dimension(dh1,dh2) :: hmatrix
    integer(kind=it),intent(in),dimension(n2)   :: species

    ! Sanity checks
    if (d/=3) stop 'Error in compute_local_real_energy: position array must be 3xn'
    if ((dh1/=3).or.(dh2/=3)) stop 'Error in compute_local_real_energy: cell array must be 3x3'

    if (initialised) then
       if (n/=natoms) then
          write(0,'("WARNING in compute_local_real_energy: Number of entries in position array has changed!")')
          write(0,'("Re-initialising module")')
          initialised = .false.
          call energy_init(n,pos,hmatrix,species)
       end if
    else
       call energy_init(n,pos,hmatrix,species)
    end if

    select case (image_flag)
    case (flag_iv)
       call compute_neighbour_list_iv(n,d,pos,dh2,dh1,hmatrix)
    case (flag_mi)       
       call compute_neighbour_list_mi(n,d,pos,dh2,dh1,hmatrix)
    case default
       stop 'Error in compute_neighbour_list- unknown image flag'       
    end select

  end subroutine compute_neighbour_list
  
  subroutine compute_neighbour_list_iv(n,d,pos,dh2,dh1,hmatrix) 
    !------------------------------------------------------------------------------!
    ! Computes a Verlet neighbour list for use in calculating energy               !
    ! This version uses stored vectors to consider multiple images.                !
    !------------------------------------------------------------------------------!
    ! D.Quigley February 2018                                                      !
    !------------------------------------------------------------------------------!
    implicit none
    integer,value,intent(in) :: n,d,dh2,dh1
    real(kind=dp),dimension(d,n),intent(in) :: pos
    real(kind=dp),dimension(dh1,dh2),intent(in) :: hmatrix

    integer :: imol,jmol,k,ni
    real(kind=dp) :: r2_ij,rn
    real(kind=dp),dimension(3) :: v_ij,ilj,jlj,tmpvect

    ! Hard coded neighbour list cutoff with skin width 
    rn = r_cut + skin
    
    call compute_ivects(hmatrix)

    do imol = 1,n

       ilj(:) = pos(:,imol)

       nn(imol) = 0
       do jmol = 1,n

          jlj(:) = pos(:,jmol)

          v_ij(:)   = jlj(:) - ilj(:)

          do k = 1,nivect
             if ( (k==1).and.(jmol==imol) )cycle

             tmpvect(:) = v_ij(:) + ivect(:,k) ! apply image             
             r2_ij = dot_product(tmpvect,tmpvect)

             if (r2_ij<rn*rn) then
                nn(imol) = nn(imol) + 1
                ni = nn(imol)        ! number of neighbours of imol
                jn(ni,imol) = jmol   ! jmol is a neighbour of imol
                vn(ni,imol) = k      ! in image k
             end if

          end do
          if (nn(imol) == maxneigh) stop 'Error in compute_neighbour _list : maximum neighbours exceeded!'
       end do
       
       !write(0,'("Molecule ",I5," has ",I5," neighbours")')imol,nn(imol,ils)

    end do


  end subroutine compute_neighbour_list_iv


  subroutine compute_neighbour_list_mi(n,d,pos,dh2,dh1,hmatrix) 
    !------------------------------------------------------------------------------!
    ! Computes a Verlet neighbour list for use in calculating energy               !
    ! This version uses the minimum image convention.                              !
    !------------------------------------------------------------------------------!
    ! D.Quigley February 2018                                                      !
    !------------------------------------------------------------------------------!
    use util, only : util_recipmatrix,util_images
    implicit none
    integer,value,intent(in) :: n,d,dh2,dh1
    real(kind=dp),dimension(d,n),intent(in) :: pos
    real(kind=dp),dimension(dh1,dh2),intent(in) :: hmatrix

    integer :: imol,jmol,ni
    real(kind=dp) :: r2_ij,rn
    real(kind=dp),dimension(3) :: v_ij,ilj,jlj

    real(kind=dp),dimension(3,3) :: recip_matrix

    call util_recipmatrix(hmatrix,recip_matrix)

    ! Hard coded neighbour list cutoff with skin width 
    rn = r_cut + skin
    
    do imol = 1,n

       ilj(:) = pos(:,imol)

       nn(imol) = 0

       do jmol = 1,n

          if (imol==jmol) cycle

          jlj(:) = pos(:,jmol)

          v_ij(:)   = jlj(:) - ilj(:)

          call util_images(v_ij,hmatrix,recip_matrix)
          
          r2_ij = dot_product(v_ij,v_ij)
          
          if (r2_ij<rn*rn) then
             nn(imol) = nn(imol) + 1
             ni = nn(imol)        ! number of neighbours of imol
             jn(ni,imol) = jmol   ! jmol is a neighbour of imol
          end if
             
          if (nn(imol) == maxneigh) stop 'Error in compute_neighbour _list : maximum neighbours exceeded!'

       end do
       
    end do


  end subroutine compute_neighbour_list_mi

  real(kind=dp) function fc(ispc,jspc,r) bind(c)
    !------------------------------------------------------------------------------!
    ! Implments the Tersoff-style tapering function used in the Zipoli et al       !
    ! potential. The inner and outer cuf off distances are a funciton of species   !
    ! and so the two species involved need to be specified as arguments.           !
    !------------------------------------------------------------------------------!
    ! D.Quigley February 2018                                                      !
    !------------------------------------------------------------------------------!
    use util, only : Pi
    implicit none
    integer(kind=it),intent(in),value :: ispc, jspc   ! Species types
    real(kind=dp),intent(in),value    :: r            ! Radial separation
    real(kind=dp) :: Rcut,Scut
    
    Scut = zipS(ispc,jspc)
    Rcut = zipR(ispc,jspc)

    ! Might want to tabulate this for efficiency
    if ( r <= Rcut ) then
       fc = 1.0_dp
    elseif ( r <= Scut ) then
       ! Constant per-species factor in cosine argument...
       fc =  0.5_dp*(1.0_dp+cos(Pi*(r-Rcut)/(Scut-Rcut))) 
    else
       fc = 0.0_dp
    end if

    return 
    
  end function fc


  real(kind=dp) function fs(ispc,zi) bind(c)
    !------------------------------------------------------------------------------!
    ! Implments the function fs of Zipoli et al equation (10) which computes a     !
    ! species-specific deviation from ideal coordination.                          !
    !------------------------------------------------------------------------------!
    ! D.Quigley February 2018                                                      !
    !------------------------------------------------------------------------------!
    use util, only : Pi
    implicit none
    integer(kind=it),intent(in),value :: ispc   ! Species type
    real(kind=dp),intent(in),value :: zi        ! Coordination number

    real(kind=dp) :: z,zd,zint
    
    zd = abs(zi-z0(ispc))
    zint = int(zd,kind=dp)
    z = zd - zint

    if ( z <= zT - zB ) then
       fs = zint
    elseif ( z <= zT + zB ) then
       fs = sign(zint + 0.5_dp*(1.0_dp+sin(0.5_dp*Pi*(z-zT)/zB)),zi-z0(ispc))
       !fs = zint + 0.5_dp*(1.0_dp+sin(0.5_dp*Pi*(z-zT)/zB))
    else
       fs = 1.0_dp
    end if

    return
    
  end function fs
  
  real(kind=dp) function tijk_fn(ctheta,isp,jsp,ksp) bind(c,name='tijk')
    !------------------------------------------------------------------------------!
    ! Implments the function tijk of Zipoli et al equation (10) which computers a  !
    ! species-specific three body interation as function of the apex angle.        !
    !------------------------------------------------------------------------------!
    ! D.Quigley February 2018                                                      !
    !------------------------------------------------------------------------------!
    use util, only : Pi
    implicit none

    real(kind=dp),intent(in),value    :: ctheta
    integer(kind=it),intent(in),value :: isp, jsp, ksp
    real(kind=dp) :: arg,denom,theta
    
    ! Modified form of angle function - eq (11) of Zipoli paper
    theta  = acos(ctheta)                
    arg    = Pi*(theta - zipTheta0(jsp))/(Pi-zipTheta0(jsp)) ! constant denominator here

    tijk_fn = 1.0_dp + zipC(isp)*zipC(isp)*invzipd(isp,jsp,ksp)*invzipd(isp,jsp,ksp)
    denom = zipd(isp,jsp,ksp)*zipd(isp,jsp,ksp) + (ziph(isp,jsp,ksp)-cos(arg))**2
    tijk_fn = tijk_fn - zipC(isp)*zipC(isp)/denom

    return
    
  end function tijk_fn

  real(kind=dp) function tijk_fast(ctheta,isp,jsp,ksp) bind(c)
    !------------------------------------------------------------------------------!
    ! Implments the function tijk of Zipoli et al equation (10) which computers a  !
    ! species-specific three body interation as function of the apex angle.        !
    !------------------------------------------------------------------------------!
    ! D.Quigley February 2018                                                      !
    !------------------------------------------------------------------------------!
    implicit none

    real(kind=dp),intent(in),value    :: ctheta
    integer(kind=it),intent(in),value :: isp, jsp, ksp
    integer :: ic,jk

    real(kind=dp),dimension(0:7,0:2,0:1) :: coeff

    coeff(:,0,0) = (/745495.10150827_dp,783327.91492384_dp,-470664.83156142_dp, &
                    -320388.95410986_dp,615431.18188016_dp,1632699.51029188_dp, &
                    414030.54164159_dp,27592.25394563_dp/) ! Ge-Ge-Ge

    coeff(:,1,0) = (/483955.06757535_dp,509088.4305563_dp,-304057.82564713_dp, &
                    -205511.66972008_dp,406518.31653139_dp, 1111489.85938716_dp, &
                    -40955.83587062_dp,-2109.42745721_dp/)  ! Ge-Te-Ge or Ge-Te-Ge

    coeff(:,2,0) = (/4.98944061e+05_dp, 5.24644655e+05_dp,-3.14021413e+05_dp, &
                    -2.12783786e+05_dp, 4.16451534e+05_dp,1.12644228e+06_dp, &
                     7.30012171e+04_dp, -8.60019592e+02_dp/) ! Ge-Te-Te
    

    coeff(:,0,1) = (/0.58025665_dp, 0.61026363_dp,-0.36489184_dp, -0.24695249_dp, &
         0.48577056_dp, 1.32073735_dp, 0.02092988_dp,  0.99734004_dp/) ! Te-Ge-Ge

    coeff(:,1,1) = (/0.61520718_dp, 0.64689648_dp, -0.38719416_dp, -0.26236631_dp, &
                     0.51349238_dp, 1.38892399_dp,  0.09001184_dp,  0.99893835_dp/)


    coeff(:,2,1) = (/0.24545735_dp, 0.25789401_dp, -0.15501728_dp, -0.10555383_dp, &
                     0.20205017_dp, 0.53372805_dp,  0.15200562_dp,  1.01167768_dp/) ! Te-Te-Te
    
    ! Polynomial fit
    jk = ksp+jsp ! zero (x-Ge-Ge), one (x-Te-Ge or x-Ge-Te), two (x-Te-Te)

    !write(*,'("isp = ",I5," jk = ",I5)')isp,jk
    
    ! This should be unrolled, do so manually if it helps
    ic = 0
    tijk_fast = coeff(7,jk,isp)
    do ic = 0,6
       tijk_fast = tijk_fast + coeff(ic,jk,isp)*(ctheta**(7-ic))
    end do

    
    return
    
  end function tijk_fast
  
                   
end module energy
