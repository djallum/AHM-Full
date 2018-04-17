! Log:
! 171020 ----- Created
! 171105 ----- Looking to add cluster data (see research update)
! 171115 ----- Corrected issues with cluster counting, started distribution of bond strengths
! 171117 ----- Looking to collect and save the distribution of bond strengths.
! 171120 ----- Changed how i calculate the distribution of cluster sizes by weighting each contribution by the number of sites in that cluster. ie; for each 5 site cluster counted, add 5 to the element that corresponds to clusters of 5 sites
! 171204 ----- Added a new function which increments the file names by 1 if there exists a file of the same name already in the folder.




  ! This is the first version of my code of an RG routine on the Anderson-Hubbard model
  ! Created: 171018; Last Edited: 171018
  ! Purpose: Extract the DoS from the AHM with strong disorder (t/W < 1) and strong interactions (t/U < 1) using a simplified RG similar to the method implemented by Eamonn Campbell in his Thesis and Johri and Bhatt (2014). Steps are below
  ! (1) Set up system of 'nSite' sites, each with a site potential stored in 'potSite'. The values for the site potentials are selected randomly from a uniform distribution of width 'DELTA' from [-DELTA/2,DELTA/2]. Between each site, a hopping ampitude is defined which governs the overlap between neighbouring orbitals, these are stored in 'hopSite' where the i'th element is the hopping amplitude from site i to i+1. Initially, the hopping amplitude is the same for each site but is renomalized during the RG process. An interaction energy is defined which is the energy cost to having two electrons on the same site. This will be stored in 'uSite' and is the same for each site.
  ! (2) Identify all the weak bonds in the cluster. A bond between sites is defined by abs(t/E_ij) (E_ij is the difference between site potentials of neighbouring sites). E_ij is chosen as the smallest of |E_i-E_j|, |(E_i+U)-E_j| and |E_i-(E_j+U)| where E_i/j+U is the upper hubbard orbital for site i/j
  ! (3) The sites between these bonds are called clusters which for the purpose of this code are independent from the lattice and are diagonalized as their own systems. This code cannot diagonalize clusters above 8 sites, so the disorder strength needs to be large so that cluster sizes are small.
  ! (4) store the DoS contributions from each cluster
  ! Will do this for many different systems to make the DoS an ensemble average.

  ! Inputs:
  ! Physical Parameters:
  !          Disorder Strength: DELTA
  !          Hopping Amplitude: hop
  !          interaction strength: uSite
  !          system size: dim
  ! Programming Parameters:
  !          number of systems: systemn
  !          bond_cutoff: the value both bonds need to be below before a site or a cluster is considered for elimination


module SHARED_MODULE
  USE Diag
  IMPLICIT NONE                    
  SAVE
  integer, parameter :: systemn = 100000 !number of systems
  integer, parameter :: dim = 1000 !system size 
  REAL, parameter :: DELTA = 10 !disorder strength
  REAL, parameter :: hop = 0 !hopping amplitude
  REAL, parameter :: bond_cutoff = 0.001 !maximum bond strength allowed for a site or cluster to be removed
  REAL, parameter :: uSite = 6 !Interaction strength
  REAL, parameter :: ChemPot = uSite/2 !Chemical Potential
  integer, parameter :: n_d = systemn*dim !number of systems*number of sites in a single system
  !-------------DoS binning parameters---------------------
  integer, parameter :: bins = 1000 !number of bins in the density of states
  real, parameter :: DoSMax = 2*uSite   ! maxval((/ DELTA/2 + 1, uSite + 1 /))
  real, parameter :: DoSMin = -DoSMax
  integer, parameter :: DOSCl = 0 !I want to collect DOS from each cluster size. I will collect DOS contributions from clusters of size 1 to DOSCl 
  integer, parameter :: ClusterMax = 6
  real DoS(bins,1+DOSCl)
  real :: DroppedDos, TotalDos
  integer iii
  real, parameter :: BinsDoS(bins) = [(DoSMin + (DoSMax-DoSMin)*iii/real(bins), iii = 1,bins)]   
  !-------------Run This part? -------------------------
  logical, parameter :: DoClusterCount = .false. !if = 1,do cluster count/distribution stuff. If = 0, don't
  logical, parameter :: DoBondStrength = .true.
  logical, parameter :: DoEnergyDiff = .false.
  !-------------Potential binning parameters---------------
  real Potential(bins)
  real, parameter :: PotMax = Delta/2
  real, parameter :: PotMin = -PotMax
  real DroppedPot 
  real :: BinsPot(bins)
contains
!
! This subroutine is from Eamonn's code and it defines a start point for the seed of the random number generator based on the time that the code is run
!
  subroutine init_random_seed()
    
    implicit none    
    integer,dimension(:), allocatable :: seed !seed for random number generator
    integer  e,clock,i5  !clock is the system time in seconds, e is the size of seed of the random number generator
    call random_seed (size = e) !find out how big the seed is to be. This size of seed is e
    allocate (seed(e)) !allocate space for seed
    call system_clock (count= clock) !obtain the system time, in seconds
    seed = clock + 37*(/(i5 - 1,i5 = 1, e)/) !set the seed using the system time
    call random_seed (put=seed) !initialize the random number generator
    deallocate (seed) !deallocate the seed array
   return
  end subroutine init_random_seed

!
  ! This subroutine creates the system, generates the site potentials, calculates bond strengths, stores hopping in array, might make it find all the clusters too
  ! Test002 shows that the site potentials are being selected correctly and the strongest bond has been chosen.
!

  subroutine create_AHM( SitePotential, Hopping, Bonds )
    implicit none
                !---------------------------Outputs------------------------------------------
    real, dimension(:), allocatable, intent(out) :: SitePotential 
    real, dimension(:), allocatable, intent(out) :: Hopping
    real, dimension(:), allocatable, intent(out) :: Bonds

                !---------------------------Coding Tools-------------------------------------
    real RandomNumber                                    ! random number for site potential calculation
    integer loop1                                        ! loop integer 
    integer IndexOfNeighbour                             ! this stores the site number of the right neighbour site in a bond. Used for
                                                         ! boundary conditions

                !---------------------------Variable Allocations-----------------------------
    allocate(SitePotential(dim))
    allocate(Hopping(dim))
    allocate(Bonds(dim))

    
                !---------------------------Initializing-------------------------------------
    RandomNumber = 0
    Hopping(1:dim) = 0 
    SitePotential(1:dim) = 0
    Bonds(1:dim) = 0

                !---------------------------Calculate Site Potentials/Hopping Amplitudes-----
    !This loop uses the RNG to determine the site potentials and initializes all the hopping potentials as the same.
    do loop1 = 1,dim
       call RANDOM_NUMBER(RandomNumber)                  ! Sets 'RandomNumber' to a random number between [0,1]
       if ( mod(loop1,2) .eq. 0 ) then
          Hopping(loop1) = 0                              ! Sets the hopping amplitude equal to the parameter 'hop'
       else
          Hopping(loop1) = hop
       end if
       SitePotential(loop1) = (0.5 - RandomNumber)*DELTA ! Shifts 'RandomNumber' from interval [0,1] to [-DELTA,DELTA]
    end do
    !print*, Hopping
                !---------------------------Calculate the 'Bonds' array----------------------
    do loop1 = 1,dim 
       IndexOfNeighbour = loop1+1                        ! Labels the neighbour site for the purposes of the boundary conditions

       if (loop1 .eq. dim) then                          ! If loop1 == number of sites then the right neighbour site is site 1
          IndexOfNeighbour = 1
       end if
                !---------------------------Calls the subroutine-----------------------------
       CALL Strongest_Bond( Bonds, SitePotential, Hopping, loop1, IndexOfNeighbour )
       

    end do
   
    return
  end subroutine create_AHM

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  ! Strongest_Bond routine. Inputs: Site potentials, Hopping Amplitude, adjacent site labels. Output: Bonds array
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine Strongest_Bond( Bonds, SitePotential, Hopping, LeftNbr, RightNbr )
    implicit none
                !---------------------------Output-------------------------------------------
    real, dimension(:), allocatable, intent(inout) :: Bonds

                !---------------------------Inputs-------------------------------------------
    real, dimension(:), intent(in) :: SitePotential
    real, dimension(:), intent(in) :: Hopping
    integer, intent(in) :: LeftNbr                       ! Site label of the left neighbour site
    integer, intent(in) :: RightNbr                      ! Site label of the right neighbour site

                !---------------------------Programming Tools--------------------------------
    real ABond, HBond1, HBond2                           ! ABond: Anderson bond; HBond1: Hubbard bond with right neighbour UHO
                                                         ! HBond2: Hubbard bond with the left neighbour UHO
    
    !---------------------------Calculate the three possible bond strengths------
    ABond = abs(Hopping(LeftNbr)/(SitePotential(LeftNbr) - SitePotential(RightNbr)))
    HBond1 = abs(Hopping(LeftNbr)/(SitePotential(LeftNbr) - (SitePotential(RightNbr) + uSite)))
    HBond2 = abs(Hopping(LeftNbr)/((SitePotential(LeftNbr) + uSite) - SitePotential(RightNbr)))

                !---------------------------Determine which is strongest---------------------
    if ( (ABond .ge. HBond1) .and. &
         (ABond .ge. HBond2) ) then                      ! If the noninteracting bond is strongest, then store that bond
       Bonds(LeftNbr) = ABond
    else if ( (HBond1 .ge. ABond) .and. &                ! If the bond between the LHO of the left site and the UHO of the right
         (HBond1 .ge. HBond2) ) then                     ! site is the strongest, then store that bond                         
       Bonds(LeftNbr) = HBond1
    else if ( (HBond2 .ge. ABond) .and. &                ! If the bond between the UHO of the left site and the LHO of the right
         (HBond2 .ge. HBond1) ) then                     ! site is the strongest, then store that bond
       Bonds(LeftNbr) = HBond2
    else if ( (IsNAN(ABond)) .or. (IsNAN(HBond1)) .or. (IsNAN(HBond2)) ) then
       Bonds(LeftNbr) = huge(0.0)
    else if ( (.not. ( Size(SitePotential) == 1 )) ) then
       print*, 'You belong in a museum, Error:001', SitePotential(LeftNbr), SitePotential(RightNbr)   ! Prints an error to screen in case there is an option not covered
                                                         ! Doesn't print when there is only a single site left in the system (bond
                                                         ! strength ends up being infinite)
       !print*, ABond, HBond1, HBond2
       Bonds(LeftNbr) = 0
    end if
       ! Note on above: only use 'greater than or equal to' when comparing bonds. This is because if they are equal it doesn't matter
       ! which is chosen
  end subroutine Strongest_Bond
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! CalcWeakBonds routine. Inputs: Bonds array, SitesRemoved. Output: WeakBonds array
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine CalcWeakBonds( Bonds, WeakBonds, SitesRemoved )
    implicit none
                !---------------------------Output-------------------------------------------
    integer, dimension(:), allocatable, intent(out) :: WeakBonds

                !---------------------------Inputs-------------------------------------------
    real, dimension(:), intent(in) :: Bonds
    integer, intent(in) :: SitesRemoved

                !---------------------------Programming Variables----------------------------
    integer NumWeakBonds                                 ! The number of weak bonds in the system at a given point in the RG
    integer loop1

                !---------------------------Initializing-------------------------------------
    NumWeakBonds = 0

                !---------------------------Count Weak Bonds---------------------------------
    do loop1 = 1, ( dim - SitesRemoved ) 
       if ( Bonds(loop1) .lt. bond_cutoff ) then         ! Check if each bond is below the cutoff. If yes, increment NumWeakBonds
          NumWeakBonds = NumWeakBonds + 1
       end if
    end do

    !print*, "Number of WeakBonds: ", NumWeakBonds
    
    if ( NumWeakBonds .gt. 0 ) then                      ! If there are any weak bonds then calculate the array
       allocate(WeakBonds(NumWeakBonds))
       NumWeakBonds = 0
       do loop1 = 1,( dim - SitesRemoved ) 
          if ( Bonds(loop1) .lt. bond_cutoff ) then
             NumWeakBonds = NumWeakBonds + 1
             WeakBonds(NumWeakBonds) = loop1
             if ( WeakBonds(NumWeakBonds) .le. 0 ) then 
                print*, 'You belong in a museum, Error:002', loop1
                stop
             end if
             
          end if
       end do
    else if ( NumWeakBonds .eq. 0 ) then                 ! If there are no weak bonds then set WeakBonds as follows. This cannot happen
                                                         ! under normal circumstances. Used to tell the code that there are no weakbond
       allocate(WeakBonds(1))
       WeakBonds(1) = 0
    end if
    !print*, "Weak Bonds: ",WeakBonds
  end subroutine CalcWeakBonds

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! GetDos routine. Inputs: Site Potentials, Cluster size, the labels for the left/right weak bonds of the cluster
  ! Outputs: The effecitve Hoppings between neighbour sites, contributions/weights to the density of states, SitesRemoved/Ignored
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine GetDos( SitePotential, Teff, ClusterSize, Cont, Weight, weakL, weakR, SitesRemoved, SitesIgnored )
    implicit none
    
                !---------------------------Inputs-------------------------------------------
    real, dimension(:), intent(in) :: SitePotential
    integer, intent(in) :: ClusterSize
    integer, intent(in) :: weakL, weakR

                !---------------------------Outputs------------------------------------------
    real, dimension(:), allocatable, intent(out) :: Cont
    real, dimension(:), allocatable, intent(out) :: Weight
    real, intent(out) :: Teff
    integer, intent(inout) :: SitesRemoved
    integer, intent(inout) :: SitesIgnored

                !---------------------------Programming Variables----------------------------
    real w1, w2, w3                                      ! The three possible grand potentials for a single site cluster
    integer Llabel                                       ! Must change weakL, relabel it so weakL doesn't have to have intent:inout
    !real, dimension(:), allocatable :: Hold             ! Array that holds the site potentials in the site potential distribution, only use if this is being recorded
    !real, dimension(:), allocatable :: tempWeight       ! Weight for the site potential distribution, only use this if it is being recorded
    real, dimension(:), allocatable :: E
    integer :: EndSite
    real :: tempDrop
    Llabel = weakL
    if ( ClusterSize .eq. dim+1 ) then
                !---------------------------Single Site Clusters-----------------------------
       SitesRemoved = SitesRemoved + 1                   ! Increment the number of sites removed
       

       if ( ( weakL .gt. weakR) ) then                      ! If weakL > weakR the cluster goes over the boundary, therefore weakL = weakR
          Llabel = weakR - ClusterSize                                      ! - ClusterSize
       end if
       if ( Size(SitePotential) .eq. 1 ) then
          Llabel = 0
       end if
       
       
       w1 = 0                                            ! Grand potential of unoccupied state
       w2 = SitePotential(Llabel+1) - ChemPot             !   "       "     of singly occupied state
       w3 = 2*SitePotential(Llabel+1) - 2*ChemPot + uSite !   "       "     of the doubly occupied state 

       !allocate(Hold(1))
       !allocate(tempWeight(1))
       !Hold(1) = SitePotential(Llabel+1)
       !tempWeight = 1
       
       !CALL Bin_Data(Potential, Hold, tempWeight, DroppedPot, PotMax, PotMin)
       !deallocate(tempWeight, Hold)
                 !---------------------------Determine the contributions/weights--------------
       ! Cont/Weight are allocatable because number of contributions depends on which ground state we are in.
       if ( (w1 .le. w2) .and. (w1 .le. w3) ) then       ! If w1 is the ground state
          allocate(Cont(1))
          allocate(Weight(1))

          ! Adding an electron when going from unoccupied state to singly occupied
          ! Convention for adding an electron is "Final grand potential" - "Initial Grand Potential"
          Cont = w2 - w1
          Weight = 1
          
       else if ( (w2 .le. w1) .and. (w2 .le. w3) ) then
          allocate(Cont(2))
          allocate(Weight(2))
  
          ! Removing an electron when going from singly occupied state to unoccupied
          ! Convention for removing an electron is "Initial Grand Potential" - "Final grand potential" 
          Cont(1) = w2 - w1
          Weight(1) = 0.5

          ! Adding an electron when going from singly occupied state to doubly occupied
          ! Convention for adding an electron is "Final grand potential" - "Initial Grand Potential"
          Cont(2) = w3 - w2
          Weight(2) = 0.5
          
       else if ( (w3 .le. w1) .and. (w3 .le. w2) ) then
          allocate(Cont(1))
          allocate(Weight(1))

          ! Removing an electron when going from doubly occupied state to singly occupied
          ! Convention for removing an electron is "Initial Grand Potential" - "Final grand potential"
          Cont = w3 - w2
          Weight = 1
          
       else
          print*, 'Error in DOS for cluster size 1'      ! Can occur if w1, w2 or w3 are NaN
       end if
       
       if ( DOSCl .ne. 0 ) then
          CALL Bin_Data( Dos(:,1+1), Cont, Weight, &                   ! Bin the contributions based on their weights using the Bin_Data routine
                  DroppedDos,TotalDos, DoSMax, DoSMin )
       end if
                !---------------------------Assign Effective Hopping-------------------------
       Teff = 0
       
                !---------------------------Larger than single site clusters-----------------
    else if ( ( ClusterSize .ge. 1 ) .and. ( ClusterSize .le. ClusterMax ) ) then
       allocate(Cont(2*ClusterSize*(4**ClusterSize)))
       allocate(Weight(2*ClusterSize*(4**ClusterSize)))
	
       Cont = 0.0
       Weight = 0.0
       allocate(E(ClusterSize))
       if ( weakL .gt. weakR) then
          EndSite = (dim - sitesremoved) - weakL
          E(1:EndSite) = SitePotential((weakL+1):(dim - sitesremoved))
          E((EndSite+1):ClusterSize) = SitePotential(1:weakR)
             
       else
          E = SitePotential((weakL+1):weakR)
       end if
       
       CALL DiagCluster(ClusterSize,4**ClusterSize,E,Cont,Weight,ChemPot,uSite,hop)
       if ( ClusterSize .le. DOSCl ) then
          CALL Bin_Data(Dos(:,1+ClusterSize),Cont,ClusterSize*Weight, tempDrop, tempDrop, DoSMax, DoSMin)
       end if
       Teff = 0
       SitesRemoved = SitesRemoved + ClusterSize
    else if ( ClusterSize .gt. clusterMax ) then                
       SitesRemoved = SitesRemoved + ClusterSize         ! Increment sites removed, even tho currently no RG is performed, only removed
       SitesIgnored = SitesIgnored + ClusterSize         ! Increment the number of sites ignored, happens when cluster is too large
       allocate(Cont(1))
       allocate(Weight(1))

       ! Here I set the weight to zero, this is how the site is ignored. Cont is irrelevant
       Cont(1) = 0
       Weight(1) = 0

       ! Effective Hopping still 0
       Teff = 0
       
    end if
    
  end subroutine GetDos

  subroutine MinCouplingLoc( WeakBonds, Bonds, ClusterStage, SitesRemoved, weakL, weakR)
    implicit none
   
    integer, dimension(:), intent(in) :: WeakBonds
    real, dimension(:), intent(in) :: Bonds
    integer, dimension(:), allocatable :: Temp
    real, dimension(:,:), allocatable :: Order
    integer, intent(out) :: weakL, weakR
    integer, intent(in) :: ClusterStage
    integer, intent(in) :: SitesRemoved
    

    integer numCluster
    integer nextWeak
    integer loop1, loop2, row
    integer ClusterSize

    allocate(Temp(Size(WeakBonds)))
    Temp = 0

    do loop1 = 1,size(WeakBonds)-1
       Temp(loop1 + 1) = WeakBonds(loop1)
    end do

    Temp(1) = WeakBonds(Size(WeakBonds)) - (dim - SitesRemoved)

    numCluster = Count(abs( Temp - WeakBonds ) .eq. ClusterStage)

    if (numCluster .gt. 0) then
       
       allocate(Order(numCluster,3))
       
       loop2 = 1
       do loop1 = 1,Size(WeakBonds) 
          weakL = WeakBonds(Loop1)                                  ! Currently cluster has weakL as its left neighbour
          
          nextWeak = Loop1 + 1                                      ! Next bond label in WeakBonds creates a cluster with weakL
          if ( (Loop1 + 1) .gt. size(WeakBonds) ) then              ! If weakL is the last weakbond in the system then search the beginning of the WeakBonds array
             nextWeak = 1                                           
             weakR = WeakBonds(nextWeak)                            ! weakR is the first weak bond in the system
             ClusterSize = weakR - ( weakL - (dim - SitesRemoved) ) ! Have to calculate the cluster size in this way due to boundary conditions
          else                                                      ! For all other cases, the weakR label is just the next element in the weakbonds array
             weakR = WeakBonds(nextWeak)
             ClusterSize = weakR - weakL
          end if
          
          
          if ( ClusterSize .eq. ClusterStage ) then
             Order(loop2, 1) = Bonds(weakL)**2 + Bonds(weakR)**2
             Order(loop2, 2) = weakL
             Order(loop2, 3) = weakR
             loop2 = loop2 + 1
          end if
          
          
       end do

       row = minloc( Order( 1:numCluster, 1 ), dim=1 ) 
       

       do loop1 = 1,numCluster
          !print*, Order(loop1,:)
       end do
       

       weakL = Order(row,2); weakR = Order(row,3)

       deallocate(Order,Temp)
    else if ( numCluster .eq. 0 ) then
       weakL = 0; weakR = 0
    else
       print*, "Error in MinCouplingLoc", numCluster
       stop
    end if
    
    
  end subroutine MinCouplingLoc
  

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Full_DoS routine. Inputs: Site Potentials, Hopping Amplitudes, Bond Strengths, Weak bonds array, Density of states parameters
  ! Outputs: Density of states for the system, number of contributions dropped from the density of states, Number of sites not
  ! renormalized
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine Full_DoS( SitePotential, Hopping, Bonds, WeakBonds, SitesIgnored )
    implicit none

                !---------------------------Inputs-------------------------------------------
    real, dimension(:), allocatable, intent(inout) :: SitePotential    
    real, dimension(:), allocatable, intent(inout) :: Hopping          
    real, dimension(:), allocatable, intent(inout) :: Bonds            
    integer, dimension(:), allocatable, intent(inout) :: WeakBonds
    
                !---------------------------Outputs------------------------------------------                                                    
    integer, intent(inout) :: SitesIgnored
    real, dimension(:), allocatable :: Cont              ! Array that will hold contributions to the DoS for a given cluster.
                                                         ! Allocatable because # of contributions not fixed
    real, dimension(:), allocatable :: Weight            ! The weights of the contributions. Same length as 'Cont', sums to 1.
    

                !---------------------------Programming Variables----------------------------
    integer SitesRemoved                                 ! Keeps track of the number of sites removed at each stage of the RG. Done to
                                                         ! change lengths of loops over the system size
    integer weakL                                        ! The label for the weak bond on the "left" of the current cluster
    integer weakR                                        ! The label for the weak bond on the "left" of the current cluster
    real Teff                                            ! Effecitve hopping amplitude between neighbouring sites after a cluster is
                                                         ! removed between them
    integer ClusterSize                                  ! Size of the current cluster
    integer IndexOfNeighbour, Loop1, Loop2               ! Loop integers
    !integer, dimension(dim) :: ClusterCount
                !---------------------------Initializing-------------------------------------
    
    SitesRemoved = 0
    weakL = 0
    weakR = 0
    
    !print*, WeakBonds
    !print*, "haha"
    !allocate(ClusterCount(dim))
    !Call AnalyzeClusters(WeakBonds,ClusterCount)
    !print*, DoSMax, DoSMin
    Loop1 = 1                                                          ! Current Cluster size, starting from 1 going to system size
    Loop2 = 0
    
    do while ( (size(WeakBonds) .gt. 1) .and. (WeakBonds(1) .ne. 0)  ) ! Stop when WeakBonds is set to have 1 element or have the first
                                                                       ! element be 0. Either only one cluster left or no weak bonds.

       !print*, "Number of WeakBonds: ", Size(WeakBonds)
                !---------------------------Determine Cluster Size and the bond labels-------
       CALL MinCouplingLoc( WeakBonds, Bonds, Loop1, SitesRemoved, weakL, weakR)
       
       if ( weakL .gt. weakR ) then                              ! If weakL is the last weakbond in the system then search the beginning of the WeakBonds array
          ClusterSize = weakR - ( weakL - (dim - SitesRemoved) ) ! Have to calculate the cluster size in this way due to boundary conditions
       else if ( weakL .lt. weakR) then                          ! For all other cases, the weakR label is just the next element in the weakbonds array
          ClusterSize = weakR - weakL
       else if ( (weakL .eq. 0) .and. (weakR .eq. 0) ) then
          ClusterSize = 0
          Loop1 = Loop1 + 1
       else
          Print*, "weakL = weakR, but not 0 in Full_DoS"
          print*, WeakBonds
          print*, weakL, weakR
          stop
       end if

                !---------------------------Perform RG is cluster is of smallest size remaining-----------
                ! Search for smallest clusters first, labeled as Loop1, Loop1 increases when clusters of a certain size are eliminated
       if ( ClusterSize .eq. Loop1 ) then                                 ! If the current cluster is small enough then remove it                      
          Loop2 = 0

                !---------------------------Calculate DoS for cluster and bin the contributions-----------
          CALL GetDoS( SitePotential, Teff, ClusterSize &                 ! Call the Get_DoS routine to calculate the contributions and their weights to the DoS
               , Cont, Weight, weakL, weakR, SitesRemoved, SitesIgnored )  

          CALL Bin_Data( Dos(:,1), Cont, ClusterSize*Weight, &                   ! Bin the contributions based on their weights using the Bin_Data routine
               DroppedDos,TotalDos, DoSMax, DoSMin )
          
          !print*, DoS
          
                !---------------------------Removing the cluster from the relevant data-------------------
          !print*, Loop1, Loop2
          CALL resize_array( SitePotential, ClusterSize, weakL + 1 )          ! Remove the sites in the cluster used above which has ClusterSize sites and is bounded on the left by weakL
          !print*, "i see"
          CALL resize_array( Hopping, ClusterSize, weakL + 1 )                ! Same process as above, but instead for the Hopping amplitude
          !print*, "ooo ahh"   

          if ( weakL .gt. Size(Hopping) ) then                            ! weakL is used for the last cluster. If the cluster in question is around the boundary and weakL referenced
             Hopping(Size(Hopping)) = Teff                            ! a bond that has a label larger than the new system size, the effective bond between sites neighbouring the
          else                                                            ! removed cluster has a label of weakL-ClusterSize
             Hopping(weakL) = Teff                                        ! Otherwise the new bond has the same value as weakL
          end if
          !print*, "hehe"
          CALL resize_array( Bonds, ClusterSize, weakL + 1 )                  ! Resize the Bonds array
          !print*, "Cough"

         
          !do loop5 = 1,( dim - SitesRemoved )                             ! Use this loop to recalculate the Bonds array. This must be done before of Teff changing that single bond
             IndexOfNeighbour = weakL+1                                   ! In early stages of this code, teff will always be zero so could just set that same element to 0 in the bonds
             ! array. But in the future teff will be nonzero.
             if ( weakL .ge. ( dim - SitesRemoved) ) then                  ! Could also consider only recalculating the single bond but this current set up works and this is not the
                IndexOfNeighbour = 1                                      ! time consuming part of the code.
                weakL = dim - SitesRemoved
             end if                                                       ! Recalculating the whole thing as allowed me to catch an error in the past.
             
             CALL Strongest_Bond( Bonds, SitePotential, Hopping, weakL, & ! Recalculate the loop5'th strongest bond.
                  IndexOfNeighbour )
             
          !end do
          
          
          
          deallocate(WeakBonds)                                           ! deallocate the weakbonds array
          
          CALL CalcWeakBonds(Bonds,WeakBonds,SitesRemoved)                ! recalculate the weakbonds array, using the newly calculated Bonds array
          deallocate(Cont,Weight)
          
          
       else if ( ClusterSize .gt. Loop1 ) then
          print*, Loop1, weakL, weakR
          stop
       else if ( (ClusterSize .lt. Loop1) .and. (ClusterSize .ne. 0) ) then
          print*, "Small Clusters Found Late"                             ! This implicitly sets weakL to weakR and weakR to weakR+1
          print*, WeakBonds
          stop
       else if ( ClusterSize .eq. 0 ) then
          !If the cluster size is 0 then there are no more clusters of that size in the RG
       else
          print*, "Error in Full_Dos"
          stop
       end if

       
       
       
    end do
    
    ClusterSize = dim - SitesRemoved
    weakL = WeakBonds(1); weakR = WeakBonds(1)
 
    CALL GetDoS( SitePotential, Teff, ClusterSize &
         , Cont, Weight, weakL, weakR, SitesRemoved, SitesIgnored ) 
    CALL Bin_Data( DoS(:,1), Cont, ClusterSize*Weight, DroppedDos, TotalDos, DoSMax, DoSMin )
        
    
  end subroutine Full_DoS

   subroutine AnalyzeClusters(WeakBonds, ClusterCount)
    implicit none
    integer loop1, loop2!loop integer
    integer, dimension(:), allocatable, intent(in) :: WeakBonds !An array where each element is 0 or 1. 1 means weak bond, 0 means strong bond
    integer, intent(out) :: ClusterCount(dim) !Counts the number of clusters recorded for each cluster size
    integer ClusterSize !contains the size of the current cluster
    integer LbondLabel, RbondLabel !The number that labels the left and right weak bonds that neighbour each cluster
   
    ClusterCount(1:size(ClusterCount)) = 0
    
    if ( WeakBonds(1) .ne. 0 ) then
       do loop1 = 1,SIZE(WeakBonds) !cycle over the clusters in the system, number of weak bonds == number of clusters
          
          LbondLabel = WeakBonds(loop1) !Extract the label for the left bond for the current cluster
          if (loop1 .eq. SIZE(WeakBonds)) then !If the current cluster is the last cluster, then the second weak bond in this cluster is also the first weak bond in the system
             loop2 = 1
             LbondLabel = LbondLabel - dim !Also redefine the left label as an integer <= 0, the last bond on the chain (bond between last and first sites) is also the 0th bond in the cluster.
          else if (loop1 .lt. SIZE(WeakBonds)) then !If the current bond is not the last bond, then the label for the right neighbour bond is one element larger than the left neighbour bond in WeakBonds
             loop2 = loop1 + 1
          else !In case something funky happens. Can't see how but better safe than sorry
             print*, 'You belong in a museum: error 003' 
          end if
          RbondLabel = WeakBonds(loop2) !Extract the label for the right neighbour bond for the current cluster
          
          ClusterSize = abs(RbondLabel - LbondLabel) !this calculates the size of the cluster. If a cluster is between bond 3 and 4, then the cluster size is 4-3=1
          if (ClusterSize .eq. 0) then
             print*, 'Rbond', RbondLabel
             print*, 'Lbond', LbondLabel
             print*, 'Cluster Size', ClusterSize
             do loop2 = 1,5
                print*, 'WeakBonds', WeakBonds(loop2)
             end do
             print*, 'length of WeakBonds', size(WeakBonds)
             stop
          end if
          
          ClusterCount(ClusterSize) = ClusterCount(ClusterSize) + ClusterSize !Increase the number of clusters of size 'ClusterSize' by an amount equal to 'ClusterSize'
       end do
       
    else if ( WeakBonds(1) .eq. 0 ) then
       ClusterCount(Size(ClusterCount)) = ClusterCount(Size(ClusterCount)) + Size(ClusterCount)
    end if
       
    
  end subroutine AnalyzeClusters

  subroutine Bin_Data( HistoData, Data, Weights, Dropped, Total, Max, Min )
    implicit none
    real, dimension(:), allocatable,  intent(in) :: Data
    real, dimension(:), allocatable, intent(in) :: Weights
    real, intent(inout) :: HistoData(bins)
    real, intent(in) :: Max
    real, intent(in) :: Min
    real, intent(inout) :: Dropped, Total
    integer Loop1, Loop2
 
    !print*, "Max: ", Max
    
    do Loop1 = 1,size(Data)
       if ( Weights(Loop1) .ne. 0 ) then
          Total = Total + Weights(Loop1)
          if ( ( Data(Loop1) .lt. Min ) .or. ( Data(Loop1) .gt. Max ) ) then
             Dropped = Dropped + Weights(Loop1)
          else
          
             Loop2 = FLOOR(bins*(Data(Loop1) - Min)/(Max - Min) + 1)

             HistoData(Loop2) = HistoData(Loop2) + Weights(Loop1)

          end if
       end if
    end do
    


  end subroutine Bin_Data
      
!  
! This subroutine, resize_array, takes in an array - given the name 'array' - and an integer 'n' and returns the same array, but with length 'SIZE(array) + n'. This only works for 1 dimensional arrays    
! 
  
  subroutine resize_array(array, numRemove, Start)
    real, dimension(:), allocatable :: tmp_arr
    real, dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: numRemove                           ! Number of elements to remove from 'array'
    integer, intent(in) :: Start                               ! The element in which to start removing the array
    integer i, j                                               ! Looping integer
    

    allocate(tmp_arr(size(array) - numRemove))

    tmp_arr(:) = 0

    j=1

    if ( (Start + numRemove) .le. size(array) ) then


       do i = 1, Start - 1

          tmp_arr(i) = array(i)
          j = j + 1

       end do


       do i = (Start + numRemove), Size(array)

          tmp_arr(j) = array(i)
          j = j + 1

       end do



    else if ( (Start + numRemove) .gt. size(array) ) then

       do i = (Start + numRemove - size(array)) , Start-1

          tmp_arr(j) = array(i)
          j = j + 1

       end do


    end if

    deallocate(array)
    allocate(array(size(tmp_arr)))
    
    array = tmp_arr
    
   
  end subroutine resize_array
  
  character(len=20) function str(k)                  ! function to change an integer to a string, used in FileNamer subroutine
     integer, intent(in) :: k

     write (str, *) k
     str = adjustl(str)
   end function str

   character(len=25) function FileNamer(Name)                                ! Increments filename labels Blah023.dat label:023
     character(len=25), intent(in) :: Name                                   ! Initial filename. Would be 'Blah' from above example
     integer num                                                             ! Keeps track of the label increment as this continues to check for name existence
     logical Exist                                                           ! Used to store if a filename exists

     Exist = .true.
     num = 0
     FileNamer = ""//trim(Name)//"000.text"                                  ! Takes Blah and puts it into the format 'Blah000.text'
  
     do while (Exist)                                                        ! Continue to change filename if the current one exists
        inquire( file = FileNamer, exist = Exist )                           ! Checks the current folder (where this .f90 file is found) if the current filename exists
        if ( Exist ) then                                                    ! If the file name exists
           num = num + 1                                                     ! Increment the label

           if ( len(trim(str(num))) .eq. 1 ) then                            ! Labels can be from 000 to 999. To preserve filename length the leading zeros are required
              FileNamer = ""//trim(Name)//"00"//trim(str(num))//".text"      ! This if statement (and the one below) check the length of the label (1 has length 1, 23 has length 2, 450 has length 3)
           else if ( len(trim(str(num))) .eq. 2) then                        ! This determines the number of leading zeros to add to the filename label
              FileNamer = ""//trim(Name)//"0"//trim(str(num))//".text"
           else
              FileNamer = ""//trim(Name)//trim(str(num))//".text"
           end if
           
        else if ( .not. Exist ) then                                         ! If it does not exist. Note: Didn't increment 'num'
           if ( len(trim(str(num))) .eq. 1 ) then
              FileNamer = ""//trim(Name)//"00"//trim(str(num))//".text"
           else if ( len(trim(str(num))) .eq. 2) then
              FileNamer = ""//trim(Name)//"0"//trim(str(num))//".text"
           else
              FileNamer = ""//trim(Name)//trim(str(num))//".text"
           end if
           Exist = .false.
        else
           print*, "You belong in a museum. Error: FileNameErr"              ! Outputs error if something strange happens. Can't imagine what
        end if
     end do
   end function FileNamer

   

end module SHARED_MODULE

              
program AndersonHubbard
     use Diag
     use SHARED_MODULE                                   ! Makes the data in the SHARED_MODULE routine available
     implicit none
     
               !----------------------------System Parameters--------------------------------
     real, dimension(:), allocatable :: SitePotential    ! Contains the site potentials of all the sites, element i == site potential i
     real, dimension(:), allocatable :: Hopping          ! Contains the hopping amplitudes between sites, element i == hopping
                                                         ! amplitude between site i and i+1
     real, dimension(:), allocatable :: Bonds            ! Contains the strongest bond between sites, element i == bond strength
                                                         ! between site i and i+1
     
               !----------------------------Important Quantities-----------------------------
     integer, dimension(:), allocatable :: WeakBonds     ! An array of length X where each element corresponds to an element in Bonds
     ! that is below the bond cutoff where X is the number of weak bonds.
     
     
     

               !----------------------------DoS Stuff----------------------------------------
                                      ! Contains the right edge of all the bins in the density of states
     integer SitesRemoved
     integer SitesIgnored                                ! Keeps track of all the sites in clusters larger than this code can handle
     

               !----------------------------Coding Tools-------------------------------------
     real BinWidth                                       ! Width of the DoS bins for the purposes of normalization
     integer Loop1,i                                       ! Loop Integer
     character(len=25) filename                          ! Base file name to be put into FileNamer
     real :: start, finish, TIME                         ! Used to keep track of CPU time and prints to screen at the end

     
               !----------------------------Initializing-------------------------------------
     DoS = 0
     
     DroppedDos = 0
     TotalDos = 0
     DroppedPot = 0
     
     do loop1 = 1,bins
     end do
  

     !----------------------------Opening files------------------------------------
     do i = 1,DOSCl+1
	if ( DOSCl .eq. 0 ) EXIT  
        filename = "DoS"//trim(str(i))//""
        filename = FileNamer(filename)
        open (unit = 100+i, file = filename ,status ='unknown' )
        
        !Comments for the file above, says whats in the file and which parameters go with the data
        write(100+i,*) "#This file contains the density of states for the ensemble. First column is", & 
             "# right edge of bin, fraction of contributions per bin"
        write(100+i,*) "#Size of clusters = #", i
        write(100+i,*) "#Disorder Strength = ", DELTA
        write(100+1,*) "#Bond Cutoff = ", bond_cutoff
        write(100+i,*) "#Interaction strength = ", uSite
        write(100+i,*) "#Dimensions = ", dim
        write(100+i,*) "#Number of systems = ", systemn
       
     end do
  
     
     filename = "DoS"
     filename = FileNamer(filename)
     open (unit = 100, file = filename ,status ='unknown' )
     
     !Comments for the file above, says whats in the file and which parameters go with the data
     write(100,*) "#This file contains the density of states for the ensemble. First column is", & 
          "# right edge of bin,second is fraction of contributions per bin"
     write(100,*) "#This is the DOS with 0 hopping between every site ", &
           "#(aka an ensemble of two site clusters)", &
          "#This uses the Diag module to diagonalize the clusters"
     write(100,*) "#Disorder Strength = ", DELTA
     write(100,*) "#Bond Cutoff = ", bond_cutoff
     write(100,*) "#Interaction strength = ", uSite
     write(100,*) "#Dimensions = ", dim
     write(100,*) "#Number of systems = ", systemn
  
!     filename = "Potential"
!     filename = FileNamer(filename)
!     open (unit = 200, file = filename ,status ='unknown' )
     
     !Comments for the file above, says whats in the file and which parameters go with the data
!     write(200,*) "#This file contains the distribution of site potentials found in single site clusters"
!     write(200,*) "#Disorder Strength = ", DELTA
!     write(200,*) "#Bond Cutoff = ", bond_cutoff
!     write(200,*) "#Interaction strength = ", uSite
!     write(200,*) "#Dimensions = ", dim
!     write(200,*) "#Number of systems = ", systemn
     call init_random_seed()
     
     SitesIgnored = 0
     SitesRemoved = 0
                !---------------------------RG begins----------------------------------------
     CALL CPU_TIME(start)
     do Loop1 = 1,systemn
        if ( mod(real(Loop1),0.1*systemn) == 0.0 ) then
           !print*, int(Loop1/real(systemn)*100), "%"
        end if
        
        

                !---------------------------Create System------------------------------------
        call create_AHM( SitePotential, Hopping, Bonds )
        call CalcWeakBonds( Bonds, WeakBonds, SitesRemoved )
        if ( WeakBonds(1) .eq. 0 ) then
           print*, "All Bonds Strong"
           if ( systemn .gt. ClusterMax) CYCLE
        end if
                !---------------------------Calculate DoS------------------------------------
        call Full_DoS( SitePotential, Hopping, Bonds, WeakBonds, SitesIgnored )
     end do
     CALL CPU_TIME(finish)
     BinWidth = BinsDoS(2) - BinsDoS(1)                ! Area of each bin (bins have equal width)
     do i=0,DOSCl
        if ( (i .eq. 0) .and. (i .eq. DOSCl) .and. (Loop1 .eq. bins) ) CYCLE
        do Loop1=1,bins
           write(100+i,400) BinsDoS(Loop1) , DoS(Loop1,i+1)/( Sum(DoS(:,i+1)) * BinWidth )
           print*, BinsDoS(Loop1)
        end do
     end do
	print*, BinsDoS
     
     BinWidth = BinsPot(2) - BinsPot(1)
     !print*, Potential
     !print*, DroppedPot
     !do Loop1=1,bins
     !   write(200,400) BinsPot(Loop1) , Potential(Loop1)/( Sum(Potential) * BinWidth )
     !end do
     
     TIME = finish - start
     print*, "Sites in large clusters: ", SitesIgnored, " out of", n_d
     !This loop checks area under the DOS curve. Should equal 1.
     !do i = 1,bins
     !   DroppedPot = DroppedPot + DOS(i,1)/Sum(DOS(:,1))
     !end do
     print*, DroppedPot
     print*, "The sum of the weights of the DOS contributions that were dropped ", DroppedDos, " out of", TotalDos
     print '("Time of New= ", I1.1," days,", I2.1," hours,", I2.1," minutes, and ", I2.1," seconds.")', &
          int(TIME/real(86400)), mod(int(TIME/real(3600)),24), mod(int(TIME/real(60)),60), mod(int(TIME),60)
  
!300  format(  i10, g12.5)                               ! Format for ( integer , real ) file
400  format( g12.5, g12.5)                               ! Format for ( real , real ) file
	 
   end program AndersonHubbard
		
   ! Last Edited: 180208
   ! Next Step: Properly comment things, explain orders of w1/w2/w3 when calculating the contribution, set weights to zero for some 
   ! each of the 4 possible contributions to show where they appear on the DoS. Write the outline for stage 2 write up
   !
