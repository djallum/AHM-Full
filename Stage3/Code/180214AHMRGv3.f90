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
!module stuff
!  implicit none
!  save

!contains
!  subroutine saysdf()
!    implicit none
!    print*, "It worked"
!  end subroutine saysdf
!end module stuff

module SHARED_MODULE
  use stuff
  IMPLICIT NONE
  SAVE
  integer, parameter :: systemn = 100 !number of systems
  integer, parameter :: dim = 1000 !system size 
  REAL, parameter :: DELTA = 12 !disorder strength
  REAL, parameter :: hop = -1 !hopping amplitude
  REAL, parameter :: bond_cutoff = 100 !maximum bond strength allowed for a site or cluster to be removed
  REAL, parameter :: uSite = 8 !Interaction strength
  REAL, parameter :: ChemPot = uSite/2 !Chemical Potential
  integer, parameter :: n_d = systemn*dim !number of systems*number of sites in a single system
  !-------------DoS binning parameters---------------------
  integer, parameter :: bins = 100 !number of bins in the density of states
  integer, parameter :: BinMin = -10 !the minimum energy value of the smallest bin in the DoS
  !-------------Run This part? -------------------------
  logical, parameter :: DoClusterCount = .false. !if = 1,do cluster count/distribution stuff. If = 0, don't
  logical, parameter :: DoBondStrength = .true.
  logical, parameter :: DoEnergyDiff = .false.
  
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
       Hopping(loop1) = hop                              ! Sets the hopping amplitude equal to the parameter 'hop'
       SitePotential(loop1) = (0.5 - RandomNumber)*DELTA ! Shifts 'RandomNumber' from interval [0,1] to [-DELTA,DELTA]
    end do

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
    real, dimension(:), allocatable, intent(in) :: SitePotential
    real, dimension(:), allocatable, intent(in) :: Hopping
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
    else if ( .not. ( Size(SitePotential) == 1 ) ) then
       print*, 'You belong in a museum, Error:001', Size(SitePotential)      ! Prints an error to screen in case there is an option not covered
                                                         ! Doesn't print when there is only a single site left in the system (bond
                                                         ! strength ends up being infinite)
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
    real, dimension(:), allocatable, intent(in) :: Bonds
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
    real, dimension(:), allocatable, intent(in) :: SitePotential
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
    Llabel = weakL
                !---------------------------2-Site variables---------------------------------
    
    !real, dimension(16) :: grand_potential               ! grand potentials (eigenenergies - mu*number electrons), currently only works for clusters of size 2
    


    if ( ClusterSize .eq. 1 ) then
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

                !---------------------------Assign Effective Hopping-------------------------
       Teff = 0
       
                !---------------------------Larger than single site clusters-----------------
    !else if ( ClusterSize .eq. 2 ) then
    !  SitesRemoved = SitesRemoved + 2                   ! Increment the number of sites removed
       

    !   if ( ( weakL .gt. weakR) ) then                      ! If weakL > weakR the cluster goes over the boundary, therefore weakL = weakR- ClusterSize
    !      Llabel = weakR - ClusterSize
    !      if ( Llabel .lt. 1 ) then
    !         Llabel = dim - SitesRemoved - Llabel
    !      end if
    !   end if
    
     !  SitesRemoved = SitesRemoved + 2                   ! Increment the number of sites removed
       
                !---------------------------Larger than single site clusters-----------------
    else if ( ClusterSize .gt. 1 ) then                
       SitesRemoved = SitesRemoved + ClusterSize         ! Increment sites removed, even tho currently no RG is performed, only removed
       SitesIgnored = SitesIgnored + ClusterSize         ! Increment the number of sites ignored, happens when cluster is too large
       allocate(Cont(1))
       allocate(Weight(1))
       call saysdf()

       ! Here I set the weight to zero, this is how the site is ignored. Cont is irrelevant
       Cont(1) = 0
       Weight(1) = 0

       ! Effective Hopping still 0
       Teff = 0
       
    end if
    
  end subroutine GetDos

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! hamiltonian_2 routine calculates the grand potentials for a 2 site system
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Full_DoS routine. Inputs: Site Potentials, Hopping Amplitudes, Bond Strengths, Weak bonds array, Density of states parameters
  ! Outputs: Density of states for the system, number of contributions dropped from the density of states, Number of sites not
  ! renormalized
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine Full_DoS( SitePotential, Hopping, Bonds, WeakBonds, DoS, DoSMax, DoSMin, BinValue, DroppedDos, SitesIgnored)
    implicit none

                !---------------------------Inputs-------------------------------------------
    real, dimension(:), allocatable, intent(inout) :: SitePotential    
    real, dimension(:), allocatable, intent(inout) :: Hopping          
    real, dimension(:), allocatable, intent(inout) :: Bonds            
    integer, dimension(:), allocatable, intent(inout) :: WeakBonds
    integer, dimension(:), allocatable :: Test
    real, dimension(:), allocatable :: Test1
    ! DoS Parameters
    real, intent(in) :: DoSMax, DoSMin
    real, intent(in) :: BinValue(bins)
    
                !---------------------------Outputs------------------------------------------    
    real, intent(out) :: DoS(bins)                                                 
    integer, intent(inout) :: DroppedDos
    integer, intent(inout) :: SitesIgnored
    real, dimension(:), allocatable :: Cont              ! Array that will hold contributions to the DoS for a given cluster.
                                                         ! Allocatable because # of contributions not fixed
    real, dimension(:), allocatable :: Weight            ! The weights of the contributions. Same length as 'Cont', sums to 1.

                !---------------------------Programming Variables----------------------------
    integer SitesRemoved                                 ! Keeps track of the number of sites removed at each stage of the RG. Done to
                                                         ! change lengths of loops over the system size
    integer nextWeak                                     ! Used to calculate weakR, labels the "next" weak bond in the weak bond array
    integer weakL                                        ! The label for the weak bond on the "left" of the current cluster
    integer weakR                                        ! The label for the weak bond on the "left" of the current cluster
    real Teff                                            ! Effecitve hopping amplitude between neighbouring sites after a cluster is
                                                         ! removed between them
    logical BinFound                                     ! Keeps track of whether a bin is found for a contribution
    integer ClusterSize, i                                  ! Size of the current cluster
    integer IndexOfNeighbour, Loop1, Loop2, Loop3, Loop4, Loop5 ! Loop integers
    
                !---------------------------Initializing-------------------------------------
    
    SitesRemoved = 0
    DoS(1:bins) = 0
    weakL = 0
    weakR = 0
    
    
    Loop1 = 1                                                          ! Current Cluster size, starting from 1 going to system size
    Loop2 = 1                                                          ! Used to loop through weakbonds array
    do while ( (size(WeakBonds) .gt. 1) .and. (WeakBonds(1) .ne. 0)  ) ! Stop when WeakBonds is set to have 1 element or have the first
                                                                       ! element be 0. Either only one cluster left or no weak bonds.

      
                !---------------------------Determine Cluster Size and the bond labels-------
       weakL = WeakBonds(Loop2)                                  ! Currently cluster has weakL as its left neighbour

       nextWeak = Loop2 + 1                                      ! Next bond label in WeakBonds creates a cluster with weakL
       if ( (Loop2 + 1) .gt. size(WeakBonds) ) then              ! If weakL is the last weakbond in the system then search the beginning of the WeakBonds array
          nextWeak = 1                                           
          weakR = WeakBonds(nextWeak)                            ! weakR is the first weak bond in the system
          ClusterSize = weakR - ( weakL - (dim - SitesRemoved) ) ! Have to calculate the cluster size in this way due to boundary conditions
       else                                                      ! For all other cases, the weakR label is just the next element in the weakbonds array
          weakR = WeakBonds(nextWeak)
          ClusterSize = weakR - weakL
       end if
       
                !---------------------------Perform RG is cluster is of smallest size remaining-----------
                ! Search for smallest clusters first, labeled as Loop1, Loop1 increases when clusters of a certain size are eliminated
       if ( ClusterSize .eq. Loop1 ) then                                 ! If the current cluster is small enough then remove it                      

          
                !---------------------------Calculate DoS for cluster and bin the contributions-----------
          CALL GetDoS( SitePotential, Teff, ClusterSize &                 ! Call the Get_DoS routine to calculate the contributions and their weights to the DoS
               , Cont, Weight, weakL, weakR, SitesRemoved, SitesIgnored )  
          
          CALL Bin_Data( DoS, Cont, Weight, BinValue, &                   ! Bin the contributions based on their weights using the Bin_Data routine
               DroppedDos, DoSMin, DoSMax )
         
          
                !---------------------------Removing the cluster from the relevant data-------------------

          CALL resize_array( SitePotential, ClusterSize, weakL )          ! Remove the sites in the cluster used above which has ClusterSize sites and is bounded on the left by weakL
          CALL resize_array( Hopping, ClusterSize, weakL )                ! Same process as above, but instead for the Hopping amplitude
             

          if ( weakL .gt. Size(Hopping) ) then                            ! weakL is used for the last cluster. If the cluster in question is around the boundary and weakL referenced
             Hopping(weakL-ClusterSize) = Teff                            ! a bond that has a label larger than the new system size, the effective bond between sites neighbouring the
          else                                                            ! removed cluster has a label of weakL-ClusterSize
             Hopping(weakL) = Teff                                        ! Otherwise the new bond has the same value as weakL
          end if

          CALL resize_array( Bonds, ClusterSize, weakL )                  ! Resize the Bonds array

          
          do loop5 = 1,( dim - SitesRemoved )                             ! Use this loop to recalculate the Bonds array. This must be done before of Teff changing that single bond
             IndexOfNeighbour = loop5+1                                   ! In early stages of this code, teff will always be zero so could just set that same element to 0 in the bonds
                                                                          ! array. But in the future teff will be nonzero.
             if (loop5 .eq. ( dim - SitesRemoved) ) then                  ! Could also consider only recalculating the single bond but this current set up works and this is not the
                IndexOfNeighbour = 1                                      ! time consuming part of the code.
             end if                                                       ! Recalculating the whole thing as allowed me to catch an error in the past.

             CALL Strongest_Bond( Bonds, SitePotential, Hopping, loop5, & ! Recalculate the loop5'th strongest bond.
                  IndexOfNeighbour )
             
          end do
          
          
          deallocate(WeakBonds)                                           ! deallocate the weakbonds array
          
          CALL CalcWeakBonds(Bonds,WeakBonds,SitesRemoved)                ! recalculate the weakbonds array, using the newly calculated Bonds array
          

          
       else
          Loop2 = Loop2 + 1                                               ! If the current cluster bordered by weakL and weakR is larger than loop1, we increase loop2 by one
       end if                                                             ! This implicitly sets weakL to weakR and weakR to weakR+1
       
          
          
       if ( Loop2 .gt. size(WeakBonds) ) then                             ! If we've exceeded the size of the weakbonds array with Loop2 that means that all clusters of size loop1
                                                                          ! have been removed. Restart at the beginning of the weakbonds array and start searching for larger clusters
          Loop2 = 1                                                       ! Restart at beginning of weakbond array 
          Loop1 = Loop1 + 1                                               ! Search for larger clusters
          
       else if ( Size(WeakBonds) .eq. 2 ) then                            ! If there are two clusters left, loop1 can 
          Loop2 = 1
          if ( ( Loop1 .lt. ClusterSize) ) then !.or. &
               !( Loop1 .gt. ClusterSize) ) then
             Loop1 = ClusterSize
          end if
       else if ( Loop1 .gt. dim ) then
          print*, "Searching for clusters that are too large"
          print*, WeakBonds
          stop
       end if
       
       
    end do
    
    ClusterSize = dim - SitesRemoved
    CALL GetDoS( SitePotential, Teff, ClusterSize &
         , Cont, Weight, weakL, weakR, SitesRemoved, SitesIgnored )
    CALL Bin_Data( DoS, Cont, Weight, BinValue, DroppedDos &
         , DoSMin, DoSMax )
    
    
  end subroutine Full_DoS

   subroutine AnalyzeClusters(WeakBonds, ClusterCount)
    implicit none
    integer loop1, loop2!loop integer
    integer, dimension(:), intent(in) :: WeakBonds !An array where each element is 0 or 1. 1 means weak bond, 0 means strong bond
    integer, intent(out) :: ClusterCount(dim) !Counts the number of clusters recorded for each cluster size
    integer ClusterSize !contains the size of the current cluster
    integer LbondLabel, RbondLabel !The number that labels the left and right weak bonds that neighbour each cluster

    do loop1 = 1,dim
       ClusterCount(loop1) = 0 !Initializes the number of clusters per cluster size to 0
    end do

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
       ClusterCount(dim) = ClusterCount(dim) + dim
    end if
       
    
  end subroutine AnalyzeClusters

  subroutine Bin_Data( HistoData, Data, Weights, BinValue, Dropped, Min, Max )
    implicit none
    real, intent(inout) :: HistoData(bins)
    real, dimension(:), allocatable, intent(in) :: Data
    real, intent(in) :: BinValue(bins)
    real, intent(in) :: Min, Max
    real, dimension(:), allocatable, intent(in) :: Weights
    integer, intent(inout) :: Dropped
    integer Loop1, Loop2
    logical BinFound

    do Loop1 = 1,size(Data)
       
       BinFound = .false.
       Loop2 = 1
       do while ( (Loop2 .le. bins) .and. (.not. BinFound) )
          if ( (Data(Loop1) .gt. Max) .or. (Data(Loop1) .lt. Min) ) then
             Dropped = Dropped + 1
             BinFound = .true.
          else if ( Data(Loop1) .le. BinValue(loop2) ) then
             HistoData(Loop2) = HistoData(Loop2) + Weights(Loop1)
             BinFound = .true.
          else
             Loop2 = Loop2 + 1
          end if
       end do
       
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
    integer NewStart
    integer i, j                                               ! Looping integer
    
    allocate(tmp_arr(size(array) - numRemove))

    if ( ( Start .gt. (Size(array) - numRemove) ) .and. &
         ( Start .lt. numRemove ) ) then

       NewStart = Start - numRemove
       
       if ( (NewStart + numRemove) .lt. size(array) + 1 ) then

          j=1
          tmp_arr(:) = 0
          !print*, Start
          do i = 1,(NewStart)
             tmp_arr(i) = array(i)
             j=j+1
          end do
       
          do i = (NewStart + numRemove + 1 ),size(array)
             
             tmp_arr(j) = array(i)
             
             j=j+1
          end do

       
          deallocate(array)
          allocate(array(size(tmp_arr)))
          array = tmp_arr
          
       else if ( NewStart .eq. Size(array) ) then
       
          j=1
          !print*, "edge"
          tmp_arr(:) = 0
          do i = mod(NewStart + NumRemove + 1 ,size(array)),NewStart
             
             tmp_arr(j) = array(i)
             j = j + 1
          end do
          
          deallocate(array)
          allocate(array(size(tmp_arr)))
          array = tmp_arr
       
       else
          !print*, "Start: ", Start
          j=1
          !print*, "Outside"
          tmp_arr(:) = 0
          !print*, "Array: ", array
          !print*, mod(Start + NumRemove + 1 ,size(array)), Start
          do i = mod(NewStart + NumRemove + 1 ,size(array)),NewStart 
          
             tmp_arr(j) = array(i)
             !print*, tmp_arr(j)
             
             j = j + 1
          end do
          
          deallocate(array)
          allocate(array(size(tmp_arr)))
          array = tmp_arr
          
       end if
       
    else
       
       NewStart = Start
       
       if ( (NewStart + numRemove) .lt. size(array) + 1 ) then

          j=1
          tmp_arr(:) = 0
          !print*, Start
          do i = 1,(NewStart)
             tmp_arr(i) = array(i)
             j=j+1
          end do
       
          do i = (NewStart + numRemove + 1 ),size(array)
             
             tmp_arr(j) = array(i)
             
             j=j+1
          end do

       
          deallocate(array)
          allocate(array(size(tmp_arr)))
          array = tmp_arr
          
       else if ( NewStart .eq. Size(array) ) then
       
          j=1
          !print*, "edge"
          tmp_arr(:) = 0
          do i = mod(NewStart + NumRemove + 1 ,size(array)),NewStart
             
             tmp_arr(j) = array(i)
             j = j + 1
          end do
          
          deallocate(array)
          allocate(array(size(tmp_arr)))
          array = tmp_arr
       
       else
          !print*, "Start: ", Start
          j=1
          !print*, "Outside"
          tmp_arr(:) = 0
          !print*, "Array: ", array
          !print*, mod(Start + NumRemove + 1 ,size(array)), Start
          do i = mod(NewStart + NumRemove + 1 ,size(array)),NewStart 
          
             tmp_arr(j) = array(i)
             !print*, tmp_arr(j)
             
             j = j + 1
          end do
          
          deallocate(array)
          allocate(array(size(tmp_arr)))
          array = tmp_arr
          
       end if
    end if
    

    
   
  end subroutine resize_array
  
  character(len=20) function str(k)                  ! function to change an integer to a string, used in FileNamer subroutine
     integer, intent(in) :: k

     write (str, *) k
     str = adjustl(str)
   end function str

   character(len=30) function FileNamer(Name)                                ! Increments filename labels Blah023.dat label:023
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
  use stuff
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
     real DoS(bins)                                      ! Density of States per system
     real DoSToT(bins)                                   ! Total Density of States
     real DoSMax                                         ! Maximum DoS value
     real DoSMin                                         ! Minimum DoS value
     real BinValue(bins)                                 ! Contains the right edge of all the bins in the density of states
     integer SitesRemoved                                ! Keeps track of the number of sites removed in each stage of the RG process
     integer SitesIgnored                                ! Keeps track of all the sites in clusters larger than this code can handle
     integer DroppedDos                                  ! number of DoS contributions outside of range [DoSMin,DoSMax]
     

               !----------------------------Coding Tools-------------------------------------
     real BinWidth                                       ! Width of the DoS bins for the purposes of normalization
     integer Loop1                                       ! Loop Integer
     character(len=25) filename                          ! Base file name to be put into FileNamer

     
               !----------------------------Initializing-------------------------------------
     DoSToT(1:bins) = 0
     DoSMax = 1.1*uSite                                  ! DoS extends to [-uSite,uSite], so I search for contributions at a range 
                                                         ! that is 10% larger
     DoSMin = -DoSMax
     DroppedDos = 0
     SitesRemoved = 0
     
     do loop1 = 1,bins
        BinValue(loop1) = DoSMin + (DoSMax-DoSMin)*loop1/real(bins)
     end do
     
  

               !----------------------------Opening files------------------------------------
     filename = "DoS"
     filename = FileNamer(filename)
     open (unit = 100, file = filename ,status ='unknown' )
     
     !Comments for the file above, says whats in the file and which parameters go with the data
     write(100,*) "#This file contains the density of states for the ensemble. First column is", & 
          " right edge of bin, fraction of contributions per bin"
     write(100,*) "#Disorder Strength = ", DELTA
     write(100,*) "#Bond Cutoff = ", bond_cutoff
     write(100,*) "#Interaction strength = ", uSite
     write(100,*) "#Dimensions = ", dim
     write(100,*) "#Number of systems = ", systemn
  
     
     call init_random_seed()
     DroppedDos = 0
     SitesIgnored = 0

                !---------------------------RG begins----------------------------------------
     
     do Loop1 = 1,systemn
        if ( mod(real(Loop1),0.1*systemn) == 0.0 ) then
           print*, int(Loop1/real(systemn)*100), "%"
        end if
        
        

                !---------------------------Create System------------------------------------
        call create_AHM( SitePotential, Hopping, Bonds )
        call CalcWeakBonds( Bonds, WeakBonds, SitesRemoved )
        if ( WeakBonds(1) .eq. 0 ) then
           print*, "All Bonds Strong"
        end if
                !---------------------------Calculate DoS------------------------------------
        call Full_DoS( SitePotential, Hopping, Bonds, WeakBonds, DoS, DoSMax, DoSMin, BinValue, DroppedDos, SitesIgnored )
        DosToT(1:bins) = DoSToT(1:bins) + DoS(1:bins)

     end do
     
     BinWidth = BinValue(2) - BinValue(1)                ! Area of each bin (bins have equal width)
     do Loop1=1,bins
        write(100,400) BinValue(Loop1) , DoSToT(Loop1)/( Sum(DoSToT) * BinWidth )
     end do

     print*, "Sites in large clusters: ", SitesIgnored, " out of", n_d
     print*, "DoS contributions outside of range: ", DroppedDos
  
!300  format(  i10, g12.5)                               ! Format for ( integer , real ) file
400  format( g12.5, g12.5)                               ! Format for ( real , real ) file
	 
   end program AndersonHubbard
		
   ! Last Edited: 180208
   ! Next Step: Properly comment things, explain orders of w1/w2/w3 when calculating the contribution, set weights to zero for some 
   ! each of the 4 possible contributions to show where they appear on the DoS. Write the outline for stage 2 write up
   !
