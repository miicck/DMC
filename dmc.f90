! Constants used by the program
module constants
implicit none
    integer,    parameter :: prec = SELECTED_REAL_KIND(15,307)
    real(prec), parameter :: pi   = 3.141592625359
end module constants

! Utilities
module utils
use constants
implicit none

contains

    ! Return a normally distributed random number
    ! using a Box-Muller transform 
    function randNormal(var) result(res)
    implicit none
        real(prec) :: var, res, u1, u2
        call random_number(u1)
        call random_number(u2)
        res = sqrt(-2*log(u1))*sin(2*pi*u2)*sqrt(var)
    end function

end module utils

! Module containing the definition of a nucleus
module nucleus_type
use constants
implicit none

    ! Represents a nucleus
    type :: nucleus
        real(prec) :: position(3)
        real(prec) :: charge
    contains
        procedure :: copyTo => nucleusCopyTo
        procedure :: info   => nucleusInfo
    end type

contains

    ! Copy this nucleus to another nucleus
    subroutine nucleusCopyTo(this, other)
    implicit none
        class(nucleus) :: this, other
        other%position = this%position
        other%charge   = this%charge
    end subroutine

    ! Print nuclus info
    subroutine nucleusInfo(this)
    implicit none
        class(nucleus) :: this      
        write(*,"(F10.3, F10.3, F10.3, F10.3)") &
            this%charge, this%position(1), this%position(2), this%position(3)
    end subroutine

end module nucleus_type

! Module containing the definition of the walker type
module walker_type
use constants
use utils
use nucleus_type
implicit none

    ! Represents an electron configuration in the DMC scheme
    type :: walker
    private
        real(prec), allocatable :: electronPositions(:,:)
    contains 
        procedure :: energy          => walkerEnergy
        procedure :: potentialEnergy => walkerPotentialEnergy
        procedure :: kineticEnergy   => walkerKineticEnergy
        procedure :: diffuse         => walkerDiffuse
        procedure :: moveElectron    => walkerMoveElectron
        procedure :: copyTo          => walkerCopyTo
        procedure :: initialize      => walkerInitialize
    end type

contains

    ! Copy the walker this onto other
    subroutine walkerCopyTo(this, other)
    implicit none
        class(walker) :: this, other

        other%electronPositions = this%electronPositions

    end subroutine

    ! Move the i^th electron by x in the j^th direction
    subroutine walkerMoveElectron(this,i,j,x)
    implicit none
        class(walker) :: this
        integer       :: i, j
        real(prec)    :: x

        this%electronPositions(j,i) = &
            this%electronPositions(j,i) + x
        
    end subroutine

    ! Return the potential energy of this walker
    function walkerPotentialEnergy(this, nuclei) result(res)
    implicit none
        class(walker)  :: this
        class(nucleus) :: nuclei(:)
        real(prec)     :: res, disp(3)
        integer        :: i, j

        res = 0

        ! Calculate electron-nuclear interaction
        do i=1,size(this%electronPositions, 2)
            do j=1,size(nuclei)
                disp = this%electronPositions(:,i)
                disp = disp - nuclei(j)%position
                res = res - nuclei(j)%charge/norm2(disp)
            enddo
        enddo

        ! Calculate electron-electron interaction
        do i=1,size(this%electronPositions,2)
            do j=i+1,size(this%electronPositions,2)
                disp = this%electronPositions(:,i) - &
                       this%electronPositions(:,j)
                res = res + 1/norm2(disp)
            enddo
        enddo

    end function

    ! Returns the kinetic energy of this walker using
    ! the virial theorem and finite differences
    function walkerKineticEnergy(this, nuclei) result(res)
    implicit none
        real(prec), parameter :: EPS = 0.0001
        class(walker)         :: this
        class(nucleus)        :: nuclei(:)
        real(prec)            :: res, pot1, pot2
        real(prec)            :: centrePosition
        integer               :: i, j

        res = 0
        do i=1,size(this%electronPositions,2)
            do j=1,3
                centrePosition = this%electronPositions(j,i)

                call this%moveElectron(i,j,EPS/2) ! Move electron forward
                pot1 = walkerPotentialEnergy(this, nuclei)

                call this%moveElectron(i,j,-EPS)  ! Move electron back
                pot2 = walkerPotentialEnergy(this, nuclei)

                call this%moveElectron(i,j,EPS/2) ! Reset electorn
                res = res + 0.5 * centrePosition * (pot1 - pot2)/EPS ! Virial theorem
            enddo
        enddo

    end function

    ! Return the energy of this walker
    function walkerEnergy(this, nuclei) result(res)
    implicit none
        class(walker)  :: this
        class(nucleus) :: nuclei(:)        
        real(prec)     :: res
        res = walkerPotentialEnergy(this, nuclei) + walkerKineticEnergy(this, nuclei)
    end function

    ! Make a diffusion move on this walker
    subroutine walkerDiffuse(this, tau)
    implicit none
        class(walker) :: this
        real(prec)    :: tau
        integer       :: i, j
        
        do i=1,size(this%electronPositions,2)
            do j=1,3
                call this%moveElectron(i,j,randNormal(tau))
            enddo
        enddo

    end subroutine

    ! Initialize this walker
    subroutine walkerInitialize(this, electrons)
    implicit none
        class(walker) :: this
        integer :: i, j, electrons

        ! Allocate electron positions
        allocate(this%electronPositions(3, electrons)) 

        ! Initialize electrons to a guassian distribution
        do i=1,electrons
            do j=1,3
                this%electronPositions(j,i) = &
                    this%electronPositions(j,i) + &
                        randNormal(real(1,prec))
            enddo
        enddo

    end subroutine

end module walker_type

! Module to carry out dmc calculations
module dmc
use constants
use utils
use walker_type
implicit none

    ! Structure with info about a DMC step
    type :: step
        real(prec) :: averageEnergy
        real(prec) :: averageKinetic
        real(prec) :: averagePotential
        integer    :: population
    end type

    ! Structure with the best values calculated so far
    type :: bestValues
        real(prec) :: energy
        real(prec) :: kinetic
        real(prec) :: potential
    end type

    ! Program state
    class(walker), allocatable  :: walkers(:)               ! The walkers in the simulation
    class(nucleus), allocatable :: nuclei(:)                ! The nuclei in the simulation
    class(step), allocatable    :: steps(:)                 ! Information about the statistics accumulation steps
    type(bestValues)            :: bestVals                 ! The best accumulated answers so far
    real(prec)                  :: tau = 0.01               ! The DMC timestep
    integer                     :: electrons = 0            ! The number of electrons in the simulation
    integer                     :: targetPopulation = 1000  ! The target walker population
    real(prec)                  :: population_tol = 4       ! How far away from targetPopulation we are allowed to go (as a ratio)
    integer                     :: steps_stats = 1000       ! The number of DMC statistics gethering itterations to carry out
    integer                     :: steps_equil = 100        ! The number of DMC equilibriation steps to carry out
    integer                     :: equil_av_window = 10     ! The averaging window for quantities during equilibriation
    integer                     :: output_block_size = 100  ! Output every output_block_size itterations
    integer                     :: currentItteration = 0    ! The itteration we are currently on

contains

    ! Carry out a single DMC itteration, propagating the set of
    ! walkers and running the branching algorithm at their new
    ! positions. Eventually, the walkers will settle down into
    ! the ground state configuration.
    subroutine dmcItteration()
    implicit none
        integer                    :: i, j, k
        real(prec)                 :: potBeforeDefuse, potAfterDefuse
        integer                    :: surviving(size(walkers)), totalSurviving
        class(walker), allocatable :: walkersAfter(:)

        ! Initialize varaibles
        currentItteration = currentItteration + 1
        totalSurviving = 0

        do i=1,size(walkers)

            ! Calculate the potential energy of this walker
            ! before and after a diffusion move
            potBeforeDefuse = walkers(i)%potentialEnergy(nuclei)
            call walkers(i)%diffuse(tau)
            potAfterDefuse = walkers(i)%potentialEnergy(nuclei)

            ! Calculate the branching probability and how
            ! many walkers should survive here as a result
            surviving(i) = walkersSurvivingMove(potBeforeDefuse, potAfterDefuse)

            ! Keep track of how many survive
            totalSurviving = totalSurviving + surviving(i)

        enddo ! loop over walkers

        ! Create the new walker array and replace
        ! the current walker array with it
        allocate(walkersAfter(totalSurviving))
        k = 0
        do i=1,size(walkers)
            do j=1,surviving(i)
                k = k + 1
                call walkers(i)%copyTo(walkersAfter(k))
            enddo
        enddo
        deallocate(walkers)
        call move_alloc(walkersAfter, walkers)

        ! The itteration is complete
        call onCompleteItteration()
    
    end subroutine

    ! Called when an itteration completes
    subroutine onCompleteItteration()
    implicit none
        integer    :: i, stats_itter
        logical    :: equil
        type(step) :: info

        stats_itter = currentItteration - steps_equil
        if (stats_itter > 0) then
            call onCompleteStatsItter(stats_itter)
        else
            call onCompleteEquilItter()
        endif 

        if (mod(currentItteration,output_block_size) .eq. 0) then
            if (stats_itter > 0) then
                print *, "Stats itteration: ", stats_itter
            else
                print *, "Equilibriation itteration: ", currentItteration
            endif
            print *, "    Population: ", size(walkers)
            print *, "    Best estimate of energy: ", bestVals%energy
            print *, "                    kinetic: ", bestVals%kinetic
            print *, "                  potential: ", bestVals%potential
            print *, ""
        endif

    end subroutine

    ! An equilibriation ittration has completed
    subroutine onCompleteEquilItter()
    implicit none
        integer    :: i 
        real(prec) :: avPot, avKin, ratio
        avPot = 0
        avKin = 0
        do i=1,size(walkers)
            avPot = avPot + walkers(i)%potentialEnergy(nuclei)
            avKin = avKin + walkers(i)%kineticEnergy(nuclei)
        enddo
        avPot = avPot/real(size(walkers),prec)
        avKin = avKin/real(size(walkers),prec)

        ratio = 1/real(equil_av_window,prec)
        bestVals%energy = (avKin + avPot) * ratio + bestVals%energy * (1-ratio)
        bestVals%potential = avPot * ratio + bestVals%potential * (1-ratio)
        bestVals%kinetic = avPot * ratio + bestVals%kinetic * (1-ratio)
    end subroutine

    ! A statistics itteration has completed
    subroutine onCompleteStatsItter(itter)
    implicit none
        integer :: i, itter

        ! Fill the step structure for this itteration
        do i=1,size(walkers)
            steps(itter)%averagePotential = steps(itter)%averagePotential + walkers(i)%potentialEnergy(nuclei)
            steps(itter)%averageKinetic   = steps(itter)%averageKinetic   + walkers(i)%kineticEnergy(nuclei)
        enddo

        steps(itter)%averageKinetic   = steps(itter)%averageKinetic   / size(walkers)
        steps(itter)%averagePotential = steps(itter)%averagePotential / size(walkers)
        steps(itter)%averageEnergy    = steps(itter)%averageKinetic + steps(itter)%averagePotential
        steps(itter)%population       = size(walkers)

        ! Calculate the best energies as the average over all stats steps so far
        bestVals%potential = 0
        bestVals%kinetic   = 0
        do i=1,itter
            bestVals%kinetic   = bestVals%kinetic   + steps(i)%averageKinetic
            bestVals%potential = bestVals%potential + steps(i)%averagePotential
        enddo
        bestVals%kinetic   = bestVals%kinetic   / real(itter,prec)
        bestVals%potential = bestVals%potential / real(itter,prec)
        bestVals%energy = bestVals%kinetic + bestVals%potential

    end subroutine

    ! The branching algorithm:
    ! Calculates how many walkers should survive at x'
    ! after a walker moves from x -> x', v(x) -> v(x')
    function walkersSurvivingMove(potentialBefore, potentialAfter) result(res)
    implicit none
        real(prec) :: potentialBefore, potentialAfter, branchingProb
        real(prec) :: trialEnergy, averagePotential, randNumber, ratio
        integer    :: res

        ! No walkers survive for NaN potentials
        res = 0
        if (isnan(potentialBefore)) return
        if (isnan(potentialAfter))  return

        call random_number(randNumber)        

        ! Check if we're hitting max/min bounds
        ratio = size(walkers)/real(targetPopulation,prec)
        if (ratio > population_tol) then
            res = 0
            if (randNumber < 1/ratio) res = 1
            return
        else if (ratio < 1/population_tol) then
            ratio = 1/ratio
            res = floor(ratio)
            ratio = ratio - res
            if (randNumber < ratio) res = res + 1
            return
        endif

        ! Calculate the branching ratio and resulting number of surviving walkers
        trialEnergy = bestVals%energy - log(size(walkers)/real(targetPopulation,prec))
        averagePotential = 0.5 * (potentialBefore + potentialAfter)      
        branchingProb = exp(-tau * (averagePotential - trialEnergy))
        res = floor(randNumber + branchingProb)
    end function

    ! Initialize the calculation based on the
    ! input parameters
    subroutine initializeCalculation()
    implicit none
        integer :: i

        allocate(walkers(targetPopulation))
        do i=1,targetPopulation
            call walkers(i)%initialize(electrons)
        enddo

        allocate(steps(steps_stats))        

    end subroutine

    ! Run our DMC calculation
    subroutine runDMC()
    implicit none
        integer :: i
        
        call parseInput
        call printParameters
        call initializeCalculation

        write(*,*) "===== BEGIN DMC ====="

        do i=1,steps_equil + steps_stats
            call dmcItteration()
        enddo   
    
        call outputStats

    end subroutine

    ! Write out statistics info
    subroutine outputStats()
    implicit none
        integer :: i
        write(*,*) "=== WRITING OUTPUT ==="
        
        open(unit=1,file="stats")

        do i=1,steps_stats
            write(1,*) steps(i)%averageEnergy
        enddo

        close(unit=1)
        write(*,*) "    success"
        write(*,*)
    end subroutine

    ! Parse the input file
    subroutine parseInput()
    implicit none
    integer        :: ioStatus
    character(100) :: tag
    real(prec)     :: charge, x, y ,z
    integer        :: nu, nd, i

        allocate(nuclei(0))
        open(unit=1, file="input")
        do
            read(1, *, iostat=ioStatus) tag
            if (ioStatus > 0) cycle ! Error reading line
            if (ioStatus < 0) exit  ! End of file

            backspace(1)
            select case (tag)

                ! Read in target population
                case ("walkers")
                    read(1,*) tag, targetPopulation
                     
                ! Read in statistics step #
                case ("steps_stats")
                    read(1,*) tag, steps_stats

                ! Read in equilibiriation step #
                case ("steps_equil")
                    read(1,*) tag, steps_equil

                ! Read in the timestep
                case ("tau")
                    read(1,*) tag, tau

                ! Read in an atom
                case ("atom")
                    read(1,*) tag, charge, nu, nd, x, y ,z
                    call addAtom(charge, nu, nd, x, y ,z)

                ! Read in the output block size
                case ("output_block_size")
                    read(1,*) tag, output_block_size

                ! Read in the equilibriation averaging window
                case ("equil_av_window")
                    read(1,*) tag, equil_av_window

                ! Read in min walkers
                case ("population_tol")
                    read(1,*) tag, population_tol

                ! Unkown tag
                case default
                    print *, "Unkown tag: ", tag
                    call exit(-1)
                    read(1,*)

            end select
        enddo
        close(unit=1)
        
    end subroutine

    ! Print the parameters of the calculation
    subroutine printParameters()
    implicit none
        integer :: i

        ! Output general parameters
        write(*,*) "General parameters:"
        write(*,*) "    walkers:           ", targetPopulation
        write(*,*) "    electrons:         ", electrons
        write(*,*) "    steps_equil:       ", steps_equil
        write(*,*) "    steps_stats:       ", steps_stats
        write(*,*) "    tau:               ", tau        
        write(*,*) "    output_block_size: ", output_block_size

        ! Output nuclei
        write(*,*) ""        
        write(*,*) "Nuclei:"
        write(*,*) "    charge    x         y         z"
        do i=1,size(nuclei)
            call nuclei(i)%info
        enddo
        write(*,*) ""

    end subroutine

    ! Add an atom of the given charge and number of electrons
    ! at position x, y, z
    subroutine addAtom(charge, nu, nd, x, y, z)
    implicit none
        real(prec)                  :: charge, x, y, z
        class(nucleus), allocatable :: newNuclei(:)
        type(nucleus)               :: newNucleus
        integer                     :: nu, nd, i

        electrons = electrons + nu + nd;
        allocate(newNuclei(size(nuclei)+1))

        do i=1,size(nuclei)
            call nuclei(i)%copyTo(newNuclei(i))
        enddo

        newNucleus%position(1) = x
        newNucleus%position(2) = y
        newNucleus%position(3) = z
        newNucleus%charge = charge
       
        call newNucleus%copyTo(newNuclei(size(newNuclei)))

        deallocate(nuclei)
        call move_alloc(newNuclei, nuclei)
    end subroutine

end module dmc

program main
use dmc
implicit none
    call runDMC
end program main