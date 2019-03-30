module parameters
  IMPLICIT NONE

  !-- Parameters -------------------------------------------------
  integer, parameter :: D=3           !D=2 or D=3  
  integer, parameter :: ScaleNum=128    !number of q bins of the external momentum
  double precision, parameter   :: UVScale=128.0     !the upper bound of energy scale
  integer            :: PID           ! the ID of this job
  integer            :: Order
  double precision   :: Mass2         ! screening length^2
  double precision   :: BareCoupling         ! bare coupling
  integer            :: Seed          ! random-number seed

  !-- Markov Chain ----------------------------------------------
  double precision                       :: Step        ! a counter to keep track of the current step number
  integer, parameter                     :: UpdateNum=4 ! number of updates
  double precision, dimension(UpdateNum) :: PropStep
  double precision, dimension(UpdateNum) :: AcceptStep

  !-- Diagram Tables  --------------------------------------------
  integer, parameter                    :: MaxOrder=4           ! Max diagram order
  integer                               :: LoopNum=2;

  double precision, dimension(ScaleNum)       :: ScaleTable(ScaleNum)
  integer                                     :: CurrScale
  double precision, dimension(MaxOrder)       :: CurrWeight
  double precision, dimension(D, MaxOrder+1)  :: LoopMom ! values to attach to each loop basis

  !-- common parameters and variables ------------------------------
  ! THIS IS PROJECT-INDEPENDENT 
  integer, parameter          :: UP=1, DOWN=0
  double precision, parameter :: tm32   = 1.d0/(2.d0**32.d0)
  double precision, parameter :: eps    = 1.d-14            ! very small number
  integer,          parameter :: Mxint  = 2147483647        ! maximum integer
  integer,          parameter :: Mnint  =-2147483647        ! minimum integer
  double precision, parameter :: pi=3.1415926
end module

INCLUDE "rng.f90"

program main
    use mt19937
    use parameters
    implicit none
    integer :: PrintCounter, SaveCounter, o
    double precision :: TotalStep  !total steps of this MC simulation
    double precision :: x
  
    print *, 'Mass2, BareCoupling, Order, TotalStep(*1e6), Seed, PID'
    read(*,*)  Mass2, BareCoupling, Order, TotalStep, Seed, PID

    ! For a given order, the bigger factor, the more accurate result 

    call sgrnd(Seed) 
  
    print *, "Initializing ..."
    call Initialize()
    print *, "Initialized!"
  
    call Test() !call test first to make sure all subroutines work as expected
  
    TotalStep=TotalStep*1.0e6
    print *, "Start simulation..."
    do while (Step<TotalStep)
      Step=Step+1.0
      x=grnd()
      if (x<1.0/UpdateNum) then
        call ChangeMom()
      endif
      !if(mod(int(Step), 4)==0) call Measure()
      call Measure()

      ! call DynamicTest()

      PrintCounter=PrintCounter+1
      if (PrintCounter==1e6)  then
        write(*,*) 
        write(*,"(f8.2, A15)") Step/1e6, "million steps"
        write(*,"(A20)") "Accept Ratio: "
        write(*,"(A16, f8.3)") "Change Mom:", AcceptStep(1)/PropStep(1)
        write(*,"(A20)") "Order Partition Sum: "
        print *, "mom1", norm2(LoopMom(:,1)), LoopMom(:,1)
        print *, "mom2", norm2(LoopMom(:,2)), LoopMom(:,2)
      endif
      if (PrintCounter==1e7)  then
        call SaveToDisk()
        PrintCounter=0
      endif

      !print *, CurrOrder

    end do

    call SaveToDisk()
    print *, "End simulation."
  
    CONTAINS

    subroutine Initialize()
      implicit none
      integer :: i, num
      double precision :: DeltaQ

      PropStep=0.0
      AcceptStep=0.0

      Step=0.0
      PrintCounter=0
      SaveCounter=0

    end subroutine

    subroutine Test()
      implicit none
      integer :: i
      double precision :: ratio, old, new
      ! Put all tests here
      !if (cabs(Green(0.d0, -1.d0)+Green(0.d0, Beta-1.d0))>1e-6) then
        !print *, "Green's function is not anti-periodic"
        !stop
      !endif

      ! old=0.0
      !call GenerateNewFreq(old, new, ratio)
      ! print *, old, new, ratio
      return
    end subroutine

    subroutine DynamicTest()
      implicit none
      return
    end subroutine

    double precision function Green(Mom, g_type)
      implicit none
      double precision :: k2
      integer :: g_type
      double precision, dimension(D) :: Mom
      k2=sum(Mom**2)
      if(k2>UVScale) then
        Green=0.0
      else
        Green=1.0/(k2+Mass2); 
      endif
      return
    end function
    
    double precision function Ver4(VerType, Scale)
      implicit none
      integer :: VerType, Scale
      Ver4=BareCoupling
      return
    end function Ver4

    !!!!!!!!!!!!!!!!! Balance the occurance of each order in the Markov Chain !!!!!!!!
    !subroutine ReWeightEachOrder()
      !implicit none
      !integer :: o
      !double precision :: TotalWeight
      !TotalWeight=sum(OrderPartitionSum(1:Order))
      !do o=1, Order
        !ReWeightFactor(o)=TotalWeight/OrderPartitionSum(o)
      !enddo
      !return
    !end subroutine

    subroutine Measure()
      implicit none
  
      !-------------------------------------------------------
      return
    end subroutine
    
    subroutine SaveToDisk()
      implicit none
      integer :: i, ref, j, o
      double precision :: Obs
      double precision :: DeltaQ
      character*10 :: ID
      character*10 :: order_str
      character*10 :: DiagIndex
      character*20 :: filename
      !Save Polarization to disk

      ! do o=1, Order
      !   write(ID, '(i10)') PID
      !   write(order_str, '(i10)') o
      !   filename="Data/Diag"//trim(adjustl(order_str))//"_"//trim(adjustl(ID))//".dat"
      !   write(*,*) "Save to disk ..."
      !   open(100, status="replace", file=trim(filename))
      !   write(100, *) "#", Step
      !   write(100, *) "#", Polarization(1, o)
      !   ref=int(kF/DeltaQ)+1
      !   do i=1, QBinNum
      !       Obs = Polarization(i, o)
      !       write(100, *) norm2(ExtMomMesh(:, i)), Obs
      !   enddo
      !   close(100)

      return
    end subroutine

    double precision function CalcWeight()
      !calculate the weight for ALL diagrams in a given sector
      implicit none
      CalcWeight=1.0;
      return
    end function CalcWeight
    
    subroutine ChangeMom()
      !randomly choose a vertex, change the space variable
      implicit none
      double precision, dimension(D) :: NewMom, OldMom
      double precision :: prop, R
      integer :: Num, i, j, NewExtMomBin
      double precision :: Weight
  
      Num = int( LoopNum*grnd() ) + 1
  
      PropStep(1) = PropStep(1) + 1.0
      OldMom = LoopMom(:, Num)
  
      call GenerateNewMom(OldMom, NewMom, prop)
  
      if(prop<0.0) return
      LoopMom(:,Num)=NewMom
  
      ! Weight = CalcWeight()
      ! R = prop*abs(Weight)/abs(CurrWeight)
  
      if(grnd()<R) then
        AcceptStep(1) = AcceptStep(1)+1.0
        ! CurrWeight = Weight
        ! call UpdateState()
      else
        LoopMom(:,Num)=OldMom
      endif
  
      return
    end subroutine
    
    subroutine GenerateNewMom(old, new, prop)
      implicit none
      double precision, dimension(D) :: old, new
      integer :: Direction, Num
      double precision :: ratio, prop, k, k_new
      double precision :: lambda
      x=grnd()
      if(x<1.0/3.0) then
          new=old
          Num=int(grnd()*D)+1
          new(Num)=new(Num)+sign(ScaleTable(CurrScale)*grnd(), grnd()-0.5)
          prop=1.0
      else if(x<3.0/3.0) then
          k=norm2(old) !sqrt(m_x^2+m_y^2+m_z^2)
          if(k==0.0)then
            prop=-1.0
            return
          endif
          lambda=1.5
          k_new=k/lambda+grnd()*(lambda-1.0/lambda)*k
          ratio=k_new/k
          new=old*ratio
          if(D==2) then
            prop=1.0 ! prop=k_old/k_new
          else
            prop=k_new/k ! prop=k_old/k_new
          endif
      else
          new=-old
          prop=1.0
      endif
    end subroutine
    
end program main
