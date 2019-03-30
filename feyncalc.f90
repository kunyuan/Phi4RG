module parameters
  IMPLICIT NONE

  !-- Parameters -------------------------------------------------
  integer, parameter :: D=3           !D=2 or D=3  
  integer, parameter :: QBinNum=64    !number of q bins of the external momentum
  integer            :: PID           ! the ID of this job
  integer            :: Order
  double precision   :: Beta, Mu, rs, EF, kF
  double precision   :: Mass2         ! screening length^2
  double precision   :: ExtMomMax     !the upper bound of the external Momentum 
  integer            :: ObsType       !0: measure zero ferq polarization, 1: measure equal time polarization
  integer            :: Seed          ! random-number seed

  !-- Markov Chain ----------------------------------------------
  double precision                       :: Step        ! a counter to keep track of the current step number
  integer                                :: CurrOrder   !keep track of the diagram order for the current state
  double precision                       :: CurrWeight  !keep track of the weight of the current state
  integer                                :: ExtMomBin   !keep track of the q bin of the external momentum in the Markov chain
  integer, parameter                     :: UpdateNum=4 ! number of updates
  double precision, dimension(UpdateNum) :: PropStep
  double precision, dimension(UpdateNum) :: AcceptStep

  !-- Diagram Tables  --------------------------------------------
  integer, parameter                    :: MaxOrder=8           ! Max diagram order
  integer, parameter                    :: MaxDiagNum=1024      ! Max diagram number 
  integer, parameter                    :: MaxIndepdentG=2048   ! Max indepdent Green function number 
  integer, parameter                    :: MaxIndepdentVer=1024 ! Max indepdent vertex number
  !integer, parameter                   :: MaxEK=10000          ! Used to tabulate G
  !integer, parameter                   :: MaxTau=10000         ! Used to tabulate G

  integer, dimension(MaxOrder)          :: iGNum, iVerNum       !Number of independent G and Vertex
  double precision, dimension(MaxOrder) :: ReWeightFactor       !reweightfactor for each order
  double precision, dimension(MaxOrder) :: OrderPartitionSum    !store partition sum for each order
  integer, dimension(MaxOrder)          :: GNum, VerNum, LoopNum, TauNum, HugenholtzNum !Number of G, Vertex, Loop, diagram, 

  integer, dimension(MaxDiagNum, 2*MaxOrder, MaxOrder)              :: GIndex  ! Index of Green function
  integer, dimension(MaxDiagNum, 2*MaxOrder, MaxOrder)              :: VerIndex ! Index of vertex
  double precision,dimension(MaxDiagNum, MaxOrder)                  :: SymFactor  ! Symmetry Factor (includes diagram sign)
  double precision,dimension(2**(MaxOrder-1), MaxDiagNum, MaxOrder) :: SpinFactor  ! Spin Factor (includes diagram sign)
  double precision,dimension(2**(MaxOrder-1), MaxOrder-1)           :: SpinCache  ! Spin Factor (includes diagram sign)
  integer, dimension(MaxOrder+1, MaxIndepdentG, MaxOrder)           :: LoopBases ! Bases for loops
  integer, dimension(MaxOrder+1, MaxIndepdentVer, MaxOrder)         :: LoopBasesVer ! Bases for loops including vertex
  integer, dimension(MaxIndepdentG, MaxOrder)                       :: GType 
  integer, dimension(MaxIndepdentVer, MaxOrder)                     :: VerType 
  integer, dimension(2*MaxOrder, MaxIndepdentG, MaxOrder)           :: TauBases ! Permutation 
  double precision, dimension(MaxDiagNum)                           :: DiagWeight, DiagWeightABS

  double precision, dimension(MaxIndepdentG)                        :: GWeight, NewGWeight !Weight of green function
  double precision, dimension(MaxIndepdentVer)                      :: VerWeight, NewVerWeight ! Weight of vertex (interation)

  double precision, dimension(D, MaxOrder+1)                        :: LoopMom ! values to attach to each loop basis
  integer, dimension(MaxOrder+1)                                    :: LoopSpin ! values to attach to each spin
  double precision, dimension(2*MaxOrder)                           :: TauTable ! Time table for each Tau, all Tau are between [0,beta)

  !-- Measurement  ------------------------------------------------
  double precision, dimension(20)                            :: F1, F2, F3   !three weight functions
  double precision, dimension(QBinNum, MaxOrder)             :: Polarization !the accumulated total weight of the Spin-zz polarization
  double precision, dimension(MaxDiagNum, QBinNum, MaxOrder) :: PolarDiag    !the accumulated weight for each diagram of the Spin-zz polarization
  double precision, dimension(D, QBinNum)                    :: ExtMomMesh

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
  
    print *, 'Beta, rs, Mass2, Order, MaxExtMom(*kF), TotalStep(*1e6), Observable, Seed, PID'
    read(*,*)  Beta, rs, Mass2, Order, ExtMomMax, TotalStep, ObsType, Seed, PID

    ! For a given order, the bigger factor, the more accurate result 
    ReWeightFactor(1:3)=(/1.0,1.0,20.0/)

    if(D==3) then
      kF=(9.0*pi/4.0)**(1.0/3.0)/rs !3D
    else if (D==2) then
      kF=sqrt(2.0)/rs !2D
    else
      print *, "Dimension ", D, " has not yet been implemented!"
      stop
    endif
    EF=kF*kF
    Mu=EF
    ExtMomMax=ExtMomMax*kF
    Beta=Beta/EF
    print *,"Inverse Temperature:", Beta
    print *,"rs:", rs
    print *,"Fermi Mom:", kF
    print *,"Fermi Energy:", EF
  
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
        call ChangeTau()
      else if (x<2.0/UpdateNum) then
        call ChangeMom()
      else if (x<3.0/UpdateNum) then
        call IncreaseOrder()
      else if (x<4.0/UpdateNum) then
        call DecreaseOrder()
      ! else if (x<5.0/UpdateNum) then
        !!call ChangeSpin()
        !call SwapMom()
      endif
      !if(mod(int(Step), 4)==0) call Measure()
      call Measure()

      ! call DynamicTest()

      PrintCounter=PrintCounter+1
      OrderPartitionSum(CurrOrder)=OrderPartitionSum(CurrOrder)+1.0/ReWeightFactor(CurrOrder)
      if (PrintCounter==1e6)  then
        write(*,*) 
        write(*,"(f8.2, A15)") Step/1e6, "million steps"
        write(*,"(A20)") "Accept Ratio: "
        write(*,"(A16, f8.3)") "Increase Order:", AcceptStep(1)/PropStep(1)
        write(*,"(A16, f8.3)") "Decrease Order:", AcceptStep(2)/PropStep(2)
        write(*,"(A16, f8.3)") "Change Tau:", AcceptStep(3)/PropStep(3)
        write(*,"(A16, f8.3)") "Change Mom:", AcceptStep(4)/PropStep(4)
        write(*,"(A20)") "Order Partition Sum: "
        write(*,*) OrderPartitionSum(1:Order)/sum(OrderPartitionSum(1:Order))
        ! write(*,*) ReWeightFactor(1:Order)
!        write(*,"(A14, A6, f8.3)") "Accept Ratio: ", "Swap:", AcceptStep(3)/PropStep(3)
        print *, "order:", CurrOrder
        print *, "mom1", norm2(LoopMom(:,1)), LoopMom(:,1)
        print *, "mom2", norm2(LoopMom(:,2)), LoopMom(:,2)
        ! print *, "Momnorm", norm2(LoopMom(:,2)+LoopMom(:,1))**2-Mu, norm2(LoopMom(:,2))**2-Mu
        !print *, "mom3", norm2(LoopMom(:,3))
        !print *, "mom4", norm2(LoopMom(:,4))
        ! print *, "tau table", TauTable(1:TauNum(CurrOrder))
        ! print *, "green", GWeight(1), GWeight(2)
        ! print *, "weight", CurrWeight
        ! call ReWeightEachOrder()
      endif
      if (PrintCounter==1e7)  then
        call SaveToDisk()
        PrintCounter=0
      endif

      !print *, CurrOrder

    end do

    call SaveToDisk()
    call SaveToDiskAdditional()
    
    print *, "End simulation."
  
    CONTAINS

    subroutine Initialize()
      implicit none
      integer :: i, num
      double precision :: DeltaQ

      CurrOrder=Order

      print *, "Reading Diagram ..."
      call ReadDiagram()
      print *, "Read Diagram done!"
      !call CreateGTable()

      PropStep=0.0
      AcceptStep=0.0

      !ExtMomMax = 3.0*kF
      DeltaQ=ExtMomMax/QBinNum

      Polarization=0.0
      PolarDiag=0.0
      F1 = 0.0
      F2 = 0.0
      F3 = 0.0

      ExtMomMesh=0.0
      do i=1, QBinNum
        ExtMomMesh(1, i)=(i-1)*DeltaQ
      enddo

      Step=0.0
      PrintCounter=0
      SaveCounter=0

      ExtMomBin=1

      do i=0, CurrOrder-1
        num=2*i+1
        if(num==1) then
          TauTable(num)=0.0
          TauTable(num+1)=1.e-6
        else
          TauTable(num)=Beta*grnd()
          TauTable(num+1) = TauTable(num)
        endif
      enddo

      LoopMom(:,1)=ExtMomMesh(:, ExtMomBin)
      LoopMom(:,2:)=kF/sqrt(real(D))
      LoopSpin(1)=0
      LoopSpin(2)=1

      ! OrderPartitionSum=1.0
      ! call ReWeightEachOrder()

      call ResetWeightTable(CurrOrder)
      CurrWeight = CalcWeight(0, CurrOrder)
      call UpdateState()
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
      integer   :: i,j
      double precision   :: weight
      double precision :: tau,k2, Ek1, Ek2
  
      call ResetWeightTable(CurrOrder)
      weight = CalcWeight(0, CurrOrder)
      if(abs(weight-CurrWeight)>1e-6) then
        print *, "Weight is wrong at step:", Step, weight, CurrWeight
        print *, "Order: ", CurrOrder
        print *, "TauTable", TauTable(1:TauNum(CurrOrder))
        print *, "LoopMom", LoopMom(:, 1:LoopNum(CurrOrder))
        print *, "GWeight", GWeight(1:GNum(CurrOrder))
        print *, "VerWeight", VerWeight(1:VerNum(CurrOrder))

        print *, "NewGWeight", NewGWeight(1:GNum(CurrOrder))
        print *, "NewVerWeight", NewVerWeight(1:VerNum(CurrOrder))
        stop
      endif
  
      do i=1, LoopNum(CurrOrder)
        do j=1, D
          if(isnan(LoopMom(j, i))) then
            print *, "Mom is NaN",Step, j, i, LoopMom(j, i)
            stop
          endif
        enddo
      enddo
  
      do i=1, iVerNum(CurrOrder)
        if(isnan(VerWeight(i))) then
          print *, "VerWeight is NaN",Step, i, VerWeight(i)
          stop
        endif
  
      enddo
  
      if(isnan(CurrWeight)) then
        print *, "CurrWeight is NaN",Step, i, CurrWeight
        stop
      endif
  
      if(CurrWeight==0.0) then
        print *, "CurrWeight is zero",Step, i, CurrWeight
        stop
      endif
  
      if(abs(ExtMomMesh(1, ExtMomBin)-LoopMom(1,1))>1e-6) then
        print *, "ExtMom is wrong!",Step, ExtMomBin, ExtMomMesh(:, ExtMomBin), LoopMom(:,1)
      endif
    end subroutine

    subroutine ReadDiagram()
      implicit none
      character(90) :: fname
      character (len=20) :: charc
      integer :: Offset, OffVer, OffDiag
      integer :: baseNum, i, numDiagV, o
      integer, dimension(MaxOrder) :: DiagNum1H, DiagNum 
      character( len = 3 ) :: Orderstr

      do o=1, Order

        OffDiag = 0
        Offset = 0
        OffVer = 0

        write( Orderstr,'(I3)' )  o
        fname = 'DiagPolar'//trim(adjustl(Orderstr))//'.txt'
    
        open(unit=10, file=trim(fname), action='read', status="old")
          Read(10, *) charc  !the DiagNum comment
          Read(10, *) HugenholtzNum(o)

          GNum(o) = 2*o  	
          VerNum(o) = o-1 	
          LoopNum(o) = o+1	
          TauNum(o) = 2*o
          DiagNum1H(o) = 2**VerNum(o)

          DiagNum(o) = DiagNum1H(o) * HugenholtzNum(o)  
          iGNum(o) = GNum(o) * HugenholtzNum(o)
          iVerNum(o) = 2 * VerNum(o)  * HugenholtzNum(o) 

          do numDiagV=1, HugenholtzNum(o)
            GIndex(numDiagV, 1:GNum(o), o) = (/ ((numDiagV-1)*GNum(o)+i, i=1, GNum(o)) /)
            VerIndex(numDiagV, 1:2*VerNum(o), o) = (/ ( (numDiagV-1)*2*VerNum(o) + i, i=1, 2*VerNum(o) ) /)

            Read(10, *) charc   !Permuation comment
            Read(10, *) TauBases(2, Offset+1:Offset+GNum(o), o)
            TauBases(1, Offset+1:Offset+GNum(o), o) = (/ (i, i=1,GNum(o)) /) 
            TauBases(2, Offset+1:Offset+GNum(o), o) = TauBases(2, Offset+1:Offset+GNum(o), o) + 1
            Read(10, *) charc
            Read(10, *) GType(Offset+1:Offset+GNum(o), o)
            Read(10, *) charc
            Read(10, *) SymFactor(numDiagV, o)
            Read(10, *) charc
            do baseNum=1, LoopNum(o)
              Read(10, *)  LoopBases(baseNum, Offset+1:Offset+GNum(o), o)
            end do
            Read(10, *) charc

            if(o>1) then
              do baseNum=1, LoopNum(o)
                Read(10, *)  LoopBasesVer(baseNum, OffVer+1:OffVer+2*VerNum(o), o)
              end do
            endif
            Read(10, *) charc
            ! print *, charc
            if(o>1) then
              Read(10, *)  VerType(OffVer+1:OffVer+2*VerNum(o), o)
              ! print *, o, VerType(OffVer+1:OffVer+2*VerNum(o),o)
            endif
            Read(10, *) charc
            Read(10, *) SpinFactor(1:2**VerNum(o), numDiagV, o)
            ! print *, SpinFactor(1:2**VerNum(o), numDiagV, o)

            OffDiag = OffDiag + DiagNum1H(o)
            Offset = Offset + GNum(o)
            OffVer = OffVer + 2*VerNum(o)
          end do
    
        CLOSE(unit=10)
      enddo
    end subroutine

    !!!!!!!!!!!!!!!!   Fock diagram self energy   !!!!!!!!!!!!!!!!!!!!!!!!!
    double precision function SelfEnergy(Mom)
      implicit none
      double precision, dimension(D) :: Mom
      double precision :: k, l, shift
      l=sqrt(Mass2)
      !l=100.0
      k=norm2(Mom)
      SelfEnergy=1.0+l/kF*atan((k-kF)/l)
      SelfEnergy=SelfEnergy-l/kF*atan((k+kF)/l)
      SelfEnergy=SelfEnergy-(l*l-k*k+kF*kF)/4.0/k/kF*log((l*l+(k-kF)**2)/(l*l+(k+kF)**2))
      SelfEnergy=SelfEnergy*(-2.0*kF)/pi

      shift=1.0-l/kF*atan(2.0*kF/l)
      shift=shift-l*l/4.0/kF/kF*log(l*l/(l*l+4*kF**2))
      shift=shift*(-2.0*kF)/pi
      SelfEnergy=SelfEnergy-shift
      SelfEnergy=SelfEnergy+k*k
      return
    end function

    double precision function Green(tau ,Mom, spin, g_type)
      implicit none
      double precision :: tau
      integer :: spin, g_type
      double precision, dimension(D) :: Mom
      if(g_type==0) then
        Green=PhyGreen(tau, Mom)
      else
        ! Green=PhyGreen(tau, Mom)
        ! print *, tau
        ! Green=PhyGreen(tau-beta/3.0, Mom)*PhyGreen(+beta/3.0, Mom)*beta
        Green=FakeGreen(tau, Mom)
      endif
    end function

    !!!!!!!!!!!!!!!!! Green's function for free fermion   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double precision function PhyGreen(tau ,Mom)
      implicit none
      double precision :: tau, k2, s, Ek, x, y, w, r, coshv
      double precision, dimension(D) :: Mom

      ! if tau is exactly zero, set tau=0^-
      if(tau==0.0) tau=-eps

      s=1.0
      if(tau<0.0) then
        tau=beta+tau
        s=-s
      endif
      if(tau>=beta) then
        tau=tau-beta
        s=-s
      endif
       !Ek=sum(Mom**2)   ! bare propagator

      Ek=SelfEnergy(Mom)   ! Fock diagram dressed propagator

      x=Beta*(Ek-Mu)/2.0
      y=2.0*tau/Beta-1.0
      if(x>100.0) then
        PhyGreen=dexp(-x*(y+1.0))
      else if(x<-100.0) then
        PhyGreen=dexp(x*(1.0-y))
      else
        PhyGreen=dexp(-x*y)/(2.0*cosh(x))
      endif
      !if(spin==1 .or. spin==-1) then
      PhyGreen=s*PhyGreen
  
      if(isnan(PhyGreen)) then
        print *, "Green is too large!", tau, Ek, PhyGreen
        stop
      endif
      return
    end function

    !!!!!!!!!!!!!!!!! Green's function for free fermion   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double precision function FakeGreen(tau ,Mom)
      implicit none
      double precision :: tau, k2, s, Ek, x, y, w, r, coshv, Factor
      double precision, dimension(D) :: Mom

      ! if tau is exactly zero, set tau=0^-
      if(tau==0.0) tau=-eps

      s=1.0
      if(tau<0.0) then
        tau=beta+tau
        s=-s
      endif
      if(tau>=beta) then
        tau=tau-beta
        s=-s
      endif
       !Ek=sum(Mom**2)   ! bare propagator

      Ek=SelfEnergy(Mom)   ! Fock diagram dressed propagator

      x=Beta*(Ek-Mu)/2.0
      y=2.0*tau/Beta-1.0
      if(x>100.0) then
        FakeGreen=dexp(-x*(y+1.0))
        Factor=tau
      else if(x<-100.0) then
        FakeGreen=dexp(x*(1.0-y))
        Factor=-(beta-tau)
      else
        FakeGreen=dexp(-x*y)/(2.0*cosh(x))
        Factor=tau*dexp(x)/(2.0*cosh(x))-(beta-tau)*dexp(-x)/(2.0*cosh(x))
      endif
      !if(spin==1 .or. spin==-1) then
      FakeGreen=s*FakeGreen*Factor

      if(isnan(FakeGreen)) then
        print *, "Green is too large!", tau, Ek, FakeGreen
        stop
      endif
      return
    end function


    !!!!!!!!!!!!!!!!! Green's function for phi4 model   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !double precision function Green(tau ,Mom, spin)
      !!calculate Green's function
      !implicit none
      !double precision :: tau, k2, s, Ek, x, y, w, r, coshv
      !integer :: spin, i
      !double precision, dimension(D) :: Mom
  !!    print *, tau, Mom, Spin
  !!    stop
      !Ek=sum(Mom**2)   !kinetic energy
      !Green=1.0/(Ek+Mass2)
  
      !return
   !end function Green

  !!!!!!!!!!!!!!!!!!!! Tabulated Green's function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !subroutine CreateGTable()
      !implicit none
      !integer :: iEk, iTau, i
      !double precision :: dEk, dTau, Ek, Tau, x
      !dEk=5.0*EF/MaxEK
      !dTau=Beta/MaxTau
      !do iEk=0, MaxEK
        !do iTau=0, MaxTau-1 !treat the last point, namely tau==beta later
          !Ek=iEk*dEk
          !Tau=iTau*dTau
          !GreenTable(iEk, iTau)=GreenFunc(Tau, Ek, 1)
        !enddo
        !GreenTable(iEk, MaxTau)=GreenFunc(Beta-1.0e-10, Ek, 1) !set the largest table point to tau--->beta^-
      !enddo
    !end subroutine

    !double precision function Green(tau, Mom, spin)
      !implicit none
      !integer :: iEK, iTau, spin
      !double precision :: Ek, s, tau, wEk, wTau, rEk, rTau
      !double precision :: g11, g12, g21, g22
      !double precision :: ratioEk, ratioTau
      !double precision, dimension(D) :: Mom
      !Ek=sum(Mom**2)   !kinetic energy
      !s=1.0
      !if(tau<0) then
        !tau=beta+tau
        !s=-s
      !endif
      !if(tau>=beta) then
        !tau=tau-beta
        !s=-s
      !endif
      !if(Ek>=5.0*EF) then
        !Green=GreenFunc(tau, Ek, spin)
      !else
        !wEk=Ek/5.0/EF*MaxEK
        !wTau=tau/Beta*MaxTau
        !iEK=int(wEk)
        !iTau=int(wTau)
        !rEk=wEk-iEk
        !rTau=wTau-iTau
        !Green=(1-rTau)*((1-rEk)*GreenTable(iEk, iTau)+rEk*GreenTable(iEk+1, iTau))
        !Green=Green+rTau*((1-rEk)*GreenTable(iEk, iTau+1)+rEk*GreenTable(iEk+1, iTau+1))
        !Green=s*Green
        !if(abs(Green-s*GreenFunc(tau, Ek, spin))>1e-4) then
          !print *, Ek, tau
          !print *, iEk, iTau
          !print *, wEk, wtau
          !print *, rEk, rTau
          !print *, GreenTable(iEk, iTau), GreenTable(iEk+1, iTau), GreenTable(iEk, iTau+1), GreenTable(iEk+1, iTau+1)
          !print *, Green, s*GreenFunc(tau, Ek, spin) 
          !stop
        !endif
      !endif
    !end function Green
    
    double precision function Interaction(tau, Mom, spin, VerType)
      implicit none
      double precision :: tau
      double precision, dimension(D) :: Mom
      integer :: spin, VerType
      if(VerType<0) then
        print *, "VerType can not be ", VerType
        stop
      endif

      Interaction=8.0*pi/(sum(Mom**2)+Mass2)

      if(VerType>0)then
        !the interaction contains counter-terms
        Interaction=Interaction*(Mass2/(sum(Mom**2)+Mass2))**VerType
        Interaction=Interaction*(-1)**VerType
      endif

      return
    end function Interaction

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
      integer :: Num, i, HugenNum
      double precision :: AbsWeight, Q, ReWeight, Weight
      double precision :: Phase
  
      AbsWeight=abs(CurrWeight)
      Phase=CurrWeight/AbsWeight
  
      Num=ExtMomBin
      !Q=norm2(ExtMomMesh(:, ExtMomBin))
      !ReWeight=exp(Q*Q/1.d0/kF/kF)
      Reweight=1.0/ReWeightFactor(CurrOrder)
  
      call NewState()
      Weight=CalcWeight(1, CurrOrder)
  
      !------------- measure sign blessing ------------------
      HugenNum=HugenholtzNum(CurrOrder)
  
      do i=1, HugenNum
        F1(i) = F1(i) + DiagWeight(i)/AbsWeight*ReWeight
        F2(i) = F2(i) + ABS(DiagWeight(i))/AbsWeight*ReWeight
        F3(i) = F3(i) + DiagWeightABS(i)/AbsWeight*ReWeight
      enddo

      F1(HugenNum+1) = F1(HugenNum+1) + SUM(DiagWeight(1:HugenNum))/AbsWeight*ReWeight
      F2(HugenNum+1) = F2(HugenNum+1) + ABS(SUM(DiagWeight(1:HugenNum)))/AbsWeight*ReWeight
      F3(HugenNum+1) = F3(HugenNum+1) + SUM(ABS(DiagWeightABS(1:HugenNum)))/AbsWeight*ReWeight
      !-------------------------------------------------------

      do i = 1, HugenNum
        PolarDiag(i, Num, CurrOrder) = PolarDiag(i, Num, CurrOrder) + DiagWeight(i)/AbsWeight*ReWeight
        Polarization(Num, CurrOrder) = Polarization(Num, CurrOrder) + DiagWeight(i)/AbsWeight*ReWeight
      enddo
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

      DeltaQ=ExtMomMax/QBinNum

      do o=1, Order
        write(ID, '(i10)') PID
        write(order_str, '(i10)') o
        filename="Data/Diag"//trim(adjustl(order_str))//"_"//trim(adjustl(ID))//".dat"
        write(*,*) "Save to disk ..."
        open(100, status="replace", file=trim(filename))
        write(100, *) "#", Step
        write(100, *) "#", Polarization(1, o)
        ref=int(kF/DeltaQ)+1
        do i=1, QBinNum
            Obs = Polarization(i, o)
            write(100, *) norm2(ExtMomMesh(:, i)), Obs
        enddo
        close(100)

        do j=1, HugenholtzNum(o)
            write(DiagIndex, '(i10)') j
            filename="Data/Diag"//trim(adjustl(order_str))//"_"//trim(adjustl(ID))//"_"//trim(adjustl(DiagIndex))//".dat"
            open(100, status="replace", file=trim(filename))
            write(100, *) "#", Step
            write(100, *) "#", PolarDiag(j,1, o)
            ref=int(kF/DeltaQ)+1
            do i=1, QBinNum
              Obs=PolarDiag(j, i, o)
              write(100, *) norm2(ExtMomMesh(:, i)), Obs
            enddo
            close(100)
        enddo
      enddo

      return
    end subroutine

    subroutine SaveToDiskAdditional()
      implicit none
      integer :: i, ref, j, o
      double precision :: Obs
      double precision :: DeltaQ
      character*10 :: ID
      character*20 :: filename
      !Save Polarization to disk

      DeltaQ=ExtMomMax/QBinNum

      filename="output.dat"
      open(100, position="append", file=trim(filename))
      write(100, *) "#para", rs, Beta, PID
      do o=1, Order
        write(100, *) o, Polarization(1, o)
      enddo
      write(100, *)
      close(100)
    
    end subroutine
    
    !subroutine SaveToDiskF()
      !implicit none
      !character*20 :: filename
      !integer :: i
      !double precision :: res=0.0
      !do i=1, HugenholtzNum
          !res = res + F2(i)
      !enddo

      !filename="F123.dat"
      !open(100, position='append', file=trim(filename))
          !write(100, *)  "OLD BASES:  q=", LoopMom(1,1) 
          !do i=1, HugenholtzNum+1
              !write(100, *) "Hugenholtz", i,":", F1(i)/F3(HugenholtzNum+1), F2(i)/F3(HugenholtzNum+1), F3(i)/F3(HugenholtzNum+1)
          !enddo
          !write(100, *) "F2   SUM:", res/F3(HugenholtzNum+1)
      !close(100)
        
    !end subroutine

    !check if the Ver is one-interaction-line reducible or not
    integer function IsReducible(index, o)
        implicit none
        integer :: index, o
        !!!! Reducibility check   !!!!!!!!!!!!!!!!!!
        if(abs(LoopBasesVer(1, index, o))==1 .and. sum(LoopBasesVer(2:LoopNum(o), index, o)**2)==0) then
          IsReducible=1
        else
          IsReducible=0
        endif
        return
    end function
   
    !reset G and Ver weight table
    subroutine ResetWeightTable(NewOrder)
      implicit none
      integer :: i, j, Spin, NewOrder
      double precision :: Tau
      double precision, dimension(D) :: Mom
      !print *, "iGNum", iGNum
      call NewState()
      do i=1, iGNum(NewOrder)
        Tau = TauTable(TauBases(2, i, NewOrder))-TauTable(TauBases(1, i, NewOrder))
        Spin=sum(LoopBases(1:LoopNum(NewOrder), i, NewOrder)*LoopSpin(1:LoopNum(NewOrder)))
        do j=1, D
          Mom(j)=sum(LoopBases(1:LoopNum(NewOrder), i, NewOrder)*LoopMom(j, 1:LoopNum(NewOrder)))
        enddo
        NewGWeight(i)=Green(Tau, Mom, Spin, GType(i, NewOrder))
      enddo
  
  
      do i=1, iVerNum(NewOrder) 
          Spin=sum(LoopBasesVer(1:LoopNum(NewOrder), i, NewOrder)*LoopSpin(1:LoopNum(NewOrder)))
          do j=1, D
              Mom(j)=sum(LoopBasesVer(1:LoopNum(NewOrder), i, NewOrder)*LoopMom(j, 1:LoopNum(NewOrder)))
          enddo

          NewVerWeight(i) = Interaction(0.d0, Mom, Spin, VerType(i, NewOrder))
          !!!! Reducibility check   !!!!!!!!!!!!!!!!!!
          if(IsReducible(i, NewOrder)==1) then
            NewVerWeight(i)=0.0
          endif
      enddo
  
      return
    end subroutine

    ! Copy the weights of G and Ver for the current state to the new state
    subroutine NewState()
      implicit none
      NewGWeight(1:iGNum(Order))=GWeight(1:iGNum(Order))
      NewVerWeight(1:iVerNum(Order))=VerWeight(1:iVerNum(Order))
    end subroutine
    
    ! Update the weights of G and Ver for the current state from the new state
    subroutine UpdateState()
      implicit none
      GWeight(1:iGNum(Order))=NewGWeight(1:iGNum(Order))
      VerWeight(1:iVerNum(Order))=NewVerWeight(1:iVerNum(Order))
    end subroutine
    
    subroutine ChangeWeightTableMom(loopindex)
      implicit none
      integer :: i, j, Spin, loopindex
      double precision :: Tau
      double precision, dimension(D) :: Mom
      call NewState()
      do i=1, iGNum(CurrOrder) 
          if (LoopBases(loopindex, i, CurrOrder)==0) then
              cycle !if the loop with #loopindex do not affect the Gline, then do not recalculate 
          endif
  
          Tau = TauTable(TauBases(2, i, CurrOrder))-TauTable(TauBases(1, i, CurrOrder))
          !Spin = sum(LoopBases(:, i)*LoopSpin(:))
          Spin=1
          do j=1, D
              Mom(j) = sum(LoopBases(1:LoopNum(CurrOrder), i, CurrOrder)*LoopMom(j, 1:LoopNum(CurrOrder)))
          enddo
          NewGWeight(i) = Green(Tau, Mom, Spin, GType(i, CurrOrder))
      enddo
  
  
      do i=1, iVerNum(CurrOrder) 

          if (LoopBasesVer(loopindex, i, CurrOrder)==0) then
            cycle !if the loop with #loopindex do not affect the Gline, then do not recalculate 
          endif
    
          Spin=sum(LoopBasesVer(1:LoopNum(CurrOrder), i, CurrOrder)*LoopSpin(1:LoopNum(CurrOrder)))
          do j=1, D
            Mom(j)=sum(LoopBasesVer(1:LoopNum(CurrOrder), i, CurrOrder)*LoopMom(j, 1:LoopNum(CurrOrder)))
          enddo
          NewVerWeight(i) = Interaction(0.d0, Mom, Spin, VerType(i, CurrOrder))
          !!!! Reducibility check   !!!!!!!!!!!!!!!!!!
          if(IsReducible(i, CurrOrder)==1) then
            NewVerWeight(i)=0.0
          endif
      enddo
  
      return
    end subroutine
    
    subroutine ChangeWeightTableTau(Tauindex)
      implicit none
      integer :: i, j, Spin, Tauindex
      double precision :: Tau
      double precision, dimension(D) :: Mom
      call NewState()
      do i=1, iGNum(CurrOrder)
        if (TauBases(1, i, CurrOrder)==Tauindex .or. TauBases(2, i, CurrOrder)==Tauindex .or. &
          TauBases(1, i, CurrOrder)==Tauindex+1 .or. TauBases(2, i, CurrOrder)==Tauindex+1) then

          Tau=TauTable(TauBases(2, i, CurrOrder))-TauTable(TauBases(1, i, CurrOrder))
          Spin=sum(LoopBases(1:LoopNum(CurrOrder), i, CurrOrder)*LoopSpin(1:LoopNum(CurrOrder)))
          do j=1, D
            Mom(j)=sum(LoopBases(1:LoopNum(CurrOrder), i, CurrOrder)*LoopMom(j, 1:LoopNum(CurrOrder)))
          enddo
          NewGWeight(i)=Green(Tau, Mom, Spin, GType(i, CurrOrder))
        endif
      enddo
      return
    end subroutine

    double precision function CalcWeight(SaveWeight, NewOrder)
      !calculate the weight for ALL diagrams in a given sector
      implicit none
      integer :: i, j, k, BlockNum,  NewOrder
      ! double precision :: Q2, Q
      double precision :: TotalGWeight, TotalVerWeight, AbsTotalVerWeight
      integer :: SaveWeight
      CalcWeight=0.0
  
      do i=1, HugenholtzNum(NewOrder)
        TotalGWeight = SymFactor(i, NewOrder)  !initialize weight wiht the symmetry factor
        do j=1, GNum(NewOrder)
          TotalGWeight = TotalGWeight * NewGWeight(GIndex(i, j, NewOrder))
        enddo

        !!!!!!!!!!!!!! Calculate Feynman diagram weight with binary tree expansion !!!!!!!!!!!!!!!
        if(NewOrder==1) then
          TotalVerWeight=SpinFactor(1, i, NewOrder)
          AbsTotalVerWeight=abs(TotalVerWeight)
        else
          SpinCache(1, 1)=NewVerWeight(VerIndex(i,1, NewOrder))
          SpinCache(2, 1)=NewVerWeight(VerIndex(i,2, NewOrder))
          BlockNum=2
          do j=2, VerNum(NewOrder)
            do k=1, BlockNum
              SpinCache(2*k-1, j)=SpinCache(k, j-1)*NewVerWeight(VerIndex(i,2*j-1, NewOrder))
              SpinCache(2*k, j)=SpinCache(k, j-1)*NewVerWeight(VerIndex(i,2*j, NewOrder))
            enddo
            BlockNum=BlockNum*2
          enddo

          TotalVerWeight=0.0
          AbsTotalVerWeight=0.0
          do j=1, 2**(NewOrder-1)
            TotalVerWeight=TotalVerWeight+SpinCache(j, VerNum(NewOrder))*SpinFactor(j, i, NewOrder)
            AbsTotalVerWeight=AbsTotalVerWeight+abs(SpinCache(j, VerNum(NewOrder))*SpinFactor(j, i, NewOrder))
          enddo

        endif

        if(SaveWeight==1) then
          DiagWeight(i) = TotalGWeight*TotalVerWeight/(2.0*pi)**(D*NewOrder)*ReWeightFactor(NewOrder)
          DiagWeightABS(i) = abs(TotalGWeight)*abs(AbsTotalVerWeight)/(2.0*pi)**(D*NewOrder)*ReWeightFactor(NewOrder)
        endif
        CalcWeight = CalcWeight + TotalGWeight*TotalVerWeight/(2.0*pi)**(D*NewOrder)*ReWeightFactor(NewOrder)
      enddo
      return
    end function CalcWeight
    
    
    subroutine IncreaseOrder()
      !increase diagram order by one/change normalization diagram to physical diagram
      implicit none
      integer :: tNum, i
      double precision :: NewTau
      double precision :: R, Weight, Kamp, dK, theta, phi, Prop
      double precision, dimension(D) :: NewMom
      if (CurrOrder==Order) return
      PropStep(1)=PropStep(1)+1.0

      ! Generate New Tau
      NewTau=grnd()*Beta

      ! Generate New Mom
      !!!! Hard way  !!!!!!!!!!!!!!!!!!!!!!!!!!
      dK=kF/sqrt(Beta)/4.0
      Kamp=kF+(grnd()-0.5)*2.0*dK
      !kF-dK<amp<kF+dK
      if(Kamp<=0.0) then
        print *, "K amplitude can not be smaller than zero!"
        STOP
      endif
      phi=2.0*pi*grnd()
      if(D==3) then
        theta=pi*grnd()
        if(theta==0.0) return
        NewMom(1)=Kamp*sin(theta)*cos(phi)
        NewMom(2)=Kamp*sin(theta)*sin(phi)
        NewMom(D)=Kamp*cos(theta)
        Prop=Beta*2.0*dK*2.0*pi*pi*sin(theta)*Kamp**(D-1)
      else if(D==2) then
        NewMom(1)=Kamp*cos(phi)
        NewMom(2)=Kamp*sin(phi)
        Prop=Beta*2.0*dK*2.0*pi*Kamp**(D-1)
      else
        print *, "Dimension ", D," has not yet been implemented!"
      endif

      !!! Simple way  !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! do i=1, D
      !   NewMom(i)=kF*(grnd()-0.5)*2
      ! enddo
      ! Prop=Beta*(2.0*kF)**D
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      tNum=TauNum(CurrOrder+1)
      TauTable(tNum-1)=NewTau
      TauTable(tNum)=NewTau
      LoopMom(:, LoopNum(CurrOrder+1))=NewMom


      call ResetWeightTable(CurrOrder+1)
      Weight = CalcWeight(0, CurrOrder+1)

      R=abs(Weight)/abs(CurrWeight)*Prop
      if(grnd()<R) then
        AcceptStep(1)=AcceptStep(1)+1.0
        CurrWeight=Weight
        CurrOrder=CurrOrder+1
        call UpdateState()
      endif
      return
    end subroutine
  
    subroutine DecreaseOrder()
    !decrease diagram order by one/change physical diagram to normalization diagram
      implicit none
      integer :: i
      double precision :: R, Weight, Prop
      double precision :: Kamp, dK, SinTheta
      if (CurrOrder==1) return
      !if the current diagrams are already in normalization sector, then return

      !Get proper K proposed probability
      !!!!!!! Hard way !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dK=kF/sqrt(Beta)/4.0
      Kamp=norm2(LoopMom(:, LoopNum(CurrOrder)))
      if(Kamp<kF-dK .or. Kamp>kF+dK) return
      if(D==3) then
        SinTheta=norm2(LoopMom(1:2, LoopNum(CurrOrder)))/Kamp
        if(SinTheta==0.0) return
        Prop=1.0/(Beta*2.0*dK*2.0*pi*pi*SinTheta*Kamp**(D-1))
      else if(D==2) then
        Prop=1.0/(Beta*2.0*dK*2.0*pi*Kamp**(D-1))
      else
        print *, "Dimension ", D," has not yet been implemented!"
      endif

      !!!!!!! Simple way !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! do i=1, D
      !   if(abs(LoopMom(i, LoopNum(CurrOrder)))>kF) return
      ! enddo
      ! Prop=1.0/(Beta*(2.0*kF)**D)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PropStep(2)=PropStep(2)+1.0
      call ResetWeightTable(CurrOrder-1)
      Weight = CalcWeight(0, CurrOrder-1)
      R=abs(Weight)/abs(CurrWeight)*Prop
      if(grnd()<R) then
        AcceptStep(2)=AcceptStep(2)+1.0
        CurrWeight=Weight
        CurrOrder=CurrOrder-1
        call UpdateState()
      endif
      return
    end subroutine
  
    subroutine ChangeTau()
    !randomly choose a vertex, change the time variable
      implicit none
      double precision :: NewTau, prop, dw, R
      double precision :: Weight, OldTau
      integer :: Num
  
      Num=int(grnd()*CurrOrder)*2+1 !1,3,5

      if(Num==1 .and. ObsType==1) return
  
      PropStep(3)=PropStep(3)+1.0
      OldTau=TauTable(Num+1) !in the case of Num==1, TauTable(1)/=TauTable(2)
      call GenerateNewTau(OldTau, NewTau, prop)
      TauTable(Num+1)=NewTau
      if(Num/=1) TauTable(Num)=NewTau !in the case of Num==1, TauTable(1)/=TauTable(2)
      call ChangeWeightTableTau(Num)
      Weight = CalcWeight(0, CurrOrder)
      R=prop*abs(Weight)/abs(CurrWeight)
      if(grnd()<R) then
        AcceptStep(3)=AcceptStep(3)+1.0
        CurrWeight=Weight
        call UpdateState()
      else
        if(Num/=1) TauTable(Num)=OldTau
        TauTable(Num+1)=OldTau
      endif
      return
    end subroutine
    
    subroutine ChangeMom()
      !randomly choose a vertex, change the space variable
      implicit none
      double precision, dimension(D) :: NewMom, OldMom
      double precision :: prop, R
      integer :: Num, i, j, NewExtMomBin
      double precision :: Weight
  
      Num = int( (LoopNum(CurrOrder))*grnd() ) + 1
  
      PropStep(4) = PropStep(4) + 1.0
      OldMom = LoopMom(:, Num)
  
      if(Num==1) then
        call GenerateNewExtMom(ExtMomBin, NewExtMomBin, prop)
        NewMom=ExtMomMesh(:, NewExtMomBin)
      else
        call GenerateNewMom(OldMom, NewMom, prop)
      endif
  
      if(prop<0.0) return
      if(Num==1 .and. norm2(NewMom)>ExtMomMax) return
      LoopMom(:,Num)=NewMom
  
      call ChangeWeightTableMom(Num)
      Weight = CalcWeight(0, CurrOrder)
      R = prop*abs(Weight)/abs(CurrWeight)
  
      if(grnd()<R) then
        AcceptStep(4) = AcceptStep(4)+1.0
        CurrWeight = Weight
        if(Num==1) ExtMomBin=NewExtMomBin
        call UpdateState()
      else
        LoopMom(:,Num)=OldMom
      endif
  
      return
    end subroutine
    
    subroutine SwapMom()
    !randomly choose a vertex, change the space variable
      implicit none
      double precision, dimension(D) :: TempMom1, TempMom2
      double precision :: R
      integer :: Num1, Num2
      double precision :: Weight
  
      Num1=int(grnd()*LoopNum(CurrOrder))+1
      Num2=int(grnd()*LoopNum(CurrOrder))+1
      if(Num1==1 .or. Num2==1) return
      if(Num1==Num2) return
  
      ! PropStep(5)=PropStep(5)+1.0
  
      TempMom1=LoopMom(:, Num1)
      TempMom2=LoopMom(:, Num2)
  
      LoopMom(:,Num1)=TempMom1
      LoopMom(:,Num2)=TempMom2
  
      call ResetWeightTable(CurrOrder)
      Weight=CalcWeight(0, CurrOrder)
      R=abs(Weight)/abs(CurrWeight)
  
      if(grnd()<R) then
      !  AcceptStep(5)=AcceptStep(5)+1.0
        CurrWeight=Weight
        call UpdateState()
      else
        LoopMom(:,Num1)=TempMom1
        LoopMom(:,Num2)=TempMom2
        !call ResetWeightTable()
      endif
  
      return
    end subroutine
    
    subroutine ChangeSpin()
      !randomly choose a loop basis, change the spin variable
      implicit none
      return
    end subroutine
    
    subroutine GenerateNewMom(old, new, prop)
      implicit none
      double precision, dimension(D) :: old, new
      integer :: Direction, Num
      double precision, parameter :: STEP=5.0
      double precision :: ratio, prop, k, k_new
      double precision :: lambda, x
      x=grnd()
      if(x<1.0/3.0) then
          new=old
          Num=int(grnd()*D)+1
          new(Num)=new(Num)+sign(STEP/Beta*grnd(), grnd()-0.5)
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
    
    
    subroutine GenerateNewExtMom(old, new, prop)
      implicit none
      double precision :: prop, x
      integer :: old, new
  
      new=int(grnd()*QBinNum)+1
      prop=1.0
    end subroutine
    
    
    subroutine GenerateNewTau(OldTau, NewTau, Prop)
      implicit none
      double precision, intent(in) :: OldTau
      double precision, intent(out) :: NewTau
      double precision, intent(out) :: Prop
      double precision :: DeltaT, x
      x=grnd()
      if(x<1.0/3.0) then
        DeltaT=Beta/3.0
        NewTau=OldTau+DeltaT*(grnd()-0.5)
      else if(x<2.0/3.0) then
        NewTau=-OldTau
      else
        NewTau=grnd()*Beta
      endif
      if (NewTau<0) then
        NewTau=NewTau+Beta
      else if (NewTau>Beta) then
        NewTau=NewTau-Beta
      endif
      Prop=1.0
      return
    end subroutine GenerateNewTau

end program main
