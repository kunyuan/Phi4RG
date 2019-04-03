module parameters
  IMPLICIT NONE
  !-- Parameters -------------------------------------------------
  integer, parameter :: D=3           !D=2 or D=3  
  integer, parameter :: kNum=512       !k bins of Green's function
  double precision, parameter :: kMax=8.0
  double precision, parameter :: DeltaK=kMax/kNum
  double precision, parameter   :: UVScale=8.0     !the upper bound of energy scale
  integer, parameter                    :: MaxOrder=4           ! Max diagram order
  integer, parameter          :: SIGMA=1, GAMMA4=2

  integer            :: PID           ! the ID of this job
  integer            :: Order
  double precision   :: Mass2         ! square mass
  double precision   :: BareCoupling         ! bare coupling
  integer            :: Seed          ! random-number seed

  !-- Markov Chain ----------------------------------------------
  double precision                        :: Step        ! a counter to keep track of the current step number
  integer, parameter                      :: UpdateNum=6 ! number of updates
  double precision, dimension(UpdateNum)  :: PropStep, AcceptStep
  double precision, dimension(0:MaxOrder, 2) :: ReWeightFactor       !reweightfactor for each order

  integer                                 :: CurrOrder
  integer                                 :: CurrType    ! SelfEnergy: 1, Ver4: 2
  integer                                 :: CurrExtMom  !external momentum for self energy
  double precision                        :: CurrWeight, CurrIRScale
  double precision, dimension(D, MaxOrder+3)  :: LoopMom ! values to attach to each loop basis
  double precision, dimension(D, kNum)        :: ExtMomMesh ! external momentum
  integer, dimension(MaxOrder, 2)         :: LoopNum !self energy has two ext loops: Order+2, gamma4: Order
  double precision, dimension(D)          :: Mom0 ! values to attach to each loop basis

  !--- Measurement ------------------------------------------------
  double precision                            :: DiffVer
  double precision, dimension(kNum)           :: DiffSigma !the subleading correction of sigma
  double precision                            :: DiffMu  !the leading order correction of sigma
  double precision                            :: Norm

  !-- Diagram Elements  --------------------------------------------
  double precision                            :: EffVer
  double precision, dimension(kNum)           :: EffSigma
  double precision                            :: EffMu

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
    integer :: iBlck, iInner, AnnealStep, TotalStep
    double precision :: x
  
    print *, 'mass2, coupling, Order, TotalStep(*1e6), Seed, PID'
    read(*,*) Mass2, BareCoupling, Order, TotalStep, Seed, PID

    ! For a given order, the bigger factor, the more accurate result 

    call sgrnd(Seed) 
  
    print *, "Initializing ..."
    call Initialize()
    print *, "Initialized!"
  
    call Test() !call test first to make sure all subroutines work as expected
  
    AnnealStep=4
    Step=0.0
    print *, "Start simulation..."
    do iBlck=1,TotalStep
      do iInner=1,1000000
        Step=Step+1.0
        x=grnd()
        if (x<1.0/UpdateNum) then
          call ChangeOrder()
        else if (x<2.0/UpdateNum) then
          call ChangeMom()
        else if (x<3.0/UpdateNum) then
          call ChangeType()
        ! else if (x<5.0/UpdateNum) then
        !   call ChangeExtMom()
        endif
        !if(mod(int(Step), 4)==0) call Measure()
        call Measure()

      enddo

      ! call DynamicTest()
      call SolveDSE()

      if (iBlck>=AnnealStep) then
        AnnealStep=int(AnnealStep*1.5)
        if(CurrIRScale>0.1) then
          CurrIRScale=CurrIRScale/1.5
        endif
      endif

      if(mod(iBlck, 10)==0) then
        !!!!!!!!!!!!  Print Info and Save Data !!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) 
        print *, iBlck, "million steps"
        write(*,"(A20)") "Accept Ratio: "
        write(*,"(A16, f8.3)") "Increase Order:", AcceptStep(1)/PropStep(1)
        write(*,"(A16, f8.3)") "Decrease Order:", AcceptStep(2)/PropStep(2)
        write(*,"(A16, f8.3)") "Change Mom:", AcceptStep(3)/PropStep(3)
        write(*,"(A16, f8.3)") "Gamma4->Sigma:", AcceptStep(4)/PropStep(4)
        write(*,"(A16, f8.3)") "Sigma->Gamma4:", AcceptStep(5)/PropStep(5)
        ! write(*,"(A16, f8.3)") "Change ExtK:", AcceptStep(6)/PropStep(6)
        ! write(*, *) "coupling: ", DiffVer(CurrScale)/Norm
        write(*,*) "IR Scale:", CurrIRScale
        write(*, *) "coupling: ", EffVer
        ! write(*, *) "coupling: ", DiffVer/Norm
        ! write(*, *) "mu: ", DiffMu/Norm
        write(*, *) "mu: ", EffMu
      endif

      if (mod(iBlck, 100)==10)  then
        call SaveToDisk()
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end do

    call SaveToDisk()
    print *, "End simulation."
  
    CONTAINS

    subroutine Initialize()
      implicit none
      integer :: i, j
      double precision :: kamp
  ! For a given order, the bigger factor, the more accurate result 
      ReWeightFactor(0:2, GAMMA4)=(/1.0,1.0,20.0/)
      ReWeightFactor(0:2, SIGMA)=(/0.0,1.0,20.0/)
      Mom0=0.0

      PropStep=1.0e-10
      AcceptStep=0.0

      DiffVer=0.0
      DiffSigma=0.0
      DiffMu=0.0
      Norm=1.0e-10

      LoopNum(1:3, SIGMA)=(/1,2,3/)
      LoopNum(1:3, GAMMA4)=(/1,2,3/)

      ExtMomMesh=0.0
      do j=1, kNum
        !1 <--> DeltaK/2, kNum <--> kMax-DeltaK/2
        kamp=(j-0.5)*DeltaK
        ExtMomMesh(1, j)=kamp
      enddo

      EffVer=BareCoupling
      EffSigma=0.0
      EffMu=Mass2

      CurrType=GAMMA4 !start with gamma4
      CurrIRScale=UVScale/2.0
      CurrOrder=0
      CurrExtMom=1
      CurrWeight=CalcWeight(CurrOrder, CurrType)
    end subroutine

    subroutine Test()
      implicit none
      return
    end subroutine

    subroutine DynamicTest()
      implicit none
      return
    end subroutine

    subroutine Measure()
      implicit none
      double precision :: Factor

      Factor=CurrWeight/abs(CurrWeight)/ReWeightFactor(CurrOrder, CurrType)
      if(CurrOrder==0) then
          Norm=Norm+Factor
      else
        if(CurrType==GAMMA4) then
            DiffVer=DiffVer+Factor
        else
          if(CurrOrder==1) then
            DiffMu=DiffMu+Factor
          else
            DiffSigma(CurrExtMom)=DiffSigma(CurrExtMom)+Factor
          endif
        endif
      endif
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

    subroutine SolveDSE()
      implicit none
      integer :: i, start, end
      double precision :: dg, dMu

      EffVer=(BareCoupling-DiffVer/Norm)*UVScale/CurrIRScale
      EffMu=DiffMu/Norm

    end subroutine

    double precision function CalcWeight(NewOrder, Type)
      !calculate the weight for ALL diagrams in a given sector
      implicit none
      integer :: NewOrder, Type
      if(NewOrder==0) then 
        CalcWeight=1.0
      else
        if(Type==GAMMA4) then
          !4-vertex
          if(NewOrder==1) then
            CalcWeight=Ver4_OneLoop(-1, 0, Mom0, Mom0, Mom0, Mom0, 1, .true.)
          endif
        else
          ! Sigma
          CalcWeight=Sigma_OneLoop(0, ExtMomMesh(:, CurrExtMom), 1)
        endif
      endif
      ! print *, CalcWeight
      CalcWeight=CalcWeight*ReWeightFactor(NewOrder, Type)
      return
    end function CalcWeight
    
    subroutine ChangeOrder()
      !increase diagram order by one/change normalization diagram to physical diagram
      implicit none
      double precision :: R, Weight, Prop
      integer :: NewOrder, Index
      if(CurrType==SIGMA) return

      if(grnd()<0.5) then
        ! increase order
        if (CurrOrder==Order) return
        Index=1
        NewOrder=CurrOrder+1
        call CreateMom(LoopMom(:, CurrOrder+1), Prop)
        if(Prop<0.0) return
      else
        if (CurrOrder==0) return
        Index=2
        NewOrder=CurrOrder-1
        call RemoveMom(LoopMom(:, CurrOrder), Prop)
        if(Prop<0.0) return
      endif

      PropStep(Index)=PropStep(Index)+1.0
      Weight = CalcWeight(NewOrder, CurrType)
      R=abs(Weight)/abs(CurrWeight)*Prop
      if(grnd()<R) then
        AcceptStep(Index)=AcceptStep(Index)+1.0
        CurrWeight=Weight
        CurrOrder=NewOrder
        CurrType=GAMMA4
      endif
      return
    end subroutine
  
    subroutine ChangeMom()
      !randomly choose a vertex, change the space variable
      implicit none
      double precision, dimension(D) :: NewMom, OldMom
      double precision :: prop, R, Weight
      integer :: Num
  
      Num = int( CurrOrder*grnd() ) + 1
  
      PropStep(3) = PropStep(3) + 1.0
      OldMom = LoopMom(:, Num)
  
      call ShiftMom(OldMom, NewMom, prop)
  
      if(prop<0.0) return
      LoopMom(:,Num)=NewMom
  
      Weight = CalcWeight(CurrOrder, CurrType)
      R = prop*abs(Weight)/abs(CurrWeight)

      ! print *, R, Weight, CurrWeight
  
      if(grnd()<R) then
        AcceptStep(3) = AcceptStep(3)+1.0
        CurrWeight = Weight
      else
        LoopMom(:,Num)=OldMom
      endif
  
      return
    end subroutine

    subroutine ChangeType()
      implicit none
      integer :: OldType, Index
      double precision :: Weight, prop, R

      if (CurrOrder==0) return

      OldType=CurrType
      if(CurrType==GAMMA4) then
        Index=4
        CurrType=SIGMA
        ! CurrExtMom=int(grnd()*ScaleTable(CurrScale)/DeltaK)+1
        ! CurrExtMom=int(ScaleNum/2)+1
        ! prop=ScaleNum*2
        CurrExtMom=1
      else
        Index=5
        CurrType=GAMMA4
        ! if(CurrExtMom>=int(ScaleTable(CurrScale)/DeltaK)+1) return
        ! if(CurrExtMom>=int(ScaleNum/2)+1) return
        ! prop=1.0/ScaleNum/2
      endif

      prop=1.0
      PropStep(Index) = PropStep(Index) + 1.0
      Weight = CalcWeight(CurrOrder, CurrType)
      R = prop*abs(Weight)/abs(CurrWeight)

      if(grnd()<R) then
        AcceptStep(Index) = AcceptStep(Index)+1.0
        CurrWeight = Weight
      else
        CurrType=OldType
      endif
      return
    end subroutine

    subroutine ChangeExtMom()
      !randomly choose a vertex, change the space variable
      implicit none
      double precision :: prop, R, Weight
      integer :: OldExtMom
      integer, parameter :: Delta=3

      if(CurrType==GAMMA4) return
  
      OldExtMom=CurrExtMom
      if(grnd()<0.5) then
        CurrExtMom=CurrExtMom+int(grnd()*Delta)
      else
        CurrExtMom=CurrExtMom-int(grnd()*Delta)
      endif
      if(CurrExtMom>kMax .or. CurrExtMom<1) then
        CurrExtMom=OldExtMom
        return
      endif
  
      PropStep(6) = PropStep(6) + 1.0
      Weight = CalcWeight(CurrOrder, CurrType)
      R = prop*abs(Weight)/abs(CurrWeight)
  
      if(grnd()<R) then
        AcceptStep(6) = AcceptStep(6)+1.0
        CurrWeight = Weight
        ! call UpdateState()
      else
        CurrExtMom=OldExtMom
      endif
  
      return
    end subroutine

    subroutine CreateMom(New, Prop)
      implicit none 
      double precision, dimension(D) :: New
      double precision :: Prop, dK, Kamp, theta, phi
      !!!! the Hard way  !!!!!!!!!!!!!!!!!!!!!!!!!!
      dK=2.0
      Kamp=grnd()*dK
      Prop=-1.0
      if(Kamp<=0.0) return
      phi=2.0*pi*grnd()
      theta=pi*grnd()
      if(theta==0.0) return
      New(1)=Kamp*sin(theta)*cos(phi)
      New(2)=Kamp*sin(theta)*sin(phi)
      New(3)=Kamp*cos(theta)
      Prop=dK*2.0*pi*pi*sin(theta)*Kamp**(D-1)

      !!! Simple way  !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! do i=1, D
      !   NewMom(i)=kF*(grnd()-0.5)*2
      ! enddo
      ! Prop=Beta*(2.0*kF)**D
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine

    subroutine RemoveMom(Old, Prop)
      implicit none
      double precision, dimension(D) :: Old
      double precision :: Prop, dK, Kamp, SinTheta
      !Get proper K proposed probability
      !!!!!!! Hard way !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dK=2.0
      Kamp=norm2(Old)
      Prop=-1.0
      if(Kamp<=0.0 .or. Kamp>=dK) return
      SinTheta=norm2(Old(1:2))/Kamp
      if(SinTheta==0.0) return
      Prop=1.0/(dK*2.0*pi*pi*SinTheta*Kamp**(D-1))

      !!!!!!! Simple way !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! do i=1, D
      !   if(abs(LoopMom(i, LoopNum(CurrOrder)))>kF) return
      ! enddo
      ! Prop=1.0/(Beta*(2.0*kF)**D)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine

    
    subroutine ShiftMom(old, new, prop)
      implicit none
      double precision, dimension(D) :: old, new
      integer :: Num
      double precision :: ratio, prop, k, k_new
      double precision :: lambda
      x=grnd()
      if(x<1.0/3.0) then
          new=old
          Num=int(grnd()*D)+1
          new(Num)=new(Num)+sign(grnd(), grnd()-0.5)
          prop=1.0
      else if(x<2.0/3.0) then
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

    double precision function Green(Mom, IRScale, g_type)
    !dimensionless green's function
      implicit none
      double precision :: kk, gg, IRScale
      integer :: g_type, kamp
      double precision, dimension(D) :: Mom

      kk=norm2(Mom)

      ! if(kk>UVScale) then
      !   Green=0.0
      !   return
      ! endif

      if(kk<kMax) then
        !1 <--> DeltaK/2, kNum <--> kMax-DeltaK/2
        kamp=int(kk/DeltaK)+1
        ! gg=1.0/(kk*kk+EffMu(Scale)+1.0+EffSigma(kamp, Scale))
        gg=1.0/(kk*kk+1.0+EffSigma(kamp))
      else
        gg=1.0/(kk*kk+1.0+EffSigma(kamp))
      endif
      

      if(g_type==0) then
        Green=gg
      else 
        Green=-2.0*gg*gg !dG/dLn\lambda
      endif

      ! if(g_type==0) then
      !   Green=1.0/(kk*kk+1.0)
      ! else 
      !   Green=-2.0/(kk*kk+1.0)/(kk*kk+1.0) !dG/dLn\lambda
      ! endif

      return
    end function
    
    double precision function Ver4(VerType)
      implicit none
      integer :: VerType
      ! Ver4=BareCoupling
      Ver4=EffVer
      return
    end function Ver4

    double precision function Ver4_OneLoop(LOrder, ROrder, MomL1, MomL2, MomR1, MomR2, InterMomIndex, Simple)
      implicit none
      double precision, dimension(D) :: MomL1, MomL2, MomR1, MomR2, Mom
      integer :: LOrder, ROrder !order of left and right ver4
      integer :: InterMomIndex
      logical :: Simple
      double precision :: LWeight, RWeight, UWeight, SWeight, TWeight
      ! if(Simple==.true.)
      if(LOrder==-1) then
        LWeight=BareCoupling  
      else if(LOrder==0) then
        LWeight=EffVer
      endif
      if(ROrder==-1) then
        RWeight=BareCoupling  
      else if(ROrder==0) then
        RWeight=EffVer
      endif
      Mom=LoopMom(:, InterMomIndex)
      UWeight=LWeight*RWeight*Green(MomL1-MomL2+Mom, CurrIRScale, 0)*Green(Mom, CurrIRScale, 0)
      if(Simple .eqv. .true.) then
        Ver4_OneLoop=UWeight*3.0/(2.0*pi)**D/2.0
        return
      endif
    end function

    double precision function Sigma_OneLoop(VerOrder, ExtK, InterMomIndex)
      implicit none
      integer :: VerOrder, InterMomIndex
      double precision, dimension (D) :: ExtK, Mom
      double precision :: VerWeight, Weight
      Sigma_OneLoop=0.0
      ! if(VerOrder==0) VerWeight=EffVer
      ! Mom=LoopMom(:, InterMomIndex)
      ! Weight=0.5*Green(Mom, CurrIRScale, 0)*VerWeight
      ! Sigma_OneLoop=Weight/(2.0*pi)**D
      return
    end function
end program main
