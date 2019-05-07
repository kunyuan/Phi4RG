module parameters
  IMPLICIT NONE
  !-- Parameters -------------------------------------------------
  integer, parameter :: D=3           !D=2 or D=3  
  integer, parameter :: ScaleNum=8    !number of scales
  integer, parameter :: kNum=8       !k bins of Green's function
  double precision, parameter :: kMax=8.0 !the maximum external K=\Lambda*8.0
  double precision, parameter :: UVScale=8.0     !the upper bound of energy scale
  double precision, parameter :: DeltaK=kMax/kNum
  integer, parameter          :: MaxOrder=4           ! Max diagram order
  integer, parameter          :: SIGMA=1, GAMMA4=2, ETA=3
  logical, parameter          :: IsCritical=.true.

  integer            :: PID           ! the ID of this job
  integer            :: Order
  double precision   :: Mass2         ! square mass
  double precision   :: BareCoupling         ! bare coupling
  integer            :: Seed          ! random-number seed

  !-- Markov Chain ----------------------------------------------
  double precision                        :: Step        ! a counter to keep track of the current step number
  integer, parameter                      :: UpdateNum=9 ! number of updates
  double precision, dimension(UpdateNum)  :: PropStep, AcceptStep
  double precision, dimension(0:MaxOrder, 3) :: ReWeightFactor       !reweightfactor for each order

  integer                                 :: CurrOrder, CurrScale, CurrIRScale
  integer                                 :: CurrType    !SelfEnergy: 1, Ver4: 2
  integer                                 :: CurrExtMom  !external momentum for self energy
  double precision                        :: CurrWeight
  double precision, dimension(D, MaxOrder+3)  :: LoopMom ! values to attach to each loop basis
  double precision, dimension(D, kNum)  :: ExtMomMesh ! external momentum 
  integer, dimension(MaxOrder, 2)         :: LoopNum !self energy has two ext loops: Order+2, gamma4: Order
  double precision, dimension(D)          :: Mom0 ! values to attach to each loop basis

  !--- Measurement ------------------------------------------------
  double precision, dimension(ScaleNum)       :: ScaleTable, dScaleTable
  double precision, dimension(ScaleNum, MaxOrder)       :: DiffVer
  double precision, dimension(kNum, ScaleNum, MaxOrder) :: DiffSigma !the subleading correction of sigma
  double precision, dimension(kNum, ScaleNum, MaxOrder) :: DiffGamma2K 
  double precision, dimension(ScaleNum)       :: Norm

  !-- Diagram Elements  --------------------------------------------
  double precision, dimension(ScaleNum)       :: EffVer
  double precision, dimension(kNum, ScaleNum) :: EffSigma
  double precision, dimension(kNum, ScaleNum) :: EffGamma2K   !K*d \Gamma^2(K, \Lambda)/dK
  double precision, dimension(ScaleNum)       :: EffMu

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
    integer :: iBlck, iInner, AnnealStep, TotalStep, scale, i
    double precision :: x, kamp
  
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
        x=grnd()*(UpdateNum-1)
        if (x<1.0) then
          call ChangeOrder()
        else if (x<2.0) then
          call ChangeScale()
        else if (x<3.0) then
          call ChangeMom()
        else if (x<4.0) then
          call ChangeExtMom()
        else if (x<6.0) then
          call ChangeType()
        else if (x<8.0) then
          call ChangeType2()
        endif
        !if(mod(int(Step), 4)==0) call Measure()
        call Measure()

      enddo

      ! call DynamicTest()
      call SolveBetaFunc()

      if (iBlck==AnnealStep) then
        AnnealStep=AnnealStep*2
        CurrIRScale=CurrIRScale/2
      endif

      if(mod(iBlck, 10)==0) then
        !!!!!!!!!!!!  Print Info and Save Data !!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) 
        write(*, *) "Differ Order 1: ", DiffVer(:, 1)/Norm
        ! write(*, *) "Differ Order 2: ", DiffVer(:, 2)/Norm
        ! write(*, *) "coupling: ", EffVer
        ! scale=ScaleNum/2+1
        scale=1
        ! write(*,*) "DiffGamma2: "
        ! do i=1, kNum
        !     print *, DiffGamma2K(i, scale, 2)*Green(ExtMomMesh(:, i)/ScaleTable(scale), scale, 0)
        ! enddo
        write(*, *) "DiffGamma2: ", DiffGamma2K(:, scale, 2)/Norm(scale)

        write(*,*) "EffGammaK2"
        do i=1, kNum
          kamp=norm2(ExtMomMesh(:, i))
          ! write(*, *) kamp,  2.0-EffGamma2K(i, scale)*kamp*Green(ExtMomMesh(:,i), scale, 0)
        enddo 
        write(*, *) "self energy: ", DiffSigma(:, scale, 2)/Norm(scale)
        ! write(*, *) "coupling: ", DiffVer/Norm
        ! write(*, *) "mu: ", DiffMu/Norm
        ! write(*, *) "mu: ", EffMu

        print *, iBlck, "million steps"
        write(*,"(A20)") "Accept Ratio: "
        write(*,"(A16, f8.6)") "Increase Order:", AcceptStep(1)/PropStep(1)
        write(*,"(A16, f8.6)") "Decrease Order:", AcceptStep(2)/PropStep(2)
        write(*,"(A16, f8.6)") "Change Scale:", AcceptStep(3)/PropStep(3)
        write(*,"(A16, f8.6)") "Change Mom:", AcceptStep(4)/PropStep(4)
        write(*,"(A16, f8.6)") "Change ExtK:", AcceptStep(5)/PropStep(5)
        write(*,"(A16, f8.6)") "Gamma4->Sigma:", AcceptStep(6)/PropStep(6)
        write(*,"(A16, f8.6)") "Sigma->Gamma4:", AcceptStep(7)/PropStep(7)
        write(*,"(A16, f8.6)") "Gamma4->Eta:", AcceptStep(8)/PropStep(8)
        write(*,"(A16, f8.6)") "Eta->Gamma4:", AcceptStep(9)/PropStep(9)

        write(*,*) "Current IRScale: ", CurrIRScale, "; Current Scale: ", CurrScale
        write(*,*) "Current Order: ", CurrOrder
        write(*,*) "Current Type: ", CurrType
        write(*,*) "Current Weight: ", CurrWeight
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
      integer :: i, j, k
      double precision :: kamp
  ! For a given order, the bigger factor, the more accurate result 
      ReWeightFactor(0:2, GAMMA4)=(/1.0,1.0,0.20/)
      ReWeightFactor(0:2, SIGMA)=(/0.0,1.0,20.00/)
      ReWeightFactor(0:2, ETA)=(/0.0,0.0,0.20/)
      Mom0=0.0

      PropStep=0.0
      AcceptStep=0.0

      DiffVer=0.0
      DiffSigma=0.0
      DiffGamma2K=0.0
      Norm=1.0e-10

      LoopNum(1:3, SIGMA)=(/1,2,3/)
      LoopNum(1:3, GAMMA4)=(/1,2,3/)

      do i=1, ScaleNum
        ScaleTable(i)=i*1.0/ScaleNum*UVScale

      end do

      ExtMomMesh=0.0
      do j=1, kNum
        !1 <--> DeltaK/2, kNum <--> kMax-DeltaK/2
        ! kamp=(j-0.5)*DeltaK
        kamp=(j-1)*DeltaK
        ExtMomMesh(1, j)=kamp
        do k=1, ScaleNum
          DiffGamma2K(j, k, 1)=2.0*kamp
          EffGamma2K(j, k)=DiffGamma2K(j, k, 1)
        enddo
      enddo

      EffVer=BareCoupling
      EffSigma=0.0
      if(IsCritical) Mass2=0.0
      EffMu=Mass2

      do i=1, ScaleNum-1
        dScaleTable(i)=ScaleTable(i+1)-ScaleTable(i)
      end do


      CurrType=GAMMA4 !start with gamma4
      CurrScale=ScaleNum
      CurrIRScale=ScaleNum/2
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

    subroutine SaveToDisk()
      implicit none
      integer :: i, ref, j, o, scale
      double precision :: Obs
      double precision :: DeltaQ
      double precision, dimension(D) :: kk
      double precision :: kamp
      character*10 :: ID
      character*10 :: order_str
      character*10 :: DiagIndex
      character*20 :: filename
      !Save Polarization to disk

      scale=1

      do o=1, Order
        write(ID, '(i10)') PID
        write(order_str, '(i10)') o
        filename="Sigma"//trim(adjustl(order_str))//"_"//trim(adjustl(ID))//".dat"
        write(*,*) "Save to disk ..."
        open(100, status="replace", file=trim(filename))
        write(100, *) "#", Step
        ! write(100, *) "#", DiffSigma(1, o)
        do i=1, kNum
            Obs = DiffSigma(i, scale, o)
            write(100, *) norm2(ExtMomMesh(:, i)), Obs
        enddo
        close(100)
      enddo

      filename="Eta"//"_"//trim(adjustl(ID))//".dat"
      write(*,*) "Save to disk ..."
      open(100, status="replace", file=trim(filename))
      write(100, *) "#", Step, ScaleTable(scale)
      ! write(100, *) "#", DiffSigma(1, o)
      do i=1, kNum
        kk=ExtMomMesh(:, i)
        kamp=norm2(kk)
        Obs = 2.0-EffGamma2K(i, scale)*Green(kk, scale, 0)*kamp
        ! Obs = 2.0-EffGamma2K(i, scale)*kamp/kamp/kamp
        write(100, *) norm2(ExtMomMesh(:, i)), Obs
      enddo
      close(100)

      return
    end subroutine

    subroutine SolveBetaFunc()
      implicit none
      integer :: i, j, start, end, o
      double precision :: dg, dMu, kk, kamp
      double precision, dimension(kNum) :: dSigma
      do i=1, ScaleNum-1
        start=ScaleNum-i+1
        end=start-1
        ! print *, start, end, dScaleTable(end), ScaleTable(start)
        !!!!!!!!!!!!!!!!  4-Vertex Renormalization  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! EffVer(end)=(EffVer(start)*ScaleTable(start)+dScaleTable(start)*DiffVer(start)/Norm)/ScaleTable(end)
        ! dg=(EffVer(start)-sum(DiffVer(start, :))/Norm(start))*dScaleTable(end)/ScaleTable(start)
        dg=(EffVer(start)-DiffVer(start, 1)/Norm(start))*dScaleTable(end)/ScaleTable(start)
        if(dg<5.0) then
          ! EffVer(end)=EffVer(start)+dg
        ! else
        !   print *, "dg", dg, start, EffVer(start)
        endif

        !!!!!!!!!!!!!!!!!!!  Sigma Renormalization  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! dMu=(2.0*EffMu(start)+DiffMu(start)/Norm(start))*dScaleTable(end)/ScaleTable(start)
        ! print *, EffMu(start), dMu
        ! EffMu(end)=EffMu(start)+dMu

        ! do o=1, Order
        !   dSigma(:)=(2.0*EffSigma(:, start)+DiffSigma(:, start, o)/Norm(start))*dScaleTable(end)/ScaleTable(start)
        ! enddo

      enddo

      do i=1, ScaleNum
        do j=1, kNum
          ! kamp=(j-0.5)*DeltaK
          kamp=(j-1)*DeltaK
          EffGamma2K(j, i)=2.0*kamp+DiffGamma2K(j, i, 2)/Norm(i)
          ! EffGamma2K(j, start)=2.0*(kamp/ScaleTable(start))
          ! if (start==2) then 
          !   print *, EffGamma2K(j, start)-2.0*(kamp/ScaleTable(start))
          ! endif
        enddo
      enddo
      !set Mu(\Lambda)=\Sigma(k=0, \Lambda)
      ! EffMu(:)=EffSigma(1,:)
    end subroutine

    subroutine Measure()
      implicit none
      double precision :: Factor, kk
      integer :: kamp

      if(CurrIRScale>=CurrScale) return

      kk=norm2(ExtMomMesh(:, CurrExtMom))

      ! if(kk<kMax) then
      !1 <--> DeltaK/2, kNum <--> kMax-DeltaK/2
      ! kamp=int(kk/DeltaK)+1

      Factor=CurrWeight/abs(CurrWeight)/ReWeightFactor(CurrOrder, CurrType)
      if(CurrOrder==0) then
          Norm(CurrScale)=Norm(CurrScale)+Factor
      else
        if(CurrType==GAMMA4) then
            DiffVer(CurrScale, CurrOrder)=DiffVer(CurrScale, CurrOrder)+Factor
        else if(CurrType==SIGMA) then
            DiffSigma(CurrExtMom, CurrScale, CurrOrder)=DiffSigma(CurrExtMom, CurrScale, CurrOrder)+Factor
        else if(CurrOrder>=2 .and. CurrType==ETA) then
            DiffGamma2K(CurrExtMom, CurrScale, CurrOrder)=DiffGamma2K(CurrExtMom, CurrScale, CurrOrder)+Factor
        endif

      endif

      return
    end subroutine
    

    double precision function CalcWeight(NewOrder, Type)
      !calculate the weight for ALL diagrams in a given sector
      implicit none
      integer :: NewOrder, Type
      double precision, dimension(D) :: ExtKScaled
      if(NewOrder==0) then 
        CalcWeight=1.0
      else
        if(Type==GAMMA4) then
          !4-vertex
          if(NewOrder==1) then
            CalcWeight=Ver4_OneLoop(0, 0, Mom0, Mom0, Mom0, Mom0, 1, .true.)
          else if(NewOrder==2) then
            CalcWeight=Ver4_TwoLoop(1)
          endif
        else if(Type==SIGMA) then
          ExtKScaled=ExtMomMesh(:, CurrExtMom)
          ! Sigma
          if(NewOrder==1) then
            CalcWeight=Sigma_OneLoop(0, ExtKScaled, 1)
          else if(NewOrder==2) then
            CalcWeight=Sigma_TwoLoop(ExtKScaled, 1)
          endif
        else if(Type==ETA) then
          ExtKScaled=ExtMomMesh(:, CurrExtMom)
          if(NewOrder==2) then
            CalcWeight=SigmaK_TwoLoop(ExtKScaled, 1)
          else
            CalcWeight=0.0
          endif
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
      if(CurrType/=GAMMA4) return

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
  
    subroutine ChangeScale()
      implicit none
      !TODO: don't forget to change all LoopMom with scale!
      integer :: OldScale
      double precision :: Weight
      double precision :: R

      OldScale=CurrScale
      if(grnd()<0.5) then
        CurrScale=CurrScale-1
      else
        CurrScale=CurrScale+1
      endif

      if(CurrScale<1 .or. CurrScale>ScaleNum) then
        CurrScale=OldScale
        return
      endif

      PropStep(3) = PropStep(3) + 1.0

      Weight = CalcWeight(CurrOrder, CurrType)
      R = abs(Weight)/abs(CurrWeight)
      if(grnd()<R) then
        AcceptStep(3) = AcceptStep(3)+1.0
        CurrWeight = Weight
      else
        CurrScale=OldScale
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
  
      PropStep(4) = PropStep(4) + 1.0
      OldMom = LoopMom(:, Num)
  
      call ShiftMom(OldMom, NewMom, prop)
  
      if(prop<0.0) return
      LoopMom(:,Num)=NewMom
  
      Weight = CalcWeight(CurrOrder, CurrType)
      R = prop*abs(Weight)/abs(CurrWeight)
  
      if(grnd()<R) then
        AcceptStep(4) = AcceptStep(4)+1.0
        CurrWeight = Weight
      else
        LoopMom(:,Num)=OldMom
      endif
  
      return
    end subroutine

    subroutine ChangeExtMom()
      !randomly choose a vertex, change the space variable
      implicit none
      double precision :: prop, R, Weight
      integer :: OldExtMom
      integer, parameter :: Delta=1

      if(CurrType==GAMMA4) return

      prop=1.0
  
      OldExtMom=CurrExtMom
      if(grnd()<0.5) then
        ! CurrExtMom=CurrExtMom+int(grnd()*Delta)
        CurrExtMom=CurrExtMom+Delta
      else
        ! CurrExtMom=CurrExtMom-int(grnd()*Delta)
        CurrExtMom=CurrExtMom-Delta
      endif
      if(CurrExtMom>kNum .or. CurrExtMom<1) then
        CurrExtMom=OldExtMom
        return
      endif
  
      PropStep(5) = PropStep(5) + 1.0
      Weight = CalcWeight(CurrOrder, CurrType)
      R = prop*abs(Weight)/abs(CurrWeight)

      ! print *, CurrExtMom
      ! print *, R, prop, Weight, CurrWeight
  
      if(grnd()<R) then
        AcceptStep(5) = AcceptStep(5)+1.0
        CurrWeight = Weight
        ! call UpdateState()
      else
        CurrExtMom=OldExtMom
      endif
  
      return
    end subroutine

    subroutine ChangeType()
      implicit none
      integer :: OldType, Index
      double precision :: Weight, prop, R

      if (CurrOrder==0) return

      if(CurrType/=GAMMA4 .and. CurrType/=SIGMA) return

      OldType=CurrType
      if(CurrType==GAMMA4) then
        Index=6
        CurrType=SIGMA
        CurrExtMom=int(grnd()*kNum)+1
        prop=kNum
      else
        Index=7
        CurrType=GAMMA4
        prop=1.0/kNum
      endif

      prop=1.0
      PropStep(Index) = PropStep(Index) + 1.0
      Weight = CalcWeight(CurrOrder, CurrType)
      R = prop*abs(Weight)/abs(CurrWeight)

      ! if(Index==5)  then
      !   print *, Weight, CurrWeight, R
      !   print *, 
      !   stop
      ! endif
  
      if(grnd()<R) then
        AcceptStep(Index) = AcceptStep(Index)+1.0
        CurrWeight = Weight
      else
        CurrType=OldType
      endif
      return
    end subroutine

    subroutine ChangeType2()
      implicit none
      integer :: OldType, Index
      double precision :: Weight, prop, R

      if (CurrOrder<2) return
      if(CurrType/=GAMMA4 .and. CurrType/=ETA) return

      OldType=CurrType
      if(CurrType==GAMMA4) then
        Index=8
        CurrType=ETA

        CurrExtMom=int(grnd()*kNum)+1
        prop=kNum
        ! CurrExtMom=1
        ! prop=1.0
      else
        Index=9
        CurrType=GAMMA4
        prop=1.0/kNum
        ! prop=1.0
      endif

      PropStep(Index) = PropStep(Index) + 1.0
      Weight = CalcWeight(CurrOrder, CurrType)
      R = prop*abs(Weight)/abs(CurrWeight)
      ! print *, "Step:", Step
      ! print *, "Weight: ", Weight, CurrWeight
      ! print *, "R: ", R
      ! print *, "Type: ", CurrType
      ! print *, "ExtMom: ", CurrExtMom, norm2(ExtMomMesh(:, CurrExtMom))
      ! print *, "Scale: ", CurrScale
      ! print *, " "
      ! if(Step>31840) stop

      if(grnd()<R) then
        AcceptStep(Index) = AcceptStep(Index)+1.0
        CurrWeight = Weight
      else
        CurrType=OldType
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

    double precision function Green(Mom, Scale, g_type)
    !dimensionless green's function
      implicit none
      double precision :: kk, gg, realKK
      integer :: g_type, Scale, kamp
      double precision, dimension(D) :: Mom
      kk=norm2(Mom)
      ! realKK=kk*ScaleTable(Scale)

      ! if(realKK>UVScale) then
      !   Green=0.0
      !   return
      ! endif

      ! if(kk<kMax) then
      !1 <--> DeltaK/2, kNum <--> kMax-DeltaK/2
      ! kamp=int(realKK/DeltaK)+1
      ! gg=1.0/(kk*kk+EffMu(Scale)+1.0+EffSigma(kamp, Scale))
      ! gg=1.0/(kk*kk+1.0+EffSigma(kamp, Scale)-EffMu(Scale))
      gg=1.0/(kk*kk+1.0)

      ! else
        ! gg=1.0/(kk*kk+1.0+EffSigma(kamp, ScaleNum))
      ! endif
      

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
    
    double precision function Ver4(VerType, Scale)
      implicit none
      integer :: VerType, Scale
      ! Ver4=BareCoupling
      Ver4=EffVer(Scale)
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
      if(LOrder==0) LWeight=EffVer(CurrScale)
      if(ROrder==0) RWeight=EffVer(CurrScale)
      Mom=LoopMom(:, InterMomIndex)
      UWeight=LWeight*RWeight*Green(MomL1-MomL2+Mom, CurrScale, 1)*Green(Mom, CurrScale, 0)
      if(Simple .eqv. .true.) then
        Ver4_OneLoop=-3.0*UWeight/(2.0*pi)**D
        return
      endif
    end function

    double precision function Ver4_TwoLoop(InterMomIndex)
      implicit none
      double precision, dimension(D) :: Mom1, Mom2
      integer :: InterMomIndex
      double precision :: UWeight, G1p, G2, G3, G4
      Mom1=LoopMom(:, InterMomIndex)
      Mom2=LoopMom(:, InterMomIndex+1)
      G1p=Green(Mom1, CurrScale, 1)  !dG/d\Lambda
      G2=Green(Mom1, CurrScale, 0)
      G3=Green(Mom2, CurrScale, 0)
      G4=Green(Mom1+Mom2, CurrScale, 0)
      UWeight=G1p*G3*(G4*G2+G4*G3-G2*G3)*EffVer(CurrScale)**3
      ! UWeight=G1p*G3*(-G2*G3)*EffVer(CurrScale)**3
      Ver4_TwoLoop=6.0*UWeight/(2.0*pi)**(2.0*D)
      if(abs(Ver4_TwoLoop)>1.0e4) then
        print *, G1p, G2, G3, G4
        print *, Step, Ver4_TwoLoop
        stop
      endif
    end function

    double precision function Sigma_OneLoop(VerOrder, ExtK, InterMomIndex)
      implicit none
      integer :: VerOrder, InterMomIndex
      double precision, dimension (D) :: ExtK, Mom
      double precision :: VerWeight, Weight
      if(IsCritical .eqv. .true.)  then
        Sigma_OneLoop=0.0
      else
        if(VerOrder==0) VerWeight=EffVer(CurrScale)
        Mom=LoopMom(:, InterMomIndex)
        Weight=0.5*Green(Mom, CurrScale, 1)*VerWeight
        Sigma_OneLoop=Weight/(2.0*pi)**D
      endif
      return
    end function

    ! dSigma/d\Lambda
    ! double precision function Sigma_TwoLoop(ExtK, InterMomIndex)
    !   implicit none
    !   integer :: InterMomIndex
    !   double precision, dimension (D) :: ExtK, Mom1, Mom2
    !   double precision :: Weight, G1, G2, G3
    !   Mom1=LoopMom(:, InterMomIndex)
    !   Mom2=LoopMom(:, InterMomIndex+1)
    !   G1=Green(Mom1, CurrScale, 1)
    !   G2=Green(Mom2, CurrScale, 0)
    !   G3=Green(Mom1+Mom2-ExtK, CurrScale, 0)
    !   Weight=-0.5*G1*(G2*G3-G2*G2)
    !   !if IsCritical is true, then one should subtract Sigma(ExtK=0).
    !   !however, it is easy to see Sigma_TwoLoop(ExtK=0) is automatically zero!
    !   Sigma_TwoLoop=Weight/(2.0*pi)**(2*D)
    ! end function

    double precision function Sigma_TwoLoop(ExtK, InterMomIndex)
      implicit none
      integer :: InterMomIndex
      double precision, dimension (D) :: ExtK, Mom1, Mom2
      double precision :: Weight, G1, G2, G3, G4
      Mom1=LoopMom(:, InterMomIndex)
      Mom2=LoopMom(:, InterMomIndex+1)
      G1=Green(Mom1, CurrScale, 0)
      G2=Green(Mom2, CurrScale, 0)
      G3=Green(Mom1+Mom2+ExtK, CurrScale, 0)
      G4=Green(Mom1+Mom2, CurrScale, 0)
      Weight=-G1*G2*G3
      if(IsCritical) then
      !if IsCritical is true, then one should subtract Sigma(ExtK=0).
        Weight=Weight+G1*G2*G4
      endif
      Sigma_TwoLoop=Weight/(2.0*pi)**(2*D)/6.0*EffVer(CurrScale)**2
    end function

    double precision function SigmaK_TwoLoop(ExtK, InterMomIndex)
      implicit none
      integer :: InterMomIndex, kamp
      double precision, dimension (D) :: ExtK, Mom1, Mom2
      double precision :: Weight, G1, G2, G3, G4, kk, realKK
      Mom1=LoopMom(:, InterMomIndex)
      Mom2=LoopMom(:, InterMomIndex+1)
      G1=Green(Mom1, CurrScale, 0)
      G2=Green(Mom1+Mom2, CurrScale, 0)
      G3=Green(Mom2+ExtK, CurrScale, 0)
      G4=Green(Mom2, CurrScale, 0)
      ! G3=Green(Mom2, CurrScale, 0)
      ! G4=Green(ExtK, CurrScale, 0)
      ! G4=Green(Mom2, CurrScale, 0)

      Weight=G1*G2*G3*G3/6.0*EffVer(CurrScale)**2
      ! Weight=G1*G2*G3*G3/6.0

      !if IsCritical is true, then one should subtract Sigma(ExtK=0).
      !however, it is easy to see Sigma_TwoLoop(ExtK=0) is automatically zero!
      SigmaK_TwoLoop=Weight/(2.0*pi)**(2*D)
      ! if(abs(Step-31287)<1.0e-3) then
      ! print *, "Sigma Step:", Step
      ! print *, "Sigma Mom:", norm2(Mom1), norm2(Mom2)
      ! print *, "G: ", G1, G2, G3
      ! endif

      kk=ExtK(1)+Mom2(1)
      !1 <--> DeltaK/2, kNum <--> kMax-DeltaK/2
      ! kamp=int(abs(kk)/DeltaK)+1
      ! if(kamp<=kNum .and. abs(kk)>1.0e-10) then
      !   SigmaK_TwoLoop=SigmaK_TwoLoop*EffGamma2K(kamp, CurrScale)*kk/abs(kk)
      ! else
      !   SigmaK_TwoLoop=SigmaK_TwoLoop*2.0*kk
      ! endif
      SigmaK_TwoLoop=SigmaK_TwoLoop*2.0*kk*ExtK(1)

      ! SigmaK_TwoLoop=0.0
      Weight=G1*G2*(G3*G3-G4*G4)/2.0*EffVer(CurrScale)**2
      ! Weight=G1*(G2*G2)*(G3-G4)/2.0*EffVer(CurrScale)**2
      SigmaK_TwoLoop=SigmaK_TwoLoop+Weight/(2.0*pi)**(2*D)*2.0
    end function
end program main
