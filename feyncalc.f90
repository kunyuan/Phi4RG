module parameters
  IMPLICIT NONE

  !-- Parameters -------------------------------------------------
  integer, parameter :: D=3           !D=2 or D=3  
  integer, parameter :: ScaleNum=64    !number of q bins of the external momentum
  double precision, parameter   :: UVScale=8.0     !the upper bound of energy scale
  integer, parameter                    :: MaxOrder=4           ! Max diagram order

  integer            :: PID           ! the ID of this job
  integer            :: Order
  double precision   :: BareCoupling         ! bare coupling
  integer            :: Seed          ! random-number seed

  !-- Markov Chain ----------------------------------------------
  double precision                       :: Step        ! a counter to keep track of the current step number
  integer, parameter                     :: UpdateNum=4 ! number of updates
  double precision, dimension(UpdateNum) :: PropStep
  double precision, dimension(UpdateNum) :: AcceptStep
  double precision, dimension(0:MaxOrder) :: ReWeightFactor       !reweightfactor for each order

  integer                               :: CurrOrder
  integer                               :: CurrScale
  double precision                      :: CurrWeight
  integer                               :: CurrIRScale
  double precision, dimension(D, MaxOrder+1)  :: LoopMom ! values to attach to each loop basis
  double precision, dimension(D)  :: Mom0 ! values to attach to each loop basis

  !--- Measurement ------------------------------------------------
  double precision, dimension(ScaleNum)       :: DiffVer(ScaleNum)
  double precision, dimension(ScaleNum)       :: Norm(ScaleNum)

  !-- Diagram Tables  --------------------------------------------
  double precision, dimension(ScaleNum)       :: EffVer(ScaleNum)
  double precision, dimension(ScaleNum)       :: ScaleTable(ScaleNum)
  double precision, dimension(ScaleNum)       :: dScaleTable(ScaleNum)

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
    double precision :: x, ClearStep, iBlck
  
    print *, 'BareCoupling, Order, TotalStep(*1e6), Seed, PID'
    read(*,*)  BareCoupling, Order, TotalStep, Seed, PID

    ! For a given order, the bigger factor, the more accurate result 

    call sgrnd(Seed) 
  
    print *, "Initializing ..."
    call Initialize()
    print *, "Initialized!"
  
    call Test() !call test first to make sure all subroutines work as expected
  
    TotalStep=TotalStep*1.0e6
    ClearStep=4
    iBlck=0
    print *, "Start simulation..."
    do while (Step<TotalStep)
      Step=Step+1.0
      x=grnd()
      if (x<1.0/UpdateNum) then
        call IncreaseOrder()
      else if (x<2.0/UpdateNum) then
        call DecreaseOrder()
      else if (x<3.0/UpdateNum) then
        call ChangeScale()
      else if (x<4.0/UpdateNum) then
        call ChangeMom()
      endif
      !if(mod(int(Step), 4)==0) call Measure()
      call Measure()

      ! call DynamicTest()

      PrintCounter=PrintCounter+1
      if (PrintCounter==1e6)  then
        iBlck=iBlck+1
        call SolveBetaFunc()
        write(*,*) 
        write(*,"(f8.2, A15)") Step/1e6, "million steps"
        write(*,"(A20)") "Accept Ratio: "
        write(*,"(A16, f8.3)") "Increase Order:", AcceptStep(1)/PropStep(1)
        write(*,"(A16, f8.3)") "Decrease Order:", AcceptStep(2)/PropStep(2)
        write(*,"(A16, f8.3)") "Change Scale:", AcceptStep(3)/PropStep(3)
        write(*,"(A16, f8.3)") "Change Mom:", AcceptStep(4)/PropStep(4)
        ! write(*, *) "coupling: ", DiffVer(CurrScale)/Norm
        write(*, *) "coupling: ", DiffVer/Norm
        write(*, *) "coupling: ", EffVer
      endif
      if (PrintCounter==1e7)  then
        ! DiffVer=0.0
        ! Norm=1.0e-10
        call SaveToDisk()
        PrintCounter=0
      endif

      if (abs(iBlck-ClearStep)<1.0e-4) then
        ! DiffVer=0.0
        ! Norm=1.0e-10
        ClearStep=ClearStep*2
        CurrIRScale=CurrIRScale/2
      endif

    end do

    call SaveToDisk()
    print *, "End simulation."
  
    CONTAINS

    subroutine Initialize()
      implicit none
      integer :: i
  ! For a given order, the bigger factor, the more accurate result 
      ReWeightFactor(0:2)=(/1.0,10.0,20.0/)
      Mom0=0.0

      PropStep=0.0
      AcceptStep=0.0

      Step=0.0
      PrintCounter=0
      SaveCounter=0

      DiffVer=0.0
      Norm=1.0e-10

      do i=1, ScaleNum
        ScaleTable(i)=i*1.0/ScaleNum*UVScale
        EffVer(i)=BareCoupling
      end do

      do i=1, ScaleNum-1
        dScaleTable(i)=ScaleTable(i+1)-ScaleTable(i)
      end do

      CurrScale=ScaleNum
      CurrOrder=0
      CurrWeight=CalcWeight(CurrOrder)
      CurrIRScale=ScaleNum/2
    end subroutine

    subroutine Test()
      implicit none
      return
    end subroutine

    subroutine DynamicTest()
      implicit none
      return
    end subroutine

    double precision function Green(Mom, Scale, g_type)
    !dimensionless green's function
      implicit none
      double precision :: k2
      integer :: g_type, Scale
      double precision, dimension(D) :: Mom
      k2=sum(Mom**2)+1.0
      
      ! if(k2>UVScale/ScaleTable(Scale)) then
      !   Green=0.0
      !   return
      ! endif

      if(g_type==0) then
        Green=1.0/k2 
      else 
        Green=-2.0/k2/k2 !dG/dLn\lambda
      endif
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

    subroutine Measure()
      implicit none
      double precision :: Factor

      Factor=CurrWeight/abs(CurrWeight)/ReWeightFactor(CurrOrder)
  
      if(CurrIRScale<CurrScale) then
        if(CurrOrder==0) then
          Norm(CurrScale)=Norm(CurrScale)+Factor
        else
          DiffVer(CurrScale)=DiffVer(CurrScale)+Factor
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

    subroutine SolveBetaFunc()
      implicit none
      integer :: i, start, end
      double precision :: dg
      do i=1, ScaleNum-1
        start=ScaleNum-i+1
        end=start-1
        ! print *, start, end, dScaleTable(end), ScaleTable(start)
        ! EffVer(end)=(EffVer(start)*ScaleTable(start)+dScaleTable(start)*DiffVer(start)/Norm)/ScaleTable(end)
        dg=(-EffVer(start)-DiffVer(start)/Norm(start))*dScaleTable(end)/ScaleTable(start)
        ! dg=(-EffVer(start)+3.0/16/pi*EffVer(start)**2)*dScaleTable(end)/ScaleTable(start)
        EffVer(end)=EffVer(start)-dg
      enddo
    end subroutine

    double precision function CalcWeight(NewOrder)
      !calculate the weight for ALL diagrams in a given sector
      implicit none
      integer :: NewOrder
      if(NewOrder==0) then 
        CalcWeight=1.0
      else if(NewOrder==1) then
        CalcWeight=Ver4_One(0, 0, Mom0, Mom0, Mom0, Mom0, 1, .true.)
      endif
      ! print *, CalcWeight
      CalcWeight=CalcWeight*ReWeightFactor(NewOrder)
      return
    end function CalcWeight
    
    subroutine IncreaseOrder()
      !increase diagram order by one/change normalization diagram to physical diagram
      implicit none
      double precision :: R, Weight, Kamp, dK, theta, phi, Prop
      if (CurrOrder==Order) return
      PropStep(1)=PropStep(1)+1.0

      ! Generate New Mom
      !!!! Hard way  !!!!!!!!!!!!!!!!!!!!!!!!!!
      dK=2.0
      Kamp=grnd()*dK
      if(Kamp<=0.0) return
      phi=2.0*pi*grnd()
      theta=pi*grnd()
      if(theta==0.0) return
      LoopMom(1, CurrOrder+1)=Kamp*sin(theta)*cos(phi)
      LoopMom(2, CurrOrder+1)=Kamp*sin(theta)*sin(phi)
      LoopMom(3, CurrOrder+1)=Kamp*cos(theta)
      Prop=dK*2.0*pi*pi*sin(theta)*Kamp**(D-1)

      !!! Simple way  !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! do i=1, D
      !   NewMom(i)=kF*(grnd()-0.5)*2
      ! enddo
      ! Prop=Beta*(2.0*kF)**D
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Weight = CalcWeight(CurrOrder+1)

      R=abs(Weight)/abs(CurrWeight)*Prop
      if(grnd()<R) then
        AcceptStep(1)=AcceptStep(1)+1.0
        CurrWeight=Weight
        CurrOrder=CurrOrder+1
      endif
      return
    end subroutine
  
    subroutine DecreaseOrder()
    !decrease diagram order by one/change physical diagram to normalization diagram
      implicit none
      double precision :: R, Weight, Prop
      double precision :: Kamp, dK, SinTheta
      if (CurrOrder==0) return
      !if the current diagrams are already in normalization sector, then return

      !Get proper K proposed probability
      !!!!!!! Hard way !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dK=2.0
      Kamp=norm2(LoopMom(:, CurrOrder))
      if(Kamp<=0.0 .or. Kamp>=dK) return
      SinTheta=norm2(LoopMom(1:2, CurrOrder))/Kamp
      if(SinTheta==0.0) return
      Prop=1.0/(dK*2.0*pi*pi*SinTheta*Kamp**(D-1))

      !!!!!!! Simple way !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! do i=1, D
      !   if(abs(LoopMom(i, LoopNum(CurrOrder)))>kF) return
      ! enddo
      ! Prop=1.0/(Beta*(2.0*kF)**D)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PropStep(2)=PropStep(2)+1.0
      Weight = CalcWeight(CurrOrder-1)
      R=abs(Weight)/abs(CurrWeight)*Prop
      ! print *, R, Weight, CurrWeight, Prop
      if(grnd()<R) then
        AcceptStep(2)=AcceptStep(2)+1.0
        CurrWeight=Weight
        CurrOrder=CurrOrder-1
      endif
      return
    end subroutine

    subroutine ChangeScale()
      implicit none
      !TODO: don't forget to change all LoopMom with scale!
      integer :: OldScale
      double precision :: Weight
      double precision :: R

      ! if(CurrOrder==0) return
      
      OldScale=CurrScale
      if(grnd()<0.5) then
        CurrScale=CurrScale-1
      else
        CurrScale=CurrScale+1
      endif

      ! print *, CurrScale

      if(CurrScale<1 .or. CurrScale>ScaleNum) then
        CurrScale=OldScale
        return
      endif

      PropStep(3) = PropStep(3) + 1.0

      Weight = CalcWeight(CurrOrder)
      R = abs(Weight)/abs(CurrWeight)
      ! print *, R, Weight, CurrWeight
  
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
      double precision :: prop, R
      integer :: Num, i, j
      double precision :: Weight
  
      Num = int( CurrOrder*grnd() ) + 1
  
      PropStep(4) = PropStep(4) + 1.0
      OldMom = LoopMom(:, Num)
  
      call GenerateNewMom(OldMom, NewMom, prop)
  
      if(prop<0.0) return
      LoopMom(:,Num)=NewMom
  
      Weight = CalcWeight(CurrOrder)
      R = prop*abs(Weight)/abs(CurrWeight)
  
      if(grnd()<R) then
        AcceptStep(4) = AcceptStep(4)+1.0
        CurrWeight = Weight
        ! call UpdateState()
      else
        LoopMom(:,Num)=OldMom
      endif
  
      return
    end subroutine

    
    subroutine GenerateNewMom(old, new, prop)
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


    double precision function Ver4_One(LOrder, ROrder, MomL1, MomL2, MomR1, MomR2, InterMomIndex, Simple)
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
        Ver4_One=UWeight*3.0/(2.0*pi)**D
        return
      endif
    end function
    
end program main
