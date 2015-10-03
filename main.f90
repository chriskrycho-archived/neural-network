! Uses stepwise iteration with Euler method to model neural network function.
! Based on the work of the Zochowski group at The University of Michigan.
!   Originated by Chris Krycho, September 23, 2008
!   Current progress as of January 21, 2009

PROGRAM zoch_method

! Modules referenced:
USE connection_module
USE output_module
USE gamma_ksi_time_module
USE values_module
USE integrate_fire_module

! Variable definitions:
implicit none
integer, parameter    :: highprec=kind(0.0d0), loprec=kind(0.0)
integer, allocatable  :: conn(:,:), fire(:), sendcount(:), receivecount(:)
integer               :: i, gmod, n, hcount, neur, ranas, min, max, net, wr, getvals
integer               :: ncount, gridcount, htotal, total, value(8), values(8), totaltime, timecount, step, tts, ttl, tinit
logical               :: Vhigh, Vlow, nrec, rec, doesFire
real, allocatable     :: alpha(:), ksi(:), ksid(:), J(:), Jc(:), Jj(:), jjd(:), Jk(:), jkd(:), gamma(:), gamma_hi(:), gamma_lo(:)
real, allocatable     :: V(:), Vstore(:), Jjtotal(:), jjdtotal(:), Jktotal(:), jkdtotal(:)
real, allocatable     :: ts(:), tstore(:,:), sumpart(:)
real, allocatable     :: W(:), X(:), Y(:), Z(:), add1(:), add2(:)
real                  :: A, B, h, t, tauf, taus, gamma_t, gamma_reset, prop, sumva, time, timeReal, hReal
real                  :: sm_ran, smra, vsm_ran, vsmra, vran1, vran2, check, tas, tas1, tas2, totalvoltage, Vstorecheck


! Initial value statements for ordinary variables:
A = 4.0d0
B = 0.40d0
t = 1.0												! Initialize at non-0 to avoid complications.
timeReal = t - 1.0d0          ! This is the amount of time since the simulation has started running.
tauf = 0.030d0
taus = 0.30d0

! Initial value statements for variables that may be assigned below.
gamma_t = 1.0d0
gamma_reset = 6.0d0


! Get variables from user to customize program:
write (6,*) 'Give the number of neurons per grid (1 - 500):'
read (5,*) ncount
write (6,*) 'Give the number of grids (1-10):'
read (5,*) gridcount
total = ncount*gridcount

! Get the numbers and proportions of neurons sending and receiving connections:
aif: if (gridcount > 1) then
  
  allocate(sendcount(gridcount),receivecount(gridcount),gamma(gridcount),gamma_hi(gridcount),gamma_lo(gridcount))
  write (6,*) 'Give the proportion (0-1) of neurons receiving internetwork input:'
  ado: do getvals=1,gridcount
    
    write (6,*) 'Network', getvals
    read (5,*) prop
    receivecount(getvals) = nint(prop*ncount)   ! receivecount thus holds the actual number of neurons receiving a connection for each network.
  
  end do ado
  write (6,*) 'Give the number of inter-grid connection sources (integer <= number of neurons per grid):'
  bdo: do getvals=1,gridcount
    
    write (6,*) 'Network', getvals
    read (5,*) sendcount(getvals)    ! sendcount thus holds the actual number of neurons sending a connection for each network.
  
  end do bdo
  write (6,*) 'Give the value of gamma_hi and gamma_lo for each network (if constant, supply the same value for hi and lo):'
  cdo: do getvals=1,gridcount
    
    write (6,*) 'Network', getvals
    write (6,*) 'gamma_hi:'
    read (5,*) gamma_hi(getvals)
    write (6,*) 'gamma_lo:'
    read (5,*) gamma_lo(getvals)
  
  end do cdo

else if (gridcount == 1) then
  
  allocate(gamma(gridcount),gamma_hi(gridcount),gamma_lo(gridcount))
  write (6,*) 'Give the value of gamma_hi:'
  read (5,*) gamma_hi(1)
  write (6,*) 'Give the value of gamma_lo:'
  read (5,*) gamma_lo(1)  

else if (gridcount < 1) then
  
  write (6,*) 'You must initialize at least one grid.'
  stop

end if aif


! Get the stepsize and time from the user:
write (6,*) 'Give the stepsize (0.1 - 0.0001):'
read (5,*) h
write (6,*) 'Give the desired number of iterations (2 - 100000):'
read (5,*) htotal
time = h * htotal
totaltime = nint(time)


! Transmission lag times: long and short respectively, specified by division of the time (8ms or 6ms) divided by the stepsize (h):
hReal = real(h)
ttl = int(0.008 / hReal)    ! long - for internetwork
tts = int(0.006 / hReal)    ! short - for intranetwork


! Account for initial t = 1.0d0:
gamma_t = gamma_t + 1.0d0
gamma_reset = gamma_reset + 1.0d0


! Allocate the arrays to hold the number of neurons in the system as specified above:
allocate(V(total), Vstore(total), alpha(total), fire(total), sumpart(ncount), ksi(total), ksid(total))
allocate(J(total), Jc(total), Jj(total), jjd(total), Jk(total), jkd(total))
allocate(Jjtotal(total), jjdtotal(total), Jktotal(total), jkdtotal(total))
allocate(ts(total), tstore(ttl,total))
allocate(W(total), X(total), Y(total), Z(total), conn(total,total))


! Reinitialize the random number generator with a new seed each time the program runs.
call date_and_time(values = values)
value = values(7) + values(8)
call random_seed(put = value)


! Initial value assignments for variables on ranges: set all to 0 to clear any previous data, then assign as below.
Vstore = 0.0d0

Jj = 0.0d0
Jjtotal = 0.0d0
jjd = 0.0d0
jjdtotal = 0.0d0
Jk = 0.0d0
Jktotal = 0.0d0
jkd = 0.0d0
jkdtotal = 0.0d0
Jc = 0.0d0

fire = 0

ts = 0.0d0
tstore = 0.0d0

alpha = 0.0d0
ksi = 0.0d0

W = 0.0d0
X = 0.0d0
Y = 0.0d0
Z = 0.0d0

conn = 0
sumpart = 0.0d0


! Use uniform random number generator for values on alpha:
do ranas=1,total
  vsmra = vsmranres(vsm_ran)
  alpha(ranas) = 1.250d0 + vsmra
end do


! Assign the initial value of the white noise factor ksi with the same random generator:
do ranas=1,total
  smra = smranres(sm_ran)
  ksi(ranas) = 1.20d0 + smra
end do


! Call random number for values on V.
do ranas=1,total
  call random_number(vran1)
  call random_number(vran2)
  call random_number(check)
  initcheck: if (check <= .05) then   ! barrier chosen to represent number of neurons in post-firing state at any "ordinary" time.
    V(ranas) = vran1 + vran2
    if ((vran1+vran2) >= 1.0d0) fire(ranas) = 1
  else initcheck
    V(ranas) = 0
  end if initcheck
end do


! Assign the time since a given element has fired, accounting for recent/non-recent firing.
do ranas=1,total
  if (V(ranas) .eq. 0) then
    call random_number(tas)
    tstore(1,ranas) = (tas*.1)
  else
    call random_number(tas1)
    call random_number(tas2)
    tstore(1,ranas) = (tas1*5) + (tas2*5)
  end if
  do tinit=2,ttl
    tstore(tinit,ranas) = tstore(1,ranas) + h * tinit
  end do
end do


! Assign proper initial values of internetwork current.
call vals(Jj,jjd,Jk,jkd,tstore,tts,ttl,taus,tauf,ncount,gridcount,total)
call sumvals(J,Jj,jjd,Jk,jkd,conn,ncount,gridcount,total,gamma,Jjtotal,jjdtotal,Jktotal,jkdtotal)


! Clear all the old file information:
call clear()


! Now set the initial state of the connections and perform a rewiring (see the corresponding module for description), writing the files for analysis.
! Initialize the network (that is, array) states.
call init(conn,ncount,gridcount,total)

! Handle intranetwork wiring.
call intra(conn,ncount,gridcount,total)

! If there is more than a single grid, add connections between the grids with subroutine inter:
if (gridcount > 1) call inter(sendcount,receivecount,gridcount,ncount,conn,total)


! Iterate with the Euler method, writing the results to a file.
!   hcount is the number of times we iterate the approximation, with stepsize h.
iteration: do hcount=1,htotal

  ! Change the value of white noise in each network for each iteration.
  call ksi_change(ksi,ksid,total,h)
  
  ! Check how gamma should now behave.
  call gamma_change(t,gamma,gamma_t,gamma_reset,gamma_lo,gamma_hi,gridcount)

  ! Now call the values, current values, and sum on values for the system with the changed state.
  call vals(Jj,jjd,Jk,jkd,tstore,tts,ttl,taus,tauf,ncount,gridcount,total)
  call sumvals(J,Jj,jjd,Jk,jkd,conn,ncount,gridcount,total,gamma,Jjtotal,jjdtotal,Jktotal,jkdtotal)

  ! Integrate incoming current across each neuron.
  call integrate(V,Vstore,tstore,tts,ttl,ncount,gridcount,total,alpha,ksi,ksid,h,A,B,Jjtotal,jjdtotal,Jktotal,jkdtotal)
    
  ! Check if a given neuron fires, and if so reset the time. Also account for boundary condition on V (must be >= 0):
  call firecheck(V,Vstore,fire,doesFire,total,tstore,ttl,ncount,h)
  
  call curvals(Jc,tstore,ttl,taus,tauf,ncount,gridcount,total)
  
  !J = Jjtotal + Jktotal
  
  ! Write the results to a file for analysis.
  call writer(ncount, gridcount, total, timeReal, V, Vstore, Jc, fire, tstore, ttl)
  
  ! Reset fire array:
  fire = 0

  ! For various subparts of the total time, print the time elapsed to the screen for user feedback.
  if ((hcount == (htotal / 10)).OR.(hcount == (htotal / 5)).OR.(hcount == (htotal / 2)).OR.(hcount == htotal)) then
    
    write (6,*) timeReal
  
  end if
  
  ! Clear all variable values prior to starting the next iteration:
  Jj = 0
  Jjtotal = 0
  jjd = 0
  jjdtotal = 0
  Jk = 0
  Jktotal = 0
  jkd = 0
  jkdtotal = 0
  Jc = 0
  J = 0
  V = 0

  ! Increase the overall time and display time by h as well:
  t = t + h
  timeReal = t - 1.0d0

end do iteration


! ---> Subroutines and Functions --->

contains

! ---> Functions --->
! Uniform random number generator:
FUNCTION  smranres( sm_ran )
! Originally created by Alan Miller. (http://users.bigpond.net.au/amiller/)
!   Modified by Chris Krycho and Dr. Kieran Mullen (University of Oklahoma).
!   Generate a random normal deviate using the polar method.
!   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!              normal variables', Siam Rev., vol.6, 260-264, 1964.

IMPLICIT none
REAL  :: sm_ran, smranres

! Local variables
REAL            :: cc,sum
REAL, SAVE      :: cd,sln
LOGICAL, SAVE   :: second = .FALSE.
REAL, PARAMETER :: one = 1.0, vsmall = TINY( one )

IF (second) THEN
! If second, use the second random number generated on last call
  second = .false.
  sm_ran = cd*sln

ELSE
! First call; generate a pair of random normals
  second = .true.
  DO
    CALL RANDOM_NUMBER( cc )
    CALL RANDOM_NUMBER( cd )
    cc = SCALE( cc, 1 ) - one
    cd = SCALE( cd, 1 ) - one
    sum = cc*cc + cd*cd + vsmall         ! vsmall added to prevent LOG(zero) / zero
    IF(sum < one) EXIT
  END DO
  sln = (SQRT(- SCALE( LOG(sum), 1 ) / sum)) / 6    ! div by 6 to ensure range
  sm_ran = cc*sln
END IF

  smranres=sm_ran

RETURN
END FUNCTION smranres


! ---> As above: an even smaller scale --->
FUNCTION  vsmranres( vsm_ran )

IMPLICIT none
REAL  :: vsm_ran, vsmranres

! Local variables
REAL            :: c,sum
REAL, SAVE      :: d,sln
LOGICAL, SAVE   :: second = .FALSE.
REAL, PARAMETER :: one = 1.0, vsmall = TINY( one )

IF (second) THEN
! If second, use the second random number generated on last call
  second = .false.
  vsm_ran = d*sln

ELSE
! First call; generate a pair of random normals
  second = .true.
  DO
    CALL RANDOM_NUMBER( c )
    CALL RANDOM_NUMBER( d )
    c = SCALE( c, 1 ) - one
    d = SCALE( d, 1 ) - one
    sum = c*c + d*d + vsmall         ! vsmall added to prevent LOG(zero) / zero
    IF(sum < one) EXIT
  END DO
  sln = (SQRT(- SCALE( LOG(sum), 1 ) / sum)) / 10    ! div by 10 to ensure range
  vsm_ran = c*sln
END IF

  vsmranres = vsm_ran

RETURN
END FUNCTION vsmranres



! ---> End subroutines and functions --->

END PROGRAM zoch_method

! ---> End <--- !