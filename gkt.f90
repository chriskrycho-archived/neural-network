module gamma_ksi_time_module

contains
SUBROUTINE gamma_change(t, gamma, gamma_t, gamma_reset, gamma_lo, gamma_hi, gridcount)

! Variables passed by reference:
real    :: t, gamma(gridcount), gamma_t, gamma_reset, gamma_lo(gridcount), gamma_hi(gridcount)
integer :: gridcount

gammatime: if ((t >= gamma_t).AND.(t <= gamma_reset)) then     ! Used if gamma, rather than E, varies
    
  gamma = gamma_lo
  
else
    
  gamma = gamma_hi
  
end if gammatime
  
END SUBROUTINE gamma_change


! --------- !


SUBROUTINE ksi_change(ksi,ksid,total,h)
! Changes the value of the white noise constant for each neuron in every iteration.

! Local variables:
implicit none
integer   :: count
logical   :: random
real      :: ksi_ran, ksi_randomizer, ksi_start

! Variables passed by argument:
integer   :: total
real      :: ksi(total), ksid(total), h

! Set the initial value for random so the function behaves properly.
random = .true.

! Change the value of the white noise on each neuron every iteration, causing only small changes in ksi.  
ksi_count: do count=1,total
  
  ! For each neuron, store the previous value of ksi, as it will be necessary later to compute the derivative.
  ksi_start = ksi(count)
  random_call: do while (random .eqv. .true.)
    
    call random_number(ksi_ran)
    random_if: if (ksi_ran > 0.1) then  ! in build c was 0.05
      
      random = .true.
    
    else
      
      random = .false.
    
    end if random_if
  
  end do random_call
  
  ! Choose which direction to alter the particular value of ksi this iteration. In general, if above 0.5 add, below 0.5 subtract.
  call random_number(ksi_randomizer)
  which_way: if (ksi_randomizer > 0.5) then
    
    ! To account for boundary conditions, add as long as the addition does not cause ksi to meet or exceed 1.4. If it does, subtract instead.
    high: if (ksi(count) + ksi_ran < 1.4)
    
      ksi(count) = ksi(count) + ksi_ran
    
    else
      
      ksi(count) = ksi(count) - ksi_ran
      
    end if high
  
  else if (ksi_randomizer  <= 0.5) then
    
    ! Likewise, only subtract if the subtraction does not cause ksi to fall to or below 1.0; add if it would.
    if (ksi(count) - ksi_ran > 1.0)
    
      ksi(count) = ksi(count) - ksi_ran
      
    else
    
      ksi(count) = ksi(count) + ksi_ran
  
  end if which_way
  
  ! The derivative of ksi is simple d(ksi)/dt, but dt is simply h, and d(ksi) is the difference in the new and previous values of ksi.
  ksid(count) = (ksi(count) - ksi_start) / h

end do ksi_count
  
END SUBROUTINE ksi_change


! --------- !


SUBROUTINE timestep(tstore,ttl,ncount,n,h,doesFire)
! Steps the timestore in the proper order.

! local variables:
integer   :: step

! variables passed by reference:
integer   :: n, ttl, ncount
real      :: tstore(ttl,ncount)
logical   :: doesFire

! Set the values of lagged time (for the lagged current values generated each cycle).
! Iterate ttl (8 msec / h) times so that fully 8 msec (the lag time for separated networks) are accounted for:
do step=ttl,2,-1
      
  tstore(step,n) = tstore(step - 1,n)
      
end do

! If the neuron fires (i.e. doesFire = true) then set tstore(1,n) = 0; otherwise increase tstore(1,n) by h.
fire: if (doesFire) then
  
  tstore(1,n) = 0

else
  
  tstore(1,n) = tstore(1,n) + h

end if fire
      
END SUBROUTINE timestep


end module gamma_ksi_time_module