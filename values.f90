module values_module

contains
SUBROUTINE curvals(jc,tstore,ttl,taus,tauf,n,grids,total)
! Calculates the current leaving a given element at a given time.

! Local variables:
implicit none
integer   :: gr, min, max, m
! Variables passed by argument:
real      :: Jc(total), tstore(ttl,total), taus, tauf
integer   :: n, grids, total, ttl

! Assign the current values in the network at the exact moment specified.

voltages: do gr=1,grids

  min = (gr - 1) * n + 1     ! min = 0n+1, n+1, 2n+1, ... , (grids-1)*n + 1
  max = n * gr             ! max = n, 2n, ... , grids*n
  mdo: do m=min,max

    Jc(m) = exp(-tstore(1 , m) / taus) - exp(-tstore(1 , m) / tauf)

  end do mdo

end do voltages

END SUBROUTINE curvals


! --------- !


SUBROUTINE vals(Jj,jjd,Jk,jkd,tstore,tts,ttl,taus,tauf,n,grids,total)
! Calculates values on a given element at a given time, accounting for lag between neurons.


! Local variables:
implicit none
integer   :: gr, min, max, m

! Variables passed by argument:
real      :: Jj(total), jjd(total), Jk(total), jkd(total), tstore(ttl , total), taus, tauf
integer   :: n, grids, total, tts, ttl

! Assign the values that each incoming current should have (and the way that each is changing, which is important
!   for an Euler approximation); J and K are identical in form but represent currents coming from different sectors.
!   Jj, Jk --> actual currents; jjd, jkd --> rate of change (derivative) of those currents.
! For both J and Jk, we chose the corresponding tstore parameter (tstore(6 msec,n) and tstore(8 msec, n) respectively)
!   rather than tstore(present, n) to account for the delay between individual neurons and between networks.

currents: do gr=1,grids

  min = (gr - 1) * n + 1     ! min = 0n+1, n+1, 2n+1, ... , (grids-1)*n + 1
  max = n * gr             ! max = n, 2n, ... , grids*n

  mdo: do m = 1, total

    if ((m >= min) .AND. (m <= max)) then

      Jj(m) = exp(-tstore(tts , m) / taus) - exp(-tstore(tts , m) / tauf)
      jjd(m) = (-1 / taus) * exp(-tstore(tts , m) / taus) + (1 / tauf) * exp(-tstore(tts , m) / tauf)

    else

      Jk(m) = exp(-tstore(ttl , m) / taus) - exp(-tstore(ttl , m) / tauf)
      jkd(m) = (-1 / taus) * exp(-tstore(ttl , m) / taus) + (1 / tauf) * exp(-tstore(ttl , m) / tauf)

    end if

  end do mdo

end do currents

end subroutine vals


! --------- !


subroutine sumvals(J,Jj,jjd,Jk,jkd,conn,n,grids,total,gamma, Jjtotal, jjdtotal, Jktotal, jkdtotal)
! Sums the values of current in preparation for integration.

! Overall method: sum the incoming current to each neuron, dealing with both internetwork and intranetwork current separately.
!   After finishing the summation, add the present current on each neuron to this total current. Then compare the total current
!   on the network to a predefined gamma_cut value, here specified as gamma. If the value achieved by summation is less than
!   that value, reset the incoming current to 0. Then pass the total values back to the main program.

! Local variables:
implicit none
integer   :: g, network, min, max, m, o, p, send, sendj, sendk, neuron, receive, receivej, receivek, ra

! Variables passed by argument:
integer   :: conn(total,total), n, grids, total
real      :: J(total), Jj(total), jjd(total), Jk(total), jkd(total), gamma(grids)
real      :: Jjtotal(total), jjdtotal(total), Jktotal(total), jkdtotal(total) ! used to sum the variables

! Initialize the values on the sum arrays as 0:
Jjtotal = 0
jjdtotal = 0
Jktotal = 0
jkdtotal = 0

! Now sum all incoming voltages, accounting for the different lags intra- and inter-network:
! If there are multiple interacting grids, do as follows:
outer: if (grids > 1) then

  ! Specify which network we are summing neurons from. Min and max will change as grid changes.
  gdo: do g = 1, grids

    ! min = 0n + 1, n + 1, 2n + 1, ... , ((grids - 1) * n) + 1
    min = (g - 1) * n + 1

    ! max = n, 2n, ... , grids * n
    max = n * g

    ! Check connections in the same network for Jj - examine both sending and receiving connections in the same network:
    receivejdo: do receivej = min, max

      sendjdo: do sendj = min, max

        ! since conn is structured conn(send,receive) and we are checking for current incoming on sendj:
        jconncheck: if (conn(sendj,receivej) == 1) then

          ! adds present current on sendj (the receiving neuron) with lagged current from receivej:
          Jjtotal(receivej) = Jjtotal(receivej) + Jj(sendj)
          jjdtotal(receivej) = jjdtotal(receivej) + jjd(sendj)

        end if jconncheck

      end do sendjdo

    end do receivejdo

    ! examine possible RECEIVING neurons in this network only.
    receivekdo: do receivek = min, max

      ! examine all possible SENDING neurons; eliminate intranetwork connections in next line.
      sendkdo: do sendk = 1, total

        rbounds: if ((sendk < min) .OR. (sendk > max)) then

          ! as in jconncheck
          kconncheck: if (conn(sendk,receivek) == 1) then

            ! adds present current on receivek with lagged current from sendk.
            Jktotal(receivek) = Jktotal(receivek) + Jk(sendk)
            jkdtotal(receivek) = jkdtotal(receivek) + jkd(sendk)

          end if kconncheck

        end if rbounds

      end do sendkdo

    end do receivekdo

  end do gdo

  ! The total current on a neuron is the sum of all incoming currents.
  J = Jjtotal + Jktotal

  ! Specifies which network (and therefore which threshold gamma_cut) we are examining.
  networkdo: do network = 1, grids

    ! min = 0n+1, n+1, 2n+1, ... , ((grids - 1) * n) + 1
    min = (network - 1) * n + 1

    ! max = n, 2n, ... , grids * n
    max = n * network
    neurondo: do neuron=min,max

      ! if the total of all incoming current is below this threshold level, it resets to 0.
      gammacheck: if (J(neuron) < gamma(network)) then

        ! Reset the value for the next outer iteration.
        J(neuron) = 0

        ! Set all those values which will be integrated to 0.
        Jjtotal(neuron) = 0
        Jktotal(neuron) = 0
        !jjdtotal(neuron) = 0
        !jkdtotal(neuron) = 0

      end if gammacheck

    end do neurondo

  end do networkdo

! If there is a single grid, instead proceed thus:
elseif (grids == 1) then

  receivedo: do receive=1,n

    senddo: do send=1,n

      ! as in jconncheck, kconncheck
      conncheck: if (conn(send,receive) == 1) then

        Jjtotal(receive) = Jjtotal(receive) + Jj(send)
        jjdtotal(receive) = jjdtotal(receive) + jjd(send)

      end if conncheck

    end do senddo

  end do receivedo

  gammacycle: do neuron=1,n

    gammacheck2: if (Jjtotal(neuron) <= gamma(1)) then

      ! Reset the sum for the next iteration:
      Jjtotal(neuron) = 0

      ! Set the value of Jj to 0 since the barrier was not passed.
      Jj(neuron) = 0

    end if gammacheck2

  end do gammacycle

! If the number of grids is less than 1, the whole thing is wrong.
else

  write(6,*) 'Error, wrong number of grids supplied.'
  stop

end if outer

END SUBROUTINE sumvals


end module values_module
