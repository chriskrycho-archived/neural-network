module connection_module

contains
SUBROUTINE init(conn,ncount,gridcount,total)
! Initialize the connection matrix.

! Variable declarations:
implicit none
integer               :: conn(total,total)
integer               :: ncount, length, gridcount, total
integer               :: send, receive, g, min, max, wr, i
logical               :: a, b, c, d, e, f, k, l, m, n, o, p

! Initialize the connection matrix with a set of standard grid values.
!   The grid is a 5x5 grid, with nodes numbered as they would be read: the first length is 1, 2, 3, 4, 5;
!   the last length is 21, 22, 23, 24, 25.
! The "i" iteration represents rows, the "j" iteration represents columns. (This is more important when,
!   later, connections become directional.) Note: we must be careful here at multiples of the number of
!   nodes, since they do not connect in the same pattern. This boundary condition is accounted for below.

! Below is the initial assignment code, written as a subroutine to be called by the main program.

length = int(sqrt(real(ncount))) ! number of neurons per length in grid

conn = 0    ! all connections are initially set to 0 so no memory bits float into the operation.

! Now assign all possible connections, ignoring boundary condition issues.
! The outer loop accounts for the various networks:
gdo: do g=1,gridcount

  ! Now we use the grid information to specify the range assigned, which allows us to move one network at a time.
  min = (g - 1) * ncount + 1     ! min = 0n+1, n+1, 2n+1, ... , (gridcount-1)*n + 1
  max = ncount * g             ! max = n, 2n, ... , gridcount*n
  
  senddo: do send=min,max
    
    receivedo: do receive=min,max
      
      a = (send .EQ. (receive - 1))
      b = (send .EQ. (receive + 1))
      c = (send .EQ. (receive - 2))
      d = (send .EQ. (receive + 2))
      e = (send .EQ. (receive - length))
      f = (send .EQ. (receive - length - 1))
      k = (send .EQ. (receive - length + 1))
      l = (send .EQ. (receive + length))
      m = (send .EQ. (receive + length - 1))
      n = (send .EQ. (receive + length + 1))
      o = (send .EQ. (receive + 2 * length))
      p = (send .EQ. (receive - 2 * length))
      ifnorm: if (a.OR.b.OR.c.OR.d.OR.e.OR.f.OR.k.OR.l.OR.m.OR.n.OR.o.OR.p) then
        
        conn(send,receive) = 1
      
      else
        
        conn(send,receive) = 0
      
      end if ifnorm
    
    end do receivedo
  
  end do senddo

end do gdo

! Account now for the "boundaries" - that is, the jump from the oth to the (length+1)th position.
!   NOTE: Will not work properly for a 2x2 (n = 4) grid because of logical construction. All
!   larger gridcount will function properly.

! Start at the first corner, and jump upwards to the end of the grid (at total) by length of rows:
corners: do i = length, total, length
  
  ! Eliminating corner problems from range 1:
  if (i + 1 <= total) then
  
    conn(i , i + 1) = 0
    conn(i + 1 , i) = 0
  
  end if
  
  ! Eliminating diagonals of range 1 (i.e. in a 3x3 grid, neuron 3 to neuron 7)
  if (i + length + 1 <= total) then
    
    conn(i , i + length + 1) = 0
    conn(i + length + 1 , i) = 0
  
  end if
  
  ! Eliminating corner problems at range 2 (i.e in a 4x4 grid, neuron 4 to neuron 6: more than range of 2 connections between)
  if (i + 2 <= total) then
    
    conn(i , i + 2) = 0
    conn(i + 2 , i) = 0
  
  end if
  
  ! Eliminating edge gaps at range 2 (i.e. in a 4x4 grid, neuron 3 to neuron 5: more than range of 2 connections between)
  outer: if (total > 9) then
    
    inner: if (i + 1 <= total) then
    
      conn(i - 1 , i + 1) = 0
      conn(i + 1 , i - 1) = 0
    
    end if inner
  
  end if outer
  
end do corners


! Eliminating the other set of diagonals on corners (for gridcount larger than 3x3):
if (total > 9) then
  
  diagonals: do i = 1, total, length
  
    conn(i , i + length - 1) = 0
    conn(i + length - 1 , i) = 0
  
  end do diagonals

end if

! Then write the files for analysis.
open(3,file="./data/conninit.txt",status='old')
  
  do wr=1,total
    
    write(3,*) conn(wr,1:total)
  
  end do
  
close(3,status='keep')

END SUBROUTINE init


! --------- !


SUBROUTINE intra(conn,ncount,gridcount,total)
! Rewire inside each network.

! Local variable declarations:
implicit none
integer               :: o, ran2i, values(8), value(8), raneck, raneck1
integer               :: k, k1, g, g2, min, max, ch, eck, wr, ii, q, r, xint, yint, blah
logical               :: gbound1, gbound2
real                  :: nre, ran1, ran2, ran3, xre, yre
! Variables passed by argument:
integer               :: conn(ncount * gridcount , ncount * gridcount), ncount, gridcount, total

! Reinitialize the random number generator with a new seed each time the program runs.
call date_and_time(values=values)
value = values(7) + values(8)
call random_seed(put=value)

nre = real(ncount)

! Rewire the system by checking whether connections exist, and if they do, possibly rewiring them.
!   As in the case of the initial wiring, do this for each grid.
outer: if (gridcount > 1) then
  
  g2do: do g2 = 1 , gridcount
    
    ! min = 0n+1, n+1, 2n+1, ... , (gridcount-1) * ncount + 1
    min = (g2 - 1) * ncount + 1
    
    ! max = ncount, 2 * ncount, ... , gridcount * ncount
    max = ncount * g2
    kdo: do k = min , max
      
      check: do ch = min , max
        
        ! if there is a connection, then possibly rewire.
        conncheck: if (conn(k , ch) == 1) then
          
          ! call a random number between 0 and 1:
          call random_number(ran1)
          
          ! probability barrier chosen to be consistent with Zoch. group: only 0.3 of the existing connections are rewired
          prob: if (ran1 >= 0.7) then
            
            ! get a random integer between the bounds of the current grid, adding the integer to min - 1 (so that 1 is a possible choice): 
            call random_number(ran2)
            ran2i = int(ran2 * nre) + (min - 1)
            
            ! Ensures we avoid the self-connection issue initially: instead of adding 0 to the position examined, add either 1 or n.
            ranbound: if (ran2i == 0) then
              
              call random_number(ran3)
              
              ! Pushes 0-values for ran2i either to the bottom or the top of the grid.
              rbin: if (ran3 < 0.5) then
                
                ran2i = min
              
              else rbin
                
                ran2i = max
              
              end if rbin
            
            end if ranbound
            
            ! Check conn(k , ran2i); if that's connected then cycle upward till it is *not* connected.
            checkloop: do eck = 0 , ncount - 1
              
              ! If eck + ran2i > ncount, loop around so as to prevent wiring to other sectors! This is the test variable.
              raneck = ran2i + eck
              
              ! Checks if raneck is in the proper grid set. Cycle up to fix if not.
              eckbound: do ii = 1 , gridcount
                
                jbih: if ((raneck >= min) .AND. (raneck <= max)) then
                  
                  exit eckbound
                
                elseif (raneck < min) then
                  
                  raneck = raneck + ncount
                  cycle eckbound
                
                elseif (raneck > max) then
                  
                  raneck = raneck - ncount
                  cycle eckbound
                
                endif jbih
              
              end do eckbound
              
              ! Then, having ensured that raneck is in the proper grid and accounted for boundaries:
              ! Is k already sending a connection to raneck?
              connif: if (conn(k , raneck) == 0) then
                
                ! Is k the same as raneck?
                self: if (k /= raneck) then
                  
                  ! If no to both, then send from k to raneck and exit the checkloop, cycling to the next connection k.
                  conn(k , raneck) = 1
                  exit checkloop
                
                ! If k is the same as raneck, cycle the checkloop to look for another value of raneck with the same value of k.
                elseif (k == raneck) then
                  
                  cycle checkloop
                
                end if self
              
              ! If k and raneck are already connected, cycle the checkloop to look for another value of raneck with the same value of k.
              elseif (conn(k,raneck) == 1) then
                
                cycle checkloop
              
              end if connif
            
            end do checkloop
            
            ! Having satisfied these conditions and generated a new connection, eliminate the old connection.
            conn(k , ch) = 0
          
          end if prob
        
        ! if no connection, cycle the loop to check the next possible connection
        else conncheck
          
          cycle check
        
        end if conncheck
      
      end do check
    
    end do kdo
  
  end do g2do

! The logic for the case that there is only one grid is the same, but the conditions relating to multiple gridcount have been left out.
elseif (gridcount == 1) then
  
  kdoa: do k = 1 , ncount
    
    checka: do ch = 1 , ncount
      
      connchecka: if (conn(k , ch) == 1) then
        
        call random_number(ran1)
        
        ! probability barrier chosen to be consistent with Zoch. group: only 0.3 of the existing connections are rewired
        proba: if (ran1 > 0.7) then
          
          ! get a random integer between the bounds of the current grid, adding the integer to min - 1 (so that 1 is a possible choice): 
          call random_number(ran2)
          ran2i = int(ran2 * nre)
          
          ! Ensures we avoid the self-connection issue initially: instead of adding 0 to the position examined, add either 1 or n.
          ranbounda: if (ran2i == 0) then
            
            call random_number(ran3)
            
            ! Pushes 0-values for ran2i either to the bottom or the top of the grid.
            rbina: if (ran3 > 0.5) then
              
              ran2i = 1
            
            else rbina
              
              ran2i = 0
            
            end if rbina
          
          end if ranbounda
          
          ! Check current ran2i, if that's connected then cycle upward till it is *not* connected.
          checkloopa: do eck = 0, ncount - 1
            
            ! If eck + ran2i > ncount, loop around so as to prevent wiring to other sectors! This is the test variable.
            raneck = ran2i + eck
            jbas: if (raneck > ncount) then
              
              raneck = raneck - ncount
            
            elseif (raneck <= 0) then
              
              raneck = raneck + ncount
            
            endif jbas
            
            ! Then, having ensured that raneck is in the proper range and accounted for boundaries:
            ! Is k already sending a connection to raneck?
            connifa: if (conn(k , raneck) == 0) then
              
              ! Is k the same as raneck?
              selfa: if (k /= raneck) then
                
                ! If no to both, then send from k to raneck and exit the checkloop, cycling to the next connection k.
                conn(k , raneck) = 1
                exit checkloopa
              
              ! If k is the same as raneck, cycle the checkloop to look for another value of raneck with the same value of k.
              elseif (k == raneck) then
                
                cycle checkloopa
              
              end if selfa
            
            ! If k and raneck are already connected, cycle the checkloop to look for another value of raneck with the same value of k.
            elseif (conn(k , raneck) == 1) then
              
              cycle checkloopa
            
            else connifa
              
              write(6,*) 'Error in conn array. Please check code.', k, ch
              stop
            
            end if connifa
          
          end do checkloopa
          
          ! Having satisfied these conditions and generated a new connection, eliminate the old connection.
          conn(k , ch) = 0
        
        end if proba
      
      else connchecka
        
        cycle checka
      
      end if connchecka
    
    end do checka
  
  end do kdoa

! for the case that gridcount <= 0.
else
  write (6,*) 'Error, please enter an integer value greater than 0 for the number of networks.'
  stop
end if outer

! Write the file after intranetwork wiring.
open(3,file="./data/connintra.txt",status='old')

  do wr=1, ncount * gridcount
  
    write(3,*) conn(wr,1:total)
  
  end do
  
close(3,status='keep')

END SUBROUTINE intra


! --------- !


SUBROUTINE inter(send,rec,gridcount,ncount,conn,total)
! Generate connections between networks. All of a user-specified proportion of the neurons in each network receive
!   connections from a user-specified number of random neurons in the other network. (Note: each neuron in the proportion
!   receives that number of inputs from the other network).

! Local variables:
implicit none
real                    :: r1, r2
integer, allocatable    :: sendarr(:,:), recarr(:,:)
integer                 :: g, h, i, j, k, l, m, maxsend, maxrec, n, ncount, x, sendstore, y, receivestore, blah, sendmax, sendmin
integer                 :: value(8), values(8), sumrec, sumsend, a, b, dup, receivemax, receivemin, wr
! Variables passed by argument:
integer                 :: gridcount, send(gridcount), rec(gridcount), conn(ncount * gridcount , ncount * gridcount), total

! Reinitialize the random number generator with a new seed each time the program runs.
call date_and_time(values=values)
value = values(7) + values(8)
call random_seed(put=value)

! Define and allocate:
maxsend = maxval(send)   ! Use max to ensure that the arrays are large enough.
maxrec = maxval(rec)

! sendarr specifies the neurons sending a connection from each grid.
! recarr specifies the neurons receiving a connection in each grid.
! Format: array(grid,number), so neuron 15 in grid 1 would be (1,15)
allocate(sendarr(gridcount,maxsend))
allocate(recarr(gridcount,maxrec))

! Set initial values of sending and receiving arrays to 0:
sendarr = 0
recarr = 0

! Choose the neurons that are the receivers of intergrid connections:
ido: do i=1,gridcount
  
  j = 1
  
  ! if 0, the grid has no incoming connections from other gridcount
  if (rec(i) > 0) then
    
    ! Where j > rec(i), the value remains 0. This keeps there from being extra connections made.
    jdo: do while (j <= rec(i))
      
      call random_number(r2)
      r2 = r2*(real(ncount))
      
      if (nint(r2) == 0) r2 = r2 + 1.0d0
      
      recarr(i,j) = nint(r2) + (i-1)*ncount
      
      ! Check for duplicates: if the value stored at (grid,dup) is the same as (grid,j), cycle j.
      !   That way each neuron only receives connections once.
      ! Check: for duplicates
      dup2do: do dup=1,ncount
        
        if (dup == j) cycle dup2do    ! Prevents the system from looping once it hits the selfsame term.
        
        ! if a duplicate is found, call a new random number by cycling jdo without increasing the count on j.
        cif: if ( recarr(i,dup) == recarr(i,j) ) then
          
          cycle jdo
        
        ! if no duplicate is found, check through the rest of the possibles, exiting normally at the end.
        else
          
          cycle dup2do
        
        end if cif
      
      enddo dup2do
      j = j + 1
    
    end do jdo
  
  else
    
    recarr(i,:) = 0
    exit ido
  
  end if

end do ido

sendstore = 1
receivestore = 1

kdo: do k=1,gridcount
  
  ldo: do l=1,gridcount
    
    if (k == l) cycle ldo
    
    ! Specify the maximum neuron for the sending grid in a given network.
    sendmax = ncount*k
    
    ! Do likewise for the minimum neuron.
    sendmin = ncount*(k-1) + 1
    
    ! Do the same but for the maximum neuron number in the receiving grid.
    receivemax = ncount*l
    
    ! Do the same for the minimum neuron in receiving grid.
    receivemin = ncount*(l-1) + 1
    
    ! Sets the values for the initial constant for the m and n cycles. Start at 1 initially.
    sendstore = 1
    receivestore = 1
    ndo: do n=1,maxrec

      ! Choose the neurons from each grid that are the source of intergrid connections
      !   (repeated for each neuron receiving a connection):
      gdo: do g=1,gridcount
        
        ! if 0, the grid sends no connections to other gridcount
        if (send(g) > 0) then
          
          h = 1
          
          ! Where h > send(g), the value remains 0.
          hdo: do while (h <= send(g))
            
            call random_number(r1)
            r1 = r1*(real(ncount))
            
            if (nint(r1) == 0) r1 = r1 + 1.0d0
            
            sendarr(g,h) = nint(r1) + (g-1)*ncount
            
            ! Check for duplicates: if the value stored at (grid,dup) is the same as (grid,j), cycle j.
            !   That way each neuron only sends connections once.
            dup1do: do dup=1,ncount
              
              ! Prevents the system from looping once it hits the selfsame term.
              if (dup == h) cycle dup1do
              
              aif: if ( sendarr(g,dup) == sendarr(g,h) ) then
                
                ! if a duplicate is found, call a new random number by cycling hdo without increasing the count on h.
                cycle hdo
              
              else
                
                ! if no duplicate is found, check through the rest of the possibles, exiting normally at the end.
                cycle dup1do
              
              end if aif
            
            enddo dup1do
            h = h + 1
          
          end do hdo
        
        else
          
          sendarr(g,:) = 0
          exit gdo
        
        end if
      
      end do gdo

      ! With this information, generate connections between networks. Each neuron in recarr should receive connections
      !   from the current random set of neurons in sendarr.
      mdo: do m=1,maxsend
        x = sendarr(k,m)
        y = recarr(l,n)
        
        if (x > sendmax) x = x - ncount
        
        if (y > receivemax) y = y - ncount
        
        ! If there is already a connection present, move to the next point to receive a connection.
        connectioncheck: if (conn(x,y) == 1) then
          
          cycle mdo
        
        ! If it is not already connected, connect.
        elseif (conn(x,y) == 0) then
            
            ! only make an input, not a bidirectional connection:
            conn(x,y) = 1
            
            ! Then continue checking the rest of the possible connections.
            cycle mdo
        
        end if connectioncheck
      
      end do mdo
    
    end do ndo
  
  end do ldo

end do kdo

! Print the results for analysis.
open(3,file="./data/connfinal.txt",status='old')
  
  do wr=1, ncount * gridcount
    
    write(3,*) conn(wr,1:total)
  
  end do
  
close(3,status='keep')

END SUBROUTINE inter


end module connection_module