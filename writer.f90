module output_module

contains
SUBROUTINE clear
! The subroutine that clears old versions of the files:

! Eliminate old versions of files. (When working properly, integrate a query here.)

open(3,file="./data/fire.txt",status='unknown')
close(3,status='delete')
open(3,file="./data/fire.txt",status='new')
close(3,status='keep')

open(3,file="./data/fireIndi.txt",status='unknown')
close(3,status='delete')
open(3,file="./data/fireIndi.txt",status='new')
close(3,status='keep')

open(3,file="./data/current.txt",status='unknown')
close(3,status='delete')
open(3,file="./data/current.txt",status='new')
close(3,status='keep')

open(3,file="./data/indicurrent.txt",status='unknown')
close(3,status='delete')
open(3,file="./data/indicurrent.txt",status='new')
close(3,status='keep')

open(3,file="./data/voltage.txt",status='unknown')
close(3,status='delete')
open(3,file="./data/voltage.txt",status='new')
close(3,status='keep')

open(3,file="./data/neuronVoltages.txt",status='unknown')
close(3,status='delete')
open(3,file="./data/neuronVoltages.txt",status='new')
close(3,status='keep')

open(3,file="./data/neuronVStore.txt",status='unknown')
close(3,status='delete')
open(3,file="./data/neuronVStore.txt",status='new')
close(3,status='keep')

open(3,file="./data/conninit.txt",status='unknown')
close(3,status='delete')
open(3,file="./data/conninit.txt",status='new')
close(3,status='keep')

open(3,file="./data/connintra.txt",status='unknown')
close(3,status='delete')
open(3,file="./data/connintra.txt",status='new')
close(3,status='keep')

open(3,file="./data/connfinal.txt",status='unknown')
close(3,status='delete')
open(3,file="./data/connfinal.txt",status='new')
close(3,status='keep')

open(3,file="./data/times.txt",status='unknown')
close(3,status='delete')
open(3,file="./data/times.txt",status='new')
close(3,status='keep')

END SUBROUTINE clear


! --------- !


SUBROUTINE writer(ncount, gridcount, total, timeReal, V, Vstore, Jc, fire, tstore, ttl)
! The function that writes the values achieved at the end of each cycle:

implicit none

! Variables passed by reference:
integer   :: gridcount, ncount, total, fire(total), ttl
real      :: timeReal, V(total), Jc(total), Vstore(ncount), tstore(ttl, ncount)

! local variables:
integer   :: net, gmod, i, firesum(gridcount), firetotal(ncount), num, neuron
real      :: sumcu(gridcount), current(ncount), totalVoltage(gridcount), voltage(ncount)
real      :: indiVoltages(3 * gridcount), indiJ(3 * gridcount)!, indiFire(3 * gridcount), indiTime(3 * gridcount), indiStore(3 * gridcount)


! Reset all local variables to 0 each loop iteration:
voltage = 0
totalVoltage = 0
current = 0
sumcu = 0
indiJ = 0
firetotal = 0
firesum = 0
indiVoltages = 0


! num is the number of representative neurons needed in total, 3 for each network.
num = 3 * gridcount


! Write to a file the total number of neurons firing in a given network at a given time:
open (3,file="./data/fire.txt",position='append')
  
  nettotal1: do net=1,gridcount
      
    gmod = ncount * (net - 1)
    ido1: do i=1,ncount
      
      firetotal(i) = fire(i + gmod)
      firesum(net) = sum(firetotal)
    
    end do ido1
  
  end do nettotal1
  write (3,*) timeReal, firesum

close (3,status='keep')


! Write to a file the firing status of 3 representative neurons from each network:
open (3,file="./data/indicurrent.txt",position='append')

  do i=1,num
    
    neuron = (total / num) * i
    indiJ(i) = Jc(neuron)
  
  end do
  write(3,*) timeReal, indiJ

close (3,status='keep')


! Write to a file the total current in a given network at a given time:
open (3,file="./data/current.txt",position='append')
  
  nettotal2: do net=1,gridcount
    
    gmod = ncount * (net - 1)
    ido2: do i=1,ncount
      
      current(i) = Jc(i + gmod)
      sumcu(net) = sum(current)
    
    end do ido2
  
  end do nettotal2
  write (3,*) timeReal, sumcu

close (3,status='keep')


! Write to a file the total voltage in a given network at a given time:
open (3, file="./data/voltage.txt",position='append')

  nettotal3: do net=1,gridcount
  
    gmod = ncount * (net - 1)
    ido3: do i=1,ncount
      
      voltage(i) = V(i + gmod)
      totalVoltage(net) = sum(voltage)
    
    end do ido3
  
  end do nettotal3
  write (3,*) timeReal, totalVoltage

close (3, status = 'keep')


! Write to a file the voltage ACHIEVED on 3 representative neurons from each network:
open (3,file="./data/neuronVoltages.txt",position='append')

  do i=1,num
    
    neuron = (total / num) * i
    indiVoltages(i) = V(neuron)
  
  end do
  write(3,*) timeReal, indiVoltages

close (3,status='keep')


! Write to a file the voltage STORED on 3 representative neurons from each network:
!open (3,file="./data/neuronVStore.txt",position='append')

!  do i=1,num
    
!    neuron = (total / num) * i
!    indiStore(i) = Vstore(neuron)
  
!  end do
!  write(3,*) timeReal, indiStore

!close (3,status='keep')


! Write to a file the time since the neuron fired on 3 representative neurons from each network:
!open (3,file="./data/times.txt",position='append')
!
!  do i=1,num
!    
!    neuron = (total / num) * i
!    indiTime(i) = tstore(1,neuron)
!  
!  end do
!  write(3,*) timeReal, indiTime
!
!close (3,status='keep')

END SUBROUTINE writer

end module output_module