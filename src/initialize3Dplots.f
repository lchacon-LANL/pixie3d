c initialize3Dplots
c####################################################################
      subroutine initialize3Dplots

c--------------------------------------------------------------------
c     Set graphics files and dumping intervals
c--------------------------------------------------------------------

      use parameters

      use variables

      use graphics

      use nlfunction_setup

      implicit none

c Call variables

c Local variables

      integer  :: ieq

c Begin program

      do ieq = 1,neqd
        array_graph(ieq)%array => u_np%array_var(ieq)%array
        array_graph(ieq)%descr =  u_np%array_var(ieq)%descr
      enddo
cc      do ieq = 1,neqd
cc        array_graph(ieq)%array => utmp%array_var(ieq)%array
cc        array_graph(ieq)%descr =  utmp%array_var(ieq)%descr
cc      enddo

      array_graph(neqd+1)%array => jx
      array_graph(neqd+1)%descr = 'Jx'

      array_graph(neqd+2)%array => jy
      array_graph(neqd+2)%descr = 'Jy'

      array_graph(neqd+3)%array => jz
      array_graph(neqd+3)%descr = 'Jz'

      array_graph(neqd+4)%array => nuu
      array_graph(neqd+4)%descr = 'nu'

      array_graph(neqd+5)%array => eeta
      array_graph(neqd+5)%descr = 'eta'

      array_graph(neqd+6)%array => divrgB
      array_graph(neqd+6)%descr = 'div(B)'

      array_graph(neqd+7)%array => vx
      array_graph(neqd+7)%descr = 'Vx (car)'

      array_graph(neqd+8)%array => vy
      array_graph(neqd+8)%descr = 'Vy (car)'

      array_graph(neqd+9)%array => vz
      array_graph(neqd+9)%descr = 'Vz (car)'

      array_graph(neqd+10)%array => bx_cov
      array_graph(neqd+10)%descr = 'Bx (car)'

      array_graph(neqd+11)%array => by_cov
      array_graph(neqd+11)%descr = 'By (car)'

      array_graph(neqd+12)%array => bz_cov
      array_graph(neqd+12)%descr = 'Bz (car)'

cc      array_graph(neqd+13:20)%descr = ''

c End program

      end subroutine
