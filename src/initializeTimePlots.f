c initializeTimePlots
c####################################################################
      subroutine initializeTimePlots

c--------------------------------------------------------------------
c     Initializes time history plots
c--------------------------------------------------------------------

      use parameters

      use equilibrium

      use graphics

      implicit none

c Call variables

c Local variables

      integer  :: ieq

c Begin program

c Define diagnostics

      diag_desc(IRHO)   = 'ln(drho)'       
      diag_desc(IVX)    = 'ln(dvx)'        
      diag_desc(IVY)    = 'ln(dvy)'        
      diag_desc(IVZ)    = 'ln(dvz)'        
      diag_desc(IBX)    = 'ln(dbx)'        
      diag_desc(IBY)    = 'ln(dby)'        
      diag_desc(IBZ)    = 'ln(dbz)'        
      diag_desc(ITMP)   = 'ln(dtmp)'       
      diag_desc(neqd+1) = 'Magnetic energy'
      diag_desc(neqd+2) = 'Kinetic energy' 
      diag_desc(neqd+3) = 'Thermal energy' 
      diag_desc(neqd+4) = 'Total energy'   
      diag_desc(neqd+5) = 'Time step'      
      diag_desc(neqd+6) = 'Growth rate'    
      diag_desc(neqd+7) = 'div(B)'         
      diag_desc(neqd+8) = 'Total particles'
      diag_desc(neqd+9) = 'Total X momentum'
      diag_desc(neqd+10)= 'Total Y momentum'
      diag_desc(neqd+11)= 'Total Z momentum'

      diag_desc(neqd+12:20) = ''

c Define corresponding independent variables

      diag_ivar(:) = 0  !This means that the independent variable is time
                        !for all variables

c End program

      end subroutine
