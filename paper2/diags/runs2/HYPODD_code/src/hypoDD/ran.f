c Pseudo-random number generator
c Numerical Recipes, p. 274

	subroutine ran(jlo, jhi, j)

	implicit none

c	Parameters:
	real	jlo, jhi	! Limit values (changed)
	real	j

c	Local variables:
	integer	im, ia, ic
	real	jran

	parameter (im= 714025, ia= 4096, ic= 150889)

	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/ran.f,v 1.4 2001/02/06 03:04:35 julian Exp julian $"/
	save rcsid

      jhi=jhi-1.0
      jran= mod(jran*ia+ic,im)         ! generator
      j= jlo+((jhi-jlo+1)*jran)/im
      return
      end
