        program rad

c---    this routine prompts for information necessary to estimate
c       the radial velocity of a star given from its prominenet
c       absorption features


C--     new edition to be used with CCD spectra instead of the
c       old reticon data.   2/17/93




c*****************************************************************************

        character*32 name,outfile,iname,fname
        character*10 name2(100),id(100),name3(100),iid(100)
	character*10 tname
	character*2 hnm1,hnm2
	character*6 hnm3
        byte object(16),nite(6)
        character*3 name1
	character*9 sfname
        character*9 stname
	character*1 junk
        character*80 line
        integer ipix(40),iindex(500),ibad(100),imark(100),iimark(100)
        integer mdepth,f1depth,f2depth
        real aa(3000,2),arrx(2,4096),arrsub(500),tot(500),cross(100),
     +  cz(100),ccz(100),ms,rlpix(100)
        real arrerr(500),a(6),deltaa(6),sigmaa(6),yfit(500),
     +  rlwidth(100),rlz(100)
        real rindex(500),sn(100),w(100),arrc(4096),arrin(4096),
     +  xp(25),yp(25)
	real czave,czsig,czmean,czcor,rvcor
        real sigmay(12),ai(100),arr(500),dummy(100),xx(100),
     +  xxline(100),yy(100),yyy(100)
        real wline(500),xline(100),xcross(100),eqw(100),dxcross(100)
        real xline3(100),xxline1(100),flg(100)
	byte head(512)
        logical exist
        common /wid/ head,arrx

c*************************************************************************** 
c     0.0  READ IN THE VALUES OF THE LINES TO SEARCH FOR	
c***************************************************************************

        range = 10.		! search range in angstroms
        sncrit=2.		! signal-to-noise cutoff (after continuum fit)
        aicrit=0.3		! asymmetry index cutoff
	ycut = 90.		! minimum value for line to be accepted 
        cee = 2.997929e5	! speed of light
        np = 20                 ! pixel range for features

c-    get beyond header in rvc.out


	open(unit=115,file='rvc.out',status='old')

	do jk = 1,5

	read(115,934) junk

	enddo


	write(*,*) 'Spectral File Name?'
	read(*,934) sfname

	open (unit=111,file=sfname,status='old')

        do ij = 1,500

934     format(a)

        read(111,934,end=936) stname


c***************************************************************************
c     1.0  PROMPT FOR FILES, OPEN THEM, AND SET UP INITIAL CONDITIONS
c***************************************************************************

c---	clear out any leftover values

        do i=1,100
	xx(i) = 0
	yy(i) = 0
	enddo

        num = 0
        nnum = 0
	inum = 0
	keep = 0
	sncrit=2.0


c***************************************************************************
c     2.0  READ IN DATA, CLIP OUTSIDE VALUES 
c***************************************************************************


        open(unit=112,file=stname//'.dad',status='old')
        open(unit=113,file=stname//'.vel',status='unknown')
        open(unit=114,file='rdvel.out',status='unknown')



c-	read in data

	mm = 0

	do i = 1,3000

        read(112,*,end=937) aa(i,1),aa(i,2)
        write(*,*) aa(i,1),aa(i,2)
	arrc(i) = aa(i,2)	

        mm = mm + 1

	enddo

937     continue


	rlmin = aa(1,1)
	rlmax = aa(mm,1)

        write(*,*) rlmin,rlmax

        lst=1

        call features (lst,rlmin,rlmax,nnum,name2,xxline1)

        num = nnum - 1
c        num = nnum

c        write(*,*) 'this is num', num

        do i=1,num
        w(i)=1.
c        xline(i) = xxline1(i+1)
        xline(i) = xxline1(i)
c        write(*,*) 'this is xline', xline(i)
        enddo



c***************************************************************************
c     3.0  MASTER LINE LOOP -- LOOK AT EACH LINE, PERFORM CALCULATIONS
c***************************************************************************


c        write(*,*) '3.0 master line loop'
        do i=1,num
        imark(i)=1        
        enddo

        do 101 i = 1,num
c        ii=2*i-1
c        wline(ii)=xline(i)-range
c        wline(ii+1)=xline(i)+range
c        ipix(ii)=sdpixl (head,wline(ii))
c        ipix(ii+1)=sdpixl (head,wline(ii+1))   
c        np=ipix(ii)-ipix(ii+1)+1
c        jp=ipix(ii+1)


        np = 20

	do j = 1,mm

	if ((aa(j,1).le.xline(i)).and.(aa(j+1,1).ge.xline(i))) then

        jp = j - 10
        goto 987
        endif
        enddo
987     continue

c         write(*,*) 'this be jp', jp
c         write(*,*) aa(j,1)
c---   	obtain position of peak, then reset window centered on peak

c	write(*,*) 'this is np', np

	do kk = 1,2

c---    store data in sub-array starting at index 1

        do j=1,np
        arrsub(j)=arrc(jp+j-1)
        arr(j)=arrsub(j)

c        write(*,*) arr(j)

c        arrerr(j)=arrx(2,jp+j-1)
        rindex(j)=j
        enddo

c--     find value of median data point

        davg = arr(10)

c---    find approximate peak

        ymin=100000.
        do j=1,np
        if (arr(j).lt.ymin) then
        ymin = arr(j)
        imin = j
        endif
        enddo

	if(kk.eq.1) jp = jp+imin-np/2		! reset jp so line is centered

c        write(*,*) 'this ymin. '
c        write(*,*) ymin

	enddo

c***************************************************************************
C     4.0  PERFORM PRELIMINARY GAUSSIAN FIT TO FEATURE
c***************************************************************************

c        write(*,*) 'preforming prelim gauss fit'

c---    use first few pixels on either edge to approximate a continuum

        sum = 0.
        summ = 0.
        do j=1,np/4
        sum = sum +arr(j)
        enddo

        do j=np+1-np/4,np
        summ = summ +arr(j)
        enddo

        cont = (sum + summ)/(np/2)
        del = cont-ymin

c        write(*,*) 'this is cont and del', cont,del

c---    measure s/n in line

        del = abs(del)
        if(cont.le.0.) cont=1			! avoid divide by zero error
        signoise = del/(cont**0.5)
        sn(i)=signoise

c---    if signoise is awful (i.e. less than 2.0) avoid line entirely

        if(signoise.le.2.0) then
        ibad(i)=0
        goto 101
        endif

c---    search for x position where y has fallen to 40% (1-sigma) of the peak

        k=0
        do j=imin,np
        if(arr(j).gt.cont-.6*(del)) k=k+1
        if(k.eq.3) then
        isig = j     
        goto 35
        endif
        enddo

35      xsig = isig-imin

        a(1)=del
        a(2)=imin
        a(3)=xsig
        a(4)=cont

        do j=1,4
        deltaa(j)=a(j)*0.10
        enddo
        deltaa(2)=a(2)*0.01		! pixel location of feature step

        nterms = 4		! continuum plus gaussian (in absorption)

        call gridls (rindex,arr,arrerr,np,nterms,0,a,deltaa,
     +               sigmaa,yfit,chisqr)

        xpix = a(2)+ jp
        sig = a(3)
        cont = a(4)

c***************************************************************************
C     5.0  DETERMINE ASYMMETRY INDEX OF THE LINE
c***************************************************************************

c        write(*,*) 'asymmetry'


c---    find point where value is 80% of peak on either side of line

        k=0
	iright= np
        do j=imin,np
        if(arr(j).gt.cont-.2*(del)) k=k+1
        if(k.eq.3) then
        iright = j     
        endif
        enddo

        k=0
	ileft= 1
        do j=imin,1,-1
        if(arr(j).gt.cont-.2*(del)) k=k+1
        if(k.eq.3) then
        ileft = j     
        endif
        enddo
       
        wl=abs(ileft-imin)
        wr=abs(iright-imin)

c        write(*,*) 'this is asym wl,wr', wl,wr

        if((wl+wr).ne.0) ai(i) = (wl-wr)/(wl+wr)

        if(abs(ai(i)).gt.aicrit) imark(i)=0

c***************************************************************************
C      6.0  OBTAIN SOME ADDITIONAL INFORMATION ABOUT THE LINE
c***************************************************************************


c***************************************************************************
C     7.0  OBTAIN POSITION OF LINE 
c***************************************************************************

c---    obtain sums for analysis

	ms = 0.5	! peak multiplier for gaussian

c        write(*,*) 'this is ms', ms

        call summer (ms,jp,np,arrsub,davg,aa,xpix,sig,iindex,tot,isize)


c        write(*,*) 'this is sig from summer', sig

c---    obtain zero point crossings

        call zero (iindex,tot,isize,xcross(i))

c        write(*,*) 'this is xcross from zero'

101     continue


c--     linearly interpolate in wavelength space to get cross


        do i=1,num

         idex= 0
         idex2= 0


         do j = 1,3000

         if((j.le.xcross(i)).and.(j+1.ge.xcross(i))) then
         idex = j
         idex2 = j+1
         goto 954
         endif
         enddo
954      continue

         tmp1 = aa(idex,1)
         tmp2 = aa(idex2,1)
         dxcross(i) = xcross(i) - idex

         cross(i) = (((tmp2-tmp1)/(idex2-idex))*dxcross(i)) + tmp1


c         write(*,*) 'this is xcross', xcross(i)
c         write(*,*) 'this is tmp1 and tmp2', tmp1,tmp2
c         write(*,*) 'this is cross', cross(i)
939      continue
         enddo


c***************************************************************************
C     8.0  ESTIMATE RADIAL VELOCITY
c***************************************************************************

c        write(*,*) 'estimate radial velocity'

c---    convert crossing positions to wavelength
c       and obtain cz

        sum=0.
        wsum=0.
        wnum=0.

c---    check each velocity one by one and mark those which are bad

2000    do i=1,num
        ibad(i)=1
        if(sn(i).lt.sncrit) ibad(i)=0
        if(abs(ai(i)).gt.aicrit) imark(i)=0
        enddo

        do i=1,num
        if(ibad(i).eq.1) then

        cz(i) = (cross(i)-xline(i))/xline(i)*cee
c        rlz(i) = (rlpix(i)-xline(i))/xline(i)*cee
        wnum=wnum+1
        endif
        enddo

        if(wnum.eq.0) then
        write(*,*)  ' '
        write(*,*) 'NO LINES FOUND !  TRYING AGAIN WITH
     +  LOWER S/N LIMIT '
        write(*,*) ' '
        SNCRIT=SNCRIT-1.0
        goto 2000       
        endif

c***************************************************************************
C     9.0  OMIT LINES ONE AT A TIME TO FIND BADDIES
c***************************************************************************

1500    do j=1,num

        if(ibad(j).eq.0) goto 2222
        sigsum=0.
        sum=0.
        wsum=0.

                k=0
	        do i=1,num
                if(i.ne.j.and.imark(j).eq.1) then
                sum=sum+w(i)*cz(i)
                wsum=wsum+w(i)
                k=k+1
	        endif
	        enddo

        avg=sum
        if(k.gt.0.) avg=sum/(wsum)

        	do i=1,num
	        if(imark(j).eq.1.and.i.ne.j) then
	        sigsum=sigsum+w(i)*(cz(i)-avg)**2
	        endif
	        enddo

        sigest=sigsum**0.5

        if(k-1.gt.0.)sigest=(k*sigsum/((k-1)*wsum))**0.5

        if(abs(cz(j)-avg).gt.3.*sigest) imark(j)=0

2222    enddo

c***************************************************************************
C     10.0  ESTIMATE THE MEAN AND STANDARD DEVIATION
c***************************************************************************


1600    sum=0.
        wnum=0.
        wsum=0.

        do i=1,num
        if(ibad(i).eq.1.and.imark(i).eq.1)then
        sum=sum+w(i)*cz(i)
        wnum=wnum+1
        wsum=wsum+w(i)
        endif
        enddo
        czavg=sum
        if(wnum.ne.0)czavg = sum/wsum

c---    determine scatter
     
        sum2=0.
        do i=1,num
        if(ibad(i).eq.1.and.imark(i).eq.1) then
        sum2 = sum2+w(i)*(cz(i)-czavg)**2
        endif
        enddo

        czsig=sum2**0.5
        if(wnum-1.ne.0..and.wsum.ne.0.)czsig = 
     +  (wnum*sum2/(wsum*(wnum-1.)))**0.5        

        czmean=czsig
        if(wnum.ne.0.)czmean = czsig/(wnum)**0.5

c---    type out results for preliminary inspection

4096    write(*,600)
600     format(/' #    line     position     s/n
     +  velocity    delta    ai   i'/
     +          ' --   ----     --------     ---      --- 
     +    ---    
     +  --------    -----  -----  -')

	sum1 = 0
	mn = 0
	hilim = (abs(czavg)) + czsig
	lolim = (abs(czavg)) - czsig



	do i=1,num
	if((abs(cz(i)).gt.(lolim)).and.(abs(cz(i)).lt.(hilim))) then
	sum1 = sum1 + cz(i)
	mn = mn + 1
	imark(i) = 2
	endif
	enddo


        do i=1,num
        if(ibad(i).eq.1) then
        write(*,601) i,name2(i),cross(i),sn(i),
     +  cz(i),cz(i)-czavg,ai(i),imark(i)
601     format(/1x,i2,2x,
     +  a,1x,f7.1,4x,f5.1,3x,f8.1,2x,f7.1,
     +  3x,f5.2,2x,i1)
        write(113,601) i,name2(i),cross(i),sn(i),
     +  cz(i),cz(i)-czavg,ai(i),imark(i)

        endif
        enddo

503     format (/1x,'mean velocity:  ',f7.1,' km/s',/
     +           1x,'one-sigma scatter:  ',f7.1,' km/s',/
     +           1x,'error in mean:  ',f7.1,' km/s')


c---  through out the values outside of 1 sigma




	czavg = sum1/mn



666     format(25x,f5.1,48x)

	read(115,666) rvcor


	czcor = czavg + rvcor


        write(114,504) stname,mn,czavg,rvcor,czcor,czsig,czmean
504     format(1x,a9,3x,i1,3x,f7.1,3x,f5.1,3x,f7.1,3x,f7.1,3x,f7.1)

518     format(19x,a2,11x,a2,8x,a6)
517     format(1x,a9,6x,f7.1,6x,f7.1,6x,f7.1)

	hnm1 = 'Vo'
	hnm2 = 'Vs'
	hnm3 = 'Vsigma'

	write(113,*)
	write(113,*)
	write(113,*)
	write(113,518) hnm1,hnm2,hnm3
	write(113,*)
	write(113,517) stname,czavg,czcor,czsig
        close(112)
        close(113)

	enddo

936     continue

        close(111)
        close(114)
	close(115)

	end


c***************************************************************************
C     11.0  PROMPT THE USER FOR ACTION
c***************************************************************************

c     RADLIB

c---  This represents a repository of auxilliary routines used for
c     obtaining radial velocities with either RAD.FOR or RAD2.for
c     and should be linked before running.

c*******************************************************************************
c      DETERM
c*******************************************************************************

       function determ (array,norder)

c---   calculates determinant of a square matrix
 
       real*8 array(10,10)

10     determ=1.
11     do 50 k=1,norder

       if(array(k,k)) 41,21,41

21     do 23 j=k,norder
       if(array(k,j)) 31,23,31
23     continue

       determ=0.
       goto 60

31     do 34 i=k,norder
       save=array(i,j)
       array(i,j)=array(i,k)
34     array(i,k)=save

       determ=-determ

41     continue
       determ=determ*array(k,k)
       if(k-norder) 43,50,50
43     k1=k+1

       do 46 i=k1,norder
       do 46 j=k1,norder
46     array(i,j)=array(i,j)-array(i,k)*array(k,j)/array(k,k)
50     continue

60     return
       end

c******************************************************************************(
c	FCHISQ
c**********************************************************************************

       function fchisq (y,sigmay,npts,nfree,mode,yfit)

       real*4 chisq,weight
       dimension y(1),sigmay(1),yfit(1)

11     chisq = 0.
12     if (nfree) 13,13,20
13     fchisq = 0.
       goto 40

c---   accumulate chi square

20     do 30 i=1,npts
21     if (mode) 22,27,29
22     if(y(i)) 25,27,23
23     weight = 1./y(i)
       goto 30
25     weight = 1./(-y(i))
       goto 30
27     weight = 1.
       goto 30
29     weight = 1./sigmay(i)**2
30     chisq = chisq + weight*(y(i)-yfit(i))**2

c---   divide by number of degrees of freedom

31     free = nfree
32     fchisq = chisq/free
40     return
       end

c*******************************************************************************
c       FEATURES
c*******************************************************************************

	subroutine features (lst,rlmin,rlmax,num,name2,xline)

c--- This routine fills up a table of wavelenths and line IDs
c    between the input minumum and maximum wavelengths (RLMIN an
c    RLMAX).  The W array is an arrray of optional line weights
c
c	NOTE:    LST = 0   -- only returns velocity list
c		 LST = 1   -- returns reg absorption list
c                LST = 2   -- returns hot absorption list
c		 LST = 3   -- returns emission list
c                LST = 4   -- returns 0, 1, 2, and extras
c		 LST = 5   -- returns lines for width measurement


	real xline(100),ab(100),em(100),xline3(100)
	real range(100),rrr(100)
	integer ip(100),iw(100),inout(100)
	character*10 name2(100),namea(100),namee(100)
	character*10 name3(100),tname

	common/rr/range

c--- open the data files and read in the data

	num = 0

	if(lst.eq.0.or.lst.eq.1.or.lst.eq.2.or.lst.eq.4) then
	open(unit=10,file='ablines.dat',status='old')
	do i=1,100
	read (10,10,end=100) ab(i),namea(i),iw(i),rrr(i)
	write(*,10) ab(i),namea(i),iw(i),rrr(i)
	num = num + 1	
	enddo
100	CLOSE (UNIT=10)
	else if (lst.eq.3) then
	open(unit=11,file='emlines.dat',status='old')
	do i=1,100
	read (11,10,end=101) em(i),namee(i),iw(i)	
	num = num + 1	
	enddo
101	CLOSE (UNIT=11)
	else if (lst.eq.5) then
	open(unit=12,file='widlines.dat',status='old')
	do i=1,100
	read (12,10,end=102) ab(i),namea(i),iw(i),rrr(i)	
	num = num + 1	
	enddo
102	CLOSE (UNIT=12)
	endif

10	format(f8.3,8x,a10,2x,i2,2x,f4.0)	

c---	sort the line list in increasing order

	do i=1,num
	ip(i) = i
	enddo

	if(lst.ne.3) then
		call sort2 (num,ab,ip)
	else if(lst.eq.3)then
		call sort2 (num,em,ip)
	endif

c--- fill up arrays to return

	if(lst.eq.0) then		! velocity list	
	
	kv=0
	do i=1,num
	if(iw(ip(i)).eq.1.and.(ab(i).ge.rlmin+10
     +    .and.ab(i).le.rlmax-10)) then
	kv = kv + 1
	name2(kv) = namea(ip(i))
        xline(kv) = ab(i)
	endif
	enddo

	num = kv

	else if (lst.eq.1) then		! show absorption

	ka = 0
	do i=1,num
	if(IW(IP(I)).EQ.0.or. iw(ip(i)).eq.1
     +  .AND.(ab(i).ge.rlmin .and. ab(i).le.rlmax)) then
	ka = ka + 1
	name2(ka) = namea(ip(i))
        xline(ka) = ab(i)
	endif
	enddo

	num = ka

	else if (lst.eq.2) then		! show hot absorption

	kh = 0
	do i=1,num
	if(IW(IP(I)).eq.2.AND.(ab(i).ge.rlmin.and.ab(i).le.rlmax)) then
	kh = kh + 1
	name2(kh) = namea(ip(i))
        xline(kh) = ab(i)
	endif
	enddo

	num = kh

	else if (lst.eq.3) then		! show emission

	ke = 0
	do i=1,num
	if(IW(IP(I)).eq.1.AND.(em(i).ge.rlmin.and.em(i).le.rlmax)) then
	ke = ke + 1
	name2(ke) = namee(ip(i))
        xline(ke) = em(i)
	endif
	enddo

	num = ke

	else if (lst.eq.4) then		! show all except emission

	kall = 0
	do i=1,num
	if((IW(IP(I)).EQ.0.or.iw(ip(i)).eq.1.or.iw(ip(i)).eq.2 
     +  .or.iw(ip(i)).eq.3) .AND.(ab(i).ge.rlmin) .and.
     +  (ab(i).le.rlmax)) then
	kall = kall + 1
	name3(kall) = namea(ip(i))
        xline3(kall) = ab(i)
	endif
	enddo

	num = kall

c--- clear out repeats from the complete list

	inout(1) = 1

	do i=1,num
	tname = name3(i)
	do j=i+1,num
	if(tname.eq.name3(j)) then
	inout(j) = -1
	else
	inout(j) = +1
	endif
	enddo
	enddo

	kall = 0
	do i = 1,num
	if (inout(i).eq.1) then
	kall = kall + 1
	name2(kall) = name3(i)
	xline(kall)= xline3(i)
	endif
	enddo

	num = kall

	else if(lst.eq.5) then		! width list	
	
	kw=0
	do i=1,num
	if(iw(ip(i)).eq.1.and.(ab(i).ge.rlmin+10
     +    .and.ab(i).le.rlmax-10)) then
	kw = kw + 1
	name2(kw) = namea(ip(i))
        xline(kw) = ab(i)
	range(kw) = rrr(i)
	endif
	enddo

	num = kw
	endif

	return
	end


c*******************************************************************************
c	FGAUSS
c*******************************************************************************

	function FGAUSS (x,a,y,dyda,ma)

c--- function fits a gaussian in absorption plus a linear

	dimension a(ma),dyda(ma)

        y = a(4)+a(5)*x
        z = (x-a(2))/a(3)
        z2 = z*z
        if(z2-50.) 16,20,20
16      y = y - a(1)*exp(-z2/2.)

	dyda(1) = exp(-z2/2.)	
	dyda(2) = a(1)/a(3)**2*(x-a(2))*exp(-z2/2.)
	dyda(3) = -dyda(2)
	dyda(4) = 1.
	dyda(5) = x

c	type *, function'
c	type *, 'y = ',y
c	type *, 'derivatives'
c	do i=1,6
c       type *, i,dyda(i)
c	enddo

20     return		
	end

c*******************************************************************************
c	FUNCTN
c*******************************************************************************

       function functn (x,i,a)

c---   subroutine called by gridls.for to supply values of the
c      fitting function


       dimension x(4),a(4)

11     x1 = x(i)
12     functn = a(4)
13     z = (x1-a(2))/a(3)
       z2 = z**2
       if(z2-50.) 16,20,20
16     functn = functn-a(1)*exp(-z2/2.)
20     return
       end

c*******************************************************************************
c	FUNCTN2
c*******************************************************************************

       function functn2 (x,i,a)

c---   subroutine called by gridls.for to supply values of the
c      fitting function (Emission features)

c*********************************************************************

       dimension x(4),a(4)

11     x1 = x(i)
12     functn2 = a(4)
13     z = (x1-a(2))/a(3)
       z2 = z**2
       if(z2-50.) 16,20,20
16     functn2 = functn2+a(1)*exp(-z2/2.)
20     return
       end

c*******************************************************************************
c	FUNCTN3
c*******************************************************************************

       function functn3 (x,i,a)

c---   subroutine called by gridls.for to supply values of the
c      fitting function

c*********************************************************************

       dimension x(5),a(5)

11     x1 = x(i)
12     functn3 = a(4)+a(5)*x1
13     z = (x1-a(2))/a(3)
       z2 = z*z
       if(z2-50.) 16,20,20
16     functn3 = functn3-a(1)*exp(-z2/2.)
20     return
       end

c*******************************************************************************c
c	GRIDLS.FOR
c*******************************************************************************

       subroutine gridls (x,y,sigmay,npts,nterms,mode,a,deltaa,
     +                    sigmaa,yfit,chisqr)

c---   copied from Bevington to do a composite least squares fit
c      to spectral data

c************************************************************************

       dimension x(1),y(1),sigmay(1),a(1),deltaa(1),sigmaa(1),yfit(1)


11     nfree = npts-nterms
       free = nfree
       chisqr = 0.
       if (nfree) 100,100,20
20     do 90 j=1,nterms

c---   evaluate chi-square at first search points

21     do 22 i=1,npts
22     yfit(i) = functn(x,i,a)
23     chisq1 = fchisq (y,sigmay,npts,nfree,mode,yfit)
       fn = 0.
       delta = deltaa(j)
41     a(j) = a(j)+delta
       do 43 i=1,npts
43     yfit(i) = functn (x,i,a)
44     chisq2 = fchisq (y,sigmay,npts,nfree,mode,yfit)

c---   exit from program if deltaa/a is too small

       if(deltaa(j)/a(j).lt.0.5e-7) goto 90

45     if (chisq1 - chisq2) 51,41,61

c---   reverse direction of search if chi square is increasing

51     delta = -delta
       a(j) = a(j) + delta
       do 54 i = 1,npts
54     yfit(i) = functn(x,i,a)
       save = chisq1
       chisq1 = chisq2
57     chisq2 = save

c---   increment a(j) until chi-square increases

61     fn = fn + 1.
       a(j) = a(j) + delta
       do 64 i=1,npts
64     yfit(i) = functn(x,i,a)
       chisq3 = fchisq(y,sigmay,npts,nfree,mode,yfit)
66     if (chisq3-chisq2) 71,81,81
71     chisq1 = chisq2
       chisq2 = chisq3
       goto 61

c---   find minimum pf parabola defined by last three points

81     if(chisq3-chisq2.le..5e-7) goto 90
       delta = delta * (1./(1.+(chisq1-chisq2)/
     +         (chisq3-chisq2))+0.5)
82     a(j) = a(j) -delta
83     sigmaa(j) = deltaa(j)*sqrt(2./(free*(chisq3-2.*chisq2+chisq1)))
84     deltaa(j) = deltaa(j)*fn/3.

90     continue

c---   evaluate fit and chi-square for final parameters

91     do 92 i=1,npts
92     yfit(i) = functn(x,i,a)
93     chisqr = fchisq(y,sigmay,npts,nfree,mode,yfit)
100    return
       end            

c********************************************************************************
c	GRIDLS2
c******************************************************************************

       subroutine gridls2 (x,y,sigmay,npts,nterms,mode,a,deltaa,
     +                    sigmaa,yfit,chisqr)

c---   copied from Bevington to do a composite least squares fit
c      to spectral data

c************************************************************************

       dimension x(1),y(1),sigmay(1),a(1),deltaa(1),sigmaa(1),yfit(1)


11     nfree = npts-nterms
       free = nfree
       chisqr = 0.
       if (nfree) 100,100,20
20     do 90 j=1,nterms

c---   evaluate chi-square at first search points

21     do 22 i=1,npts
22     yfit(i) = functn2(x,i,a)
23     chisq1 = fchisq (y,sigmay,npts,nfree,mode,yfit)
       fn = 0.
       delta = deltaa(j)
41     a(j) = a(j)+delta
       do 43 i=1,npts
43     yfit(i) = functn2 (x,i,a)
44     chisq2 = fchisq (y,sigmay,npts,nfree,mode,yfit)

c---   exit from program if deltaa/a is too small

       if(deltaa(j)/a(j).lt..5e-7) goto 90

45     if (chisq1 - chisq2) 51,41,61

c---   reverse direction of search if chi square is increasing

51     delta = -delta
       a(j) = a(j) + delta
       do 54 i = 1,npts
54     yfit(i) = functn2(x,i,a)
       save = chisq1
       chisq1 = chisq2
57     chisq2 = save

c---   increment a(j) until chi-square increases

61     fn = fn + 1.
       a(j) = a(j) + delta
       do 64 i=1,npts
64     yfit(i) = functn2(x,i,a)
       chisq3 = fchisq(y,sigmay,npts,nfree,mode,yfit)
66     if (chisq3-chisq2) 71,81,81
71     chisq1 = chisq2
       chisq2 = chisq3
       goto 61

c---   find minimum pf parabola defined by last three points

81     if(chisq3-chisq2.le..5e-7) goto 90
       delta = delta * (1./(1.+(chisq1-chisq2)/
     +         (chisq3-chisq2))+0.5)
82     a(j) = a(j) -delta
83     sigmaa(j) = deltaa(j)*sqrt(2./(free*(chisq3-2.*chisq2+chisq1)))
84     deltaa(j) = deltaa(j)*fn/3.
90     continue

c---   evaluate fit and chi-square for final parameters

91     do 92 i=1,npts
92     yfit(i) = functn2(x,i,a)
93     chisqr = fchisq(y,sigmay,npts,nfree,mode,yfit)
100    return
       end            

c*******************************************************************************c
c	GRIDLS3.FOR
c*******************************************************************************

       subroutine gridls3 (x,y,sigmay,npts,nterms,mode,a,deltaa,
     +                    sigmaa,yfit,chisqr)

c---   copied from Bevington to do a composite least squares fit
c      to spectral data  -- uses a linear function plus gaussian

c************************************************************************

       dimension x(1),y(1),sigmay(1),a(1),deltaa(1),sigmaa(1),yfit(1)


11     nfree = npts-nterms
       free = nfree
       chisqr = 0.
       if (nfree) 100,100,20
20     do 90 j=1,nterms

c---   evaluate chi-square at first search points

21     do 22 i=1,npts
22     yfit(i) = functn3(x,i,a)
23     chisq1 = fchisq (y,sigmay,npts,nfree,mode,yfit)
       fn = 0.
       delta = deltaa(j)
41     a(j) = a(j)+delta
       do 43 i=1,npts
43     yfit(i) = functn3 (x,i,a)
44     chisq2 = fchisq (y,sigmay,npts,nfree,mode,yfit)

c---   exit from program if deltaa/a is too small

c       if(deltaa(j)/a(j).lt.0.5e-7) goto 90

45     if (chisq1 - chisq2) 51,41,61

c---   reverse direction of search if chi square is increasing

51     delta = -delta
       a(j) = a(j) + delta
       do 54 i = 1,npts
54     yfit(i) = functn3(x,i,a)
       save = chisq1
       chisq1 = chisq2
57     chisq2 = save

c---   increment a(j) until chi-square increases

61     fn = fn + 1.
       a(j) = a(j) + delta
       do 64 i=1,npts
64     yfit(i) = functn3(x,i,a)
       chisq3 = fchisq(y,sigmay,npts,nfree,mode,yfit)
66     if (chisq3-chisq2) 71,81,81
71     chisq1 = chisq2
       chisq2 = chisq3
       goto 61

c---   find minimum pf parabola defined by last three points

81     if(chisq3-chisq2.le..5e-7) goto 90
       delta = delta * (1./(1.+(chisq1-chisq2)/
     +         (chisq3-chisq2))+0.5)
82     a(j) = a(j) -delta
83     sigmaa(j) = deltaa(j)*sqrt(2./(free*(chisq3-2.*chisq2+chisq1)))
84     deltaa(j) = deltaa(j)*fn/3.

90     continue

c---   evaluate fit and chi-square for final parameters

91     do 92 i=1,npts
92     yfit(i) = functn3(x,i,a)
93     chisqr = fchisq(y,sigmay,npts,nfree,mode,yfit)
100    return
       end            

c*******************************************************************************
c	INTERP
c*******************************************************************************

       subroutine interp (x,y,npts,nterms,xin,yout)

c---   Bevington routine to fit and interpolate

       real x(1),y(1)
       real*8 deltax,delta(10),a(10),prod,sum

c---   search for appropriate value of x(1)

11     do 19 i=1,npts
       if(xin-x(i)) 13,17,19
13     i1 = i- nterms/2
       if(i1) 15,15,21
15     i1=1
       goto 21
17     yout = y(i)
18     goto 61
19     continue
       i1=npts-nterms+1
21     i2=i1+nterms-1
       if(npts-i2) 23,31,31
23     i2=npts
       i1 = i2-nterms+1
25     if(i1) 26,26,31
26     i1=1
27     nterms = i2-i1+1

c---   evaluate deviations delta

31     denom = x(i1+1)-x(i1)
       if (denom.le.1.e-8)  denom = 1.e-8
       deltax = (xin-x(i1))/denom
       do 35 i=1,nterms
       ix = i1+i-1
35     delta(i) = (x(ix)-x(i1))/denom

c---   accumulate coefficients a

40     a(1) = y(i1)
41     do 50 k=2,nterms
       prod =1.
       sum = 0.
       imax = k-1
       ixmax = i1+imax
       do 49 i=1,imax
       j=k-i
       prod = prod*(delta(k)-delta(j))
49     sum = sum-a(j)/prod
50     a(k) = sum + y(ixmax)/prod

c---   accumulate sum of expansion

51     sum=a(1)
       do 57, j=2,nterms
       prod=1.
       imax=j-1
       do 56 i=1,imax
56     prod = prod*(deltax-delta(i))
57     sum=sum+a(j)*prod
60     yout=sum
61     return
       end

c******************************************************************************
c	INVERT
c******************************************************************************

c*******************************************************************************
c	POLFIT
c*******************************************************************************

       subroutine polfit (x,y,sigmay,npts,nterms,mode,a,chisqr)

       dimension x(1),y(1),sigmay(1),a(1),sumx(19),sumy(10)
       real*8 array(10,10)

c---   accumulate weighted sums

11     nmax=2*nterms-1
       do 13 n=1,nmax
13     sumx(n)=0.
       do 15 j=1,nterms
15     sumy(j)=0.
       chisq=0.
21     do 50 i=1,npts
       xi=x(i)
       yi=y(i)
31     if (mode) 32,37,39
32     if(yi) 35,37,33
33     weight = 1./yi
       goto 41
35     weight = 1./(-yi)
       goto 41
37     weight = 1.
       goto 41
39     weight = 1./ sigmay(i)**2

41     xterm=weight
       do 44 n=1,nmax
       sumx(n)= sumx(n)+xterm
44     xterm=xterm*xi
45     yterm=weight*yi
       do 48 n=1,nterms
       sumy(n)=sumy(n)+yterm
48     yterm=yterm*xi
49     chisq=chisq+weight*yi**2
50     continue

c---   construct matrices and obtain coefficients

51     do 54 j=1,nterms
       do 54 k=1,nterms
       n=j+k-1
54     array(j,k)=sumx(n)
       delta=determ(array,nterms)
       if(delta) 61,57,61
57     chisqr=0.
       do 59 j=1,nterms
59     a(j)=0.
       go to 80
61     do 70 l=1,nterms
62     do 66 j=1,nterms
       do 65 k=1,nterms
       n=j+k-1
65     array(j,k)=sumx(n)
66     array(j,l)=sumy(j)
       a(l)=determ(array,nterms)/delta
70     continue

c---   calculate chisqr

71     do 75 j=1,nterms
       chisq=chisq - 2.*a(j)*sumy(j)
       do 75 k=1,nterms
       n=j+k-1
75     chisq=chisq+a(j)*a(k)*sumx(n)
76     free=npts-nterms
77     chisqr=chisq/free
80     return
       end

c*******************************************************************************
c       SDFLAT 
c*******************************************************************************

c*******************************************************************************
c	SHOW
c********************************************************************************

c*******************************************************************************
c	SHOW2
c********************************************************************************

C
C.............................................................................
C
c        SUBROUTINE PLOT_VALUES(lst,WLOW,WHIGH,Z,HEAD,BUF)
C
c        REAL*4 W,Z,LINE(100),line2(100),FUDGE,X,WLOW,WHIGH
c        REAL*4 BUF(2,4096),FLUX
c        BYTE HEAD(512),SYMBOL
c        CHARACTER *10 name(100), name2(100)
c        SYMBOL=30
C
C...NIGHT SKY LINES. DO NOT MULTIPLY BY Z!
C
c        LINE(1)=3125.668
c        NAME(1)='HG'
c        LINE(2)=3650.144
c        NAME(2)='HG'
c        LINE(3)=4046.557
c        NAME(3)='HG'
c        LINE(4)=4358.343
c        NAME(4)='HG'
c        LINE(5)=5460.742
c        NAME(5)='HG'
c        LINE(6)=5577.350
c        NAME(6)='[OI]'
c        LINE(7)=6300.23
c        NAME(7)='[OI]'
c        LINE(8)=6363.88
c        NAME(8)='[OI]'
cC
cC...VARIOUS OTHER LINES
c
cC
c        call features (lst,wlow,whigh,num,name2,line2)
c
c        k=0
c        do i=9,num+8
c        k = k + 1
c        name(i)=name2(k)
c        line(i)=line2(k)
c        enddo
c
Cc
Cc...MULTIPLY THE STELLAR LINES BY 1+Z
Cc
c       DO 10 I=9,num+8
c         LINE(I)=LINE(I)*(1.+Z)
c10     CONTINUE
c
C
c        DO 20 I=1,num+8
c          IF ((LINE(I).LT.WLOW).OR.(LINE(I).GT.WHIGH)) GO TO 20
c          X=LINE(I)
c          IPIXEL=SDPIXL(HEAD,X)
c          FLUX=BUF(1,IPIXEL)
cC
cC...LOWER THE Y VALUE BY 20% FOR CRT
cC
c          FLUX=FLUX*0.8
c          CALL PGPOINT(1,X,FLUX,SYMBOL)
cC
cC...LOWER THE Y VALUE BY 10% TO TYPE NAME
cC
c          FLUX=FLUX*0.9
c          CALL PGTEXT(X,FLUX,NAME(I))
c 20     CONTINUE
c        RETURN
c        END


c******************************************************************************
c	SMOOTH
c*******************************************************************************

       subroutine smooth (y,npts)

c---   smooths a set of data by averaging adjacent channels

       real y(1)

11     imax = npts -1
       y1 = y(1)
21     do 24 i=1,imax
       ynew = (y1 + 2*y(i)+y(i+1))/4.
       y1 = y(i)
24     y(i) = ynew
25     y(npts) = (y1 + 3.*y(npts))/4.

       return
       end

c******************************************************************************
c	SMOOTHER
c******************************************************************************

      SUBROUTINE SMOOTHER (M,Y,S,H,R)
C
C SMOOTH DATA: SEE HAYDEN, COMPUTERS IN PHYSICS NOV/DEC 1987, P. 74
C EACH CALL DOES ONE PASS; SET S(I)=0 BEFORE FIRST CALL.
C SUBROUTINE CODED IN FORTRAN BY
C BERNARD GOTTSCHALK, HARVARD CYCLOTRON LAB NOV 87
C
C M   NUMBER OF POINTS
C Y   RAW Y VALUES
C S   SMOOTHED Y VALUES
C H   SCRATCH ARRAY
C R   RMS DEV OF INSTRUMENT FUNCTION (CHANNELS)
C
      PARAMETER (C=.39894228)  !  1/SQRT(2*PI)
C
      DIMENSION Y(M),S(M),H(M)

C

      DO 2400 I=1,M
      JMIN=I-5.*R

C  USE GAUSSIAN ONLY OVER +/- 5 SIGMA

      IF(JMIN .LT. 1) JMIN=1
      JMAX=I+5.*R
      IF(JMAX .GT. M) JMAX=M
      H(I)=S(I)

      DO 2400 J=JMIN,JMAX
      U=I-J
      H(I)=H(I)+(Y(J)-S(J))*C*EXP(-.5*(U/R)**2)/R
2400	continue

      DO 2600 I=1,M
 2600 S(I)=H(I)

      RETURN
C
      END

c******************************************************************************
c	SORT
c******************************************************************************

        subroutine sort (x,n)

c--     performs a simple low to high sort       

        dimension x(n)
        common/extras/id(1000)

        do 10 i=n,2,-1
        do 10 j=2,i
        if(x(j).gt.x(j-1)) goto 10
        temp=x(j)
        itemp=id(j)
        x(j)=x(j-1)
        x(j-1)=temp
        id(j-1)=itemp
10      continue

        return
        end       
 
      SUBROUTINE QSORT(N,RA)

c--- routine to do a heapsort of a data array RA
c--- stolen (unabashedly) from NUMERICAL RECIPES

      DIMENSION RA(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END

c******************************************************************************
c       SORT2
c******************************************************************************

      SUBROUTINE SORT2(N,RA,RB)
        real  RA(N)
	integer rb(n),rrb

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
      END	

c*******************************************************************************
c	SUMMER
c********************************************************************************

      subroutine summer(ms,jp,np,xarr,xavg,aa,xpix,sig,index,tot,isize)

c---    this routine obtains the summation of spectrum convolved
c       with a gaussian derivative which is then solved later
c       for the zero-point crossing
c
c       jp -- starting pixel index
c       np -- number of pixels in data array
c       xarr -- array containing the spectral data
c       xavg -- value of median data point
c       aa -- the total array of spectral values
c       xpix -- mean pixel position
c       sig -- gaussian sigma
c       index -- array containing the pixel positions with sums
c       tot  -- array containing the sums
c       isize -- final size of the tot array
c       iflg -- flags ridiculously bad points

c**************************************************************************

        real xarr (500), tot(500),fit(4000)
        real arrx(2,4096),aa(3000,2),ms
        integer index(500)

        isig = nint(sig)
        ipix = nint(xpix)
        sigstar = ms*sig/(2.**0.5)
        isize = 0

c---    obtain a new data array centered on the approximate
c       line center pixel (within +- 3-sigma)

        mult=3

        jp = ipix -mult*isig
        np = isig*2*mult

        mm = 0

        do i = jp,jp+np
        mm = mm + 1

        arrx(1,mm) = aa(i,2)
        enddo

c---    fill up sub-array and normalize to numbers around 100

        do i=1,np
        xarr(i)=arrx(1,i)/xavg*.01
        enddo

c---    determine maximum and change to an emission feature

        xmax=-1000.
        do i=1,np
        xmax=amax1(xmax,xarr(i))
        enddo

        do i=1,np
        xarr(i)=xmax-xarr(i)
        if(xarr(i).le.0.) xarr(i)=1.e-8
        enddo

c---    set xpix to the lower limit

        ixlo=jp
        ixhi=jp+np


c---    perform summation between ixhi and ixlo for value of xpix
c       ranging over xpix +- mult*sig (to get peak)

        ihi=int(ipix+isig)-int(ipix-isig)

        ipix = ipix - isig -1

        do j=1,ihi
        ipix = ipix+1
        sum=0.
             do i=1,500
             ipos = ixlo+i-1
             delta = ipos-ipix
             f1 = delta*exp(-delta**2./2./sigstar**2.)
             f2 = xarr(i)
             sum = sum +f1*f2
             if(ipos.eq.ixhi) goto 10
             enddo
10      continue

        index(j)=ipix
        tot(j)=sum
        isize=isize+1

        write(*,*) 'this is index',index(j)
        write(*,*) 'this is tot', tot(j)

        enddo

100     return
        end

c******************************************************************************
c	SUMMER2
c********************************************************************************

        subroutine summer2 (jp,np,xarr,xpix,sig,index,tot,isize)

c---    this routine obtains the summation of spectrum convolved
c       with a gaussian derivative which is then solved later
c       for the zero-point crossing
c
c       jp -- starting pixel index
c       np -- number of pixels in data array
c       xarr -- array containing the spectral data
c       pmax -- maximum data value
c       xpix -- mean pixel position
c       sig -- gaussian sigma
c       index -- array containing the pixel positions with sums
c       tot  -- array containing the sums
c       isize -- final size of the tot array

c**************************************************************************

        real xarr (500), tot(500)
        integer index(500)
        real arrx (2,4096)
        common arrx,pmax

        isig = nint(sig)
        ipix = nint(xpix)
        sigstar = sig/(2.**0.5)
        isize = 0

c---    obtain a new data array centered on the approximate
c       line center pixel (within +- 3-sigma)

        mult=3

        jp = ipix -3*mult*isig
        np = 6*isig*mult

        if(jp+np.ge.3744-200..or.jp+np.le.200.) goto 100

c---    fill up sub-array and renormalize to 100

        do i=1,np
        xarr(i) = arrx(1,jp+i-1)/pmax*100.
        enddo

c---    set xpix to the lower limit

        ixlo=jp
        ixhi=jp+np

c---    perform summation between ixhi and ixlo for value of xpix
c       ranging over xpix +- mult*sig (to get peak)

        ihi=int(ipix+mult*isig)-int(ipix-mult*isig)

        ipix = ipix - mult*isig

        do j=1,ihi
        ipix = ipix+1
        sum=0.
             do i=1,500
             ipos = ixlo+i-1
             delta = ipos-ipix
             f1 = delta*exp(-delta**2./2./sigstar**2.)
             f2 = xarr(i)
             sum = sum +f1*f2
             if(ipos.eq.ixhi) goto 10
             enddo
10      continue

        index(j)=ipix
        tot(j)=sum
        isize=isize+1
        enddo

100     return
        end

c*******************************************************************************
c	ZERO
c*******************************************************************************

       subroutine zero (index,tot,isize,xcross)

c---   this routine finds the zero point crossing of
c      a set of data assumed to fit a cubic relation ;
c      also by interpolation to compare results

c***********************************************************

       real tot(500),rindex(500),mult
       integer index(500)

c---   estimate by interpolation

c---   start with best guess as location of gaussian peak
c      and narrow guess region by a binary search criteria

c---   find pixel indices of zero-point straddle position

c----  check if location is between first two pixels

       if(tot(1)*tot(2).lt.0.) then
       k=index(1)
       goto 101
       endif

       do i=1,isize
       mult=tot(i)*tot(i+1)
       if(mult.gt.0.) then
       j=i
       goto 100
       endif
       enddo

100    continue

c---   now find first one with the opposite sign

       do i=j,isize-1
       mult=tot(i)*tot(i+1)       
       if(mult.lt.0.) then
       k=index(i)
       goto 101
       endif
       enddo

101    do i=1,isize
       rindex(i)=index(i)
       enddo

       x1 = k
       x2 = k+1
       x3 = k+0.5
       nterms = 4		! cubic fit

       do i=1,20
       call interp (rindex,tot,isize,nterms,x3,yout)
       if(yout.lt.0.) then
       s=-1.
       else if(yout.gt.0.) then
       s=1.
       endif
       if (i.eq.20) goto 30
       dx = s*abs(x2-x3)
       x2=x3
       x3=x3+dx/2.
       enddo

30     xcross = x3

c       type *, 'accuracy:  ',yout

       return
       end
  

      SUBROUTINE MRQCOF(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NALP,CH
     *ISQ,FUN)
      PARAMETER (MMAX=20)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),ALPHA(NALP,NALP),BETA(MA),
     *    DYDA(MMAX),LISTA(MFIT)
      DO 12 J=1,MFIT
        DO 11 K=1,J
          ALPHA(J,K)=0.
11      CONTINUE
        BETA(J)=0.
12    CONTINUE
      CHISQ=0.
      DO 15 I=1,NDATA
c        CALL FGauss(X(I),A,YMOD,DYDA,MA)
        SIG2I=1./(SIG(I)*SIG(I))
        DY=Y(I)-YMOD
        DO 14 J=1,MFIT
          WT=DYDA(LISTA(J))*SIG2I
          DO 13 K=1,J
            ALPHA(J,K)=ALPHA(J,K)+WT*DYDA(LISTA(K))
13        CONTINUE
          BETA(J)=BETA(J)+DY*WT
14      CONTINUE
        CHISQ=CHISQ+DY*DY*SIG2I
15    CONTINUE
      DO 17 J=2,MFIT
        DO 16 K=1,J-1
          ALPHA(K,J)=ALPHA(J,K)
16      CONTINUE
17    CONTINUE
      RETURN
      END

      SUBROUTINE MRQMIN(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,
     *    COVAR,ALPHA,NCA,CHISQ,FUN,ALAMDA)
      PARAMETER (MMAX=20)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),LISTA(MFIT),
     *  COVAR(NCA,NCA),ALPHA(NCA,NCA),ATRY(MMAX),BETA(MMAX),DA(MMAX)
      IF(ALAMDA.LT.0.)THEN
        KK=MFIT+1
        DO 12 J=1,MA
          IHIT=0
          DO 11 K=1,MFIT
            IF(LISTA(K).EQ.J)IHIT=IHIT+1
11        CONTINUE
          IF (IHIT.EQ.0) THEN
            LISTA(KK)=J
            KK=KK+1
          ELSE IF (IHIT.GT.1) THEN
            PAUSE 'Improper permutation in LISTA'
          ENDIF
12      CONTINUE
        IF (KK.NE.(MA+1)) PAUSE 'Improper permutation in LISTA'
        ALAMDA=0.001
        CALL MRQCOF(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NCA,CHISQ,F
     *UN)
        OCHISQ=CHISQ
        DO 13 J=1,MA
          ATRY(J)=A(J)
13      CONTINUE
      ENDIF
      DO 15 J=1,MFIT
        DO 14 K=1,MFIT
          COVAR(J,K)=ALPHA(J,K)
14      CONTINUE
        COVAR(J,J)=ALPHA(J,J)*(1.+ALAMDA)
        DA(J)=BETA(J)
15    CONTINUE
      CALL GAUSSJ(COVAR,MFIT,NCA,DA,1,1)
      IF(ALAMDA.EQ.0.)THEN
        CALL COVSRT(COVAR,NCA,MA,LISTA,MFIT)
        RETURN
      ENDIF
      DO 16 J=1,MFIT
        ATRY(LISTA(J))=ATRY(LISTA(J))+DA(J)
16    CONTINUE
      CALL MRQCOF(X,Y,SIG,NDATA,ATRY,MA,LISTA,MFIT,COVAR,DA,NCA,CHISQ,FU
     *N)
      IF(CHISQ.LT.OCHISQ)THEN
        ALAMDA=0.1*ALAMDA
        OCHISQ=CHISQ
        DO 18 J=1,MFIT
          DO 17 K=1,MFIT
            ALPHA(J,K)=COVAR(J,K)
17        CONTINUE
          BETA(J)=DA(J)
          A(LISTA(J))=ATRY(LISTA(J))
18      CONTINUE
      ELSE
        ALAMDA=10.*ALAMDA
        CHISQ=OCHISQ
      ENDIF
      RETURN
      END

      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END

      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      DIMENSION COVAR(NCVM,NCVM),LISTA(MFIT)
      DO 12 J=1,MA-1
        DO 11 I=J+1,MA
          COVAR(I,J)=0.
11      CONTINUE
12    CONTINUE
      DO 14 I=1,MFIT-1
        DO 13 J=I+1,MFIT
          IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
          ENDIF
13      CONTINUE
14    CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
        COVAR(1,J)=COVAR(J,J)
        COVAR(J,J)=0.
15    CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
        COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16    CONTINUE
      DO 18 J=2,MA
        DO 17 I=1,J-1
          COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18    CONTINUE
      RETURN
      END

	

      FUNCTION BRENT(AX,BX,CX,NPOL,AA,TOL,XMIN)
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0E-10)
      EXTERNAL POLY
      REAL AA(NPOL)

      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
      FX=POLY(AA,NPOL,X)
      FV=FX
      FW=FX
      DO 11 ITER=1,ITMAX
        XM=0.5*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.*TOL1
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3
        IF(ABS(E).GT.TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.*(Q-R)
          IF(Q.GT.0.) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR. 
     *        P.GE.Q*(B-X)) GOTO 1
          D=P/Q
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
1       IF(X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
2       IF(ABS(D).GE.TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=POLY(AA,NPOL,U)
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
11    CONTINUE
      PAUSE 'Brent exceed maximum iterations.'
3     XMIN=X
      BRENT=FX
      RETURN
      END

      SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)
      PARAMETER (NMAX=100)
      DIMENSION U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
      DO 12 J=1,N
        S=0.
        IF(W(J).NE.0.)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=0.
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      RETURN
      END
      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
      PARAMETER (NMAX=100)
      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
      G=0.0
      SCALE=0.0
      ANORM=0.0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0
        S=0.0
        SCALE=0.0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0
        S=0.0
        SCALE=0.0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0
            V(J,I)=0.0
31        CONTINUE
        ENDIF
        V(I,I)=1.0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0) THEN
          G=1.0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0
          S=1.0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=SQRT(F*F+G*G)
              W(I)=H
              H=1.0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.30) PAUSE 'No convergence in 30 iterations'
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G=SQRT(F*F+1.0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=1.0
          S=1.0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=SQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=SQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.0.0) THEN
              Z=1.0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END

      SUBROUTINE SVDFIT(X,Y,SIG,NDATA,A,MA,U,V,W,MP,NP,CHISQ,FUNCS)
      PARAMETER(NMAX=1000,MMAX=50,TOL=1.E-5)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),V(NP,NP),
     *    U(MP,NP),W(NP),B(NMAX),AFUNC(MMAX)
      DO 12 I=1,NDATA
c        CALL FUNCS(X(I),AFUNC,MA)
        TMP=1./SIG(I)
        DO 11 J=1,MA
          U(I,J)=AFUNC(J)*TMP
11      CONTINUE
        B(I)=Y(I)*TMP
12    CONTINUE
      CALL SVDCMP(U,NDATA,MA,MP,NP,W,V)
      WMAX=0.
      DO 13 J=1,MA
        IF(W(J).GT.WMAX)WMAX=W(J)
13    CONTINUE
      THRESH=TOL*WMAX
      DO 14 J=1,MA
        IF(W(J).LT.THRESH)W(J)=0.
14    CONTINUE
      CALL SVBKSB(U,W,V,NDATA,MA,MP,NP,B,A)
      CHISQ=0.
      DO 16 I=1,NDATA
c        CALL FUNCS(X(I),AFUNC,MA)
        SUM=0.
        DO 15 J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
15      CONTINUE
        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
16    CONTINUE
      RETURN
      END

      SUBROUTINE SVDVAR(V,MA,NP,W,CVM,NCVM)
      PARAMETER (MMAX=20)
      DIMENSION V(NP,NP),W(NP),CVM(NCVM,NCVM),WTI(MMAX)
      DO 11 I=1,MA
        WTI(I)=0.
        IF(W(I).NE.0.) WTI(I)=1./(W(I)*W(I))
11    CONTINUE
      DO 14 I=1,MA
        DO 13 J=1,I
          SUM=0.
          DO 12 K=1,MA
            SUM=SUM+V(I,K)*V(J,K)*WTI(K)
12        CONTINUE
          CVM(I,J)=SUM
          CVM(J,I)=SUM
13      CONTINUE
14    CONTINUE
      RETURN
      END
 
        FUNCTION POLY (A,N,X)
 
        REAL A(N)

	POLY = A(N) 
        DO 10 J=N-1,1,-1
        POLY = POLY*X+A(J)
10      CONTINUE
 
        RETURN
        END
 
      SUBROUTINE FPOLY(X,P,NP)
      DIMENSION P(NP)
      P(1)=1.
      DO 11 J=2,NP
        P(J)=P(J-1)*X
11    CONTINUE
      RETURN
      END

      SUBROUTINE FLEG(X,PL,NL)
      DIMENSION PL(NL)
      PL(1)=1.
      PL(2)=X
      IF(NL.GT.2) THEN
        TWOX=2.*X
        F2=X
        D=1.
        DO 11 J=3,NL
          F1=D
          F2=F2+TWOX
          D=D+1.
          PL(J)=(F2*PL(J-1)-F1*PL(J-2))/D
11      CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE SMOOFT(Y,N,PTS)
      PARAMETER(MMAX=1024)
      DIMENSION Y(MMAX)
      M=2
      NMIN=N+2.*PTS
1     IF(M.LT.NMIN)THEN
        M=2*M
      GO TO 1
      ENDIF
      CONST=(PTS/M)**2
      Y1=Y(1)
      YN=Y(N)
      RN1=1./(N-1.)
      DO 11 J=1,N
        Y(J)=Y(J)-RN1*(Y1*(N-J)+YN*(J-1))
11    CONTINUE
      IF(N+1.LE.M)THEN
        DO 12 J=N+1,M
          Y(J)=0.
12      CONTINUE
      ENDIF
      MO2=M/2
      CALL REALFT(Y,MO2,1)
      Y(1)=Y(1)/MO2
      FAC=1.
      DO 13 J=1,MO2-1
        K=2*J+1
        IF(FAC.NE.0.)THEN
          FAC=AMAX1(0.,(1.-CONST*J**2)/MO2)
          Y(K)=FAC*Y(K)
          Y(K+1)=FAC*Y(K+1)
        ELSE
          Y(K)=0.
          Y(K+1)=0.
        ENDIF
13    CONTINUE
      FAC=AMAX1(0.,(1.-0.25*PTS**2)/MO2)
      Y(2)=FAC*Y(2)
      CALL REALFT(Y,MO2,-1)
      DO 14 J=1,N
        Y(J)=RN1*(Y1*(N-J)+YN*(J-1))+Y(J)
14    CONTINUE
      RETURN
      END
      SUBROUTINE SNCNDN(UU,EMMC,SN,CN,DN)
      PARAMETER (CA=.0003)
      LOGICAL BO
      DIMENSION EM(13),EN(13)
      EMC=EMMC
      U=UU
      IF(EMC.NE.0.)THEN
        BO=(EMC.LT.0.)
        IF(BO)THEN
          D=1.-EMC
          EMC=-EMC/D
          D=SQRT(D)
          U=D*U
        ENDIF
        A=1.
        DN=1.
        DO 11 I=1,13
          L=I
          EM(I)=A
          EMC=SQRT(EMC)
          EN(I)=EMC
          C=0.5*(A+EMC)
          IF(ABS(A-EMC).LE.CA*A)GO TO 1
          EMC=A*EMC
          A=C
11      CONTINUE
1       U=C*U
        SN=SIN(U)
        CN=COS(U)
        IF(SN.EQ.0.)GO TO 2
        A=CN/SN
        C=A*C
        DO 12 II=L,1,-1
          B=EM(II)
          A=C*A
          C=DN*C
          DN=(EN(II)+A)/(B+A)
          A=C/B
12      CONTINUE
        A=1./SQRT(C**2+1.)
        IF(SN.LT.0.)THEN
          SN=-A
        ELSE
          SN=A
        ENDIF
        CN=C*SN
2       IF(BO)THEN
          A=DN
          DN=CN
          CN=A
          SN=SN/D
        ENDIF
      ELSE
        CN=1./COSH(U)
        DN=CN
        SN=TANH(U)
      ENDIF
      RETURN

      END
      SUBROUTINE REALFT(DATA,N,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      THETA=6.28318530717959D0/2.0D0/DBLE(N)
      WR=1.0D0
      WI=0.0D0
      C1=0.5
      IF (ISIGN.EQ.1) THEN
        C2=-0.5
        CALL FOUR1(DATA,N,+1)
        DATA(2*N+1)=DATA(1)
        DATA(2*N+2)=DATA(2)
      ELSE
        C2=0.5
        THETA=-THETA
        DATA(2*N+1)=DATA(2)
        DATA(2*N+2)=0.0
        DATA(2)=0.0
      ENDIF
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=SIN(THETA)
      N2P3=2*N+3
      DO 11 I=1,N/2+1
        I1=2*I-1
        I2=I1+1
        I3=N2P3-I2
        I4=I3+1
        WRS=SNGL(WR)
        WIS=SNGL(WI)
        H1R=C1*(DATA(I1)+DATA(I3))
        H1I=C1*(DATA(I2)-DATA(I4))
        H2R=-C2*(DATA(I2)+DATA(I4))
        H2I=C2*(DATA(I1)-DATA(I3))
        DATA(I1)=H1R+WRS*H2R-WIS*H2I
        DATA(I2)=H1I+WRS*H2I+WIS*H2R
        DATA(I3)=H1R-WRS*H2R+WIS*H2I
        DATA(I4)=-H1I+WRS*H2I+WIS*H2R
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
11    CONTINUE
      IF (ISIGN.EQ.1) THEN
        DATA(2)=DATA(2*N+1)
      ELSE
        CALL FOUR1(DATA,N,-1)
      ENDIF
      RETURN
      END

      SUBROUTINE FOUR1(DATA,NN,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*SIN(0.5D0*THETA)**2
        WPI=SIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END

