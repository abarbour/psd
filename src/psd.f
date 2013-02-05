                      program psd
c
c  Performs the PSD estimation by the sine multitaper method of
c  Riedel & Sidorenko, IEEE Tr Sig Pr,43, 188, 1995
c  See ascan for a synopsis of commands
c
c$$$$ calls autoco getdat getone getint getpsd plot prinz ascan stats
c
c
c  Common blocks of ascan:
      parameter (inmx=500)
      character *80 input(inmx)
      common /store/ input
      common /ndict/ iecho,nin
c
c  Common blocks of getdat, state:
      common /state/ xbar,varx,smin,smax,ebar,kmin,kmax
c
c  Notice: to increase the length of series that can be treated,
c  simply raise mxx in every routine.
      parameter (mxx=500 000)
      common /series/ nx,dt,t1,x(mxx)
c  Common block from getpsd:
      parameter (nfmx=mxx + 2)
      common /result/ nt,nf,fNyq,df,sx(nfmx),kopt(nfmx)
c
      parameter (pi=3.14159265359)
c
ce
      dt=1.0
      nx=0
c
c  Read in the commands
      do 1000 kase=1, 10000
        call ascan(kase)
        call getint('quit', ignor, kwit)
        if (kwit.ge. 0) then
           dt=0.0
c  Flush all plots on 'hold'
           call getint('hold', ignor, kold)
           if (kold.ge. 0) call plot
           write(*,'(/a)') 'Normal termination'
           stop
        endif
c
c
c  Replot on another scale spectra from previous calculations
        call getone('replot',  dum, jagain)
        if (jagain.eq. 1 .and. nx.gt.0) then
          call plot
          goto 1000
        endif
c
c
c  Input the data series; find and print summary statistics
        call getdat
        call stats
c
c
c  Compute the PSD
        call getpsd
c
c
c  Print the results and save them to disk
        call prinz
c
c  Make a plotfile for the results and plot to screen if required
        call plot
c
c
c  Find and plot autocorrelation function
        call autoco
c
c
        write(*,'(/a)')'Enter further commands or "quit" to finish'
 1000 continue
      end
c_______________________________________________________________
      subroutine prinz
c$$$$ calls units getchr ljust
c  Prints summary of results and writes full tables to a diskfile
c
      common /state/ xbar,varx,smin,smax,ebar,kmin,kmax
      parameter (mxx=500 000)
      common /series/ nx,dt,t1,x(mxx)
      parameter (nfmx=mxx + 2)
      common /result/ nt,nf,fNyq,df,sx(nfmx),kopt(nfmx)
c
      parameter (pi=3.14159265359)
      character*64 unit,spout,fruq*20
      data unit/' '/, spout/'fort.1'/
      data i1,i2,i3/1,1,1/
c
ce
      if (nx.le. 1) return
c  Print summary of useful results
      call units(i1,i2,i3, unit)
      fruq=unit(1:i1)
      if (fruq.eq. '/s ' .or. fruq.eq.'sec') fruq='Hz'
      write(*,'(/(a,g12.4,1x,a))')
     $'           Nyquist frequency:',fNyq,fruq,
     $'           Frequency spacing:',df,  fruq
      write(*,'(/(a,i7))')
     $'Length of time series analysed:',2*nf-2,
     $'         Number of frequencies:',nf
c
c  Write results to a file ONLY when output command is issued
      call getchr('output', spout, named)
      if (spout.eq.'-none') named=-1
      if (named.ge. 0) open (unit=1, file=spout)
c
c  Write to unit 1: f, S(f), e, res, Kopt
c  where e=root variance, res=frequency resolution
      ttime=dt*(2*nf - 2)
      ebar=0.0
      kbar=0.0
      kmax=0
      kmin=nf
      do 1300 j=1, nf
        e=sx(j) / sqrt(real(kopt(j))/1.2)
        ebar=ebar + e/sx(j)
        kbar=kbar + kopt(j)
        kmax=max(kopt(j), kmax)
        kmin=min(kopt(j), kmin)
        if (named.ge. 0)
     $  write(1,*) df*(j-1),sx(j),e,kopt(j)/ttime,kopt(j)
 1300 continue
      ebar=ebar/nf
      kbar=kbar/nf
c
      write(*,'(a,3(i5,a)/a,f6.3,a)')
     $'      Average number of tapers: ',kbar,
     $'   [min',kmin,'     max',kmax,']',
     $'         Average 1-sigma error:',ebar,' * Sx(f)'
      write(*,'(a,g12.4,1x,a)')'  Typical frequency resolution:',
     $  fNyq /(nf*ebar**2), fruq
c
      if (named.ge. 0) then
        call ljust(spout,spout, np)
        write(*,'(/2a/a)')'One-sided PSD written to: ',spout(1:np),
     $ 'The file contains:  freq, PSD, err, resolution, # tapers,'
      else
        write(*,'(/a)')'This spectrum was NOT written to disk'
      endif
      return
      end
c_______________________________________________________________
      subroutine getpsd
c$$$$ calls getint getone getchr plot prolat quick yule adapt nearn later
c$$$$ ljust period
c  Performs the PSD estimation by the sine multitaper method of
c  Riedel & Sidorenko, IEEE Tr Sig Pr,43, 188, 1995
c  May invoke prewhitening via yule - Yule-Walker equations
      parameter (mxx=500 000, mpool=8*mxx+4)
      parameter (pi=3.141592653589793d0)
      parameter (nfmx=mxx + 2)
      common /series/ nx,dt,t1,x(mxx)
      common /state/ xbar,varx,smin,smax,ebar,kmin,kmax
      common /result/ nt,nf,fNyq,df,sx(nfmx),kopt(nfmx)
c
      character*72 taped,white,line, signal*64
      data  ntimes/4/,ktap/0/, taped/'*'/, white/'fort.13'/
c
c  Adjust series length if necessary to be 2**i * 3**j * 5**k
c  with i > 0, j,k >= 0.  Find other constants
ce
      if (nx.le. 1) then
        write(*,'(a)')
     $ '>>> Too few terms in the series to estimate a spectrum'
        return
      endif
      nt=2*nearn(nx/2)
      nf=1 + nt/2
      nt=2*nf - 2
      fNyq=0.5/dt
      df=fNyq/(nf - 1)
c
      nord=0
      call getint('prewhiten', nord, ignor)
c  Prewhiten series with FIR filter; but first find
c  a pilot spectral estimate.  On returning from yule x() has
c  been modified.
      if (nord.gt. 1) then
        write(*,'(/5x,a,i4)')
     $ 'Series prewhitened by FIR filter, order ',nord
c
c  Pilot estimate of PSD then prewhiten time series
        call quick(1, -nint(5.0 + 0.2*sqrt(real(nx))))
        call yule(1, nord)
      elseif (nord.eq. -1) then
        call yule(3, nord)
        write(*,'(a,i4)')
     $ 'Series prewhitened by prior FIR filter, order ',nord
      endif
c
c  Save prewhitened time series to diskfile
      white='fort.13'
      call getchr('save', white, kule)
      if (kule.ge. 0) then
        open(unit=13, file=white)
        write(13,'(1p,2g14.5)') (t1+(j-1)*dt,x(j),j=1,nx)
        nwhi=max(kule,8)
        write(*,'(a,a)')'Prewhitened series saved to: ',white(1:nwhi)
      endif
c
c       --------Periodogram Estimate--------------
c
      call getint('tapers', ntap, nu)
      if (nu.le. 0) goto 1000
      if (ntap.le. 0) then
        call period
        write(*,*)' Zero tapers requested: periodogram returned'
        goto 3000
      endif
c
c       --------Prolate Spheroidal multitapers--------------
c
 1000 call getchr('prolate',  taped, nprole)
      if (nprole.lt. 0) goto 2000
c
      call getone('bandwidth', width, nband)
      call getone('timebp',    timebp, ntime)
      if (max(nband, ntime).lt. 1) then
        write(*,'(a)')
     $ '>>> Bandwidth not provided - PSD will use sine tapers'
        goto 2000
      endif
      if (nband.gt. 0) timebp=0.5*nt*width/fNyq
      write(*,'(/a,g12.4)')
     $'Time-bandwidth product for prolate tapers:',timebp
      ntap=1.1*timebp
      call getint('tapers', ntap, ignor)
c  Save tapers to a diskfile if named in command
      ktap=0
      if (nprole.gt. 0) then
        open(unit=16, file=taped)
        ktap=16
      endif
c
      call prolat(ntap, timebp, ktap)
c
      call ljust(taped,taped, npr)
      if (ktap.ne.0)
     $write(*,'(3x,2a)') 'Tapers written to file: ',taped(1:npr)
      goto  3000
c
c       --------Sine multitapers--------
c
c  Number of tapers is requested for uniform tapering or to
c  place a bound on their number (less >0)
 2000 call getchr('tapers', line, ignor)
      call getint('tapers', ntap, iftap)
      less = index(line,'<')
      nbound=0
      if (iftap.gt.0 .and. less.gt.0) nbound=ntap
      call getint('adapt', ntimes, ignor)
c
c  Decide whether to make a simple estimate or an adaptive one
      if (later('tapers', 'adapt').gt. 0 .and. less.eq. 0
     $     .and. iftap.gt.0) then
c  Estimate S with uniform number of tapers across spectrum
        write(*,'(/a,i6,a/a,f6.1)')
     $ 'Uniform spectral estimation with',ntap,' tapers',
     $ '     viz Time-bandwidth product:',0.5*ntap
        call quick(1, -ntap)
c
      else
c  Make adaptive estimate
        fact=1.0
        call getone('smooth', fact, ignor)
        initap = 3.0 + sqrt(fact*nx)/5.0
        write(*,'(/a,i5,a/(a,i5))')
     $ 'Adaptive spectral estimation with',ntimes,' iterations',
     $ '         Initial number of tapers',initap
          if (fact .ne. 1.0) write(*,'(a,f5.2)')
     $ '   Smooth penalizes var/bias^2 by ',fact**2
c
        call adapt(ntimes, initap, nbound)
c
      endif
c
c
c  Undo the prewhitening, if any; but first, plot prewhiten spectrum
 3000 if (nord.gt. 1) then
        call getchr('prewhiten', signal, ignor)
        if (index(signal, 'plot').gt. 0) call plot
        call yule(2, nord)
      endif
c
c
c  Normalize the spectrum
c  Various constants calculated spectrum scaled to correct units
      const=dt/nx
      smax=sx(1)*const
      smin=smax
      do 5000 k=1, nf
        sx(k)=const*sx(k)
        smax=max(sx(k), smax)
        smin=min(sx(k), smin)
 5000 continue
c
      return
c
      end
c______________________________________________________________________
      subroutine plot
c$$$$ calls units getchr getint clear getarr getone system ljust
c
c  Creates plotxy plotfile of current spectrum, or of all
c  spectra.  Invokes plotxy via 'system' call
c
c  Common blocks of ascan:
      parameter (inmx=500)
      character *80 input(inmx)
      common /store/ input
      common /ndict/ iecho,nin
c
c  Common blocks of getdat, state
      parameter (mxx=500 000)
      common /series/ nx,dt,t1,x(mxx)
      common /state/ xbar,varx,smin,smax,ebar,kmin,kmax
c
c  Common block from PSD
      parameter (nfmx=mxx + 2)
      common /result/ nt,nf,fNyq,df,sx(nfmx),kopt(nfmx)
c
      character*64 unit,splot,title,detl, code*4,bks*1
      dimension caylim(2)
      data unit/' '/,  splot/'fort.2'/, title/' '/
      data i1,i2,i3/1,1,1/, nplot/0/, ncolor/0/
      data caylim /0.0, 0.0/, dsh/0.0/,  llog/0/
      save nplot, slim
c
c  Because of special status of \ in Unix, generate it from ASCII
ce
      bks=char(92)
c  Special flag to flush the 'hold' plot buffer before quitting
      if (dt.eq. 0.0) goto 2000
      if (nx.le. 1) return
c  Get units of the series and frequency
      call units(i1,i2,i3, unit)
c
c  Make a plotfile for results
c
      call getchr('plot', splot, newplt)
c
c  Omit if splot = 'off' or '-none'
      if (splot .eq. 'off' .or. splot.eq.'-none') return
c
      open (unit=2, file=splot)
      call clear ('plot')
      if (newplt.ge. 0) nplot=0
      nplot=nplot + 1
c  User wants this plot on top of previous one
      ntop =-1
      if (nplot.gt. 1) call getone('super', dum, ntop)
c
c  Rewrite the plotfile
      rewind (unit=2)
      if (ntop.ge. 0) then
         do 1500 lin=1, 100000
           read (2, '(a4)', err=1510, end=1520) code
           if (code .eq. 'plot') then
             backspace (unit=2)
             goto 1520
           endif
 1500    continue
 1510    write(*,'(2a)')' ','>>> Cannot superimpose spectra - ',
     $   'No prior graph'
 1520    continue
         ncolor=1 + ncolor
      else
         ncolor=0
         nplot=1
      endif
      dsh=0.03*(ncolor/9)
c
c  See if a subinterval of freq axis is requested
      caylim(1)=-fNyq
      caylim(2)= fNyq
      llog=0
      call getarr('detail', caylim, 2, ndet)
      call getchr('detail', detl, kog)
      llog=index(detl, 'log')*kog
      call getarr('replot', caylim, 2, ndet)
      call clear('replot')
      nf1= max(1.0, caylim(1)/df + 1.0)
      nf2= min(real(nf), caylim(2)/df+ 1.0)
      if (nf1 .ge. nf2) then
         nf1=1
         nf2=nf
      endif
c  Log frequency axis?
      call getint('logfreq', ignor, logf)
      if (logf.ge. 0) llog=1
      if (llog.gt. 0) nf1=max(nf1,2)
      if (nplot.le. 1) then
        slim=0.0
        do  1600 j=nf1, nf2
           slim=max(slim, sx(j))
 1600   continue
      endif
c
c  Get the title if there is one
      call getchr('title', title, nti)
      if (nti.eq. -99) call getchr('file', title, nti)
c
c  Write plotfile to unit 2
      if (nti.gt. 0) write(2,'(2a)') 'title ',title(1:nti)
      write(2,'(a)')'logxy linlog','frame gridsol','char 0.11'
      if (unit(1:1).eq. ' ' .and. dt.eq. 1.0) then
        write(2,'(5a)') 'xlabel Frequency in 1/',bks,'D',bks,'t'
      elseif (unit(1:1).eq. ' ') then
        write(2, '(a)') 'xlabel Frequency '
      elseif (unit(2:i1).eq. 's' .or. unit(2:i1).eq.'sec') then
        write(2, '(a)') 'xlabel Frequency Hz'
      else
        write(2,'(4a)')'xlabel Frequency ',unit(2:i1-1),bks,'sup{-1}'
      endif
      if (i2 .le. i1+1) then
        write(2,'(a)') 'ylabel PSD in arbitrary units'
      else
        write(2,'(5a)')
     $    'ylabel PSD ',unit(i1:i2-1),bks,'sup{2 }',unit(2:i1)
      endif
c  Write out the spectrum
      write(2,'(a,2f7.3)') 'dash ',dsh,0.05
      if (nplot.ge. 2) write(2,'(a,i4)') 'color ',1+mod(ncolor,9)
      write(2,'(a/a,4g12.4)') 'mode 1',
     $  'affine ',df,(nf1-2)*df,1.0,0.0
      write(2,'(a,i7/(1p,7e10.3))') 'read ',nf2-nf1+1,
     $  (sx(k), k=nf1, nf2)
      write(2,'(a)')'affine 1 0 1 0','ylim 5.2','xlim 4.6'
c  Draw error bars
      psig=1.03*slim
      pmid=psig +  max(smin, psig-psig/sqrt(real(kmin)))
      fdel= -df*(nf2 - nf1)*0.05
      fbar= fdel*(nplot+1) + df*nf2
      if (llog.eq.1 ) fbar= df*nf1*real(nf2/nf1)**(0.95-0.05*nplot)
      if (llog.eq.1 ) fdel= -0.05*fbar*log(real(nf2/nf1))
      write (2,'(a/a/a/3g15.7/a)') 'symb dot 0.05','mode 3','read 1',
     $  fbar+fdel*.2,1.05*pmid,psig/sqrt(real(kmax)),'mode 2'
      write (2,'(a/a/a/3g15.7/a)') 'symb dot 0.05','mode 3','read 1',
     $  fbar,1.05*pmid,psig/sqrt(real(kmin)),'mode 2','symb off'
      if (nplot.eq. 1) write(2,'(a,g12.4,a,g12.4,10a)')
     $'note (',fbar-fdel/2,' ',pmid*(1-ebar/3),
     $')',bks,'s',bks,bks,'sub{',bks,'0007',bks,'min,max}'
      if (ndet.eq.2) write(2,'(a,1p,2g14.7)')
     $'xlim 4.6 ',caylim(1),caylim(2)
      if (llog.gt. 0) write(2,'(a)') 'logxy 3'
c  Finish up
      write(2,'(a)') 'color 0','plot 1 4.6','stop'
      close (unit=2)
c
      call ljust(splot,splot, np)
      write(*,'(/2a)')'Plotxy file written to ',splot(1:np)
c
c
c  If 'hold' is on do not plot to screen.
      call getint('hold', ignor, kold)
      if (kold.ge. 0) return
 2000 if (nx.eq. 0) return
c  WARNING: see notes at end of subroutine autoc
      call system('sleep 1')
      call system('(plotxy < '//splot//' > /dev/null;mv mypost my$$;
     $gv my$$ ; sleep 5; rm my$$)&')
       return
       end
c______________________________________________________________________
      subroutine autoco
c$$$$ calls  fft getchr ljust system units
c  Estimates the sutocorrelation function from the PSD by taking
c  the Fourier cosine transform.
c
c  nf  number of frequencies in spectrum
c  sx() nf-vector of spectral estimates
c
      parameter (pi=3.141592653589793)
      parameter (mxx=500 000)
      parameter (nfmx=mxx + 2)
      character*64 rfile,title,unit,psfile*10
      common /series/ nx,dt,t1,x(mxx)
      common /state/ xbar,varx,smin,smax,ebar,kmin,kmax
      common /result/ nt,nf,fNyq,df,sx(nfmx),kopt(nfmx)
      common /pool/ fxre(8*mxx + 4)
      complex fx(4*mxx+2)
      dimension rx(mxx)
      equivalence (fxre(1),fx(1)),(fxre(mxx+1),rx(1))
      data rfile/'fort.15'/, nauto/0/
      save nauto
c
ce
      if (nx.le. 1) return
c
c  See if the autorrelation if needed, open a file
      call getchr('autocorrelation', rfile, needr)
      if (needr.lt. 0) return
      nauto=nauto + 1
      open (unit=15, file=rfile)
c
c  Adjust series length
      nf1=nf-1
c
c  Perform a complex FFT on sx; adjust 1st term for trapezoidal
c  quadrature; pad with zeros to get proper dt in autocovar fn
      fx(1)=sx(1)/2.0
      fx(nf1+1)=0.0
      do 1100 j=2, nf1
          fx(j)=sx(j)
          fx(j+nf1)=0.0
 1100 continue
c
c  Perform mixed-radix FFT
      nf2=2*nf1
      call fft(fxre, fxre(2), nf2,nf2,nf2, 2, ierr)
      if (ierr .ne. 0) then
        write(*,'(a,i7)')'>>> In autoco, FFT routine fails for n =',nf1,
     $    '    This error should never arise: autocorrelation not found'
        return
      endif
c
c  Normalize and copy into real array
      dtr=dt
      const=1.0/real(fx(1))
      do 1200 j=1, nf1
        rx(j)=const * real(fx(j))
 1200 continue
      amp=0.05*rx(1)
      do 1250 j=nf1, 1, -1
        if (abs(rx(j)).gt. amp) goto 1260
 1250 continue
 1260 limrx=min(nf1, nint(1.2*j+10))
c
      call getchr('title', title, nti)
c
c  Write plotfile to unit 15
      call units(i1,i2,i3, unit)
      if (nti.gt. 0) write(15,'(2a)') 'title ',title(1:nti)
      write (15,'(a/a/(2g12.4))')
     $'dash 0.02 0.04','read 2',0.0,0.0,limrx*dtr,0.0
      write(psfile,'(a,I2.2,a)')'auto',nauto,'.ps'
      write (15,'(a)') 'file '//rfile,'xlim 5.5', 'ylim 3.9',
     $'char 0.11','output '//psfile,'xlab Lag, '//unit(2:i1),
     $'ylab Autocorrelation','smooth','skip 26','dash 0 0'
      write (15,'(a,i7)') 'color ',3,'read',limrx
      write (15,'(a)') 'symb dot .14','file','skip 26','mode 2'
      write (15,'(a,i7)') 'color ',2,'read',limrx
      write (15,'(a)') 'color black','plot 1 6.8','stop',' '
      write (15,'(1p,2g13.6)') ((j-1)*dtr,rx(j),j=1, nf1)
      close (unit=15)
c
      call ljust(rfile,rfile, nr)
      write(*,'(/a,a)')
     $'Autocorrelation series and plot commands written to ',rfile(1:nr)
      if (limrx.lt. nf1) write(*,'(a,i7,a)')
     $'Only the first ',limrx,' terms are plotted'
c
c  WARNING: 'system' is a command not generally available in FORTRAN.
c  It launches a shell script that invokes plotxy and screen preview
c  command.  You may need to replace gv (the ghostview PostScript
c  viewer) with your favorite previewer.
      call system('(plotxy < '//rfile//' > /dev/null; sleep 1;
     $gv '//psfile//')& ')
c
      return
      end
c__________________________________________________________________
      subroutine adapt(ntimes, initap, nbound)
c$$$$ calls quick qorth getone kurb
c
c  Perform adaptive multitaper spectral estimation.
c  From pilot estimate of spectrum, computes estimates of S"
c  to be used in (13) of R & s (1995) for the Min Square Err
c  spectrum.
c
c  ntimes  number of iterations
c  initap    initial number of tapers
c  nbound    upper bound on number tapers applied
c  Returns in /result/:
c  nf      number of frequencies in PSD
c  sx()    nf-vector of PSD estimates
c  kopt()  integer nf-vector of number of tapers at each freq
c
      parameter (mxx=500 000)
      parameter (nfmx=mxx + 2)
      parameter (eps=1e-35)
      common /series/ nx,dt,t1,x(mxx)
      common /result/ nt,nf,fNyq,df,sx(nfmx),kopt(nfmx)
      parameter (my=mxx/2+1)
      dimension optk(my),y(my)
      equivalence (optk(1), sx(1)), (y(1), sx(my+1))
c
c  c2=480**0.2 = constant for parabolic weighting of
c  tapered spectra in quick; for uniform weighting c2=12.0**0.4=2.702
      data c2/3.437/
c
ce
c
c  Set the degree of smoothing by the spline
      kmax=1
      fact=1.0
      call getone('smooth', fact, ignor)
c
c  Get pilot estimate of spectrum
      initap = 3 + fact*sqrt(nx/25.0)
      call quick(1, -initap)
c
c
c  Adaptive iteration for MSE spectrum
      do 1600 iter=1, ntimes
c
        do 1100 j=1, nf
           y(j)= log(sx(j))
 1100   continue
c
c
c  Estimate K, number of tapers at each freq for MSE spectrum
c  R = S"/S -- use R = Y" + (Y')**2 , Y=ln S.
c  Note  c2=3.437 = 480*0.2
        do 1200 j=1, nf
c  Note smoothing for Y" is over twice nominal kernel width
          ispan = kopt(j)
          call qorth(nf, j-ispan, j+ispan, y, dy, ddy)
c
          optk(j)=c2/((eps + abs(ddy + dy**2)) /fact)**0.4
c
 1200   continue
c
c
c  Curb runaway growth of Kopt near zeros of R
        call kurb(nf, optk)
c
c  Always use 3 or more tapers; never average over more than length of psd
c  and apply preset bound
        bound = nf/3.0
        if (nbound .gt. 0) bound = min(bound, real(nbound))
        do 1550 j=1, nf
          kopt(j) = min(bound, max(3.0, optk(j)))
 1550   continue
c
c  Recompute spectrum with optimal variable taper numbers
        call quick(0, kmax)
 1600 continue
c
      return
      end
c_______________________________________________________________________
      subroutine qorth(n, i1, i2, s, ds, dds)
c$$$$ calls nothing
c  Performs LS fit to s on the interval (i1,i2)  to s by
c  a quadratic polynomial in an orthogonal basis.
c  Returns ds = estimate of 1st derivative  ds/dn  at center of record
c  Returns dds = estimate of 2nd derivative
c
      dimension s(*)
c
      L = i2 - i1 + 1
      el=L
      gamma = (el**2 - 1.0)/12.0
      u1sq = el*(el**2 - 1.0)/12.0
      u2sq = (el*(el**2 - 1.0)*(el**2- 4.0))/180.0
      amid= 0.5*(el + 1.0)
      dot1=0.0
      dot2=0.0
      do 1100 kk=1, L
        i=kk + i1 - 1
c  Negative or excessive index uses even function assumption
        if (i.le. 0) i=2 - i
        if (i.gt. n) i=2*n - i
        dot1 = dot1 + (kk - amid) * s(i)
        dot2 = dot2 + ((kk - amid)**2 - gamma)*s(i)
 1100 continue
      ds = dot1/u1sq
      dds = 2.0*dot2/u2sq
c
      return
      end
c_______________________________________________________________________
      subroutine kurb(nf, optk)
c$$$$ calls nothing
c  Curbs runaway growth of Kopt caused by zeros of S".  Limits
c  slopes to less than 1 in magnitude, this avoiding greedy
c  averaging
      dimension optk(*)
c  Scan forward and create new series where slopes <= 1
      istate=0;
      do 1500 j = 2,nf
        if (istate .eq. 0) then
           slope=optk(j)-optk(j-1);
           if (slope .ge. 1.0) then 
              istate=1;
              optk(j)=optk(j-1)+1.0;
           endif
        else
           if (optk(j) .ge. optk(j-1)+1.0) then
              optk(j) = optk(j-1)+1.0;
           else
              istate=0;
           endif
        endif 
 1500 continue
c
c  Scan backward to bound slopes >= -1
      istate=0;
      do 1600 j = nf, 2, -1
        if (istate .eq. 0) then
           slope=optk(j-1)-optk(j);
           if (slope .ge. 1.0) then 
              istate=1;
              optk(j-1)=optk(j)+1.0;
           endif
        else
           if (optk(j-1) .ge. optk(j)+1.0) then
              optk(j-1) = optk(j)+1.0;
           else
              istate=0;
           endif
        endif
 1600 continue
      return
      end
c________________________________________________________________________
      subroutine quick(new, ktop)
c$$$$ calls fft nearn
c  Sine multitaper routine.  Performs only 1 FFT then
c  constructs FT[sin(q*n)*x(n)] from FT[x(n)].
c
c  if new > 0 finds FFT of data series and saves it in /pool/
c       otherwise, uses uses previously saved transform
c  ktop if > 0, largest number in array kopt
c       if < 0, -ktop = constant value to be used; insert
c       -ktop into the array kopt.
c  kopt() nf-vector giving number of tapers to be averaged
c    into each estimate a frequency number j, if ktop > 0.
c  Works on:
c  x() n-vector of the time series in /series/
c  Returns into /results/:
c  nf  number of frequencies in spectrum
c  sx() nf-vector of spectral estimates
c  if ktop < 0 fills kopt with -ktop
c
c  The spectrum is based on parabolic weighting of the individual
c  estimates.
      parameter (pi=3.141592653589793)
      parameter (mxx=500 000)
      parameter (nfmx=mxx + 2)
      complex fx,zz
      common /series/ nx,dt,t1,x(mxx)
      common /result/ nt,nf,fNyq,df,sx(nfmx),kopt(nfmx)
      common /pool/ fxre(8*mxx + 4)
      dimension fx(4*mxx+2)
      equivalence (fxre(1),fx(1))
c
c
c  Adjust series length if necessary to be 2**i * 3**j * 5**k
c  with i > 0, j,k >= 0.  Find other constants
ce
      npower=2*nearn(nx/2)
      nptwo=2*npower
      nf=1 + npower/2
      nt=2*nf - 2
      fNyq=0.5/dt
      df=fNyq/(nf - 1)
c
c  Perform a double-length complex FFT on data and save it.
c  Done only once per time series!
      if (new.ge. 1) then
        do 1100 j=1, npower
          fx(j)=x(j)
 1100   continue
        do 1150 j=npower+1, nptwo
          fx(j)=0.0
 1150   continue
c
c  Perform mixed-radix FFT
        call fft(fxre, fxre(2), nptwo,nptwo,nptwo, 2, ierr)
        if (ierr .ne. 0) then
          write(*,'(a,i7)')'>>> FFT routine fails for n =',nptwo,
     $    '    This error should never arise'
          stop
         endif
      endif
c
c  Fill kopt vector if necessary
      if (ktop.lt. 0) then
        do 1250 j=1, nf
          kopt(j)=-ktop
 1250   continue
      endif
c
c
c  Loop over frequency
      do 1500 m=0, nf-1
        m2=2*m
        sx(m+1)=0.0
        klim=kopt(m+1)
        ck=1.0/klim**2
c  Average over tapers, parabolic weighting wk
        do 1300 k=1, klim
            j1=mod(m2+nptwo-k, nptwo)
            j2=mod(m2+k, nptwo)
            zz=fx(j1+1) - fx(j2+1)
            wk=1.0 - ck*(k-1)**2
            sx(m+1)=sx(m+1) + (real(zz)**2 + aimag(zz)**2)*wk
 1300   continue
*       sx(m+1)=sx(m+1)/klim
c  Exact normalization for parabolic factor
        sx(m+1)=sx(m+1) * (6.0*klim)/(4*klim**2+3*klim-1)
 1500 continue
c
      return
      end
c__________________________________________________________________
      function nearn(n)
c$$$$ calls nothing
c  Finds the largest integer less than or equal to n whose
c  prime factors are in the set {2, 3, 5}.
c
      parameter (c=0.693147181,b=1.098612289,a=1.609437912)
ce
      p=log(real(n)+0.5)
      ka=p/a
      kb=p/b
      kc=p/c
c
      dmin=p
      do 1500 ja=0, ka
        do 1400 jb=0, kb
          do 1300 jc=0, kc
            d=p -ja*a -jb*b -jc*c
            if (d.lt. 0.0) goto 1400
            if (d.lt. dmin) dmin=d
 1300     continue
 1400   continue
 1500 continue
c
      nearn=nint(exp(p-dmin))
      return
      end
c_______________________________________________________________
c=======================================================================
c                Unit D:  Get the Data
c=======================================================================
      subroutine getdat
c$$$$ calls getarr getchr getone getint spline eval evlak akima ljust
c  Gets the file name and reads data into vector x().
c  Also spline interpolates onto equal interval if so requested; natural
c  or Akima splines available.
c  Decimates input file if required.
c  Issues error messages.
      character*64 name, even, spln
      parameter (mxx=500 000,  pi=3.14159265)
      dimension tab(100),column(2), decn(2),terms(2)
      dimension t(mxx),u(mxx),ddu(mxx),wrk(mxx)
      common /series/ nx,dt,t1,x(mxx)
      common /pool/ work(8*mxx+4)
      equivalence (t(1),work(1)),(u(1),work(mxx+1)),
     $ (ddu(1),work(2*mxx+1)),(wrk(1),work(3*mxx+1))
      data name(1:1)/'  '/,nterm/mxx/, column/1,2/, decn/0,0/
      data spln(1:3)/'nat'/, even/'fort.14'/, dtnew/1/, kon/0/
      data nfilt/1/, ieof/0/
      save nterm, dtnew, kon
ce
      call getchr('file', name, nbyte)
      if (nbyte.eq. 0) goto 3200
c
c  Continue reading from the same file if name='-continue'
c  otherwise close previous unit, reopen new file
      if (name(1:4) .ne. '-con')  then
        close(unit=11)
        open (unit=11, file=name, status='OLD', err=3100)
        kon=1
      elseif (kon.eq. 0) then
        write(*,'(a)')
     $ '>>> No data file currently open; PSD cannot continue'
        stop
      endif
c
c  Skip records before reading
      iskip=0
      call getint('skip', iskip, none)
      call getarr('terms',terms,2, ntrms)
      if (ntrms.eq. 2) then
        iskip=nint(terms(1))-1
        rewind(unit=11)
      endif
      if (iskip.gt. 0) then
        do 1010 j=1,  iskip
          read (11,'(a1)', end=3600) dum
 1010   continue
        write(*,'(a,i7)')'   Skipped records before data:',iskip
      endif
c
c  Decimate the series?  Compute antialias filter weights u
      nd=1
      call getarr('decimate', decn, 2, idecid)
      if (decn(1).lt. 2.0) idecid=0
      if (idecid.le. 0) goto 1100
      nd=decn(1)
      nfilt=max(10*nd*(2-idecid) + decn(2)*(idecid-1), nd+1.0)
      usum=0.0
      do 1020 j=1, nfilt
        thet=pi*(0.5*nfilt+0.5-j) + 1e-6
        u(j)=sin(thet/decn(1)) * cos(thet/(nfilt+1.0)) / thet
        usum=usum + u(j)
 1020 continue
      do 1021 j=1, nfilt
        u(j)=u(j)/usum
 1021 continue
c
c
c  Get number of terms to be read
 1100 nt=mxx
      call getint('nterms', nterm, nte)
      call getarr('terms', terms, 2, ntrms)
      if (ntrms.eq. 2) nterm=nint(terms(2)) - iskip
      if (ntrms.eq. 1) nterm=nint(terms(1))
      if (nte.gt. 0 .or. ntrms.ge. 1) nt=min(mxx, nterm)
c
c  Get the sampling interval
      call getone('interval', dtnew, newdt)
      dt=dtnew * nd
      t1=0.0
      if (nd.gt. 1)write(*,'(a,g12.4)')
     $' Sampling interval before decimation: ',dt/nd
c
c  Check for spline parameters
      call getchr('spline', spln, isplin)
c
      if ((newdt.le.0 .or. dt.eq.0) .and. isplin.ge.0) goto 3700
      if (isplin.gt.0 .and. idecid.gt.1) write(*,'(a)')
     $'>>> Spline and decimation incompatible: PSD applies spline only'
c
c  Get column(s) to be read
      call getarr('column', column, 2, kols)
c
      if (isplin.lt. 0) then
c  Interpolation not requested
        if (kols.eq. 2) write(*,'(2a)') '>>> Warning: column takes ',
     $ '1 number, unless spline is specified - First number accepted'
c
c  Read from a single column; no decimation
        koln=column(1)
        if (koln .ge. 1 .and. idecid.le. 0) then
          do 1600 j=1, nt
             read (11,*,err=3000,end=1650) (tab(l),l=1,koln)
             x(j)=tab(koln)
 1600     continue
          j=nt + 1
 1650     nx=j - 1
          nr=nx
          ieof=nt - j
c
c  Decimate series
        elseif (koln .ge. 1 .and. idecid.gt.0) then
c
c  If 'terms' command is used n2 is term number in the file so:
          if (ntrms.eq. 2) nt=nt/nd
          do 1710 k=1, nfilt-nd
            read(11,*,err=3000,end=1760) (tab(l),l=1,koln)
            t(k)=tab(koln)
 1710     continue
          do 1750 j=1, nt
            do 1715 k=1, nd
              read(11,*,err=3000,end=1760) (tab(l),l=1,koln)
              t(nfilt+k-nd)=tab(koln)
 1715       continue
            x(j)=0.0
            do 1725 k=1, nfilt
              x(j)=x(j) + t(k)*u(k)
 1725       continue
            do 1735 k=1, nfilt-nd
              t(k)=t(k+nd)
 1735       continue
 1750     continue
          j=nt+1
          k=1
 1760     nx=j-1
          nr=(j-1)*nd + k - 1
          ieof=nd*nt - nr
c
c  Read every number on each line (column 0)
        elseif (column(1).eq.0) then
          if (idecid.gt. 0) write(*,'(a)')
     $ 'Decimation allowed only when column >= 1'
c  Fill x with a junk value
          do 1800 j=1, nt
            x(j)=-0.4342944e+15
 1800     continue
          read (11, *, err=3000, end=1810) (x(j),j=1, nt)
 1810     do 1850 j=1, nt
            if (x(j).eq. -0.4342944e+15) goto 1855
 1850     continue
          j=nt + 1
 1855     nx=j - 1
          nr=nx
          ieof=nt - j
        endif
c
        if (ieof.gt. 0) then
          write(*,'(a)')'EOF detected in data file'
          rewind(unit=11)
        endif
        write(*,'(a,i7)')'Number of terms read from file:',nr
        if (nx.eq. mxx) write(*,'(2a,i7,a)')'>>> Array space filled:',
     $  ' series truncated to',mxx,' terms'
c
c
c  Spline interpolation requested
      else
        nt=min(nt, mxx)
        if (kols.lt. 2) then
           kol1=1
           kol2=2
           write(*,'(2a)')'>>> Columns undefined for spline: ',
     $    ' reading from cols 1 & 2'
        else
           kol1=column(1)
           kol2=column(2)
        endif
        koln=max(kol1, kol2)
        list=0
        gape=0.0
        gap=gape
        do 2050 j=1, nt
          read (11,*,err=3000, end=2100) (tab(l),l=1,koln)
          t(j)=tab(kol1)
          u(j)=tab(kol2)
c  Check for improper time and large gaps
          if (j.ge. 2 ) then
            if (t(j).le. t(j-1)) then
              list=list + 1
              wrk(list)=j
            endif
            gap=max(gap, t(j)-t(j-1))
            if (gap.gt. gape) jump=j
            gape=gap
          endif
 2050   continue
        if (nt.eq. mxx) write(*,'(a,i7)')
     $  '>>> Array space filled: series truncated at ',mxx
 2100   nx=j - 1
        if (list.gt. 0) goto 3500
        write(*,'(a,i7)')'  Length of input time series:',nx
        write(*,'(2(a,g12.4))')
     $  ' Time variable runs from',t(1),' to ',t(nx)
        if (gap/dt.ge. 10.0) write(*,'(/a,i7/a,f9.1)')
     $  '>>> Warning: A data gap of > 10*dt was encountered at n=',
     $  jump,'    gap/dt=',gap/dt
c
        if (spln(1:3).eq. 'nat') then
          call spline(nx, t, u, ddu, wrk)
        else
          call akima(nx, t, u, wrk)
        endif
c
        t1=t(1)
        do 2200 j=1, mxx
          tx=t1 + (j-1)*dt
          if (tx.gt. t(nx)) goto 2300
          if (spln(1:3).eq. 'nat') x(j)=eval(tx, nx, t, u, ddu)
          if (spln(1:3).ne. 'nat') x(j)=evlak(tx, nx, t, u, wrk)
 2200   continue
        j=mxx+1
        write(*,'(a,g12.4)')'>>> Interpolated series truncated at t =',
     $  tx
 2300   nx=j - 1
        write(*,'(a,i7)')'   Length after interpolation:',nx
c
      endif
c
c  Scale data
      scale=0.0
      call getone('scale', scale, ignor)
      if (scale.ne. 0.0) then
        do 2500 j=1, nx
          x(j)=scale*x(j)
 2500   continue
      endif
c
c  If requested write interpolated data to a diskfile
      call getchr('validate', even, nval)
      if (nval.ge. 0 .and. (idecid.gt. 0 .or. isplin.gt. 0)) then
        if (nval.gt. 0) open(unit=14, file=even)
        t1=t(1)
        if (idecid.gt. 0) t1=0.5*dt*(nfilt+1.0)/nd
        write(14,'(2g13.5)') (t1+(j-1)*dt,x(j),j=1,nx)
        call ljust(even,even, nev)
        write(*,'(/a,a)')
     $ 'Interpolated/decimated series written to: ',even(1:nev)
      endif
c
      return
c
c  Error out
 3000 write(*,'(a,i7)')'>>> Unreadable number in data file at point',j
      stop
 3100 write(*,'(2a)')'>>> Unable to open the file: ',name
      stop
 3200 write(*,'(a)')'>>> PSD cannot continue without a file name'
      stop
 3500 write(*,'(a,/a/(10i6))')'>>> Time does not increase uniformly',
     $ '    Bad points at these term numbers:',
     $(int(wrk(j)),j=1, min(list,100))
      if (list.gt. 100) write(*,'(a)')'List terminated after 100 terms'
      stop
 3600 write(*,'(a)')'>>> End-of-file encountered while skipping'
      stop
 3700 write(*,'(a)')'>>> You must set the interval dt to use spline'
      stop
      end
c_______________________________________________________________________
      subroutine stats
c$$$$ calls getone units
c  Finds some statistics of the series - removes mean and trend
c  Prints results in the proper units
      parameter (mxx=500 000)
      common /state/ xbar,varx,smin,smax,ebar,kmin,kmax
      common /series/ nx,dt,t1,x(mxx)
      double precision sum,xmean
      character*64 unit
      data i1,i2,i3/1,1,1/
ce
      if (nx.le. 1) return
      sum=0.0
c  Find mean
      xup=x(1)
      xlo=xup
      do 1100 j=1, nx
        sum=x(j) + sum
        xup=max(x(j), xup)
        xlo=min(x(j), xlo)
 1100 continue
      xmean=sum/nx
c  Find variance
      sum=0.0
      do 1500 j=1, nx
        sum=sum + (x(j)-xmean)**2
 1500 continue
      varx=sum/nx
c
c  Removes mean from the series
      do 1600 j=1, nx
        x(j)=x(j) - xmean
 1600 continue
      xbar=xmean
c
c  Estimates slope of LS line
      tx=0.0
      tbar=0.5+ 0.5*nx
      do 1700 j=1, nx
        tx=tx + (j-tbar)*x(j)
 1700 continue
      xn=nx
      slope=12.0*tx/(xn*(xn**2+6.0*xn-4.0))
      write(*,'(/a,g12.4)')
     $'      Offset produced by trend:  ',slope*nx
c  Remove the slope if requested
      call getone('detrend', dummy, mark)
      if (mark.ge. 0) then
        sum=0.0
        do 1800 j=1, nx
           x(j)=x(j) - (j-tbar)*slope
           sum=sum + x(j)**2
 1800   continue
        varx=sum/nx
        write(*,'(6x,a)')'The series has been detrended'
      endif
      sig=sqrt(varx)
c
c  Get units for the statistical values
      call units(i1,i2,i3, unit)
c
      call getone('decimate', dec, idecid)
      if (idecid.gt. 0) write(*,'(/a,i4)')
     $'Statistics after decimation by',nint(dec)
      write(*,'(/(a,1p,g12.4,1x,a))')
     $'          Mean value of series:',xbar,unit(i1:i2),
     $'        Variance of the series:',varx,unit(i2:i3),
     $'                    RMS signal:',sig, unit(i1:i2),
     $'   Sampling interval of series:',dt,' '//unit(2:i1),
     $'            Duration of series:',dt*(nx-1),' '//unit(2:i1)
c
      write(*,'(/a)')'The mean value has been removed from the series'
c
      if (max(xup-xbar, xbar-xlo).gt. 10.0*sig) write(*,'(a)')' ',
     $'>>> Warning: Series contains outliers > 10 sigma in magnitude'
c
      return
      end
c_______________________________________________________________________
      subroutine units(i1,i2,i3, unit)
c$$$$ calls getchr
c  Reads in the names of the units of time variable and signal variable.
c  On return unit(2:i1) holds time variable unit + a space
c            unit(i1+1,i2) holds signal variable unit name + space
c            unit(i2+1,i3) holds signal-squared unit name + space
c            unit(1:1)='/'
      character*64 unit
ce
      call getchr('units', unit(2:64), inut)
      if (inut.le. 0) then
        i1=1
        i2=1
        i3=1
        unit=' '
      elseif (inut.ge. 1) then
        unit(1:1)='/'
        i1=index(unit, ' ')
        i2=index(unit(i1+1:64), ' ')+i1
        i3=2*i2 - i1 + 2
        if (i2.gt.i1+1) unit(i2+1:i3)=unit(i1+1:i2-1)//'^2'
      endif
      return
      end
c_______________________________________________________________________
      subroutine yule(iphase, nord)
c$$$$ calls toepl slowft
c  Prewhitening of the spectrum.  See Percival and Walden,
c  "Spectral Analaysis ...", Chap 9, 1993.
c
c  nord = number of FIR coeffs to be found.
c  iphase=1, uses spectral estimates  sx()  to calculate
c    autocovariances ar() for Yule-Walker equations.  Solves
c    Y-W equations then prewhitens time series x().
c  iphase=2, modifies sx() - undoes prewhitening by dividing
c    out the squared FT of the FIR filter.
c  iphase=3, skip generation and use prior filter already in memory
      parameter (pi=3.141592653589793)
      parameter (mxx=500 000)
      parameter (nfmx=mxx + 2)
      common /series/ nx,dt,t1,x(mxx)
      common /result/ nt,nf,fNyq,df,sx(nfmx),kopt(nfmx)
      parameter (ndmx=20)
      dimension row(2*ndmx),rhs(ndmx),ar(ndmx),wrk(2*ndmx)
      save ar,nrd
      data nrd/0/
ce
      goto (1000, 2000, 1200), iphase
 1000 nrd=min(ndmx-1, nord)
      if (nrd .ne. nord) write(*,'(a,i3)')
     $'>>> Order of prewhitening filter reduced to ',nrd
c
c  Phase 1: prewhiten time series.
c
      nff=nf-1
c  Fourier transform PSD estimate to find autocovariance function.
c  Load values into Toeplitz matrix and rhs vector
      row(nrd+2)=0.0
      do 1100 i=nrd+1, 1 , -1
         omi=pi*(i-1)/nff
         call slowft(nff, sx, omi, row(i), st)
         row(nrd+i)=row(i)
         rhs(i)=row(i+1)
 1100 continue
c
c  Solve symmetric Toeplitz system for filter coeffs ar()
      ar(1)=-1.0
      call toepl(nrd, row, rhs, wrk, ar(2), iok)
c
c  Convolve filter ar() into series x(), overwriting
 1200 if (nrd.gt. 0) then
        nord=nrd
        do 1500 j=nx, 1, -1
          sum=x(j)
          do 1400 jf=2, nrd+1
            jx=max(1, j-jf+1)
            sum=sum - ar(jf)*x(jx)
 1400     continue
          x(j)=sum
 1500   continue
c
        write(*,'(a/(1p,5g12.5))') 'Prewhitening FIR filter:',
     $ (ar(j),j=1,nrd+1)
      else
        write(*,'(a)')'>>>> No prior prewhitening filter'
      endif
      return
c
c  Phase 2: Correct spectral estimates for prewhitening.
c  Uses ar() found in phase 1.
 2000 if (nrd.eq. 0) return
      dom=pi/nf
c  Divide out from sx() the FIR spectrum of filter coeffs ar()
      do 2100 k=1, nf
        call slowft(nrd+1, ar, (k-1)*dom, ct,st)
        ftsq= st**2 + ct**2
        sx(k)=sx(k) / ftsq
 2100 continue
c
      return
      end
c________________________________________________________________
      subroutine toepl(n, r, y, ff, x, iok)
c$$$$ calls nothing
c  Solves linear system in Toeplitz form.
c  n    order of matrix.
c  r(1) ... r(n)     1st row of Toeplitz array.  Overwritten.
c  r(n+1) ... r(2*n) 1st col of Toeplitz array.  Overwritten.
c  y()    the rhs vector.
c  ff()   working array of size at least  2*n.
c  x()    the solution vector.
c  iok  1 if everything ok, 0 if singular principal minor.
c
      dimension r(*),x(*),y(*),ff(2,*)
c
c  Re-order  r  into crazy pattern:
c  r(1),r(2), ... r(n):       the first row of the matrix BACKWARDS!
c  r(n),r(n+1), ... r(2*n-1): the first column of the matrix.
ce
      do 1100 i=1, n/2
        ri=r(i)
        r(i)=r(n-i+1)
        r(n-i+1)=ri
 1100 continue
      do 1150 i=n+1, 2*n
        r(i-1)=r(i)
 1150 continue
c
      iok=1
      if(r(n).eq.0.) goto 99
      x(1)=y(1)/r(n)
      if(n.eq.1)return
      ff(1,1)=r(n-1)/r(n)
      ff(2,1)=r(n+1)/r(n)
      do 15 m=1,n
        m1=m+1
        sxn=-y(m1)
        sd=-r(n)
        do 11 j=1,m
          sxn=sxn+r(n+m1-j)*x(j)
          sd=sd+r(n+m1-j)*ff(1,m-j+1)
11      continue
        if(sd.eq.0.)goto 99
        x(m1)=sxn/sd
        do 12 j=1,m
          x(j)=x(j)-x(m1)*ff(1,m-j+1)
12      continue
        if(m1.eq.n) return
        sgn=-r(n-m1)
        shn=-r(n+m1)
        sgd=-r(n)
        do 13 j=1,m
          sgn=sgn+r(n+j-m1)*ff(1,j)
          shn=shn+r(n+m1-j)*ff(2,j)
          sgd=sgd+r(n+j-m1)*ff(2,m-j+1)
13      continue
        if(sd.eq.0..or.sgd.eq.0.)goto 99
        ff(1,m1)=sgn/sgd
        ff(2,m1)=shn/sd
        k=m
        m2=(m+1)/2
        pp=ff(1,m1)
        qq=ff(2,m1)
        do 14 j=1,m2
          pt1=ff(1,j)
          pt2=ff(1,k)
          qt1=ff(2,j)
          qt2=ff(2,k)
          ff(1,j)=pt1-pp*qt2
          ff(1,k)=pt2-pp*qt1
          ff(2,j)=qt1-qq*pt2
          ff(2,k)=qt2-qq*pt1
          k=k-1
14      continue
15    continue
c
99    iok=0
      return
      end
c_______________________________________________________________________
      subroutine slowft(n, x, om, ct, st)
c$$$$ calls nothing
c
c  Calculates Fourier transform of real sequence x(i),i=1,...n
c  at angular frequency om normalized so that nyquist=pi:
c     cmplx(ct, st) = sum {k=0 to n-1} x(k+1)*exp(i*k*om)
c  The sine transform is returned in st and the cosine transform in
c  ct.  Algorithm is that of Goertzal with modifications by Gentleman,
c  Comp.J. 1969.  The transform is not normalized.  To normalize
c  one-sided ft, divide by sqrt(data length).  For positive om, the ft
c  is defined as ct-(0.,1.)st or like slatec cfftf
c
c  Slightly modified from Slatec library
c
      implicit double precision (a-h,o-z)
      real x, om, ct, st
      parameter (pi=3.141592653589793238d0,tp=2.d0*pi)
      dimension x(n)
ce
      omega=om
      np1=n+1
      l=6.d0*omega/tp
      s=sin(omega)
      a=0.d0
      c=0.d0
      d=0.d0
      e=0.d0
      if(l.eq.0)then
c  Recursion for low frequencies (.lt. nyq/3)
        b=-4.d0*sin(omega/2.d0)**2
        do 10 k=1,n
          c=a
          d=e
          a=x(np1-k)+b*d+c
   10     e=a+d
      elseif(l.eq.1)then
c  Regular Goertzal algorithm for intermediate frequencies
        b=2.d0*cos(omega)
        do 20 k=1,n
          a=x(np1-k)+b*e-d
          d=e
   20     e=a
      else
c  Recursion for high frequencies (.gt. 2*nyq/3)
        b=4.d0*cos(omega/2.d0)**2
        do 30 k=1,n
          c=a
          d=e
          a=x(np1-k)+b*d-c
   30     e=a-d
      endif
      st=-s*d
      ct=a-b*d/2.d0
      return
      end
c_______________________________________________________________________
      subroutine fft(a, b, ntot, n, nspan, isn, ierr)
c$$$$ calls nothing
c  Multivariate complex Fourier transform, computed in place using
c  mixed-radix fast Fourier transform algorithm.  By R. C. Singleton,
c  Stanford Research Institute, Sept. 1968.  Arrays a and b originally
c  hold the real and imaginary components of the data, and return the
c  real and imaginary components of the resulting Fourier
c  coefficients.  Multivariate data is indexed according to the Fortran
c  array element successor function, without limit on the number of
c  implied multiple subscripts.  The subroutine is called once for each
c  variate.  The calls for a multivariate transform may be in any
c  order.
c   ntot is the total number of complex data values.
c   n is the dimension of the current variable.
c   nspan/n is the spacing of consecutive data values while indexing
c   the current variable.
c  The integer ierr is an error return indicator. It is normally zero,
c  but is set to 1 if the number of terms cannot be factored in the
c  space available. If it is permissible the appropriate action at this
c  stage is to  enter fft again after having reduced the length of the
c  series by one term.  The sign of isn determines the sign of the
c  complex exponential, and the magnitude of isn is normally one.  A
c  tri-variate transform with a(n1,n2,n3), b(n1,n2,n3) is computed by
c      call fft(a, b, n1*n2*n3, n1, n1, 1)
c      call fft(a, b, n1*n2*n3, n2, n1*n2, 1)
c      call fft(a, b, n1*n2*n3, n3, n1*n2*n3, 1)
c   For a single-variate
c   transform, ntot = n = nspan = (number of complex data values),
c   e.g.    call fft(a, b, n, n, n, 1)
c  With most Fortran compilers the data
c  can alternatively be stored in a single complex array a, then the
c  magnitude of isn changed to 2 to give the correct indexing
c  increment and a(2) used to pass the initial address for the sequence
c  of imaginary values, e.g.
c      call fft(a, a(2), ntot, n, nspan, 2)
c  Arrays at(maxf), ck(maxf),
c   bt(maxf), sk(maxf), and np(maxp) are used for temporary storage.
c  if the available storage is insufficient, the program is terminated
c  by the error return option
c   maxf must be .ge. the maximum prime factor of n.  maxp must be .gt.
c   the number of prime factors of n.  In addition, if the square-free
c  portion k of n has two or more prime factors, then maxp must be .ge.
c  k-1.
c  Array storage in nfac for a maximum of 15 prime factors of n.
c  If n has more than one square-free factor, the product of the
c    square-free factors must be .le. 210
      parameter (mmaxf=23, maxp=209)
c  Array storage for maximum prime factor of 23
      dimension at(mmaxf),ck(mmaxf),bt(mmaxf),sk(mmaxf)
      dimension nfac(11),np(maxp)
      dimension a(*),b(*)
      equivalence (i,ii)
      data k3/0/, c2,c3,s2,s3/4*0/
c
      maxf=mmaxf
      ierr=0
      if(n.lt. 2) return
      inc=isn
      c72=0.30901699437494742
      s72=0.95105651629515357
      s120=0.86602540378443865
      rad=6.2831853071796
      if(isn.ge. 0) goto 10
      s72=-s72
      s120=-s120
      rad=-rad
      inc=-inc
   10 nt=inc*ntot
      ks=inc*nspan
      kspan=ks
      nn=nt-inc
      jc=ks/n
      radf=rad*float(jc)*0.5
      i=0
      jf=0
c  determine the factors of n
      m=0
      k=n
      goto 20
   15 m=m+1
      nfac(m)=4
      k=k/16
   20 if(k-(k/16)*16 .eq. 0) goto 15
      j=3
      jj=9
      goto 30
   25 m=m+1
      nfac(m)=j
      k=k/jj
   30 if(mod(k,jj).eq. 0) goto 25
      j=j+2
      jj=j**2
      if(jj .le. k) goto 30
      if(k.gt. 4) goto 40
      kt=m
      nfac(m+1)=k
      if(k .ne. 1) m=m+1
      goto 80
   40 if(k-(k/4)*4 .ne. 0) goto 50
      m=m+1
      nfac(m)=2
      k=k/4
   50 kt=m
      j=2
   60 if(mod(k,j) .ne. 0) goto 70
      m=m+1
      nfac(m)=j
      k=k/j
   70 j=((j+1)/2)*2+1
      if(j .le. k) goto 60
   80 if(kt .eq. 0) goto 100
      j=kt
   90 m=m+1
      nfac(m)=nfac(j)
      j=j-1
      if(j .ne. 0) goto 90
c  compute Fourier transform
  100 sd=radf/float(kspan)
      cd=2.0*sin(sd)**2
      sd=sin(sd+sd)
      kk=1
      i=i+1
      if(nfac(i) .ne. 2) goto 400
c  transform for factor of 2 (including rotation factor)
      kspan=kspan/2
      k1=kspan+2
  210 k2=kk+kspan
      ak=a(k2)
      bk=b(k2)
      a(k2)=a(kk)-ak
      b(k2)=b(kk)-bk
      a(kk)=a(kk)+ak
      b(kk)=b(kk)+bk
      kk=k2+kspan
      if(kk .le. nn) goto 210
      kk=kk-nn
      if(kk .le. jc) goto 210
      if(kk .gt. kspan) goto 800
  220 c1=1.0-cd
      s1=sd
  230 k2=kk+kspan
      ak=a(kk)-a(k2)
      bk=b(kk)-b(k2)
      a(kk)=a(kk)+a(k2)
      b(kk)=b(kk)+b(k2)
      a(k2)=c1*ak-s1*bk
      b(k2)=s1*ak+c1*bk
      kk=k2+kspan
      if(kk .lt. nt) goto 230
      k2=kk-nt
      c1=-c1
      kk=k1-k2
      if(kk .gt. k2) goto 230
      ak=cd*c1+sd*s1
      s1=(sd*c1-cd*s1)+s1
      c1=c1-ak
      kk=kk+jc
      if(kk .lt. k2) goto 230
      k1=k1+inc+inc
      kk=(k1-kspan)/2+jc
      if(kk .le. jc+jc) goto 220
      goto 100
c  transform for factor of 3 (optional code)
  320 k1=kk+kspan
      k2=k1+kspan
      ak=a(kk)
      bk=b(kk)
      aj=a(k1)+a(k2)
      bj=b(k1)+b(k2)
      a(kk)=ak+aj
      b(kk)=bk+bj
      ak=-0.5*aj+ak
      bk=-0.5*bj+bk
      aj=(a(k1)-a(k2))*s120
      bj=(b(k1)-b(k2))*s120
      a(k1)=ak-bj
      b(k1)=bk+aj
      a(k2)=ak+bj
      b(k2)=bk-aj
      kk=k2+kspan
      if(kk .lt. nn) goto 320
      kk=kk-nn
      if(kk .le. kspan) goto 320
      goto 700
c  transform for factor of 4
  400 if(nfac(i) .ne. 4) goto 600
      kspnn=kspan
      kspan=kspan/4
  410 c1=1.0
      s1=0
  420 k1=kk+kspan
      k2=k1+kspan
      k3=k2+kspan
      akp=a(kk)+a(k2)
      akm=a(kk)-a(k2)
      ajp=a(k1)+a(k3)
      ajm=a(k1)-a(k3)
      a(kk)=akp+ajp
      ajp=akp-ajp
      bkp=b(kk)+b(k2)
      bkm=b(kk)-b(k2)
      bjp=b(k1)+b(k3)
      bjm=b(k1)-b(k3)
      b(kk)=bkp+bjp
      bjp=bkp-bjp
      if(isn .lt. 0) goto 450
      akp=akm-bjm
      akm=akm+bjm
      bkp=bkm+ajm
      bkm=bkm-ajm
      if(s1 .eq. 0) goto 460
  430 a(k1)=akp*c1-bkp*s1
      b(k1)=akp*s1+bkp*c1
      a(k2)=ajp*c2-bjp*s2
      b(k2)=ajp*s2+bjp*c2
      a(k3)=akm*c3-bkm*s3
      b(k3)=akm*s3+bkm*c3
      kk=k3+kspan
      if(kk .le. nt) goto 420
  440 c2=cd*c1+sd*s1
      s1=(sd*c1-cd*s1)+s1
      c1=c1-c2
      c2=c1**2-s1**2
      s2=2.0*c1*s1
      c3=c2*c1-s2*s1
      s3=c2*s1+s2*c1
      kk=kk-nt+jc
      if(kk .le. kspan) goto 420
      kk=kk-kspan+inc
      if(kk .le. jc) goto 410
      if(kspan .eq. jc) goto 800
      goto 100
  450 akp=akm+bjm
      akm=akm-bjm
      bkp=bkm-ajm
      bkm=bkm+ajm
      if(s1 .ne. 0) goto 430
  460 a(k1)=akp
      b(k1)=bkp
      a(k2)=ajp
      b(k2)=bjp
      a(k3)=akm
      b(k3)=bkm
      kk=k3+kspan
      if(kk .le. nt) goto 420
      goto 440
c  transform for factor of 5 (optional code)
  510 c2=c72**2-s72**2
      s2=2.0*c72*s72
  520 k1=kk+kspan
      k2=k1+kspan
      k3=k2+kspan
      k4=k3+kspan
      akp=a(k1)+a(k4)
      akm=a(k1)-a(k4)
      bkp=b(k1)+b(k4)
      bkm=b(k1)-b(k4)
      ajp=a(k2)+a(k3)
      ajm=a(k2)-a(k3)
      bjp=b(k2)+b(k3)
      bjm=b(k2)-b(k3)
      aa=a(kk)
      bb=b(kk)
      a(kk)=aa+akp+ajp
      b(kk)=bb+bkp+bjp
      ak=akp*c72+ajp*c2+aa
      bk=bkp*c72+bjp*c2+bb
      aj=akm*s72+ajm*s2
      bj=bkm*s72+bjm*s2
      a(k1)=ak-bj
      a(k4)=ak+bj
      b(k1)=bk+aj
      b(k4)=bk-aj
      ak=akp*c2+ajp*c72+aa
      bk=bkp*c2+bjp*c72+bb
      aj=akm*s2-ajm*s72
      bj=bkm*s2-bjm*s72
      a(k2)=ak-bj
      a(k3)=ak+bj
      b(k2)=bk+aj
      b(k3)=bk-aj
      kk=k4+kspan
      if(kk .lt. nn) goto 520
      kk=kk-nn
      if(kk .le. kspan) goto 520
      goto 700
c  transform for odd factors
  600 k=nfac(i)
      kspnn=kspan
      kspan=kspan/k
      if(k .eq. 3) goto 320
      if(k .eq. 5) goto 510
      if(k .eq. jf) goto 640
      jf=k
      s1=rad/float(k)
      c1=cos(s1)
      s1=sin(s1)
      if(jf .gt. maxf) goto 998
      ck(jf)=1.0
      sk(jf)=0.0
      j=1
  630 ck(j)=ck(k)*c1+sk(k)*s1
      sk(j)=ck(k)*s1-sk(k)*c1
      k=k-1
      ck(k)=ck(j)
      sk(k)=-sk(j)
      j=j+1
      if(j .lt. k) goto 630
  640 k1=kk
      k2=kk+kspnn
      aa=a(kk)
      bb=b(kk)
      ak=aa
      bk=bb
      j=1
      k1=k1+kspan
  650 k2=k2-kspan
      j=j+1
      at(j)=a(k1)+a(k2)
      ak=at(j)+ak
      bt(j)=b(k1)+b(k2)
      bk=bt(j)+bk
      j=j+1
      at(j)=a(k1)-a(k2)
      bt(j)=b(k1)-b(k2)
      k1=k1+kspan
      if(k1 .lt. k2) goto 650
      a(kk)=ak
      b(kk)=bk
      k1=kk
      k2=kk+kspnn
      j=1
  660 k1=k1+kspan
      k2=k2-kspan
      jj=j
      ak=aa
      bk=bb
      aj=0.0
      bj=0.0
      k=1
  670 k=k+1
      ak=at(k)*ck(jj)+ak
      bk=bt(k)*ck(jj)+bk
      k=k+1
      aj=at(k)*sk(jj)+aj
      bj=bt(k)*sk(jj)+bj
      jj=jj+j
      if(jj .gt. jf) jj=jj-jf
      if(k .lt. jf) goto 670
      k=jf-j
      a(k1)=ak-bj
      b(k1)=bk+aj
      a(k2)=ak+bj
      b(k2)=bk-aj
      j=j+1
      if(j .lt. k) goto 660
      kk=kk+kspnn
      if(kk .le. nn) goto 640
      kk=kk-nn
      if(kk .le. kspan) goto 640
c  multiply by rotation factor (except for factors of 2 and 4)
  700 if(i .eq. m) goto 800
      kk=jc+1
  710 c2=1.0-cd
      s1=sd
  720 c1=c2
      s2=s1
      kk=kk+kspan
  730 ak=a(kk)
      a(kk)=c2*ak-s2*b(kk)
      b(kk)=s2*ak+c2*b(kk)
      kk=kk+kspnn
      if(kk .le. nt) goto 730
      ak=s1*s2
      s2=s1*c2+c1*s2
      c2=c1*c2-ak
      kk=kk-nt+kspan
      if(kk .le. kspnn) goto 730
      c2=c1-(cd*c1+sd*s1)
      s1=s1+(sd*c1-cd*s1)
      kk=kk-kspnn+jc
      if(kk .le. kspan) goto 720
      kk=kk-kspan+jc+inc
      if(kk .le. jc+jc) goto 710
      goto 100
c  permute the results to normal order---done in two stages
c  permutation for square factors of n
  800 np(1)=ks
      if(kt .eq. 0) goto 890
      k=kt+kt+1
      if(m .lt. k) k=k-1
      j=1
      np(k+1)=jc
  810 np(j+1)=np(j)/nfac(j)
      np(k)=np(k+1)*nfac(j)
      j=j+1
      k=k-1
      if(j .lt. k) goto 810
      k3=np(k+1)
      kspan=np(2)
      kk=jc+1
      k2=kspan+1
      j=1
      if(n .ne. ntot) goto 850
c  permutation for single-variate transform (optional code)
  820 ak=a(kk)
      a(kk)=a(k2)
      a(k2)=ak
      bk=b(kk)
      b(kk)=b(k2)
      b(k2)=bk
      kk=kk+inc
      k2=kspan+k2
      if(k2 .lt. ks) goto 820
  830 k2=k2-np(j)
      j=j+1
      k2=np(j+1)+k2
      if(k2 .gt. np(j)) goto 830
      j=1
  840 if(kk .lt. k2) goto 820
      kk=kk+inc
      k2=kspan+k2
      if(k2 .lt. ks) goto 840
      if(kk .lt. ks) goto 830
      jc=k3
      goto 890
c  permutation for multivariate transform
  850 k=kk+jc
  860 ak=a(kk)
      a(kk)=a(k2)
      a(k2)=ak
      bk=b(kk)
      b(kk)=b(k2)
      b(k2)=bk
      kk=kk+inc
      k2=k2+inc
      if(kk .lt. k) goto 860
      kk=kk+ks-jc
      k2=k2+ks-jc
      if(kk .lt. nt) goto 850
      k2=k2-nt+kspan
      kk=kk-nt+jc
      if(k2 .lt. ks) goto 850
  870 k2=k2-np(j)
      j=j+1
      k2=np(j+1)+k2
      if(k2 .gt. np(j)) goto 870
      j=1
  880 if(kk .lt. k2) goto 850
      kk=kk+jc
      k2=kspan+k2
      if(k2 .lt. ks) goto 880
      if(kk .lt. ks) goto 870
      jc=k3
  890 if(2*kt+1 .ge. m) return
      kspnn=np(kt+1)
c  permutation for square-free factors of n
      j=m-kt
      nfac(j+1)=1
  900 nfac(j)=nfac(j)*nfac(j+1)
      j=j-1
      if(j .ne. kt) goto 900
      kt=kt+1
      nn=nfac(kt)-1
      if(nn .gt. maxp) goto 998
      jj=0
      j=0
      goto 906
  902 jj=jj-k2
      k2=kk
      k=k+1
      kk=nfac(k)
  904 jj=kk+jj
      if(jj .ge. k2) goto 902
      np(j)=jj
  906 k2=nfac(kt)
      k=kt+1
      kk=nfac(k)
      j=j+1
      if(j .le. nn) goto 904
c  determine the permutation cycles of length greater than 1
      j=0
      goto 914
  910 k=kk
      kk=np(k)
      np(k)=-kk
      if(kk .ne. j) goto 910
      k3=kk
  914 j=j+1
      kk=np(j)
      if(kk .lt. 0) goto 914
      if(kk .ne. j) goto 910
      np(j)=-j
      if(j .ne. nn) goto 914
      maxf=inc*maxf
c  reorder a and b, following the permutation cycles
      goto 950
  924 j=j-1
      if(np(j) .lt. 0) goto 924
      jj=jc
  926 kspan=jj
      if(jj .gt. maxf) kspan=maxf
      jj=jj-kspan
      k=np(j)
      kk=jc*k+ii+jj
      k1=kk+kspan
      k2=0
  928 k2=k2+1
      at(k2)=a(k1)
      bt(k2)=b(k1)
      k1=k1-inc
      if(k1 .ne. kk) goto 928
  932 k1=kk+kspan
      k2=k1-jc*(k+np(k))
      k=-np(k)
  936 a(k1)=a(k2)
      b(k1)=b(k2)
      k1=k1-inc
      k2=k2-inc
      if(k1 .ne. kk) goto 936
      kk=k2
      if(k .ne. j) goto 932
      k1=kk+kspan
      k2=0
  940 k2=k2+1
      a(k1)=at(k2)
      b(k1)=bt(k2)
      k1=k1-inc
      if(k1 .ne. kk) goto 940
      if(jj .ne. 0) goto 926
      if(j .ne. 1) goto 924
  950 j=k3+1
      nt=nt-kspnn
      ii=nt-inc+1
      if(nt .ge. 0) goto 924
      return
c  Error finish, insufficient array storage
998   write(*,*)' Array bounds exceeded within subroutine fft'
      ierr=1
      end
c_______________________________________________________________________
c=======================================================================
c           Unit S: Splines
c=======================================================================
      subroutine spline(nn, x, u, s, a)
c$$$$ calls nothing
c  Finds array  s  for spline interpolator  eval.
c  Based on an alogorithm by Richard F. Thompson in 'Spline
c  interpolation on a digital computer' (x-692-70-261) 1970 Goddard
c  Space Flight Center.  Note these splines are not 'natural' since the
c  2nd derivative is not zero x(1) and x(n).
c  nn  number of data points supplied (may be negative, see below)
c  x()  nn-array of x-coordinates where function is sampled.  The
c      sequence xx(1),xx(2)...  must be a strictly increasing sequence
c      but THIS IS NOT CHECKED by spline.
c  u() nn-vector of function values that are to be interpolated.
c  s() nn-vector output of 2nd derivative at sample points.
c  a() nn-vector, working space.
c
c  If the user wishes derivatives at the ends of series to assume
c  specified values, du(1)/dx, du(nn)/dx  must be placed in  s1, s2
c  and in the call nn=-number of terms in series.  Normally a parabola
c  is fitted through the 1st and last 3 points to find the slopes.
c
c  If 2 points are supplied, a straight lines is fitted.
c  If 3 points are supplied, a parabola is fitted.
c
c  At evaluation time, straight-line extrapolation is provided
c  outside the interval (x(1), x(n)).
      common /startx/ istart
      dimension x(*),u(*),s(*),a(*)
ce
      q(u1, x1, u2, x2)=(u1/x1**2-u2/x2**2)/(1.0/x1-1.0/x2)
c
c  Assign derivatives at the ends.
      istart=0
      n=abs(nn)
c  Series too short for cubic spline - use straight lines.
      if (n.le. 2) then
        s(1)=0.0
        s(2)=0.0
        return
      endif
c
      q1=q(u(2)-u(1), x(2)-x(1), u(3)-u(1), x(3)-x(1))
      qn=q(u(n-1)-u(n), x(n-1)-x(n), u(n-2)-u(n), x(n-2)-x(n))
c  Three points - fit a parabola
      if (n.eq. 3) then
        s(1)=(qn - q1)/(x(3) - x(1))
        s(2)=s(1)
        s(3)=s(1)
        return
      endif
c  Force the derivatives to take user-defined values
      if (nn.le.0) then
        q1=s(1)
        qn=s(2)
      endif
      s(1)=6.0*((u(2)-u(1))/(x(2)-x(1)) - q1)
      n1= n - 1
      do 2000 i=2,n1
        s(i)= (u(i-1)/(x(i)-x(i-1)) - u(i)*(1.0/(x(i)-x(i-1))+
     +  1.0/(x(i+1)-x(i))) + u(i+1)/(x(i+1)-x(i)))*6.0
 2000 continue
      s(n)=6.0*(qn + (u(n1)-u(n))/(x(n)-x(n1)))
      a(1)=2.0*(x(2)-x(1))
      a(2)=1.5*(x(2)-x(1)) + 2.0*(x(3)-x(2))
      s(2)=s(2) - 0.5*s(1)
      do 3000 i=3,n1
        c=(x(i)-x(i-1))/a(i-1)
        a(i)=2.0*(x(i+1)-x(i-1)) - c*(x(i)-x(i-1))
        s(i)=s(i) - c*s(i-1)
 3000 continue
      c=(x(n)-x(n1))/a(n1)
      a(n)=(2.0-c)*(x(n)-x(n1))
      s(n)=s(n) - c*s(n1)
c  Back substitiute.
      s(n)= s(n)/a(n)
      do 4000 j=1,n1
        i=n-j
        s(i) =(s(i) - (x(i+1)-x(i))*s(i+1))/a(i)
 4000 continue
      return
      end
c_____________________________________________________________________
      function eval(y, nn, x, u, s)
c$$$$ calls nothing
c  Performs cubic spline interpolation of a function sampled unequally
c  in  x.  The routine  spline  must be called to set up the array s.
c  See comments in  spline.
c  y   the coordinate at which function value is desired.
c  nn  number of samples of original function.
c  x() nn-vector of x values in original data sample.
c  u() nn-vector of function values.
c  s() nn-vector of 2nd derivatives at sample points.  Found by  spline
c      which must be called once before beginning interpolation.
c  If  y  falls outside  (x(1), x(nn))  extrapolation is based upon
c  the tangent straight line at the endpoint.
      common /startx/ istart
      dimension x(*),u(*),s(*)
      data l0/1/
ce
      n=abs(nn)
      if (istart.gt.n .or. istart.lt.1)
     $    istart=1+(n-1)*(y-x(1))/(x(n)-x(1))
      if (y .le. x(1))  go to 3000
      if (y .ge. x(n)) go to 3100
c  Locate interval (x(l0),x(l1))  of y.
      if (y .lt. x(istart)) goto 1200
c  Scan up the x array.
      do 1100 l1=istart, n
        if (x(l1) .gt. y) go to 1150
 1100 continue
 1150 l0=l1 - 1
      go to 1500
c  Scan downwards in  x  array.
 1200 do 1300 l1=1, istart
        l0=istart - l1
        if (x(l0) .le. y) go to 1350
 1300 continue
 1350 l1=l0 + 1
 1500 istart=l0
c  Evaluate interpolation.
      xi=(y - x(l0))/(x(l1) - x(l0))
      h6=(x(l1) - x(l0))**2/6.0
      eval=u(l0) + xi*(u(l1)-u(l0) - h6*(1.0-xi)*(s(l1)*(1.0+xi)
     $     + s(l0)*(2.0-xi)))
      return
c  Out of range.  Substitute straight-line extrapolant.
 3000 h=x(2) - x(1)
      eval=u(1) + (y-x(1))*((u(2)-u(1))/h - h*(s(2)-2.0*s(1))/6.0)
      return
 3100 h=x(n) - x(n-1)
      eval=u(n)+(y-x(n))*((u(n)-u(n-1))/h+h*(2.0*s(n)-s(n-1))/6.0)
      return
      end
c_____________________________________________________________________
      subroutine akima(n, x,f, df)
c$$$$ calls nothing
c  Given the array of samples of a function f and sample points x,
c  assumed to be montone, generates slopes for an interpolating rule
c  according to Akima's algorithm (Lancaster and Salkauskas,
c  Curve and surface fitting, 1986 academic press, p 82).
c
      dimension x(n),f(n),df(n),S(4)
c
      q(u1,x1,u2,x2)=(u1/x1**2-u2/x2**2)/(1.0/x1-1.0/x2)
c
ce
      S(1)=(f(2) - f(1))/(x(2) - x(1))
      S(2)=S(1)
      S(3)=S(1)
      S(4)=(f(3) - f(2))/(x(3) - x(2))
      Sn=(f(n) - f(n-1))/(x(n) - x(n-1))
      eps=1.0e-6*abs(x(n) - x(1))
      do 1200 i=1, n
        D1=abs(S(2) - S(1))
        D2=abs(S(4) - S(3))
        df(i)=(D2*S(2) + D1*S(3))/(D1 + D2 + eps)
c
        S(1)=S(2)
        S(2)=S(3)
        S(3)=S(4)
        if (i+3 .le. n) then
          S(4)=(f(i+3) - f(i+2))/(x(i+3)-x(i+2))
        else
          S(4)=Sn
        endif
 1200 continue
c
      if (n.eq. 2) return
c  If 3 or more points use gradient from a parabola for 1st & last df
      df(1)=q(f(2)-f(1),x(2)-x(1),f(3)-f(1),x(3)-x(1))
      df(n)=q(f(n-1)-f(n),x(n-1)-x(n),f(n-2)-f(n),x(n-2)-x(n))
      return
      end
c____________________________________________________________________
      function evlak(y, n,x,f,df)
c$$$$ calls nothing
c  Given a montone increasing array of x values (Note: routine does not
c  check this), samples of values of the function f and 1st derivative
c  df at the xs, interpolates by cubic polynomial.  The array df can be
c  found by calling subroutine akima.
c
c  If y falls outside the interval, the polynomial on the
c  nearest interval is used to extrapolate.
      dimension x(n), f(n), df(n)
c
c  Search for proper interval initiated at previous call if possible.
      save init
      data init/1/
c
c  Locate sample interval containing y:  after y lies in [x(init),
c  x(init+1)), unless y lies outside [x(1), x(n)] when the interval
c  is the intervals containing the apprpriate end point.
ce
      init=min(init, n)
      if (y .gt. x(init)) then
         do 1100 k=init, n
           if (x(k) .gt. y) then
             init=k-1
             goto  1300
           endif
 1100    continue
         init=n-1
       else
         do 1200 k=init, 1, -1
           if (x(k) .le. y)  then
             init=k
             goto 1300
           endif
 1200   continue
        init=1
      endif
 1300 dx=x(init+1) - x(init)
c
c  Evaluate the cubic interpolator
      t=(y - x(init))/dx
      s=1.0 - t
      evlak=s**2*((1.0 + 2.0*t)*f(init) + t*dx*df(init)) +
     $      t**2*((1.0 + 2.0*s)*f(init+1) - s*dx*df(init+1))
      return
      end
c____________________________________________________________________
      subroutine prolat(ntap, tbp, kpt)
c$$$$ calls staper, fft
c  Finds and applies the Thomson-Slepian prolate spheroidal tapers
c  ntap is the number of tapers to be used
c  tbp  is the time-bandwidth product of the frequency interval
c  for averaging.
c  kpt if > 0 is the unit to write out tapers
      parameter (mxx=500 000)
      parameter (nfmx=mxx + 2)
      common /series/ nx,dt,t1,x(mxx)
      common /result/ nt,nf,fNyq,df,sx(nfmx),kopt(nfmx)
c
      double precision vn
      common /pool/ vn(4*mxx + 2)
      dimension vnx(2,nfmx)
c
c  Find the set of tapers  - results in vn
ce
      mxtap=(4*mxx+2)/(nt+1)-2
      if (ntap .gt. mxtap) then
        write(*,'(a,i5)')
     $'>>> Memory limitations reduce number of tapers to:',mxtap
        ntap=mxtap
      else
        write(*,'(12x,a,i5)')'Number of prolate tapers used:',ntap
      endif
      ntn1=(nt+1)*ntap + 1
      ntn2=ntn1 + nt + 2
      ndim=nt+1
      call staper(nt,dble(tbp),ntap, vn,ndim, vn(ntn1),vn(ntn2))
c
      write(*,'(3x,a/(4f18.14))')'Prolate spheroidal eigenvalues:',
     $(vn(ntn1+j-1),j=1,ntap)
c
c  Zero out the target array and fill kopt()
      do 1100 j=1, nf
         sx(j)=0.0
         kopt(j)=ntap
 1100 continue
c
c  Loop over individual tapers
      do 1500 ktap=1, ntap
        kj=(ktap-1)*ndim
        do 1200 j=1, nt
          vnx(1,j)=vn(j+kj)*x(j)
          vnx(2,j)=0.0
 1200   continue
        if (kpt.gt.0)write(kpt,'(2i5,g14.6)')(j,ktap,vn(j+kj),j=1,nt)
c
c  In-place mixed-radix FFT of (current taper * x)
        call fft(vnx,vnx(2,1), nt,nt,nt, 2, ierr)
c
c  Accumulate the spectrum in sx
        do 1300 j=1, nf
          sx(j)=sx(j) + vnx(1,j)**2 + vnx(2,j)**2
 1300   continue
 1500 continue
c
c  Normalize
      sc=2.0/ntap
      do 1600 j=1, nf
        sx(j)=sc*sx(j)
 1600 continue
c
      return
      end
c______________________________________________________________________
      subroutine period
c$$$$ calls fft
c  Finds primitive periodogram estimate of PSD
      parameter (mxx=500 000)
      parameter (nfmx=mxx + 2)
      common /series/ nx,dt,t1,x(mxx)
      common /result/ nt,nf,fNyq,df,sx(nfmx),kopt(nfmx)
c
      common /pool/ vn(2,4*mxx + 2)
c
c  Find the set of tapers  - results in vn
ce
      do 1100 j=1, nt
         vn(1,j)=x(j)
         vn(2,j)=0.0
         kopt(j)=1
 1100 continue
c
c  In-place mixed-radix FFT of x
      call fft(vn,vn(2,1), nt,nt,nt, 2, ierr)
c
c  Accumulate 1-sided spectrum in sx
      do 1300 j=1, nf
        sx(j)=2.0*(vn(1,j)**2 + vn(2,j)**2)
 1300 continue
c
      return
      end
c______________________________________________________________________
c=======================================================================
c           Unit P: Prolates
c=======================================================================
      subroutine staper(nt, fw, nev, v, ndim, a, w)
c$$$$ calls tsturm
c  Slepian - Thomson multi-taper procedure
c  Slepian, D.     1978  Bell Sys Tech J v57 n5 1371-1430
c  Thomson, D. J.  1982  Proc IEEE v70 n9 1055-1096
c    nt    the number of points in the series
c    fw    the time-bandwidth product (number of Rayleigh bins)
c    nev   the desired number of tapers
c    v     the eigenvectors (tapers) are returned in v(.,nev)
c    a, w  work arrays dimensioned at least nt long (nt+1, nt odd)
c    a(1..nev) contains bandwidth retention factors on output.
c  The tapers are the eigenvectors of the tridiagonal matrix sigma(i,j)
c  [see Slepian(1978) eq 14 and 25.] They are also the eigenvectors of
c  the Toeplitz matrix eq. 18. We solve the tridiagonal system in
c  tsturm for the tapers and use them in Slepians eq 18 to get the
c  bandwidth retention factors (i.e. the eigenvalues) Thomson's
c  normalisation is used with no attention to sign.
      implicit real*8(a-h,o-z)
      dimension a(*),w(*),v(ndim,*)
      parameter (pi=3.14159265358979d0,r2=1.414213562373095d0)
ce
      if(nt.lt.2) return
      nx=mod(nt,2)
      lh=(nt/2)+nx
      lp1=nt+1
      om=2.d0*pi*fw/nt
      com=cos(om)
      hn=0.5d0*lp1
      do 10 i=1,lh
        a(i)=com*(i-hn)**2
   10   w(i)=0.5d0*float(i*(nt-i))
      if(nx.eq.0) then
        asav=a(lh)-w(lh)
        a(lh)=a(lh)+w(lh)
        rbd=1.d0/(a(lh)+w(lh-1))
      else
        asav=w(lh-1)
        rbd=1.d0/(w(lh)+w(lh-1))
        w(lh-1)=r2*w(lh-1)
      endif
      do 15 i=1,lh
        a(i+lh)=w(i)*rbd
        w(i)=a(i+lh)**2
   15   a(i)=a(i)*rbd
      neven=max0((nev+1)/2,1)
      nodd=nev-neven
c  Do the even tapers
      call tsturm(nt,lh,a,a(lh+1),w,neven,v,ndim,w(lh+1),0)
      do 20 i=1,neven
        k=2*i-1
        if(nx.eq.1) v(lh,k)=r2*v(lh,k)
          do 20 j=1,lh
   20     v(lp1-j,k)=v(j,k)
      if(nodd.le.0) goto 34
c  Do the odd tapers
      if(nx.eq.0) then
        a(lh)=asav*rbd
      else
        a(nt)=asav*rbd
        w(lh-1)=asav*asav
      endif
      call tsturm(nt,lh-nx,a,a(lh+1),w,nodd,v,ndim,w(lh+1),1)
      do 30 i=1,nodd
        k=2*i
        if(nx.eq.1) v(lh,k)=0.d0
          do 30 j=1,lh
   30     v(lp1-j,k)=-v(j,k)
   34 ntot=neven+nodd
c  Calculate bandwidth retention parameters
      dc=2.d0*com
      sm=0.d0
      s=sin(om)
      w(1)=om/pi
      w(2)=s/pi
      do 35 j=3,nt
        sn=dc*s-sm
        sm=s
        s=sn
   35   w(j)=s/(pi*(j-1))
      do 55 m=1,ntot
        vmax=abs(v(1,m))
        kmax=1
        do 40 kk=2,lh
          if(abs(v(kk,m)).le.vmax) goto 40
          kmax=kk
          vmax=abs(v(kk,m))
   40     continue
        a(m)=0.d0
        nlow=kmax-1
          do 45 j=1,nlow
   45     a(m)=a(m)+w(j+1)*v(nlow+1-j,m)
        nup=nt-nlow
          do 50 j=1,nup
   50     a(m)=a(m)+w(j)*v(nlow+j,m)
   55 a(m)=a(m)/v(kmax,m)
      return
      end
c______________________________________________________________________
      subroutine tsturm(nt,n,a,b,w,nev,r,ndim,ev,ipar)
c$$$$ calls root
c  Uses bisection and Sturm counting to isolate the eigenvalues of the
c  symmetric tridiagonal matrix with main diagonal a(.) and sub/super
c  diagonal b(.).  Newton's method is used to refine the eigenvalue in
c  subroutine root then direct recursion is used to get the eigenvector
c  as this is always stable.  Note  ipar=0 for even tapers   =1 for odd
c  tapers
      implicit real*8(a-h,o-z)
      parameter (eps=1.d-15,eps1=5.d-15)
      dimension a(*),b(*),ev(*),w(*),r(ndim,*)
ce
      if(n.le.0.or.nev.le.0) return
      umeps=1.d0-eps
      do 5 i=1,nev
    5 ev(i)=-1.d0
      u=1.d0
      do 1000 ik=1,nev
      if(ik.gt.1) u=ev(ik-1)*umeps
      el=min(ev(ik),u)
   10 elam=.5d0*(u+el)
      if(abs(u-el).le.eps1) goto 35
      iag=0
      q=a(1)-elam
      if(q.ge.0.d0) iag=iag+1
      do 15 i=2,n
      if(q.eq.0.d0) x=abs(b(i-1))/eps
      if(q.ne.0.d0) x=w(i-1)/q
      q=a(i)-elam-x
      if(q.ge.0.d0) iag=iag+1
      if(iag.gt.nev) goto 20
   15 continue
      if(iag.ge.ik) go to 20
      u=elam
      go to 10
   20 if(iag.eq.ik) go to 30
      m=ik+1
      do 25 i=m,iag
   25 ev(i)=elam
      el=elam
      go to 10
   30 el=elam
      call root(u,el,elam,a,b,w,n,ik)
   35 ev(ik)=elam
      jk=2*ik+ipar-1
      r(1,jk)=1.d0
      r(2,jk)=-(a(1)-ev(ik))/b(1)
      ddot=1.d0+r(2,jk)*r(2,jk)
      jm1=2
      do 45 j=3,n
      r(j,jk)=-((a(jm1)-ev(ik))*r(jm1,jk)+b(j-2)*r(j-2,jk))/b(jm1)
      ddot=ddot+r(j,jk)*r(j,jk)
   45 jm1=j
      rnorm=sqrt(nt/(2.d0*ddot))
      do 50 j=1,n
   50 r(j,jk)=r(j,jk)*rnorm
 1000 continue
      return
      end
c______________________________________________________________________
      subroutine root(u,el,elam,a,bb,w,n,ik)
c$$$$ calls nothing
      implicit real*8(a-h,o-z)
      parameter (eps=1.d-15,eps1=5.d-15)
      dimension a(*),bb(*),w(*)
ce
    5 elam=.5d0*(u+el)
   10 if(abs(u-el).le.1.5d0*eps1) return
      an=a(1)-elam
      b=0.d0
      bn=-1.d0/an
      iag=0
      if(an.ge.0.d0) iag=iag+1
      do 20 i=2,n
      if(an.eq.0.d0) x=abs(bb(i-1))/eps
      if(an.ne.0.d0) x=w(i-1)/an
      an=a(i)-elam-x
      if(an.eq.0.d0) an=eps
      bm=b
      b=bn
      bn=((a(i)-elam)*b-bm*x-1.d0)/an
      if(an.ge.0.d0) iag=iag+1
   20 continue
      if(iag.eq.ik) goto 25
      u=elam
      goto 30
   25 el=elam
   30 del=1.d0/bn
      if(abs(del).le.eps1) del=sign(eps1,del)
      elam=elam-del
      if(elam.ge.u.or.elam.le.el) goto 5
      goto 10
      end
c______________________________________________________________________
c=======================================================================
c        Unit K:    CLIP Command Line Interface Programs
c=======================================================================
      subroutine ascan(kase)
c$$$$ calls clear getchr icheck ljust
c  Command input routine
c
c  Reads from the standard input until eof or code "execute " or "quit".
c  Saves lines in /store/ for later retrieval by getarr, getone or getchr.
c
c  Prints a glossary upon request
c
      parameter (inmx=500)
      character *80 input(inmx),line, code*4
      common /store/ input
      common /ndict/ iecho,nin
      common /long/ lg(inmx)
c
      if (nin.eq. 0) write(*,'(a)') ' ',
     $'Enter commands for spectral analysis (? for help)'
c
      do 1500 l=nin+1, inmx
        read(*,'(80a)', end=3000) line
        if (line(1:1).eq. ' ' .or. line(1:1).eq. '%') goto 1500
        if (line(1:4).eq.'exec') goto 2000
c  List a glossary of codes
        if (line(1:1).eq. '?') then
          write(*, '(a/a/(2x,a))')
     $'Enter commands from the following list:',
     $'?:        Remind me of the command list again',
     $'execute:           Quit reading commands and begin calculations',
     $'file filename:        Enter filename of time series (mandatory)',
     $'file -continue:             Read further data from current file',
     $'output filename:               Provide the filename for results',
     $'plot   filename:               Provide filename for plotxy file',
     $'interval dt:                          Sampling interval of data',
     $'units time signal:        Give measurement units of independent',
     $'                         variable and of signal (eg: sec volts)',
     $'scale fac:             Multiply signal values by the factor fac',
     $'terms n1 n2:               Read from term n1 to n2 in the file',
     $'nterms n:            Number of data points to be read from file',
     $'skip k:                   Skip over k lines before reading data',
     $'column k [n]:        Read data from column k [and n] in a table',
     $'                            (2 columns only when interpolating)',
     $'spline [kind]:    Input t, x pairs & spline interpolate to even',
     $'                     time interval dt;  kind = natural or Akima',
     $'validate filename:          Write interpolated series to a file'
          write(*, '(2x,a)')
     $'decimate m [nw]:        Decimate input series by the factor m,',
     $'                            using nw anti-alias filter weights',
     $'detrend:                        If present, series is detrended',
     $'superimpose:               Plot next spectrum on previous graph',
     $'hold:           While on, save up plots of superimposed spectra',
     $'detail f1 f2:               Plot spectra only in f1 < freq < f2',
     $'replot f1 f2:         Replot spectrum from previous calculation',
     $'title text:                         Identify plot with the text',
     $'review:                           Display current command stack',
     $'adapt ntimes:                     Apply adaptive process ntimes',
     $'smooth s:         Increase s>1, or decrease s<1, PSD smoothness',
     $'tapers nt:               Use exactly nt tapers (overides adapt)',
     $'tapers < nt:        Use no more than nt tapers in adaptive code',
     $'prewhiten nd:                    Prewhiten with filter order nd',
     $'save file:                 Save prewhitened time series to file',
     $'autocorr [file]:     Calculate autocorrelation fn, save in file',
     $'clear command:            Delete all occurrences of the command',
     $'prolate [file]           Use prolate tapers, write them to file',
     $'time-bp p:    Specify time-bandwidth product for prolate tapers',
     $' '
c
c  Review the command stack
        elseif (line(1:4).eq. 'revi') then
          write(*,'(5x,a)')' ', '=================== ',
     $    (input(i)(1:lg(i)),i=1,nin),'=================== ',' '
c
        else
c  Translate homonyms
          if (line(1:4).eq. 'data') line(1:4)='file'
          if (line(1:4).eq. 'cols') line(1:4)='colu'
          if (line(1:4).eq. 'echo') iecho=-iecho
          nin=nin + 1
          if (nin.gt. inmx) then
            write(*,'(a)') '>>>> Too many commands - memory full'
            stop
          endif
c  Check for recognized commands
          call icheck(line)
c  Special treatment of '<' in 'tapers': puts '<' after number
          less = index(line, '<')
          if (line(1:4).eq. 'tape' .and. less.ne. 0) then
              call ljust(line,line, last)
              line(less:less)=' '
              line(last+2:last+2) = '<'
          endif
c  Save current command line in input store
          call ljust(line,input(nin), lg(nin))
c
          if (line(1:4).eq. 'quit') return
c  Clear a command and clear clear itself
          if (line(1:4).eq. 'clea') then
            call getchr('clear', code, ignor)
            call clear (code)
            nin=nin - 1
          endif
        endif
 1500 continue
c
 2000 continue
      n1=max(1, nin-24)
      if (kase.gt. 1) write(*, '(a,i4)') '--------- Spectrum',kase
      write(*,'(5x,a)') '=================== ',
     $(input(i)(1:lg(i)),i=n1,nin),'=================== ',' '
      return
c
 3000 stop
      end
c______________________________________________________________________
      subroutine icheck(line)
c$$$$ calls nothing
c  Checks the command fragment com against the catalog; appends a
c  warning if com is not present in the list.
      character*80 line, com*4
ce
      com=line(1:4)
      if (0 .eq.
     $ index('adap auto clea colu deci deta detr file hold inte ',com)
     $+index('nter outp plot prew repl skip smoo spli scal supe',com)
     $+index('tape term titl unit vali logf prol time band save quit',
     $ com))
     $ line=line(1:20)//' %<<<<<<< Unrecognized command'
      return
      end
c______________________________________________________________________
      subroutine getarr(code, values, nwant, nfound)
c$$$$ calls ljust
c
c  Extracts an array of numbers from the input store.  That store is
c  a large array in common /store/ which has been filled earlier.
c
c  code    A 4-bye identifying code.  Only lines in the input store
c          beginning with this code are scanned for numbers.
c  values  the real output array of values found.
c  nwant   the maximum number of numbers expected .
c  nfound  the number of numbers actually found in the input.
c          If the line contains fewer than nwant  values, this is the
c          value returned in  nfound.  If an error is discovered
c          nfound=-n, where  n  is the number of numbers properly
c          decoded.  If there are no numbers after the codeword
c          nfound=0.  Finally, if the code is absent from the store
c          nfound=-99 and the array is left undisturbed.
c
      parameter (inmx=500)
      character *80 input(inmx),line,local,char, code*4
      common /store/ input
      common /ndict/ iecho,nin
      dimension values(*)
c
c  Read the store in reverse order (Thus last entry is obeyed)
ce
      do 1010 lin=nin, 1, -1
        line=input(lin)
c  Check for code
        if (code .eq. line(1:4)) goto 1020
 1010 continue
c  Code word not found
      nfound=-99
      return
c
 1020 if (iecho.ge. 1) write(*,'(2a)')'==> ',line
      n1=index(line, ' ')+1
      n2=index(line, '%')
      n2=80 + min(n2,1)*(n2 - 81)
      char=line(n1:n2)
c
      do 1500 n=1, nwant
        call ljust(char, local, ignor)
        lbl=index(local, ' ')
        if (lbl.eq. 1) goto 1510
        read (local, *, err=2000) values(n)
        char=local(lbl:80)
 1500 continue
      n=nwant+1
 1510 nfound=n - 1
      return
c
 2000 write(*,'(a)')' ',
     $'>>> Unreadable numbers in this input line:',line
      nfound=1 - n
      return
      end
c______________________________________________________________
      subroutine getone(code, value, nfound)
c$$$$ calls getarr
c
c  Extracts a single number from the input store.  That store is
c  a large array in common /store/ which has been filled earlier.
c
c  code    A 4-bye identifying code.  Only lines in the input store
c          beginning with this code are scanned for numbers.
c  value   the real output variable containing the desired number.
c  nfound  is 1 if a number is successfully read in; it is zero
c          the number is absent or unreadable.  nfound = -99 if the
c          code is absent from the input store.
c
      character*4 code
      dimension v(1)
ce
      call getarr(code, v, 1, nfound)
      if (nfound.eq. 1) value=v(1)
      return
      end
c______________________________________________________________
      subroutine getint(code, number, nfound)
c$$$$ calls getarr
c
c  Extracts a single integer from the input store.
c  See getone for details.
c
      character*4 code
      dimension v(1)
ce
      call getarr(code, v, 1, nfound)
      if (nfound.eq. 1) then
        number=nint(v(1))
        if (number .ne. v(1)) write(*,'(4a,1p,g16.8/a,i10)')
     $  '>>> Warning: ',code,' expects an integer argument, but',
     $  ' found: ',v(1),'    Returns: ',number,' '
      endif
      return
      end
c______________________________________________________________
      subroutine getchr(code, char, nbytes)
c$$$$ calls ljust
c  Extracts a single character variable from the input store.  That
c  store is a large array in common /store/ which has been filled
c  earlier.
c
c  code    A 4-bye identifying code.  Only lines in the input store
c          beginning with this code are scanned.
c  char    the character output variable containing the desired string.
c  nbytes  is the length of the string read; it is zero if
c          the line is blank after the code.  nbytes = -99 if the
c          code is absent from the input store.
c
      parameter (inmx=500)
      character *80 input(inmx),line, char*(*),  code*4
      common /store/ input
      common /ndict/ iecho,nin
c  Inspect the store in reverse order (thus reading latest entry)
ce
      do 1010 lin=nin, 1, -1
        line=input(lin)
c  Check for code word
        if (code .eq. line(1:4)) goto 1020
 1010 continue
c  Code word not found
      nbytes=-99
      return
c
 1020 if (iecho.ge. 1) write(*,'(2a)')'==> ',line
c
      n1=index(line, ' ')+1
      n2=index(line, '%')
      n2=80 + min(n2,1)*(n2 - 81)
c  Check for blank string
      nbytes=0
      if (line(n1:n2).eq. ' ') return
c  Save in char everything from 1st non-blank to last non-blank
      call ljust(line(n1:n2), char, nbytes)
      return
      end
c______________________________________________________________
      subroutine clear(code)
c$$$$ calls nothing
c
c  Removes the last command entry identified by code.
c  code    A 4-byte identifying code.  Only lines in the input store
c          beginning with this code are scanned.
c
      parameter (inmx=500)
      character *80 input(inmx),line, code*4
      common /store/ input
      common /ndict/ iecho,nin
c
c  Inspect the store in normal order
ce
      do 1010 lin= nin, 1, -1
        line=input(lin)
c  Check code  and clear the line if present, for every instance
        if (code.eq. line(1:4)) then
          input(lin)=' --'
c         return
        endif
 1010 continue
      return
      end
c______________________________________________________________
      subroutine ljust(str1, str2, len2)
c$$$$ calls nothing
c  Copy str1 (up to 80 bytes in length) into str2, left justified; len2
c  is length of str2 with trailing blanks deleted.
c  str1 and str2 may be the same variable in the call.
      character*(*) str1, str2, str3*80
c
ce
      str2=str1
      len2=1
      if (str1.eq. ' ') return
c
      l1=len(str1)
      str3=str1
      do 1010 j=1, l1
        if (str1(j:j) .ne. ' ') goto 1011
 1010 continue
 1011 str2=str3(j:l1)
      do 2100 len2 = len(str2),1,-1
        if (str2(len2:len2).ne. ' ') return
 2100 continue
      end
c______________________________________________________________
      blockdata iounit
      common /ndict/ iecho,nin
c  Don't echo command lines when read
      data iecho/-1/
c  Initial command line count
      data nin/0/
      end
c______________________________________________________________________
      function later(code1, code2)
c$$$$ calls nothing
c  Returns the difference in order number in the stack of the
c  most recent occurrence of the commands code1, code2.  If the
c  command isn't present assigns it the order number zero.
c  Thus later > 0 if code1 occurs after code 2, < 0 if before,
c  and later=0 if both are absent and code1 and code2 are different.
      parameter (inmx=500)
      character *80 input(inmx),line, code1*4,code2*4
      common /store/ input
      common /ndict/ iecho,nin
c  Inspect the store
ce
      kode1=0
      kode2=0
      do 1010 lin=1, nin
        line=input(lin)
        if (code1 .eq. line(1:4)) kode1=lin
        if (code2 .eq. line(1:4)) kode2=lin
 1010 continue
c
      later=kode1 - kode2
      return
      end
c______________________________________________________________
