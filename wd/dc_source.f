c  This is the Differential Corrections Main Program.                      
c                                                                          
c  Version of April 11, 2016
c                                             
c     This version can iterate subsets and also the base set.
C                                                                          
C     PARAMETER NUMBER 9 IS A, THE RELATIVE ORBITAL SEMI-MAJOR AXIS, IF    
C     SIMULTANEOUS LIGHT AND VELOCITY SOLUTIONS ARE BEING DONE. HOWEVER,   
C     IF ONLY VELOCITY CURVES ARE BEING SOLVED, PARAMETER 9 WILL           
C     EFFECTIVELY BE A*SIN(I), PROVIDED AN INCLINATION OF 90 DEG. IS       
C     ENTERED. CONCEIVABLY,IN SOME RARE SITUATIONS IT MAY BE POSSIBLE TO   
C     FIND A AND I SEPARATELY, USING ONLY VELOCITY CURVES. THIS COULD BE   
C     THE CASE IF THE VELOCITY PROXIMITY EFFECTS ARE IMPORTANT.            
C                                                                          
C     OTHER PROGRAM UNITS: ORBITAL S-M AXIS IN SOLAR RADII (6.96d5 KM),    
C     PERIOD IN DAYS, PHASE IN 2 PI RADIANS, SYSTEMIC VELOCITY AND         
C     THIRD LIGHT IN SAME UNITS AS VELOCITY AND LIGHT OBSERVATIONS,        
C     INCLINATION IN DEGREES, TEMPERATURES IN 10000K., SPOT LATITUDES      
C     IN radians (0=NORTH POLE, pi=SOUTH POLE), SPOT LONGITUDES      
C     IN radians (0=LINE OF CENTERS MERIDIAN, INCREASING COUNTER-          
C     CLOCKWISE AS SEEN FROM NORTH POLE TO 2*pi),                   
C     SPOT ANGULAR RADII IN radians. SPOT TEMPERATURE FACTOR IS            
C     DIMENSIONLESS.                                                       
c     
c  Integer MAGLITE is 1 for magnitude input and 0 for flux input.
c   
      program dc   
      implicit real*8 (a-h,o-z)                                            
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      ARRAY DIMENSIONING WRAPPER
c                            March 6, 2007
c
c     The following parameters are maximum array sizes.
c     There is no need to change any dimension numbers except these
c     in order to accomodate finer grids.
c
c        Nmax    ..    maximum grid fineness parameters N1 and N2
c                        default:   Nmax =    100
c      igsmax    ..    maximum number of surface elements (depends on N1 and N2)
c                        e.g. igsmax=762 for N=30, 3011 for N=60, etc.
c                        default: igsmax =   8331
c      ifoumax   ''    maximum number of points on the eclipsing horizon
c      ispmax    ..    maximum number of spots
c                        default: ispmax =    100
c      iclmax    ..    maximum number of clouds
c                        default: iclmax =    100
c      iptmax    ..    maximum number of observed data points, including
c                        blank points on last lines of the velocity and light curve
c                        data sets and on stop lines
c                        default: iptmax =  20000
c       ncmax    ..    maximum number of input data 'curves' (velocity +light + timings)
c                        default: 50
c       ipmax    ..    maximum number of parameters that are actually
c            adjusted, with band-independent parameters counted once each and
c            band-dependent parameters counted N_band times each.
c                        default: ipmax= 50
c
c       nbmax    ..    number of defined photometric bands
c
      parameter (Nmax=     100)
      parameter (igsmax=  8331)
      parameter (ifoumax=  5600)
      parameter (ispmax=   100)
      parameter (iclmax=   100)
      parameter (iptmax= 90000)
      parameter (ncmax=     50)
      parameter (ipmax=     50)
      parameter (nbmax=100)
      parameter (ntmax=80)
      parameter (ngmax=11)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Other array dimensions that are set automatically are listed
c     below and should not be changed, as the above parameter statements
c     determine their values.
c
c        MMmax    ..    maximum dimension of the MMSAVE array
c       ifrmax    ..    maximum dimension of the horizon polar coordinate arrays
c       istmax    ..    maximum dimension of storage arrays OBS and HOLD
c                         (iptmax * (no. of adjusted parameters + 1).
c
      parameter (MMmax=2*Nmax+4)
      parameter (ifrmax=4*Nmax)
      parameter (istmax=iptmax*ipmax)
      parameter (nplcofmax=50*nbmax)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Finally, the following dimensions are considered static and
c     their size does not depend on parameters.
c
c       ichno    ..    number of parameter channels (currently 60)
c
      parameter (ichno=60)
      character*5 arad
      dimension xtha(4),xfia(4),arad(4),po(2),omcr(2)
      dimension iband(50),xx(nbmax),
     $calib(50),refflux(50),xunit(50),aextinc(50),sphas1(50),
     $sphas2(50),sphas3(50),sphas4(50),extarray(nbmax),ffband(nbmax)
      dimension aa(20),bb(20),hotr(10),coolr(10)
c
c     Nothing needs to be changed beyond this point to accomodate finer grids.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      dimension rv(igsmax),grx(igsmax),gry(igsmax),grz(igsmax),
     $rvq(igsmax),grxq(igsmax),gryq(igsmax),grzq(igsmax),slump1(igsmax),
     $slump2(igsmax),srv(igsmax),sgrx(igsmax),sgry(igsmax),sgrz(igsmax),
     $srvq(igsmax),sgrxq(igsmax),sgryq(igsmax),sgrzq(igsmax),
     $srvl(igsmax),sgrxl(igsmax),sgryl(igsmax),sgrzl(igsmax),
     $srvql(igsmax),sgrxql(igsmax),sgryql(igsmax),sgrzql(igsmax),
     $slmp1(igsmax),slmp2(igsmax),slmp1l(igsmax),slmp2l(igsmax),
     $fr1(igsmax),fr2(igsmax),glump1(igsmax),glump2(igsmax),
     $grv1(igsmax),grv2(igsmax),xx1(igsmax),xx2(igsmax),yy1(igsmax),
     $yy2(igsmax),zz1(igsmax),zz2(igsmax),gmag1(igsmax),gmag2(igsmax),
     $csbt1(igsmax),csbt2(igsmax),rf1(igsmax),rf2(igsmax),
     $rftemp(igsmax),sxx1(igsmax),sxx2(igsmax),syy1(igsmax),
     $syy2(igsmax),szz1(igsmax),szz2(igsmax),sgmg1(igsmax),
     $sgmg2(igsmax),sgrv1(igsmax),sgrv2(igsmax),sglm1(igsmax),
     $sglm2(igsmax),scsb1(igsmax),scsb2(igsmax),srf1(igsmax),
     $srf2(igsmax),sglm1l(igsmax),sglm2l(igsmax),sgrv1l(igsmax),
     $sgrv2l(igsmax),sxx1l(igsmax),sxx2l(igsmax),syy1l(igsmax),
     $syy2l(igsmax),szz1l(igsmax),szz2l(igsmax),sgmg1l(igsmax),
     $sgmg2l(igsmax),scsb1l(igsmax),scsb2l(igsmax),srf1l(igsmax),
     $srf2l(igsmax),glog1(igsmax),glog2(igsmax),erv(igsmax),
     $egrx(igsmax),egry(igsmax),egrz(igsmax),elmp1(igsmax),
     $eglm1(igsmax),egrv1(igsmax),exx1(igsmax),eyy1(igsmax),
     $ezz1(igsmax),egmg1(igsmax),ecsb1(igsmax),erf1(igsmax),
     $ervq(igsmax),egrxq(igsmax),egryq(igsmax),egrzq(igsmax),
     $elmp2(igsmax),eglm2(igsmax),egrv2(igsmax),exx2(igsmax),
     $eyy2(igsmax),ezz2(igsmax),egmg2(igsmax),ecsb2(igsmax),
     $erf2(igsmax),ervl(igsmax),egrxl(igsmax),egryl(igsmax),
     $egrzl(igsmax),elmp1l(igsmax),eglm1l(igsmax),egrv1l(igsmax),
     $exx1l(igsmax),eyy1l(igsmax),ezz1l(igsmax),egmg1l(igsmax),
     $ecsb1l(igsmax),erf1l(igsmax),ervql(igsmax),egrxql(igsmax),
     $egryql(igsmax),egrzql(igsmax),elmp2l(igsmax),eglm2l(igsmax),
     $egrv2l(igsmax),exx2l(igsmax),eyy2l(igsmax),ezz2l(igsmax),
     $egmg2l(igsmax),ecsb2l(igsmax),erf2l(igsmax),sfr1(igsmax),
     $sfr1l(igsmax),efr1(igsmax),efr1l(igsmax),sfr2(igsmax),
     $sfr2l(igsmax),efr2(igsmax),efr2l(igsmax)
      dimension stldh(2*igsmax),stldl(2*igsmax),etldh(2*igsmax),
     $etldl(2*igsmax)
      dimension obs(istmax),hold(istmax),obser(ifoumax)
      dimension xcl(iclmax),ycl(iclmax),zcl(iclmax),rcl(iclmax),
     $op1(iclmax),fcl(iclmax),dens(iclmax),encl(iclmax),
     $edens(iclmax),xmue(iclmax)
      dimension xlat(2,ispmax),xlong(2,ispmax),radsp(2,ispmax),
     $temsp(2,ispmax),tstart(2,ispmax),tmax1(2,ispmax),tmax2(2,ispmax),
     $tfinal(2,ispmax)
      dimension mmsavh(MMmax),mmsavl(MMmax),mntype(iptmax)
      dimension phas(iptmax),flux(iptmax),wt(iptmax),br(iptmax),
     $bl(iptmax),phjd(iptmax),dfdph(iptmax),dfdap(iptmax),xmag(iptmax),
     $cgsta(iptmax),dfdrs(iptmax),timecomp(iptmax)
      dimension hla(ncmax),cla(ncmax),x1a(ncmax),x2a(ncmax),y1a(ncmax),
     $y2a(ncmax),el3a(ncmax),wla(ncmax),noise(ncmax),sigma(ncmax),
     $sigmainit(ncmax),opsfa(ncmax),knobs(ncmax+2),ksd(ncmax),
     $ratlum(ncmax)
      dimension snthh(2*Nmax),csthh(2*Nmax),snthl(2*Nmax),csthl(2*Nmax),
     $snfih(2*igsmax),csfih(2*igsmax),snfil(2*igsmax),csfil(2*igsmax)
      dimension hld(igsmax),tldh(2*igsmax),tldl(2*igsmax)
      dimension theta(ifrmax),rho(ifrmax)
      dimension del(ichno),keep(ichno+1),kep(ichno),nshift(ichno+1),
     $low(ichno),dlif(ichno)
      dimension clc(ipmax),out(ipmax),sd(ipmax),ccl(ipmax),ll(ipmax),
     $mm(ipmax),cnc(ipmax**2),cn(ipmax**2),cnn(ipmax**2)
      dimension para(30+5*ncmax),v(ipmax),cnout(ipmax**2),
     $parainit(30+5*ncmax)
c                                                                          
c The following dimensioned variables are not used by DC. They are         
c    dimensioned only for compatibility with usage of subroutine           
c    LIGHT by program LC.                                                  
c                                                                          
      dimension fbin1(1),fbin2(1),delv1(1),delv2(1),count1(1),             
     $count2(1),delwl1(1),delwl2(1),resf1(1),resf2(1),wl1(1),wl2(1),       
     $dvks1(1),dvks2(1),tau1(1),tau2(1),taug(1),hbarw1(1),hbarw2(1),       
     $emm1(1),emm2(1),emmg(1),yskp(1),zskp(1)   
      common /abscissa/ aba(29),ind0(10)
      common /klsp/ tratio,lsp
      common /xybol/ xbol1,xbol2,ybol1,ybol2
      common /abung/ abun(19),glog(ngmax)
      common /effwave/ effwvl(nbmax)
      common /planckleg/ plcof(nplcofmax)
      common /nummod/ nbp,ntemp,ngl
      common /atmmessages/ message(2,4),komp
      common /ramprange/ tlowtol,thightol,glowtol,ghightol   
      COMMON /ECCEN/ EC,aeb,pzero,VGA,SINI,VF,VFAC,VGAM,VOL1,VOL2,IFC       
      COMMON /DPDX/ DPDX1,DPDX2,POT1,POT2                                  
      COMMON /SUMM/ SUMM1,SUMM2                                            
      common /ardot/ dperdt,hjd,hjd0,perdum                                
      COMMON /FLVAR/ PSH,DP,EF,EFC,ECOS,pert,PHPER,pconsc,pconic,            
     $PHPERI,VSUM1,VSUM2,VRA1,VRA2,VKM1,VKM2,VUNIT,vfvu,trc,qfacd          
      COMMON /INVAR/ KH,IPB,IRTE,NREF,IRVOL1,IRVOL2,mref,ifsmv1,ifsmv2,    
     $icor1,icor2,ld1,ld2,ncl,jdphs,ipc,nr                                         
      COMMON /SPOTS/ SNLAT(2,ispmax),CSLAT(2,ispmax),SNLNG(2,ispmax),
     $CSLNG(2,ispmax),rdsp(2,ispmax),tmsp(2,ispmax),xlng(2,ispmax),
     $kks(2,ispmax),Lspot(2,ispmax)
      common /spev/ tstrt(2,ispmax),tmx1(2,ispmax),tmx2(2,ispmax),
     $tfnal(2,ispmax),amax(2,ispmax),radinc(2,ispmax),norscalc
      common /nspt/ nsp1,nsp2,kspot,kspev
      common /cld/ acm,opsf                                                
      common /ipro/ nbins,nl,inmax,inmin,nf1,nf2                           
      common /prof2/ duma,dumb,dumc,dumd,du1,du2,du3,du4,du5,du6,du7   
      common /inprof/ in1min,in1max,in2min,in2max,mpage,nl1,nl2   
      common /con3b/ xm03b,t03b,e3b,p3b,perr3b,vfac3b,ecos3b,if3b 
      common /poleint/ polin1,polin2
      data ARAD(1),ARAD(2),ARAD(3),ARAD(4)/' pole','point',' side',' bac
     $k'/   
c     
c           Bandpass Label Assignments for Stellar Atmospheres   
c    Label   Bandpass   Reference for Response Function   
c    -----   --------   -------------------------------   
c       1        u      Crawford, D.L. and Barnes, J.V. 1970, AJ, 75, 978   
c       2        v          "                "           "   
c       3        b          "                "           "   
c       4        y          "                "           "   
c       5        U      Buser, R. 1978, Ang, 62, 411   
c       6        B      Azusienis and Straizys 1969, Sov. Astron., 13, 316   
c       7        V          "             "                "   
c       8        R      Johnson, H.L. 1965, ApJ, 141, 923   
c       9        I         "            "    "   
c      10        J         "            "    "   
c      11        K         "            "    "   
c      12        L         "            "    "   
c      13        M         "            "    "   
c      14        N         "            "    "   
c      15        R_c    Bessell, M.S. 1983, PASP, 95, 480   
c      16        I_c       "            "    "   
c      17      230      Kallrath, J., Milone, E.F., Terrell, D., Young, A.T.   
c                          1998, ApJ, 508, 308   
c      18      250         "             "             "           "   
c      19      270         "             "             "           "   
c      20      290         "             "             "           "   
c      21      310         "             "             "           "   
c      22      330         "             "             "           "   
c      23     'TyB'    Tycho catalog B
c      24     'TyV'    Tycho catalog V
c      25     'HIP'    Hipparcos catalog
c   
c  **************************************************************************
c   SOME PUBLISHED ABSOLUTE CALIBRATIONS OF STANDARD PHOTOMETRIC SYSTEMS
c      (All for stars of magnitude 0.00 in the given system)
c
c   Bessell (1979) gives (converted to erg/(cm^2*sec*cm):
c     calib=4.19d-01  (U)
c     calib=6.60d-01  (B)
c     calib=3.61d-01  (V)
c     calib=2.25d-01  (R Cousins)
c     calib=1.22d-01  (I Cousins)
c     calib=4.01d-03  (K)
c
c   Johnson (1966) gives (again converted to erg/(cm^2*sec*cm):
c     calib=4.35d-01  (U)
c     calib=6.88d-01  (B)
c     calib=3.78d-01  (V)
c     calib=1.85d-01  (R)
c     calib=8.99d-02   (I)
c     calib=3.40d-02   (J)
c     calib=3.90d-03   (K)
c     calib=8.04d-04   (L)
c     calib=2.16d-04   (M)
c     calib=1.24d-05  (N)
c
c  ****************************************************************************
   67 FORMAT(20A4)                                                         
   21 FORMAT(d27.16,d36.16,d22.6)                                          
   55 FORMAT(11(3X,d8.1))                                                  
   56 FORMAT(11(1X,d7.1))                                                  
   20 format(1x,2(4i1,1x),7i1,1x,9(5i1,1x),i2,d10.3,f6.3)
    1 FORMAT(I3,I6,I6,I7,I7,I5,I5,I5,f15.6,d13.5,f10.5,f12.3,f14.4)        
  701 FORMAT(4I2,4I4,f13.6,d13.5,F8.5,F9.3)                                
  131 FORMAT(5(F14.5,f11.6,F8.3))                                           
  129 format(5(F14.5,i6,F13.3))
  130 format(5(F14.5,i6,F13.3))
  132 format(5(F14.5,i6,' *',F11.3))
  218 FORMAT(i3,2d13.6,4F7.3,d10.3,d12.5,4f8.5,f10.6,i2)
  318 FORMAT(i3,2d13.6,4F7.3,d10.3,d12.5,4f9.5,f10.6,i4)
   18 format(i3,2d13.6,4f7.3,d12.4,d10.3,i2,d12.5,4f8.5,i2)
   85 FORMAT(i3,2d13.6,4(1X,F6.3),d13.4,d10.3,i6,d13.5,4f9.5,i4)
   86 format(i3,f10.6,f10.4,d12.4,d13.5)
   99 format(f9.6,f8.4,d11.4,d12.5)
   37 FORMAT(1X,11F12.7)                                                   
  137 FORMAT(1X,F11.7)                                                     
  715 format(22x,'Input-Output in F Format')                               
  716 format(22x,'Input-Output in D Format')                               
   43 format('No.',2x,'Curve',5x,'Input Param.',8x,'Correction',6x,        
     $'Output Param.',5x,'Standard Error',5x,'L1/(L1+L2)',5x,'S.E.')                             
  615 format(i2,i7,2(f19.10,f18.10),2f13.5)                                                
  616 format(i2,i7,2(d19.10,d18.10),2x,2d13.5)                                                
  138 FORMAT(1X,'SUM OF ABSOLUTE VALUES OF CHECKS IS',1X,D12.6)            
  181 FORMAT(7X,'Normal Equations')                                        
  183 FORMAT (7X,'Correlation Coefficients')                               
  184 FORMAT(7X,'Normal Equations Times Inverse')                          
  185 FORMAT(1X,'Check of Computed Dependent Variables from Normal Equat   
     $ions')                                                               
   82 FORMAT(7X,'Unweighted Observational Equations')                      
   83 FORMAT(7X,'Weighted Observational Equations')                        
    9 FORMAT(33X,'Observations')                                           
  955 FORMAT(3(9x,'phase   V rad     wt'))                                   
   10 FORMAT(3(9x,'phase   light     wt'))                                   
  755 FORMAT(3(9x,'JD      V rad     wt'))                                   
  756 FORMAT(3(9x,'JD      light     wt'))                                   
  110 FORMAT(3(9x,'phase     mag       wt  '))
  757 FORMAT(3(9x,'JD        mag       wt  '))
  758 FORMAT(3(4x,'min. time    type',9x,'wt  '))
   40 FORMAT('   Mean residual for input values       Mean residual pred   
     $icted          determinant')                                         
   11 format('band',6x,'L1',11x,'L2',8x,'x1     x2     y1     y2',9x, 
     $'3rdlt',4x,'opsf    NOISE   Sigma',6x,'phase 1',2x,'phase 2',2x,
     $'phase 3',2x,'phase 4',2x,'KSD')
   88 format('band',3x,'Wave L',4x,'extinc',6x,'xunit',6x,'calib')
  111 FORMAT('band',6x,'L1',11x,'L2',8x,'x1',5x,'x2',5x,'y1',5x,'y2',   
     $5x,'opsf',7x,'Sigma',4x,'phase 1',2x,'phase 2',2x,'phase 3',2x,'ph
     $ase 4','   Wave L','   KSD')                                                       
   12 FORMAT('MODE   IPB  IFAT1  IFAT2   N1   N2  N1L  N2L',7x,'Arg Per'   
     $,5x,'dperdt',7x,'TH e',4x,'V unit(km/s)     V FAC')                     
  206 FORMAT(F6.5,d13.6,2F11.4,F11.4,f10.3,2f8.3,f7.2,2f11.4,i6,i7,i7)                    
  205 FORMAT('  ecc',5x,'s-m axis',7x,'F1',9x,'F2',7x,'Vgam',7x,'Incl',          
     $6x,'g1',6x,'g2',4x,'[M/H]',5x,'Fspot1',5x,'Fspot2',2x,
     $'Nspot1 Nspot 2  iab')
  402 FORMAT(6x,'del a      del e     del PER    del F1     del F2     d      
     $el PHS     del i    del g1     del g2     del T1     del T2')               
  403 FORMAT('    del ALB1   del ALB2   del POT1   del POT2   del Q        
     $ del L1     del L2     del X1     del X2')                           
  406 FORMAT(' Adjustment control integers; 1 suppresses adjustment, 0 a   
     $llows adjustment.')                                                  
  702 FORMAT(F6.5,d13.6,2F10.4,F10.4,f9.3,2f7.3,f7.2,2f10.4)                           
  706 FORMAT(F7.4,f8.4,2F7.3,3d13.6,4f7.3,f9.5)                           
  408 FORMAT(2F8.4,2F9.3,2d15.6,d13.6,4f9.3,f12.5)                               
  705 FORMAT(I1,1X,I1,1X,9I2)                                              
   23 format('     T1      T2',5x,'Alb 1',4x,'Alb 2',9x,'Pot 1',10x,       
     $'Pot 2',8x,'M2/M1  x1(bolo) x2(bolo) y1(bolo) y2(bolo)')
   54 FORMAT('     T1      T2',5x,'Alb 1',4x,'Alb 2',9x,'Pot 1',10x,       
     $'Pot 2',8x,'M2/M1  x1(bolo) x2(bolo) y1(bolo) y2(bolo)',2x,
     $'log(d) [pc]')            
  707 FORMAT('    IFVC1   IFVC2   NLC IFTIMES  KO   KDISK   ISYM   NPPL   
     $  IFDER  IFLCIN  IFOC')         
  917 format('nref',3x,'mref',3x,'ifsmv1',3x,'ifsmv2',3x,'icor1',3x,       
     $'icor2',3x,'if3b',4x,'ld1',5x,'ld2',4x,'kspev',2x,'kspot',2x,
     $'nomax',2x,'ifcgs',2x,'maglite',2x,'linkext',2x,'desig. ext.')
  708 FORMAT(11(4X,I3))                                                     
  912 format(i3,i7,i8,i9,i8,i8,i8,i8,5i7,i10,i8,f13.4)                                         
  650 FORMAT(20X,'radii and related quantities (from iteration output)')              
  651 FORMAT(5X,'dOm1/dq',5X,'dOm2/dq',5X,'Om1-q corr.',5X,'Om2-q corr.'   
     $,5X,'Om1 s.e.',4X,'Om2 s.e.',4X,'q  s.e.')                           
  652 FORMAT(1X,3F12.6,4X,F12.6,4X,3F12.6)                                 
  653 FORMAT(' component',7X,'r=R/a',9X,'dr/dOm',8X,'dr/dq',11X,'s.e.')       
  654 FORMAT(2X,I1,1X,A6,4F14.6)                                           
  684 FORMAT(i2,4F13.5,4f14.5)                                                    
  985 FORMAT(4f9.5,4f14.5)                                                        
  983 FORMAT(1X,'star  co-latitude  longitude  spot radius temp factor',   
     $2x,'onset time',5x,'begin max',5x,'end max',5x,'end time')
  399 FORMAT('    del lat    del long   del rad    del tempf  del lat      
     $ del long   del rad    del tempf')                                   
   60 FORMAT(4I3)                                                          
   61 FORMAT(1X,4I6)                                                       
   66 FORMAT('   STAR  SPOT   STAR  SPOT')                                 
  166 FORMAT(' SPOTS TO BE ADJUSTED')                                      
  440 FORMAT(' AS1=FIRST ADJUSTED SPOT')                                   
  441 FORMAT(' AS2=SECOND ADJUSTED SPOT')                                  
  405 FORMAT(' ORDER OF PARAMETERS IS AS FOLLOWS:')                        
 1440 FORMAT('  (1) - AS1 LATITUDE')                                       
 1441 FORMAT('  (2) - AS1 LONGITUDE')                                      
 1442 FORMAT('  (3) - AS1 ANGULAR RADIUS')                                 
 1443 FORMAT('  (4) - AS1 TEMPERATURE FACTOR')                             
 1444 FORMAT('  (5) - AS2 LATITUDE')                                       
 1445 FORMAT('  (6) - AS2 LONGITUDE')                                      
 1446 FORMAT('  (7) - AS2 ANGULAR RADIUS')                                 
 1447 FORMAT('  (8) - AS2 TEMPERATURE FACTOR')                             
 1448 FORMAT('  (9) - A=ORBITAL SEMI-MAJOR AXIS in Solar radii')                          
 1449 FORMAT(' (10) - E=ORBITAL ECCENTRICITY')                             
 1450 FORMAT(' (11) - PERR0=ARGUMENT of PERIASTRON at time HJD0 in radi 
     $ans')          
 1451 FORMAT(' (12) - F1=STAR 1 ROTATION PARAMETER')                       
 1452 FORMAT(' (13) - F2=STAR 2 ROTATION PARAMETER')                       
 1453 FORMAT(' (14) - PHASE SHIFT= PHASE OF PRIMARY CONJUNCTION')          
 1454 FORMAT(' (15) - VGAM=SYSTEMIC RADIAL VELOCITY')                      
 1455 FORMAT(' (16) - INCL=orbital inclination in degrees')                           
 1456 FORMAT(' (17) - g1=STAR 1 GRAVITY DARKENING EXPONENT')               
 1457 FORMAT(' (18) - g2=STAR 2 GRAVITY DARKENING EXPONENT')               
 1458 FORMAT(' (19) - T1=STAR 1 AVERAGE SURFACE TEMPERATURE in 10000 K')              
 1459 FORMAT(' (20) - T2=STAR 2 AVERAGE SURFACE TEMPERATURE in 10000 K')              
 1460 FORMAT(' (21) - ALB1=STAR 1 BOLOMETRIC ALBEDO')                      
 1461 FORMAT(' (22) - ALB2=STAR 2 BOLOMETRIC ALBEDO')                      
 1462 FORMAT(' (23) - POT1=STAR 1 modified surface potential')                      
 1463 FORMAT(' (24) - POT2=STAR 2 modified surface potential')                      
 1464 FORMAT(' (25) - Q=MASS RATIO (STAR 2/STAR 1)')                       
 1470 FORMAT(' (26) - HJD0= Hel. JD reference time')                       
 1471 FORMAT(' (27) - PZERO= orbital period in days at time HJD0')                             
 1472 FORMAT(' (28) - DPDT= dimensionless time derivative of orbital pe 
     $riod (e.g. days/day)')            
 1473 FORMAT(' (29) - DPERDT= time derivative of argument of periastron    
     $in radians per day')                                                                    
 1474 FORMAT(' (30) - a3b=relative semi-major axis for 3rd body orbit in 
     $ Solar radii')       
 1465 FORMAT(' (31) - P3b=3rd body orbit period in days')        
 1466 FORMAT(' (32) - i3b=inclination of 3rd body orbit to plane of sky  
     $in degrees')        
 1467 FORMAT(' (33) - e3b=eccentriicity of 3rd body orbit')               
 1468 FORMAT(' (34) - omega3b=argument of periastron of 3rd body orbit i 
     $n radians')               
 1469 FORMAT(' (35) - T03b=HJD of superior conjunction of inner system c 
     $enter of mass with regard to the outer orbit')                                    
 1475 FORMAT(' (41) - log d= log (base 10) of distance in parsecs')        
 1476 FORMAT(' (42) - A_desig = interstellar extinction in designated ph
     $otometric band')        
 1484 format(' (43) - onset time of spot 1')
 1485 format(' (44) - start time of constant maximum for spot 1')
 1486 format(' (45) - end time of constant maximum for spot 1')
 1487 format(' (46) - end time of spot 1')
 1495 format(' (47) - onset time of spot 2')
 1489 format(' (48) - start time of constant maximum for spot 2')
 1490 format(' (49) - end time of constant maximum for spot 2')
 1496 format(' (50) - end time of spot 2')
 1479 FORMAT(' (56) - L1=STAR 1 RELATIVE MONOCHROMATIC LUMINOSITY')        
 1480 FORMAT(' (57) - L2=STAR 2 RELATIVE MONOCHROMATIC LUMINOSITY')        
 1481 FORMAT(' (58) - X1=STAR 1 LIMB DARKENING COEFFICIENT')               
 1482 FORMAT(' (59) - X2=STAR 2 LIMB DARKENING COEFFICIENT')               
 1483 FORMAT(' (60) - el3=third light')                                    
  119 format(1x,i6,i13,d22.8,i13,d26.8,d19.4,f16.4)                                              
  159 format(' Standard Deviations for Computation of Curve-dependent We
     $ights')   
  169 format('    Curve   No. of obs.        Std. dev.',7x,'No. of obs
     $',2x,'Std. dev. in phase range',3x,'reference flux',3x,'reference 
     $mag')
 1063 format(3f9.4,f7.4,d11.4,f9.4,d11.3,f9.4,f7.3)                        
   64 format(3f10.4,f9.4,d12.4,f10.4,d12.4,f9.4,f9.3,d12.4)                
   69 format('      xcl       ycl       zcl      rcl       op1         f   
     $cl        ne       mu e      encl     dens')                         
  170 format(i3,f17.6,d18.10,d14.6,f10.4,f11.5,i5)                                  
  649 format(i1,f15.6,d17.10,d14.6,f10.4,f8.5,i3)                                  
  171 format('JDPHS',5x,'J.D. zero',8x,'P zero',11x,'dPdt',                
     $6x,'Ph. shift',5x,'delph',4x,'NGA')
  911 format(7(i1,1x),2(i2,1x),6(i1,1x),f7.4)                                                     
  840 format('Do not try to adjust the ephemeris or any time derivative    
     $parameters','when JDPHS = 2')                                        
  839 format('Ordinarily one should not try to adjust both PSHIFT and',    
     $' HJD0. They are perfectly correlated if the period is constant',    
     $' and extremely highly correlated if the period is not constant')    
  283 format('log g below ramp range for at least one point',  
     $' on star',i2,', black body applied locally.')  
  284 format('log g above ramp range for at least one point',  
     $' on star',i2,', black body applied locally.')  
  285 format('T above ramp range for at least one',  
     $' point on star',i2,', black body applied locally.')  
  286 format('T below ramp range for at least one point',  
     $' on star',i2,', black body applied locally.')  
  287 format('Input [M/H] = ',f6.3,' is not a value recognized by ',  
     $'the program. Replaced by ',f5.2)  
  101 format(d12.6,d14.7,f11.5,f9.6,f10.7,f17.8) 
  102 format(d12.6,d14.7,f11.5,f9.6,f13.7,f17.8,d15.7) 
  103 format(5x,'a3b',9x,'P 3b',8x,'incl 3b',5x,'e 3b',3x, 
     $'arg. perr. 3b',4x,'T conj 3b',6x,'m3/(m1+m2)') 
  104 format('Computed third star mass is negative. Input parameters  
     $unrealistic.') 
  105 format('With assumption of a dynamical third body (control 
     $integer if3b=1), the program works only with time as the 
     $independent variable (JDPHS=1). Program $stopping.')
  149 format(20f20.10)
  147 format(f13.5,22f20.10)
  356 format('The distance parameter cannot be adjusted with IFCGS=0')
  354 format('Interstellar extinction cannot be adjusted with IFCGS=0')
  182 format(5x,'time',10x,'phase',7x,'observed',7x,'computed')
  175 format(5x,'time',10x,'phase')
  192 format(5x,'phase',7x,'observed',7x,'computed')
  176 format(5x,'phase')
  829 format('Third body parameters (numbers 30 to 35) cannot be adjuste
     $d with if3b=0.')
  830 format('Their adjustment has been turned off.')
    2 FORMAT(F6.5,d13.6,2F10.4,F10.4,f9.3,2f7.3,f7.2,2f10.4)
    6 FORMAT(F7.4,1X,f7.4,2f7.3,3d13.6,4F7.3,f8.5)
  568 FORMAT(4I2,2I4,f13.6,d14.6,f8.5,F8.2)
    4 FORMAT(i3,2d13.6,4F7.3,d12.4,d11.4,F8.3,F8.4,f10.6,f8.4,d12.5)
  567 FORMAT(4f9.5,4f14.5)
   22 format(8(i1,1x),2(i2,1x),4(i1,1x))
  811 format(i1,f15.6,d17.10,d14.6,f10.4,f8.5,i3,d11.4,i2,f11.0)
  217 format(f14.6,f15.6,f13.6,4f12.6,f10.4,i2,f8.4)
  566 format(d12.6,d14.7,f11.5,f9.6,f10.7,f17.8)
  210 format(f4.0)
  211 format(i1)
 3333 format('The iteration above is number',i4,' of',i4)
  832 format('Note: tmax2 is set equal to tmax1 and cannot be adjusted i
     $f nomax=1. Adjustment control integers 45 and 49 have been re-set 
     $accordingly.')
  212 format('The largest NGA allowed in this program version is 10.')
 6066 format('nbpsav, nbp mismatch, program stopped.',2i5)
   77 format(82x,d11.5,32x,i2)
   76 format(88x,d11.5,38x,i2)
  820 format('NGA cannot be less than 1. Reset to 1')
  139 format(F14.5,i6,6x)
  216 format(f9.0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
c  Ramp ranges are set here. The values below seem to work well, however,   
c  they can be changed.   
      tlowtol=1500.d0   
      thightol=50000.d0   
      glowtol=4.0d0   
      ghightol=4.0d0   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
      abun(1)=1.d0 
      abun(2)=0.5d0 
      abun(3)=0.3d0 
      abun(4)=0.2d0 
      abun(5)=0.1d0 
      abun(6)=0.0d0 
      abun(7)=-0.1d0 
      abun(8)=-0.2d0 
      abun(9)=-0.3d0 
      abun(10)=-0.5d0 
      abun(11)=-1.0d0 
      abun(12)=-1.5d0 
      abun(13)=-2.0d0 
      abun(14)=-2.5d0 
      abun(15)=-3.0d0 
      abun(16)=-3.5d0 
      abun(17)=-4.0d0 
      abun(18)=-4.5d0 
      abun(19)=-5.0d0 
      glog(1)=0.0d0 
      glog(2)=0.5d0 
      glog(3)=1.0d0 
      glog(4)=1.5d0 
      glog(5)=2.0d0 
      glog(6)=2.5d0 
      glog(7)=3.0d0 
      glog(8)=3.5d0 
      glog(9)=4.0d0 
      glog(10)=4.5d0 
      glog(11)=5.0d0 
      message(1,1)=0  
      message(1,2)=0  
      message(2,1)=0  
      message(2,2)=0  
      message(1,3)=0  
      message(1,4)=0  
      message(2,3)=0  
      message(2,4)=0  
      open(unit=24,file='effwvl.dat',status='old')
      do 2928 j=1,nbmax
      read(24,*,end=2929) jdum,effwvl(j)
 2928 continue
 2929 nbpsav=j-1
      close (24)
      nplcof=50*nbpsav
      open(unit=23,file='atmcofplanck.dat',status='old')   
      read(23,*) (plcof(j),j=1,nplcof)
      close(23)
      open(unit=5,file='dcin.active',status='old')  
      open(unit=6,file='dcout.active')  
c
c  lsp is used by LC but not needed in DC. IT is assigned a value to 
c    avoid triggering diagnostic warnings.
c
      lsp=0
      tt=2.d0/3.d0                                                         
      nbins=1                                                              
      nl=1                                                                 
      inmax=1                                                              
      inmin=1                                                              
      nf1=1                                                                
      nf2=1                                                                
      toldis=1.d-5                                                         
      enap=dexp(1.d0)
      xlog10e=dlog10(enap)
      pi=dacos(-1.d0)                                                      
      xlg10=dlog(10.d0)
      twopi=pi+pi 
      pih=.5d0*pi 
      dtr=pi/180.d0 
      threepi=twopi+pi 
      en0=6.02214199d23                                                        
      clight=2.99792458d5 
      rsunkm=6.95660d5 
      rsuncm=rsunkm*1.d5
      aupc=6.48d5/pi
      aucm=1.49597870700d13
      pcrsun=rsuncm/(aucm*aupc)
      dtcon=rsunkm/(8.64d4*clight) 
      XTHA(1)=0.d0                                                         
      XTHA(2)=.5d0*PI                                                      
      XTHA(3)=.5d0*PI                                                      
      XTHA(4)=.5d0*PI                                                      
      XFIA(1)=0.d0                                                         
      XFIA(2)=0.d0                                                         
      XFIA(3)=.5d0*PI                                                      
      XFIA(4)=PI                                                           
      niter=0
c   
c  The initializations in the 886 and 887 loops are just to avoid   
c    triggering error messages from some compilers. The quantities do   
c    not otherwise need initialization. Same for du1, ...., du7.   
c    Same for mpage, nl1, nl2.   
c   
      key=0
      do 886 ikks=1,2   
      do 886 jkks=1,ispmax   
  886 kks(ikks,jkks)=0   
      do 887 immsav=1,MMmax   
      mmsavh(immsav)=0   
  887 mmsavl(immsav)=0   
      du1=0.d0   
      du2=0.d0   
      du3=0.d0   
      du4=0.d0   
      du5=0.d0   
      du6=0.d0   
      du7=0.d0   
      mpage=0   
      nl1=0   
      nl2=0   
      IBEF=0                                                               
      KH=25                                                                
      NS=1                                                                 
      NI=0                                                                 
      NY=0                                                                 
      KNOBS(1)=0                                                           
      WRITE(6,405)                                                         
      write(6,*) 
      WRITE(6,440)                                                         
      WRITE(6,441)                                                         
      write(6,*) 
      WRITE(6,1440)                                                        
      WRITE(6,1441)                                                        
      WRITE(6,1442)                                                        
      WRITE(6,1443)                                                        
      WRITE(6,1444)                                                        
      WRITE(6,1445)                                                        
      WRITE(6,1446)                                                        
      WRITE(6,1447)                                                        
      WRITE(6,1448)                                                        
      WRITE(6,1449)                                                        
      WRITE(6,1450)                                                        
      WRITE(6,1451)                                                        
      WRITE(6,1452)                                                        
      WRITE(6,1453)                                                        
      WRITE(6,1454)                                                        
      WRITE(6,1455)                                                        
      WRITE(6,1456)                                                        
      WRITE(6,1457)                                                        
      WRITE(6,1458)                                                        
      WRITE(6,1459)                                                        
      WRITE(6,1460)                                                        
      WRITE(6,1461)                                                        
      WRITE(6,1462)                                                        
      WRITE(6,1463)                                                        
      WRITE(6,1464)                                                        
      WRITE(6,1470)                                                        
      WRITE(6,1471)                                                        
      WRITE(6,1472)                                                        
      WRITE(6,1473)                                                        
      WRITE(6,1474)                                                        
      WRITE(6,1465)                                                        
      WRITE(6,1466)                                                        
      WRITE(6,1467)                                                        
      WRITE(6,1468)                                                        
      WRITE(6,1469)                                                        
      write(6,1475)
      write(6,1476)
      write(6,1484)
      write(6,1485)
      write(6,1486)
      write(6,1487)
      write(6,1495)
      write(6,1489)
      write(6,1490)
      write(6,1496)
      WRITE(6,1479)                                                        
      WRITE(6,1480)                                                        
      WRITE(6,1481)                                                        
      WRITE(6,1482)                                                        
      WRITE(6,1483)                                                        
      READ(5,56)(DEL(I),I=1,8)                                             
      READ(5,56)(DEL(I),I=9,14),(DEL(I),I=16,20)                          
      READ(5,56)(DEL(I),I=21,25),(del(i),i=56,59) 
      READ(5,20)(KEP(I),I=1,60),niters,xlamda,vlr
      READ(5,60) KSPA,NSPA,KSPB,NSPB                                       
      READ(5,705)IFVC1,IFVC2,NLC,iftime,KO,KDISK,ISYM,nppl,ifder,iflcin,
     $ifoc                        
      read(5,911)nref,mref,ifsmv1,ifsmv2,icor1,icor2,if3b,ld1,ld2,kspev,
     $kspot,nomax,ifcgs,maglite,linkext,desextinc         
      if(kep(45)*kep(49).eq.1.or.nomax.eq.0) goto 833
      kep(45)=1
      kep(49)=1
      write(6,*)
      write(6,832)
  833 continue
      ld1abs=iabs(ld1)
      ld2abs=iabs(ld2)
      ifadj3b=1-kep(30)*kep(31)*kep(32)*kep(33)*kep(34)*kep(35)
      if(ifadj3b.eq.0.or.if3b.eq.1) goto 827
      write(6,*)
      write(6,829)
      write(6,830)
      write(6,*)
      do 828 ikep=30,35
      kep(ikep)=1
  828 continue
  827 continue
      read(5,649) jdphs,hjd0,pzero,dpdt,pshift,delph,nga
      if(nga.gt.10) write(6,212)
      if(nga.gt.0) goto 821
      write(6,820)
      nga=1
  821 continue
      nhalf=nga/2
      emmm=.5d0*delph
      if(jdphs.eq.2.and.kep(26).eq.0) write(6,840)                         
      if(jdphs.eq.2.and.kep(27).eq.0) write(6,840)                         
      if(jdphs.eq.2.and.kep(28).eq.0) write(6,840)                         
      if(kep(14).eq.0.and.kep(26).eq.0) write(6,839)                       
      READ(5,701)MODE,IPB,IFAT1,IFAT2,N1,N2,N1L,N2L,perr0,dperdt,THE,      
     $vunit                                                                
      READ(5,702)E,A,F1,F2,VGA,XINCL,GR1,GR2,abunin,fspot1,fspot2   
      READ(5,706) TAVH,TAVC,ALB1,ALB2,poth,potc,RM,xbol1,
     $xbol2,ybol1,ybol2,dpclog      
      phsv=poth
      pcsv=potc
      read(5,101) a3b,p3b,xinc3b,e3b,perr3b,t03b 
      if(nga.le.1) goto 228
      do 227 iff=1,nga
      hotr(iff)=0.d0
  227 coolr(iff)=0.d0
      call gaussquad(1,nga,emmm,hotr,outt)
  228 continue
c 
c  Compute constants related to third star if necessary 
c 
      if(if3b.eq.0) goto 870 
      if(jdphs.eq.1) goto 871 
      write(6,*) 
      write(6,105) 
      write(6,*) 
      stop 
  871 continue 
      rvfac=rsunkm/(8.64d4*vunit) 
      efac3b=(1.d0+e3b)/(1.d0-e3b) 
      rtefac3b=dsqrt(efac3b) 
      twopioppr=twopi/p3b 
      arat=a/a3b 
      prat=p3b/pzero 
      ertfac=dsqrt(1.d0-e3b**2) 
      rm3b=1.d0/(arat**3*prat**2)-1.d0 
      if(rm3b.lt.0.d0) write(6,104) 
      xinc3brad=dtr*xinc3b 
      sini3b=dsin(xinc3brad) 
      cosi3b=dcos(xinc3brad) 
      tr3b0=pih-perr3b 
c  Next line is a default in case tr3b0/2 is too close to pi/2 
      ecan3b0=pi 
      parecande3b0=0.d0
      parecandtr3b0=rtefac3b
      trtest1=dabs(tr3b0-pi) 
      trtest2=dabs(tr3b0-threepi) 
c
c  Next line is a default in case tr3b0/2 is too close to 3/2 pi 
c
c     if(trtest2.lt.1.d-8) ecan3b0=pi32 
      if(trtest2.lt.1.d-9) ecan3b0=threepi 
      trtest3=dabs(tr3b0+pi) 
      trtest4=dabs(tr3b0+threepi) 
      if(trtest1.lt.1.d-9.or.trtest2.lt.1.d-9) goto 767 
      if(trtest3.lt.1.d-9.or.trtest4.lt.1.d-9) goto 767 
      tanhtr3b0=dtan(0.5d0*tr3b0)
      tanhtr3b0sq=tanhtr3b0**2
      ecan3b0=2.d0*datan(tanhtr3b0/rtefac3b) 
      parecande3b0=-2.d0*rtefac3b*tanhtr3b0/(1.d0+e3b)**2/
     $(1.d0+tanhtr3b0sq/efac3b)
      parecandtr3b0=(1.d0+tanhtr3b0sq)/(1.d0+
     $tanhtr3b0sq/efac3b)
      parecandtr3b0=parecandtr3b0/rtefac3b
  767 continue 
      xm03b=ecan3b0-e3b*dsin(ecan3b0) 
      dxmde3b0=(1.d0-e3b*dcos(ecan3b0))*parecande3b0
     $-dsin(ecan3b0)
      apfac=arat**3*prat**2 
      apfacm=1.d0-apfac 
      sinperr3b=dsin(perr3b) 
      cosperr3b=dcos(perr3b) 
      ecos3b=e3b*cosperr3b 
      vfac3b=twopioppr*rvfac*a3b*sini3b*apfacm/ertfac 
  870 continue 
      acm=rsuncm*a                                                       
      nn1=n1   
      CALL SINCOS(1,N1,N1,SNTHH,CSTHH,SNFIH,CSFIH,MMSAVH)                  
      CALL SINCOS(2,N2,N1,SNTHH,CSTHH,SNFIH,CSFIH,MMSAVH)                  
      CALL SINCOS(1,N1L,N1L,SNTHL,CSTHL,SNFIL,CSFIL,MMSAVL)                
      CALL SINCOS(2,N2L,N1L,SNTHL,CSTHL,SNFIL,CSFIL,MMSAVL)                
      dint1=pi*(1.d0-xbol1/3.d0)                                           
      if(ld1abs.eq.2) dint1=dint1+pi*ybol1*2.d0/9.d0                           
      if(ld1abs.eq.3) dint1=dint1-pi*ybol1*.2d0                                
      dint2=pi*(1.d0-xbol2/3.d0)                                           
      if(ld2abs.eq.2) dint2=dint2+pi*ybol2*2.d0/9.d0                           
      if(ld2abs.eq.3) dint2=dint2-pi*ybol2*.2d0                                
      IS=ISYM+1                                                            
      KEEP(61)=0                                                           
      MM1=N1+1                                                             
      MM2=N1+N2+2                                                          
      MM3=N1L+1                                                            
      MM4=N1L+N2L+2                                                        
      M1H=MMSAVH(MM1)                                                      
      M2H=MMSAVH(MM2)                                                      
      M1L=MMSAVL(MM3)                                                      
      M2L=MMSAVL(MM4)                                                      
      MTLH=M1H+M2H                                                         
      MTLL=M1L+M2L                                                         
      NVC=IFVC1+IFVC2                                                      
      NLVC=NLC+NVC                                                         
      nlvct=nlvc+iftime
      NVCP=NVC+1                                                           
      IF(NVC.NE.0) GOTO 288                                                
      KEP(9)=1                                                             
      KEP(15)=1                                                            
  288 CONTINUE                                                             
      DO 84 I=1,ichno                                                         
      KEEP(I)=KEP(I)                                                       
   84 LOW(I)=1                                                             
      LOW(1)=0                                                             
      LOW(2)=0                                                             
      LOW(3)=0                                                             
      LOW(5)=0                                                             
      LOW(6)=0                                                             
      LOW(7)=0                                                             
      LOW(9)=0                                                             
      LOW(10)=0                                                            
      LOW(11)=0                                                            
      LOW(12)=0                                                            
      LOW(13)=0                                                            
      LOW(14)=0                                                            
      LOW(16)=0                                                            
      LOW(23)=0                                                            
      LOW(24)=0                                                            
      LOW(25)=0                                                            
      ifap=1-keep(29)                                                      
      ifphi=1-keep(26)*keep(27)*keep(28)                                   
      KOSQ=(KO-2)*(KO-2)                                                   
      do 1135 i35=1,50
 1135 aextinc(i35)=0.d0
      IF(NVC.EQ.0) GOTO 195                                                
      DO 90 I=1,NVC                                                        
   90 READ(5,218) iband(i),HLA(I),CLA(I),X1A(I),X2A(I),y1a(i),y2a(i),        
     $opsfa(i),sigma(i),sphas1(i),sphas2(i),sphas3(i),sphas4(i),wla(i),
     $ksd(i)                                           
  195 CONTINUE                                                             
      IF(NLC.EQ.0) GOTO 194                                             
      DO 190 I=nvcp,nlvc                                                   
  190 read(5,18) iband(i),hla(i),cla(i),x1a(i),x2a(i),y1a(i),y2a(i),  
     $el3a(i),opsfa(i),noise(i),sigma(i),sphas1(i),sphas2(i),
     $sphas3(i),sphas4(i),ksd(i)
  194 if(iftime.eq.1) read(5,77) sigma(nlvct),ksd(nlvct)
      if(linkext.ne.0)call bandext(desextinc,linkext,xx,ffband,extarray)
      if(nlc.eq.0) goto 172
      do 112 i=nvcp,nlvc
      read(5,99) wla(i),aextinc(i),xunit(i),calib(i)
      if(linkext.eq.0) goto 112
      iban=iband(i)
      aextinc(i)=extarray(iban)
  112 continue
  172 continue
      NSP1=0                                                               
      NSP2=0                                                               
      DO 988 KP=1,2                                                        
      DO 987 I=1,ispmax                                                  
      READ(5,985)XLAT(KP,I),XLONG(KP,I),RADSP(KP,I),TEMSP(KP,I),
     $tstart(kp,i),tmax1(kp,i),tmax2(kp,i),tfinal(kp,i)
      if(nomax.ne.0) tmax2(kp,i)=tmax1(kp,i)
      xlng(kp,i)=xlong(kp,i)                                               
      IF(XLAT(KP,I).GE.200.d0) GOTO 988                                    
      SNLAT(KP,I)=dsin(XLAT(KP,I))                                         
      CSLAT(KP,I)=dcos(XLAT(KP,I))                                         
      SNLNG(KP,I)=dsin(XLONG(KP,I))                                        
      CSLNG(KP,I)=dcos(XLONG(KP,I))                                        
      RDSP(KP,I)=RADSP(KP,I)                                               
      amax(kp,i)=twopi*(1.d0-dcos(rdsp(kp,i)))
      TMSP(KP,I)=TEMSP(KP,I)                                               
      tstrt(kp,i)=tstart(kp,i)
      tmx1(kp,i)=tmax1(kp,i)
      tmx2(kp,i)=tmax2(kp,i)
      tfnal(kp,i)=tfinal(kp,i)
      IF(KP.EQ.1) NSP1=NSP1+1                                              
  987 IF(KP.EQ.2) NSP2=NSP2+1                                              
  988 CONTINUE                                                             
      NSTOT=NSP1+NSP2                                                      
      ncl=0                                                                
      do 1062 i=1,iclmax                                                      
      read(5,1063) xcl(i),ycl(i),zcl(i),rcl(i),op1(i),fcl(i),edens(i),     
     $xmue(i),encl(i)                                                      
      if(xcl(i).gt.100.d0) goto 1066                                       
      ncl=ncl+1                                                            
      dens(i)=edens(i)*xmue(i)/en0                                         
 1062 continue                                                             
 1066 continue                                                             
      para(1)=xlat(kspa,nspa)                                              
      para(2)=xlong(kspa,nspa)                                             
      para(3)=radsp(kspa,nspa)                                             
      para(4)=temsp(kspa,nspa)                                             
      para(5)=xlat(kspb,nspb)                                              
      para(6)=xlong(kspb,nspb)                                             
      para(7)=radsp(kspb,nspb)                                             
      para(8)=temsp(kspb,nspb)                                             
      para(9)=a                                                            
      para(10)=e                                                           
      para(11)=perr0                                                       
      para(12)=f1                                                          
      para(13)=f2                                                          
      para(14)=pshift                                                      
      para(15)=vga                                                         
      para(16)=xincl                                                       
      para(17)=gr1                                                         
      para(18)=gr2                                                         
      para(19)=tavh                                                        
      para(20)=tavc                                                        
      para(21)=alb1                                                        
      para(22)=alb2                                                        
      para(23)=phsv                                                        
      para(24)=pcsv                                                        
      para(25)=rm                                                          
      para(26)=hjd0                                                        
      para(27)=pzero                                                      
      para(28)=dpdt                                                        
      para(29)=dperdt                                                      
      if(ifadj3b.eq.0) goto 878 
      para(30)=a3b                                                        
      para(31)=p3b                                                        
      para(32)=xinc3b                                                        
      xinc3brad=dtr*xinc3b
      sini3b=dsin(xinc3brad)
      cosi3b=dcos(xinc3brad)
      para(33)=e3b                                                        
      para(34)=perr3b                                                        
      para(35)=t03b                                                        
  878 continue 
      para(41)=dpclog
      para(42)=desextinc
      para(43)=tstart(kspa,nspa)
      para(44)=tmax1(kspa,nspa)
      para(45)=tmax2(kspa,nspa)
      para(46)=tfinal(kspa,nspa)
      para(47)=tstart(kspb,nspb)
      para(48)=tmax1(kspb,nspb)
      para(49)=tmax2(kspb,nspb)
      para(50)=tfinal(kspb,nspb)
      do 341 ipara=1,55
  341 parainit(ipara)=para(ipara)
      ib=nvc                                                               
      do 191 irx=56,55+nlc                                                 
      ib=ib+1                                                              
      parainit(irx)=hla(ib)
  191 para(irx)=hla(ib)                                                    
      ib=nvc                                                               
      do 186 irx=56+nlc,55+2*nlc                                           
      ib=ib+1                                                              
      parainit(irx)=cla(ib)
  186 para(irx)=cla(ib)                                                    
      ib=nvc                                                               
      do 187 irx=56+2*nlc,55+3*nlc                                         
      ib=ib+1                                                              
      parainit(irx)=x1a(ib)
  187 para(irx)=x1a(ib)                                                    
      ib=nvc                                                               
      do 188 irx=56+3*nlc,55+4*nlc                                         
      ib=ib+1                                                              
      parainit(irx)=x2a(ib)
  188 para(irx)=x2a(ib)                                                    
      ib=nvc                                                               
      do 189 irx=56+4*nlc,55+5*nlc                                         
      ib=ib+1                                                              
      parainit(irx)=el3a(ib)
  189 para(irx)=el3a(ib)                                                   
      PERT=perr0                                                           
      aeb=a
      EC=E                                                                 
      hjd=hjd0                                                             
      PSH=PSHIFT                                                           
      IRTE=0                                                               
      IRVOL1=0                                                             
      IRVOL2=0                                                             
c   
c  The following lines take care of abundances that may not be among   
c  the 19 Kurucz values (see abun array). abunin is reset at the 
c  allowed value nearest the input value.   
c
      call binnum(abun,19,abunin,iab)   
      dif1=abunin-abun(iab)   
      if(iab.eq.19) goto 7702   
      dif2=abun(iab+1)-abun(iab)   
      dif=dif1/dif2   
      if((dif.ge.0.d0).and.(dif.le.0.5d0)) goto 7702   
      iab=iab+1   
 7702 continue   
      call cofprep(iab)
      call atmcof(iab)
      if(dif1.ne.0.d0) write(6,287) abunin,abun(iab)   
      abunin=abun(iab)   
      if(nbp.eq.nbpsav) goto 7703
      write(6,6066) nbpsav,nbp
      stop
 7703 continue
      CALL MODLOG(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVH,FR1,FR2,HLD,    
     $RM,PHSV,PCSV,GR1,GR2,ALB1,ALB2,N1,N2,F1,F2,MOD,XINCL,THE,MODE,       
     $SNTHH,CSTHH,SNFIH,CSFIH,GRV1,GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,GLUMP1,    
     $GLUMP2,CSBT1,CSBT2,GMAG1,GMAG2,glog1,glog2)                                      
      call ellone(f1,dp,rm,xld,omcr(1),xldd,omd)                           
      rmr=1.d0/rm                                                          
      call ellone(f2,dp,rmr,xld,omcr(2),xldd,omd)                          
      omcr(2)=rm*omcr(2)+(1.d0-rm)*.5d0                                    
      po(1)=phsv                                                           
      po(2)=pcsv                                                           
      CALL VOLUME(VOL1,RM,PHSV,DP,F1,nn1,N1,1,RV,GRX,GRY,GRZ,RVQ,           
     $GRXQ,GRYQ,GRZQ,MMSAVH,FR1,FR2,HLD,SNTHH,CSTHH,SNFIH,CSFIH,SUMMD      
     $,SMD,GRV1,GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2     
     $,GMAG1,GMAG2,glog1,glog2,GR1,1)                                                  
      CALL VOLUME(VOL2,RM,PCSV,DP,F2,N2,N1,2,RV,GRX,GRY,GRZ,RVQ,           
     $GRXQ,GRYQ,GRZQ,MMSAVH,FR1,FR2,HLD,SNTHH,CSTHH,SNFIH,CSFIH,SUMMD      
     $,SMD,GRV1,GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2     
     $,GMAG1,GMAG2,glog1,glog2,GR2,1)                                                  
      IF(E.EQ.0.d0) GOTO 134                                               
      DAP=1.d0+E                                                           
      P1AP=PHSV-2.d0*E*RM/(1.d0-E*E)                                       
      VL1=VOL1                                                             
      CALL VOLUME(VL1,RM,P1AP,DAP,F1,nn1,N1,1,RV,GRX,GRY,GRZ,RVQ,           
     $GRXQ,GRYQ,GRZQ,MMSAVH,FR1,FR2,HLD,SNTHH,CSTHH,SNFIH,CSFIH,SUMMD      
     $,SMD,GRV1,GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2     
     $,GMAG1,GMAG2,glog1,glog2,GR1,2)                                                  
      DPDX1=(PHSV-P1AP)*(1.d0-E*E)*.5d0/E                                  
      P2AP=PCSV-2.d0*E/(1.d0-E*E)                                          
      VL2=VOL2                                                             
      CALL VOLUME(VL2,RM,P2AP,DAP,F2,N2,N1,2,RV,GRX,GRY,GRZ,RVQ,           
     $GRXQ,GRYQ,GRZQ,MMSAVH,FR1,FR2,HLD,SNTHH,CSTHH,SNFIH,CSFIH,SUMMD      
     $,SMD,GRV1,GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2     
     $,GMAG1,GMAG2,glog1,glog2,GR2,2)                                                  
      DPDX2=(PCSV-P2AP)*(1.d0-E*E)*.5d0/E                                  
  134 CONTINUE                                                             
      PHP=PHPER                                                            
      POTH=PHSV                                                            
      POTC=PCSV                                                            
      POT1=PHSV                                                            
      POT2=PCSV                                                            
      if(nlvc.eq.0) goto 728
      DO 24 I=1,NLVC                                                       
      opsf=opsfa(i)                                                        
   24 CALL BBL(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVH,FR1,FR2,HLD,       
     $SLUMP1,SLUMP2,THETA,RHO,AA,BB,POTH,POTC,N1,N2,F1,F2,d,hla(i),        
     $cla(i),x1a(i),x2a(i),y1a(i),y2a(i),gr1,gr2,wla(i),sm1,sm2,tpolh,     
     $tpolc,sbrh,sbrc,tavh,tavc,alb1,alb2,xbol1,xbol2,   
     $ybol1,ybol2,fspot1,fspot2,php,rm,xincl,hot,cool,snthh,csthh,snfih,
     $csfih,tldh,glump1,glump2,xx1,xx2,yy1,yy2,zz1,zz2,dint1,dint2,grv1,
     $grv2,rftemp,rf1,rf2,csbt1,csbt2,gmag1,gmag2,glog1,glog2,obser,
     $fbin1,fbin2,delv1,delv2,      
     $count1,count2,delwl1,delwl2,resf1,resf2,wl1,wl2,dvks1,dvks2,         
     $tau1,tau2,emm1,emm2,hbarw1,hbarw2,xcl,ycl,zcl,rcl,op1,fcl,dens,          
     $encl,edens,taug,emmg,yskp,zskp,mode,iband(i),ifat1,ifat2,1)
  728 continue                                         
      DEL(15)=0.d0                                                         
      del(26)=0.d0                                                         
      del(27)=0.d0                                                         
      del(28)=0.d0                                                         
      del(29)=0.d0                                                         
      del(30)=0.d0                                                         
      del(31)=0.d0                                                         
      del(32)=0.d0                                                         
      del(33)=0.d0                                                         
      del(34)=0.d0                                                         
      del(35)=0.d0                                                         
      DEL(60)=0.d0                                                         
      WRITE(6,*)                                                         
      WRITE(6,399)                                                         
      WRITE(6,55)(DEL(I),I=1,8)                                            
      WRITE(6,*)                                                         
      WRITE(6,402)                                                         
      WRITE(6,55)(DEL(I),I=9,14),(DEL(I),I=16,20)                         
      WRITE(6,*)                                                         
      WRITE(6,403)                                                         
      WRITE(6,55)(DEL(I),I=21,25),(del(i),i=56,59)                         
      WRITE(6,*)                                                         
      WRITE(6,406)                                                         
      WRITE(6,20)(KEP(I),I=1,60),niters,xlamda,vlr
      WRITE(6,*)                                                         
      WRITE(6,166)                                                         
      WRITE(6,66)                                                          
      WRITE(6,61) KSPA,NSPA,KSPB,NSPB                                      
      WRITE(6,*)                                                         
      WRITE(6,707)                                                         
      WRITE(6,708)IFVC1,IFVC2,NLC,iftime,KO,KDISK,ISYM,nppl,ifder,iflcin                       
     $,ifoc
      WRITE(6,*)                                                         
      write(6,917)                                                         
      write(6,912)nref,mref,ifsmv1,ifsmv2,icor1,icor2,if3b,ld1,ld2,
     $kspev,kspot,nomax,ifcgs,maglite,linkext,desextinc                   
      WRITE(6,*)                                                         
      write(6,171)                                                         
      write(6,170) jdphs,hjd0,pzero,dpdt,pshift,delph,nga
      WRITE(6,*)                                                         
      WRITE(6,12)                                                          
      WRITE(6,1)MODE,IPB,IFAT1,IFAT2,N1,N2,N1L,N2L,perr0,dperdt,THE,       
     $vunit,vfac                                                           
      WRITE(6,*)                                                         
      WRITE(6,205)                                                         
      WRITE(6,206) E,A,F1,F2,VGA,XINCL,GR1,GR2,abunin,fspot1,fspot2,
     $nsp1,nsp2,iab                   
      WRITE(6,*)                                                         
      if(ifcgs.ne.0) WRITE(6,54)                                                          
      if(ifcgs.eq.0) write(6,23)                                                          
      if(ifcgs.eq.0) WRITE(6,408) TAVH,TAVC,ALB1,ALB2,PHSV,PCSV,RM,
     $xbol1,xbol2,ybol1,ybol2     
      if(ifcgs.ne.0) WRITE(6,408) TAVH,TAVC,ALB1,ALB2,PHSV,PCSV,RM,
     $xbol1,xbol2,ybol1,ybol2,dpclog     
      write(6,*) 
      write(6,103) 
      write(6,102) a3b,p3b,xinc3b,e3b,perr3b,t03b,rm3b 
      IF(NVC.EQ.0) GOTO 196                                                
      WRITE(6,*)                                                         
      WRITE(6,111)                                                         
      DO 91 I=1,NVC                                                        
   91 WRITE(6,318) iband(I),HLA(I),CLA(I),X1A(I),X2A(I),y1a(i),y2a(i),  
     $opsfa(i),sigma(i),sphas1(i),sphas2(i),sphas3(i),sphas4(i),wla(i),
     $ksd(i)
      if(nlc.eq.0.and.iftime.eq.1) write(6,76) sigma(nlvct),ksd(nlvct)
  196 CONTINUE                                                             
      if(nlvc.eq.0.and.iftime.eq.1) goto 193
      IF(NLC.EQ.0) GOTO 197
      WRITE(6,*)                                                         
      WRITE(6,11)                                                          
      DO 92 I=NVCP,NLVC                                                    
   92 write(6,85)iband(i),hla(i),cla(i),x1a(i),x2a(i),y1a(i),y2a(i),  
     $el3a(i),opsfa(i),noise(i),sigma(i),sphas1(i),sphas2(i),
     $sphas3(i),sphas4(i),ksd(i)
  193 WRITE(6,*)
      if(iftime.eq.1) write(6,76) sigma(nlvct),ksd(nlvct)
      if(nlc.eq.0) goto 197
      write(6,*)
      write(6,88)
      do 95 i=nvcp,nlvc
   95 write(6,86) iband(i),wla(i),aextinc(i),xunit(i),calib(i)
  197 CONTINUE                                                             
      WRITE(6,*)                                                         
      IF(NSTOT.GT.0) WRITE(6,983)                                          
      DO 688 KP=1,2                                                        
      IF((NSP1+KP-1).EQ.0) GOTO 688                                        
      IF((NSP2+(KP-2)**2).EQ.0) GOTO 688                                   
      NSPOT=NSP1                                                           
      IF(KP.EQ.2) NSPOT=NSP2                                               
      DO 687 I=1,NSPOT                                                     
  687 WRITE(6,684)KP,XLAT(KP,I),XLONG(KP,I),RADSP(KP,I),TEMSP(KP,I),
     $tstart(kp,i),tmax1(kp,i),tmax2(kp,i),tfinal(kp,i)        
  688 WRITE(6,*)                                                         
      if(ncl.eq.0) goto 1067                                               
      write(6,69)                                                          
      do 68 i=1,ncl                                                        
   68 write(6,64) xcl(i),ycl(i),zcl(i),rcl(i),op1(i),fcl(i),edens(i),      
     $xmue(i),encl(i),dens(i)                                              
      write(6,*)                                                         
 1067 continue                                                             
      WRITE(6,*)                                                         
      WRITE(6,9)                                                           
c
c   Read all RV curve, light curve, and eclipse timing data in the
c   DO 75 loop
c
      DO 75 LCV=1,nlvct                                                     
      WRITE(6,*)                                                         
      if(lcv.gt.nlvc) goto 759
      IF(LCV.LE.NVC.and.jdphs.eq.2) WRITE(6,955)                           
      IF(LCV.LE.NVC.and.jdphs.eq.1) WRITE(6,755)                           
      IF(LCV.GT.NVC.and.jdphs.eq.2.and.maglite.eq.0) WRITE(6,10)
      IF(LCV.GT.NVC.and.jdphs.eq.2.and.maglite.eq.1) WRITE(6,110)
      IF(LCV.GT.NVC.and.jdphs.eq.1.and.maglite.eq.0) WRITE(6,756)
      IF(LCV.GT.NVC.and.jdphs.eq.1.and.maglite.eq.1) WRITE(6,757)
  759 if(lcv.gt.nlvc) WRITE(6,758)
      DO 74 I=NS,900000                                                      
      ifirst=nppl*(i-1)+NY+1                                               
      last=ifirst+nppl-1                                                   
      if(LCV.gt.nlvc) goto 816
      if(maglite.eq.1.and.lcv.gt.nvc) goto 814
c
c   Read and write a line of RV curve or light curve data
c
      READ(5,131) (phjd(in),flux(in),wt(in),in=ifirst,last)                  
      WRITE(6,131) (phjd(in),flux(in),wt(in),in=ifirst,last)                 
c
c  After writing the input fluxes, multiply by the unit in each band (415 loop).
c  So after 415 the fluxes will have been changed to cgs for internal computations, but have already been
c    printed with their original values that are in terms of the arbitrary unit 'xunit'.
c
      xmult=1.d0
      if(ifcgs.eq.1) xmult=calib(lcv)*xunit(lcv)
      do 415 in=ifirst,last
      if(lcv.gt.nvc.and.ifcgs.eq.1) flux(in)=xmult*flux(in)
  415 continue
      goto 815
c
c  Read and write a line of up to 5 light curve points in magnitudes
c    and then convert to fluxes in DO 127 loop
c
  814 READ(5,131) (phjd(in),xmag(in),wt(in),in=ifirst,last)                  
      WRITE(6,131) (phjd(in),xmag(in),wt(in),in=ifirst,last)                 
      do 127 inn=ifirst,last
      xmult=1.d0
      if(ifcgs.eq.1) xmult=calib(lcv)
      flux(inn)=xmult*10.d0**(-.4d0*xmag(inn))
  127 continue
      goto 815
c
c  Read a line of up to 5 eclipse timings, timing types, and weights
c    (a timing type is 1 for superior conjunction of star 1, 2 for 
c     inferior conjunction of star 1, and 0 if DC is to decide on the
c     minimum type.
c
  816 read(5,129) (phjd(in),mntype(in),wt(in),in=ifirst,last)
c
c  Determine whether the input timing is for min I or min II in the do 812
c    loop (mntype=1 for star 1 superior conjunction, mntype=2 for its
c     inferior conjunction) 
c
      ino=0
      do 812 in=ifirst,last
      mnsv=mntype(in)
      argper=perr0+dperdt*(phjd(in)-hjd0)
      call mintype(ec,argper,pshift,hjd0,pzero,dpdt,phjd(in),newtype)
      if(mnsv.ne.newtype) ino=ino+1
      if(mnsv.eq.0) mntype(in)=newtype
  812 continue
      if(ino.eq.0) write(6,130) (phjd(in),mntype(in),wt(in),in=ifirst,
     $last)
c
c  Write an asterisk at the end of the I/O line if any minimum types on
c    that line have been changed by DC
c
      if(ino.ne.0) write(6,132) (phjd(in),mntype(in),wt(in),in=ifirst,
     $last)
  815 IF(phjd(ifirst).gt.-10000.d0) GOTO 74                                
      NI=-(phjd(ifirst)+10000.d0)                                          
      NY=NY+NI                                                             
      NOBS=nppl*(I-NS-1)+NI                                                
      goto 150
   74 continue
  150 NS=I-1                                                               
      LC1=LCV+1                                                            
   75 KNOBS(LC1)=NOBS+KNOBS(LCV)                                           
c
c  DC loops back to label 209 if a given subset is solved iteratively
c
  209 continue
      xlat(kspa,nspa)=para(1)
      xlong(kspa,nspa)=para(2)
      radsp(kspa,nspa)=para(3)
      temsp(kspa,nspa)=para(4)
      xlat(kspb,nspb)=para(5)
      xlong(kspb,nspb)=para(6)
      radsp(kspb,nspb)=para(7)
      temsp(kspb,nspb)=para(8)
      a=para(9)
      e=para(10)
      perr0=para(11)
      f1=para(12)
      f2=para(13)
      pshift=para(14)
      vga=para(15)
      xincl=para(16)
      gr1=para(17)
      gr2=para(18)
      tavh=para(19)
      tavc=para(20)
      alb1=para(21)
      alb2=para(22)
      if(mode.eq.1) alb2=alb1
      phsv=para(23)
      if(mod.eq.1) para(24)=phsv
      pcsv=para(24)
      rm=para(25)
      dp=1.d0-e
      if(mode.eq.4) call ellone(f1,dp,rm,xld,phsv,xldd,omd)
      if(mode.eq.5) call ellone(f2,dp,rm,xld,pcsv,xldd,omd)
      hjd0=para(26)
      pzero=para(27)
      dpdt=para(28)
      dperdt=para(29)
      if(ifadj3b.eq.0) goto 978
      a3b=para(30)
      arat=a/a3b
      p3b=para(31)
      prat=p3b/pzero
      twopioppr=twopi/p3b
      xinc3b=para(32)
      xinc3brad=dtr*xinc3b 
      sini3b=dsin(xinc3brad) 
      cosi3b=dcos(xinc3brad) 
      e3b=para(33)
      ertfac=dsqrt(1.d0-e3b**2)
      perr3b=para(34)
      sinperr3b=dsin(perr3b)
      cosperr3b=dcos(perr3b)
      t03b=para(35)
  978 continue
      dpclog=para(41)
      desextinc=para(42)
      tstart(kspa,nspa)=para(43)
      tmax1(kspa,nspa)=para(44)
      tmax2(kspa,nspa)=para(45)
      tfinal(kspa,nspa)=para(46)
      tstart(kspb,nspb)=para(47)
      tmax1(kspb,nspb)=para(48)
      tmax2(kspb,nspb)=para(49)
      tfinal(kspb,nspb)=para(50)
      ib=nvc
      do 991 irx=56,55+nlc
      ib=ib+1
  991 hla(ib)=para(irx)
      ib=nvc
      do 986 irx=56+nlc,55+2*nlc
      ib=ib+1
      cla(ib)=para(irx)
  986 continue
      ib=nvc
      do 587 irx=56+2*nlc,55+3*nlc
      ib=ib+1
  587 x1a(ib)=para(irx)
      ib=nvc
      do 888 irx=56+3*nlc,55+4*nlc
      ib=ib+1
      x2a(ib)=para(irx)
      if(mode.eq.1) x2a(ib)=x1a(ib)
  888 continue
      ib=nvc
      do 789 irx=56+4*nlc,55+5*nlc
      ib=ib+1
  789 el3a(ib)=para(irx)
      if(linkext.ne.0)call bandext(desextinc,linkext,xx,ffband,extarray)
      do 114 i=nvcp,nlvc
      if(linkext.eq.0) goto 114
      iban=iband(i)
      aextinc(i)=extarray(iban)
  114 continue
      if(nlc.eq.0) goto 717
      inx=55+nlc
      DO 634 k=nvcp,nlvc                                                       
      opsf=opsfa(k)                                                        
      irterem=irte
      irvol1rem=irvol1
      irvol2rem=irvol2
      irte=0
      irvol1=0
      irvol2=0
      CALL MODLOG(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVH,FR1,FR2,HLD,    
     $RM,PHSV,PCSV,GR1,GR2,ALB1,ALB2,N1,N2,F1,F2,MOD,XINCL,THE,MODE,       
     $SNTHH,CSTHH,SNFIH,CSFIH,GRV1,GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,GLUMP1,    
     $GLUMP2,CSBT1,CSBT2,GMAG1,GMAG2,glog1,glog2)                                      
      CALL BBL(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVH,FR1,FR2,HLD,
     $SLUMP1,SLUMP2,THETA,RHO,AA,BB,PHSV,PCSV,N1,N2,F1,F2,d,hla(k),
     $cla(k),x1a(k),x2a(k),y1a(k),y2a(k),gr1,gr2,wla(k),sm1,sm2,tpolh,
     $tpolc,sbrh,sbrc,tavh,tavc,alb1,alb2,xbol1,xbol2,
     $ybol1,ybol2,fspot1,fspot2,php,rm,xincl,hot,cool,snthh,csthh,snfih,
     $csfih,tldh,glump1,glump2,xx1,xx2,yy1,yy2,zz1,zz2,dint1,dint2,grv1,
     $grv2,rftemp,rf1,rf2,csbt1,csbt2,gmag1,gmag2,glog1,glog2,obser,
     $fbin1,fbin2,delv1,delv2,
     $count1,count2,delwl1,delwl2,resf1,resf2,wl1,wl2,dvks1,dvks2,
     $tau1,tau2,emm1,emm2,hbarw1,hbarw2,xcl,ycl,zcl,rcl,op1,fcl,dens,
     $encl,edens,taug,emmg,yskp,zskp,mode,iband(k),ifat1,ifat2,1)
      ratlum(k)=hla(k)/(hla(k)+cla(k))
      irte=irterem
      irvol1=irvol1rem
      irvol2=irvol2rem
      inx=inx+1
      para(inx)=cla(k)
  634 continue
  717 continue                                         
      do 275 ijp=1,knobs(lc1)                                              
      phas(ijp)=phjd(ijp)                                                  
c 
c  The commented lines here are for the future episodic period change version 
c 
c     if(jdphs.eq.1.and.if3b.gt.0) call litetime(phjd(ijp),0.d0,hjd0, 
c    $pzero,dpdt,xjddum,phas(ijp))     
c     if(jdphs.eq.1.and.if3b.gt.0) call jdph(hjdnew,0.d0,hjd0,pzero, 
c    $dpdt,xjddum,phas(ijp)) 
c
      if(jdphs.eq.1) call jdph(phjd(ijp),0.d0,hjd0,pzero,dpdt,xjddum, 
     $phas(ijp)) 
  275 continue                                                             
      matrix=56                                                            
      do 427 ima=1,55                                                      
  427 matrix=matrix-keep(ima)                                              
      matrix=matrix+nlc*(5-keep(56)-keep(57)-keep(58)-keep(59)-keep(60))   
      MAT=MATRIX-1                                                         
      NCOEFF=MATRIX*KNOBS(LC1)                                             
      NMAT=MAT*KNOBS(LC1)                                                  
      NCOF=NCOEFF                                                          
      if(ibef.eq.1) goto 208
      DO 63 J=1,NCOF                                                       
   63 HOLD(J)=0.d0                                                         
  208 continue
      IF(KOSQ.EQ.1) GOTO 71                                                
      DO 416 IOBS=1,NCOEFF                                                 
  416 OBS(IOBS)=0.d0                                                       
      KSMAX=62                                                             
      KSSMAX=62                                                            
      IF(E.EQ.0.d0) KSMAX=1                                                
      IF(E.NE.0.d0) KSSMAX=1                                               
c
c  Start of outer parameter loop
c
      DO 419 KSS=1,KSSMAX                                                  
c
c  Start of dataset (velocity, light curve, or eclipse timing) loop
c
      DO 417 IB=1,nlvct                                                     
      VC1=0.d0                                                             
      VC2=0.d0                                                             
      ELIT=0.d0                                                            
      IF(IB.GT.NVC.and.ib.le.nlvc) ELIT=1.d0                                              
      IF(IB.EQ.IFVC1) VC1=1.d0                                             
      IF(IB.EQ.(IFVC2*(1+IFVC1))) VC2=1.d0                                 
      vc12=0.d0 
      if(vc1.gt.0.d0.or.vc2.gt.0.d0) vc12=1.d0 
      IST=KNOBS(IB)+1                                                      
      IB1=IB+1                                                             
      ISP=KNOBS(IB1)                                                       
      atten=1.d0
      if(ifcgs.eq.1.and.ib.le.nlvc) atten=10.d0**(-.4d0*aextinc(ib))
c
c Start of observation point loop
c
      DO 418 IX=IST,ISP                                                    
      hjd=phjd(ix)                                                         
      dtime=hjd-hjd0                                                       
      argper=perr0+dperdt*dtime
      if(ib.gt.nlvc) call mintype(ec,argper,pshift,hjd0,pzero,dpdt,
     $phjd(ix),ntype)
      if(ib.gt.nlvc.and.mntype(ix).ne.ntype) goto 418
      period=pzero+dpdt*dtime
      if(if3b.eq.0) goto 266
      xm3b=twopioppr*(hjd-t03b) 
      xmean3b=xm03b+xm3b 
      dmdp3b=-xm3b/p3b
      call kepler(xmean3b,e3b,ecan3b,tr3b) 
c
c  Next line is the default in case E/2 is too close to pi/2 or 3*pi/2 
c
      partr3bdecan3b=dsqrt(1.d0/efac3b) 
      partr3bde3b=0.d0
      ecantest1=dabs(ecan3b-pi) 
      ecantest2=dabs(ecan3b-threepi) 
      ecantest3=dabs(ecan3b+pi) 
      ecantest4=dabs(ecan3b+threepi) 
      if(ecantest1.lt.1.d-9.or.ecantest2.lt.1.d-9) goto 764 
      if(ecantest3.lt.1.d-9.or.ecantest4.lt.1.d-9) goto 764 
      tanhecan=dtan(.5d0*ecan3b) 
      tanhecansq=tanhecan*tanhecan 
      partr3bdecan3b=rtefac3b*(1.d0+tanhecansq)/(1.d0+efac3b*tanhecansq) 
      xnumpartrde=tanhecan*(1.d0+efac3b)/(1.d0-e3b) 
      denompartrde=rtefac3b*(1.d0+efac3b*tanhecan**2) 
      partr3bde3b=xnumpartrde/denompartrde 
  764 continue 
      decande3b=(dsin(ecan3b)+dxmde3b0)/(1.d0-e3b*dcos(ecan3b))
      dtrde3b=partr3bde3b+partr3bdecan3b*decande3b
      dtr3bdp3b=partr3bdecan3b*dmdp3b/(1.d0-e3b*dcos(ecan3b))
      decan3bdperr3b=-parecandtr3b0*(1.d0-e3b*dcos(ecan3b0))/
     $(1.d0-e3b*dcos(ecan3b))
      dtr3bdperr3b=partr3bdecan3b*decan3bdperr3b
      sintr3b=dsin(tr3b) 
      costr3b=dcos(tr3b) 
      ecosp1=1.d0+e3b*costr3b 
      eupsfac=(1.d0-e3b**2)/ecosp1 
      sinsum=dsin(tr3b+perr3b) 
      cossum=dcos(tr3b+perr3b) 
      eccfac=ecos3b+cossum 
      prat=p3b/period
      apfac=arat**3*prat**2
      apfacm=1.d0-apfac
      vfac3b=twopioppr*rvfac*a3b*sini3b*apfacm/ertfac
      vfac3bsav=vfac3b
      deltat=dtcon*a3b*sini3b*apfacm*eupsfac*sinsum 
      pardeltdtr3b=dtcon*a3b*sini3b*apfacm*(1.d0-e3b**2)*(ecosp1*cossum 
     $+sinsum*e3b*sintr3b)/ecosp1**2 
      parvdtr3b=-twopioppr*rvfac*a3b*sini3b*apfacm*sinsum/ertfac 
  266 continue 
c
c Next 12 lines check whether the df/drs numerical derivatives have
c   already been computed or not for spots A and B.
c
      idida=0
      k4344=keep(43)*keep(44)
      k434445=k4344*keep(45)
      if(kh.eq.44.and.keep(43).eq.0) idida=1
      if(kh.eq.45.and.k4344.eq.0) idida=1
      if(kh.eq.46.and.k434445.eq.0) idida=1
      ididb=0
      k4748=keep(47)*keep(48)
      k474849=k4748*keep(49)
      if(kh.eq.48.and.keep(47).eq.0) ididb=1
      if(kh.eq.49.and.k4748.eq.0) ididb=1
      if(kh.eq.50.and.k474849.eq.0) ididb=1
c
c Start of inner parameter loop
c
      DO 420 KS=1,KSMAX                                                    
      norscalc=0
      IRTE=0                                                               
      IRVOL1=0                                                             
      IRVOL2=0                                                             
      IF(E.NE.0.d0) GOTO 297                                               
      IF(IX.NE.IST) IRTE=1                                                 
      IF(IX.EQ.ist) GOTO 297                                                 
      IRVOL1=1                                                             
      IRVOL2=1                                                             
  297 CONTINUE                                                             
      KSR=KS                                                               
      IF(E.EQ.0.d0) KSR=KSS                                                
      if(e.ne.0.d0) goto 1110                                              
      IF(ISYM.NE.1) GOTO 1110                                              
      IF(KSR.EQ.1) GOTO 420                                                
 1110 CONTINUE                                                             
      KH=KSR-2                                                             
      IF(KH.GT.0) GOTO 740                                                 
      KH=1                                                                 
      GOTO 941                                                             
  740 CONTINUE                                                             
      if(ifap*kh.eq.11) goto 842                                           
      if(ifphi*kh.eq.14) goto 842                                          
      if(ifadj3b*kh.eq.14) goto 842 
      IF(KEEP(KH).EQ.1) GOTO 420                                           
  842 IF(E.EQ.0.d0) GOTO 889                                               
      IF(KSR.LE.2) GOTO 889                                                
      IF(KH.LE.8) IRTE=1                                                   
      IF(KH.LE.8) IRVOL1=1                                                 
      IF(KH.LE.8) IRVOL2=1                                                 
      if(kh.eq.9.and.if3b.eq.0) IRTE=1
      if(kh.eq.9.and.if3b.eq.0) IRVOL1=1
      if(kh.eq.9.and.if3b.eq.0) IRVOL2=1
      IF(KH.EQ.12) IRVOL2=1                                                
      IF(KH.EQ.13) IRVOL1=1                                                
      IF(KH.EQ.15) IRTE=1                                                  
      IF(KH.EQ.15) IRVOL1=1                                                
      IF(KH.EQ.15) IRVOL2=1                                                
      IF(KH.EQ.16) IRTE=1                                                  
      IF(KH.EQ.16) IRVOL1=1                                                
      IF(KH.EQ.16) IRVOL2=1                                                
      IF(KH.EQ.17) IRVOL2=1                                                
      IF(KH.EQ.18) IRVOL1=1                                                
      IF(KH.EQ.19) IRVOL1=1                                                
      IF(KH.EQ.19) IRVOL2=1                                                
      IF(KH.EQ.20) IRVOL1=1                                                
      IF(KH.EQ.20) IRVOL2=1                                                
      IF(KH.EQ.21) IRVOL1=1                                                
      IF(KH.EQ.21) IRVOL2=1                                                
      IF(KH.EQ.22) IRVOL1=1                                                
      IF(KH.EQ.22) IRVOL2=1                                                
      IF(KH.EQ.23) IRVOL2=1                                                
      IF(KH.EQ.24) IRVOL1=1                                                
      IF(KH.GE.56) IRVOL1=1                                                
      IF(KH.GE.56) IRVOL2=1                                                
  889 CONTINUE                                                             
      LCF=0                                                                
      IF(KH.GT.55) LCF=IB-NVC                                              
      KPCT1=0                                                              
      KPCT2=0                                                              
      KSP=KH                                                               
      IF(KH.GT.55) KSP=55                                                  
      IF(KH.LT.2) GOTO 808                                                 
      DO 804 ICT=1,KSP                                                     
  804 KPCT1=KPCT1+1-KEEP(ICT)                                              
      GOTO 809                                                             
  808 KPCT1=1                                                              
  809 CONTINUE                                                             
      IF(KH.LT.56) GOTO 806                                                
      DO 805 ICT=56,KH                                                     
  805 KPCT2=KPCT2+1-KEEP(ICT)                                              
      GOTO 807                                                             
  806 KPCT2=1                                                              
  807 CONTINUE                                                             
      II=(KPCT1+NLC*(KPCT2-1)+LCF-1)*KNOBS(LC1)+IX                         
      IF(KH.EQ.9.and.if3b.eq.0) GOTO 300                                                 
      if(KH.eq.9.and.ib.gt.nlvc) goto 765
      if(kh.eq.10.and.ib.gt.nlvc) goto 761
      if(kh.eq.11.and.ib.gt.nlvc) goto 762
      IF(KH.EQ.15) GOTO 308                                                
      if(kh.eq.26) goto 844                                                
      if(kh.eq.27) goto 845                                                
      if(kh.eq.28) goto 846                                                
      if(kh.eq.29) goto 847                                                
      if(kh.eq.30) goto 861                                                
      if(kh.eq.31) goto 862                                                
      if(kh.eq.32) goto 863                                                
      if(kh.eq.33) goto 864                                                
      if(kh.eq.34) goto 865                                                
      if(kh.eq.35) goto 867                                                
      if(kh.eq.41) goto 357                                                
      if(kh.eq.42) goto 355                                                
      if(kh.gt.42.and.kh.lt.47) goto 743
      if(kh.gt.46.and.kh.lt.51) goto 747
      IF(KH.EQ.60) GOTO 301                                                
      IF(KH.NE.56) GOTO 941                                                
      IF(MODE.LE.0) GOTO 941                                               
      IF(IPB.EQ.1) GOTO 941                                                
      if(nlvc.eq.0) goto 420
c
c  Compute analytic derivative for L1
c
      el3g=el3a(ib)
      if(ifcgs.ne.0) el3g=el3a(ib)*atten
      IF(IB.GT.nvc.and.ib.le.nlvc) OBS(II)=(BR(IX)-el3g)/HLA(IB)                      
      GOTO 420                                                             
c
  941 CONTINUE                                                             
      DL=DEL(KH)                                                           
      IF(ISYM.EQ.1) DL=.5d0*DEL(KH)                                        
      SIGN=1.d0                                                            
      ISS=1                                                                
      DO 421 IH=1,ichno                                                       
  421 DLIF(IH)=0.d0                                                        
      IF(KSR.LE.2) GOTO 749                                                
      ISS=IS                                                               
      DLIF(KH)=1.d0                                                        
  749 CONTINUE                                                             
      DO 319 IL=1,ISS                                                      
      IF(E.NE.0.d0) GOTO 4011                                              
      IF(ISYM.EQ.1.and.ix.ne.ist) GOTO 4012                                              
      GOTO 940                                                             
 4011 IF(KSR.LE.2) GOTO 940                                                
      goto 4014                                                            
 4012 IF(IL.EQ.2) GOTO 4016                                                
 4014 IF(LOW(KH).EQ.1) GOTO 314                                            
      VOL1=SVOL1                                                           
      VOL2=SVOL2                                                           
      SUMM1=SSUM1                                                          
      SUMM2=SSUM2                                                          
      SM1=SSM1                                                             
      SM2=SSM2                                                             
      sbrh=ssbrh
      sbrc=ssbrc
      polin1=spolin1
      polin2=spolin2
      DO 851 IRE=1,MTLH                                                    
  851 TLDH(IRE)=STLDH(IRE)                                                 
      DO 508 IRE=1,M1H                                                     
      RV(IRE)=SRV(IRE)                                                     
      GRX(IRE)=SGRX(IRE)                                                   
      GRY(IRE)=SGRY(IRE)                                                   
      GRZ(IRE)=SGRZ(IRE)                                                   
      GLUMP1(IRE)=SGLM1(IRE)                                               
      GRV1(IRE)=SGRV1(IRE)                                                 
      XX1(IRE)=SXX1(IRE)                                                   
      YY1(IRE)=SYY1(IRE)                                                   
      ZZ1(IRE)=SZZ1(IRE)                                                   
      GMAG1(IRE)=SGMG1(IRE)                                                
      CSBT1(IRE)=SCSB1(IRE)                                                
      RF1(IRE)=SRF1(IRE)                                                   
      FR1(IRE)=SFR1(IRE)                                                   
  508 SLUMP1(IRE)=SLMP1(IRE)                                               
      DO 309 IRE=1,M2H                                                     
      RVQ(IRE)=SRVQ(IRE)                                                   
      GRXQ(IRE)=SGRXQ(IRE)                                                 
      GRYQ(IRE)=SGRYQ(IRE)                                                 
      GRZQ(IRE)=SGRZQ(IRE)                                                 
      GLUMP2(IRE)=SGLM2(IRE)                                               
      GRV2(IRE)=SGRV2(IRE)                                                 
      XX2(IRE)=SXX2(IRE)                                                   
      YY2(IRE)=SYY2(IRE)                                                   
      ZZ2(IRE)=SZZ2(IRE)                                                   
      GMAG2(IRE)=SGMG2(IRE)                                                
      CSBT2(IRE)=SCSB2(IRE)                                                
      RF2(IRE)=SRF2(IRE)                                                   
      FR2(IRE)=SFR2(IRE)                                                   
  309 SLUMP2(IRE)=SLMP2(IRE)                                               
      GOTO 940                                                             
 4016 IF(LOW(KH).EQ.1) GOTO 4018                                           
      VOL1=EVOL1                                                           
      VOL2=EVOL2                                                           
      SUMM1=ESUM1                                                          
      SUMM2=ESUM2                                                          
      SM1=ESM1                                                             
      SM2=ESM2                                                             
      sbrh=esbrh
      sbrc=esbrc
      polin1=epolin1
      polin2=epolin2
      DO 852 IRE=1,MTLH                                                    
  852 TLDH(IRE)=ETLDH(IRE)                                                 
      DO 1508 IRE=1,M1H                                                    
      RV(IRE)=ERV(IRE)                                                     
      GRX(IRE)=EGRX(IRE)                                                   
      GRY(IRE)=EGRY(IRE)                                                   
      GRZ(IRE)=EGRZ(IRE)                                                   
      GLUMP1(IRE)=EGLM1(IRE)                                               
      GRV1(IRE)=EGRV1(IRE)                                                 
      XX1(IRE)=EXX1(IRE)                                                   
      YY1(IRE)=EYY1(IRE)                                                   
      ZZ1(IRE)=EZZ1(IRE)                                                   
      GMAG1(IRE)=EGMG1(IRE)                                                
      CSBT1(IRE)=ECSB1(IRE)                                                
      RF1(IRE)=ERF1(IRE)                                                   
      FR1(IRE)=EFR1(IRE)                                                   
 1508 SLUMP1(IRE)=ELMP1(IRE)                                               
      DO 1309 IRE=1,M2H                                                    
      RVQ(IRE)=ERVQ(IRE)                                                   
      GRXQ(IRE)=EGRXQ(IRE)                                                 
      GRYQ(IRE)=EGRYQ(IRE)                                                 
      GRZQ(IRE)=EGRZQ(IRE)                                                 
      GLUMP2(IRE)=EGLM2(IRE)                                               
      GRV2(IRE)=EGRV2(IRE)                                                 
      XX2(IRE)=EXX2(IRE)                                                   
      YY2(IRE)=EYY2(IRE)                                                   
      ZZ2(IRE)=EZZ2(IRE)                                                   
      GMAG2(IRE)=EGMG2(IRE)                                                
      CSBT2(IRE)=ECSB2(IRE)                                                
      RF2(IRE)=ERF2(IRE)                                                   
      FR2(IRE)=EFR2(IRE)                                                   
 1309 SLUMP2(IRE)=ELMP2(IRE)                                               
      GOTO 940                                                             
 4018 CONTINUE                                                             
      VOL1=EVOL1L                                                          
      VOL2=EVOL2L                                                          
      SUMM1=ESUM1L                                                         
      SUMM2=ESUM2L                                                         
      SM1=ESM1L                                                            
      SM2=ESM2L                                                            
      sbrh=esbrhl
      sbrc=esbrcl
      polin1=epolin1L
      polin2=epolin2L
      DO 853 IRE=1,MTLL                                                    
  853 TLDL(IRE)=ETLDL(IRE)                                                 
      DO 310 IRE=1,M1L                                                     
      RV(IRE)=ERVL(IRE)                                                    
      GRX(IRE)=EGRXL(IRE)                                                  
      GRY(IRE)=EGRYL(IRE)                                                  
      GRZ(IRE)=EGRZL(IRE)                                                  
      GLUMP1(IRE)=EGLM1L(IRE)                                              
      GRV1(IRE)=EGRV1L(IRE)                                                
      XX1(IRE)=EXX1L(IRE)                                                  
      YY1(IRE)=EYY1L(IRE)                                                  
      ZZ1(IRE)=EZZ1L(IRE)                                                  
      GMAG1(IRE)=EGMG1L(IRE)                                               
      CSBT1(IRE)=ECSB1L(IRE)                                               
      RF1(IRE)=ERF1L(IRE)                                                  
      FR1(IRE)=EFR1L(IRE)                                                  
  310 SLUMP1(IRE)=ELMP1L(IRE)                                              
      DO 311 IRE=1,M2L                                                     
      RVQ(IRE)=ERVQL(IRE)                                                  
      GRXQ(IRE)=EGRXQL(IRE)                                                
      GRYQ(IRE)=EGRYQL(IRE)                                                
      GRZQ(IRE)=EGRZQL(IRE)                                                
      GLUMP2(IRE)=EGLM2L(IRE)                                              
      GRV2(IRE)=EGRV2L(IRE)                                                
      XX2(IRE)=EXX2L(IRE)                                                  
      YY2(IRE)=EYY2L(IRE)                                                  
      ZZ2(IRE)=EZZ2L(IRE)                                                  
      GMAG2(IRE)=EGMG2L(IRE)                                               
      CSBT2(IRE)=ECSB2L(IRE)                                               
      RF2(IRE)=ERF2L(IRE)                                                  
      FR2(IRE)=EFR2L(IRE)                                                  
  311 SLUMP2(IRE)=ELMP2L(IRE)                                              
      GOTO 940                                                             
  314 CONTINUE                                                             
      VOL1=SVOL1L                                                          
      VOL2=SVOL2L                                                          
      SUMM1=SSUM1L                                                         
      SUMM2=SSUM2L                                                         
      SM1=SSM1L                                                            
      SM2=SSM2L                                                            
      sbrh=ssbrhl
      sbrc=ssbrcl
      polin1=spolin1L
      polin2=spolin2L
      DO 854 IRE=1,MTLL                                                    
  854 TLDL(IRE)=STLDL(IRE)                                                 
      DO 1310 IRE=1,M1L                                                    
      RV(IRE)=SRVL(IRE)                                                    
      GRX(IRE)=SGRXL(IRE)                                                  
      GRY(IRE)=SGRYL(IRE)                                                  
      GRZ(IRE)=SGRZL(IRE)                                                  
      GLUMP1(IRE)=SGLM1L(IRE)                                              
      GRV1(IRE)=SGRV1L(IRE)                                                
      XX1(IRE)=SXX1L(IRE)                                                  
      YY1(IRE)=SYY1L(IRE)                                                  
      ZZ1(IRE)=SZZ1L(IRE)                                                  
      GMAG1(IRE)=SGMG1L(IRE)                                               
      CSBT1(IRE)=SCSB1L(IRE)                                               
      RF1(IRE)=SRF1L(IRE)                                                  
      FR1(IRE)=SFR1L(IRE)                                                  
 1310 SLUMP1(IRE)=SLMP1L(IRE)                                              
      DO 1311 IRE=1,M2L                                                    
      RVQ(IRE)=SRVQL(IRE)                                                  
      GRXQ(IRE)=SGRXQL(IRE)                                                
      GRYQ(IRE)=SGRYQL(IRE)                                                
      GRZQ(IRE)=SGRZQL(IRE)                                                
      GLUMP2(IRE)=SGLM2L(IRE)                                              
      GRV2(IRE)=SGRV2L(IRE)                                                
      XX2(IRE)=SXX2L(IRE)                                                  
      YY2(IRE)=SYY2L(IRE)                                                  
      ZZ2(IRE)=SZZ2L(IRE)                                                  
      GMAG2(IRE)=SGMG2L(IRE)                                               
      CSBT2(IRE)=SCSB2L(IRE)                                               
      RF2(IRE)=SRF2L(IRE)                                                  
      FR2(IRE)=SFR2L(IRE)                                                  
 1311 SLUMP2(IRE)=SLMP2L(IRE)                                              
  940 CONTINUE                                                             
      DELS=DL*SIGN                                                         
      SIGN=-1.d0                                                           
      IF(NSPA.EQ.0) GOTO 470                                               
      xlt=xlat(kspa,nspa)+dels*dlif(1)                                     
      xlng(kspa,nspa)=xlong(kspa,nspa)+dels*dlif(2)                        
      snlat(kspa,nspa)=dsin(xlt)                                           
      cslat(kspa,nspa)=dcos(xlt)                                           
      snlng(kspa,nspa)=dsin(xlng(kspa,nspa))                               
      cslng(kspa,nspa)=dcos(xlng(kspa,nspa))                               
      RDSP(KSPA,NSPA)=RADSP(KSPA,NSPA)+DELS*DLIF(3)                        
      amax(kspa,nspa)=twopi*(1.d0-dcos(rdsp(kspa,nspa)))
      TMSP(KSPA,NSPA)=TEMSP(KSPA,NSPA)+DELS*DLIF(4)                        
      tstrt(kspa,nspa)=tstart(kspa,nspa)+dels*dlif(43)
      tmx1(kspa,nspa)=tmax1(kspa,nspa)+dels*dlif(44)
      tmx2(kspa,nspa)=tmax2(kspa,nspa)+dels*dlif(45)
      tfnal(kspa,nspa)=tfinal(kspa,nspa)+dels*dlif(46)
  470 CONTINUE                                                             
      IF(NSPB.EQ.0) GOTO 471                                               
      xlt=xlat(kspb,nspb)+dels*dlif(5)                                     
      xlng(kspb,nspb)=xlong(kspb,nspb)+dels*dlif(6)                        
      snlat(kspb,nspb)=dsin(xlt)                                           
      cslat(kspb,nspb)=dcos(xlt)                                           
      snlng(kspb,nspb)=dsin(xlng(kspb,nspb))                               
      cslng(kspb,nspb)=dcos(xlng(kspb,nspb))                               
      RDSP(KSPB,NSPB)=RADSP(KSPB,NSPB)+DELS*DLIF(7)                        
      amax(kspb,nspb)=twopi*(1.d0-dcos(rdsp(kspb,nspb)))
      TMSP(KSPB,NSPB)=TEMSP(KSPB,NSPB)+DELS*DLIF(8)                        
      tstrt(kspb,nspb)=tstart(kspb,nspb)+dels*dlif(47)
      tmx1(kspb,nspb)=tmax1(kspb,nspb)+dels*dlif(48)
      tmx2(kspb,nspb)=tmax2(kspb,nspb)+dels*dlif(49)
      tfnal(kspb,nspb)=tfinal(kspb,nspb)+dels*dlif(50)
  471 CONTINUE                                                             
      aeb=a+dels*dlif(9)
      EC=E+DELS*DLIF(10)                                                   
      pconsc=0.d0
      pconic=0.5d0
      argper=perr0+dperdt*dtime
      if(ec.gt.0.d0.and.ib.gt.nlvc) call conjph(ec,argper,pshift,trsc,
     $tric,econsc,econic,xmsc,xmic,pconsc,pconic)
      PERT=perr0+DELS*DLIF(11)                                             
      FF1=F1+DELS*DLIF(12)                                                 
      FF2=F2+DELS*DLIF(13)                                                 
      PSH=PSHIFT+DELS*DLIF(14)                                             
      XINC=XINCL+DELS*DLIF(16)                                             
      G1=GR1+DELS*DLIF(17)                                                 
      G2=GR2+DELS*DLIF(18)                                                 
      T1=TAVH+DELS*DLIF(19)                                                
      T2=TAVC+DELS*DLIF(20)                                                
      A1=ALB1+DELS*DLIF(21)                                                
      A2=ALB2+DELS*DLIF(22)                                                
      POT1=PHSV+DELS*DLIF(23)                                              
      POT2=PCSV+DELS*DLIF(24)                                              
      RMASS=RM+DELS*DLIF(25)                                               
      HL=HLA(IB)+DELS*DLIF(56)                                             
      CL=CLA(IB)+DELS*DLIF(57)                                             
      X1=X1A(IB)+DELS*DLIF(58)                                             
      X2=X2A(IB)+DELS*DLIF(59)                                             
      y1=y1a(ib)                                                           
      y2=y2a(ib)                                                           
      opsf=opsfa(ib)                                                       
      phasin0=phas(ix) 
      if(if3b.eq.0) goto 8722
c
c  Note that parameter 9 ('a', orbit size of the inner 2-body system)
c    is treated differently from the other parameters in regard to 
c    the light time effect (next dozen or so lines). The reason is
c    that 'a' is the only parameter (p) that 
c    (A) affects light times and 
c    (B) has its derivatives [partial(f)/partial(p)] computed 
c    NUMERICALLY in DC.
c
      if(kh.ne.9) goto 8721
      prat=p3b/period
      arateb=aeb/a3b
      apfaceb=arateb**3*prat**2
      apfacmeb=1.d0-apfaceb
      deltateb=dtcon*a3b*sini3b*apfacmeb*eupsfac*sinsum
      vfac3b=twopioppr*rvfac*a3b*sini3b*apfacmeb/ertfac
      phasin0=phas(ix)-deltateb/period
      goto 8722
 8721 continue
      phasin0=phas(ix)-deltat/period
      vfac3b=vfac3bsav
 8722 continue
      phasin=phasin0
      if(ib.gt.nlvc) goto 801
      IF(KSR.EQ.1) GOTO 802                                                
      IF(KSR.EQ.2) GOTO 872                                                
      IF(LOW(KH).EQ.1) GOTO 802                                            
  872 continue
      do 550 iga=1,nga
      if(nga.eq.1) goto 552
      plusmin=-1.d0
      if(iga.le.nhalf) ind=ind0(nga)+iga-1
      if(iga.gt.nhalf) ind=ind0(nga)+nga-iga
      if(iga.gt.nhalf) plusmin=1.d0
      abscis=plusmin*aba(ind)
      phasin=phasin0+emmm*abscis
  552 CALL MODLOG(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVH,FR1,FR2,HLD,
     $RMASS,POT1,POT2,G1,G2,A1,A2,N1,N2,FF1,FF2,MOD,XINC,THE,MODE,         
     $SNTHH,CSTHH,SNFIH,CSFIH,GRV1,GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,GLUMP1,    
     $GLUMP2,CSBT1,CSBT2,GMAG1,GMAG2,glog1,glog2)                                      
      CALL BBL(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVH,FR1,FR2,HLD,
     $SLUMP1,SLUMP2,THETA,RHO,AA,BB,POT1,POT2,N1,N2,FF1,FF2,d,hl,cl,x1,
     $x2,y1,y2,g1,g2,wla(ib),sm1,sm2,tph,tpc,sbrh,sbrc,t1,
     $t2,a1,a2,xbol1,xbol2,ybol1,ybol2,fspot1,fspot2,phasin,rmass,xinc,
     $hotr(iga),coolr(iga),snthh,csthh,snfih,csfih,tldh,glump1,glump2,
     $xx1,xx2,yy1,yy2,zz1,zz2,dint1,dint2,grv1,grv2,rftemp,rf1,rf2,
     $csbt1,csbt2,gmag1,gmag2,glog1,glog2,obser,
     $fbin1,fbin2,delv1,delv2,count1,count2,delwl1,delwl2,resf1,resf2,
     $wl1,wl2,dvks1,dvks2,tau1,tau2,emm1,emm2,hbarw1,hbarw2,xcl,ycl,zcl,
     $rcl,op1,fcl,dens,encl,edens,taug,emmg,yskp,zskp,mode,iband(ib),
     $ifat1,ifat2,1)
  550 continue
      hot=hotr(1)
      cool=coolr(1)
      if(nga.eq.1) goto 801
      call gaussquad(2,nga,emmm,hotr,hotg)
      call gaussquad(2,nga,emmm,coolr,coolg)
      hot=hotg/delph
      cool=coolg/delph
      GOTO 801                                                             
  802 CONTINUE                                                             
      do 551 iga=1,nga
      if(nga.eq.1) goto 553
      plusmin=-1.d0
      if(iga.le.nhalf) ind=ind0(nga)+iga-1
      if(iga.gt.nhalf) ind=ind0(nga)+nga-iga
      if(iga.gt.nhalf) plusmin=1.d0
      abscis=plusmin*aba(ind)
      phasin=phasin0+.5d0*delph*abscis
  553 CALL MODLOG(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVL,FR1,FR2,HLD,
     $RMASS,POT1,POT2,G1,G2,A1,A2,N1L,N2L,FF1,FF2,MOD,XINC,THE,MODE,       
     $SNTHL,CSTHL,SNFIL,CSFIL,GRV1,GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,GLUMP1,    
     $GLUMP2,CSBT1,CSBT2,GMAG1,GMAG2,glog1,glog2)                                      
      CALL BBL(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVL,FR1,FR2,HLD,
     $SLUMP1,SLUMP2,THETA,RHO,AA,BB,POT1,POT2,N1L,N2L,FF1,FF2,d,hl,cl,
     $x1,x2,y1,y2,g1,g2,wla(ib),sm1,sm2,tph,tpc,sbrh,sbrc,
     $t1,t2,a1,a2,xbol1,xbol2,ybol1,ybol2,fspot1,fspot2,phasin,rmass,
     $xinc,hotr(iga),coolr(iga),snthl,csthl,snfil,csfil,tldl,glump1,
     $glump2,xx1,xx2,yy1,yy2,zz1,zz2,dint1,dint2,grv1,grv2,rftemp,rf1,
     $rf2,csbt1,csbt2,gmag1,gmag2,glog1,glog2,obser,
     $fbin1,fbin2,delv1,delv2,count1,count2,delwl1,delwl2,resf1,resf2,
     $wl1,wl2,dvks1,dvks2,tau1,tau2,emm1,emm2,hbarw1,hbarw2,xcl,ycl,zcl,
     $rcl,op1,fcl,dens,encl,edens,taug,emmg,yskp,zskp,mode,iband(ib),
     $ifat1,ifat2,1)
  551 continue
      hot=hotr(1)
      cool=coolr(1)
      if(nga.eq.1) goto 801
      call gaussquad(2,nga,emmm,hotr,hotg)
      call gaussquad(2,nga,emmm,coolr,coolg)
      hot=hotg/delph
      cool=coolg/delph
  801 CONTINUE                                                               
      vfac3b=vfac3bsav
      IF(E.NE.0.d0) GOTO 4111                                              
      IF(ISYM.EQ.0) GOTO 602                                               
      IF(IX.NE.IST) GOTO 602                                               
      if(ksr.eq.1) goto 602
      if(ksr.le.2) goto 4119
      IF(LOW(KH).EQ.1) GOTO 4112                                           
      GOTO 4119                                                            
 4111 IF(IRTE.EQ.1) GOTO 602                                               
      IF(KSR.GE.3) GOTO 602                                                
      IF(KSR.EQ.2) GOTO 4119                                               
 4112 IF(IL.EQ.2) GOTO 4116                                                
      SVOL1L=VOL1                                                          
      SVOL2L=VOL2                                                          
      SSUM1L=SUMM1                                                         
      SSUM2L=SUMM2                                                         
      SSM1L=SM1                                                            
      SSM2L=SM2                                                            
      ssbrhl=sbrh
      ssbrcl=sbrc
      spolin1L=polin1
      spolin2L=polin2
      DO 855 IHL=1,MTLL                                                    
  855 STLDL(IHL)=TLDL(IHL)                                                 
      DO 603 IHL=1,M1L                                                     
      SRVL(IHL)=RV(IHL)                                                    
      SGRXL(IHL)=GRX(IHL)                                                  
      SGRYL(IHL)=GRY(IHL)                                                  
      SGRZL(IHL)=GRZ(IHL)                                                  
      SLMP1L(IHL)=SLUMP1(IHL)                                              
      SGLM1L(IHL)=GLUMP1(IHL)                                              
      SGRV1L(IHL)=GRV1(IHL)                                                
      SXX1L(IHL)=XX1(IHL)                                                  
      SYY1L(IHL)=YY1(IHL)                                                  
      SZZ1L(IHL)=ZZ1(IHL)                                                  
      SGMG1L(IHL)=GMAG1(IHL)                                               
      SCSB1L(IHL)=CSBT1(IHL)                                               
      SRF1L(IHL)=RF1(IHL)                                                  
      SFR1L(IHL)=FR1(IHL)                                                  
  603 CONTINUE                                                             
      DO 606 IHL=1,M2L                                                     
      SRVQL(IHL)=RVQ(IHL)                                                  
      SGRXQL(IHL)=GRXQ(IHL)                                                
      SGRYQL(IHL)=GRYQ(IHL)                                                
      SGRZQL(IHL)=GRZQ(IHL)                                                
      SLMP2L(IHL)=SLUMP2(IHL)                                              
      SGLM2L(IHL)=GLUMP2(IHL)                                              
      SGRV2L(IHL)=GRV2(IHL)                                                
      SXX2L(IHL)=XX2(IHL)                                                  
      SYY2L(IHL)=YY2(IHL)                                                  
      SZZ2L(IHL)=ZZ2(IHL)                                                  
      SGMG2L(IHL)=GMAG2(IHL)                                               
      SCSB2L(IHL)=CSBT2(IHL)                                               
      SRF2L(IHL)=RF2(IHL)                                                  
      SFR2L(IHL)=FR2(IHL)                                                  
  606 CONTINUE                                                             
      GOTO 602                                                             
 4116 CONTINUE                                                             
      EVOL1L=VOL1                                                          
      EVOL2L=VOL2                                                          
      ESUM1L=SUMM1                                                         
      ESUM2L=SUMM2                                                         
      ESM1L=SM1                                                            
      ESM2L=SM2                                                            
      esbrhl=sbrh
      esbrcl=sbrc
      epolin1L=polin1
      epolin2L=polin2
      DO 856 IHL=1,MTLL                                                    
  856 ETLDL(IHL)=TLDL(IHL)                                                 
      DO 1603 IHL=1,M1L                                                    
      ERVL(IHL)=RV(IHL)                                                    
      EGRXL(IHL)=GRX(IHL)                                                  
      EGRYL(IHL)=GRY(IHL)                                                  
      EGRZL(IHL)=GRZ(IHL)                                                  
      ELMP1L(IHL)=SLUMP1(IHL)                                              
      EGLM1L(IHL)=GLUMP1(IHL)                                              
      EGRV1L(IHL)=GRV1(IHL)                                                
      EXX1L(IHL)=XX1(IHL)                                                  
      EYY1L(IHL)=YY1(IHL)                                                  
      EZZ1L(IHL)=ZZ1(IHL)                                                  
      EGMG1L(IHL)=GMAG1(IHL)                                               
      ECSB1L(IHL)=CSBT1(IHL)                                               
      ERF1L(IHL)=RF1(IHL)                                                  
      EFR1L(IHL)=FR1(IHL)                                                  
 1603 CONTINUE                                                             
      DO 1606 IHL=1,M2L                                                    
      ERVQL(IHL)=RVQ(IHL)                                                  
      EGRXQL(IHL)=GRXQ(IHL)                                                
      EGRYQL(IHL)=GRYQ(IHL)                                                
      EGRZQL(IHL)=GRZQ(IHL)                                                
      ELMP2L(IHL)=SLUMP2(IHL)                                              
      EGLM2L(IHL)=GLUMP2(IHL)                                              
      EGRV2L(IHL)=GRV2(IHL)                                                
      EXX2L(IHL)=XX2(IHL)                                                  
      EYY2L(IHL)=YY2(IHL)                                                  
      EZZ2L(IHL)=ZZ2(IHL)                                                  
      EGMG2L(IHL)=GMAG2(IHL)                                               
      ECSB2L(IHL)=CSBT2(IHL)                                               
      ERF2L(IHL)=RF2(IHL)                                                  
      EFR2L(IHL)=FR2(IHL)                                                  
 1606 CONTINUE                                                             
      GOTO 602                                                             
 4119 IF(IL.EQ.2) GOTO 4120                                                
      SVOL1=VOL1                                                           
      SVOL2=VOL2                                                           
      SSUM1=SUMM1                                                          
      SSUM2=SUMM2                                                          
      SSM1=SM1                                                             
      SSM2=SM2                                                             
      ssbrh=sbrh
      ssbrc=sbrc
      spolin1=polin1
      spolin2=polin2
      DO 857 IHH=1,MTLH                                                    
  857 STLDH(IHH)=TLDH(IHH)                                                 
      DO 601 IHH=1,M1H                                                     
      SRV(IHH)=RV(IHH)                                                     
      SGRX(IHH)=GRX(IHH)                                                   
      SGRY(IHH)=GRY(IHH)                                                   
      SGRZ(IHH)=GRZ(IHH)                                                   
      SLMP1(IHH)=SLUMP1(IHH)                                               
      SGLM1(IHH)=GLUMP1(IHH)                                               
      SGRV1(IHH)=GRV1(IHH)                                                 
      SXX1(IHH)=XX1(IHH)                                                   
      SYY1(IHH)=YY1(IHH)                                                   
      SZZ1(IHH)=ZZ1(IHH)                                                   
      SGMG1(IHH)=GMAG1(IHH)                                                
      SCSB1(IHH)=CSBT1(IHH)                                                
      SRF1(IHH)=RF1(IHH)                                                   
      SFR1(IHH)=FR1(IHH)                                                   
  601 CONTINUE                                                             
      DO 605 IHH=1,M2H                                                     
      SRVQ(IHH)=RVQ(IHH)                                                   
      SGRXQ(IHH)=GRXQ(IHH)                                                 
      SGRYQ(IHH)=GRYQ(IHH)                                                 
      SGRZQ(IHH)=GRZQ(IHH)                                                 
      SLMP2(IHH)=SLUMP2(IHH)                                               
      SGLM2(IHH)=GLUMP2(IHH)                                               
      SGRV2(IHH)=GRV2(IHH)                                                 
      SXX2(IHH)=XX2(IHH)                                                   
      SYY2(IHH)=YY2(IHH)                                                   
      SZZ2(IHH)=ZZ2(IHH)                                                   
      SGMG2(IHH)=GMAG2(IHH)                                                
      SCSB2(IHH)=CSBT2(IHH)                                                
      SRF2(IHH)=RF2(IHH)                                                   
      SFR2(IHH)=FR2(IHH)                                                   
  605 CONTINUE                                                             
      GOTO 602                                                             
 4120 CONTINUE                                                             
      EVOL1=VOL1                                                           
      EVOL2=VOL2                                                           
      ESUM1=SUMM1                                                          
      ESUM2=SUMM2                                                          
      ESM1=SM1                                                             
      ESM2=SM2                                                             
      esbrh=sbrh
      esbrc=sbrc
      epolin1=polin1
      epolin2=polin2
      DO 858 IHH=1,MTLH                                                    
  858 ETLDH(IHH)=TLDH(IHH)                                                 
      DO 1601 IHH=1,M1H                                                    
      ERV(IHH)=RV(IHH)                                                     
      EGRX(IHH)=GRX(IHH)                                                   
      EGRY(IHH)=GRY(IHH)                                                   
      EGRZ(IHH)=GRZ(IHH)                                                   
      ELMP1(IHH)=SLUMP1(IHH)                                               
      EGLM1(IHH)=GLUMP1(IHH)                                               
      EGRV1(IHH)=GRV1(IHH)                                                 
      EXX1(IHH)=XX1(IHH)                                                   
      EYY1(IHH)=YY1(IHH)                                                   
      EZZ1(IHH)=ZZ1(IHH)                                                   
      EGMG1(IHH)=GMAG1(IHH)                                                
      ECSB1(IHH)=CSBT1(IHH)                                                
      ERF1(IHH)=RF1(IHH)                                                   
      EFR1(IHH)=FR1(IHH)                                                   
 1601 CONTINUE                                                             
      DO 1605 IHH=1,M2H                                                    
      ERVQ(IHH)=RVQ(IHH)                                                   
      EGRXQ(IHH)=GRXQ(IHH)                                                 
      EGRYQ(IHH)=GRYQ(IHH)                                                 
      EGRZQ(IHH)=GRZQ(IHH)                                                 
      ELMP2(IHH)=SLUMP2(IHH)                                               
      EGLM2(IHH)=GLUMP2(IHH)                                               
      EGRV2(IHH)=GRV2(IHH)                                                 
      EXX2(IHH)=XX2(IHH)                                                   
      EYY2(IHH)=YY2(IHH)                                                   
      EZZ2(IHH)=ZZ2(IHH)                                                   
      EGMG2(IHH)=GMAG2(IHH)                                                
      ECSB2(IHH)=CSBT2(IHH)                                                
      ERF2(IHH)=RF2(IHH)                                                   
      EFR2(IHH)=FR2(IHH)                                                   
 1605 CONTINUE                                                             
  602 CONTINUE                                                             
      HTT=HOT                                                              
      IF(MODE.EQ.-1) HTT=0.d0                                              
      XR=(HTT+COOL+el3a(ib))*ELIT+VKM1*VC1+VKM2*VC2
      phasn=phas(ix)
 9488 continue
      if(phasn.gt.1.d0) phasn=phasn-1.d0
      if(phasn.gt.1.d0) goto 9488
 9493 continue
      if(phasn.lt.0.d0) phasn=phasn+1.d0
      if(phasn.lt.0.d0) goto 9493
c
c  Compute observable flux in cgs units
c
      cgs1=1.d0
      cgs2=1.d0
      if(ifcgs.eq.0) goto 358
      dpc=10.d0**dpclog
      aodsq=(pcrsun*a/dpc)**2
      cgs1=aodsq*htt*polin1/sbrh
      cgs2=aodsq*cool*polin2/sbrc
      cgsta(ix)=cgs1+cgs2+el3a(ib)
      cgstot=cgsta(ix)*atten
      xr=cgstot*elit+vkm1*vc1+vkm2*vc2
  358 continue
      IF(KSR.NE.1) GOTO 710                                                
      BL(IX)=XR                                                            
      GOTO 420                                                             
  710 CONTINUE                                                             
      IF(KSR.NE.2) GOTO 711                                                
      BR(IX)=XR                                                            
      II=NMAT+IX                                                           
      OBS(II)=FLUX(IX)-XR                                                  
      if(ib.le.nlvc) goto 822
c
c  Compute timing residuals (down to label 822)
c
c  Compute phase corresponding to observed eclipse time.
c
      call jdph(phjd(ix),0.d0,hjd0,pzero,dpdt,xjddum,phaz) 
      xmp1=2.d0-dfloat(mntype(ix))
      xmp2=dfloat(mntype(ix))-1.d0
      phazconj=xmp1*pconsc+xmp2*pconic
c
c  Shift observed phase in accordance with eccentricity and eclipse type 
c       (primary vs. secondary eclipse).
c
      phaze=phaz-phazconj-deltat/period
c
c  Round observed conjunction to nearest whole cycle and restore eccentricity
c        and eclipse type shift.
c
      phazcomp=dfloat(nint(phaze))+phazconj
c
c  Compute time corresponding to conjunction phase and ephemeris.
c 
      call jdph(0.d0,phazcomp,hjd0,pzero,dpdt,hjdcomp,phdum) 
c  Difference observed and computed timing to get timing residual.
c
      timecomp(ix)=hjdcomp+deltat
      obs(ii)=phjd(ix)-timecomp(ix)
  822 continue
      if(ib.gt.nlvc) goto 420
      if(iss.eq.2) goto 319                                                
      GOTO 420                                                             
  711 CONTINUE                                                             
      XLR=BR(IX)                                                           
      IF(LOW(KH).EQ.1) XLR=BL(IX)                                          
      IF(IL.NE.2) GOTO 388                                                 
      XLR=XR                                                               
      GOTO 87                                                              
  388 XNR=XR                                                               
   87 khkeep=kh*keep(kh)                                                   
      if(khkeep.ne.3.and.khkeep.ne.7.and.khkeep.ne.11.and.khkeep.ne.14)
     $OBS(II)=(XNR-XLR)/DEL(KH)
      if(kh.eq.11) dfdap(ix)=(xnr-xlr)/del(kh)                             
      if(kh.eq.14) dfdph(ix)=(xnr-xlr)/del(kh)   
  319 CONTINUE                                                             
      GOTO 420                                                             
  300 IF(IB.LE.NVC) OBS(II)=(BR(IX)-VGA)/A                                 
      if(ib.gt.nvc.and.ib.le.nlvc.and.ifcgs.eq.1) obs(ii)=2.d0*br(ix)/a
      GOTO 420                                                             
  308 IF(IB.LE.NVC) OBS(II)=1.d0                                           
      GOTO 420                                                             
  357 if(ifcgs.eq.0) write(6,356)
      if(ib.gt.nvc.and.ib.le.nlvc) obs(ii)=-2.d0*br(ix)/xlog10e
      goto 420
  355 if(ifcgs.eq.0) write(6,354)
      iban=iband(ib)
      if(ib.gt.nvc.and.ib.le.nlvc) obs(ii)=-cgsta(ix)*.4d0*xlg10*atten*
     $ffband(iban)
      goto 420
  743 continue
c
c  After label 743, DC computes numeric-analytic hybrid derivatives,
c    df/d(tstart), df/d(tmax1), df/d(tmax2), or df/dt(final), for 
c    spot A.
c
      t1a=tstart(kspa,nspa)
      t2a=tmax1(kspa,nspa)
      t3a=tmax2(kspa,nspa)
      t4a=tfinal(kspa,nspa)
      obs(ii)=0.d0
      if(hjd.le.t1a.or.hjd.ge.t4a) goto 420
      if(hjd.ge.t2a.and.hjd.le.t3a) goto 420
c
c  NORSCALC=1 turns off the spot radius calculation in subroutine BBL.
c
      norscalc=1
c
c  Subroutine SPOTEV computes the spot radius at the time of observation.
c
      call spotev(hjd,t1a,t2a,t3a,t4a,amax(kspa,nspa),area,radspot)
      if(radspot.le.0.d0) goto 420
c
c  IDIDA can cause DC to skip computation of df/drs if that derivative
c    has already been computed (IDIDA=1).
c
      if(idida.eq.1.or.ib.gt.nlvc) goto 238
c
c  Now compute derivative df/d(rspot) numerically in two passes
c    thru subroutine BBL in the DO 221 loop.
c
      phasin=phas(ix)
      if(if3b.gt.0) phasin=phasin-deltat/period
      abdel=dabs(del(3))
c
c    Need to make sure the calls to BBL in the DO 221 and DO 222 loops
c    are carried out with the correct vfac3b value, transferred to BBL
c    via the common con3b block. If vfac3b was changed earlier
c    when calculating numerical derivatives with respect to a,
c    it may have acquired a value corresponding to an incremented value
c    of a. The next line resets vfac3b.
c
      vfac3b=vfac3bsav
      do 221 i221=1,2
      signdel=dfloat(2*i221-3)
      delrad=dfloat(i221-1)*abdel
      if(radspot.gt.abdel) delrad=.5d0*abdel*signdel 
      radinc(kspa,nspa)=radspot+delrad
      CALL BBL(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVH,FR1,FR2,HLD,
     $SLUMP1,SLUMP2,THETA,RHO,AA,BB,POT1,POT2,N1,N2,FF1,FF2,d,hla(ib),
     $cla(ib),x1,x2,y1,y2,g1,g2,wla(ib),sm1,sm2,tph,tpc,sbrh,sbrc,t1,
     $t2,a1,a2,xbol1,xbol2,ybol1,ybol2,fspot1,fspot2,phasin,rmass,xinc,
     $hot,cool,snthh,csthh,snfih,csfih,tldh,glump1,glump2,xx1,xx2,yy1,
     $yy2,zz1,zz2,dint1,dint2,grv1,grv2,rftemp,rf1,rf2,csbt1,csbt2,
     $gmag1,gmag2,glog1,glog2,obser,
     $fbin1,fbin2,delv1,delv2,count1,count2,delwl1,delwl2,resf1,resf2,
     $wl1,wl2,dvks1,dvks2,tau1,tau2,emm1,emm2,hbarw1,hbarw2,xcl,ycl,zcl,
     $rcl,op1,fcl,dens,encl,edens,taug,emmg,yskp,zskp,mode,iband(ib),
     $ifat1,ifat2,1)
      htt=hot
      if(mode.eq.-1) htt=0.d0
      total=(htt+cool+el3a(ib))*atten
      if(i221.eq.1) xr0=total*elit+vkm1*vc1+vkm2*vc2
      xrinc=total*elit+vkm1*vc1+vkm2*vc2
      if(ifcgs.eq.0) goto 221
      cgs1=aodsq*htt*polin1/sbrh
      cgs2=aodsq*cool*polin2/sbrc
      cgstot=(cgs1+cgs2+el3a(ib))*atten
      if(i221.eq.1) xr0=cgstot*elit+vkm1*vc1+vkm2*vc2
      xrinc=cgstot*elit+vkm1*vc1+vkm2*vc2
  221 continue
      dfdrs(ix)=(xrinc-xr0)/abdel
  238 continue
      call daspdt (hjd,amax(kspa,nspa),t1a,t2a,t3a,t4a,dadt1,dadt2,
     $dadt3,dadt4)
      drsdas=1.d0/(twopi*dsin(radspot))
      if(kh.eq.43) obs(ii)=dfdrs(ix)*drsdas*dadt1
      if(kh.eq.44) obs(ii)=dfdrs(ix)*drsdas*dadt2
      if(kh.eq.45) obs(ii)=dfdrs(ix)*drsdas*dadt3
      if(kh.eq.46) obs(ii)=dfdrs(ix)*drsdas*dadt4
      goto 420
  747 continue
c
c  After label 747, DC computes numeric-analytic hybrid derivatives,
c    df/d(tstart), df/d(tmax1), df/d(tmax2), or df/dt(final), for 
c    spot B.
c
      t1b=tstart(kspb,nspb)
      t2b=tmax1(kspb,nspb)
      t3b=tmax2(kspb,nspb)
      t4b=tfinal(kspb,nspb)
      obs(ii)=0.d0
      if(hjd.le.t1b.or.hjd.ge.t4b) goto 420
      if(hjd.ge.t2b.and.hjd.le.t3b) goto 420
c
c  NORSCALC=1 turns off the spot radius calculation in subroutine BBL.
c
      norscalc=1
c
c  Subroutine SPOTEV computes the spot radius at the time of observation.
c
      call spotev(hjd,t1b,t2b,t3b,t4b,amax(kspb,nspb),area,radspot)
      if(radspot.le.0.d0) goto 420
c
c  IDIDB can cause DC to skip computation of df/drs if that derivative
c    has already been computed (IDIDB=1).
c
      if(ididb.eq.1.or.ib.gt.nlvc) goto 239
c
c  Compute derivative df/drs numerically in two passes
c    thru subroutine BBL in the DO 222 loop.
c
      phasin=phas(ix)
      if(if3b.gt.0) phasin=phasin-deltat/period
      abdel=dabs(del(7))
      vfac3b=vfac3bsav
      do 222 i222=1,2
      signdel=dfloat(2*i222-3)
      delrad=dfloat(i222-1)*abdel
      if(radspot.gt.abdel) delrad=.5d0*abdel*signdel
      radinc(kspb,nspb)=radspot+delrad
      CALL BBL(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVH,FR1,FR2,HLD,
     $SLUMP1,SLUMP2,THETA,RHO,AA,BB,POT1,POT2,N1,N2,FF1,FF2,d,hla(ib),
     $cla(ib),x1,x2,y1,y2,g1,g2,wla(ib),sm1,sm2,tph,tpc,sbrh,sbrc,t1,
     $t2,a1,a2,xbol1,xbol2,ybol1,ybol2,fspot1,fspot2,phasin,rmass,xinc,
     $hot,cool,snthh,csthh,snfih,csfih,tldh,glump1,glump2,xx1,xx2,yy1,
     $yy2,zz1,zz2,dint1,dint2,grv1,grv2,rftemp,rf1,rf2,csbt1,csbt2,
     $gmag1,gmag2,glog1,glog2,obser,
     $fbin1,fbin2,delv1,delv2,count1,count2,delwl1,delwl2,resf1,resf2,
     $wl1,wl2,dvks1,dvks2,tau1,tau2,emm1,emm2,hbarw1,hbarw2,xcl,ycl,zcl,
     $rcl,op1,fcl,dens,encl,edens,taug,emmg,yskp,zskp,mode,iband(ib),
     $ifat1,ifat2,1)
      htt=hot
      if(mode.eq.-1) htt=0.d0
      total=(htt+cool+el3a(ib))*atten
      if(i222.eq.1) xr0=total*elit+vkm1*vc1+vkm2*vc2
      xrinc=total*elit+vkm1*vc1+vkm2*vc2
      if(ifcgs.eq.0) goto 222
      cgs1=aodsq*htt*polin1/sbrh
      cgs2=aodsq*cool*polin2/sbrc
      cgstot=(cgs1+cgs2+el3a(ib))*atten
      if(i222.eq.1) xr0=cgstot*elit+vkm1*vc1+vkm2*vc2
      xrinc=cgstot*elit+vkm1*vc1+vkm2*vc2
  222 continue
      dfdrs(ix)=(xrinc-xr0)/abdel
  239 continue
      call daspdt (hjd,amax(kspb,nspb),t1b,t2b,t3b,t4b,dadt1,dadt2,
     $dadt3,dadt4)
      drsdas=1.d0/(twopi*dsin(radspot))
      if(kh.eq.47) obs(ii)=dfdrs(ix)*drsdas*dadt1
      if(kh.eq.48) obs(ii)=dfdrs(ix)*drsdas*dadt2
      if(kh.eq.49) obs(ii)=dfdrs(ix)*drsdas*dadt3
      if(kh.eq.50) obs(ii)=dfdrs(ix)*drsdas*dadt4
      goto 420
c
c   Compute derivative of superior or inferior conjunction time with
c     respect to eccentricity
c
  761 delsp=1.d-8
      ecpl=ec+delsp
      call conjph(ecpl,argper,pshift,trsc,tric,econsc,econic,xmsc,xmic,
     $pconscpl,pconicpl)
      obs(ii)=period*(pconscpl-pconsc)/delsp
      if(mntype(ix).eq.2) obs(ii)=period*(pconicpl-pconic)/delsp
      goto 420
c
c   Compute derivative of superior or inferior conjunction time with
c     respect to argument of periastron parameter
c
  762 delsp=1.d-8
      argperpl=argper+delsp
      call conjph(ec,argperpl,pshift,trsc,tric,econsc,econic,xmsc,xmic,
     $pconscpl,pconicpl)
      obs(ii)=period*(pconscpl-pconsc)/delsp
      if(mntype(ix).eq.2) obs(ii)=period*(pconicpl-pconic)/delsp
      goto 420
  844 obs(ii)=1.d0
      if(ib.le.nlvc) obs(ii)=dfdph(ix)/period                                
      if(if3b.eq.0) goto 420
      obs2=obs(ii)
      if(ib.gt.nlvc) goto 768
      goto 420                                                             
  845 if(ib.gt.nlvc) goto 859
      brac=pzero+dtime*dpdt                                               
      obs(ii)=dfdph(ix)*dtime/(brac*pzero)                                
      goto 420
  859 continue
      call jdph(phjd(ix),0.d0,hjd0,pzero,dpdt,xjddum,phaz)
      kphaz=phaz
      pcon=dfloat(kphaz)+pconsc
      if(mntype(ix).eq.2) pcon=dfloat(kphaz)+pconic
      qq=pcon*dpdt
c
c  Here's the dt_conj/dP series expansion that's applied for small dpdt
c
      obs(ii)=pcon*(1.d0+.5d0*qq+qq*qq/6.d0+qq**3/24.d0+qq**4/120.d0+
     $qq**5/720.d0)
c
c  Here's the exact dt_conj/dP expression for larger dpdt
c
      if(dabs(dpdt).gt.1.d-6) obs(ii)=(dexp(qq)-1.d0)/dpdt
      obs2=obs(ii)
      if(if3b.ne.0) goto 766
      goto 420                                                             
  846 dis=dabs(dtime*dpdt/period)
      if(ib.gt.nlvc) goto 860
      if(dis.gt.toldis) goto 848
      brac2=2.d0*pzero+dtime*dpdt
      dphpd=-2.d0*(dtime/brac2)**2
     $+tt*dtime**3*(2.d0*brac2*dpdt-3.d0*dpdt**2*dtime)/brac2**4
      goto 849
  848 brac=pzero+dtime*dpdt
      dphpd=dtime/(brac*dpdt)-(dlog(brac)-dlog(pzero))/dpdt**2
  849 obs(ii)=-dfdph(ix)*dphpd                                             
      goto 420
  860 continue
      call jdph(phjd(ix),0.d0,hjd0,pzero,dpdt,xjddum,phaz)
      kphaz=phaz
      pcon=dfloat(kphaz)+pconsc
      if(mntype(ix).eq.2) pcon=dfloat(kphaz)+pconic
      if(dis.gt.toldis) goto 837
      obs(ii)=pcon*pzero*(.5d0*pcon+pcon*pcon*dpdt/3.d0+(pcon**3*
     $dpdt**2)/8.d0+(pcon**4*dpdt**3)/30.d0+(pcon**5*dpdt**4)/144.d0)
      obs2=obs(ii)
      if(if3b.ne.0) goto 7677
      goto 420
  837 continue
      obs(ii)=pzero*pcon*dexp(pcon*dpdt)/dpdt-(dexp(pcon*dpdt)-1.d0)*
     $pzero/dpdt**2
      obs2=obs(ii)
      if(if3b.ne.0) goto 7677
      goto 420
c
c   Compute derivative of superior or inferior conjunction time with
c     respect to d(omega)/dt parameter (omega is argument of periastron)
c
  847 if(ib.le.nlvc) goto 763
      delsp=1.d-8
      argperpl=argper+delsp  
      call conjph(ec,argperpl,pshift,trsc,tric,econsc,econic,xmsc,xmic,
     $pconscpl,pconicpl)
      obs(ii)=period*dtime*(pconscpl-pconsc)/delsp
      if(mntype(ix).eq.2) obs(ii)=period*dtime*(pconicpl-pconic)/delsp
      goto 420
  763 obs(ii)=dtime*dfdap(ix)                                              
      goto 420                                                             
c 
c   compute analytic derivatives for 3rd body light time and radial velocity 
c    with respect to a3b 
c 
  861 ddeltda3b=dtcon*sini3b*eupsfac*sinsum*(1.d0+2.d0*apfac) 
      if(ib.le.nlvc) goto 771
      obs(ii)=ddeltda3b
      goto 420
  771 dvda3b=rvfac*twopioppr*sini3b*(1.d0+2.d0*apfac)*eccfac/ertfac 
      term1=dfdph(ix)*ddeltda3b/period 
      term2=vc12*dvda3b 
      obs(ii)=term1+term2 
      goto 420 
c
c   compute analytic derivatives for 3rd body light time and radial velocity
c    with respect to a
c
  765 ddeltda=-3.d0*dtcon*a3b*sini3b*eupsfac*sinsum*apfac/a
      if(ib.le.nlvc) goto 777
      obs(ii)=ddeltda
      goto 420
  777 dvda=-3.d0*(a3b/a)*rvfac*twopioppr*sini3b*apfac*eccfac/ertfac
      term1=dfdph(ix)*ddeltda/period
      term2=vc12*dvda
      xmean3b=xm03b+twopi*(hjd-t03b)/p3b
      call kepler(xmean3b,e3b,ecan3b,tr3b)
      prat=p3b/period
      arat=a/a3b
      apfac=arat**3*prat**2
      apfacm=1.d0-apfac
      vfac3b=twopioppr*rvfac*a3b*sini3b*apfacm/ertfac
      v3b=vfac3b*(ecos3b+dcos(tr3b+perr3b))
      term3=vc12*(br(ix)-v3b)/a
      obs(ii)=term1+term2+term3
      goto 420
c
c   compute analytic derivatives for 3rd body light time and radial velocity
c    with respect to P0
c
  766 ddeltdp=2.d0*dtcon*a3b*sini3b*eupsfac*sinsum*apfac/period
      if(ib.le.nlvc) goto 778
      obs(ii)=obs2+ddeltdp
      goto 420
  778 dvdp=(2.d0/period)*a3b*rvfac*twopioppr*sini3b*apfac*eccfac/ertfac
      term1=dfdph(ix)*ddeltdp/period
      term2=vc12*dvdp
      obs(ii)=term1+term2
      goto 420
c
c   compute analytic derivatives for 3rd body light time and radial velocity
c    with respect to dp/dt
c
 7677 ddeltdpdot=dtime*2.d0*dtcon*a3b*sini3b*eupsfac*sinsum*apfac/period
      if(ib.le.nlvc) goto 779
      obs(ii)=obs2+ddeltdpdot
      goto 420
  779 dvdpdot=dtime*(2.d0/period)*a3b*rvfac*twopioppr*sini3b*apfac
      dvdpdot=dvdpdot*eccfac/ertfac
      term1=dfdph(ix)*ddeltdpdot/period
      term2=vc12*dvdpdot
      obs(ii)=term1+term2
      goto 420
c
c   compute analytic derivatives for 3rd body light time and radial velocity
c    with respect to T0
c
  768 ddeltdt0=2.d0*dtcon*a3b*sini3b*eupsfac*sinsum*apfac/period
      ddeltdt0=-ddeltdt0*dpdt
      if(ib.le.nlvc) goto 780
      obs(ii)=obs2+ddeltdt0
      goto 420
  780 dvdt0=(2.d0/period)*a3b*rvfac*twopioppr*sini3b*apfac*eccfac/ertfac
      dvdt0=-dvdt0*dpdt
      term1=dfdph(ix)*ddeltdt0/period
      term2=vc12*dvdt0
      obs(ii)=term1+term2
      goto 420
c 
c   compute analytic derivatives for 3rd body light time and radial velocity 
c    with respect to P3b 
c 
  862 pardeltdp3b=-2.d0*dtcon*a3b*sini3b*eupsfac*sinsum*arat**3*p3b/ 
     $pzero**2 
      ddeltdp3b=pardeltdp3b+pardeltdtr3b*dtr3bdp3b 
      if(ib.le.nlvc) goto 772
      obs(ii)=ddeltdp3b
      goto 420
  772 parvdp3b=-rvfac*twopi*a3b*sini3b*eccfac* 
     $(arat**3/pzero**2+1.d0/p3b**2)/ertfac 
      dvdp3b=parvdp3b+parvdtr3b*dtr3bdp3b 
      term1=dfdph(ix)*ddeltdp3b/period
      term2=vc12*dvdp3b 
      obs(ii)=term1+term2 
      goto 420 
c 
c   compute analytic derivatives for 3rd body light time and radial velocity 
c    with respect to i3b 
c 
  863 ddeltdi3b=dtr*dtcon*a3b*cosi3b*apfacm*eupsfac*sinsum 
      if(ib.le.nlvc) goto 773
      obs(ii)=ddeltdi3b
      goto 420
  773 dvdi3b=dtr*rvfac*twopioppr*a3b*cosi3b*apfacm* 
     $eccfac/ertfac 
      term1=dfdph(ix)*ddeltdi3b/period
      term2=vc12*dvdi3b 
      obs(ii)=term1+term2 
      goto 420 
c 
c   compute analytic derivatives for 3rd body light time and radial velocity 
c    with respect to e3b 
c 
  864 pardeltde3b=-dtcon*a3b*sini3b*apfacm*sinsum*(costr3b* 
     $(1.d0+e3b**2)+2.d0*e3b)/ecosp1**2 
      parvde3b=twopioppr*rvfac*a3b*sini3b*apfacm*(e3b*eccfac/ertfac**3 
     $+cosperr3b/ertfac) 
      ddeltde3b=pardeltde3b+pardeltdtr3b*dtrde3b 
      if(ib.le.nlvc) goto 774
      obs(ii)=ddeltde3b
      goto 420
  774 dvde3b=parvde3b+parvdtr3b*dtrde3b 
      term1=dfdph(ix)*ddeltde3b/period
      term2=vc12*dvde3b 
      obs(ii)=term1+term2 
      goto 420 
c 
c   compute analytic derivatives for 3rd body light time and radial velocity 
c    with respect to argument of periastron of 3rd body orbit 
c 
  865 pardeltdperr3b=dtcon*a3b*sini3b*apfacm*eupsfac*cossum 
      parvdperr3b=-rvfac*twopioppr*a3b*sini3b*apfacm* 
     $(e3b*sinperr3b+sinsum)/ertfac 
      ddeltdperr3b=pardeltdperr3b+pardeltdtr3b*dtr3bdperr3b
      if(ib.le.nlvc) goto 775
      obs(ii)=ddeltdperr3b
      goto 420
  775 dvdperr3b=parvdperr3b+parvdtr3b*dtr3bdperr3b
      term1=dfdph(ix)*ddeltdperr3b/period 
      term2=vc12*dvdperr3b 
      obs(ii)=term1+term2 
      goto 420 
c 
c   compute analytic derivatives for 3rd body light time and radial velocity 
c    with respect to T0 
c 
  867 dtr3bdt0=-twopioppr*partr3bdecan3b/(1.d0-e3b*dcos(ecan3b))
      ddeltdt0=pardeltdtr3b*dtr3bdt0 
      if(ib.le.nlvc) goto 776
      obs(ii)=ddeltdt0
      goto 420
  776 dvdt3b0=parvdtr3b*dtr3bdt0 
      term1=dfdph(ix)*ddeltdt0/period
      term2=vc12*dvdt3b0 
      obs(ii)=term1+term2 
      goto 420 
  301 if(ib.gt.nvc.and.ib.le.nlvc) obs(ii)=atten
  420 CONTINUE                                                             
  418 CONTINUE                                                             
  417 CONTINUE                                                             
  419 CONTINUE                                                             
c
c  Next compute the sum of squares (SS) of residuals over two specified sub-ranges of
c   phase (spresq) and also SS over all phases (resq) in 298, 299 loops.
c
c  Note that DC actually works with INTERNAL fluxes that are the input fluxes times
c   the calibration factor (CALIB) and times a flux unit (XUNIT). Or if input is in
c   magnitudes, DC figures fluxes from the line (somewhere above)
c     flux(inn)=calib(lcv)*10.d0**(-.4d0*xmag(inn))
c   In either case, DC is consistent in
c   computing the weighted sums of squares of residuals for the various bands from
c   the internal fluxes, and subsequently applying them TO the internal fluxes. So the 
c   idea is to use the box of 
c     "Sums of squares of residuals for separate curves, including only individual weights"
c   to calculate SIGMAS. Do not calculate SIGMA's
c   from the direct input fluxes unless XUNIT and CALIB are unity.
c
      write(6,*)
      write(6,*)
      write(6,*)
      write(6,159)
      write(6,*)
      write(6,169)
      do 298 icv=1,nlvct                                                    
      icvp=icv+1                                                           
      nbs=knobs(icvp)-knobs(icv)                                           
      iw=knobs(icv)+1                                                      
      jstart=nmat+iw                                                       
      jstop=jstart+nbs-1                                                   
      resq=0.d0                                                            
      spresq=0.d0
      kspresq=0
      iww=iw-1                                                             
      obssum=0.d0
      tmag=0.d0
      do 299 jres=jstart,jstop                                             
      iww=iww+1                                                            
      if(icv.gt.nlvc) goto 645
      sphas=phas(iww)
  598 if(sphas.le.1.d0) goto 599
      sphas=sphas-1.d0
      goto 598
  599 continue
      if(sphas.ge.0.d0) goto 597
      sphas=sphas+1.d0
      goto 599
  597 continue
      if(sphas.lt.sphas1(icv).or.sphas.gt.sphas2(icv)) goto 644
      spresq=spresq+wt(iww)*obs(jres)**2
      obssum=obssum+flux(iww)
      if(icv.gt.nvc) tmag=tmag+xmag(iww)
      kspresq=kspresq+1
  644 continue
      if(sphas.lt.sphas3(icv).or.sphas.gt.sphas4(icv)) goto 645
      spresq=spresq+wt(iww)*obs(jres)**2
      obssum=obssum+flux(iww)
      if(icv.gt.nvc) tmag=tmag+xmag(iww)
      kspresq=kspresq+1
  645 continue
  299 resq=resq+wt(iww)*obs(jres)**2                                       
      if(kspresq.gt.5.or.ksd(icv).eq.0) goto 400
      write(6,*)
      write(6,401) kspresq,icv,icv
      write(6,*)
  400 continue
  401 format('Only ',i2,' data points in specified ranges for curve',
     $i3,' Probably KSD should be set to 0 (zero) for curve',i3,' or the
     $ specified phase ranges made wider')
c
c  Calculate this curve's standard deviations in next 7 lines
c
      enminus=dfloat(nbs-1)
      if(enminus.gt.0.d0) sdevall=dsqrt(resq/enminus)
      if(ksd(icv).eq.1) sigma(icv)=sdevall
      if(icv.gt.nlvc) goto 760
      enminussp=dfloat(kspresq-1)
      if(enminussp.gt.0.d0) sdevsp=dsqrt(spresq/enminussp)
      if(ksd(icv).eq.2) sigma(icv)=sdevsp
  760 continue
c
c  Compute mean flux and mean magnitude within the specified phase range
c       (next 8 lines)
c
      obsmean=0.d0
      if(kspresq.ne.0) obsmean=obssum/dfloat(kspresq)
      if(kspresq.ne.0) avmag=tmag/dfloat(kspresq)
      refflux(icv)=obsmean
      if(maglite.ne.0) refflux(icv)=10.d0**(-.4d0*avmag)*calib(icv)
      if(icv.le.nvc) write(6,119) icv,nbs,sdevall,kspresq,sdevsp
      if(icv.gt.nvc.and.icv.le.nlvc) write(6,119) icv,nbs,sdevall,
     $kspresq,sdevsp,refflux(icv),avmag                                     
      if(icv.gt.nlvc) write(6,119) icv,nbs,sdevall
      if(ibef.eq.0) sigmainit(icv)=sigma(icv)
  298 continue
      write(6,*)                                                         
      do 909 komp=1,2  
      if(message(komp,1).eq.1) write(6,283) komp  
      if(message(komp,2).eq.1) write(6,284) komp  
      if(message(komp,3).eq.1) write(6,285) komp  
      if(message(komp,4).eq.1) write(6,286) komp  
  909 continue  
      write(6,*)                                                         
      GOTO 65                                                              
c
c   DC starts the next subset at label 71
c
   71 JF=0                                                                 
      do 404 ib=1,nlvct
  404 sigma(ib)=sigmainit(ib)
      niter=0
      do 793 kpar=1,60
      imax=1
      if(kpar.gt.55) imax=nlc
      if(keep(kpar).eq.1) goto 793
      do 794 kurv=1,imax
      kcurv=kurv
      if(kpar.le.55) kcurv=0
      ipar=kpar
      if(kpar.gt.55) ipar=55+kurv+(kpar-56)*nlc
      para(ipar)=parainit(ipar)
  794 continue
  793 continue
      IF(KO.EQ.0) stop                                                     
      if(key.eq.1) goto 204
      DO 261 J=1,NCOF                                                      
  261 OBS(J)=HOLD(J)                                                       
  204 continue
      IF(KDISK.EQ.0) GOTO 72                                               
      REWIND 9                                                             
      READ(9,67)(OBS(J),J=1,NCOF)                                          
   72 READ(5,20)(KEEP(I),I=1,60),niters,xlamda,vlr
      if(keep(1).ne.2) goto 866 
      close (5)  
      close (6)  
      if(iflcin.ne.1) stop
c 
c  The following default values can be changed in the generated input file
c    for the LC program
c
      stdev=0.d0
      seed=138472375.d0
      hjdst=hjd0
      hjdsp=hjdst+1.d0
      hjdin=0.1d0
      phstrt=-0.2d0
      phstop=+1.0d0
      phin=0.01d0
      phn=0.25d0
      phobs=0.75d0
      lspin=1
      tobs=tavh
      zero=+8.00d0
      factor=1.d0
      stop150=150.d0
      stoptime=-10000.d0
      kstop9=9
c
c  Here (to label 678) DC generates an input file for the LC program, with the parameter 
c    input values of the last DC iteration.
c
      open(unit=6,file='lcin.input_from_DC')
      do 678 iwrt=1,nlc+iftime+1
      mpage=2
      if(iwrt.gt.1) mpage=1
      if(iwrt.gt.nlc+1) mpage=6
      jdphw=2
      if(if3b.eq.1) jdphw=1
      if(iwrt.gt.nlc+1) jdphw=1
      iwrp=iwrt
      if(iwrt.gt.1) iwrp=iwrt+nvc-1
      write(6,22) mpage,nref,mref,ifsmv1,ifsmv2,icor1,icor2,if3b,ld1,
     $ld2,kspev,kspot,nomax,ifcgs             
      write(6,811) jdphw,hjd0,pzero,dpdt,pshift,delph,nga,stdev,
     $noise(iwrp),seed          
      write(6,217) hjdst,hjdsp,hjdin,phstrt,phstop,phin,phn,phobs,lspin,
     $tobs                
      write(6,568) MODE,IPB,IFAT1,IFAT2,N1,N2,perr0,dperdt,the,VUNIT         
      write(6,2) E,A,F1,F2,VGA,XINCL,GR1,GR2,abunin,fspot1,fspot2                               
      write(6,6) tavh,tavc,alb1,alb2,phsv,pcsv,rm,xbol1,xbol2,ybol1,       
     $ybol2,dpclog                                                            
      write(6,566) a3b,p3b,xinc3b,e3b,perr3b,t03b
      write(6,4)iband(iwrp),hla(iwrp),cla(iwrp),x1a(iwrp),x2a(iwrp),
     $y1a(iwrp),y2a(iwrp),el3a(iwrp),opsfa(iwrp),zero,factor,wla(iwrp),
     $aextinc(iwrp),calib(iwrp) 
      DO 478 KP=1,2                                                        
      nspot=nsp1
      if(kp.eq.2) nspot=nsp2
      DO 787 I=1,nspot+1                                                
      write(6,567)XLAT(KP,I),XLONG(KP,I),RADSP(KP,I),TEMSP(KP,I),            
     $tstart(kp,i),tmax1(kp,i),tmax2(kp,i),tfinal(kp,i)
  787 continue
  478 continue
      write(6,210) stop150
      if(mpage.ne.6) goto 678
      ist=knobs(nlvct)+1
      isp=knobs(nlvct+1)
      do 125 imn=ist,isp
  125 write(6,139) phjd(imn),mntype(imn)
      write(6,216) stoptime
  678 continue
      write(6,211) kstop9
      close(6)
      stop 
  866 continue 
      DO 232 I=1,ichno                                                        
  232 IF(KEP(I).EQ.1) KEEP(I)=1                                            
      NOBS=KNOBS(LC1)                                                      
      matrix=56                                                            
      do 428 ima=1,55                                                      
  428 matrix=matrix-keep(ima)                                              
      matrix=matrix+nlc*(5-keep(56)-keep(57)-keep(58)-keep(59)-keep(60))   
      MAT=MATRIX-1                                                         
      NCOEFF=MATRIX*NOBS                                                   
      KC=1                                                                 
      NSHIFT(1)=0                                                          
      DO 59 I=2,61                                                         
      IF(I.GT.56) KC=NLC                                                   
      KE=0                                                                 
      J=I-1                                                                
      IF(KEEP(J).GT.KEP(J)) KE=1                                           
   59 NSHIFT(I)=NOBS*KE*KC+NSHIFT(J)                                       
      NOBBS=NOBS                                                           
      DO 30 I=1,61                                                         
      IF(KEEP(I).EQ.1) GOTO 30                                             
      IF(I.GT.55) NOBBS=NOBS*NLC                                           
      IF(I.EQ.61) NOBBS=NOBS                                               
      DO 32 J=1,NOBBS                                                      
      JF=JF+1                                                              
      KX=JF+NSHIFT(I)                                                      
   32 OBS(JF)=OBS(KX)                                                      
   30 CONTINUE                                                             
   65 continue                                                         
      write(6,*)                                                         
      WRITE(6,20)(KEEP(I),I=1,60),niters,xlamda,vlr
      NOBS=KNOBS(LC1)                                                      
      WRITE(6,*)                                                         
      if(ifder.eq.0) goto 817
c
c  The next line causes writing of the unweighted derivatives and 
c    residuals to be skipped except for the first and last iterations
c    of a given solution.
c
      if(niter.ne.0.and.niter.ne.(niters-1)) goto 817
c
c  Write unweighted derivatives and residuals
c
      WRITE(6,82)                                                          
      WRITE(6,*)                                                         
      if(ifoc.eq.0) write(6,175)
      if(jdphs.eq.1.and.ifoc.eq.1) WRITE(6,182)
      if(jdphs.eq.2.and.ifoc.eq.1) WRITE(6,192)
      DO 96 IB=1,nlvct                                                      
      IST=KNOBS(IB)+1                                                      
      IB1=IB+1                                                             
      ISP=KNOBS(IB1)                                                       
      DO 96 I=IST,ISP                                                      
      phasn=phas(i)
 1488 continue
      if(phasn.gt.1.d0) phasn=phasn-1.d0
      if(phasn.gt.1.d0) goto 1488
 1493 continue
      if(phasn.lt.0.d0) phasn=phasn+1.d0
      if(phasn.lt.0.d0) goto 1493
      obsrv=flux(i)
      cmptd=br(i)
      if(ib.gt.nlvc) obsrv=phjd(i)
      if(ib.gt.nlvc) cmptd=timecomp(i)
      if(jdphs.eq.1.and.ifoc.eq.0) write(6,147) phjd(i),phasn,
     $(obs(j),j=i,ncoeff,nobs)
      if(jdphs.eq.1.and.ifoc.eq.1) write(6,147) phjd(i),phasn,
     $obsrv,cmptd,(obs(j),j=i,ncoeff,nobs)
      if(jdphs.eq.2.and.ifoc.eq.0) write(6,149) phasn,
     $(obs(j),j=i,ncoeff,nobs)
      if(jdphs.eq.2.and.ifoc.eq.1) write(6,149) phasn,obsrv,cmptd,
     $(obs(j),j=i,ncoeff,nobs)
   96 CONTINUE                                                             
  817 continue
      IF(KO.LE.1) GOTO 70                                                  
      IF(IBEF.EQ.1) GOTO 70                                                
      DO 62 J=1,NCOEFF                                                     
   62 HOLD(J)=OBS(J)                                                       
      IF(KDISK.EQ.0) GOTO 73                                               
      REWIND 9                                                             
      WRITE(9,67)(OBS(J),J=1,NCOEFF)                                       
   73 CONTINUE                                                             
   70 WRITE(6,*)                                                         
c
c  Compute and apply curve-dependent, level-dependent, and assigned weights
c
      sumwt=0.d0
      DO 97 IB=1,nlvct                                                      
      IST=KNOBS(IB)+1                                                      
      IB1=IB+1                                                             
      ISP=KNOBS(IB1)                                                       
      NOIS=NOISE(IB)                                                       
      DO 97 I=IST,ISP                                                      
      IF(IB.GT.NVC.and.ib.le.nlvc) GOTO 444                                               
      ROOTWT=dsqrt(WT(I))/SIGMA(IB)                               
      GOTO 445                                                             
c
c  SIGMA for light curves applies at the reference level, so the effective SIGMA
c   at another level is SIGMA*dsqrt(FLUX/REFFLUX)**NOIS.
c
  444 ROOTWT=dsqrt(WT(I))/(SIGMA(IB)*dsqrt(FLUX(I)/refflux(ib))**NOIS)
  445 CONTINUE                                                             
      sumwt=sumwt+rootwt**2
      DO 97 LOB=I,NCOEFF,NOBS                                              
   97 OBS(LOB)=OBS(LOB)*ROOTWT                                             
      if(ibef.eq.0) sumwthold=sumwt
      if(ifder.eq.0) goto 818
c
c  The next line cause writing of the weighted derivatives and 
c    residuals to be skipped except for the first and last iterations
c    of a given solution.
c
      if(niter.ne.0.and.niter.ne.(niters-1)) goto 818
c
c  Write weighted derivatives and residuals
c
      WRITE(6,83)                                           
      WRITE(6,*)                                          
      if(jdphs.eq.1) WRITE(6,175)
      if(jdphs.eq.2) WRITE(6,176)
      DO 98 I=1,NOBS                                                       
      phasn=phas(i)
 1492 continue
      if(phasn.gt.1.d0) phasn=phasn-1.d0
      if(phasn.gt.1.d0) goto 1492
 1491 continue
      if(phasn.lt.0.d0) phasn=phasn+1.d0
      if(phasn.lt.0.d0) goto 1491
      if(jdphs.eq.1) write(6,147) phjd(i),phasn,(obs(j),j=i,ncoeff,nobs)
      if(jdphs.eq.2) write(6,149) phasn,(obs(j),j=i,ncoeff,nobs)
   98 CONTINUE                                                             
  818 continue
      CALL square (OBS,NOBS,MAT,OUT,sd,xlamda,deter,CN,CNN,cnc,clc,         
     $pred,ccl,ll,mm)                                                         
      MSQ=MAT*MAT                                                          
      WRITE(6,*)                                                         
      WRITE(6,181)                                                         
      WRITE(6,*)                                                         
      DO 38 JR=1,MAT                                                       
   38 WRITE(6,37) (CN(JX),JX=JR,MSQ,MAT),CCL(JR)                           
      WRITE(6,*)                                                         
      WRITE(6,183)                                                         
      WRITE(6,*)                                                         
      NO1=23                                                               
      NO2=24                                                               
      NRM=25                                                               
      DO 334 IRM=1,24                                                      
      IF(IRM.LE.23) NO1=NO1-KEEP(IRM)                                      
      NO2=NO2-KEEP(IRM)                                                    
  334 NRM=NRM-KEEP(IRM)                                                    
      CORO1=1-KEEP(23)                                                     
      CORO2=1-KEEP(24)                                                     
      CORQ=1-KEEP(25)                                                      
      co1q=0.d0                                                            
      co2q=0.d0                                                            
      DO 34 JM=1,MAT                                                       
      DO 33 JQ=1,MAT                                                       
      JT=JM+MAT*(JQ-1)                                                     
      IJM=(MAT+1)*(JM-1)+1                                                 
      IJQ=(MAT+1)*(JQ-1)+1                                                 
   33 V(JQ)=CNN(JT)/DSQRT(CNN(IJM)*CNN(IJQ))                               
      if(jm.eq.nrm.and.no1.gt.0) co1q=v(no1)*corq*coro1                    
      if(jm.eq.nrm.and.no2.gt.0) co2q=v(no2)*corq*coro2                    
   34 WRITE(6,37)(V(IM),IM=1,MAT)                                          
      WRITE(6,*)                                                         
      WRITE(6,184)                                                         
      WRITE(6,*)                                                         
      CALL DGMPRD(CN,CNN,CNOUT,MAT,MAT,MAT)                                
      DO 116 J8=1,MAT                                                      
  116 WRITE(6,37)(CNOUT(J7),J7=J8,MSQ,MAT)                                 
      WRITE(6,*)                                                         
      WRITE(6,185)                                                         
      WRITE(6,*)                                                         
      ANSCH=0.D0                                                           
      DO 118 J5=1,MAT                                                      
      V(J5)=0.D0                                                           
      DO 117 J6=1,MAT                                                      
      int=mat*(j6-1)                                                       
      idi=j6+int                                                           
      I9=J5+int                                                            
  117 V(J5)=OUT(J6)*CN(I9)*dsqrt(cnc(idi))+V(J5)                           
      ERR=V(J5)-CCL(J5)                                                    
  118 ANSCH=ANSCH+DABS(ERR)                                                
      WRITE(6,137)(V(J4),J4=1,MAT)                                         
      WRITE(6,*)                                                         
      WRITE(6,138) ANSCH                                                   
      WRITE(6,*)                                                         
      WRITE(6,*)                                                         
c
c  Write input-output in F format
c
      write(6,715)                                                         
      WRITE(6,*)                                                         
      write(6,43)                                                          
      iout=0                                                               
      do 93 kpar=1,60                                                      
      if(keep(kpar).eq.1) goto 93                                          
      imax=1                                                               
      if(kpar.gt.55) imax=nlc                                              
      do 94 kurv=1,imax                                                    
      kcurv=kurv                                                           
      if(kpar.le.55) kcurv=0                                               
      iout=iout+1                                                          
      ipar=kpar                                                            
      if(kpar.gt.55) ipar=55+kurv+(kpar-56)*nlc                            
      parout=para(ipar)+out(iout)*vlr                                          
      if(kpar.ne.56) write(6,615) kpar,kcurv,para(ipar),out(iout),
     $parout,sd(iout)         
      if(kpar.ne.56) goto 94
      kpr=nvc+kurv
      sdrat=sd(iout)/(hla(kpr)+cla(kpr))
      write(6,615) kpar,kcurv,para(ipar),out(iout),parout,sd(iout),
     $ratlum(kpr),sdrat
   94 continue                                                             
   93 continue                                                             
      WRITE(6,*)                                                         
      WRITE(6,*)                                                         
c
c  Write input-output in D format
c
      write(6,716)                                                         
      WRITE(6,*)                                                         
      write(6,43)                                                          
      iout=0                                                               
      do 53 kpar=1,60                                                      
      if(keep(kpar).eq.1) goto 53                                          
      imax=1                                                               
      if(kpar.gt.55) imax=nlc                                              
      do 52 kurv=1,imax                                                    
      kcurv=kurv                                                           
      if(kpar.le.55) kcurv=0                                               
      iout=iout+1                                                          
      ipar=kpar                                                            
      if(kpar.gt.55) ipar=55+kurv+(kpar-56)*nlc                            
      parout=para(ipar)+out(iout)*vlr                                          
      if(kpar.ne.56) write(6,616) kpar,kcurv,para(ipar),out(iout),
     $parout,sd(iout)         
      if(kpar.ne.56) goto 52
      kpr=nvc+kurv
      sdrat=sd(iout)/(hla(kpr)+cla(kpr))
      write(6,616) kpar,kcurv,para(ipar),out(iout),parout,sd(iout),
     $ratlum(kpr),sdrat
   52 para(ipar)=parout
   53 continue                                                             
      WRITE(6,*)                                                         
      RESSQ=0.d0                                                           
      JST=MAT*NOBS+1                                                       
      DO 199 JRES=JST,NCOEFF                                               
  199 RESSQ=RESSQ+OBS(JRES)**2                                             
      divisor=sumwt*dfloat(nobs)
      if(key.eq.0) divisor=sumwthold*dfloat(nobs)
      resmean=dsqrt(ressq/divisor)
      spred=dsqrt(pred/divisor)
      WRITE(6,*)                                                         
      WRITE(6,40)                                                          
      WRITE(6,21) resmean,spred,deter                                              
      IBEF=1                                                               
      write(6,*)
      niter=niter+1
      write(6,3333) niter,niters
      WRITE(6,*)                                                         
      WRITE(6,*)                                                         
      WRITE(6,650)                                                         
      WRITE(6,*)                                                         
      WRITE(6,*)                                                         
      WRITE(6,653)                                                         
      WRITE(6,*)                                                         
      do1=0.d0
      if(no1.ne.0) do1=sd(no1)*coro1
      do2=0.d0
      if(no2.ne.0) do2=sd(no2)*coro2
      if(mod.eq.1) do2=do1                                                 
      if(mod.eq.1) co2q=co1q                                               
      DQ=sd(NRM)*CORQ                                                      
      COQ=CO1Q                                                             
      ecc=para(10)
      dp=1.d0-ecc                                                            
      ome=para(23)
      rmas=para(25)
      fro=para(12)
      DOM=DO1                                                              
      KOMP=0                                                               
  925 CONTINUE                                                             
      KOMP=KOMP+1                                                          
      if(komp.eq.2) ome=para(24)
      if(komp.eq.2) fro=para(13)
      DO 926 KD=1,4                                                        
      if(kd.ne.2) goto 928                                                 
      if(po(komp).ge.omcr(komp)) goto 928                                  
      goto 926                                                             
  928 continue                                                             
      TH=XTHA(KD)                                                          
      FI=XFIA(KD)                                                          
      CALL ROMQ(ome,rmas,fro,dp,ecc,TH,FI,R,DRDO,DRDQ,DODQ,KOMP,MODE)            
      DR=dsqrt(DRDQ**2*DQ**2+DRDO**2*DOM**2+2.d0*COQ*DRDQ*DRDO*DQ*DOM)     
      WRITE(6,654)KOMP,ARAD(KD),R,DRDO,DRDQ,DR                             
  926 CONTINUE                                                             
      DO2DQ=DODQ                                                           
      IF(KOMP.EQ.1)DO1DQ=DODQ                                              
      COQ=CO2Q                                                             
      DOM=DO2                                                              
      WRITE(6,*)                                                         
      IF(KOMP.EQ.1) GOTO 925                                               
      WRITE(6,*)                                                         
      WRITE(6,651)                                                         
      IF(KOMP.EQ.2) WRITE(6,652)DO1DQ,DO2DQ,CO1Q,CO2Q,DO1,DO2,DQ
c           
c  Next 2 lines: key=0 means that OBS is to be read from the storage array
c   called HOLD. That's so the first iteration of the next subset
c   "comes free" (don't have to re-compute the partial derivatives and
c    residuals - they've already been computed for the base set). So
c    if the last iteration for the currently active subset was just 
c    finished, DC goes to label 71 and reads the saved initial parameter values
c    back into the para array and starts the new subset solution with the 
c    "free" derivatives and residuals.
c
      key=0
      if(niter.ge.niters) goto 71
c
c  Next 2 lines: key=1 means that the OBS array is NOT to be read from
c    HOLD, but re-computed from scratch, which is done on the 2nd, 3rd,
c    etc. iterations for the new subset. Go back to label 209 and copy
c    the current solution iteration results from the para array into
c    the individual parameter quantities to set up the next iteration.
c
      key=1
      goto 209
      END                                                                  
      SUBROUTINE SINCOS (KOMP,N,N1,SNTH,CSTH,SNFI,CSFI,MMSAVE)           
c  Version of November 9, 1995                                           
      implicit real*8 (a-h,o-z)                                          
      DIMENSION SNTH(*),CSTH(*),SNFI(*),CSFI(*),MMSAVE(*)                
      IP=(KOMP-1)*(N1+1)+1                                               
      IQ=IP-1                                                            
      IS=0                                                               
      IF(KOMP.EQ.2) IS=MMSAVE(IQ)                                        
      MMSAVE(IP)=0                                                       
      EN=N                                                               
      DO 8 I=1,N                                                         
      EYE=I                                                              
      EYE=EYE-.5d0                                                       
      TH=1.570796326794897d0*EYE/EN                                      
      IPN1=I+N1*(KOMP-1)                                                 
      SNTH(IPN1)=dsin(TH)                                                
      CSTH(IPN1)=dcos(TH)                                                
      EM=SNTH(IPN1)*EN*1.3d0                                             
      MM=EM+1.d0                                                         
      XM=MM                                                              
      IP=(KOMP-1)*(N1+1)+I+1                                             
      IQ=IP-1                                                            
      MMSAVE(IP)=MMSAVE(IQ)+MM                                           
      DO 8 J=1,MM                                                        
      IS=IS+1                                                            
      XJ=J                                                               
      FI=3.141592653589793d0*(XJ-.5d0)/XM                                
      CSFI(IS)=dcos(FI)                                                  
      SNFI(IS)=dsin(FI)                                                  
    8 CONTINUE                                                           
      RETURN                                                             
      END                                                                
      SUBROUTINE DURA(F,XINCL,RM,D,THE,OMEG,R)                           
c  Version of May 19, 1996                                               
C                                                                        
C     PARAMETER 'THE' IS THE SEMI-DURATION OF X-RAY ECLIPSE, AND SHOULD  
C     BE IN CIRCULAR MEASURE.                                            
      IMPLICIT REAL*8(A-H,O-Z)                                           
      DELX=0.D0                                                          
      FSQ=F*F                                                            
      RMD=1.d0/RM                                                        
      RMD1=RMD+1.D0                                                      
      XINC=.017453293d0*XINCL                                            
      TH=6.2831853071795865d0*THE                                        
      CI=DCOS(XINC)                                                      
      SI=DSIN(XINC)                                                      
      DSQ=D*D                                                            
      ST=DSIN(TH)                                                        
      CT=DCOS(TH)                                                        
      COTI=CI/SI                                                         
      TT=ST/CT                                                           
      C1=CT*SI                                                           
      C2=TT*ST*SI                                                        
      C3=C1+C2                                                           
      C4=COTI*CI/CT                                                      
      C5=C3+C4                                                           
      C6=C2+C4                                                           
      C7=(ST*ST+COTI*COTI)/CT**2                                         
      X=D*(SI*SI*ST*ST+CI*CI)+.00001D0                                   
   15 X=X+DELX                                                           
      PAR=X*X+C7*(D-X)**2                                                
      RPAR=DSQRT(PAR)                                                    
      PAR32=PAR*RPAR                                                     
      PAR52=PAR*PAR32                                                    
      FC=(C6*D-C5*X)/PAR32+C1**3*C5*RMD/(D-X)**2+C3*FSQ*RMD1*X-C2*FSQ*D* 
     $RMD1-C1*RMD/DSQ                                                    
      DFCDX=(-C5*PAR-3.D0*(C6*D-C5*X)*((1.D0+C7)*X-C7*D))/PAR52+2.D0*C1  
     $**3*C5*RMD/(D-X)**3+C3*FSQ*RMD1                                    
      DELX=-FC/DFCDX                                                     
      ABDELX=DABS(DELX)                                                  
      IF(ABDELX.GT.1.d-7) GOTO 15                                        
      Y=-(D-X)*TT                                                        
      Z=-(D-X)*COTI/CT                                                   
      YZ2=Y*Y+Z*Z                                                        
      OMEG=1.D0/DSQRT(X*X+YZ2)+RMD/DSQRT((D-X)**2+YZ2)+.5D0*RMD1*FSQ*    
     $(X*X+Y*Y)-RMD*X/DSQ                                                
      OMEG=RM*OMEG+.5d0*(1.d0-RM)                                        
      R=DSQRT(X*X+YZ2)                                                   
      RETURN                                                             
      END                                                                
      SUBROUTINE RING(Q,OM,KOMP,L,FR,HLD,R1,RL)                          
c   Version of September 14, 1998                                        
      IMPLICIT REAL*8(A-H,O-Z)                                           
      DIMENSION RAD(100),THET(100),AA(3),BB(3),FI(150),THA(150),FR(*),   
     $HLD(*)                                                             
      IX=0                                                               
      LR=L+1                                                             
      DO 92 I=1,LR                                                       
      THA(I)=0.D0                                                        
   92 FI(I)=-.1D0                                                        
      OMEGA=OM                                                           
      K=3                                                                
      EL=dfloat(L)                                                       
      DEL=2.D0/EL                                                        
      CALL ELLONE(1.d0,1.d0,Q,xlsv,OM1,XL2,OM2)                          
      CALL NEKMIN(Q,OM,xlsv,Z)                                           
      XL=xlsv                                                            
      QQ=Q                                                               
      XLSQ=XL*XL                                                         
      IF(Q.GT.1.D0) QQ=1.D0/Q                                            
      RMAX=DEXP(.345D0*DLOG(QQ)-1.125D0)                                 
      R=RMAX*(OM1-OMEGA)/(OM1-OM2)                                       
      DO 22 IT=1,L                                                       
      EYT=dfloat(IT)                                                     
      TH=EYT*1.570796326794897d0/EL                                      
      COSQ=DCOS(TH)**2                                                   
      DELR=0.D0                                                          
   14 R=DABS(R+DELR)                                                     
      RSQ=R*R                                                            
      X2R2=XLSQ+RSQ                                                      
      RX2R2=DSQRT(X2R2)                                                  
      XM2R2=(XL-1.D0)**2+RSQ                                             
      RXM2R2=DSQRT(XM2R2)                                                
      OM=1.D0/RX2R2+Q*(1.D0/RXM2R2-XL)+.5D0*(Q+1.D0)*(XLSQ+RSQ*COSQ)     
      DOMDR=-R/(X2R2*RX2R2)-Q*R/(XM2R2*RXM2R2)+(Q+1.D0)*COSQ*R           
      DELR=(OMEGA-OM)/DOMDR                                              
      ABDELR=DABS(DELR)                                                  
      IF(ABDELR.GT..00001D0) GOTO 14                                     
      RAD(IT)=R                                                          
   22 THET(IT)=TH*4.D0                                                   
      R1=RAD(1)                                                          
      RL=RAD(L)                                                          
      R90SQ=RL*RL                                                        
      DO 18 IJ=1,L                                                       
      EYJ=IJ                                                             
      RAD(IJ)=RAD(IJ)-(RL-R1)*(EYJ-1.D0)/EL                              
   18 CONTINUE                                                           
      DO 65 N=1,K                                                        
      AA(N)=0.D0                                                         
   65 BB(N)=0.D0                                                         
      DO 29 J=1,L                                                        
      DO 29 N=1,K                                                        
      EN=N-1                                                             
      ENTHET=EN*THET(J)                                                  
      AA(N)=AA(N)+RAD(J)*DCOS(ENTHET)*DEL                                
   29 BB(N)=BB(N)+RAD(J)*DSIN(ENTHET)*DEL                                
      AA(1)=AA(1)*.5D0                                                   
      IF(KOMP.EQ.2) XL=1.d0-xlsv                                         
      XLSQ=XL*XL                                                         
      DIS=RL/XL-.0005D0                                                  
      DO 42 IR=1,L                                                       
      LL=IR-1                                                            
      EY=dfloat(L+1-IR)                                                  
      THA(IR)=1.570796326794897D0*EY/EL                                  
      IF(THA(IR).LT.1.570796326794897D0) GOTO 82                         
      COT=0.D0                                                           
      GOTO 83                                                            
   82 COT=1.D0/DTAN(THA(IR))                                             
   83 IF(COT.GE.DIS) GOTO 50                                             
      COSSQ=DCOS(THA(IR))**2                                             
      A0=AA(1)                                                           
      A1=AA(2)                                                           
      A2=AA(3)                                                           
      B1=BB(2)                                                           
      B2=BB(3)                                                           
      DELSIN=0.D0                                                        
      KNTR=0                                                             
      SINTH=DSQRT(COSSQ*(XLSQ+R90SQ)/R90SQ)                              
   88 SINTH=SINTH+DELSIN                                                 
      KNTR=KNTR+1                                                        
      IF(SINTH.GT.1.D0) SINTH=1.D0/SINTH                                 
      CSQ=1.D0-SINTH*SINTH                                               
      COSTH=DSQRT(CSQ)                                                   
      SINSQ=SINTH*SINTH                                                  
      SIN4=8.D0*COSTH**3*SINTH-4.D0*COSTH*SINTH                          
      COS4=8.D0*CSQ*(CSQ-1.D0)+1.D0                                      
      C4SQ=COS4*COS4                                                     
      SINCOS=SIN4*COS4                                                   
      RRR=A0+A1*COS4+A2*(C4SQ+C4SQ-1.D0)+B1*SIN4+(B2+B2)*SINCOS          
      ARC=dasin(SINTH)                                                   
      RR=RRR+(RL-R1)*(2.D0*ARC/3.141592653589793D0-1.D0/EL)              
      IF(KNTR.GT.30) GOTO 42                                             
      P=RR*SINTH                                                         
      DRDSIN=-A1*SINTH/COSTH-4.D0*A2*SINTH+B1-(B2+B2)*SINSQ/COSTH+(B2+B2 
     $)*COSTH+(RL+RL-R1-R1)/(3.141592653589793D0*COSTH)                  
      DPDSIN=RR+SINTH*DRDSIN                                             
      F=P*P/COSSQ-RR*RR-XLSQ                                             
      DFDSIN=(P+P)*DPDSIN/COSSQ-(RR+RR)*DRDSIN                           
      DELSIN=-F/DFDSIN                                                   
      ABDEL=DABS(DELSIN)                                                 
      IF(ABDEL.GT..00001D0)  GOTO 88                                     
   42 FI(IR)=DATAN(RR*COSTH/XL)                                          
   50 LL1=LL+1                                                           
      DELTH=1.570796326794897D0/EL                                       
      DO 75 I=1,L                                                        
      EY=dfloat(L+1-I)-.5d0                                              
      THE=1.570796326794897D0*EY/EL                                      
      SNTH=DSIN(THE)                                                     
      EM=dsin(THE)*EL*1.3d0                                              
      MM=EM+1.d0                                                         
      XM=dfloat(MM)                                                      
      DELFI=3.141592653589793D0/XM                                       
      HDELFI=1.570796326794897D0/XM                                      
      DO 75 J=1,MM                                                       
      IX=IX+1                                                            
      IF(I.LE.LL1) GOTO 43                                               
      HLD(IX)=1.d0                                                       
      GOTO 75                                                            
   43 XJ=MM+1-J                                                          
      FE=3.141592653589793D0*(XJ-.5D0)/XM                                
      PH2=FE+HDELFI                                                      
      PHB=PH2                                                            
      IF(FI(I).GT.(FE-HDELFI)) GOTO 51                                   
      HLD(IX)=1.d0                                                       
      GOTO 75                                                            
   51 IPL=I+1                                                            
      IF(FI(IPL).GT.0.D0) GOTO 66                                        
      RR=A0+A1-A2+(RL-R1)*(1.D0-1.D0/EL)                                 
      PH1=DELFI*(XJ-1.D0)                                                
      TH1=DATAN(XL/RR)                                                   
      GOTO 56                                                            
   66 IF(FI(IPL).LT.(FE+HDELFI)) GOTO 52                                 
      HLD(IX)=0.d0                                                       
      GOTO 75                                                            
   52 IF(FI(IPL).LT.(FE-HDELFI)) GOTO 53                                 
      PH1=FI(IPL)                                                        
      TH1=THA(IPL)                                                       
      GOTO 56                                                            
   53 DELSIN=0.D0                                                        
      SINTH=DSQRT(COSSQ*(XLSQ+R90SQ)/R90SQ)                              
      TANFE=DTAN(FE-HDELFI)                                              
   77 SINTH=SINTH+DELSIN                                                 
      IF(SINTH.GT.1.D0) SINTH=1.D0/SINTH                                 
      SINSQ=SINTH*SINTH                                                  
      CSQ=1.D0-SINSQ                                                     
      COSTH=DSQRT(CSQ)                                                   
      SIN4=8.D0*COSTH**3*SINTH-4.D0*COSTH*SINTH                          
      COS4=8.D0*CSQ*(CSQ-1.D0)+1.D0                                      
      C4SQ=COS4*COS4                                                     
      SINCOS=SIN4*COS4                                                   
      RRR=A0+A1*COS4+A2*(C4SQ+C4SQ-1.D0)+B1*SIN4+(B2+B2)*SINCOS          
      ARC=dasin(SINTH)                                                   
      RR=RRR+(RL-R1)*(2.D0*ARC/3.141592653589793D0-1.D0/EL)              
      DRDSIN=-A1*SINTH/COSTH-4.D0*A2*SINTH+B1-(B2+B2)*SINSQ/COSTH+(B2+B2 
     $)*COSTH+(RL+RL-R1-R1)/(3.141592653589793D0*COSTH)                  
      F=RR*COSTH-XL*TANFE                                                
      DFDSIN=COSTH*DRDSIN-RR*SINTH/COSTH                                 
      DELSIN=-F/DFDSIN                                                   
      ABDEL=DABS(DELSIN)                                                 
      IF(ABDEL.GT..00001D0)  GOTO 77                                     
      PH1=FE-HDELFI                                                      
      TH1=DATAN(XL/(RR*SINTH*DCOS(PH1)))                                 
   56 IF(FI(I).GT.(FE+HDELFI)) GOTO 57                                   
      PHB=FI(I)                                                          
      TH2=THA(I)                                                         
      GOTO 60                                                            
   57 DELSIN=0.D0                                                        
      SINTH=DSQRT(COSSQ*(XLSQ+R90SQ)/R90SQ)                              
      TANFE=DTAN(FE+HDELFI)                                              
   78 SINTH=SINTH+DELSIN                                                 
      IF(SINTH.GT.1.D0) SINTH=1.D0/SINTH                                 
      SINSQ=SINTH*SINTH                                                  
      CSQ=1.D0-SINSQ                                                     
      COSTH=DSQRT(CSQ)                                                   
      SIN4=8.D0*COSTH**3*SINTH-4.D0*COSTH*SINTH                          
      COS4=8.D0*CSQ*(CSQ-1.D0)+1.D0                                      
      C4SQ=COS4*COS4                                                     
      SINCOS=SIN4*COS4                                                   
      RRR=A0+A1*COS4+A2*(C4SQ+C4SQ-1.D0)+B1*SIN4+(B2+B2)*SINCOS          
      ARC=dasin(SINTH)                                                   
      RR=RRR+(RL-R1)*(2.D0*ARC/3.141592653589793D0-1.D0/EL)              
      DRDSIN=-A1*SINTH/COSTH-4.D0*A2*SINTH+B1-(B2+B2)*SINSQ/COSTH+(B2+B2 
     $)*COSTH+(RL+RL-R1-R1)/(3.141592653589793D0*COSTH)                  
      F=RR*COSTH-XL*TANFE                                                
      DFDSIN=COSTH*DRDSIN-RR*SINTH/COSTH                                 
      DELSIN=-F/DFDSIN                                                   
      ABDEL=DABS(DELSIN)                                                 
      IF(ABDEL.GT..00001D0)  GOTO 78                                     
      TH2=DATAN(XL/(RR*SINTH*DCOS(PH2)))                                 
   60 CTHT=DCOS(THA(IPL))                                                
      CTH1=DCOS(TH1)                                                     
      CTH2=DCOS(TH2)                                                     
      STH1=DSIN(TH1)                                                     
      STH2=DSIN(TH2)                                                     
      DTH=TH2-TH1                                                        
      DCTH=CTH1-CTH2                                                     
      OMDP=PH2*DCTH-.5D0*(PH1*STH1+PHB*STH2)*DTH                         
      OMP=DELFI*(CTHT-CTH1)                                              
      OMN=OMP+OMDP                                                       
      HLD(IX)=OMN/(DELTH*DELFI*SNTH)                                     
   75 CONTINUE                                                           
      DO 94 JB=1,IX                                                      
      JA=IX+1-JB                                                         
   94 FR(JB)=HLD(JA)                                                     
      RETURN                                                             
      END                                                                
      SUBROUTINE DGMPRD(A,B,R,N,M,L)                                     
c  Version of April 9, 1992                                              
      DIMENSION A(*),B(*),R(*)                                           
      DOUBLE PRECISION A,B,R                                             
      IR=0                                                               
      IK=-M                                                              
      DO 10 K=1,L                                                        
      IK=IK+M                                                            
      DO 10 J=1,N                                                        
      IR=IR+1                                                            
      JI=J-N                                                             
      IB=IK                                                              
      R(IR)=0.D0                                                         
      DO 10 I=1,M                                                        
      JI=JI+N                                                            
      IB=IB+1                                                            
   10 R(IR)=R(IR)+A(JI)*B(IB)                                            
      RETURN                                                             
      END                                                                
      SUBROUTINE DMINV(A,N,D,L,M)                                        
c  Version of January 9, 2002                                              
      DIMENSION A(*),L(*),M(*)                                           
      DOUBLE PRECISION A,D,BIGA,HOLD                                     
      D=1.D0                                                             
      NK=-N                                                              
      DO 80 K=1,N                                                        
      NK=NK+N                                                            
      L(K)=K                                                             
      M(K)=K                                                             
      KK=NK+K                                                            
      BIGA=A(KK)                                                         
      DO 20 J=K,N                                                        
      IZ=N*(J-1)                                                         
      DO 20 I=K,N                                                        
      IJ=IZ+I                                                            
      IF(DABS(BIGA).GE.DABS(A(IJ))) GOTO 20                              
      BIGA=A(IJ)                                                         
      L(K)=I                                                             
      M(K)=J                                                             
   20 CONTINUE                                                           
      J=L(K)                                                             
      IF(J.LE.K) GOTO 35                                                 
      KI=K-N                                                             
      DO 30 I=1,N                                                        
      KI=KI+N                                                            
      HOLD=-A(KI)                                                        
      JI=KI-K+J                                                          
      A(KI)=A(JI)                                                        
   30 A(JI) =HOLD                                                        
   35 I=M(K)                                                             
      IF(I.LE.K) GOTO 45                                                 
      JP=N*(I-1)                                                         
      DO 40 J=1,N                                                        
      JK=NK+J                                                            
      JI=JP+J                                                            
      HOLD=-A(JK)                                                        
      A(JK)=A(JI)                                                        
   40 A(JI) =HOLD                                                        
   45 IF(BIGA.NE.0.D0) GOTO 48                                           
      D=0.D0                                                             
      RETURN                                                             
   48 DO 55 I=1,N                                                        
      IF(I.EQ.K) GOTO 55                                                 
      IK=NK+I                                                            
      A(IK)=A(IK)/(-BIGA)                                                
   55 CONTINUE                                                           
      DO 65 I=1,N                                                        
      IK=NK+I                                                            
      HOLD=A(IK)                                                         
      IJ=I-N                                                             
      DO 65 J=1,N                                                        
      IJ=IJ+N                                                            
      IF(I.EQ.K) GOTO 65                                                 
      IF(J.EQ.K) GOTO 65                                                 
      KJ=IJ-I+K                                                          
      A(IJ)=HOLD*A(KJ)+A(IJ)                                             
   65 CONTINUE                                                           
      KJ=K-N                                                             
      DO 75 J=1,N                                                        
      KJ=KJ+N                                                            
      IF(J.EQ.K) GOTO 75                                                 
      A(KJ)=A(KJ)/BIGA                                                   
   75 CONTINUE                                                           
      D=D*BIGA                                                           
      A(KK)=1.D0/BIGA                                                    
   80 CONTINUE                                                           
      K=N                                                                
  100 K=(K-1)                                                            
      IF(K.LE.0) RETURN                                                  
      I=L(K)                                                             
      IF(I.LE.K) GOTO 120                                                
      JQ=N*(K-1)                                                         
      JR=N*(I-1)                                                         
      DO 110 J=1,N                                                       
      JK=JQ+J                                                            
      HOLD=A(JK)                                                         
      JI=JR+J                                                            
      A(JK)=-A(JI)                                                       
  110 A(JI) =HOLD                                                        
  120 J=M(K)                                                             
      IF(J.LE.K) GOTO 100                                                
      KI=K-N                                                             
      DO 130 I=1,N                                                       
      KI=KI+N                                                            
      HOLD=A(KI)                                                         
      JI=KI-K+J                                                          
      A(KI)=-A(JI)                                                       
  130 A(JI) =HOLD                                                        
      GO TO 100                                                          
      END                                                                
      SUBROUTINE NEKMIN(RM,OMEG,X,Z)                                     
c  Version of October 9, 1995                                            
      IMPLICIT REAL*8(A-H,O-Z)                                           
      DIMENSION DN(4),EN(2),OUT(2),LL(2),MM(2)                           
      Z=.05d0                                                            
   15 P1=X*X+Z*Z                                                         
      RP1=DSQRT(P1)                                                      
      P115=P1*RP1                                                        
      P2=(1.d0-X)**2+Z*Z                                                 
      RP2=DSQRT(P2)                                                      
      P215=P2*RP2                                                        
      DODZ=-Z/P115-RM*Z/P215                                             
      OM=1.d0/RP1+RM/RP2+(1.d0+RM)*.5d0*X*X-RM*X                         
      DELOM=OMEG-OM                                                      
      DELZ=DELOM/DODZ                                                    
      Z=DABS(Z+DELZ)                                                     
      ABDELZ=DABS(DELZ)                                                  
      IF(ABDELZ.GT..00001d0) GOTO 15                                     
   16 P1=X*X+Z*Z                                                         
      RP1=DSQRT(P1)                                                      
      P115=P1*RP1                                                        
      P125=P1*P115                                                       
      P2=(1.d0-X)**2+Z*Z                                                 
      RP2=DSQRT(P2)                                                      
      P215=P2*RP2                                                        
      P225=P2*P215                                                       
      DN(1)=-X/P115+RM*(1.d0-X)/P215+(1.d0+RM)*X-RM                      
      DN(2)=(3.d0*X*X-P1)/P125+(3.d0*RM*(1.d0-X)**2-RM*((1.d0-X)**2      
     $+z*z))/p225+(RM+1.d0)                                              
      DN(3)=-Z/P115-RM*Z/P215                                            
      DN(4)=3.d0*X*Z/P125-3.d0*RM*Z*(1.d0-X)/P225                        
      OME=1.d0/RP1+RM/RP2+(1.d0+RM)*.5d0*X*X-RM*X                        
      EN(1)=OMEG-OME                                                     
      EN(2)=-DN(1)                                                       
      CALL DMINV(DN,2,D,LL,MM)                                           
      CALL DGMPRD(DN,EN,OUT,2,2,1)                                       
      DT1=OUT(1)                                                         
      DT2=OUT(2)                                                         
      ABDX=DABS(DT1)                                                     
      X=X+DT1                                                            
      ABDZ=DABS(DT2)                                                     
      Z=Z+DT2                                                            
      IF(ABDX.GT.1.d-8) GOTO 16                                          
      IF(ABDZ.GT.1.d-8) GOTO 16                                          
      RETURN                                                             
      END                                                                
      SUBROUTINE ROMQ(omein,Q,F,D,EC,TH,FI,R,DRDO,DRDQ,DODQ,KOMP,MODE)     
c  Version of December 5, 2003                                         
      implicit real*8 (a-h,o-z)                                          
      theq=1.570796326794897d0                                           
      MOD46=(MODE-5)**2                                                  
      MOD56=(2*MODE-11)**2                                               
      modkom=mode*(komp+komp-3)                                          
      ome=omein 
      DQ=1.d-4*Q                                                         
      QP=Q+DQ                                                            
      TOL=5.d-8                                                          
C     TH, FI SHOULD BE IN RADIANS.                                       
      sinth=dsin(th)                                                     
      XNUSQ=sinth*sinth                                                  
      XLAM=sinth*dcos(FI)                                                
      RMA=Q                                                              
      QF=1.d0                                                            
      DP=1.d0-EC                                                         
      QFM=1.d0                                                           
      IF(KOMP.NE.2) GOTO 23                                              
      RMA=1.d0/Q                                                         
      QF=1.d0/Q                                                          
      QFM=-1.d0/Q**2                                                     
   23 CONTINUE                                                           
      CALL ELLONE(F,DP,RMA,X,OMEG,XLD,OMD)                               
      OM2SAV=OMEG                                                        
      RMAP=QP                                                            
      IF(KOMP.NE.2) GOTO 92                                              
      OMEG=OMEG*Q+(1.d0-Q)*.5d0                                          
      IF(MOD56.EQ.1) OME=OMEG                                            
      RMAP=1.d0/QP                                                       
      GOTO 93                                                            
   92 CONTINUE                                                           
      IF(MOD46.EQ.1) OME=OMEG                                            
   93 CONTINUE                                                           
      POT=OME                                                            
      IF(KOMP.EQ.2) POT=OME/Q+.5d0*(Q-1.d0)/Q                            
      CALL ELLONE(F,DP,RMAP,XP,OMP,XLD,OMD)                              
      DODQ=(OMP-OM2SAV)/DQ                                               
      RM1=RMA+1.d0                                                       
      DS=D*D                                                             
      RF=F*F                                                             
      R=1.d0/POT                                                         
      KOUNT=0                                                            
      DELR=0.d0                                                          
      IF(FI.NE.0.d0) GOTO 85                                             
      IF(TH.NE.THEQ) GOTO 85                                             
      IF(MODE.EQ.6) GOTO 114                                             
      IF(MODE.NE.4) GOTO 80                                              
      IF(KOMP.EQ.1) GOTO 114                                             
      GOTO 85                                                            
   80 IF(MODE.NE.5) GOTO 85                                              
      IF(KOMP.EQ.2) GOTO 114                                             
   85 CONTINUE                                                           
   14 R=R+DELR                                                           
      KOUNT=KOUNT+1                                                      
      IF(KOUNT.LT.20) GOTO 70                                            
  217 if(mode.eq.6) goto 114                                             
      if(modkom.eq.-4) goto 114                                          
      if(modkom.eq.5) goto 114                                           
      DOMR=-1.d15                                                        
      R=-1.d0                                                            
      GOTO 116                                                           
   70 RSQ=R*R                                                            
      PAR=DS-2.d0*XLAM*R*D+RSQ                                           
      RPAR=dsqrt(PAR)                                                    
      OM=1.d0/R+RMA*(1.d0/RPAR-XLAM*R/DS)+RM1*.5d0*RSQ*XNUSQ*RF          
      DOMR=1.d0/(RF*RM1*XNUSQ*R-1.d0/RSQ-(RMA*(R-XLAM*D))/(PAR*RPAR)-    
     $RMA*XLAM/DS)                                                       
      DELR=(POT-OM)*DOMR                                                 
      ABDELR=dabs(DELR)                                                  
      IF(ABDELR.GT.TOL) GOTO 14                                          
      DOMRSV=DOMR                                                        
      IF(R.GE.1.d0) GOTO 217                                             
      IF(FI.NE.0.d0) GO TO 116                                           
      IF(TH.NE.THEQ)GO TO 116                                            
      IF(OME-OMEG) 217,114,116                                           
  114 DOMR=1.d15                                                         
      R=X                                                                
      goto 118                                                           
  116 DRDQ=(1.d0/RPAR-R*XLAM/DS+.5d0*RF*RSQ*XNUSQ)/(1.d0/RSQ+RMA*        
     $((1.d0/(PAR*RPAR))*(R-XLAM*D)+XLAM/DS)-RF*XNUSQ*RM1*R)             
      DRDQ=DRDQ*QFM                                                      
  118 drdo=domr*qf                                                       
      IF(MODE.EQ.6) GOTO 215                                             
      IF(MODE.NE.4) GOTO 180                                             
      IF(KOMP.EQ.1) GOTO 215                                             
      RETURN                                                             
  180 IF(MODE.NE.5) RETURN                                               
      IF(KOMP.EQ.2) GOTO 215                                             
      RETURN                                                             
  215 IF(FI.NE.0.d0) GOTO 230                                            
      IF(TH.NE.THEQ) GOTO 230                                            
      DRDQ=(XP-X)/DQ                                                     
      RETURN                                                             
  230 DRDQ=DRDQ+DOMRSV*DODQ                                              
      RETURN                                                             
      END                                                                
      subroutine cloud(cosa,cosb,cosg,x1,y1,z1,xc,yc,zc,rr,wl,op1,       
     $opsf,edens,acm,en,cmpd,ri,dx,dens,tau)                             
c  Version of January 9, 2002                                         
      implicit real*8 (a-h,o-z)                                          
      dx=0.d0                                                            
      tau=0.d0                                                           
      ri=1.d0                                                            
      sige=.6653d-24                                                     
      dtdxes=sige*edens                                                  
c  cosa can be zero, so an alternate path to the solution is needed      
      dabcoa=dabs(cosa)                                                  
      dabcob=dabs(cosb)                                                  
      if(dabcoa.lt.dabcob) goto 32                                       
      w=cosb/cosa                                                        
      v=cosg/cosa                                                        
      u=y1-yc-w*x1                                                       
      t=z1-zc-v*x1                                                       
      aa=1.d0+w*w+v*v                                                    
      bb=2.d0*(w*u+v*t-xc)                                               
      cc=xc*xc+u*u+t*t-rr*rr                                             
      dubaa=aa+aa                                                        
      dis=bb*bb-4.d0*aa*cc                                               
      if(dis.le.0.d0) return                                             
      sqd=dsqrt(dis)                                                     
      xx1=(-bb+sqd)/dubaa                                                
      xx2=(-bb-sqd)/dubaa                                                
      yy1=w*(xx1-x1)+y1                                                  
      yy2=w*(xx2-x1)+y1                                                  
      zz1=v*(xx1-x1)+z1                                                  
      zz2=v*(xx2-x1)+z1                                                  
      goto 39                                                            
   32 w=cosa/cosb                                                        
      v=cosg/cosb                                                        
      u=x1-xc-w*y1                                                       
      t=z1-zc-v*y1                                                       
      aa=1.d0+w*w+v*v                                                    
      bb=2.d0*(w*u+v*t-yc)                                               
      cc=yc*yc+u*u+t*t-rr*rr                                             
      dubaa=aa+aa                                                        
      dis=bb*bb-4.d0*aa*cc                                               
      if(dis.le.0.d0) return                                             
      sqd=dsqrt(dis)                                                     
      yy1=(-bb+sqd)/dubaa                                                
      yy2=(-bb-sqd)/dubaa                                                
      xx1=w*(yy1-y1)+x1                                                  
      xx2=w*(yy2-y1)+x1                                                  
      zz1=v*(yy1-y1)+z1                                                  
      zz2=v*(yy2-y1)+z1                                                  
   39 dis=bb*bb-4.d0*aa*cc                                               
      if(dis.le.0.d0) return                                             
      sqd=dsqrt(dis)                                                     
      xs1=(xx1-cmpd)*cosa+yy1*cosb+zz1*cosg                              
      xs2=(xx2-cmpd)*cosa+yy2*cosb+zz2*cosg                              
      xxnear=xx1                                                         
      yynear=yy1                                                         
      zznear=zz1                                                         
      xxfar=xx2                                                          
      yyfar=yy2                                                          
      zzfar=zz2                                                          
      xsnear=xs1                                                         
      xsfar=xs2                                                          
      if(xs1.gt.xs2) goto 38                                             
      xxnear=xx2                                                         
      yynear=yy2                                                         
      zznear=zz2                                                         
      xxfar=xx1                                                          
      yyfar=yy1                                                          
      zzfar=zz1                                                          
      xsnear=xs2                                                         
      xsfar=xs1                                                          
   38 continue                                                           
      xss=(x1-cmpd)*cosa+y1*cosb+z1*cosg                                 
      if(xss.ge.xsnear) return                                           
      if(xss.le.xsfar) goto 20                                           
      xxfar=x1                                                           
      yyfar=y1                                                           
      zzfar=z1                                                           
   20 continue                                                           
      dtaudx=dtdxes+(op1*wl**en+opsf)*dens                               
      dx=dsqrt((xxnear-xxfar)**2+(yynear-yyfar)**2+(zznear-zzfar)**2)    
      tau=dx*dtaudx*acm                                                  
      ri=dexp(-tau)                                                      
      return                                                             
      end                                                                
      subroutine linpro(komp,dvks,hbarw,tau,emm,count,taug,emmg,fbin,       
     $delv) 
c  Version of November 3, 2000                                              
      implicit real*8(a-h,o-z)                                           
      dimension dvks(*),hbarw(*),tau(*),emm(*),count(*),fbin(*),delv(*), 
     $taug(*),emmg(*) 
      common /flpro/ vks,binc,binw,difp,dum1,dum2                        
      common /ipro/ nbins,nl,inmax,inmin,idum1,idum2                     
      COMMON /SPOTS/ SINLAT(2,100),COSLAT(2,100),SINLNG(2,100),COSLNG    
     $(2,100),RAD(2,100),TEMSP(2,100),xlng(2,100),kks(2,100),            
     $Lspot(2,100)                                                       
      inmin=300000                                                        
      inmax=0                                                            
c 
c  The 83 loop pre-computes the limiting bin numbers, encompassing all lines 
c 
      do 83 iln=1,nl                                                     
      if(Lspot(komp,iln).eq.0) goto 83                                   
      vksg=vks+dvks(iln)                                                 
      vksp=vksg+hbarw(iln)                                               
      vksm=vksg-hbarw(iln)                                               
      indp=vksp/binw+binc                                                
      indm=vksm/binw+binc                                                
      if(indm.lt.inmin) inmin=indm                                       
      if(indp.gt.inmax) inmax=indp                                       
   83 continue                                                           
      do 82 i=inmin,inmax                                                
      emmg(i)=0.d0 
   82 taug(i)=0.d0                                                       
c 
c  The 84 loop puts fractional contributions into the two end bins 
c    (first part, up to 28 continue) and puts full contributions 
c    into the middle bins (the 26 loop). 
c 
      do 84 iln=1,nl                                                     
      if(Lspot(komp,iln).eq.0) goto 84                                   
      vksg=vks+dvks(iln)                                                 
      vksp=vksg+hbarw(iln)                                               
      vksm=vksg-hbarw(iln)                                               
      indp=vksp/binw+binc                                                
      indm=vksm/binw+binc                                                
      vks1=(dfloat(indm+1)-binc)*binw                                    
      vks2=(dfloat(indp)-binc)*binw                                      
      fr1=(vks1-vksm)/binw                                               
      fr2=(vksp-vks2)/binw                                               
      taug(indm)=taug(indm)+fr1*tau(iln)                                 
      emmg(indm)=emmg(indm)+fr1*emm(iln) 
      delv(indm)=delv(indm)+fr1*vksm                                     
      count(indm)=count(indm)+fr1                                        
      taug(indp)=taug(indp)+fr2*tau(iln)                                 
      emmg(indp)=emmg(indp)+fr2*emm(iln) 
      delv(indp)=delv(indp)+fr2*vksp                                     
      count(indp)=count(indp)+fr2                                        
      if(indp.ne.indm) goto 28                                           
      taug(indp)=taug(indp)-tau(iln)                                     
      emmg(indp)=emmg(indp)-emm(iln) 
      delv(indp)=delv(indp)-.5d0*(vksm+vksp)                             
      count(indp)=count(indp)-1.d0                                       
   28 continue                                                           
      ind=indm                                                           
      idmax=indp-indm-1                                                  
      if(idmax.le.0) goto 84                                             
      do 26 id=1,idmax                                                   
      ind=ind+1                                                          
      vksb=(dfloat(ind)-binc)*binw                                       
      taug(ind)=taug(ind)+tau(iln)                                       
      emmg(ind)=emmg(ind)+emm(iln) 
      delv(ind)=delv(ind)+vksb                                           
      count(ind)=count(ind)+1.d0                                         
   26 continue                                                           
   84 continue                                                           
c 
c  The 85 loop collects the absorption and emission contributions to all 
c     active bins, with the absorption lines summed via an optical thickness 
c     treatment and the emission lines summed directly according to contributed 
c     flux. The sign on emmg is negative because emmg is negative. 
c 
      do 85 i=inmin,inmax                                                
   85 fbin(i)=fbin(i)+(1.d0-dexp(-taug(i))+emmg(i))*difp                         
      return                                                             
      end                                                                
      subroutine jdph(xjdin,phin,t0,p0,dpdt,xjdout,phout)                
c  Version of February 2, 1999                                           
c                                                                        
c  Subroutine jdph computes a phase (phout) based on an input            
c   JD (xjdin), reference epoch (t0), period (p0), and dP/dt (dpdt).     
c   It also computes a JD (xjdout) from an input phase (phin) and the    
c   same ephemeris. So jdph can be used either to get phase from         
c   JD or JD from phase.                                                 
c                                                                        
      implicit real*8(a-h,o-z)                                           
      tol=1.d-6                                                          
      abdpdt=dabs(dpdt)                                                  
      deltop=(xjdin-t0)/p0                                               
      fcsq=0.d0                                                          
      fccb=0.d0                                                          
      fc4th=0.d0                                                         
      fc5th=0.d0                                                         
      fc=deltop*dpdt                                                     
      if(dabs(fc).lt.1.d-18) goto 25                                     
      fcsq=fc*fc                                                         
      if(dabs(fcsq).lt.1.d-24) goto 25                                   
      fccb=fc*fcsq                                                       
      if(dabs(fccb).lt.1.d-27) goto 25                                   
      fc4th=fc*fccb                                                      
      if(dabs(fc4th).lt.1.d-28) goto 25                                  
      fc5th=fc*fc4th                                                     
   25 phout=deltop*(1.d0-.5d0*fc+fcsq/3.d0-.25d0*fccb+.2d0*fc4th-        
     $fc5th/6.d0)                                                        
      pddph=dpdt*phin                                                    
      xjdout=p0*phin*(1.d0+.5d0*pddph+pddph**2/6.d0+pddph**3/24.d0       
     $+pddph**4/120.d0+pddph**5/720.d0)+t0                               
      if(abdpdt.lt.tol) return                                           
      phout=dlog(1.d0+deltop*dpdt)/dpdt                                  
      xjdout=(dexp(dpdt*phin)-1.d0)*p0/dpdt+t0                           
      return                                                             
      end                                                                
      subroutine legendre(x,pleg,n) 
c  Version of January 7, 2002 
      implicit real*8 (a-h,o-z) 
      dimension pleg(n) 
      pleg(1)=1.d0 
      pleg(2)=x 
      if(n.le.2) return 
      denom=1.d0 
      do 1 i=3,n 
      fac1=x*(2.d0*denom+1.d0) 
      fac2=denom 
      denom=denom+1.d0 
      pleg(i)=(fac1*pleg(i-1)-fac2*pleg(i-2))/denom 
   1  continue 
      return 
      end 
      subroutine binnum(x,n,y,j) 
c  Version of January 7, 2002 
      implicit real*8(a-h,o-z) 
      dimension x(n) 
      mon=1 
      if(x(1).gt.x(2)) mon=-1 
      do 1 i=1,n 
      if(mon.eq.-1) goto 3 
      if(y.le.x(i)) goto 2 
      goto 1 
   3  if(y.gt.x(i)) goto 2 
   1  continue 
   2  continue 
      j=i-1 
      return 
      end 
      SUBROUTINE square (OBS,NOBS,ML,OUT,sd,xlamda,D,CN,CNN,cnc,          
     $clc,ss,cl,ll,mm)                                                   
c  Version of January 16, 2002                                               
      implicit real*8 (a-h,o-z)                                          
      dimension obs(*),out(*),sd(*),cn(*),cnn(*),cnc(*),clc(*),          
     $cl(*),ll(*),mm(*)                                                  
c                                                                        
c  cnc ("cn copy") is the original normal equation matrix                
c  cn is the re-scaled version of cnc                                    
c  cnn comes in as the original n.e. matrix, then is copied              
c    from cn to become the re-scaled n.e.'s, and finally is              
c    inverted by DMINV to become the inverse of the re-scaled            
c    n.e. matrix.                                                        
c                                                                        
      S=0.D0                                                             
      CLL=0.D0                                                           
      CAY=NOBS-ML                                                        
      JMAX=ML*ML                                                         
      DO 20 J=1,JMAX                                                     
   20 CN(J)=0.D0                                                         
      DO 21 J=1,ML                                                       
   21 CL(J)=0.D0                                                         
      DO 25 NOB=1,NOBS                                                   
      III=NOB+NOBS*ML                                                    
      OBSQQ=OBS(III)                                                     
      DO 23 K=1,ML                                                       
      DO 23 I=1,ML                                                       
      II=NOB+NOBS*(I-1)                                                  
      KK=NOB+NOBS*(K-1)                                                  
      J=I+(K-1)*ML                                                       
      OBSII=OBS(II)                                                      
      OBSKK=OBS(KK)                                                      
      CN(J)=CN(J)+OBSII*OBSKK                                            
   23 cnc(j)=cn(j)                                                       
      DO 24 I=1,ML                                                       
      II=NOB+NOBS*(I-1)                                                  
      OBSII=OBS(II)                                                      
   24 CL(I)=CL(I)+OBSQQ*OBSII                                            
   25 CLL=CLL+OBSQQ*OBSQQ                                                
      do 123 k=1,ml                                                      
      do 123 i=1,ml                                                      
      xlf=0.d0                                                           
      if(i.eq.k) xlf=xlamda                                              
      j=i+(k-1)*ml                                                       
      ji=i+(i-1)*ml                                                      
      ki=k+(k-1)*ml                                                      
  123 cn(j)=cn(j)/dsqrt(cnc(ji)*cnc(ki))+xlf                             
      do 124 i=1,ml                                                      
      ji=i+(i-1)*ml                                                      
      clc(i)=cl(i)                                                       
  124 cl(i)=cl(i)/dsqrt(cnc(ji))                                         
      DO 50 J=1,JMAX                                                     
   50 CNN(J)=CN(J)                                                       
      CALL DMINV(CNN,ML,D,LL,MM)                                         
      CALL DGMPRD(CNN,CL,OUT,ML,ML,1)                                    
      do 125 i=1,ml                                                      
      ji=i+(i-1)*ml                                                      
  125 out(i)=out(i)/dsqrt(cnc(ji))                                       
      DO 26 I=1,ML                                                       
   26 S=S+clc(I)*OUT(I)                                                  
      S=CLL-S                                                            
      SS=S                                                               
      SIGSQ=S/CAY                                                        
      DO 27 J=1,ML                                                       
      JJ=J*ML+J-ML                                                       
      CNJJ=CNN(JJ)                                                       
      ARG=SIGSQ*CNJJ                                                     
   27 sd(J)=dsqrt(arg/cnc(jj))                                           
      RETURN                                                             
      END                                                                
      subroutine conjph(ecc,argper,phzero,trsc,tric,econsc,econic,
     $xmsc,xmic,pconsc,pconic)
      implicit real*8(a-h,o-z)
c  Version of December 15, 2003
c
c  Subroutine conjph computes the phases of superior and inferior conjunction
c    (pconsc and pconic) of star 1
c
      pi=dacos(-1.d0)
      pih=.5d0*pi
      pi32=1.5d0*pi
      twopi=pi+pi
      ecfac=dsqrt((1.d0-ecc)/(1.d0+ecc))
c
c  sc in variable names (like trsc) means superior conjunction, and
c  ic means inferior conjunction (always for star 1).
c
      trsc=pih-argper
      tric=pi32-argper
      econsc=2.d0*datan(ecfac*dtan(.5d0*trsc))
      econic=2.d0*datan(ecfac*dtan(.5d0*tric))
      xmsc=econsc-ecc*dsin(econsc)
      xmic=econic-ecc*dsin(econic)
      pconsc=(xmsc+argper)/twopi-.25d0+phzero
      pconic=(xmic+argper)/twopi-.25d0+phzero
      return
      end
      SUBROUTINE KEPLER(XM,EC,ECAN,TR)                                   
c  Version of July 22, 2004                                            
      IMPLICIT REAL*8(A-H,O-Z)                                           
      TOL=1.d-10                                                          
      DLECAN=0.D0                                                        
      ECAN=XM                                                            
   18 ECAN=ECAN+DLECAN                                                   
      XMC=ECAN-EC*DSIN(ECAN)                                             
      DEDM=1.D0/(1.D0-EC*DCOS(ECAN))                                     
      DLECAN=(XM-XMC)*DEDM                                               
      ABDLEC=DABS(DLECAN)                                                
      IF(ABDLEC.GT.TOL) GOTO 18                                          
      TR=2.D0*DATAN(DSQRT((1.D0+EC)/(1.D0-EC))*DTAN(.5D0*ECAN))          
      IF(TR.LT.0.) TR=TR+6.2831853071795865d0                            
      RETURN                                                             
      END                                                                
      SUBROUTINE BBL(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVE,FR1,FR2,   
     $HLD,SLUMP1,SLUMP2,THETA,RHO,AA,BB,PHSV,PCSV,N1,N2,F1,F2,d,hlum,    
     $clum,xh,xc,yh,yc,gr1,gr2,wl,sm1,sm2,tpolh,tpolc,sbrh,sbrc,   
     $tavh,tavc,alb1,alb2,xbol1,xbol2,ybol1,ybol2,fspot1,fspot2,phas,rm,        
     $xincl,hot,cool,snth,csth,snfi,csfi,tld,glump1,glump2,xx1,xx2,      
     $yy1,yy2,zz1,zz2,dint1,dint2,grv1,grv2,rftemp,rf1,rf2,csbt1,        
     $csbt2,gmag1,gmag2,glog1,glog2,obser,fbin1,fbin2,delv1,delv2, 
     $count1,count2,delwl1,delwl2,resf1,resf2,wl1,wl2,dvks1,dvks2,tau1,
     $tau2,emm1,emm2,hbarw1,hbarw2,xcl,ycl,zcl,rcl,op1,fcl,dens,encl, 
     $edens,taug,emmg,yskp,zskp,mode,iband,ifat1,ifat2,ifphn)      
c  Version of February 28, 2012                                         
      implicit real*8 (a-h,o-z)                                          
      parameter (ispmax=100)
      DIMENSION RV(*),GRX(*),GRY(*),GRZ(*),RVQ(*),GRXQ(*),GRYQ(*),       
     $GRZQ(*),MMSAVE(*),FR1(*),FR2(*),HLD(*),SLUMP1(*),SLUMP2(*),        
     $THETA(*),RHO(*),AA(*),BB(*),SNTH(*),CSTH(*),SNFI(*),CSFI(*),TLD(*) 
     $,GLUMP1(*),GLUMP2(*),XX1(*),XX2(*),YY1(*),YY2(*),ZZ1(*),ZZ2(*),    
     $GRV1(*),GRV2(*),RFTEMP(*),RF1(*),RF2(*),CSBT1(*),CSBT2(*)          
     $,GMAG1(*),GMAG2(*),glog1(*),glog2(*),obser(*)                                                 
      dimension fbin1(*),fbin2(*),delv1(*),delv2(*),count1(*),count2(*), 
     $delwl1(*),delwl2(*),resf1(*),resf2(*),wl1(*),wl2(*),dvks1(*),      
     $dvks2(*),tau1(*),tau2(*),hbarw1(*),hbarw2(*),taug(*),emm1(*), 
     $emm2(*),emmg(*)               
      dimension xcl(*),ycl(*),zcl(*),rcl(*),op1(*),fcl(*),dens(*),       
     $edens(*),encl(*),yskp(*),zskp(*)                                   
      common /klsp/ tratio,lsp
      COMMON /INVAR/ KH,ipbdum,IRTE,NREF,IRVOL1,IRVOL2,mref,ifsmv1,
     $ifsmv2,icor1,icor2,ld1,ld2,ncl,jdphs,ipc,nr
      COMMON /FLVAR/ PSHIFT,DP,EF,EFC,ECOS,perr0,PHPER,pconsc,pconic,         
     $PHPERI,VSUM1,VSUM2,VRA1,VRA2,VKM1,VKM2,VUNIT,vfvu,trc,qfacd        
      common /nspt/ nsp1,nsp2,kspot,kspev                             
      common /ardot/ dperdt,hjd,hjd0,perr                                
      COMMON /spots/ snlat(2,ispmax),cslat(2,ispmax),snlng(2,ispmax),
     $cslng(2,ispmax),rdsp(2,ispmax),tmsp(2,ispmax),xlng(2,ispmax),
     $kks(2,ispmax),Lspot(2,ispmax)
      COMMON /ECCEN/ E,A,PERIOD,VGA,SINI,VF,VFAC,VGAM,VOL1,VOL2,IFC      
      common /prof2/ vo1,vo2,ff1,ff2,du1,du2,du3,du4,du5,du6,du7         
      common /con3b/ xm03b,t03b,e3b,p3b,perr3b,vfac3b,ecos3b,ifadj3b
      common /spev/ tstart(2,ispmax),tmaxa(2,ispmax),tmaxb(2,ispmax),
     $tfinal(2,ispmax),amax(2,ispmax),sprad(2,ispmax),norscalc
      pi=dacos(-1.d0)                                             
      twopi=pi+pi                                                        
      ff1=f1                                                             
      ff2=f2                                                             
      qfac1=1.d0/(1.d0+rm)                                               
      qfac=rm*qfac1                                                      
      if(mode.ne.1.or.ld1.lt.0.or.ld2.lt.0) goto 84
      xc=xh
      yc=yh
   84 continue
      PSFT=PHAS-PHPERI                                                   
   29 if(PSFT.GT.1.d0) PSFT=PSFT-1.d0                                    
      if(psft.gt.1.d0) goto 29                                           
   30 if(PSFT.LT.0.d0) PSFT=PSFT+1.d0                                    
      if(psft.lt.0.d0) goto 30                                           
      XMEAN=PSFT*twopi                                                   
      tr=xmean                                                           
      if(norscalc.eq.1) goto 63
      do 60 kp=1,2                                                       
      nsp=nsp1*(2-kp)+nsp2*(kp-1)      
      do 62 iev=1,nsp                                                      
c
c   The next line sets the spot radius default for the case of no spot evolution
c
      sprad(kp,iev)=rdsp(kp,iev)
      if(kspev.eq.0) goto 62
      call spotev(hjd,tstart(kp,iev),tmaxa(kp,iev),tmaxb(kp,iev),
     $tfinal(kp,iev),amax(kp,iev),area,sprad(kp,iev))
   62 continue
      ifsmv=ifsmv1*(2-kp)+ifsmv2*(kp-1)                                  
      if(ifsmv.eq.0) goto 60                                             
      do 61 i=1,nsp                                                      
      ff=fspot1*dfloat(2-kp)+fspot2*dfloat(kp-1)                                 
      xlg=xlng(kp,i)+twopi*ff*(phas-pconsc)-(tr-trc)                      
      snlng(kp,i)=dsin(xlg)                                              
      cslng(kp,i)=dcos(xlg)                                             
   61 continue                                                           
   60 continue                                                           
   63 continue
      if(e.ne.0.d0) call KEPLER(XMEAN,E,DUM,TR)                          
      U=TR+PERR                                                          
      COSU=dcos(U)                                                       
      GPHA=U*.1591549d0-.25d0                                            
   40 if(GPHA.lt.0.d0) GPHA=GPHA+1.d0                                    
      if(gpha.lt.0.d0) goto 40                                           
   50 if(GPHA.GE.1.d0) GPHA=GPHA-1.d0                                    
      if(gpha.ge.1.d0) goto 50                                           
      D=EF/(1.d0+E*dcos(TR))                                             
      qfacd=qfac*d                                                       
      IF(IRTE.EQ.1) GOTO 19                                              
      CALL LCR(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVE,FR1,FR2,HLD, 
     $slump1,SLUMP2,RM,PHSV,PCSV,N1,N2,F1,F2,D,HLUM,CLUM,xh,xc,yh,yc, 
     $gr1,gr2,SM1,SM2,TPOLH,TPOLC,SBRH,SBRC,IFAT1,IFAT2,TAVH,TAVC, 
     $alb1,alb2,xbol1,xbol2,ybol1,ybol2,vol1,vol2,snth,csth,snfi,csfi, 
     $tld,glump1,glump2,xx1,xx2,yy1,yy2,zz1,zz2,dint1,dint2,grv1,grv2, 
     $csbt1,csbt2,rftemp,rf1,rf2,gmag1,gmag2,glog1,glog2,mode,iband) 
   19 CONTINUE                                                           
      VO1=qfac*SINI*(ECOS+COSU)/EFC+VGAM                                 
      VO2=-qfac1*SINI*(ECOS+COSU)/EFC+VGAM                               
c
c  Here's the special call to compute TPOLE/TOBS
c
      if(lsp.eq.0) goto 111
      call lightsp(gpha,xincl,xh,xc,yh,yc,n1,n2,hot,cool,rv,grx,gry,grz,   
     $rvq,grxq,gryq,grzq,mmsave,theta,rho,aa,bb,somhot,
     $somkul,d,snth,csth,snfi,csfi,tld,gmag1,gmag2,obser,
     $glump1,glump2)
  111 continue
      call light(gpha,xincl,xh,xc,yh,yc,n1,n2,hot,cool,rv,grx,gry,grz,   
     $rvq,grxq,gryq,grzq,mmsave,theta,rho,aa,bb,slump1,slump2,somhot,    
     $somkul,d,wl,snth,csth,snfi,csfi,tld,gmag1,gmag2,glog1,glog2,obser,
     $fbin1,fbin2,delv1,delv2,count1,count2,delwl1,delwl2,resf1,resf2,
     $wl1,wl2,dvks1,dvks2,tau1,tau2,emm1,emm2,hbarw1,hbarw2,xcl,ycl,zcl,
     $rcl,op1,fcl,edens,encl,dens,taug,emmg,yskp,zskp,iband,ifat1,ifat2,
     $ifphn)                                    
      VRA1=0.d0                                                          
      VRA2=0.d0                                                          
      IF(HOT.GT.0.d0) VRA1=F1*SOMHOT/HOT                                 
      IF(COOL.GT.0.d0) VRA2=F2*SOMKUL/COOL                               
      vsum1=vo1                                                          
      vsum2=vo2                                                          
      if(icor1.eq.1) vsum1=vo1+vra1                                      
      if(icor2.eq.1) vsum2=vo2+vra2                                    
c
c  Compute velocity due to 3rd body if necessary ******************
c
      v3b=0.d0
      if(ifadj3b.eq.0) goto 104
      xmean3b=xm03b+twopi*(hjd-t03b)/p3b
      call kepler(xmean3b,e3b,ecan3b,tr3b)
      v3b=vfac3b*(ecos3b+dcos(tr3b+perr3b))
  104 continue
      vfcc=vfac/vunit 
      VKM1=VSUM1*vfcc+v3b                                                    
      VKM2=VSUM2*vfcc+v3b                                                    
      RETURN                                                             
      END                                                                
      Subroutine lum (xlum,xlimb,ylimb,tpoll,n,n1,komp,sbr,rv,rvq,
     $glump1,glump2,glog1,glog2,grv1,grv2,mmsave,summ,fr,sm,ifat,vol,
     $rm,om,f,d,snth,iband)
c
c   Version of October 18, 2017                                       
c
      implicit real*8 (a-h,o-z)                                          
      dimension rv(*),rvq(*),mmsave(*),fr(*),snth(*),glump1(*),          
     $glump2(*),glog1(*),glog2(*),grv1(*),grv2(*)                                          
      dimension message(2,4) 
      common /atmmessages/ message,kompcom 
      common /radi/ R1H,RLH,R1C,RLC                                      
      common /invar/ khdum,ipbdum,irtedm,nrefdm,irv1dm,irv2dm,mrefdm,    
     $is1dm,is2dm,ic1dm,ic2dm,ld1,ld2,ncl,jdphs,ipc,nr                           
      common /gpoles/ gplog1,gplog2 
      common /poleint/ polin1,polin2
      common /coflimbdark/ xlocal,ylocal,xldmean1,xldmean2,yldmean1,
     $yldmean2
      kompcom=komp 
      ld=ld1
      if(komp.eq.2) ld=ld2
      ldabs=iabs(ld)
      pi=dacos(-1.d0)
      pih=.5d0*pi
      xlocal=xlimb
      ylocal=ylimb
      TPOLE=10000.d0*TPOLL                                               
      KR=0                                                               
      gplog=gplog1
      if(komp.eq.2) gplog=gplog2
      if(ifat.ne.0) goto 142
      if(ld.lt.0) call limdark(iband,ldabs,tpole,gplog,xlocal,ylocal)
      call planckint(tpole,iband,pollog,polin) 
  142 continue
      IF(IFAT.NE.0) call atmx(tpole,gplog,iband,pollog,polin)                            
      if(komp.eq.1) polin1=polin
      if(komp.eq.2) polin2=polin
      EN=dfloat(N)                                                       
      delth=pih/en
      SUM=0.d0                                                           
      SUMM=0.d0                                                          
      SM=0.d0                                                            
      sumxld=0.d0
      sumyld=0.d0
      VOL=0.d0                                                           
      DO 36 I=1,N                                                        
      IPN1=I+N1*(komp-1)                                                
      SINTH=SNTH(IPN1)                                                   
      EM=SINTH*EN*1.3d0                                                  
      MM=int(EM+1.d0)                                                         
      XM=dfloat(MM)                                                      
      delfi=pi/xm
      DFST=DELFI*SINTH                                                   
      SUMJ=0.d0                                                          
      SUMMJ=0.d0                                                         
      SMJ=0.d0                                                           
      sumjxld=0.d0
      sumjyld=0.d0
      VOLJ=0.d0                                                          
      DO 26 J=1,MM                                                       
      IP=(komp-1)*(N1+1)+I                                              
      IX=MMSAVE(IP)+J                                                    
      IF(komp.EQ.1) GOTO 39                                             
      IF(RVQ(IX).EQ.-1.d0) GOTO 25                                       
      R=RVQ(IX)                                                          
      grav=grv2(ix)
      glogg=glog2(ix)
      di=glump2(ix)
      GOTO 49                                                            
   39 IF(RV(IX).EQ.-1.d0) GOTO 25                                        
      R=RV(IX)                                                           
      grav=grv1(ix)
      glogg=glog1(ix)
      di=glump1(ix)
   49 continue
      TLOCAL=TPOLE*dsqrt(dsqrt(GRAV))                                    
      if(ifat.eq.0) call planckint(tlocal,iband,xinlog,xint) 
      IF(IFAT.NE.0) CALL atmx(tlocal,glogg,iband,xinlog,xint)                              
      GRAVM=xint/polin                         
      if(ld.lt.0) call limdark(iband,ldabs,tlocal,glogg,xlocal,ylocal)
      dif=di*gravm
      if(ld.ge.0) goto 112
      darkin=pi*(1.d0-xlocal/3.d0)
      if(ldabs.eq.2) darkin=darkin+.6981317d0*ylocal
      if(ldabs.eq.3) darkin=darkin-.6283185d0*ylocal
      dif=dif*darkin
  112 continue
      DIFF=DI*GRAV                                                       
      SMJ=SMJ+DI                                                         
      SUMJ=SUMJ+DIF                                                      
      SUMMJ=SUMMJ+DIFF                                                   
      sumjxld=sumjxld+xlocal*dif
      sumjyld=sumjyld+ylocal*dif
      VOLJ=VOLJ+R*R*R*FR(IX)                                             
      GOTO 26                                                            
   25 KR=1                                                               
   26 CONTINUE                                                           
      SMJ=SMJ*DELFI                                                      
      SUMJ=SUMJ*DELFI                                                    
      SUMMJ=SUMMJ*DELFI                                                  
      sumjxld=sumjxld*delfi
      sumjyld=sumjyld*delfi
      SM=SM+SMJ                                                          
      SUMM=SUMM+SUMMJ                                                    
      sumxld=sumxld+sumjxld
      sumyld=sumyld+sumjyld
      VOL=VOL+VOLJ*DFST                                                  
   36 SUM=SUM+SUMJ                                                       
      sbr=.25d0*xlum/(sum*delth)
      if(ld.lt.0) goto 111
      darkin=pi*(1.d0-xlocal/3.d0)
      if(ldabs.eq.2) darkin=darkin+.6981317d0*ylocal
      if(ldabs.eq.3) darkin=darkin-.6283185d0*ylocal
      sumxld=sumxld*darkin
      sumyld=sumyld*darkin
      sbr=sbr/darkin
  111 continue
      xldmean=sumxld/sum
      yldmean=sumyld/sum
      if(ld.ge.0) xldmean=xldmean/darkin
      if(ld.ge.0) yldmean=yldmean/darkin
      if(komp.eq.1) xldmean1=xldmean
      if(komp.eq.1) yldmean1=yldmean
      if(komp.eq.2) xldmean2=xldmean
      if(komp.eq.2) yldmean2=yldmean
      SM=SM*DELTH*4.d0                                                   
      SUMM=SUMM*DELTH*4.d0                                               
      VOL=VOL*1.3333333333333d0*DELTH                                    
      IF(KR.EQ.0) RETURN                                                 
      CALL ELLONE(F,D,RM,XL1,OMD,XLD,omdum)                                
      CALL NEKMIN(RM,OM,XL1,ZD)                                          
      IF(komp.EQ.2) XL1=D-XL1                                           
      R1=R1H
      if(komp.eq.2) R1=R1C
      RL=RLH
      if(komp.eq.2) RL=RLC
      VOL=VOL+1.047198d0*XL1*R1*RL                                       
      RETURN                                                             
      END                                                                
      Subroutine light(phs,xincl,xh,xc,yh,yc,n1,n2,sumhot,sumkul,rv,grx, 
     $gry,grz,rvq,grxq,gryq,grzq,mmsave,theta,rho,aa,bb,slump1,slump2,   
     $somhot,somkul,d,wl,snth,csth,snfi,csfi,tld,gmag1,gmag2,glog1,  
     $glog2,obser,fbin1,fbin2,delv1,delv2,count1,count2,delwl1,delwl2,   
     $resf1,resf2,wl1,wl2,dvks1,dvks2,tau1,tau2,emm1,emm2,hbarw1,hbarw2,
     $xcl,ycl,zcl,rcl,op1,fcl,edens,encl,dens,taug,emmg,yskp,zskp,iband,                             
     $ifat1,ifat2,ifphn)                              
c   Version of March 2, 2012                                        
      implicit real*8 (a-h,o-z)                                          
      parameter (ispmax=100)
      DIMENSION RV(*),GRX(*),GRY(*),GRZ(*),RVQ(*),GRXQ(*),GRYQ(*),GRZQ(* 
     $),SLUMP1(*),SLUMP2(*),MMSAVE(*),THETA(*),RHO(*),AA(*),BB(*)        
      DIMENSION SNTH(*),CSTH(*),SNFI(*),CSFI(*),tld(*),gmag1(*), 
     $gmag2(*),glog1(*),glog2(*),obser(*) 
      dimension xcl(*),ycl(*),zcl(*),rcl(*),op1(*),fcl(*),dens(*),       
     $encl(*),edens(*),yskp(*),zskp(*)                                   
      dimension fbin1(*),fbin2(*),delv1(*),delv2(*),count1(*),count2(*), 
     $delwl1(*),delwl2(*),resf1(*),resf2(*),wl1(*),wl2(*),dvks1(*),      
     $dvks2(*),tau1(*),tau2(*),hbarw1(*),hbarw2(*),taug(*),emm1(*), 
     $emm2(*),emmg(*)               
      dimension message(2,4) 
      common /kfac/ kff1,kff2,kfo1,kfo2
      common /atmmessages/ message,komp 
      common /coflimbdark/ xlimb,ylimb,dm1,dm2,dm3,dm4 
      COMMON /misc/ X1                                                          
      COMMON /NSPT/ NSP1,NSP2,kspot,kspev                               
      common /invar/ khdum,ipbdum,irtedm,nrefdm,irv1dm,irv2dm,mrefdm     
     $,ifs1dm,ifs2dm,icr1dm,icr2dm,ld1,ld2,ncl,jdphs,ipc,nr                      
      common /flvar/ du2,du3,du4,du5,du6,du7,du8,du9,du10,du11,          
     $du12,du13,du14,du15,du16,du17,vunit,vfvu,du20,qfacd                
      common /prof2/ vo1,vo2,ff1,ff2,binw1,binw2,sc1,sc2,sl1,sl2,        
     $clight                                                             
      common /cld/ acm,opsf                                              
      common /inprof/ in1min,in1max,in2min,in2max,mpage,nl1,nl2          
      common /setest/ sefac                                              
      common /flpro/ vksf,binc,binw,difp,deldel,renfsq                   
      common /ipro/ nbins,nl,inmax,inmin,nf1,nf2                         
      COMMON /spots/ snlat(2,ispmax),cslat(2,ispmax),snlng(2,ispmax),
     $cslng(2,ispmax),rdsp(2,ispmax),tmsp(2,ispmax),xlng(2,ispmax),
     $kks(2,ispmax),Lspot(2,ispmax)
      ot=1.d0/3.d0
      pi=dacos(-1.d0)                                                    
      twopi=pi+pi                                                        
      pih=.5d0*pi                                                        
      dtr=pi/180.d0                                                      
      kstp=4                                                             
      cirf=.002d0                                                        
      if(ifphn.eq.1) goto 16                                             
      if(mpage.ne.3) goto 16                                             
      nbins=90000                                                        
      binc1=.5d0*dfloat(nbins)                                           
      binc2=binc1                                                        
      in1max=0                                                           
      in2max=0                                                           
      in1min=300000                                                       
      in2min=300000                                                       
      marm1=10                                                           
      marp1=10                                                           
      marm2=10                                                           
      marp2=10                                                           
      do 916 i=1,nbins                                                   
      fbin1(i)=0.d0                                                      
      fbin2(i)=0.d0                                                      
      count1(i)=0.d0                                                     
      count2(i)=0.d0                                                     
      delv1(i)=0.d0                                                      
      delv2(i)=0.d0                                                      
  916 continue                                                           
   16 continue                                                           
      PHA=PHS*twopi                                                      
      K=6                                                                
      KK=K+1                                                             
      XINC=XINCL*dtr                                                     
      L=0                                                                
      TEST=(PHS-.5d0)**2                                                 
      TESTS=(TEST-.071525d0)**2                                          
      SINI=dsin(XINC)                                                    
      COSPH=dcos(PHA)                                                    
      SINPH=dsin(PHA)                                                    
      SINSQ=SINPH**2                                                     
      COSI=dcos(XINC)                                                    
      NP1=N1+1                                                           
      NP2=N1+N2+2                                                        
      LLL1=MMSAVE(NP1)                                                   
      LLL2=MMSAVE(NP2)                                                   
      NPP2=NP2-1                                                         
      LL1=MMSAVE(N1)+1                                                   
      LL2=MMSAVE(NPP2)+1                                                 
      LLLL1=(LL1+LLL1)/2                                                 
      LLLL2=(LL2+LLL2)/2                                                 
      SINSQE=0.d0                                                        
      IF(SINI.GT.0.d0) SINSQE=((1.02d0*(RV(LLL1)+RVQ(LLL2))/D)**2        
     $-cosi**2)/SINI**2                                                  
      if(sinsqe.gt..96d0) sinsqe=.96d0
      CICP=COSI*COSPH                                                    
      CISP=COSI*SINPH                                                    
      XLOS=COSPH*SINI                                                    
      YLOS=-SINPH*SINI                                                   
      ZLOS=COSI                                                          
      SUM=0.d0                                                           
      SOM=0.d0                                                           
      IF(TEST.LE..0625d0) GOTO 18                                        
      COMP=-1.d0                                                         
      CMP=1.d0                                                           
      COMPP=1.d0                                                         
      KOMP=2                                                             
      kff=kff2
      ld=ld2
      nl=nl2                                                             
      ffc=ff2                                                            
      voc=vo2*vfvu                                                       
      NSPOT=NSP2                                                         
      IFAT=IFAT2                                                         
      CMPP=0.d0                                                          
      xlimb=xc
      ylimb=yc
      EN=N2                                                              
      NPH=N2                                                             
      NP=2*N2                                                            
      nf=nf2                                                             
      GOTO 28                                                            
   18 xlimb=xh                                                               
      ylimb=yh                                                               
      COMP=1.d0                                                          
      KOMP=1                                                             
      kff=kff1
      ld=ld1
      nl=nl1                                                             
      ffc=ff1                                                            
      voc=vo1*vfvu                                                       
      NSPOT=NSP1                                                         
      IFAT=IFAT1                                                         
      CMP=0.d0                                                           
      COMPP=-1.d0                                                        
      CMPP=1.d0                                                          
      EN=N1                                                              
      NPH=N1                                                             
      NP=2*N1                                                            
      nf=nf1                                                             
   28 DELTH=pih/EN                                                       
      ldabs=iabs(ld)
      enf=dfloat(nf)                                                     
      renfsq=1.d0/(enf*enf)                                              
      nfm1=nf-1                                                          
      r2nfdt=0.5d0*delth/enf                                             
      vfvuff=vfvu*ffc                                                    
      AR=CMPP*RV(LLLL1)+CMP*RVQ(LLLL2)                                   
      BR=CMPP*RV(1)+CMP*RVQ(1)                                           
      ASQ=AR*AR                                                          
      BSQ=BR*BR                                                          
      AB=AR*BR                                                           
      absq=ab*ab                                                         
      ASBS=ASQ-BSQ                                                       
      CMPPD=CMPP*D                                                       
      CMPD=CMP*D                                                         
      NPP=NP+1                                                           
      TEMF=1.d0                                                          
      ipc=0                                                              
      DO 36 I=1,NP                                                       
      nmi=np-i
      if(i.eq.1) nmi=0
      IF(I.GT.NPH) GOTO 54                                                
      UPDOWN=1.d0                                                        
      IK=I                                                               
      GOTO 55                                                            
   54 UPDOWN=-1.d0                                                       
      IK=NPP-I                                                           
   55 CONTINUE                                                           
      IPN1=IK+(KOMP-1)*N1                                                
      SINTH=SNTH(IPN1)                                                   
      COSTH=CSTH(IPN1)*UPDOWN                                            
      tanth=sinth/costh                                                  
      EM=SINTH*EN*1.3d0                                                  
      MM=int(EM+1.d0)                                                         
      XM=MM                                                              
      MH=MM                                                              
      MM=2*MM                                                            
      DELFI=pi/XM                                                        
      r2nfdf=.5d0/enf                                                    
      deldel=delth*delfi                                                 
      IP=(KOMP-1)*NP1+IK                                                 
      IY=MMSAVE(IP)+1                                                    
      IF(TEST.LE..0625d0)GOTO 19                                         
      GX=GRXQ(IY)                                                        
      GY=-GRYQ(IY)                                                       
      GZ=UPDOWN*GRZQ(IY)                                                 
      grmag=gmag2(iy)                                                    
      GOTO 29                                                            
   19 GX=GRX(IY)                                                         
      GY=-GRY(IY)                                                        
      GZ=UPDOWN*GRZ(IY)                                                  
      grmag=gmag1(iy)                                                    
   29 COSSAV=(XLOS*GX+YLOS*GY+ZLOS*GZ)/GRMAG                             
      SUMJ=0.d0                                                          
      SOMJ=0.d0                                                          
      MPP=MM+1                                                           
      IY=IY-1                                                            
      DO 26 J=1,MM                                                       
      IF(J.GT.MH) GOTO 58                                                
      RTLEFT=1.d0                                                        
      JK=J                                                               
      GOTO 59                                                            
   58 RTLEFT=-1.d0                                                       
      JK=MPP-J                                                           
   59 CONTINUE                                                           
      IX=IY+JK                                                           
      IS=IX+(KOMP-1)*LLL1                                                
      SINFI=SNFI(IS)*RTLEFT                                              
      COSFI=CSFI(IS)                                                     
      STSF=SINTH*SINFI                                                   
      STCF=SINTH*COSFI                                                   
      IF(TEST.LE..0625d0)GOTO 39                                         
      IF(RVQ(IX).EQ.-1.d0) GOTO 26                                       
      GX=GRXQ(IX)                                                        
      GY=RTLEFT*GRYQ(IX)                                                 
      GZ=UPDOWN*GRZQ(IX)                                                 
      R=RVQ(IX)                                                          
      grmag=gmag2(ix)                                                    
      GOTO 49                                                            
   39 IF(RV(IX).EQ.-1.d0) GOTO 26                                        
      GX=GRX(IX)                                                         
      GY=RTLEFT*GRY(IX)                                                  
      GZ=UPDOWN*GRZ(IX)                                                  
      R=RV(IX)           
      grmag=gmag1(ix)                                                    
   49 COSGAM=(XLOS*GX+YLOS*GY+ZLOS*GZ)/GRMAG                             
      abscosgam=dabs(cosgam)
      xcoordr=stcf*r
      xneck=sefac*(cmp+comp*x1)
      ZZ=R*COSTH
      YY=R*COMP*STSF
      XX=CMPD+COMP*STCF*R                                                
      YSKY=XX*SINPH+YY*COSPH-cmpd*SINPH                                  
      ZSKY=-XX*CICP+YY*CISP+ZZ*SINI+CMPD*CICP                            
      rptsq=YSKY**2+ZSKY**2                                              
      rtstsq=absq/(BSQ+ASBS*(ZSKY**2/rptsq))                             
      if(mpage.ne.5) goto 174                                            
c
c  Next 2 lines delete self-eclipsed surface elements from images
c
      if(kff.eq.0.or.tests.lt.2.2562d-3.or.xcoordr.le.xneck) goto 194
      if(rptsq.le.rtstsq) GOTO 174                                        
  194 continue
      if(cosgam.gt.0.d0) goto 174                                        
      ipc=ipc+1                                                          
      yskp(ipc)=(xx-qfacd)*sinph+yy*cosph
      zskp(ipc)=(-xx+qfacd)*cicp+yy*cisp+zz*sini
      if(nspot.eq.0) goto 174                                            
      call frspot(komp,kspot,nmi,nspot,sinth,costh,sinfi,cosfi,delth,
     $delfi,temf)
      if(temf.eq.1.d0) goto 174                                          
      yskr=yskp(ipc)                                                     
      zskr=zskp(ipc)                                                     
      kstp=4                                                             
      cirf=.002d0                                                        
      do 179 kang=1,kstp                                           
      ang=twopi*dfloat(kang)/dfloat(kstp)                                             
      ipc=ipc+1                                                          
      yskp(ipc)=yskr+dsin(ang)*cirf                                      
      zskp(ipc)=zskr+dcos(ang)*cirf                                      
  179 continue                                                           
  174 continue                                                           
      if(sinsq.ge.sinsqe) goto 27
      PROD=COSSAV*COSGAM                                                 
      cossav=cosgam
c
c  Next line: Points in the neck region are not allowed to be horizon points.
c  
      if(kff.eq.0) goto 95
      if(xcoordr.gt.xneck) goto 27
      if(rptsq.le.0.98d0*rtstsq) GOTO 27                                        
   95 continue
c
c   Put the point in the horizon array if it's a horizon point.
c
      if(prod.gt.0.d0) goto 27
      ysky=xx*sinph+yy*cosph-cmpd*sinph
      zsky=-xx*cicp+yy*cisp+zz*sini+cmpd*cicp
      xrho=dsqrt(ysky**2+zsky**2)
      xtheta=dasin(zsky/xrho)
      if(ysky.lt.0.d0) goto 92
      xtheta=twopi+xtheta
      goto 93
   92 xtheta=pi-xtheta
   93 if(xtheta.ge.twopi) xtheta=xtheta-twopi
      L=L+1
      rho(L)=xrho
      theta(L)=xtheta
   27 continue
      IF(COSGAM.GE.0.d0) GOTO 26                                         
c
c  Surface elements for OC components that sky-project within the elliptical
c    horizon approximation are self-eclipsed. Skip to next element.
c
      if(kff.eq.0) goto 94
      if(tests.ge.2.2562d-3.and.rptsq.le.rtstsq.and.xcoordr.gt.xneck)
     $goto 26
   94 continue
      glogg=cmpp*glog1(ix)+cmp*glog2(ix) 
      if(ld.lt.0) call limdark(iband,ldabs,tld(is),glogg,xlimb,ylimb)
      CORFAC=1.d0                                                        
      if(mpage.ne.3) goto 933
      do 923 jn=1,nl                                                     
      Lspot(komp,jn)=0                                                   
  923 if(kks(komp,jn).eq.0) Lspot(komp,jn)=1                             
  933 continue
      IF(NSPOT.EQ.0) GOTO 640                                            
      call frspot(komp,kspot,nmi,nspot,sinth,costh,sinfi,cosfi,delth,
     $delfi,temf)
      IF(TEMF.EQ.1.d0) GOTO 640                                          
      TSP=TLD(IS)*TEMF                                                   
      if(ifat.eq.0) call planckint(tld(is),iband,xintlog,xintbase)
      if(ld.lt.0) call limdark(iband,ldabs,tsp,glogg,xlimb,ylimb)
      if(ifat.eq.0) call planckint(tsp,iband,xintlog,xintspot)
      IF(IFAT.EQ.0) GOTO 941                                             
      CALL atmx(TLD(IS),glogg,iband,xintlog,xintbase)                                        
      CALL atmx(TSP,glogg,iband,xintlog,xintspot)                                            
  941 CORFAC=xintspot/xintbase                  
  640 CONTINUE                                                           
      rit=1.d0                                                           
      if(ncl.eq.0) goto 818                                              
      do 815 icl=1,ncl                                                   
      opsfcl=opsf*fcl(icl)                                               
      call cloud(xlos,ylos,zlos,xx,yy,zz,xcl(icl),ycl(icl),zcl(icl),     
     $rcl(icl),wl,op1(icl),opsfcl,edens(icl),acm,encl(icl),cmpd,         
     $ri,dx,dens(icl),tau)                                               
      rit=rit*ri                                                         
  815 continue                                                           
  818 continue                                                           
      DARKEN=1.d0-xlimb+xlimb*abscosgam                                            
      if(ldabs.eq.2.and.cosgam.ne.0) darken=darken-ylimb*abscosgam*
     $dlog(abscosgam)                                
      if(ldabs.eq.3) darken=darken-ylimb*(1.d0-dsqrt(abscosgam))                   
      if(darken.lt.0.d0) darken=0.d0                                     
      DIF=rit*abscosgam*DARKEN*CORFAC*(CMP*SLUMP2(IX)+CMPP*SLUMP1(IX))      
      v=-r*(STCF*YLOS-stsf*XLOS)*COMP                                    
      if(ifphn.eq.1) goto 423                                            
      if(mpage.ne.3) goto 423                                            
      vflump=vfvuff*r*comp*costh                                         
      vcks=v*vfvuff                                                      
      vks=vcks+voc                                                       
      vksf=vks                                                           
      dvdr=vcks/r                                                        
      dvdth=vcks/tanth                                                   
      dvdfib=vfvuff*r*comp*(sinfi*ylos+cosfi*xlos)                       
c     dvdfic=dvdfib*sinth                                                
      difp=dif*deldel*renfsq                                             
c  dvdth and dvdfi (below) each need another term involving dr/d(theta)  
c    or dr/d(fi), that I will put in later. There will be a small loss   
c    of accuracy for distorted stars without those terms. See notes.     
      if(komp.eq.2) goto 422                                             
      binc=binc1                                                         
      binw=binw1                                                         
      do 1045 ifn=-nfm1,nfm1,2                                           
      dthf=dfloat(ifn)*r2nfdt                                            
      dvdfi=dvdfib*(sinth+costh*dthf)                                    
      do 1046 jfn=-nfm1,nfm1,2                                           
      if(nf.eq.1) goto 1047                                              
      dfif=dfloat(jfn)*r2nfdf*delfi                                      
      dvdth=-vflump*((cosfi-sinfi*dfif)*ylos-(sinfi+cosfi*dfif)*xlos)    
      dlr=0.d0                                                           
      vksf=vks+dvdr*dlr+dvdth*dthf+dvdfi*dfif                            
 1047 call linpro(komp,dvks1,hbarw1,tau1,emm1,count1,taug,emmg,fbin1,        
     $delv1) 
      if(inmin.lt.in1min) in1min=inmin                                   
      if(inmax.gt.in1max) in1max=inmax                                   
 1046 continue                                                           
 1045 continue                                                           
      goto 423                                                           
  422 continue                                                           
      binc=binc2                                                         
      binw=binw2                                                         
      do 1145 ifn=-nfm1,nfm1,2                                           
      dthf=dfloat(ifn)*r2nfdt                                            
      dvdfi=dvdfib*(sinth+costh*dthf)                                    
      do 1146 jfn=-nfm1,nfm1,2                                           
      if(nf.eq.1) goto 1147                                              
      dfif=dfloat(jfn)*r2nfdf*delfi                                      
      dvdth=-vflump*((cosfi-sinfi*dfif)*ylos-(sinfi+cosfi*dfif)*xlos)    
      dlr=0.d0                                                           
      vksf=vks+dvdr*dlr+dvdth*dthf+dvdfi*dfif                            
      ffi=dacos(cosfi)                                                   
      if(sinfi.lt.0.d0) ffi=twopi-ffi                                    
 1147 call linpro(komp,dvks2,hbarw2,tau2,emm2,count2,taug,emmg,fbin2,       
     $delv2) 
      if(inmin.lt.in2min) in2min=inmin                                   
      if(inmax.gt.in2max) in2max=inmax                                   
 1146 continue                                                           
 1145 continue                                                           
  423 continue                                                           
      DIFF=DIF*V                                                         
      SOMJ=SOMJ+DIFF                                                     
      SUMJ=SUMJ+DIF                                                      
   26 CONTINUE                                                           
      SOMJ=SOMJ*DELFI                                                    
      SUMJ=SUMJ*DELFI                                                    
      SOM=SOM+SOMJ                                                       
   36 SUM=SUM+SUMJ                                                       
      IF(SINSQ.GE.SINSQE) GOTO 75                                        
      LK=k                                                               
      if(L.lt.14) LK=L/2-1                                               
      CALL fourls(theta,rho,obser,L,LK,aa,bb)                                  
   75 continue
      IF(TEST.LE..0625d0) GOTO 118                                       
      SUMKUL=SUM*DELTH                                                   
      SOMKUL=SOM*DELTH                                                   
      xlimb=xh
      ylimb=yh
      KOMP=1                                                             
      ld=ld1
      nl=nl1                                                             
      ffc=ff1                                                            
      voc=vo1*vfvu                                                       
      NSPOT=NSP1                                                         
      IFAT=IFAT1                                                         
      EN=N1                                                              
      SAFTY=2.6d0*RV(LLL1)/EN                                            
      RMAX=RVQ(LLL2)+SAFTY                                               
      RMIN=RVQ(1)-SAFTY                                                  
      NPH=N1                                                             
      NP=2*N1                                                            
      nf=nf1                                                             
      GOTO 128                                                           
  118 xlimb=xc                                                               
      ylimb=yc                                                               
      KOMP=2                                                             
      ld=ld2
      nl=nl2                                                             
      ffc=ff2                                                            
      voc=vo2*vfvu                                                       
      NSPOT=NSP2                                                         
      IFAT=IFAT2                                                         
      SUMHOT=SUM*DELTH                                                   
      SOMHOT=SOM*DELTH                                                   
      if(inmax.gt.in1max) in1max=inmax                                   
      if(inmin.lt.in1min) in1min=inmin                                   
      EN=N2                                                              
      SAFTY=2.6d0*RVQ(LLL2)/EN                                           
      RMAX=RV(LLL1)+SAFTY                                                
      RMIN=RV(1)-SAFTY                                                   
      NPH=N2                                                             
      NP=2*N2                                                            
      nf=nf2                                                             
  128 DELTH=pih/EN                                                       
      ldabs=iabs(ld)
      enf=dfloat(nf)                                                     
      nfm1=nf-1                                                          
      renfsq=1.d0/(enf*enf)                                              
      r2nfdt=.5d0*delth/enf                                              
      vfvuff=vfvu*ffc                                                    
      SOM=0.d0                                                           
      SUM=0.d0                                                           
      NPP=NP+1                                                           
      TEMF=1.d0                                                          
      inmin=300000                                                        
      inmax=0                                                            
      DO 136 I=1,NP                                                      
      nmi=np-i
      if(i.eq.1) nmi=0
      IF(I.GT.NPH) GOTO 154                                              
      UPDOWN=1.d0                                                        
      IK=I                                                               
      GOTO 155                                                           
  154 UPDOWN=-1.d0                                                       
      IK=NPP-I                                                           
  155 CONTINUE                                                           
      IPN1=IK+(KOMP-1)*N1                                                
      SINTH=SNTH(IPN1)                                                   
      COSTH=CSTH(IPN1)*UPDOWN                                            
      tanth=sinth/costh                                                  
      EM=SINTH*EN*1.3d0                                                  
      MM=int(EM+1.d0)                                                         
      XM=MM                                                              
      MH=MM                                                              
      MM=2*MM                                                            
      DELFI=pi/XM                                                        
      deldel=delth*delfi                                                 
      SOMJ=0.d0                                                          
      SUMJ=0.d0                                                          
      SIGN=0.d0                                                          
      DRHO=1.d0                                                          
      MPP=MM+1                                                           
      DO 126 J=1,MM                                                      
      IF(J.GT.MH) GOTO 158                                               
      RTLEFT=1.d0                                                        
      JK=J                                                               
      GOTO 159                                                           
  158 RTLEFT=-1.d0                                                       
      JK=MPP-J                                                           
  159 CONTINUE                                                           
      IP=(KOMP-1)*NP1+IK                                                 
      IX=MMSAVE(IP)+JK                                                   
      IS=IX+LLL1*(KOMP-1)                                                
      SINFI=SNFI(IS)*RTLEFT                                              
      COSFI=CSFI(IS)                                                     
      STSF=SINTH*SINFI                                                   
      STCF=SINTH*COSFI                                                   
      IF(TEST.LE..0625d0)GOTO 139                                        
      IF(RV(IX).EQ.-1.d0) GOTO 126                                       
      GX=GRX(IX)                                                         
      GY=RTLEFT*GRY(IX)                                                  
      GZ=UPDOWN*GRZ(IX)                                                  
      R=RV(IX)                                                           
      grmag=gmag1(ix)                                                    
      GOTO 149                                                           
  139 IF(RVQ(IX).EQ.-1.d0) GOTO 126                                      
      GX=GRXQ(IX)                                                        
      GY=RTLEFT*GRYQ(IX)                                                 
      GZ=UPDOWN*GRZQ(IX)                                                 
      R=RVQ(IX)                                                          
      grmag=gmag2(ix)                                                    
  149 COSGAM=(XLOS*GX+YLOS*GY+ZLOS*GZ)/GRMAG                             
      abscosgam=dabs(cosgam)
      IF(COSGAM.LT.0.d0) GOTO 104                                        
      SIGN=0.d0                                                          
      OLSIGN=0.d0                                                        
      GOTO 126                                                           
  104 continue                                                     
      ZZ=R*COSTH                                                         
      YY=R*COMPP*STSF                                                    
      XX=CMPPD+COMPP*STCF*R                                              
      OLDIF=DIF                                                          
      CORFAC=1.d0                                                        
      glogg=cmp*glog1(ix)+cmpp*glog2(ix)
      if(ld.lt.0) call limdark(iband,ldabs,tld(is),glogg,xlimb,ylimb)
      if(mpage.ne.3) goto 833
      do 823 jn=1,nl                                                     
      Lspot(komp,jn)=0                                                   
  823 if(kks(komp,jn).eq.0) Lspot(komp,jn)=1 
  833 continue                            
      IF(NSPOT.EQ.0) GOTO 660                                            
      call frspot(komp,kspot,nmi,nspot,sinth,costh,sinfi,cosfi,delth,
     $delfi,temf)
      IF(TEMF.EQ.1.d0) GOTO 660                                      
      TSP=TLD(IS)*TEMF                                                   
      if(ifat.eq.0) call planckint(tld(is),iband,xintlog,xintbase)
      if(ld.lt.0) call limdark(iband,ldabs,tsp,glogg,xlimb,ylimb)
      if(ifat.eq.0) call planckint(tsp,iband,xintlog,xintspot)
      IF(IFAT.EQ.0) GOTO 661                                             
      CALL atmx(TLD(IS),glogg,iband,xintlog,xintbase)                                        
      CALL atmx(TSP,glogg,iband,xintlog,xintspot)                                            
  661 CORFAC=xintspot/xintbase                  
  660 CONTINUE                                                           
      rit=1.d0                                                           
      if(ncl.eq.0) goto 718                                              
      do 715 icl=1,ncl                                                   
      opsfcl=opsf*fcl(icl)                                               
      call cloud(xlos,ylos,zlos,xx,yy,zz,xcl(icl),ycl(icl),zcl(icl),     
     $rcl(icl),wl,op1(icl),opsfcl,edens(icl),acm,encl(icl),cmppd,        
     $ri,dx,dens(icl),tau)                                               
      rit=rit*ri                                                         
  715 continue                                                           
  718 continue                                                           
      DARKEN=1.d0-xlimb+xlimb*abscosgam
      if(ldabs.eq.2.and.cosgam.ne.0.d0) darken=darken-ylimb*abscosgam*
     $dlog(abscosgam)
      if(ldabs.eq.3) darken=darken-ylimb*(1.d0-dsqrt(abscosgam))
      if(darken.lt.0.d0) darken=0.d0
      DIF=rit*abscosgam*DARKEN*CORFAC*(CMPP*SLUMP2(IX)+CMP*SLUMP1(IX))      
      v=R*(STCF*YLOS-STSF*XLOS)*COMP                                     
      DIFF=DIF*V                                                         
      IF(SINSQ.GT.SINSQE) GOTO 63                                        
      OLSIGN=SIGN                                                        
      OLDRHO=DRHO                                                        
      YSKY=XX*SINPH+YY*COSPH-cmpd*SINPH                                  
      ZSKY=-XX*CICP+yy*CISP+ZZ*SINI+CMPD*CICP                            
      RRHO=dsqrt(ysky*ysky+zsky*zsky)                                    
      IF(RRHO.GT.RMAX) GOTO 63                                            
      IF(RRHO.LT.RMIN) GOTO 126                                           
      THET=dasin(ZSKY/RRHO)                                              
      IF(YSKY.LT.0.d0) GOTO 192                                          
      THET=twopi+THET                                                    
      GOTO 193                                                           
  192 THET=pi-THET                                                       
  193 IF(THET.GE.twopi) THET=THET-twopi                                  
      RHHO=0.d0                                                          
      DO 52 N=1,KK                                                       
      ENNN=N-1                                                           
      ENTHET=ENNN*THET                                                   
   52 RHHO=RHHO+AA(N)*dcos(ENTHET)+BB(N)*dsin(ENTHET)                    
      if(kff.le.0) goto 869
      rgap=dsqrt(absq/(BSQ+ASBS*(ZSKY/rrho)**2))+ot*safty
      if(rhho.gt.rgap) rhho=rgap
  869 continue
      SIGN=1.d0                                                          
      IF(RRHO.LE.RHHO) sign=-1.d0                                        
      if(mpage.eq.3) goto 861                                            
      DRHO=dabs(RRHO-RHHO)                                               
      IF((SIGN*OLSIGN).GE.0.d0) GOTO 60                                  
      SUMDR=DRHO+OLDRHO                                                  
      FACT=-(.5d0-DRHO/SUMDR)                                            
      IF(FACT.LT.0.d0) GOTO 198                                          
      RDIF=OLDIF                                                         
      GOTO 199                                                           
  198 RDIF=DIF                                                           
  199 CORR=FACT*RDIF*SIGN                                                
      CORRR=CORR*V                                                       
      SUMJ=SUMJ+CORR                                                     
      SOMJ=SOMJ+CORRR                                                    
   60 IF(SIGN.LT.0.d0) GOTO 126                                          
   63 SUMJ=SUMJ+DIF                                                      
      SOMJ=SOMJ+DIFF                                                     
      if(mpage.ne.5) goto 127                                            
      ipc=ipc+1                                                          
      yskp(ipc)=(xx-qfacd)*sinph+yy*cosph                                
      zskp(ipc)=(-xx+qfacd)*cicp+yy*cisp+zz*sini                         
      if(nspot.eq.0) goto 126                                            
      call frspot(komp,kspot,nmi,nspot,sinth,costh,sinfi,cosfi,delth,
     $delfi,temf)
      if(temf.eq.1.d0) goto 126                                          
      yskr=yskp(ipc)                                                     
      zskr=zskp(ipc)                                                     
      do 189 kang=1,kstp
      ang=twopi*dfloat(kang)/dfloat(kstp)
      ipc=ipc+1
      yskp(ipc)=yskr+dsin(ang)*cirf
      zskp(ipc)=zskr+dcos(ang)*cirf
  189 continue
      goto 126                                                           
  127 continue                                                           
      if(mpage.ne.3) goto 126                                            
      if(ifphn.eq.1) goto 126                                            
  861 vflump=vfvuff*r*comp*costh                                         
      vcks=v*vfvuff                                                      
      vks=vcks+voc                                                       
      vksf=vks                                                           
      dvdr=vcks/r                                                        
      dvdth=vcks/tanth                                                   
      dvdfib=vfvuff*r*comp*(sinfi*ylos+cosfi*xlos)                       
      difp=dif*deldel*renfsq                                             
      if(komp.eq.2) goto 452                                             
      binc=binc1                                                         
      binw=binw1                                                         
      do 1245 ifn=-nfm1,nfm1,2                                           
      dthf=dfloat(ifn)*r2nfdt                                            
      snthl=costh*dthf                                                   
      zz=r*(costh-sinth*dthf)                                            
      dvdfi=dvdfib*(sinth+costh*dthf)                                    
      do 1246 jfn=-nfm1,nfm1,2                                           
      if(nf.eq.1) goto 1247                                              
      dfif=dfloat(jfn)*r2nfdf*delfi                                      
      dlr=0.d0                                                           
      xx=cmppd+compp*r*snthl*(cosfi-sinfi*dfif)                          
      yy=r*compp*snthl*(sinfi+cosfi*dfif)                                
      ysky=(xx-cmpd)*sinph+yy*cosph                                      
      zsky=(cmpd-xx)*cicp+yy*cisp+zz*sini                                
      rrho=dsqrt(ysky*ysky+zsky*zsky)                                    
      if(rrho.lt.rhho) goto 1246                                         
      dvdth=-vflump*((cosfi-sinfi*dfif)*ylos-(sinfi+cosfi*dfif)*xlos)    
      vksf=vks+dvdr*dlr+dvdth*dthf+dvdfi*dfif                            
 1247 call linpro(komp,dvks1,hbarw1,tau1,emm1,count1,taug,emmg,fbin1,      
     $delv1) 
      if(inmax.gt.in1max) in1max=inmax                                   
      if(inmin.lt.in1min) in1min=inmin                                   
 1246 continue                                                           
 1245 continue                                                           
      goto 126                                                           
  452 continue                                                           
      binc=binc2                                                         
      binw=binw2                                                         
      do 1445 ifn=-nfm1,nfm1,2                                           
      dthf=dfloat(ifn)*r2nfdt                                            
      snthl=costh*dthf                                                   
      zz=r*(costh-sinth*dthf)                                            
      dvdfi=dvdfib*(sinth+costh*dthf)                                    
      do 1446 jfn=-nfm1,nfm1,2                                           
      if(nf.eq.1) goto 1447                                              
      dfif=dfloat(jfn)*r2nfdf*delfi                                      
      dvdth=-vflump*((cosfi-sinfi*dfif)*ylos-(sinfi+cosfi*dfif)*xlos)    
      dlr=0.d0                                                           
      xx=cmppd+compp*r*snthl*(cosfi-sinfi*dfif)                          
      yy=r*compp*snthl*(sinfi+cosfi*dfif)                                
      ysky=(xx-cmpd)*sinph+yy*cosph                                      
      zsky=(cmpd-xx)*cicp+yy*cisp+zz*sini                                
      rrho=dsqrt(ysky*ysky+zsky*zsky)                                    
      if(rrho.lt.rhho) goto 1446                                         
      vksf=vks+dvdr*dlr+dvdth*dthf+dvdfi*dfif                            
 1447 call linpro(komp,dvks2,hbarw2,tau2,emm2,count2,taug,emmg,fbin2,       
     $delv2) 
      if(inmax.gt.in2max) in2max=inmax                                   
      if(inmin.lt.in2min) in2min=inmin                                   
 1446 continue                                                           
 1445 continue                                                           
  126 CONTINUE                                                           
      SOMJ=SOMJ*DELFI                                                    
      SUMJ=SUMJ*DELFI                                                    
      SOM=SOM+SOMJ                                                       
  136 SUM=SUM+SUMJ                                                       
      if(mpage.eq.5) return                                              
      IF(TEST.LE..0625d0) GOTO 120                                       
      SOMHOT=SOM*DELTH                                                   
      SUMHOT=SUM*DELTH                                                   
      GOTO 121                                                           
  120 SUMKUL=SUM*DELTH                                                   
      SOMKUL=SOM*DELTH                                                   
  121 continue                                                           
      if(ifphn.eq.1) return                                              
      if(mpage.ne.3) return                                              
      in1min=in1min-marm1                                                
      in1max=in1max+marp1                                                
      in2min=in2min-marm2                                                
      in2max=in2max+marp2                                                
      if(nl1.eq.0) goto 3115                                             
      do 2912 i=in1min,in1max                                            
      fbin1(i)=1.d0-fbin1(i)/sumhot                                        
      if(count1(i).eq.0.d0) goto 2918                                    
      delv1(i)=delv1(i)/count1(i)                                        
      goto 2919                                                          
 2918 delv1(i)=binw1*(dfloat(i)-binc1)                                   
 2919 vdc=delv1(i)/clight                                               
      vfc=dsqrt((1.d0+vdc)/(1.d0-vdc))                                  
      delwl1(i)=wl*(vfc-1.d0)                                           
      wl1(i)=wl*vfc                                                     
      resf1(i)=(sl1*delwl1(i)+sc1)*fbin1(i)                              
 2912 continue                                                           
 3115 if(nl2.eq.0) return                                                
      do 2914 i=in2min,in2max                                            
      fbin2(i)=1.d0-fbin2(i)/sumkul                                      
      if(count2(i).eq.0.d0) goto 2917                                    
      delv2(i)=delv2(i)/count2(i)                                        
      goto 2920                                                          
 2917 delv2(i)=binw2*(dfloat(i)-binc2)                                   
 2920 vdc=delv2(i)/clight                                             
      vfc=dsqrt((1.d0+vdc)/(1.d0-vdc))                                  
      delwl2(i)=wl*(vfc-1.d0)                                           
      wl2(i)=wl*vfc                                                     
      resf2(i)=(sl2*delwl2(i)+sc2)*fbin2(i)                              
 2914 continue                                                           
      return                                                             
      END                                                                
      SUBROUTINE LUMP(GRX,GRY,GRZ,GRXQ,GRYQ,GRZQ,SLUMP1,SLUMP2,          
     $MMSAVE,ALB,tav,tavo,TPOLL,SBR,N1,N2,KOMP,IFAT,fr,snth,                   
     $TLD,GLUMP1,GLUMP2,XX1,XX2,YY1,YY2,ZZ1,ZZ2,xbol,ybol                
     $,GRV1,GRV2,SBR1B,SBR2B,RF,RFO,GMAG1,GMAG2,glog1,glog2,DINT,iband)                    
c   Version of June 8, 2017                                        
      implicit real*8 (a-h,o-z)                                          
      DIMENSION GRX(*),GRY(*),GRZ(*),GRXQ(*),GRYQ(*),grzq(*),            
     $SLUMP1(*),SLUMP2(*),MMSAVE(*),FR(*),SNTH(*),                       
     $TLD(*),GLUMP1(*),GLUMP2(*),XX1(*),XX2(*),YY1(*)                    
     $,YY2(*),ZZ1(*),ZZ2(*),GRV1(*),GRV2(*),RF(*),RFO(*),                
     $GMAG1(*),GMAG2(*),glog1(*),glog2(*)                                                  
      dimension message(2,4) 
      common /atmmessages/ message,kompcom 
      common /invar/ khdum,ipbdum,irtedm,nrefdm,irv1dm,irv2dm,mrefdm,     
     $ifs1dm,ifs2dm,icr1dm,icr2dm,ld1,ld2,ncl,jdphs,ipc,nr                      
      common /gpoles/ gplog1,gplog2 
      common /coflimbdark/ xld,yld,dm1,dm2,dm3,dm4
      kompcom=komp 
      ot=1.d0/3.d0
      IQ=(KOMP-1)*(N1+1)                                                 
      is=0
      if(komp.eq.2) is=mmsave(iq)
      xbolr=xbol
      ybolr=ybol
      dintl=dint
      pi=dacos(-1.d0)
      PIH=.5d0*PI                                                        
      TPOLE=10000.d0*TPOLL                                               
      cmp=dfloat(komp-1) 
      cmpp=dfloat(2-komp) 
      gplog=cmpp*gplog1+cmp*gplog2 
      ld=(2-komp)*ld1+(komp-1)*ld2
      ldabs=iabs(ld)
      ldo=(komp-1)*ld1+(2-komp)*ld2
      ldoabs=iabs(ldo)
      if(ld.lt.0) call limdark(iband,ldabs,tpole,gplog,xld,yld)
      if(ifat.eq.0) call planckint(tpole,iband,pollog,pint) 
      IF(IFAT.NE.0) CALL atmx(tpole,gplog,iband,pollog,pint)                            
      COMPP=dfloat(2*KOMP-3)                                             
      COMP=-COMPP                                                        
      N=(2-KOMP)*N1+(KOMP-1)*N2                                          
      NO=(2-KOMP)*N2+(KOMP-1)*N1                                         
      NOD=2*NO                                                           
      EN=dfloat(N)                                                       
      ENO=dfloat(NO)                                                     
      DELTHO=PIH/ENO                                                     
      CNST=ALB*DELTHO*SBR2B/SBR1B                                 
      IF(KOMP.EQ.2) CNST=ALB*DELTHO*SBR1B/SBR2B                   
      isss=0
      if(komp.eq.2) isss=mmsave(iq)
      DO 191 I=1,N                                                       
      IPN1=I+N1*(KOMP-1)                                                 
      SINTH=SNTH(IPN1)                                                   
      EM=SINTH*EN*1.3d0                                                  
      MM=int(EM+1.d0)                                                         
      IP=(KOMP-1)*(N1+1)+I                                               
      IY=MMSAVE(IP)                                                      
      DO 193 J=1,MM                                                      
      IX=IY+J                                                            
      iss=0
      if(komp.eq.1) iss=mmsave(n1+1)
      SUM=0.d0                                                           
      IF(FR(IX).EQ.0.d0) GOTO 193                                        
      DO 190 IOTH=1,NOD                                                  
      IOTHS=IOTH                                                         
      IF(IOTH.GT.NO) IOTHS=NOD-IOTH+1                                    
      IPNO=IOTHS+N1*(2-KOMP)                                             
      SINTHO=SNTH(IPNO)                                                  
      EMO=SINTHO*ENO*1.3d0                                               
      MMO=int(EMO+1.d0)                                                       
      MMOD=2*MMO                                                         
      IPO=(2-KOMP)*(N1+1)+IOTHS                                          
      IYO=MMSAVE(IPO)                                                    
      XMO=MMO                                                            
      DELFIO=PI/XMO                                                      
      DO 190 JOFI=1,MMOD                                                 
      JOFU=JOFI                                                          
      IF(JOFI.GT.MMO) JOFU=MMOD-JOFI+1                                   
      IXO=IYO+JOFU                                                       
c
c  Re-compute local bolometric limb darkenings on irradiating star if 
c   its 'LD' integer is negative (up to 297 continue).
c
      if(ldo.gt.0) goto 297
      tother=tavo*1.d4
      iss=iss+1
      if(nr.gt.1) tother=tld(iss)
      glog=cmp*glog1(ixo)+cmpp*glog2(ixo) 
      call limdark(0,ldo,tother,glog,xbolr,ybolr)
  297 continue
      IX1=IX                                                             
      IX2=IXO                                                            
      IF(KOMP.EQ.1) GOTO 200                                             
      IF(GLUMP1(IXO).EQ.0.d0) GOTO 184                                   
      IX1=IXO                                                            
      IX2=IX                                                             
      GOTO 201                                                           
  200 CONTINUE                                                           
      IF(GLUMP2(IXO).EQ.0.d0) GOTO 179                                   
  201 RTL1=1.d0                                                          
      RTL2=1.d0                                                          
      UPD1=1.d0                                                          
      UPD2=1.d0                                                          
      IF(KOMP.EQ.2) GOTO 22                                              
      IF(JOFI.GT.MMO) RTL2=-1.d0                                         
      IF(IOTH.GT.NO) UPD2=-1.d0                                          
      GOTO 23                                                            
   22 IF(JOFI.GT.MMO) RTL1=-1.d0                                         
      IF(IOTH.GT.NO) UPD1=-1.d0                                          
   23 CONTINUE                                                           
      GX2=GRXQ(IX2)                                                      
      GY2=GRYQ(IX2)*RTL2                                                 
      GZ2=GRZQ(IX2)*UPD2                                                 
      X1C=XX1(IX1)                                                       
      X2C=XX2(IX2)                                                       
      Y1C=YY1(IX1)*RTL1                                                  
      Y2C=YY2(IX2)*RTL2                                                  
      Z1C=ZZ1(IX1)*UPD1                                                  
      Z2C=ZZ2(IX2)*UPD2                                                  
      DX=(X2C-X1C)*COMP                                                  
      DY=(Y2C-Y1C)*COMP                                                  
      DZ=(Z2C-Z1C)*COMP                                                  
      DLRSQ=DX*DX+DY*DY+DZ*DZ                                            
      CSNUM2=(DX*GX2+DY*GY2+DZ*GZ2)*COMPP                                
      IF(CSNUM2.GE.0.d0) GOTO 190                                        
      GX1=GRX(IX1)                                                       
      GY1=GRY(IX1)*RTL1                                                  
      GZ1=GRZ(IX1)*UPD1                                                  
      CSNUM1=(DX*GX1+DY*GY1+DZ*GZ1)*COMP                                 
      IF(CSNUM1.GE.0.d0) GOTO 190                                        
      DMAG=dsqrt(DLRSQ)                                                  
      CSGM1=-CSNUM1/(DMAG*GMAG1(IX1))                                    
      CSGM2=-CSNUM2/(DMAG*GMAG2(IX2))                                    
      IF(KOMP.EQ.2) GOTO 181                                             
      DGAM2=1.d0-xbolr+xbolr*CSGM2                                         
      if(ldoabs.ne.2) goto 179                                               
      if(csgm2.eq.0.d0) goto 179                                         
      dgam2=dgam2-ybolr*csgm2*dlog(csgm2)                                 
      goto 147                                                           
  179 continue                                                           
      if(ldoabs.eq.3) dgam2=dgam2-ybolr*(1.d0-dsqrt(csgm2))                   
  147 if(dgam2.lt.0.d0) dgam2=0.d0                                       
      DSUM=GRV2(IXO)*GLUMP2(IXO)*RFO(IXO)*CSGM1*CSGM2*DGAM2/DLRSQ        
      GOTO 182                                                           
  181 continue
      DGAM1=1.d0-xbolr+xbolr*CSGM1                                         
      if(ldoabs.ne.2) goto 184                                               
      if(csgm1.eq.0.d0) goto 184                                         
      dgam1=dgam1-ybolr*csgm1*dlog(csgm1)                                 
      goto 148                                                           
  184 continue                                                           
      if(ldoabs.eq.3) dgam1=dgam1-ybolr*(1.d0-dsqrt(csgm1))                   
  148 if(dgam1.lt.0.d0) dgam1=0.d0                                       
      DSUM=GRV1(IXO)*GLUMP1(IXO)*RFO(IXO)*CSGM2*CSGM1*DGAM1/DLRSQ        
  182 CONTINUE                                                           
      SUM=SUM+DSUM*DELFIO                                                
  190 CONTINUE                                                           
c
c  If input integer LD is negative, interpolate irradiated star's local
c    bolometric limb darkening coefficients, xboll and yboll from log g
c    and T_eff. Then compute corresponding local bolometric integrated limb 
c    darkening factor, dintl (up to 197 continue). For positive LD,
c    apply fixed input limb darkening.
c
      if(ld.ge.0) goto 197
      isss=isss+1
      glog=cmpp*glog1(ix)+cmp*glog2(ix)
      tlocal=tav*1.d4
      if(nr.gt.1) tlocal=tld(isss)
      call limdark(0,ld,tlocal,glog,xboll,yboll)
      dintl=pi*(1.d0-ot*xboll)
      if(ldabs.eq.2) dintl=dintl+pi*2.d0*yboll/9.d0
      if(ldabs.eq.3) dintl=dintl-.2d0*pi*yboll
  197 continue
      RF(IX)=(CNST*SUM/(CMPP*GRV1(IX)+CMP*GRV2(IX)))/dintl+1.d0                
  193 CONTINUE                                                           
  191 CONTINUE                                                           
      DO 8 I=1,N                                                         
      IPN1=I+N1*(KOMP-1)                                                 
      SINTH=SNTH(IPN1)                                                   
      EM=SINTH*EN*1.3d0                                                  
      MM=int(EM+1.d0)                                                         
      IP=(KOMP-1)*(N1+1)+I                                               
      IY=MMSAVE(IP)                                                      
      DO 8 J=1,MM                                                        
      IS=IS+1                                                            
      IX=IY+J                                                            
      IF(FR(IX).EQ.0.d0) GOTO 8                                          
      glogg=cmpp*glog1(ix)+cmp*glog2(ix) 
      grv=cmpp*grv1(ix)+cmp*grv2(ix) 
      TNEW=TPOLE*dsqrt(dsqrt(GRV*RF(IX)))                                
      TLD(IS)=TNEW                                                       
      if(ld.lt.0) call limdark(iband,ldabs,tnew,glogg,xld,yld)
      if(ifat.eq.0) call planckint(tnew,iband,xintlog,xint) 
      if(ifat.ne.0) call atmx(tnew,glogg,iband,xintlog,xint)                            
      GRREFL=xint/pint                        
      IF(KOMP.EQ.1) GOTO 77                                              
      slump2(ix)=glump2(ix)*grrefl*sbr                                   
      GOTO 8                                                             
   77 slump1(ix)=glump1(ix)*grrefl*sbr                                   
    8 CONTINUE                                                           
      RETURN                                                             
      END                                                                
      Subroutine lightsp(phs,xincl,xh,xc,yh,yc,n1,n2,sumhot,sumkul,rv,
     $grx,gry,grz,rvq,grxq,gryq,grzq,mmsave,theta,rho,aa,bb,
     $somhot,somkul,d,snth,csth,snfi,csfi,tld,gmag1,gmag2,
     $obser,glump1,glump2)
c   Version of March 2, 2012                                        
      implicit real*8 (a-h,o-z)                                          
      DIMENSION RV(*),GRX(*),GRY(*),GRZ(*),RVQ(*),GRXQ(*),GRYQ(*),GRZQ(* 
     $),glump1(*),glump2(*),MMSAVE(*),THETA(*),
     $RHO(*),AA(*),BB(*)
      DIMENSION SNTH(*),CSTH(*),SNFI(*),CSFI(*),tld(*),gmag1(*), 
     $gmag2(*),obser(*) 
      dimension message(2,4) 
      common /kfac/ kff1,kff2,kfo1,kfo2
      common /klsp/ tratio,lsp
      common /xybol/ xbol1,xbol2,ybol1,ybol2
      common /atmmessages/ message,komp 
      common /coflimbdark/ xlimb,ylimb,dm1,dm2,dm3,dm4 
      COMMON /misc/ X1                                                          
      COMMON /NSPT/ nsp1,nsp2,kspot,kspev                               
      common /invar/ khdum,ipbdum,irtedm,nrefdm,irv1dm,irv2dm,mrefdm     
     $,ifs1dm,ifs2dm,icr1dm,icr2dm,ld1,ld2,ncl,jdphs,ipc,nr                      
      common /setest/ sefac                                              
      common /flpro/ vksf,binc,binw,difp,deldel,renfsq                   
      common /ipro/ nbins,nl,inmax,inmin,nf1,nf2                         
      pi=dacos(-1.d0)                                                    
      twopi=pi+pi                                                        
      pih=.5d0*pi                                                        
      dtr=pi/180.d0                                                      
      PHA=PHS*twopi                                                      
      K=6                                                                
      KK=K+1                                                             
      XINC=XINCL*dtr                                                     
      L=0                                                                
      TEST=(PHS-.5d0)**2                                                 
      SINI=dsin(XINC)                                                    
      COSPH=dcos(PHA)                                                    
      SINPH=dsin(PHA)                                                    
      SINSQ=SINPH**2                                                     
      COSI=dcos(XINC)                                                    
      NP1=N1+1                                                           
      NP2=N1+N2+2                                                        
      LLL1=MMSAVE(NP1)                                                   
      LLL2=MMSAVE(NP2)                                                   
      NPP2=NP2-1                                                         
      LL1=MMSAVE(N1)+1                                                   
      LL2=MMSAVE(NPP2)+1                                                 
      LLLL1=(LL1+LLL1)/2                                                 
      LLLL2=(LL2+LLL2)/2                                                 
      SINSQE=0.d0                                                        
      IF(SINI.GT.0.d0) SINSQE=((1.02d0*(RV(LLL1)+RVQ(LLL2))/D)**2        
     $-cosi**2)/SINI**2                                                  
      if(sinsqe.gt..96d0) sinsqe=.96d0
      CICP=COSI*COSPH                                                    
      CISP=COSI*SINPH                                                    
      XLOS=COSPH*SINI                                                    
      YLOS=-SINPH*SINI                                                   
      ZLOS=COSI                                                          
      suma=0.d0
      sumd=0.d0
      SUM=0.d0                                                           
      SOM=0.d0                                                           
      IF(TEST.LE..0625d0) GOTO 18                                        
      COMP=-1.d0                                                         
      CMP=1.d0                                                           
      COMPP=1.d0                                                         
      KOMP=2                                                             
      kff=kff2
      ld=ld2
      NSPOT=NSP2                                                         
      CMPP=0.d0                                                          
      xlimb=xc
      ylimb=yc
      xbol=xbol2
      ybol=ybol2
      EN=N2                                                              
      NPH=N2                                                             
      NP=2*N2                                                            
      nf=nf2                                                             
      GOTO 28                                                            
   18 xlimb=xh                                                               
      ylimb=yh                                                               
      xbol=xbol1
      ybol=ybol1
      COMP=1.d0                                                          
      KOMP=1                                                             
      kff=kff1
      ld=ld1
      NSPOT=NSP1                                                         
      CMP=0.d0                                                           
      COMPP=-1.d0                                                        
      CMPP=1.d0                                                          
      EN=N1                                                              
      NPH=N1                                                             
      NP=2*N1                                                            
      nf=nf1                                                             
   28 DELTH=pih/EN                                                       
      ldabs=iabs(ld)
      enf=dfloat(nf)                                                     
      renfsq=1.d0/(enf*enf)                                              
      AR=CMPP*RV(LLLL1)+CMP*RVQ(LLLL2)                                   
      BR=CMPP*RV(1)+CMP*RVQ(1)                                           
      ASQ=AR*AR                                                          
      BSQ=BR*BR                                                          
      AB=AR*BR                                                           
      absq=ab*ab                                                         
      ASBS=ASQ-BSQ                                                       
      CMPPD=CMPP*D                                                       
      CMPD=CMP*D                                                         
      NPP=NP+1                                                           
      TEMF=1.d0                                                          
      ipc=0                                                              
      DO 36 I=1,NP                                                       
      nmi=np-i
      if(i.eq.1) nmi=0
      IF(I.GT.NPH)GOTO 54                                                
      UPDOWN=1.d0                                                        
      IK=I                                                               
      GOTO 55                                                            
   54 UPDOWN=-1.d0                                                       
      IK=NPP-I                                                           
   55 CONTINUE                                                           
      IPN1=IK+(KOMP-1)*N1                                                
      SINTH=SNTH(IPN1)                                                   
      COSTH=CSTH(IPN1)*UPDOWN                                            
      EM=SINTH*EN*1.3d0                                                  
      MM=int(EM+1.d0)                                                         
      XM=MM                                                              
      MH=MM                                                              
      MM=2*MM                                                            
      DELFI=pi/XM                                                        
      deldel=delth*delfi                                                 
      IP=(KOMP-1)*NP1+IK                                                 
      IY=MMSAVE(IP)+1                                                    
      IF(TEST.LE..0625d0)GOTO 19                                         
      GX=GRXQ(IY)                                                        
      GY=-GRYQ(IY)                                                       
      GZ=UPDOWN*GRZQ(IY)                                                 
      grmag=gmag2(iy)                                                    
      GOTO 29                                                            
   19 GX=GRX(IY)                                                         
      GY=-GRY(IY)                                                        
      GZ=UPDOWN*GRZ(IY)                                                  
      grmag=gmag1(iy)                                                    
   29 COSSAV=(XLOS*GX+YLOS*GY+ZLOS*GZ)/GRMAG                             
      sumaj=0.d0
      sumdj=0.d0
      MPP=MM+1                                                           
      IY=IY-1                                                            
      DO 26 J=1,MM                                                       
      IF(J.GT.MH) GOTO 58                                                
      RTLEFT=1.d0                                                        
      JK=J                                                               
      GOTO 59                                                            
   58 RTLEFT=-1.d0                                                       
      JK=MPP-J                                                           
   59 CONTINUE                                                           
      IX=IY+JK                                                           
      IS=IX+(KOMP-1)*LLL1                                                
      SINFI=SNFI(IS)*RTLEFT                                              
      COSFI=CSFI(IS)                                                     
      STSF=SINTH*SINFI                                                   
      STCF=SINTH*COSFI                                                   
      IF(TEST.LE..0625d0)GOTO 39                                         
      IF(RVQ(IX).EQ.-1.d0) GOTO 26                                       
      GX=GRXQ(IX)                                                        
      GY=RTLEFT*GRYQ(IX)                                                 
      GZ=UPDOWN*GRZQ(IX)                                                 
      R=RVQ(IX)                                                          
      grmag=gmag2(ix)                                                    
      GOTO 49                                                            
   39 IF(RV(IX).EQ.-1.d0) GOTO 26                                        
      GX=GRX(IX)                                                         
      GY=RTLEFT*GRY(IX)                                                  
      GZ=UPDOWN*GRZ(IX)                                                  
      R=RV(IX)           
      grmag=gmag1(ix)                                                    
   49 COSGAM=(XLOS*GX+YLOS*GY+ZLOS*GZ)/GRMAG                             
      abscosgam=dabs(cosgam)
      xcoordr=stcf*r
      xneck=sefac*(cmp+comp*x1)
      ZZ=R*COSTH
      YY=R*COMP*STSF
      XX=CMPD+COMP*STCF*R                                                
      YSKY=XX*SINPH+YY*COSPH-cmpd*SINPH                                  
      ZSKY=-XX*CICP+YY*CISP+ZZ*SINI+CMPD*CICP                            
      rptsq=YSKY**2+ZSKY**2                                              
      rtstsq=absq/(BSQ+ASBS*(ZSKY**2/rptsq))                             
      if(sinsq.ge.sinsqe) goto 27
      PROD=COSSAV*COSGAM                                                 
      cossav=cosgam
c
c  Next line: Points in the neck region are not allowed to be horizon points.
c
      if(kff.eq.0) goto 95
      if(xcoordr.gt.xneck) goto 27
      if(rptsq.le.0.98d0*rtstsq) GOTO 27
   95 continue
c
c   Put the point in the horizon array if it's a horizon point.
c
c     if(rptsq.le.rtstsq.and.xcoordr.gt.xneck) goto 26
      if(prod.gt.0.d0) goto 27
      ysky=xx*sinph+yy*cosph-cmpd*sinph
      zsky=-xx*cicp+yy*cisp+zz*sini+cmpd*cicp
      xrho=dsqrt(ysky**2+zsky**2)
      xtheta=dasin(zsky/xrho)
      if(ysky.lt.0.d0) goto 92
      xtheta=twopi+xtheta
      goto 93
   92 xtheta=pi-xtheta
   93 if(xtheta.ge.twopi) xtheta=xtheta-twopi
      L=L+1
      rho(L)=xrho
      theta(L)=xtheta
   27 continue
      IF(COSGAM.GE.0.d0) GOTO 26                                         
c
c  Surface elements for OC components that sky-project within the elliptical
c    horizon approximation are self-eclipsed. Skip to next element.
c
      if(kff.eq.0) goto 94
      if(rptsq.le.rtstsq.and.xcoordr.gt.xneck)
     $goto 26
   94 continue
      if(nspot.gt.0) call frspot(komp,kspot,nmi,nspot,sinth,costh,sinfi,
     $cosfi,delth,delfi,temf)
      darkbol=1.d0-xbol+xbol*abscosgam
      if(ldabs.ne.2) goto 242
      if(abscosgam.eq.0.d0) goto 242
      darkbol=darkbol-ybol*abscosgam*dlog(abscosgam)
      goto 248
  242 continue
      if(ldabs.ne.3) goto 248
      darkbol=darkbol-ybol*(1.d0-dsqrt(abscosgam))
  248 if(darkbol.lt.0.d0) darkbol=0.d0
      difa=abscosgam*darkbol*temf**4*(cmp*glump2(ix)+cmpp*glump1(ix))
      ispol=1+(komp-1)*lll1
      tfac=(tld(is)/tld(ispol))**4
      difden=tfac*difa
      sumaj=sumaj+difa
      sumdj=sumdj+difden
   26 CONTINUE                                                           
      sumaj=sumaj*delfi
      sumdj=sumdj*delfi
      suma=suma+sumaj
      sumd=sumd+sumdj                                                  
   36 continue
      IF(SINSQ.GE.SINSQE) GOTO 75                                        
      LK=k                                                               
      if(L.lt.14) LK=L/2-1                                               
      CALL fourls(theta,rho,obser,L,LK,aa,bb)                                  
   75 IF(TEST.LE..0625d0) GOTO 118                                       
      sumakul=suma*delth
      sumdkul=sumd*delth
      SUMKUL=SUM*DELTH                                                   
      SOMKUL=SOM*DELTH                                                   
      xlimb=xh
      ylimb=yh
      xbol=xbol1
      ybol=ybol1
      KOMP=1                                                             
      ld=ld1
      NSPOT=NSP1                                                         
      EN=N1                                                              
      SAFTY=2.6d0*RV(LLL1)/EN                                            
      RMAX=RVQ(LLL2)+SAFTY                                               
      RMIN=RVQ(1)-SAFTY                                                  
      NPH=N1                                                             
      NP=2*N1                                                            
      nf=nf1                                                             
      GOTO 128                                                           
  118 xlimb=xc                                                               
      ylimb=yc                                                               
      xbol=xbol2
      ybol=ybol2
      KOMP=2                                                             
      ld=ld2
      NSPOT=NSP2                                                         
      sumahot=suma*delth
      sumdhot=sumd*delth
      SUMHOT=SUM*DELTH                                                   
      SOMHOT=SOM*DELTH                                                   
      EN=N2                                                              
      SAFTY=2.6d0*RVQ(LLL2)/EN                                           
      RMAX=RV(LLL1)+SAFTY                                                
      RMIN=RV(1)-SAFTY                                                   
      NPH=N2                                                             
      NP=2*N2                                                            
      nf=nf2                                                             
  128 DELTH=pih/EN                                                       
      ldabs=iabs(ld)
      enf=dfloat(nf)                                                     
      renfsq=1.d0/(enf*enf)                                              
      suma=0.d0
      sumd=0.d0
      SOM=0.d0                                                           
      SUM=0.d0                                                           
      NPP=NP+1                                                           
      TEMF=1.d0                                                          
      DO 136 I=1,NP                                                      
      nmi=np-i
      if(i.eq.1) nmi=0
      IF(I.GT.NPH) GOTO 154                                              
      UPDOWN=1.d0                                                        
      IK=I                                                               
      GOTO 155                                                           
  154 UPDOWN=-1.d0                                                       
      IK=NPP-I                                                           
  155 CONTINUE                                                           
      IPN1=IK+(KOMP-1)*N1                                                
      SINTH=SNTH(IPN1)                                                   
      COSTH=CSTH(IPN1)*UPDOWN                                            
      EM=SINTH*EN*1.3d0                                                  
      MM=int(EM+1.d0)                                                         
      XM=MM                                                              
      MH=MM                                                              
      MM=2*MM                                                            
      DELFI=pi/XM                                                        
      deldel=delth*delfi                                                 
      sumaj=0.d0
      sumdj=0.d0
      SIGN=0.d0                                                          
      DRHO=1.d0                                                          
      MPP=MM+1                                                           
      DO 126 J=1,MM                                                      
      IF(J.GT.MH) GOTO 158                                               
      RTLEFT=1.d0                                                        
      JK=J                                                               
      GOTO 159                                                           
  158 RTLEFT=-1.d0                                                       
      JK=MPP-J                                                           
  159 CONTINUE                                                           
      IP=(KOMP-1)*NP1+IK                                                 
      IX=MMSAVE(IP)+JK                                                   
      IS=IX+LLL1*(KOMP-1)                                                
      SINFI=SNFI(IS)*RTLEFT                                              
      COSFI=CSFI(IS)                                                     
      STSF=SINTH*SINFI                                                   
      STCF=SINTH*COSFI                                                   
      IF(TEST.LE..0625d0)GOTO 139                                        
      IF(RV(IX).EQ.-1.d0) GOTO 126                                       
      GX=GRX(IX)                                                         
      GY=RTLEFT*GRY(IX)                                                  
      GZ=UPDOWN*GRZ(IX)                                                  
      R=RV(IX)                                                           
      grmag=gmag1(ix)                                                    
      GOTO 149                                                           
  139 IF(RVQ(IX).EQ.-1.d0) GOTO 126                                      
      GX=GRXQ(IX)                                                        
      GY=RTLEFT*GRYQ(IX)                                                 
      GZ=UPDOWN*GRZQ(IX)                                                 
      R=RVQ(IX)                                                          
      grmag=gmag2(ix)                                                    
  149 COSGAM=(XLOS*GX+YLOS*GY+ZLOS*GZ)/GRMAG                             
      abscosgam=dabs(cosgam)
      IF(COSGAM.LT.0.d0) GOTO 104                                        
      SIGN=0.d0                                                          
      OLSIGN=0.d0                                                        
      GOTO 126                                                           
  104 continue                                                     
      ZZ=R*COSTH                                                         
      YY=R*COMPP*STSF                                                    
      XX=CMPPD+COMPP*STCF*R                                              
      darkbol=1.d0-xbol+xbol*abscosgam
      if(ldabs.ne.2) goto 142
      if(abscosgam.eq.0.d0) goto 142
      darkbol=darkbol-ybol*abscosgam*dlog(abscosgam)
      goto 148
  142 continue
      if(ldabs.eq.3) darkbol=darkbol-ybol*(1.d0-dsqrt(abscosgam))
  148 if(darkbol.lt.0.d0) darkbol=0.d0
      oldifa=difa
      oldifden=difden
      if(nspot.gt.0) call frspot(komp,kspot,nmi,nspot,sinth,costh,sinfi,
     $cosfi,delth,delfi,TEMF)
      darkbol=1.d0-xbol+xbol*abscosgam
      if(ldabs.ne.2) goto 342
      if(abscosgam.eq.0.d0) goto 342
      darkbol=darkbol-ybol*abscosgam*dlog(abscosgam)
      goto 348
  342 continue
      if(ldabs.ne.3) goto 348
      darkbol=darkbol-ybol*(1.d0-dsqrt(abscosgam))
  348 if(darkbol.lt.0.d0) darkbol=0.d0
      difa=abscosgam*darkbol*temf**4*(cmpp*glump2(ix)+cmp*glump1(ix))
      ispol=1+(komp-1)*lll1
      tfac=(tld(is)/tld(ispol))**4
      difden=tfac*difa
      IF(SINSQ.GT.SINSQE) GOTO 63                                        
      OLSIGN=SIGN                                                        
      OLDRHO=DRHO                                                        
      YSKY=XX*SINPH+YY*COSPH-cmpd*SINPH                                  
      ZSKY=-XX*CICP+yy*CISP+ZZ*SINI+CMPD*CICP                            
      RRHO=dsqrt(ysky*ysky+zsky*zsky)                                    
      IF(RRHO.GT.RMAX)GOTO 63                                            
      IF(RRHO.LT.RMIN)GOTO 126                                           
      THET=dasin(ZSKY/RRHO)                                              
      IF(YSKY.LT.0.d0) GOTO 192                                          
      THET=twopi+THET                                                    
      GOTO 193                                                           
  192 THET=pi-THET                                                       
  193 IF(THET.GE.twopi) THET=THET-twopi                                  
      RHHO=0.d0                                                          
      DO 52 N=1,KK                                                       
      ENNN=N-1                                                           
      ENTHET=ENNN*THET                                                   
   52 RHHO=RHHO+AA(N)*dcos(ENTHET)+BB(N)*dsin(ENTHET)                    
      if(kff.le.0) goto 869
      rgap=dsqrt(absq/(BSQ+ASBS*(ZSKY/rrho)**2))+safty
      if(rhho.gt.rgap) rhho=rgap
  869 continue
      SIGN=1.d0                                                          
      IF(RRHO.LE.RHHO) sign=-1.d0                                        
      DRHO=dabs(RRHO-RHHO)                                               
      IF((SIGN*OLSIGN).GE.0.d0) GOTO 60                                  
      SUMDR=DRHO+OLDRHO                                                  
      FACT=-(.5d0-DRHO/SUMDR)                                            
      IF(FACT.LT.0.d0) GOTO 198                                          
      rdifa=oldifa
      rdifden=oldifden
      GOTO 199                                                           
  198 continue
      rdifa=difa
      rdifden=difden
  199 continue
      corra=fact*rdifa*sign
      corrd=fact*rdifden*sign
      sumaj=sumaj+corra
      sumdj=sumdj+corrd
   60 IF(SIGN.LT.0.d0) GOTO 126                                          
   63 continue
      sumaj=sumaj+difa
      sumdj=sumdj+difden
  126 CONTINUE                                                           
      sumaj=sumaj*delfi
      sumdj=sumdj*delfi
      suma=suma+sumaj
      sumd=sumd+sumdj
  136 continue
      IF(TEST.LE..0625d0) GOTO 120                                       
      sumahot=suma*delth
      sumdhot=sumd*delth
      SOMHOT=SOM*DELTH                                                   
      SUMHOT=SUM*DELTH                                                   
      GOTO 121                                                           
  120 SUMKUL=SUM*DELTH                                                   
      SOMKUL=SOM*DELTH                                                   
      sumakul=suma*delth
      sumdkul=sumd*delth
  121 continue                                                           
      if(lsp.eq.1) tratio=dsqrt(dsqrt(sumahot/sumdhot))
      if(lsp.eq.2) tratio=dsqrt(dsqrt(sumakul/sumdkul))
      return                                                             
      END                                                                
      subroutine bandext (extin,iband,xx,ff,extout)
c  Version of June 7, 2013
      implicit real*8(a-h,o-z)
      dimension xw(15),xx(*),ff(*),extout(*)
      common /effwave/ effwvl(100)
      data xw(5),xw(6),xw(7),xw(8),xw(9),xw(10),xw(11),xw(12),xw(13)/
     $2.78d0,2.27d0,1.82d0,1.43d0,1.11d0,0.80d0,0.46d0,0.29d0,0.63d0/
c
c  The approximation functions applied here are from Cardelli, Clayton, & Mathis (1989)
c  BANDEXT computes ratios of extinctions in each of 93 photometric
c    bands to extinction in the Johnson V  band.  
c    It outputs extinction in stellar magnitudes for all the bands in the 
c    EXTOUT array. Adopted effective wavelengths were computed by W. Van Hamme 
c    from response curves in the literature except for bands 5 thtough 12
c    which are from Cardelli, Clayton, & Mathis. Those bands are: 
c 
c       5        U      Buser, R. 1978, Ang, 62, 411
c       6        B      Azusienis and Straizys 1969, Sov. Astron., 13, 316
c       7        V          "             "                "
c       8        R      Johnson, H.L. 1965, ApJ, 141, 923
c       9        I         "            "    "
c      10        J         "            "    "
c      11        K         "            "    "
c      12        L         "            "    "

c  Input 'extin' is extinction in band 'iband' and dimensioned variable 'ff' contains
c     the extinction ratios. 
c
c
      do 23 iyn=1,93
      xx(iyn)=1.d3/effwvl(iyn)
      if(iyn.ge.5.and.iyn.le.12) xx(iyn)=xw(iyn)
   23 continue
      rrrv=3.1d0
      do 18 ix=1,93
      xm=xx(ix)
      if(xm.ge.3.30d0) goto 14
      if(xm.le.1.10d0) goto 11
      yyy=xm-1.82d0
      ysq=yyy*yyy
      ycb=ysq*yyy
      y4=ycb*yyy
      y5=y4*yyy
      y6=y5*yyy
      y7=y6*yyy
      aofx=1.d0+.17699d0*yyy-.50447d0*ysq-.02427d0*ycb+.72085d0*y4
     $+.01979d0*y5-.77530d0*y6+.32999d0*y7
      bofx=1.41338d0*yyy+2.28305d0*ysq+1.07233d0*ycb-5.38434d0*y4
     $-.62251d0*y5+5.30260d0*y6-2.09002d0*y7
      goto 12
   11 continue
      aofx=.574d0*xm**1.61d0
      bofx=-.527d0*xm**1.61d0
      goto 12
   14 fa=0.d0
      fb=0.d0
      if(xm.lt.5.9d0) goto 15
      xmf2=(xm-5.9d0)**2
      xmf3=(xm-5.9d0)**3
      fa=-0.04473d0*xmf2-0.009779d0*xmf3
      fb=0.2130d0*xmf2+0.1207d0*xmf3
   15 continue
      aofx=1.752d0-0.316d0*xm-0.104d0/((xm-4.67d0)**2+0.341d0)+fa
      bofx=-3.090d0+1.825d0*xm+1.206d0/((xm-4.62d0)**2+0.263d0)+fb
   12 ff(ix)=aofx+bofx/rrrv
   18 continue
c
c  Calculate V-band extinction (as an intermediary scaling quantity)
c
      ev=extin/ff(iband)
c
c  Calculate extinction in all bands
c
      do 25 iyn=1,93
   25 extout(iyn)=ff(iyn)*ev
      return
      end
      subroutine cofprep(iab)
c   Version of June 22, 2015
      implicit real*8 (a-h,o-z)
      character*80 a
      parameter (nbmax=100,ntmax=80,ngmax=11)
      common /limco/ glog(ngmax),te(ntmax),nt(ngmax),kgrid(ngmax,ntmax),
     $imin(ngmax),imax(ngmax)
      common /coflim/ cof(ntmax,ngmax,5,nbmax)
      common /nummod/ nbp,ntemp,ngl
      goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),iab
 1    open(unit=20,file='limcof_bp_p10.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_p10.dat',status='old')
      goto 98
 2    open(unit=20,file='limcof_bp_p05.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_p05.dat',status='old')
      goto 98
 3    open(unit=20,file='limcof_bp_p03.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_p03.dat',status='old')
      goto 98
 4    open(unit=20,file='limcof_bp_p02.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_p02.dat',status='old')
      goto 98
 5    open(unit=20,file='limcof_bp_p01.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_p01.dat',status='old')
      goto 98
 6    open(unit=20,file='limcof_bp_p00.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_p00.dat',status='old')
      goto 98
 7    open(unit=20,file='limcof_bp_m01.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m01.dat',status='old')
      goto 98
 8    open(unit=20,file='limcof_bp_m02.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m02.dat',status='old')
      goto 98
 9    open(unit=20,file='limcof_bp_m03.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m03.dat',status='old')
      goto 98
 10   open(unit=20,file='limcof_bp_m05.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m05.dat',status='old')
      goto 98
 11   open(unit=20,file='limcof_bp_m10.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m10.dat',status='old')
      goto 98
 12   open(unit=20,file='limcof_bp_m25.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m25.dat',status='old')
      goto 98
 13   open(unit=20,file='limcof_bp_m20.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m20.dat',status='old')
      goto 98
 14   open(unit=20,file='limcof_bp_m25.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m25.dat',status='old')
      goto 98
 15   open(unit=20,file='limcof_bp_m30.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m30.dat',status='old')
      goto 98
 16   open(unit=20,file='limcof_bp_m35.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m35.dat',status='old')
      goto 98
 17   open(unit=20,file='limcof_bp_m40.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m40.dat',status='old')
      goto 98
 18   open(unit=20,file='limcof_bp_m45.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m45.dat',status='old')
      goto 98
 19   open(unit=20,file='limcof_bp_m50.dat',status='old')
      open(unit=21,file='limcof_bp_preamble_m50.dat',status='old')
 98   continue
      read(21,*) nmod,nmodtr,nbp,ntemp,ngl
      do 20 i=1,ngl
      read(21,*) j,te(i),glog(i)
 20   continue
      ngl1=ngl+1
      do 21 i=ngl1,ntemp
      read(21,*) j,te(i)
 21   continue
      close (21)
      do 23 j=1,ngl
      nt(j)=0
      do 23 i=1,ntemp
      kgrid(j,i)=0
 23   continue
      do 24 ii=1,500
      if(ii.ge.nmodtr) goto 25
      read(20,80,end=99) t,g
      goto 26
 25   read(20,81,end=99) t,g
 26   continue
  80  format(16X,f5.0,13X,F3.1)
  81  format(16X,f6.0,13X,f3.1)
      read(20,77) a
  77  format(80A)
      do 27 j=1,ngl
      if(dabs(g-glog(j)).lt.1.0d-04) goto 28
 27   continue
 28   continue
      do 29 i=1,ntemp
      if(dabs(t-te(i)).lt.1.d-04) goto 30
 29   continue
 30   continue
      kgrid(j,i)=kgrid(j,i)+1
      nt(j)=nt(j)+1
      nb1=nbp+1
      do 31 k=1,nb1
      read(20,82) cof(i,j,1,k),cof(i,j,2,k),cof(i,j,3,k),
     $cof(i,j,4,k),cof(i,j,5,k)
 31   continue
      do 333 j=1,ngl
      do 334 i=1,ntemp
 334  if(kgrid(j,i).ne.0) goto 335
 335  imin(j)=i
      do 336 i=ntemp,1,-1
 336  if(kgrid(j,i).ne.0) goto 337
 337  imax(j)=i
 333  continue
  82  format(25X,f7.3,12X,f6.3,1X,f6.3,12X,f6.3,1X,f6.3)
      read(20,77,end=99) a
      read(20,77,end=99) a
  24  continue
  99  continue
      close (20)
      return
      end
      subroutine limdark(iband,ld,t,g,x,y)
      implicit real*8 (a-h,o-z)
      parameter (nbmax=100, ntmax=80, ngmax=11)
      dimension xx(2)
      common /limco/ glog(ngmax),te(ntmax),nt(ngmax),kgrid(ngmax,ntmax),
     $imin(ngmax),imax(ngmax)
      common /coflim/ cof(ntmax,ngmax,5,nbmax)
      common /nummod/ nbp,ntemp,ngl
      ifil=iband+1
      ico1=4
      if(ld.eq.1) ico1=1
      if(ld.eq.2) ico1=2
      ico2=ico1+1
      if(ld.eq.1) ico2=1
      call binnum(te,ntemp,t,i)
      call binnum(glog,ngl,g,j)
      if(j.eq.0) j=1
      if(j.eq.ngl) j=ngl-1
      if(i.eq.0) i=1
      if(i.eq.ntemp) i=ntemp-1
      if(kgrid(j,i).eq.0) goto 11
      if(kgrid(j,i+1).eq.0) goto 11
      if(kgrid(j+1,i+1).eq.0) goto 11
      if(kgrid(j+1,i).eq.0) goto 11
      b=(g-glog(j))/(glog(j+1)-glog(j))
      a=(t-te(i))/(te(i+1)-te(i))
      do 10 kk=1,2
      if(kk.eq.1) kj=ico1
      if(kk.eq.2) kj=ico2
      xx(kk)=(1.d0-a)*(1.d0-b)*cof(i,j,kj,ifil)+
     $a*(1.d0-b)*cof(i+1,j,kj,ifil)
      xx(kk)=xx(kk)+a*b*cof(i+1,j+1,kj,ifil)+
     $(1.d0-a)*b*cof(i,j+1,kj,ifil)
  10  continue
      if(ld.eq.1) xx(2)=0.d0
      goto 16
  11  continue
      diftsav=1.0d+09
      do 12 jj=j,11
      if(kgrid(jj,i).ne.0) goto 15
  12  continue
  15  continue
      ib=imin(jj)
      ie=imax(jj)
      do 13 ii=ib,ie
      dift=dabs(t-te(ii))
      if(dift.gt.diftsav) goto 13
      ik=ii
      diftsav=dift
  13  continue
      do 14 kk=1,2
      if(kk.eq.1) kj=ico1
      if(kk.eq.2) kj=ico2
      xx(kk)=cof(ik,jj,kj,ifil)
  14  continue
      if(ld.eq.1) xx(2)=0.d0
  16  continue
      x=xx(1)
      y=xx(2)
      if(x.eq.0.d0) stop
      return 
      end
      subroutine planckint(t,ifil,ylog,y) 
      implicit real*8 (a-h,o-z) 
c  Made by WVH
c  Version of April 7, 2012 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c  IMPORTANT README 
c  This subroutine returns the log10 (ylog) of a Planck central 
c  intensity (y), as well as the Planck central intensity (y) itself. 
c  The subroutine ONLY WORKS FOR TEMPERATURES GREATER THAN OR EQUAL 
c  501 K (log10(Teff) = 2.7) OR LOWER THAN 501,187 K (log10(Teff) = 5.7). 
c  For teperatures outside this range, the program stops and prints a message. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      parameter (nbmax=100)
      parameter (nplcofmax=50*nbmax)
      dimension pl(10) 
      common /invar/ id1,id2,id3,id4,id5,id6,id7,id8,id9, 
     $id10,id11,ld1,ld2,id13,id14,id15,nr 
      common /planckleg/ plcof(nplcofmax)
      common /coflimbdark/ xld,yld,dm1,dm2,dm3,dm4 
      common /atmmessages/ message(2,4),komp
      ld=(2-komp)*ld1+(komp-1)*ld2
      ldabs=iabs(ld)
      tlog=dlog10(t)
      if(tlog.lt.2.7d0) goto 11 
      if(tlog.ge.3.279d0) goto 1 
      tb=2.7d0 
      te=3.279d0 
      ibin=1 
      goto 5 
   1  if(tlog.ge.3.740d0) goto 2 
      tb=3.279d0 
      te=3.74d0 
      ibin=2 
      goto 5 
   2  if(tlog.ge.4.301d0) goto 3 
      tb=3.74d0 
      te=4.301d0 
      ibin=3 
      goto 5 
   3  if(tlog.ge.5.d0) goto 4 
      tb=4.301d0 
      te=5.0d0 
      ibin=4 
      goto 5 
   4  if(tlog.gt.5.7d0) goto 11 
      tb=5.0d0 
      te=5.7d0 
      ibin=5 
   5  continue 
      ib=(ifil-1)*50+(ibin-1)*10 
      phas=(tlog-tb)/(te-tb) 
      call legendre(phas,pl,10) 
      y=0.d0 
      do 6 j=1,10 
      jj=j+ib 
   6  y=y+pl(j)*plcof(jj) 
      dark=1.d0-xld/3.d0 
      if(ldabs.eq.2) dark=dark+yld/4.5d0 
      if(ldabs.eq.3) dark=dark-0.2d0*yld 
      ylog=y-dlog10(dark)-0.49714987269413d0 
      y=10.d0**ylog 
      return 
  11  continue 
      write(6,80) 
      stop 
  80  format('Program stopped in PLANCKINT, 
     $T outside 501 - 501,187 K range.') 
      end 
      SUBROUTINE OLUMP(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,SLUMP1,SLUMP2,  
     $MMSAVE,GREXP,ALB,RB,TPOLL,SBR,SUMM,N1,N2,KOMP,IFAT,xoth,yoth,D,       
     $SNTH,CSTH,SNFI,CSFI,tld,glump1,glump2,glog1,glog2,grv1,grv2,iband)                             
c   Version of September 8, 2009                                        
      implicit real*8 (a-h,o-z)                                          
      DIMENSION RV(*),GRX(*),GRY(*),GRZ(*),RVQ(*),GRXQ(*),GRYQ(*),GRZQ(* 
     $),SLUMP1(*),SLUMP2(*),MMSAVE(*),F(3),W(3),SNTH(*),CSTH(*),         
     $SNFI(*),CSFI(*),tld(*),glump1(*),glump2(*),glog1(*),glog2(*),                         
     $grv1(*),grv2(*) 
      dimension message(2,4) 
      common /atmmessages/ message,kompcom 
      common /invar/ khdum,ipbdum,irtedm,nrefdm,irv1dm,irv2dm,mrefdm,
     $ifs1dm,ifs2dm,icr1dm,icr2dm,ld1,ld2,ncl,jdphs,ipc,nr
      common /gpoles/ gplog1,gplog2 
      kompcom=komp 
      ld=(2-komp)*ld1+(komp-1)*ld2
      ldabs=iabs(ld)
      IQ=(KOMP-1)*(N1+1)                                                 
      is=0
      if(iq.gt.0) IS=(KOMP-1)*MMSAVE(IQ)                                             
      FP=7.957747d-2                                                     
      pi=dacos(-1.d0)                                            
      pih=.5d0*pi
      pi32=1.5d0*pi
      F(1)=.1127017d0                                                    
      F(2)=.5d0                                                          
      F(3)=.8872983d0                                                    
      W(1)=.277777777777777d0                                            
      W(2)=.444444444444444d0                                            
      W(3)=.277777777777777d0                                            
      TPOLE=10000.d0*TPOLL                                               
      cmp=dfloat(komp-1) 
      cmpp=dfloat(2-komp) 
      gplog=cmpp*gplog1+cmp*gplog2 
      if(ifat.eq.0) call planckint(tpole,iband,pollog,pint) 
      IF(IFAT.NE.0) CALL atmx(tpole,gplog,iband,pollog,pint)                            
      COMPP=dfloat(2*KOMP-3)                                             
      COMP=-COMPP                                                        
      CMPD=CMP*D                                                         
      CMPPD=CMPP*D                                                       
      N=(2-KOMP)*N1+(KOMP-1)*N2                                          
      ENN=(15.d0+xoth)*(1.d0+GREXP)/(15.d0-5.d0*xoth)                          
      NP=N1+1+(2-KOMP)*(N2+1)                                            
      NPP=N1*(KOMP-1)+(NP-1)*(2-KOMP)                                    
      LL=MMSAVE(NPP)+1                                                   
      LLL=MMSAVE(NP)                                                     
      LLLL=(LL+LLL)/2                                                    
      AR=RV(LLL)*CMP+RVQ(LLL)*CMPP                                       
      BR=RV(LLLL)*CMP+RVQ(LLLL)*CMPP                                     
      CR=RV(1)*CMP+RVQ(1)*CMPP                                           
      BOA=BR/AR                                                          
      BOAL=1.d0-BOA*BOA                                                  
      BOC2=(BR/CR)**2                                                    
      CC=1.d0/(1.d0-.25d0*ENN*(1.d0-BOA**2)*(.9675d0-.3008d0*BOA))       
      HCN=.5d0*CC*ENN                                                    
      DF=1.d0-xoth/3.d0                                                     
      if(ldabs.eq.2) df=df+2.d0*yoth/9.d0
      if(ldabs.eq.3) df=df-.2d0*yoth
      EN=dfloat(N)                                                       
      DO 8 I=1,N                                                         
      IPN1=I+N1*(KOMP-1)                                                 
      SINTH=SNTH(IPN1)                                                   
      COSTH=CSTH(IPN1)                                                   
      EM=SINTH*EN*1.3d0                                                  
      MM=EM+1.d0                                                         
      IP=(KOMP-1)*(N1+1)+I                                               
      IY=MMSAVE(IP)                                                      
      DO 8 J=1,MM                                                        
      IS=IS+1                                                            
      ix=iy+j
      STCF=SINTH*CSFI(IS)                                                
      STSF=SINTH*SNFI(IS)                                                
      IF(KOMP.EQ.1) GOTO 39                                              
      IF(RVQ(IX).EQ.-1.d0) GOTO 8                                        
      GX=GRXQ(IX)                                                        
      GY=GRYQ(IX)                                                        
      GZ=GRZQ(IX)                                                        
      R=RVQ(IX)                                                          
      GOTO 49                                                            
   39 IF(RV(IX).EQ.-1.d0)GOTO 8                                          
      GX=GRX(IX)                                                         
      GY=GRY(IX)                                                         
      GZ=GRZ(IX)                                                         
      R=RV(IX)                                                           
   49 GRMAG=dsqrt(GX*GX+GY*GY+GZ*GZ)                                     
      ZZ=R*COSTH                                                         
      YY=R*COMP*STSF                                                     
      XX=CMPD+COMP*STCF*R                                                
      XXREF=(CMPPD+COMPP*XX)*COMPP                                       
      GRAV=cmpp*grv1(ix)+cmp*grv2(ix)                                         
      TLOCAL=TPOLE*dsqrt(dsqrt(GRAV))                                    
      DIST=dsqrt(XXREF*XXREF+YY*YY+ZZ*ZZ)                                
      RMX=dasin(.5d0*(BR+CR)/DIST)                                       
      XCOS=XXREF/DIST                                                    
      YCOS=YY/DIST                                                       
      ZCOS=ZZ/DIST                                                       
      COSINE=(XCOS*GX+YCOS*GY+ZCOS*GZ)/GRMAG                             
      RC=PIH-dacos(COSINE)                                               
      AH=RC/RMX                                                          
      RP=dabs(AH)                                                        
      IF(AH.LE..99999d0) GOTO 22                                         
      P=1.d0                                                             
      GOTO 16                                                            
   22 IF(AH.GE.-.99999d0) GOTO 24                                        
      ALBEP=0.d0                                                         
      GOTO 19                                                            
   24 SUM=0.d0                                                           
      FIST=dasin(RP)                                                     
      FII=PIH-FIST                                                       
      DO 15 IT=1,3                                                       
      FE=FII*F(IT)+FIST                                                  
      PAR=1.d0-(RP/dsin(FE))**2                                          
      RPAR=dsqrt(PAR)                                                    
      SUM=PAR*RPAR*W(IT)+SUM                                             
   15 CONTINUE                                                           
      FTRI=(1.d0-xoth)*RP*dsqrt(1.d0-RP**2)+.666666666666666d0*xoth*        
     $fii-.666666666666667d0*xoth*sum*fii                                      
      FSEC=(PIH+FIST)*DF                                                 
      P=(FTRI+FSEC)/(PI*DF)                                              
      IF(COSINE.LT.0.d0) P=1.d0-P                                        
      RTF=dsqrt(1.d0-AH**2)                                              
      DENO=PI32-3.d0*(AH*RTF+dasin(AH))                                  
      IF(DENO.NE.0.d0) GOTO 117                                          
      ABAR=1.d0                                                          
      GOTO 116                                                           
  117 ABAR=2.d0*RTF**3/DENO                                              
  116 COSINE=dcos(PIH-RMX*ABAR)                                          
   16 COSQ=1.d0/(1.d0+(YY/XXREF)**2)                                     
      COT2=(ZZ/XXREF)**2                                                 
      Z=BOAL/(1.d0+BOC2*COT2)                                            
      E=CC-HCN*COSQ*Z                                                    
      ALBEP=ALB*E*P                                                      
   19 IF(COSINE.LE.0.d0) ALBEP=0.d0                                      
      TNEW=TLOCAL*dsqrt(dsqrt(1.d0+(FP*SUMM/(DIST*DIST*GRAV))*           
     $cosine*rb*ALBEP))                                                  
      TLD(IS)=TNEW                                                       
      glogg=cmpp*glog1(ix)+cmp*glog2(ix)
      if(ifat.eq.0) call planckint(tnew,iband,xintlog,xint) 
      if(ifat.ne.0) CALL atmx(TNEW,glogg,iband,xintlog,xint) 
      grrefl=xint/pint                                                
      IF(KOMP.EQ.1) GOTO 77                                              
      slump2(ix)=glump2(ix)*grrefl*sbr                                   
      GOTO 8                                                             
   77 slump1(ix)=glump1(ix)*grrefl*sbr                                   
    8 CONTINUE                                                           
      RETURN                                                             
      END                                                                
      SUBROUTINE MODLOG(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVE,FR1,FR2 
     $,HLD,RM,POTH,POTC,GR1,GR2,ALB1,ALB2,N1,N2,F1,F2,MOD,XINCL,THE,     
     $MODE,SNTH,CSTH,SNFI,CSFI,GRV1,GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,GLUMP1  
     $,GLUMP2,CSBT1,CSBT2,GMAG1,GMAG2,glog1,glog2)                                   
c    Version of July 23, 2009                                        
      implicit real*8 (a-h,o-z)                                          
      DIMENSION RV(*),GRX(*),GRY(*),GRZ(*),RVQ(*),GRXQ(*),GRYQ(*),GRZQ   
     $(*),MMSAVE(*),FR1(*),FR2(*),HLD(*),GRV1(*),GRV2(*),XX1(*),YY1(*),  
     $ZZ1(*),XX2(*),YY2(*),ZZ2(*),GLUMP1(*),GLUMP2(*),CSBT1(*),CSBT2(*)  
     $,GMAG1(*),GMAG2(*),glog1(*),glog2(*)                                                 
      DIMENSION DRR(4),RES(2),ANS(2),LX(2),MX(2)                         
      DIMENSION SNTH(*),CSTH(*),SNFI(*),CSFI(*)                          
      common /kfac/ kff1,kff2,kfo1,kfo2                                  
      common /setest/ sefac                                              
      common /ardot/ dperdt,hjd,hjd0,perr                                
      COMMON /FLVAR/ PSHIFT,DP,EF,EFC,ECOS,perr0,PHPER,pconsc,pconic,
     $PHPERI,VSUM1,VSUM2,VRA1,VRA2,VKM1,VKM2,VUNIT,vfvu,trc,qfacd        
      COMMON /ECCEN/ E,A,PERIOD,VGA,SINI,VF,VFAC,VGAM,VOL1,VOL2,IFC      
      COMMON /INVAR/ KH,IPBDUM,IRTE,NREF,IRVOL1,IRVOL2,mref,ifsmv1,      
     $ifsmv2,icor1,icor2,ld1,ld2,ncl,jdphs,ipc,nr                                
   95 FORMAT(' WARNING: ALTHOUGH COMPONENT 2 DOES NOT EXCEED ITS LIMITIN 
     $G LOBE AT THE END OF ECLIPSE, IT DOES EXCEED THE LOBE AT PERIASTRO 
     $N')                                                                
   99 FORMAT(' SPECIFIED ECLIPSE DURATION INCONSISTENT WITH OTHER PARAME 
     $TERS')                                                             
      perr=perr0+dperdt*(hjd-hjd0)                        
      DP=1.d0-E                                                          
      MOD=(MODE-2)**2                                                    
      IF(MODE.EQ.1) GR2=GR1                                               
      IF(MODE.EQ.1) ALB2=ALB1                                             
      IF(MOD.EQ.1) POTC=POTH                                             
      MD4=(MODE-5)**2                                                    
      MD5=(2*MODE-11)**2                                                 
      call ellone(f1,dp,rm,xl1,po1cr,xl2,omo1)                           
      sefac=.8712d0                                                      
      doc=(po1cr-poth)/(po1cr-omo1)                                      
      if(doc.gt.0.d0) sefac=.201d0*doc*doc-.386d0*doc+.8712d0            
      RMR=1.d0/RM                                                        
      CALL ELLONE(F2,DP,RMR,XL1,po2c,XL2,omo2)                           
      po2cr=rm*po2c+(1.d0-rm)*.5d0                                       
      if(md4.eq.1) poth=po1cr                                            
      if(md5.eq.1) potc=po2cr                                            
      kff1=0                                                             
      kff2=0                                                             
      if(poth.lt.po1cr) kff1=1                                           
      if(potc.lt.po2cr) kff2=1                                           
      kfo1=0                                                             
      kfo2=0                                                             
      if(e.ne.0.d0) goto 100                                             
      if(f1.ne.1.d0) goto 105                                            
      if(poth.lt.omo1) kfo1=1                                            
  105 if(f2.ne.1.d0) goto 100                                            
      if(potc.lt.omo1) kfo2=1                                            
  100 continue                                                           
      SINI=dsin(.017453292519943d0*XINCL)                                
      VF=50.61455d0/PERIOD                                               
      VFAC=VF*A                                                          
      VGAM=VGA*VUNIT/VFAC                                                
      VFVU=VFAC                                                    
      IFC=2                                                              
      IF(e.eq.0.d0) IFC=1
      TRC=1.570796326794897d0-perr                                       
   39 if(TRC.LT.0.d0) TRC=TRC+6.283185307179586d0                        
      if(trc.lt.0.d0) goto 39                                            
   40 if(trc.ge.6.283185307179586d0) trc=trc-6.283185307179586d0         
      if(trc.ge.6.283185307179586d0) goto 40                             
      HTRC=.5d0*TRC                                                      
      IF(dabs(1.570796326794897d0-HTRC).LT.7.d-6) GOTO 101               
      IF(dabs(4.712388980384690d0-HTRC).LT.7.d-6) GOTO 101               
      ECAN=2.d0*datan(dsqrt((1.d0-E)/(1.d0+E))*dtan(HTRC))               
      GOTO 103                                                           
  101 ECAN=3.141592653589793d0                                           
  103 XMC=ECAN-E*dsin(ECAN)                                              
      IF(XMC.LT.0.d0) XMC=XMC+6.283185307179586d0                        
      PHPER=1.d0-XMC/6.283185307179586d0                                 
      call conjph(e,perr,pshift,trsc,tric,econsc,econic,xmsc,xmic,
     $pconsc,pconic)
   38 if(pconsc.ge.1.d0) pconsc=pconsc-1.d0                                 
      if(pconsc.ge.1.d0) goto 38                                          
   41 if(pconsc.lt.0.d0) pconsc=pconsc+1.d0                                 
      if(pconsc.lt.0.d0) goto 41                                          
   68 if(pconic.ge.1.d0) pconic=pconic-1.d0                                 
      if(pconic.ge.1.d0) goto 68                                          
   71 if(pconic.lt.0.d0) pconic=pconic+1.d0                                 
      if(pconic.lt.0.d0) goto 71                                          
      PHPERI=PHPER+pconsc                                                 
      EF=1.d0-E*E                                                        
      EFC=dsqrt(EF)                                                      
      ECOS=E*dcos(perr)                                                  
      IF(MODE.NE.-1) RETURN                                              
      if(kh.eq.17) goto 241                                              
      if((kh-12)**2.eq.1) goto 241                                       
      if((kh-12)**2.eq.4) goto 241                                       
      IF((KH-11)**2.LE.1) GOTO 241                                       
      IF((2*KH-41)**2.EQ.81) GOTO 241                                    
      RETURN                                                             
  241 CONTINUE                                                           
      EFCC=dsqrt((1.d0-E)/(1.d0+E))                                      
      THER=THE*6.283185307179586d0                                       
      DELTR=.001d0                                                       
      DTR1=0.d0                                                          
      DTR2=0.d0                                                          
      VOLTOL=5.d-6                                                       
      DXMTOL=5.d-6                                                       
      TR0=1.570796326794897d0-perr                                       
      HTR0=.5d0*TR0                                                      
      IF((1.570796326794897d0-dabs(HTR0)).LT.7.d-6) GOTO 201             
      IF((4.712388980384690d0-dabs(HTR0)).LT.7.d-6) GOTO 201             
      ECAN0=2.d0*datan(dsqrt((1.d0-E)/(1.d0+E))*dtan(HTR0))              
      GOTO 203                                                           
  201 ECAN0=3.141592653589793d0                                          
  203 XM0=ECAN0-E*dsin(ECAN0)                                            
      XM1=XM0-THER*(1.d0-.2d0*E)                                         
      XM2=XM0+THER*(1.d0-.2d0*E)                                         
      CALL KEPLER(XM1,E,DUM,TRR1)                                        
      CALL KEPLER(XM2,E,DUM,TRR2)                                        
  160 TRR1=TRR1+DTR1                                                     
      TRR2=TRR2+DTR2                                                     
      DO 161 IB=1,3                                                      
      TR1=TRR1                                                           
      TR2=TRR2                                                           
      IF(IB.EQ.2) TR1=TRR1+DELTR                                         
      IF(IB.EQ.3) TR2=TRR2+DELTR                                         
      IF(TR1.GT.TR0) TR0=TR0+6.283185307179586d0                         
      IF(TR0.GT.TR2) TR2=TR2+6.283185307179586d0                         
      DS1=EF/(1.d0+E*dcos(TR1))                                          
      DS2=EF/(1.d0+E*dcos(TR2))                                          
      TRE1=(TR0-TR1)/6.283185307179586d0                                 
      TRE2=(TR2-TR0)/6.283185307179586d0                                 
      CALL DURA(F2,XINCL,RM,DS1,TRE1,POTR,RA)                            
      CALL VOLUME(VS1,RM,POTR,DS1,F2,N2,N1,2,RV,GRX,GRY,GRZ,RVQ,GRXQ     
     $,GRYQ,GRZQ,MMSAVE,FR1,FR2,HLD,SNTH,CSTH,SNFI,CSFI,SUMMD,SMD,GRV1,  
     $GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2,GMAG1,      
     $GMAG2,glog1,glog2,GR1,1)                                                       
      CALL DURA(F2,XINCL,RM,DS2,TRE2,POTR,RA)                            
      CALL VOLUME(VS2,RM,POTR,DS2,F2,N2,N1,2,RV,GRX,GRY,GRZ,RVQ,GRXQ     
     $,GRYQ,GRZQ,MMSAVE,FR1,FR2,HLD,SNTH,CSTH,SNFI,CSFI,SUMMD,SMD,GRV1,  
     $GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2,GMAG1,      
     $GMAG2,glog1,glog2,GR2,1)                                                       
      IF(IB.NE.1) GOTO 185                                               
      ECAN1=2.d0*datan(dsqrt((1.d0-E)/(1.d0+E))*dtan(.5d0*TR1))          
      ECAN2=2.d0*datan(dsqrt((1.d0-E)/(1.d0+E))*dtan(.5d0*TR2))          
      POTC=POTR                                                          
      DTHE=DS2                                                           
      DVOL=VS2-VS1                                                       
      XM1=ECAN1-E*dsin(ECAN1)                                            
      XM2=ECAN2-E*dsin(ECAN2)                                            
      IF(XM1.LT.0.d0) XM1=XM1+6.283185307179586d0                        
      IF(XM2.LT.0.d0) XM2=XM2+6.283185307179586d0                        
      DXM=XM2-XM1-2.d0*THER                                              
      DDMDN1=-EFCC*(1.d0-E*dcos(ECAN1))*dcos(.5d0*ECAN1)**2/             
     $dcos(.5d0*tr1)**2                                                  
      DDMDN2=EFCC*(1.d0-E*dcos(ECAN2))*dcos(.5d0*ECAN2)**2/              
     $dcos(.5d0*tr2)**2                                                  
  185 CONTINUE                                                           
      IF(IB.NE.2) GOTO 162                                               
      DRR(1)=(VS2-VS1-DVOL)/DELTR                                        
      DRR(2)=DDMDN1                                                      
  162 CONTINUE                                                           
      IF(IB.NE.3) GOTO 161                                               
      DRR(3)=(VS2-VS1-DVOL)/DELTR                                        
      DRR(4)=DDMDN2                                                      
  161 CONTINUE                                                           
      RES(1)=-DVOL                                                       
      RES(2)=-DXM                                                        
      CALL DMINV(DRR,2,DUMM,LX,MX)                                       
      CALL DGMPRD(DRR,RES,ANS,2,2,1)                                     
      DTR1=ANS(1)                                                        
      DTR2=ANS(2)                                                        
      IF(dabs(DTR1).GT.VOLTOL) GOTO 160                                  
      IF(dabs(DTR2).GT.DXMTOL) GOTO 160                                  
      POTH=9999.99d0                                                     
      RMR=1.d0/RM                                                        
      CALL ELLONE(F2,DTHE,RMR,XLA,OM1,XL2,OM2)                           
      OM1=RM*OM1+(1.d0-RM)*.5d0                                          
      IF(POTC.LT.OM1) GOTO 22                                            
      IF(RA.LE.XLA) GOTO 28                                              
   22 WRITE(6,99)                                                        
      RETURN                                                             
   28 CONTINUE                                                           
      IF(E.NE.0.d0) CALL VOLUME(VTHE,RM,POTC,DTHE,F2,N2,N1,2,RV,GRX,     
     $GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVE,FR1,FR2,HLD,SNTH,CSTH,SNFI,CSFI, 
     $SUMMD,SMD,GRV1,GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,    
     $GLUMP2,GMAG1,GMAG2,glog1,glog2,GR2,1)                                          
      IF(E.NE.0.d0) CALL VOLUME(VTHE,RM,POTC,DP,F2,N2,N1,2,RV,GRX,       
     $GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVE,FR1,FR2,HLD,SNTH,CSTH,SNFI,CSFI, 
     $SUMMD,SMD,GRV1,GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,    
     $GLUMP2,GMAG1,GMAG2,glog1,glog2,GR2,2)                                          
      CALL ELLONE(F2,DP,RMR,XLD,OMP,XL2,OM2)                             
      OMP=RM*OMP+(1.d0-RM)*.5d0                                          
      IF(POTC.LT.OMP) WRITE(6,95)                                        
      RETURN                                                             
      END                                                                
      SUBROUTINE LCR(RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,MMSAVE,FR1,FR2, 
     $hld,SLUMP1,SLUMP2,RM,POTH,POTC,N1,N2,F1,F2,D,HLUM,CLUM,xh,xc,yh, 
     $yc,GR1,GR2,SM1,SM2,TPOLH,TPOLC,SBRH,SBRC,IFAT1,IFAT2,TAVH,TAVC,    
     $alb1,alb2,xbol1,xbol2,ybol1,ybol2,vol1,vol2,snth,csth,snfi,csfi,   
     $tld,glump1,glump2,xx1,xx2,yy1,yy2,zz1,zz2,dint1,dint2,grv1,grv2,   
     $csbt1,csbt2,rftemp,rf1,rf2,gmag1,gmag2,glog1,glog2,mode,iband)                   
c  Version of March 19, 2008                                         
      implicit real*8 (a-h,o-z)                                          
      DIMENSION RV(*),GRX(*),GRY(*),GRZ(*),RVQ(*),GRXQ(*),GRYQ(*),GRZQ(* 
     $),SLUMP1(*),SLUMP2(*),MMSAVE(*),FR1(*),FR2(*),HLD(*),SNTH(*),      
     $CSTH(*),SNFI(*),CSFI(*),TLD(*),GLUMP1(*),GLUMP2(*),XX1(*),XX2(*)   
     $,YY1(*),YY2(*),ZZ1(*),ZZ2(*),GRV1(*),GRV2(*),RFTEMP(*),RF1(*),     
     $RF2(*),CSBT1(*),CSBT2(*),GMAG1(*),GMAG2(*),glog1(*),glog2(*)                         
      dimension message(2,4) 
      common /atmmessages/ message,komp 
      common /coflimbdark/ xld,yld,xldmean1,xldmean2,yldmean1,yldmean2 
      COMMON /DPDX/ DPDX1,DPDX2,PHSV,PCSV                                
      COMMON /ECCEN/ E,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,ifc 
      COMMON /SUMM/ SUMM1,SUMM2                                          
      COMMON /INVAR/ KHDUM,IPB,IRTE,NREF,IRVOL1,IRVOL2,mref,ifsmv1,      
     $ifsmv2,icor1,icor2,ld1,ld2,ncl,jdphs,ipc,nr                                
      common /gpoles/ gplog1,gplog2 
      nn1=n1 
      VL1=VOL1                                                           
      VL2=VOL2                                                           
      DP=1.d0-E                                                          
      IF(IRVOL1.EQ.1) GOTO 88                                            
      CALL VOLUME(VL1,RM,POTH,DP,F1,nn1,N1,1,RV,GRX,GRY,GRZ,RVQ,GRXQ,     
     $GRYQ,GRZQ,MMSAVE,FR1,FR2,HLD,SNTH,CSTH,SNFI,CSFI,SUMM1,SM1,GRV1,   
     $GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2,GMAG1,GMAG2 
     $,glog1,glog2,GR1,1)                                                            
      IF(E.EQ.0.d0) GOTO 88                                              
      POTHD=PHSV                                                         
      IF(IFC.EQ.2) POTHD=PHSV+DPDX1*(1.d0/D-1.d0/(1.d0-E))               
      CALL VOLUME(VL1,RM,POTHD,D,F1,nn1,N1,1,RV,GRX,GRY,GRZ,RVQ,GRXQ,     
     $GRYQ,GRZQ,MMSAVE,FR1,FR2,HLD,SNTH,CSTH,SNFI,CSFI,SUMM1,SM1,GRV1,   
     $GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2,GMAG1,GMAG2 
     $,glog1,glog2,GR1,IFC)                                                          
   88 CONTINUE                                                           
      IF(IRVOL2.EQ.1) GOTO 86                                            
      CALL VOLUME(VL2,RM,POTC,DP,F2,N2,N1,2,RV,GRX,GRY,GRZ,RVQ,GRXQ,     
     $GRYQ,GRZQ,MMSAVE,FR1,FR2,HLD,SNTH,CSTH,SNFI,CSFI,SUMM2,SM2,GRV1,   
     $GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2,GMAG1,GMAG2 
     $,glog1,glog2,GR2,1)                                                            
      IF(E.EQ.0.d0) GOTO 86                                              
      POTCD=PCSV                                                         
      IF(IFC.EQ.2) POTCD=PCSV+DPDX2*(1.d0/D-1.d0/(1.d0-E))               
      CALL VOLUME(VL2,RM,POTCD,D,F2,N2,N1,2,RV,GRX,GRY,GRZ,RVQ,GRXQ,     
     $GRYQ,GRZQ,MMSAVE,FR1,FR2,HLD,SNTH,CSTH,SNFI,CSFI,SUMM2,SM2,GRV1,   
     $GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2,GMAG1,GMAG2 
     $,glog1,glog2,GR2,IFC)                                                          
   86 CONTINUE                                                           
      TPOLH=TAVH*dsqrt(dsqrt(SM1/SUMM1))                                 
      TPOLC=TAVC*dsqrt(dsqrt(SM2/SUMM2))                                 
      g1=gmag1(1)                                                        
      g2=gmag2(1)                                                        
      IF(MODE.EQ.1)TPOLC=TPOLH*dsqrt(dsqrt((G2/G1)**GR1))                
      IF(MODE.EQ.1)TAVC=TPOLC/dsqrt(dsqrt(SM2/SUMM2))                    
      tph=10000.d0*tpolh                                                 
      tpc=10000.d0*tpolc                                                 
      komp=1 
      xld=xh 
      yld=yh 
      ldabs=iabs(ld1)
      if(ld1.lt.0.and.ifat1.eq.0) call limdark(iband,ldabs,tph,gplog1,
     $xld,yld)
      if(ifat1.eq.0) call planckint(tph,iband,xintlog1,xint1) 
      IF(IFAT1.NE.0) CALL atmx(tph,gplog1,iband,xintlog1,xint1)                               
      call lum(hlum,xh,yh,tpolh,n1,n1,1,sbrh,rv,rvq,glump1,           
     $glump2,glog1,glog2,grv1,grv2,mmsave,summ1d,fr1,sm1d,ifat1,vold,rm,       
     $poth,f1,d,snth,iband)                                                         
      komp=2 
      xld=xc 
      yld=yc 
      ldabs=iabs(ld2)
      if(ld2.lt.0.and.ifat2.eq.0) call limdark(iband,ldabs,tpc,gplog2,
     $xld,yld)
      if(ifat2.eq.0) call planckint(tpc,iband,xintlog2,xint2) 
      IF(IFAT2.NE.0) CALL atmx(tpc,gplog2,iband,xintlog2,xint2)                               
      sbrc=sbrh*xint2/xint1            
      call lum(clum,xc,yc,tpolc,n2,n1,2,sbrt,rv,rvq,glump1,           
     $glump2,glog1,glog2,grv1,grv2,mmsave,summ2d,fr2,sm2d,ifat2,vold,rm,       
     $potc,f2,d,snth,iband)                                                         
      IF(IPB.EQ.1) SBRC=SBRT                                             
      IF(MODE.GT.0)CLUM=CLUM*SBRC/SBRT                                   
      IF(MODE.LE.0)SBRC=SBRT                                             
      if(mref.eq.2) goto 30                                              
      radrat=(vol1/vol2)**(1.d0/3.d0) 
      ratbol=radrat**2*(tavh/tavc)**4 
      rb=1.d0/ratbol                                                     
      xoth=xc
      yoth=yc
      if(ld2.lt.0) xoth=xldmean2
      if(ld2.lt.0) yoth=yldmean2
      call olump(rv,grx,gry,grz,rvq,grxq,gryq,grzq,slump1,slump2,mmsave  
     $,gr1,alb1,rb,tpolh,sbrh,summ1,n1,n2,1,ifat1,xoth,yoth,d,snth,        
     $csth,snfi,csfi,tld,glump1,glump2,glog1,glog2,grv1,grv2,iband)                                 
      rb=ratbol                                                          
      xoth=xh
      yoth=yh
      if(ld1.lt.0) xoth=xldmean1
      if(ld1.lt.0) yoth=yldmean1
      call olump(rv,grx,gry,grz,rvq,grxq,gryq,grzq,slump1,slump2,mmsave  
     $,gr2,alb2,rb,tpolc,sbrc,summ2,n1,n2,2,ifat2,xoth,yoth,d,snth,        
     $csth,snfi,csfi,tld,glump1,glump2,glog1,glog2,grv1,grv2,iband)                                 
      return                                                             
   30 continue                                                           
      sbr1b=tpolh**4/dint1                                               
      sbr2b=tpolc**4/dint2                                               
      LT=N1+1                                                            
      IMAX1=MMSAVE(LT)                                                   
      DO 80 I=1,IMAX1                                                    
      RFTEMP(I)=1.d0                                                     
   80 RF1(I)=1.d0                                                        
      LT=N1+N2+2                                                         
      IMAX2=MMSAVE(LT)                                                   
      DO 81 I=1,IMAX2                                                    
   81 RF2(I)=1.d0                                                        
      DO 93 NR=1,NREF                                                    
      CALL LUMP(GRX,GRY,GRZ,GRXQ,GRYQ,GRZQ,SLUMP1,SLUMP2,MMSAVE,         
     $alb1,tavh,tavc,tpolh,sbrh,n1,n2,1,ifat1,fr1,snth,                         
     $tld,glump1,glump2,xx1,xx2,yy1,yy2,zz1,zz2,xbol2,ybol2,grv1,        
     $grv2,sbr1b,sbr2b,rftemp,rf2,gmag1,gmag2,glog1,glog2,dint1,iband)                     
      CALL LUMP(GRX,GRY,GRZ,GRXQ,GRYQ,GRZQ,SLUMP1,SLUMP2,MMSAVE,         
     $ALB2,tavc,tavh,TPOLC,SBRC,N1,N2,2,IFAT2,fr2,snth,                         
     $tld,glump1,glump2,xx1,xx2,yy1,yy2,zz1,zz2,xbol1,ybol1,             
     $grv1,grv2,sbr1b,sbr2b,rf2,rf1,gmag1,gmag2,glog1,glog2,dint2,iband)                   
      DO 70 I=1,IMAX1                                                    
   70 RF1(I)=RFTEMP(I)                                                   
   93 CONTINUE                                                           
      RETURN                                                             
      END                                                                
      subroutine atmx(t,g,ifil,xintlog,xint)   
      implicit real*8 (a-h,o-z)   
c Version of April 7, 2012  
      parameter (nbmax=100,ntmax=80,ngmax=11)
      parameter (ntmax2=2*ntmax)
      dimension yy(4),pha(4),tte(2) 
      dimension y2(ntmax),tef(ntmax),y(ntmax)
      common /abung/abun(19),glog(ngmax)   
      common /arrayspline/ grand(ntmax2,ngmax,nbmax),
     $tem(ntmax2,ngmax,nbmax),yf(ntmax2,ngmax,nbmax)
      common /effwave/ effwvl(nbmax)
      common /numtemp/ kend1(ngmax,nbmax),kend2(ngmax,nbmax)
      common /ramprange/ tlowtol,thightol,glowtol,ghightol   
      common /atmmessages/ message(2,4),komp   
      common /coflimbdark/ xld,yld,dm1,dm2,dm3,dm4 
      common /invar/ id1,id2,id3,id4,id5,id6,id7,id8,id9, 
     $id10,id11,ld1,ld2,id13,id14,id15,nr 
      ld=(2-komp)*ld1+(komp-1)*ld2
      ldabs=iabs(ld)
      tlog=dlog10(t)   
      trec=1.d0/t   
      tlow=3500.d0-tlowtol   
      if(t.le.tlow) goto 66   
      thigh=50000.d0+thightol   
      fractol=thightol/50000.d0   
      glow=0.d0-glowtol   
      if(g.le.glow) goto 77   
      ghigh=5.d0+ghightol   
      if(g.ge.ghigh) goto 78   
      tt=t   
      gg=g   
      if(g.ge.0.d0) goto 11   
      gg=0.d0   
      goto 12   
  11  if(g.le.5.d0) goto 12   
      gg=5.d0   
  12  continue   
ccccccccccccccccccccccccccccccccccccccccccccccccccccc   
c The following is for 4-point interpolation in log g.   
ccccccccccccccccccccccccccccccccccccccccccccccccccccc   
      m=4   
      ifreturn=0   
      call binnum(glog,11,g,j)   
      k=min(max(j-(m-1)/2,1),12-m)   
      if(g.le.0.d0) j=1   
  10  continue   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
      do 4 ii=1,m
      nbin=1
      njj=k+ii-1   
      ken1=kend1(njj,ifil)
      nsp=ken1
      if (tem(81,njj,ifil).ne.0.d0) nbin=2
      tb=tem(1,njj,ifil)
      te=tem(ken1,njj,ifil)
      telog=dlog10(te)
      tblog=dlog10(tb)
      do 551 iii=1,ken1
      tef(iii)=tem(iii,njj,ifil)
      y2(iii)=grand(iii,njj,ifil)
      y(iii)=yf(iii,njj,ifil)
 551  continue
      if(tt.lt.te.or.nbin.eq.1) goto 55
      ken2=kend2(njj,ifil)
      nsp=ken2
      tb=tem(81,njj,ifil)
      te=tem(80+ken2,njj,ifil)
      do 552 iii=1,ken2
      iik=80+iii
      tef(iii)=tem(iik,njj,ifil)
      y2(iii)=grand(iik,njj,ifil)
      y(iii)=yf(iik,njj,ifil)
 552  continue
 55   continue
      thigh=te+fractol*te   
      pha(ii)=(tt-tb)/(te-tb)   
      yy(ii)=0.d0   
      call splinterpol(tef,y,y2,nsp,tt,yy(ii))
      if(pha(ii).lt.0.d0) call splinterpol(tef,y,y2,nsp,tb,yy(ii))
      if(pha(ii).ge.0.d0) goto 4   
      tlow=tb-tlowtol   
      if(ld.lt.0) call limdark(ifil,ldabs,tlow,g,xld,yld)
      call planckint(tlow,ifil,yylow,dum)   
      if(t.ge.tlow) goto 424   
      if(ld.lt.0) call limdark(ifil,ldabs,t,g,xld,yld)
      call planckint(t,ifil,yy(ii),dum)   
      goto 4   
  424 continue   
      tlowmidlog=0.5d0*dlog10(tb*tlow)   
      wvlmax=10.d0**(6.4624d0-tlowmidlog)   
      if(effwvl(ifil).lt.wvlmax) goto 425   
      tblog=dlog10(tb)   
      tlowlog=dlog10(tlow)   
      slope=(yy(ii)-yylow)/(tblog-tlowlog)   
      yy(ii)=yylow+slope*(tlog-tlowlog)   
      goto 4   
  425 continue   
      tbrec=1.d0/tb   
      tlowrec=1.d0/tlow   
      slope=(yy(ii)-yylow)/(tbrec-tlowrec)   
      yy(ii)=yylow+slope*(trec-tlowrec)   
   4  continue   
c Next, do a m-point Lagrange interpolation.   
      xintlog=0.d0   
      do 501 ii=1,m   
      xnum=1.d0   
      denom=1.d0   
      nj=k+ii-1   
      do 500 iij=1,m   
      njj=k+iij-1   
      if(ii.eq.iij) goto 500   
      xnum=xnum*(gg-glog(njj))   
      denom=denom*(glog(nj)-glog(njj))   
 500  continue   
      xintlog=xintlog+yy(ii)*xnum/denom   
 501  continue   
cccccccccccccccccccccccccccccccccccccccccccccccc   
cccccccccccccccccccccccccccccccccccccccccccccccc   
c  Check if a ramp function will be needed, or if we are   
c  close to the border and need to interpolate between less   
c  than 4 points.   
ccccccccccccccccccccccccccccccccccccccccccccccccc   
      if(g.lt.0.d0) goto 7   
      if(g.gt.5.d0) goto 9   
      if(t.lt.3500.d0) goto 99   
      if(pha(1).le.1.d0) goto 99   
      if(ifreturn.eq.1) goto 99  
      if(j.eq.1) goto 5   
      if(pha(3).gt.1.d0) goto 5   
      k=k+1   
      if(pha(2).gt.1.d0) goto 41    
  42  continue   
      if(k.gt.8) m=12-k   
      ifreturn=1   
      goto 10   
  41  continue   
      if(j.lt.10) goto 5   
      k=k+1   
      goto 42   
ccccccccccccccccccccccccccccccccccccccccccccccccc   
   5  continue   
      do 61 kik=1,2   
      kjj=j+kik-1
      nbin=1
      ken1=kend1(kjj,ifil)
      nsp=ken1
      if (tem(81,kjj,ifil).ne.0.d0) nbin=2
      tb=tem(1,kjj,ifil)
      te=tem(ken1,kjj,ifil)
      do 553 iii=1,ken1
      tef(iii)=tem(iii,kjj,ifil)
      y2(iii)=grand(iii,kjj,ifil)
      y(iii)=yf(iii,kjj,ifil)
 553  continue
      if(tt.lt.te.or.nbin.eq.1) goto 56
      ken2=kend2(kjj,ifil)
      nsp=ken2
      tb=tem(81,kjj,ifil)
      te=tem(80+ken2,kjj,ifil)
      do 554 iii=1,ken2
      iik=80+iii
      tef(iii)=tem(iik,kjj,ifil)
      y2(iii)=grand(iik,kjj,ifil)
      y(iii)=yf(iik,kjj,ifil)
 554  continue
 56   continue
      tte(kik)=t   
      if(t.gt.te) tte(kik)=te   
      pha(kik)=(tte(kik)-tb)/(te-tb)   
      call splinterpol(tef,y,y2,nsp,tte(kik),yy(kik))
  61  continue   
      if(g.gt.5.d0) goto 43   
      if(g.lt.0.d0) goto 47   
      slope=(yy(2)-yy(1))*2.d0   
      yy(1)=yy(2)+slope*(g-glog(j+1))   
      slope=(tte(2)-tte(1))*2.d0   
      te=tte(1)+slope*(g-glog(j))   
      thigh=te*(1.d0+fractol)   
      if(t.gt.thigh) goto 79   
      if(ld.lt.0) call limdark(ifil,ldabs,thigh,g,xld,yld)
      call planckint(thigh,ifil,yyhigh,dum)   
      thighmidlog=0.5d0*dlog10(te*thigh)   
      wvlmax=10.d0**(6.4624d0-thighmidlog)   
      if(effwvl(ifil).lt.wvlmax) goto 426   
      thighlog=dlog10(thigh)   
      telog=dlog10(te)   
      slope=(yyhigh-yy(1))/(thighlog-telog)   
      xintlog=yyhigh+slope*(tlog-thighlog)   
      goto 99   
  426 continue   
      thighrec=1.d0/thigh   
      terec=1.d0/te   
      slope=(yyhigh-yy(1))/(thighrec-terec)   
      xintlog=yyhigh+slope*(trec-thighrec)   
      goto 99   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
  43  yy(1)=yy(2)   
      te=tte(2)   
      if(ld.lt.0) call limdark(ifil,ldabs,thigh,g,xld,yld)
      call planckint(thigh,ifil,yyhigh,dum)   
      thighmidlog=0.5d0*dlog10(te*thigh)   
      wvlmax=10.d0**(6.4624d0-thighmidlog)   
      if(effwvl(ifil).lt.wvlmax) goto 427   
      thighlog=dlog10(thigh)   
      telog=dlog10(te)   
      slope=(yyhigh-yy(1))/(thighlog-telog)   
      xintlog=yyhigh+slope*(tlog-thighlog)   
      goto 44   
  427 continue   
      thighrec=1.d0/thigh   
      terec=1.d0/te   
      slope=(yyhigh-yy(1))/(thighrec-terec)   
      xintlog=yyhigh+slope*(trec-thighrec)   
      goto 44   
  47  continue   
      te=tte(1)   
      if(ld.lt.0) call limdark(ifil,ldabs,thigh,g,xld,yld)
      call planckint(thigh,ifil,yyhigh,dum)   
      thighmidlog=0.5d0*dlog10(te*thigh)   
      wvlmax=10.d0**(6.4624d0-thighmidlog)   
      if(effwvl(ifil).lt.wvlmax) goto 428   
      thighlog=dlog10(thigh)   
      telog=dlog10(te)   
      slope=(yyhigh-yy(1))/(thighlog-telog)   
      xintlog=yyhigh+slope*(tlog-thighlog)   
      goto 63   
  428 continue   
      thighrec=1.d0/thigh   
      terec=1.d0/te   
      slope=(yyhigh-yy(1))/(thighrec-terec)   
      xintlog=yyhigh+slope*(trec-thighrec)   
      goto 63   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
   7  continue   
      thigh=6000.d0*(1.d0+fractol)   
      if(t.gt.thigh) goto 79   
      if(pha(1).le.1.d0) goto 63   
      goto 5   
  63  continue   
      if(ld.lt.0) call limdark(ifil,ldabs,t,g,xld,yld)
      call planckint(t,ifil,yylow,dum)   
      slope=(yylow-xintlog)/glow   
      xintlog=yylow+slope*(g-glow)   
      goto 99   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
   9  continue   
      thigh=50000.d0*(1.d0+fractol)   
      if(t.gt.thigh) goto 79   
      if(t.gt.50000.d0) goto 52   
  44  continue   
      if(ld.lt.0) call limdark(ifil,ldabs,t,g,xld,yld)
      call planckint(t,ifil,yyhigh,dum)   
      slope=(yyhigh-xintlog)/(ghigh-5.d0)   
      xintlog=yyhigh+slope*(g-ghigh)   
      goto 99   
  52  continue   
      j=10   
      goto 5   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
  66  continue   
      message(komp,4)=1   
      if(ld.lt.0) call limdark(ifil,ldabs,t,g,xld,yld)
      call planckint(t,ifil,xintlog,xint)   
      return   
  77  continue   
      message(komp,1)=1   
      if(ld.lt.0) call limdark(ifil,ldabs,t,g,xld,yld)
      call planckint(t,ifil,xintlog,xint)   
      return   
  78  continue   
      message(komp,2)=1   
      if(ld.lt.0) call limdark(ifil,ldabs,t,g,xld,yld)
      call planckint(t,ifil,xintlog,xint)   
      return   
  79  continue   
      message(komp,3)=1   
      if(ld.lt.0) call limdark(ifil,ldabs,t,g,xld,yld)
      call planckint(t,ifil,xintlog,xint)   
      return   
  99  continue  
      xint=10.d0**xintlog   
      return  
      end   
      subroutine gabs(komp,smaxis,qq,ecc,period,dd,rad,xm,xmo,absgr, 
     $glog) 
      implicit real*8(a-h,o-z) 
c  Version of August 6, 2009 
c 
c  Input definitions: 
c   smaxis is the length of the orbital semi-major axis in solar radii. 
c   qq is the mass ratio in the sense m2/m1. Stars 1 and 2 are as defined 
c     in the external program (star 1 is near superior conjunction at 
c     phase zero). 
c   ecc is orbital eccentricity 
c   period is orbit period in days 
c   dd is the instantaneous separation of the star centers in unit of th 
c     orbital semi-major axis 
c   rad is the polar radius if the star at issue in unit of the orbital  
c     semi-major axis 
c  Output definitions: 
c   absgr is the polar acceleration due to effective gravity in cm^2/sec 
c   glog is log_10 of absgr 
c 
      twopi=6.2831853072d0 
      gbig=6.674215d-8 
      sunmas=1.988435d33 
      sunrad=6.9600d10 
      psec=8.64d4*period 
      acm=sunrad*smaxis 
      pyears=period/365.2422d0 
      aau=smaxis/214.93947d0 
      tmass=aau**3/pyears**2 
      qf=1.d0/(1.d0+qq) 
      qfm=qq*qf 
      sign=-1.d0 
      if(komp.eq.2) goto 10 
      qfm=qf 
      qf=qq*qf 
      sign=1.d0 
   10 continue 
      xm=tmass*qfm 
      xmo=tmass*qf 
      gbigm=gbig*xm*sunmas 
      gbigmo=gbig*xmo*sunmas 
      rcm=rad*acm 
      dcm=dd*acm 
      dcmsq=dcm*dcm 
      efac=dsqrt((1.d0+ecc)*(1.d0-ecc)) 
      av=twopi*efac/(psec*dd*dd) 
      avsq=av*av 
      rcmsq=rcm*rcm 
      hypsq=rcmsq+dcmsq 
      hyp=dsqrt(hypsq) 
      snalf=rcm/hyp 
      csalf=dcm/hyp 
      gz=-gbigm/rcmsq 
      gzo=-snalf*gbigmo/hypsq 
      gxo=sign*csalf*gbigmo/hypsq 
      gxcf=-sign*avsq*dcm*qf 
      gxs=gxo+gxcf 
      gzs=gz+gzo 
      absgr=dsqrt(gxs*gxs+gzs*gzs) 
      glog=dlog10(absgr) 
      return 
      end 
      SUBROUTINE SURFAS(RMASS,POTENT,N,N1,KOMP,RV,GRX,GRY,GRZ,RVQ,       
     $GRXQ,GRYQ,GRZQ,MMSAVE,FR1,FR2,HLD,FF,D,SNTH,CSTH,SNFI,CSFI,GRV1,   
     $GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2,GMAG1,      
     $GMAG2,glog1,glog2,GREXP)                                                       
c  Version of November 22, 2011                                            
      implicit real*8 (a-h,o-z)                                          
      DIMENSION RV(*),GRX(*),GRY(*),GRZ(*),RVQ(*),GRXQ(*),GRYQ(*),GRZQ(* 
     $),MMSAVE(*),FR1(*),FR2(*),HLD(*),SNTH(*),CSTH(*),SNFI(*),CSFI(*)   
     $,GRV1(*),GRV2(*),XX1(*),YY1(*),ZZ1(*),XX2(*),YY2(*),ZZ2(*),GLUMP1  
     $(*),GLUMP2(*),CSBT1(*),CSBT2(*),GMAG1(*),GMAG2(*),glog1(*),                 
     $glog2(*) 
      common /gpoles/ gplog1,gplog2 
      common /radi/ R1H,RLH,R1C,RLC                                      
      COMMON /misc/ X1                                                        
      COMMON /ECCEN/e,smaxis,period,vgadum,sindum,vfdum,vfadum,vgmdum,     
     $v1dum,v2dum,ifcdum                                                 
      pi=dacos(-1.d0)
      DSQ=D*D                                                            
      RMAS=RMASS                                                         
      IF(KOMP.EQ.2) RMAS=1.d0/RMASS                                      
      RF=FF**2                                                           
      RTEST=0.d0                                                         
      IP=(KOMP-1)*(N1+1)+1                                               
      IQ=IP-1                                                            
      IS=0                                                               
      isx=0
      if(iq.gt.0) ISX=(KOMP-1)*MMSAVE(IQ)                                            
      MMSAVE(IP)=0                                                       
      KFLAG=0                                                            
      CALL ELLONE (FF,D,RMAS,X1,OMEGA,XL2,OM2)                           
      IF(KOMP.EQ.2) OMEGA=RMASS*OMEGA+.5d0*(1.d0-RMASS)                  
      X2=X1                                                              
      IF(KOMP.EQ.2) X1=1.d0-X1                                           
      IF(E.NE.0.d0) GOTO 43                                              
      IF(POTENT.LT.OMEGA) CALL NEKMIN(RMASS,POTENT,X1,ZZ)                
      IF(POTENT.LT.OMEGA) X2=1.d0-X1                                     
   43 COMP=dfloat(3-2*KOMP)                                              
      CMP=dfloat(KOMP-1)                                                 
      CMPD=CMP*D                                                         
      TESTER=CMPD+COMP*X1
      RM1=RMASS+1.d0                                                     
      RMS=RMASS                                                          
      RM1S=RM1                                                           
      IF(KOMP.NE.2) GOTO 15                                              
      POT=POTENT/RMASS+.5d0*(RMASS-1.d0)/RMASS                           
      RM=1.d0/RMASS                                                      
      RM1=RM+1.d0                                                        
      GOTO 20                                                            
   15 POT=POTENT                                                         
      RM=RMASS                                                           
   20 EN=N                                                               
c  Find the relative polar radius, R/a 
      DELR=0.d0                                                          
      R=1.d0/pot                                                            
      knt=0
  714 R=R+DELR                                                           
      knt=knt+1
      tolr=1.d-6*dabs(r)
      RSQ=R*R                                                            
      PAR=DSQ+RSQ                                          
      RPAR=dsqrt(PAR)                                                    
      OM=1.d0/R+RM/RPAR       
      DOMR=1.d0/(-1.d0/RSQ-RM*R/(PAR*RPAR))     
      DELR=(POT-OM)*DOMR                                                 
      ABDELR=dabs(DELR)                                                  
      IF(ABDELR.GT.tolr) GOTO 714                                     
      rpole=r 
      rsave=r
c  Now compute GRPOLE (exactly at the pole) 
      x=cmpd 
      zsq=rpole*rpole 
      PAR1=x*x+zsq                                                    
      RPAR1=dsqrt(PAR1)                                                  
      XNUM1=1.d0/(PAR1*RPAR1)                                            
      XL=D-X                                                             
      PAR2=XL**2+zsq                                                   
      RPAR2=dsqrt(PAR2)                                                  
      XNUM2=1.d0/(PAR2*RPAR2)                                            
      OMZ=-rpole*(XNUM1+RMS*XNUM2)                                           
      OMX=RMS*XL*XNUM2-X*XNUM1+RM1S*X*RF-RMS/DSQ                         
      IF(KOMP.EQ.2) OMX=RMS*XL*XNUM2-X*XNUM1-RM1S*XL*RF+1.d0/DSQ         
      grpole=dsqrt(OMX*OMX+OMZ*OMZ)                               
      call gabs(komp,smaxis,rmass,e,period,d,rpole,xmas,xmaso,absgr, 
     $glogg) 
      if(komp.eq.1) gplog1=glogg 
      if(komp.eq.2) gplog2=glogg 
      oldmu=0.d0
      par32=0.d0
      DO 8 I=1,N                                                         
      IF(I.NE.2) GOTO 82                                                 
      IF(KOMP.EQ.1) RTEST=.3d0*RV(1)                                     
      IF(KOMP.EQ.2) RTEST=.3d0*RVQ(1)                                    
   82 CONTINUE                                                           
      IPN1=I+N1*(KOMP-1)                                                 
      SINTH=SNTH(IPN1)                                                   
      XNU=CSTH(IPN1)                                                     
      XNUSQ=XNU**2                                                       
      EM=SINTH*EN*1.3d0                                                  
      XLUMP=1.d0-XNUSQ                                                   
      MM=EM+1.d0                                                         
      afac=rf*rm1*xlump 
      DO 8 J=1,MM                                                        
      KOUNT=0                                                            
      IS=IS+1                                                            
      ISX=ISX+1                                                          
      DELR=0.d0                                                          
      COSFI=CSFI(ISX)                                                    
      XMU=SNFI(ISX)*SINTH                                                
      XLAM=SINTH*COSFI                                                   
      bfac=xlam*d 
      efac=rm*xlam/dsq 
      R=RSAVE                                                            
      oldr=r
      knth=0
c
c   Compute starting estimate for r based on r and d(omega)/d(phi)
c     times d(omega)/dr at previous point (only for J = 2 thru 6)
c
      if(j.eq.1.or.j.gt.6) goto 717
      XM=dfloat(MM)
      delfi=pi/xm
      dodfi=-rm*oldmu*((r*d/par32-r/dsq))
      delr=-dodfi*delfi*domr
  717 continue
   14 R=R+DELR                                                           
      tolr=1.d-6*dabs(r)
      if(kount.lt.1) goto 170
      if(knth.gt.20) goto 170
      if(r.gt.0.d0.and.r.lt.tester) goto 170
      knth=knth+1
      delr=0.5d0*delr
      r=oldr
      goto 14
  170 continue
      KOUNT=KOUNT+1                                                      
      IF(KOUNT.LT.80) GOTO 70                                            
      KFLAG=1                                                            
      R=-1.d0                                                            
      GOTO 86                                                            
   70 continue
      RSQ=R*R                                                            
      rcube=r*rsq 
      PAR=DSQ-2.d0*XLAM*R*D+RSQ                                          
      RPAR=dsqrt(PAR)                                                    
      par32=par*rpar 
      par52=par*par32 
      OM=1.d0/R+RM*((1.d0/RPAR)-XLAM*R/DSQ)+RM1*.5d0*RSQ*XLUMP*RF        
      denom=RF*RM1*XLUMP*R-1.d0/RSQ-(RM*(R-XLAM*D))/par32-efac      
      domr=1.d0/denom 
      d2rdo2=-domr*(afac+2.d0/rcube-rm*(1.d0/par32-3.d0*(r-bfac)**2/
     $par52))/denom**2 
      DELR=(POT-OM)*DOMR+.5d0*(pot-om)**2*d2rdo2
      oldr=r
      ABDELR=dabs(DELR)                                                  
      IF(ABDELR.GT.tolr) GOTO 14                                     
      ABR=dabs(R)                                                        
      IF(R.GT.RTEST) GOTO 74                                             
      KFLAG=1                                                            
      R=-1.d0                                                            
      IF(KOMP.EQ.2) GOTO 98                                              
      GOTO 97                                                            
   74 IF(ABR.LT.TESTER) RSAVE=R                                          
      Z=R*XNU                                                            
      Y=COMP*R*XMU                                                       
      X2T=ABR*XLAM                                                       
      X=CMPD+COMP*X2T                                                    
      IF(KOMP.EQ.2) GOTO 62                                              
      IF(X.LT.X1) GOTO 65                                                
      KFLAG=1                                                            
      R=-1.d0                                                            
      GOTO 97                                                            
   62 IF(X2T.LT.X2) GOTO 65                                              
      KFLAG=1                                                            
      R=-1.d0                                                            
      GOTO 98                                                            
   65 SUMSQ=Y**2+Z**2                                                    
      PAR1=X**2+SUMSQ                                                    
      RPAR1=dsqrt(PAR1)                                                  
      XNUM1=1.d0/(PAR1*RPAR1)                                            
      XL=D-X                                                             
      PAR2=XL**2+SUMSQ                                                   
      RPAR2=dsqrt(PAR2)                                                  
      XNUM2=1.d0/(PAR2*RPAR2)                                            
      OMZ=-Z*(XNUM1+RMS*XNUM2)                                           
      OMY=Y*(RM1S*RF-XNUM1-RMS*XNUM2)                                    
      OMX=RMS*XL*XNUM2-X*XNUM1+RM1S*X*RF-RMS/DSQ                         
      IF(KOMP.EQ.2) OMX=RMS*XL*XNUM2-X*XNUM1-RM1S*XL*RF+1.d0/DSQ         
      GRMAG=dsqrt(OMX*OMX+OMY*OMY+OMZ*OMZ)                               
      grvrat=grmag/grpole 
      GRAV=grvrat**GREXP                                         
      A=COMP*XLAM*OMX                                                    
      B=COMP*XMU*OMY                                                     
      C=XNU*OMZ                                                          
      COSBET=-(A+B+C)/GRMAG                                              
      IF(COSBET.LT..7d0) COSBET=.7d0                                     
   86 IF(KOMP.EQ.2) GOTO 98                                              
   97 RV(IS)=R                                                           
      GRX(IS)=OMX                                                        
      GRY(IS)=OMY                                                        
      GRZ(IS)=OMZ                                                        
      GMAG1(IS)=dsqrt(OMX*OMX+OMY*OMY+OMZ*OMZ)                           
      glog1(is)=dlog10(grvrat*absgr) 
      FR1(IS)=1.d0                                                       
      GLUMP1(IS)=R*R*SINTH/COSBET                                        
      GRV1(IS)=GRAV                                                      
      XX1(IS)=X                                                          
      YY1(IS)=Y                                                          
      ZZ1(IS)=Z                                                          
      CSBT1(IS)=COSBET                                                   
      GOTO 8                                                             
   98 RVQ(IS)=R                                                          
      GRXQ(IS)=OMX                                                       
      GRYQ(IS)=OMY                                                       
      GRZQ(IS)=OMZ                                                       
      GMAG2(IS)=dsqrt(OMX*OMX+OMY*OMY+OMZ*OMZ)                           
      glog2(is)=dlog10(grvrat*absgr) 
      FR2(IS)=1.d0                                                       
      GLUMP2(IS)=R*R*SINTH/COSBET                                        
      GRV2(IS)=GRAV                                                      
      XX2(IS)=X                                                          
      YY2(IS)=Y                                                          
      ZZ2(IS)=Z                                                          
      CSBT2(IS)=COSBET                                                   
      oldmu=xmu
    8 CONTINUE                                                           
      if(e.ne.0.d0.or.ff.ne.1.d0) goto 53
      IF(KFLAG.EQ.0) GOTO 53                                             
      ISS=IS-1                                                           
      IF(KOMP.NE.1) GOTO 50                                              
      CALL RING(RMASS,POTENT,1,N,FR1,HLD,R1H,RLH)                        
      DO 55 I=1,ISS                                                      
      IPL=I+1                                                            
      IF(RV(I).GE.0.d0)GOTO 55                                           
      FR1(IPL)=FR1(IPL)+FR1(I)                                           
      FR1(I)=0.d0                                                        
   55 CONTINUE                                                           
   53 IF(KOMP.EQ.2) GOTO 54                                              
      IS=0                                                               
      DO 208 I=1,N                                                       
      IPN1=I+N1*(KOMP-1)                                                 
      EM=SNTH(IPN1)*EN*1.3d0                                             
      MM=EM+1.d0                                                         
      DO 208 J=1,MM                                                      
      IS=IS+1                                                            
      GLUMP1(IS)=FR1(IS)*GLUMP1(IS)                                      
  208 CONTINUE                                                           
      RETURN                                                             
   50 if(e.ne.0.d0.or.ff.ne.1.d0) goto 54
      CALL RING(RMASS,POTENT,2,N,FR2,HLD,R1C,RLC)                        
      DO 56 I=1,IS                                                       
      IPL=I+1                                                            
      IF(RVQ(I).GE.0.d0) GOTO 56                                         
      FR2(IPL)=FR2(IPL)+FR2(I)                                           
      FR2(I)=0.d0                                                        
   56 CONTINUE                                                           
   54 CONTINUE                                                           
      IS=0                                                               
      DO 108 I=1,N                                                       
      IPN1=I+N1*(KOMP-1)                                                 
      EM=SNTH(IPN1)*EN*1.3d0                                             
      MM=EM+1.d0                                                         
      DO 108 J=1,MM                                                      
      IS=IS+1                                                            
      GLUMP2(IS)=FR2(IS)*GLUMP2(IS)                                      
  108 CONTINUE                                                           
      RETURN                                                             
      END                                                                
      subroutine fourls(th,ro,obs,nobs,nth,aa,bb)                            
      implicit real*8(a-h,o-z)                                           
c
c   version of June 2, 2011                                        
c                                                                        
c    Input integer nth is the largest Fourier term fitted (e.g.          
c       for nth=6, terms up to sine & cosine of 6 theta are              
c       evaluated).                                                      
c    This subroutine can handle nth only up to 6. Additional             
c      programming is needed for larger values.                          
c                                                                        
      dimension aa(*),bb(*),th(*),ro(*),obs(*),ll(14),mm(14),         
     $cn(196),cl(14),out(14)                                             
      mpl=nth+1                                                          
      ml=mpl+nth                                                         
      jjmax=ml*ml                                                        
      nobsml=nobs*ml                                                     
      nobmpl=nobs*mpl                                                    
      do 90 i=1,nobs                                                     
      obs(i)=1.d0                                                        
      iz=nobsml+i                                                        
      obs(iz)=ro(i)                                                      
      if(nth.eq.0) goto 90                                               
      ic=i+nobs                                                          
      is=i+nobmpl                                                        
      sint=dsin(th(i))                                                   
      cost=dcos(th(i))                                                   
      obs(ic)=cost                                                       
      obs(is)=sint                                                       
      if(nth.eq.1) goto 90                                               
      ic=ic+nobs                                                         
      is=is+nobs                                                         
      sncs=sint*cost                                                     
      cs2=cost*cost                                                      
      obs(ic)=cs2+cs2-1.d0                                               
      obs(is)=sncs+sncs                                                  
      if(nth.eq.2) goto 90                                               
      ic=ic+nobs                                                         
      is=is+nobs                                                         
      sn3=sint*sint*sint                                                 
      cs3=cs2*cost                                                       
      obs(ic)=4.d0*cs3-3.d0*cost                                         
      obs(is)=3.d0*sint-4.d0*sn3                                         
      if(nth.eq.3) goto 90                                               
      ic=ic+nobs                                                         
      is=is+nobs                                                         
      cs4=cs2*cs2                                                        
      obs(ic)=8.d0*(cs4-cs2)+1.d0                                        
      obs(is)=4.d0*(2.d0*cs3*sint-sncs)                                  
      if(nth.eq.4) goto 90                                               
      ic=ic+nobs                                                         
      is=is+nobs                                                         
      cs5=cs3*cs2                                                        
      obs(ic)=16.d0*cs5-20.d0*cs3+5.d0*cost                              
      obs(is)=16.d0*sn3*sint*sint-20.d0*sn3+5.d0*sint                    
      if(nth.eq.5) goto 90                                               
      ic=ic+nobs                                                         
      is=is+nobs                                                         
      obs(ic)=32.d0*cs3*cs3-48.d0*cs4+18.d0*cs2-1.d0                     
      obs(is)=32.d0*sint*(cs5-cs3)+6.d0*sncs                             
   90 continue                                                           
      do 20 jj=1,jjmax                                                   
   20 cn(jj)=0.d0                                                        
      do 21 j=1,ml                                                       
   21 cl(j)=0.d0                                                         
      do 24 nob=1,nobs                                                   
      iii=nob+nobsml                                                     
      do 23 k=1,ml                                                       
      do 23 i=1,ml                                                       
      ii=nob+nobs*(i-1)                                                  
      kk=nob+nobs*(k-1)                                                  
      j=i+(k-1)*ml                                                       
   23 cn(j)=cn(j)+obs(ii)*obs(kk)                                        
      do 24 i=1,ml                                                       
      ii=nob+nobs*(i-1)                                                  
   24 cl(i)=cl(i)+obs(iii)*obs(ii)                                       
      call dminv(cn,ml,d,ll,mm)                                          
      call dgmprd(cn,cl,out,ml,ml,1)                                     
      do 51 i=2,mpl                                                      
      aa(i)=out(i)                                                       
      ipl=i+nth                                                          
   51 bb(i)=out(ipl)                                                     
      aa(1)=out(1)                                                       
      bb(1)=0.d0                                                         
      return                                                             
      end                                                                
      subroutine frspot(komp,kspot,nmi,N,sinth,costh,sinfi,cosfi,delth,
     $delfi,temf)
c                                                                        
c   Subroutine FRSPOT computes an overall temperature factor (TEMF, ratio
c      of temperature to spot-free temperature) for a given surface element
c      that may be affected by up to N spots. If an element is in more
c      than one spot, FRSPOT adopts the product of the spot temperature
c      factors. An element that is partly covered by a spot is assigned 
c      a fractional covered area that is computed by spherical geometry,
c      based on overlap between the element and the circular spot. The
c      effective temperature factor for that element-spot combo is a
c      mean of the out-of-spot factor (unity) and the in-spot factor,
c      weighted by fractional covered area.
c                                                                        
c   "Latitudes" here run from 0 at the +z pole to pi radians         
c      at the other.                                                     
C                                                                        
c   Version of October 4, 2011                                         
c                                                                        
      implicit real*8 (a-h,o-z)                                          
      parameter (ispmax=   100)
      dimension cslat(4),snlat(4),cslon(4),snlon(4)
      common /inprof/ in1min,in1max,in2min,in2max,mpage,nl1,nl2          
      COMMON /SPOTS/ SINLAT(2,ispmax),COSLAT(2,ispmax),SINLNG(2,ispmax),
     $COSLNG(2,ispmax),rdsp(2,ispmax),temsp(2,ispmax),xlng(2,ispmax),
     $kks(2,ispmax),Lspot(2,ispmax)
      common /spev/ tstart(2,ispmax),tmaxa(2,ispmax),tmaxb(2,ispmax),
     $tfinal(2,ispmax),amax(2,ispmax),sprad(2,ispmax),ndum
      kvec=3
      pi=dacos(-1.d0)
      hdelth=0.5d0*delth
      hdelfi=0.5d0*delfi
      tolsprad=1.d-2
c
c  If kspot=1, apply simple logic based on whether the center of the surface
c    element is inside the spot.
c  If kspot=2, apply fractional area logic based on the part of the element
c    that's inside the spot.
c
c  spotmax is half of the maximum effective angular size of a surface element
c    at this latitude. It's for a fast test to see if the element is clearly
c    outside a given spot.
c
      spotmax=1.6d0*sinth*hdelfi
      TEMF=1.d0                                                          
      as1=0.d0
      at1=0.d0
      as2=0.d0
      at2=0.d0
      nl=(2-komp)*nl1+(komp-1)*nl2                                       
      DO 15 I=1,N                                                        
      if(sprad(komp,i).lt.tolsprad) goto 15
      COSDFI=COSFI*COSLNG(KOMP,I)+SINFI*SINLNG(KOMP,I)
      arc=dacos(costh*coslat(komp,i)+sinth*sinlat(komp,i)*cosdfi)
      if(kspot.eq.2) goto 43
      IF(arc.gt.sprad(komp,i)) goto 15
      TEMF=TEMF*TEMSP(KOMP,I)
      goto 44
   43 dcr=sprad(komp,i)+spotmax
      if(arc.gt.dcr) GOTO 15
      do 42 j=1,nl                                                       
   42 if(kks(komp,j).eq.-i) Lspot(komp,j)=Lspot(komp,j)+1                
c
c  Copy coordinate information for vertices 1, 2, and 3 of the local
c    surface element into TRICIR input arrays (to form a triangular
c    half-element). Surface element vertices are labeled counter-clockwise
c    (seen from outside the unit sphere). Picture of the 4 vertices is below:
c
c            1    <--   4
c            |          |
c            2    -->   3
c
      theta=dacos(costh)
      phi=dacos(cosfi)
      if(sinfi.lt.0.d0) phi=2.d0*pi-phi
      cslat(1)=dcos(theta-hdelth)
      snlat(1)=dsin(theta-hdelth)
      cslat(2)=dcos(theta+hdelth)
      snlat(2)=dsin(theta+hdelth)
      cslat(3)=cslat(2)
      snlat(3)=snlat(2)
      cslon(1)=dcos(phi-hdelfi)
      snlon(1)=dsin(phi-hdelfi)
      cslon(2)=cslon(1)
      snlon(2)=snlon(1)
      cslon(3)=dcos(phi+hdelfi)
      snlon(3)=dsin(phi+hdelfi)
c
c  Copy spot center coordinate information into TRICIR input arrays.
c
      cslat(4)=coslat(komp,i)
      snlat(4)=sinlat(komp,i)
      cslon(4)=coslng(komp,i)
      snlon(4)=sinlng(komp,i)
      radsp=sprad(komp,i)
c
c   call TRICIR for a triangular "half" of the local surface element
c
      call tricir(cslat,snlat,cslon,snlon,radsp,pi,kvec,as1,at1)
      if(nmi.eq.0) goto 119
c
c  Replace vertex 2, entered above, with vertex 4 of the surface element
c   to make the other triangular half-element.
c
      cslat(2)=cslat(1)
      snlat(2)=snlat(1)
      cslon(2)=cslon(3)
      snlon(2)=snlon(3)
c
c   call TRICIR for the other half of the element.
c
      call tricir(cslat,snlat,cslon,snlon,radsp,pi,kvec,as2,at2)
  119 continue
      denom=at1+at2
      if(denom.gt.0.d0) fract=(as1+as2)/denom
c
c  The next two lines prevent fract from exceeding its logical range of
c    zero (none of the surface element is outside the spot) to unity
c    (all of the element is outside the spot).
c
      if(fract.gt.1.d0) fract=1.d0
      if(fract.lt.0.d0) fract=0.d0
      temf=temf*(fract+(1.d0-fract)*temsp(komp,i))
   44 if(mpage.ne.3) goto 15                                             
      do 24 j=1,nl                                                       
      kk=kks(komp,j)                                                     
      if(kk.eq.-i) Lspot(komp,j)=0                                       
      if(kk.eq.i) Lspot(komp,j)=Lspot(komp,j)+1                          
   24 continue                                                           
   15 continue                                                           
      RETURN                                                             
      END                                                                
      subroutine tricir(cslat,snlat,cslon,snlon,rad,pi,kvec,aclear,atri)
c
c   May 18, 2011
c
      implicit real*8(a-h,o-z)
c
c   The name TRICIR comes from: TRIangle-CIRcle intersection.
c
c   Subroutine TRICIR computes the intersection and non-intersection areas of a spherical 
c     triangle that intersects a circle on the unit sphere. Inputs are the
c     triangle's 3 sets of vertex coordinates (latitudes and longitudes) and the circle's
c     center coordinates and radius (so 4 vectors). The circle-center vector must
c     be vector 4. All arcs and coordinates are in radians. Areas are in steradians, as subtended at the
c     center of the sphere, rad is the circle's angular radius and must be in the range
c     0 to pi radians. 
c        All of this is done on a sphere of unit radius.
c
      dimension vec1(3),vec2(3),vec3(3),vec4(3),cslat(*),snlat(*),
     $cslon(*),snlon(*),side1(300),side2(300),side3(300),
     $arc1(99),arc2(99),arc3(99),crossdum(3),veco(3),veco1(3),
     $veco2(3),veco3(3),vside(3),s1fr(2),s2fr(2),s3fr(2)
      twopi=pi+pi
c
c   Operation:
c     TRICIR is called for ONE triangle at a time (call TRICIR more 
c   than once for a polygon of more than three sides):
c      First represent each triangle side by k+2 points, one at each vertex
c   and k equally spaced ones in between (k is input quantity kvec). So we have 3*(1+k)
c   vectors that represent the sides, with three of the vectors for the 
c   vertices and the rest spaced between. 
c   Input arrays cslat, snlat, cslon, and snlon are cosines and sines of the (vector)
c   points. In most realistic situations for stars, the surface elements are
c   much smaller than the (circular) spots, so TRICIR neglects the case where
c   a (triangular) element extends through a spot, with all three vertices outside
c   the circle and two sides crossing the circle. 
c      How to interpolate between two vectors to make vectors equally spaced
c   in angle, in the plane of the original two vectors? The cross product
c   of the original vectors is A*B sin(theta), theta being the angle between 
c   them. Since the original vectors have unit length in the immediate problem,
c   the cross product length will be sin(theta), and the angle increments will be
c   uniformly stepped from [1/(k+1)](theta) to [k/(k+1)](theta).
c   So how to get rectangular components of these new vectors?
c   We know the vectors' plane (via the cross product's direction)
c   and we know the angles, as computed just above, so the components can be
c   computed by rotating vector 1 in equally spaced steps about the cross product. 
c   Subroutine VECIN does all that.
c   
c   Initializations:
c
      rangemin=0.d0
      tolarc=1.d-5
      frinc=1.d-8
      do 91 iarc=1,99
      arc1(iarc)=0.d0
      arc2(iarc)=0.d0
   91 arc3(iarc)=0.d0
      kvec1=kvec+1
      kvec2=kvec+2
      cvec1=1.d0/dfloat(kvec1)
      cvec2=1.d0/dfloat(kvec2)
      knt=0
c      
c    Defaults:
c
      kv1=1
      kv2=1
      kv3=1
      ks1=1
      ks2=1
      ks3=1
c
c  rectangular components of the vectors that represent the triangle vertices and
c    circle center are:
c
      vec1(1)=cslon(1)*snlat(1)
      vec1(2)=snlon(1)*snlat(1)
      vec1(3)=cslat(1)
      vec2(1)=cslon(2)*snlat(2)
      vec2(2)=snlon(2)*snlat(2)
      vec2(3)=cslat(2)
      vec3(1)=cslon(3)*snlat(3)
      vec3(2)=snlon(3)*snlat(3)
      vec3(3)=cslat(3)
      vec4(1)=cslon(4)*snlat(4)
      vec4(2)=snlon(4)*snlat(4)
      vec4(3)=cslat(4)
c
c  Compute triangle area ('atri') from vertex vectors. 'aclear' is the part of the 
c    triangle area that's not within the circle (it's "in the clear").
c
      call arcver(vec1,vec2,vec3,arcv1,arcv2,arcv3,ver1,ver2,ver3,3,
     $iflag)
      vertotal=dabs(ver1)+dabs(ver2)+dabs(ver3)
      atri=vertotal-pi  
c
c  Default for aclear 
c
      aclear=atri
c
c  Now check whether each triangle vertex is within the circle.
c
c  Compute arc from the circle center to vertex 1:
      call arcver(vec4,vec1,vec3,dum,dum1,arcv1,dum3,dum4,dum5,1,iflag)
      if(arcv1.lt.rad) kv1=0
c  Compute arc from the circle center to vertex 2:
      call arcver(vec4,vec2,vec3,dum,dum1,arcv2,dum3,dum4,dum5,1,iflag)
      if(arcv2.lt.rad) kv2=0
c  Compute arc from the circle center to vertex 3:
      call arcver(vec4,vec3,vec2,dum,dum1,arcv3,dum3,dum4,dum5,1,iflag)
      if(arcv3.lt.rad) kv3=0
c
      kvvert=kv1+kv2+kv3
c
c  Case: Triangle entirely within the circle. Set aclear=0 and return.
c   Otherwise go on to detailed checks.
c
      if(kvvert.ne.0) goto 99
      aclear=0.d0
      return
   99 continue
c
c  Interpolate 'inside' point vectors for sides 1, 2, and 3 in calls to VECIN.
c    Output point vectors SIDE1, SIDE2, and SIDE3 have the rectangular components
c    of the triangle vertices as their first three and last three elements, with
c    the interpolated vector components in their middle elements. Because the
c    interpolated side vectors run from an "odd man" vertex (outside the circle,
c    with the other two vertices inside, or the reverse situation) toward each of 
c    the other two vertices, the VECIN calls are either done in pairs or not at all.
c
      call vecin(vec1,vec2,crossdum,side3,veco,0.5d0,kvec,1)
      call vecin(vec3,vec1,crossdum,side2,veco,0.5d0,kvec,1)
      call vecin(vec2,vec3,crossdum,side1,veco,0.5d0,kvec,1)
c
c  Find which case applies.
c  The following 6 defaults are set to small negative flags so as 
c   to be easily identified as flags (by being negative) and,
c   (being small) not be likely to produce significant errors 
c   in any rare cases that may arise.
c
      s1fr(1)=-9.d-5
      s1fr(2)=-9.d-5
      s2fr(1)=-9.d-5
      s2fr(2)=-9.d-5
      s3fr(1)=-9.d-5
      s3fr(2)=-9.d-5
c
c  Compute arc from the circle center to each point of 
c    side 1 and decide if the side is entirely beyond the circle (ks1=1).
c
      do 35 iv=1,kvec2
      do 31 iz=1,3
   31 vside(iz)=side1(3*(iv-1)+iz)
      call arcver(vec4,vside,vec3,dum,dum1,arc1(iv),dum3,dum4,dum5,1,
     $iflag)
      if(arc1(iv).ge.rad) goto 35
      ks1=0
   35 continue
c
c  Compute arc from the circle center to each point of 
c    side 2 and decide if the side is entirely beyond the circle (ks2=1).
c
      do 45 iv=1,kvec2
      do 32 iz=1,3
   32 vside(iz)=side2(3*(iv-1)+iz)
      call arcver(vec4,vside,vec3,dum,dum1,arc2(iv),dum3,dum4,dum5,1,
     $iflag)
      if(arc2(iv).ge.rad) goto 45
      ks2=0
   45 continue
c
c  Compute arc from the circle center to each point of 
c    side 3 and decide if the side is entirely beyond the circle (ks3=1).
c
      do 55 iv=1,kvec2
      do 33 iz=1,3
   33 vside(iz)=side3(3*(iv-1)+iz)
      call arcver(vec4,vside,vec2,dum,dum1,arc3(iv),dum3,dum4,dum5,1,
     $iflag)
      if(arc3(iv).ge.rad) goto 55
      ks3=0
   55 continue
      kside=ks1+ks2+ks3
c
c  Compute the three vertex angles for each of three sub-triangles made of the circle center
c    and pairs of vertices in the element triangle, to allow computation
c    of the sub-triangle areas. If subroutine ARCVER sends back IFLAG=1, thus
c    telling of a problem with at least one of the vertex angles, then keep
c    the defaults and return.
c
      call arcver(vec4,vec2,vec3,dum,dum1,dum2,ver1,ver2,ver3,3,iflag)
      if(iflag.eq.1) return
      asub1=dabs(ver1)+dabs(ver2)+dabs(ver3)-pi
      call arcver(vec4,vec3,vec1,dum,dum1,dum2,ver1,ver2,ver3,3,iflag)
      if(iflag.eq.1) return
      asub2=dabs(ver1)+dabs(ver2)+dabs(ver3)-pi
      call arcver(vec4,vec1,vec2,dum,dum1,dum2,ver1,ver2,ver3,3,iflag)
      if(iflag.eq.1) return
      asub3=dabs(ver1)+dabs(ver2)+dabs(ver3)-pi
      areasum=asub1+asub2+asub3
c
c  At this point TRICIR has determined which, if any, of the 3 triangle vertices
c    and the 3 triangle sides lie outside the circle, according to the values (0 or 1)
c    of kv1, kv2, & kv3 for vertices, ks1, ks2, & ks3 for sides, and whether the
c    angles subtended at the circle-center by the three vertex pairs add to more than 6.0d0.
c    TRICIR also remembers the angular distances of all vertices
c    and all the 3k side points from the circle center. 
c  Next is to distinguish between cases where the circle center is inside or outside
c    of the triangle. Do that according to whether areasum (combined area of three
c    sub triangles) exceeds the main triangle area.
c
c  Case: All vertices and sides outside the circle, areasum > atri.
c    (So circle entirely outside the triangle. Keep the defaults.)
c    This will usually be the most common case, so it is tested first.
c
      if(kvvert.eq.3.and.kside.eq.3.and.areasum.ge.atri) return
c
c  Case: All vertices and sides outside the circle, areasum.le.atri.
c    (So circle entirely within the triangle.)
c
      if(kvvert.ne.3.or.kside.ne.3.or.areasum.ge.atri) goto 103
c
c  Compute circle area
c
      acir=twopi*(1.d0-dcos(rad))
      aclear=atri-acir
      return
c
c  Remaining cases involve intersection of at least one
c   side with the circle, so these cases are not treated separately - just check 
c   sides 1, 2, and 3.
c
c  Check side 1:
c
  103 if(ks1.eq.1) goto 104
c
      kside=1
      oldif=arc1(1)-rad
c
c  Compute the fractional intersection location along side 1 by simple scaling
c    of triangle sides in 201 loop.
c
      do 201 is=1,kvec1
      dif=arc1(is+1)-rad
      if(dif*oldif.ge.0.d0) goto 201
      knt=knt+1
      if(knt.gt.2) goto 209
      ratdif=dabs(oldif)/(dabs(dif)+dabs(oldif))
      s1fr(knt)=dfloat(is-1)*cvec1+ratdif*cvec2
  201 oldif=dif
c
c  Save up to two intersection points for side 1, with a maximum total of two
c    for all three sides.
c
      frct=0.d0
      if(s1fr(1).gt.rangemin) frct=s1fr(1)
      if(s1fr(2).gt.rangemin) frct=s1fr(2)
      if(frct.eq.0.d0) goto 104
      delfr=0.d0
      iter=0
c
c  The 97 loop refines the side 1 intersection point, accurate within tolerance 'tolarc'
c
   97 continue
      iter=iter+1
      frct=frct+delfr
c
c   Compute numerical derivative, dardfr, and correction to the fractional
c    interval, delfr, in 96 loop and following two lines.
c
      do 96 iar=1,2
      if(iar.eq.2) frct=frct+frinc
      call vecin(vec2,vec3,crossdum,side1,veco,frct,kvec,2)
      call arcver(vec4,veco,vec3,dum1,dum2,arc,dum3,dum4,dum5,1,iflag)
      if(iar.eq.1) arcbase=arc
   96 continue
      dardfr=(arc-arcbase)/frinc
      delfr=(rad-arcbase)/dardfr
      if(dabs(delfr).gt.tolarc) goto 97
      s1fr(knt)=frct
c
c  Check side 2:
c
  104 if(ks2.eq.1) goto 105
      kside=2
      oldif=arc2(1)-rad
c
c  Compute the fractional intersection location along side 2 by simple scaling
c    of plane triangle sides in 202 loop.
c
      do 202 is=1,kvec1
      dif=arc2(is+1)-rad
      if(dif*oldif.ge.0.d0) goto 202
      knt=knt+1
      if(knt.gt.2) goto 209
      ratdif=dabs(oldif)/(dabs(dif)+dabs(oldif))
      s2fr(knt)=dfloat(is-1)*cvec1+ratdif*cvec2
  202 oldif=dif
c
c  Save up to two intersection points for side 2.
c
      frct=0.d0
      if(s2fr(1).gt.rangemin) frct=s2fr(1)
      if(s2fr(2).gt.rangemin) frct=s2fr(2)
      if(frct.eq.0.d0) goto 105
      delfr=0.d0
      iter=0
c
c  The 98 loop refines the side 2 intersection point, accurate within tolerance 'tolarc'
c
   98 continue
      iter=iter+1
      frct=frct+delfr
c
c   Compute numerical derivative, dardfr, and correction to the fractional
c    interval, delfr, in 56 loop and following two lines.
c
      do 56 iar=1,2
      if(iar.eq.2) frct=frct+frinc
      call vecin(vec3,vec1,crossdum,side2,veco,frct,kvec,2)
      call arcver(vec4,veco,vec3,dum1,dum2,arc,dum3,dum4,dum5,1,iflag)
      if(iar.eq.1) arcbase=arc
   56 continue
      dardfr=(arc-arcbase)/frinc
      delfr=(rad-arcbase)/dardfr
      if(dabs(delfr).gt.tolarc) goto 98
      s2fr(knt)=frct
c
c  Check side 3:
c
  105 if(ks3.eq.1) goto 209
      kside=3
      oldif=arc3(1)-rad
c
c  Compute the fractional intersection location along side 3 by simple scaling
c    of plane triangle sides in 203 loop.
c
      do 203 is=1,kvec1
      dif=arc3(is+1)-rad
      if(dif*oldif.ge.0.d0) goto 203
      knt=knt+1
      if(knt.gt.2) goto 209
      ratdif=dabs(oldif)/(dabs(dif)+dabs(oldif))
      s3fr(knt)=dfloat(is-1)*cvec1+ratdif*cvec2
  203 oldif=dif
c
c  Save up to two intersection points for side 3.
c
      frct=0.d0
      if(s3fr(1).gt.rangemin) frct=s3fr(1)
      if(s3fr(2).gt.rangemin) frct=s3fr(2)
      if(frct.eq.0.d0) goto 209
      delfr=0.d0
      iter=0
c
c  The 92 loop refines the side 3 intersection point, accurate within tolerance 'tolarc'
c
   92 continue
      iter=iter+1
      frct=frct+delfr
c
c   Compute numerical derivative, dardfr, and correction to the fractional
c    interval, delfr, in 88 loop and following two lines.
c
      do 88 iar=1,2
      if(iar.eq.2) frct=frct+frinc
      call vecin(vec1,vec2,crossdum,side3,veco,frct,kvec,2)
      call arcver(vec4,veco,vec3,dum1,dum2,arc,dum3,dum4,dum5,1,iflag)
      if(iar.eq.1) arcbase=arc
   88 continue
      dardfr=(arc-arcbase)/frinc
      delfr=(rad-arcbase)/dardfr
      if(dabs(delfr).gt.tolarc) goto 92
      s3fr(knt)=frct
  209 continue
c
c  Choose meaningful value (if any) for fractional arc along each triangle
c   side in next 6 lines (i.e. reject arcs represented by flags). 
c
      fr1=s1fr(1)
      if(fr1.lt.0.d0.and.s1fr(2).gt.0.d0) fr1=s1fr(2)
      fr2=s2fr(1)
      if(fr2.lt.0.d0.and.s2fr(2).gt.0.d0) fr2=s2fr(2)
      fr3=s3fr(1)
      if(fr3.lt.0.d0.and.s3fr(2).gt.0.d0) fr3=s3fr(2)
c
c  Next 6 lines: reverse the fractional arcs along the triangle sides if
c   they refer to a sense convention opposite to that applied in their
c   computation.
c
      if(kv1.eq.0.and.kv2.eq.1.and.kv3.eq.1) fr2=1.d0-fr2
      if(kv1.eq.1.and.kv2.eq.0.and.kv3.eq.1) fr3=1.d0-fr3
      if(kv1.eq.1.and.kv2.eq.1.and.kv3.eq.0) fr1=1.d0-fr1
      if(kv1.eq.1.and.kv2.eq.0.and.kv3.eq.0) fr2=1.d0-fr2
      if(kv1.eq.0.and.kv2.eq.1.and.kv3.eq.0) fr3=1.d0-fr3
      if(kv1.eq.0.and.kv2.eq.0.and.kv3.eq.1) fr1=1.d0-fr1
c
c       Compute sub-area as sub-triangle inside the circle, 
c     with eqn 2 (aclear=atri-asub) and vertex 1. Neglect 
c     departure of the circle from a great circle. Here it must be
c     that vertex 1 (i.e. vec1) is the one that is inside the circle
c     (so kv1=0, kv2=1, kv3=1).
c
      if(kv2.eq.0.or.kv3.eq.0.or.kv1.eq.1) goto 402
      call vecin(vec1,vec2,crossdum,side3,veco3,fr3,kvec,2)
      call vecin(vec1,vec3,crossdum,side2,veco2,fr2,kvec,2)
      call arcver(vec1,veco2,veco3,arcv1,arcv2,arcv3,ver1,ver2,ver3,3,
     $iflag)
      asub=dabs(ver1)+dabs(ver2)+dabs(ver3)-pi
      aclear=atri-asub
      return
c
c       Compute sub-area as sub-triangle inside the circle, 
c     with eqn 2 and vertex 2. Here it must be that vertex 2 (i.e. vec2)
c     is the one that is inside the circle (so kv1=1, kv2=0, kv3=1).
c
  402 if(kv1.eq.0.or.kv2.eq.1.or.kv3.eq.0) goto 403
      call vecin(vec2,vec1,crossdum,side3,veco3,fr3,kvec,2)
      call vecin(vec2,vec3,crossdum,side1,veco1,fr1,kvec,2)
      call arcver(vec2,veco1,veco3,arcv2,arcv1,arcv3,ver1,ver2,ver3,3,
     $iflag)
      asub=dabs(ver1)+dabs(ver2)+dabs(ver3)-pi
      aclear=atri-asub
      return
c
c       Compute sub-area as sub-triangle inside the circle, 
c     with eqn 2 and vertex 3. Here it must be that vertex 3 (i.e. vec3)
c     is the one that is inside the circle (so kv1=1, kv2=1, kv3=0).
c
  403 if(kv1.eq.0.or.kv2.eq.0.or.kv3.eq.1) goto 404
      call vecin(vec3,vec2,crossdum,side1,veco1,fr1,kvec,2)
      call vecin(vec3,vec1,crossdum,side2,veco2,fr2,kvec,2)
      call arcver(vec3,veco1,veco2,arcv3,arcv1,arcv2,ver1,ver2,ver3,3,
     $iflag)
      asub=dabs(ver1)+dabs(ver2)+dabs(ver3)-pi
      aclear=atri-asub
      return
c
c       Compute sub-area as sub-triangle outside the circle, 
c     with eqn 1 and vertex 1. Here it must be that vertex 1 (i.e. vec1)
c     is the one that is outside the circle (so kv1=1, kv2=0, kv3=0).
c
  404 if(kv1.eq.0.or.kv2.eq.1.or.kv3.eq.1) goto 405
      call vecin(vec1,vec2,crossdum,side3,veco3,fr3,kvec,2)
      call vecin(vec1,vec3,crossdum,side2,veco2,fr2,kvec,2)
      call arcver(vec1,veco2,veco3,arcv1,arcv2,arcv3,ver1,ver2,ver3,3,
     $iflag)
      asub=dabs(ver1)+dabs(ver2)+dabs(ver3)-pi
      aclear=asub
      return
c
c       Compute sub-area as sub-triangle outside the circle, 
c     with eqn 1 and vertex 2. Here it must be that vertex 2 (i.e. vec2)
c     is the one that is outside the circle (so kv1=0, kv2=1, kv3=0).
c
  405 if(kv1.eq.1.or.kv3.eq.1.or.kv2.eq.0) goto 406
      call vecin(vec2,vec1,crossdum,side1,veco3,fr3,kvec,2)
      call vecin(vec2,vec3,crossdum,side3,veco1,fr1,kvec,2)
      call arcver(vec2,veco1,veco3,arcv1,arcv2,arcv3,ver1,ver2,ver3,3,
     $iflag)
      asub=dabs(ver1)+dabs(ver2)+dabs(ver3)-pi
      aclear=asub
      return
c
c       Compute sub-area as sub-triangle outside the circle, 
c     with eqn 1 and vertex 3. Here it must be that vertex 3 (i.e. vec3)
c     is the one that is outside the circle (so kv1=0, kv2=0, kv3=1).
c
  406 if(kv1.eq.1.or.kv2.eq.1.or.kv3.eq.0) return
      call vecin(vec3,vec1,crossdum,side1,veco2,fr2,kvec,2)
      call vecin(vec3,vec2,crossdum,side3,veco1,fr1,kvec,2)
      call arcver(vec3,veco1,veco2,arcv1,arcv2,arcv3,ver1,ver2,ver3,3,
     $iflag)
      asub=dabs(ver1)+dabs(ver2)+dabs(ver3)-pi
      aclear=asub
      return
      end
      subroutine vecin(vec1,vec2,cross,vecout,veco,fra,kvec,kay)
c
c   March 16, 2011
c
c  Subroutine VECIN interpolates 'point vectors'* with uniform angular spacing 
c   along a great circle between a pair of unit input vectors (VEC1 and VEC2)
c   in 3-d. It works in either of two ways according to whether input integer
c   KAY is 1 or 2. It first computes the VEC1 X VEC2 cross product
c   so as to establish the plane and angular separation of those vectors. It
c   then rotates VEC1 about the cross product, either in equal angular steps 
c   (for KAY=1) or in one step (for KAY=2). 
c      For KAY=1 the number of interpolated
c   vectors is (input integer) kvec. Output array VECOUT contains the components
c   of VEC1 in its first 3 places, then components of all the interpolated
c   vectors, and then the components of VEC2 in its last 3 places (so it has all the
c   equally spaced vectors along a great circle arc, including the ends). 
c      For KAY=2, the fractional arc of the one interpolated vector between VEC1
c   and VEC2 is specified by input quantity FRA. Output array VECO contains
c   the interpolated vector components.
c      KVEC is ignored if KAY=2. FRA is ignored if KAY=1.
c
c   * Here 'point vector' means a vector that represents a point on the unit sphere.
c
      implicit real*8 (a-h,o-z)
      dimension vec1(3),vec2(3),veco(3),cross(3),rotmat(9),vecout(*)
      toldif=1.d-6
      cvec1=1.d0/dfloat(kvec+1)
c
c  Compute cross product components.
c
      call vcross(vec1,vec2,cross)
      knt=0
c
c  VEC1 and VEC2 have magnitude unity so the cross product's magnitude is sin(theta).
c
      sth=dsqrt(cross(1)**2+cross(2)**2+cross(3)**2)
      cth=dsqrt(dabs(1.d0-sth**2))
c
c  Compute the rotation matrix for vec1--> vec2
c
      uf=cross(1)
      vf=cross(2)
      wf=cross(3)
      uv=uf*vf
      uw=uf*wf
      vw=vf*wf
      usq=uf*uf
      vsq=vf*vf
      wsq=wf*wf
      den=usq+vsq+wsq
      denrecip=1.d0/den
      rden=dsqrt(den)
      onemcos=1.d0-cth
c
c  Rotate vector 1 about the cross product to co-incidence with vector 2.
c
      rotmat(1)=(usq+(vsq+wsq)*cth)*denrecip
      rotmat(5)=(vsq+(usq+wsq)*cth)*denrecip
      rotmat(9)=(wsq+(usq+vsq)*cth)*denrecip
   36 continue
      if(knt.gt.0) sth=-sth
      knt=1
      rotmat(2)=(uv*onemcos+wf*rden*sth)*denrecip
      rotmat(3)=(uw*onemcos-vf*rden*sth)*denrecip
      rotmat(4)=(uv*onemcos-wf*rden*sth)*denrecip
      rotmat(6)=(vw*onemcos+uf*rden*sth)*denrecip
      rotmat(7)=(uw*onemcos+vf*rden*sth)*denrecip
      rotmat(8)=(vw*onemcos-uf*rden*sth)*denrecip
      call dgmprd(rotmat,vec1,veco,3,3,1)
      theta=dasin(sth)
      vecdif=dabs(veco(1)-vec2(1))+dabs(veco(2)-vec2(2))+
     $dabs(veco(3)-vec2(3))
c
c  If sin(theta) has the wrong sign, go back and try with the opposite sign.
c
      if(vecdif.gt.toldif) goto 36
c
c  Copy the vec1 components into the first 3 elements of vecout.
c
      do 39 im=1,3
   39 vecout(im)=vec1(im)
c
c  Now rotate vector 1 about the cross product by fractional amounts
c    and put all the vector components in the middle elements of vecout.
c
      kmax=kvec
      if(kay.eq.2) kmax=1
      do 40 i40=1,kmax
      if(kay.eq.1) thet=theta*dfloat(i40)*cvec1
      if(kay.eq.2) thet=fra*theta
      sth=dsin(thet)
      cth=dcos(thet)
      onemcos=1.d0-cth
      rotmat(1)=(usq+(vsq+wsq)*cth)*denrecip
      rotmat(2)=(uv*onemcos+wf*rden*sth)*denrecip
      rotmat(3)=(uw*onemcos-vf*rden*sth)*denrecip
      rotmat(4)=(uv*onemcos-wf*rden*sth)*denrecip
      rotmat(5)=(vsq+(usq+wsq)*cth)*denrecip
      rotmat(6)=(vw*onemcos+uf*rden*sth)*denrecip
      rotmat(7)=(uw*onemcos+vf*rden*sth)*denrecip
      rotmat(8)=(vw*onemcos-uf*rden*sth)*denrecip
      rotmat(9)=(wsq+(usq+vsq)*cth)*denrecip
      call dgmprd(rotmat,vec1,veco,3,3,1)
      kv=3*(i40-1)+4
      vecout(kv)=veco(1)
      vecout(kv+1)=veco(2)
      vecout(kv+2)=veco(3)
   40 continue
c
c  Copy the vec2 components into the last 3 elements of vecout.
c
      imp=3*(kvec+1)
      do 41 im=1,3
      imp=imp+1
   41 vecout(imp)=vec2(im)
      return
      end
      subroutine arcver(vec1,vec2,vec3,arc1,arc2,arc3,vert1,vert2,
     $vert3,nsides,iflag)
      implicit real*8 (a-h,o-z)
c
c  Version of May 5, 2011
c
c   Subroutine ARCVER computes ARC lengths on a unit sphere (equivalent to
c    central angles) and VERtex angles, all in radians. Triangle
c    vertices are specified by position vectors with rectangular components.
c    Geometry is that of the unit sphere.
c    There are three cases according to whether input integer NSIDES 
c      is 1, 2, or 3:
c
c  NSIDES=1:   Compute one arc and (of course) no vertex angle.
c              Input is vec1 & vec2 (vec3 is a dummy, not used).
c              Output is arc3 (arc1 & arc2 are set to zero).
c              The three vertex angles are set to zero.  
c
c  NSIDES=2:   Compute two arcs and the included vertex angle.
c              Input is vec1, vec2, and vec3.
c              Output is arc2 and arc3 (arc1 is set to zero).
c              Vertex angle vert1 is computed (vert2 & vert3 are set to zero).
c
c  NSIDES=3:   Compute three arcs and three vertex angles of a spherical triangle.
c              Input is vec1, vec2, and vec3.
c              Output is arc1, arc2, and arc3
c              as well as vertex angles vert1, vert2, and vert3.
c
      dimension vec1(*),vec2(*),vec3(*)
      iflag=0
      arc1=0.d0
      arc2=0.d0
      vert1=0.d0
      vert2=0.d0
      vert3=0.d0
c
c  Compute the square of the 3-d straight line (i.e. not arc) distances between 
c    points on the unit sphere.
c
      side3sq=(vec2(1)-vec1(1))**2+(vec2(2)-vec1(2))**2
     $+(vec2(3)-vec1(3))**2
      if(nsides.ge.2) side2sq=(vec3(1)-vec1(1))**2+(vec3(2)-vec1(2))**2
     $+(vec3(3)-vec1(3))**2
      if(nsides.ge.2) side1sq=(vec2(1)-vec3(1))**2+(vec2(2)-vec3(2))**2
     $+(vec2(3)-vec3(3))**2
c
c  Compute cosines and sines of the corresponding arc lengths, and then the
c    arc lengths.
c
      cosang3=1.d0-.5d0*side3sq
      arc3=dacos(cosang3)
      if(nsides.eq.1) return
      sinang3=dsqrt(1.d0-cosang3**2)
      cosang2=1.d0-.5d0*side2sq
      sinang2=dsqrt(1.d0-cosang2**2)
      cosang1=1.d0-.5d0*side1sq
      sinang1=dsqrt(1.d0-cosang1**2)
      if(nsides.lt.3) goto 46
      arc2=dacos(cosang2)
      arc1=dacos(cosang1)
c
c  Compute cosines of the vertex angles from the spherical law of cosines, 
c    and then the angles. If the cosines are outside the allowable range
c    and NSIDES=3, set IFLAG=1 (meaning bad result) and return.
c
   46 continue
      cos1=(cosang1-cosang2*cosang3)/(sinang2*sinang3)
      if(dabs(cos1).lt.1.d0) goto 48
      iflag=1
      return
   48 continue
      vert1=dacos(cos1)
      cos2=(cosang2-cosang1*cosang3)/(sinang1*sinang3)
      if(dabs(cos2).lt.1.d0) goto 49
      iflag=1
      return
   49 continue
      vert2=dacos(cos2)
      cos3=(cosang3-cosang2*cosang1)/(sinang2*sinang1)
      if(dabs(cos3).lt.1.d0) goto 50
      iflag=1
      return
   50 continue
      vert3=dacos(cos3)
      return
      end
      subroutine vcross(vector1,vector2,vectorout)
      implicit real*8(a-h,o-z)
c
c   Version of February 9, 2011
c
c   Vcross computes the cross product components of a pair of
c      vectors in rectangular coordinates
c
      dimension vector1(3),vector2(3),vectorout(3)
      vectorout(1)=vector1(2)*vector2(3)-vector1(3)*vector2(2)
      vectorout(2)=vector1(3)*vector2(1)-vector1(1)*vector2(3)
      vectorout(3)=vector1(1)*vector2(2)-vector1(2)*vector2(1)
      return
      end
      subroutine spotev(hjd,tstart,tmaxst,tmaxsp,tfinal,amax,area,sprad)
      implicit real*8 (a-h,o-z)
c   Version of October 27, 2011
      pi=dacos(-1.d0)
      twopi=pi+pi
      area=0.d0
      sla=amax/(tmaxst-tstart)
      slb=amax/(tmaxsp-tfinal)
      if(hjd.gt.tstart.and.hjd.le.tmaxst) area=sla*(hjd-tstart)
      if(hjd.gt.tmaxst.and.hjd.le.tmaxsp) area=amax
      if(hjd.gt.tmaxsp.and.hjd.lt.tfinal) area=slb*(hjd-tfinal)
      sprad=dacos(1.d0-area/twopi)
      return
      end
      subroutine daspdt (time,amax,t1,t2,t3,t4,dadt1,dadt2,dadt3,dadt4)
      implicit real*8 (a-h,o-z)
c  Version of October 5, 2011
      dadt1=0.d0
      dadt2=0.d0
      dadt3=0.d0
      dadt4=0.d0
      if(time.le.t1.or.time.ge.t4) return
      if(time.ge.t2.and.time.le.t3) return
      if(time.gt.t1.and.time.lt.t2) dadt1=amax*(time-t2)/(t2-t1)**2
      if(time.gt.t1.and.time.lt.t2) dadt2=amax*(t1-time)/(t2-t1)**2
      if(time.gt.t3.and.time.lt.t4) dadt3=amax*(t4-time)/(t4-t3)**2
      if(time.gt.t3.and.time.lt.t4) dadt4=amax*(time-t3)/(t4-t3)**2
      return
      end
      SUBROUTINE ELLONE(FF,dd,rm,xl1,OM1,XL2,OM2)                         
c  Version of December 4, 2003                                             
C     XL2 AND OM2 VALUES ASSUME SYNCHRONOUS ROTATION AND CIRCULAR ORBIT.  
C     THEY ARE NOT NEEDED FOR NON-SYNCHRONOUS OR NON-CIRCULAR CASES.      
c
c  Starting on August 13, 2003, ELLONE includes a 2nd derivative term 
c    in the N-R solution for the null point of effective gravity (XL1)  
c 
      IMPLICIT REAL*8(A-H,O-Z)                                            
      COMMON /ECCEN/ecc,smaxis,period,vgadum,sindum,vfdum,vfadum,vgmdum,
     $v1dum,v2dum,ifcdum
      ot=1.d0/3.d0 
      icase=2
      if(ff.ne.1.d0.or.ecc.gt.0.d0) icase=1
      rmass=rm                                                            
      d=dd                                                                
      xl=d/(1.d0+dsqrt(rm))
      oldxl=xl
      DO 5 I=1,icase                                                          
      RFAC=ff*ff                                                          
      IF(I.EQ.2) RFAC=1.D0                                                
      IF(I.EQ.2) D=1.D0                                                   
      DSQ=D*D                                                             
      DELXL=0.D0                                                          
      RM1=RMASS+1.D0                                                      
      kount=0 
   88 XL=XL+DELXL                                                         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  The next block of lines halves the delxl step in case
c  xl were to jump beyond the value d or below value
c  0.0 during the iteration.
      if(i.eq.2) goto 170
      if(xl.lt.dd.and.xl.gt.0.d0) goto 170
      delxl=0.5d0*delxl
      xl=oldxl
      goto 88
  170 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      kount=kount+1 
      XSQ=XL*XL                                                           
      P=(D-XL)**2                                                         
      RP=DABS(D-XL)                                                       
      PRP=P*RP                                                            
      F=RFAC*RM1*XL-1.D0/XSQ-RMASS*(XL-D)/PRP-RMASS/DSQ                   
      DXLDF=1.D0/(RFAC*RM1+2.D0/(XSQ*XL)+2.D0*RMASS/PRP)                  
      d2xldf2=6.d0*dxldf**3*(1.d0/xl**4-rmass/rp**4) 
      DELXL=-F*DXLDF+.5d0*f*f*d2xldf2                                                      
      ABDEL=DABS(DELXL)                                                   
      oldxl=xl
      IF(ABDEL.GT.1.d-10) GOTO 88                                      
      IF(I.EQ.2) GOTO 8                                                   
      xl1=xl                                                              
      OM1=1.D0/XL+RMASS*((1.D0/RP)-XL/DSQ)+RM1*.5D0*XSQ*RFAC              
      IF(rm.GT.1.d0) RMASS=1.D0/RMASS                                     
      XMU3=RMASS/(3.D0*(RMASS+1.D0))                                      
      XMU3CR=XMU3**ot                                         
    5 XL=1.D0+XMU3CR+XMU3CR*XMU3CR/3.D0+XMU3/9.D0                         
    8 IF(rm.GT.1.d0) XL=D-XL                                              
      rm1=rm+1.d0 
      XL2=XL                                                              
      OM2=1.D0/DABS(XL)+rm*((1.D0/DSQRT(1.D0-XL-XL+XL*XL))-XL)+RM1*    
     $.5D0*XL*XL                                                          
      RETURN                                                              
      END                                                                 
      SUBROUTINE VOLUME(V,Q,P,D,FF,N,N1,KOMP,RV,GRX,GRY,GRZ,RVQ,         
     $GRXQ,GRYQ,GRZQ,MMSAVE,FR1,FR2,HLD,SNTH,CSTH,SNFI,CSFI,SUMM,SM,     
     $GRV1,GRV2,XX1,YY1,ZZ1,XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2,GMAG1  
     $,GMAG2,glog1,glog2,GREXP,IFC)                                                  
c  Version of December 5, 2003                                         
      implicit real*8 (a-h,o-z)                                          
      DIMENSION RV(*),GRX(*),GRY(*),GRZ(*),RVQ(*),GRXQ(*),GRYQ(*),GRZQ(* 
     $),MMSAVE(*),FR1(*),FR2(*),HLD(*),SNTH(*),CSTH(*),SNFI(*),CSFI(*)   
     $,GRV1(*),GRV2(*),GLUMP1(*),GLUMP2(*),XX1(*),YY1(*),ZZ1(*),XX2(*),  
     $YY2(*),ZZ2(*),CSBT1(*),CSBT2(*),GMAG1(*),GMAG2(*),glog1(*),                 
     $glog2(*) 
      if(ifc.eq.1) v=0.d0
      DP=1.d-5*P                                                         
      ot=1.d0/3.d0 
      IF (IFC.EQ.1) DP=0.d0                                              
      tolr=1.d-8 
      DELP=0.d0                                                          
      KNTR=0                                                             
   16 P=P+DELP                                                           
      KNTR=KNTR+1                                                        
      IF(KNTR.GE.20) tolr=tolr+tolr                                         
      PS=P                                                               
      DO 17 I=1,IFC                                                      
      P=PS                                                               
      IF(I.EQ.1) P=P+DP                                                  
      CALL SURFAS(Q,P,N,N1,KOMP,RV,GRX,GRY,GRZ,RVQ,GRXQ,GRYQ,GRZQ,       
     $MMSAVE,FR1,FR2,HLD,FF,D,SNTH,CSTH,SNFI,CSFI,GRV1,GRV2,XX1,YY1,ZZ1, 
     $XX2,YY2,ZZ2,CSBT1,CSBT2,GLUMP1,GLUMP2,GMAG1,GMAG2,glog1,glog2,          
     $grexp) 
      IF(KOMP.EQ.2) GOTO 14                                              
      call lum(1.d0,1.d0,0.d0,1.d0,n,n1,1,sbrd,rv,rvq,glump1,       
     $glump2,glog1,glog2,grv1,grv2,mmsave,summ,fr1,sm,0,vol,q,p,ff,d,         
     $snth,7) 
      GOTO 15                                                            
   14 call lum(1.d0,1.d0,0.d0,1.d0,n,n1,2,sbrd,rv,rvq,glump1,       
     $glump2,glog1,glog2,grv1,grv2,mmsave,summ,fr2,sm,0,vol,q,p,ff,d,         
     $snth,7) 
   15 CONTINUE                                                           
      IF(I.EQ.1) VOLS=VOL                                                
      VOL2=VOLS                                                          
   17 VOL1=VOL                                                           
      rmean=(.238732414d0*vol)**ot 
      rmsq=rmean**2 
c 
c  Here use a polar estimate for d(potential)/dr (absolute value). 
c 
      domdrabs=1.d0/rmsq+q*rmean/(d*d+rmsq) 
      tolp=domdrabs*tolr 
      IF(IFC.EQ.1) V=VOL                                                 
      IF(IFC.EQ.1) RETURN                                                
      DPDV=DP/(VOL2-VOL1)
      DELP=(V-VOL1)*DPDV                                                 
      ABDELP=dabs(DELP)                                                  
      IF(ABDELP.GT.tolp) GOTO 16                                          
      P=PS                                                               
      RETURN                                                             
      END                                                                
      subroutine gaussquad (kway,n,emmm,ff,out)
c  Version of February 17, 2012
c  
c  Input N is the number of Gaussian abscissas, which are from a
c    website by P. Holoborodko. 
c    (http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/)
c 
c  The present version works only for N of 2, 3, 4, 5, 6, 7, 8, 9, or 10.
c  Usually even low N (say 3 or 4) gives very good accuracy in cases where
c  the integrated function is reasonably smooth over the range (i.e. not
c  too "busy").

c  Input EMMM is a scale conversion factor (half the range, 0.5*(B-A),
c  where B and A are the upper and lower range limits.

c  The calling program needs to know the Gaussian abscissas so it can
c  compute the input function (ff) at the correct values of the 
c  independent variable. So the subroutine can be called in two ways
c  according to whether control integer KWAY is 1 or 2, as follows:
c
c  KWAY=1: Return only the Gauss abscissas. 

c  KWAY=2: Return the integrated result.
c
c  Output OUT is the integrated result. 
c
      implicit real*8 (a-h,o-z)
      dimension ff(*),wt(30)
      common /abscissa/ abs(29),ind0(10)
      data ind0(2),ind0(3),ind0(4),ind0(5),ind0(6),ind0(7),ind0(8),
     $ind0(9),ind0(10)/1,2,4,6,9,12,16,20,25/
      data abs(1),abs(2),abs(3),abs(4),abs(5),abs(6),abs(7),abs(8),
     $abs(9),abs(10),abs(11),abs(12),abs(13),abs(14),abs(15),abs(16),
     $abs(17),abs(18),abs(19),abs(20),abs(21),abs(22),abs(23),abs(24),
     $abs(25),abs(26),abs(27),abs(28),abs(29)/.5773502691896258d0,
     $.7745966692414834d0,.0d0,.8611363115940526d0,.3399810435848563d0,
     $.9061798459386640d0,.5384693101056831d0,.0d0,.9324695142031520d0,
     $.6612093864662645d0,.2386191860831969d0,.9491079123427585d0,
     $.7415311855993944d0,.4058451513773972d0,.0d0,.9602898564975362d0,
     $.7966664774136267d0,.5255324099163290d0,.1834346424956498d0,
     $.9681602395076261d0,.8360311073266358d0,.6133714327005904d0,
     $.3242534234038089d0,.0d0,.9739065285171717d0,.8650633666889845d0,
     $.6794095682990244d0,.4333953941292472d0,.1488743389816312d0/
      data wt(1),wt(2),wt(3),wt(4),wt(5),wt(6),wt(7),wt(8),wt(9),wt(10),
     $wt(11),wt(12),wt(13),wt(14),wt(15),wt(16),wt(17),wt(18),wt(19),
     $wt(20),wt(21),wt(22),wt(23),wt(24),wt(25),wt(26),wt(27),wt(28),
     $wt(29)/1.d0,.5555555555555556d0,.8888888888888889d0,
     $.3478548451374539d0,.6521451548625461d0,.2369268850561891d0,
     $.4786286704993665d0,.5688888888888889d0,.1713244923791703d0,
     $.3607615730481386d0,.4679139345726910d0,.1294849661688697d0,
     $.2797053914892767d0,.3818300505051189d0,.4179591836734694d0,
     $.1012285362903763d0,.2223810344533745d0,.3137066458778873d0,
     $.3626837833783620d0,.0812743883615477d0,.1806481606948574d0,
     $.2606106964029355d0,.3123470770400028d0,.3302393550012598d0,
     $.0666713443088688d0,.1494513491505806d0,.2190863625159820d0,
     $.2692667193099964d0,.2955242247147529d0/
      if(kway.eq.1) return
      no2=n/2
      out=0.d0
      do 88 i88=1,n
      if(i88.le.no2) ind=ind0(n)+i88-1
      if(i88.gt.no2) ind=ind0(n)+n-i88
      out=out+wt(ind)*ff(i88)
   88 continue
      out=emmm*out
      return
      end
      subroutine splinterpol(x,y,y2,n,xx,yy)
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n),y2(n)
      jlow=1
      jhigh=n
  1   jhmjl=jhigh-jlow
      if(jhmjl.le.1) goto 4
      j=(jhigh+jlow)/2
      if(x(j).le.xx) goto 2
      jhigh=j
      goto 3
  2   continue
      jlow=j
  3   continue
      go to 1
  4   continue
      dif=x(jhigh)-x(jlow)
      if(dif.eq.0.) goto 5
      a=(x(jhigh)-xx)/dif
      b=(xx-x(jlow))/dif
      yy=a*y(jlow)+b*y(jhigh)+
     $((a**3-a)*y2(jlow)+(b**3-b)*y2(jhigh))*(dif**2)/6.d0
      return
  5   continue
      write(6,80)
      stop
  80  format('Program stopped in splinterpol, problem in x array.')
      end
      subroutine atmcof(iab) 
      implicit real*8 (a-h,o-z)
      parameter (nbmax=100,ntmax=80,ngmax=11)
      parameter (ntmax2=2*ntmax)
      common /arrayspline/ grand(ntmax2,ngmax,nbmax),
     $tem(ntmax2,ngmax,nbmax),yf(ntmax2,ngmax,nbmax)
      common /numtemp/ kend1(ngmax,nbmax),kend2(ngmax,nbmax)
      common /nummod/ nbp,ntemp,ngl
      open(unit=22,file='atmcof.dat',status='old')
      do 6880 i=1,nbmax
      do 6880 j=1,ngmax
      do 6880 k=1,ntmax2
      grand(k,j,i)=0.d0
 6880 continue
      do 58 kfil=1,nbp
      if (kfil.ne.1) goto 6771
      irec=33*(iab-1)
      do 6770 krec=1,irec
      read(22,865) ndum
 6770 continue
      goto 6773
 6771 continue
      irec=18*33
      do 6772 krec=1,irec
      read(22,865) ndum
 6772 continue
 6773 continue
      do 68 ig=1,ngl
      read(22,865,advance='no') nbin,kend1(ig,kfil)
      ken1=kend1(ig,kfil)
      if(nbin.eq.1) goto 6881
      read(22,866,advance='no') (grand(k,ig,kfil),k=1,ken1)
      read(22,867,advance='no') kend2(ig,kfil)
      ken2=kend2(ig,kfil)
      read(22,866) (grand(k+80,ig,kfil),k=1,ken2)
      read(22,868,advance='no') nbin,ken1,(tem(k,ig,kfil),k=1,ken1)
      read(22,869) ken2,(tem(k+80,ig,kfil),k=1,ken2)
      read(22,868,advance='no') nbin,ken1,(yf(k,ig,kfil),k=1,ken1)
      read(22,869) ken2,(yf(k+80,ig,kfil),k=1,ken2)
      goto 6882
 6881 continue
      read(22,866) (grand(k,ig,kfil),k=1,ken1)
      read(22,868) nbin,ken1,(tem(k,ig,kfil),k=1,ken1)
      read(22,868) nbin,ken1,(yf(k,ig,kfil),k=1,ken1)
 6882 continue
  68  continue
  58  continue
      close (22)
      return
 865  format(2i3)
 866  format(80D17.9)
 867  format(i3)
 868  format(2i3,80D17.9)
 869  format(i3,80D17.9)
      end
      subroutine mintype (ecc,argper,phzero,hjd0,period,dpdt,hjd,mtype)
c  Version of October 15, 2012
      implicit real*8 (a-h,o-z)
      phin=0.d0 
      call conjph(ecc,argper,phzero,trsc,tric,econsc,econic,xmsc,xmic,
     $pconsc,pconic)
   14 if(pconsc.lt.0.d0) pconsc=pconsc+1.d0
      if(pconsc.lt.0.d0) goto 14
   16 if(pconsc.gt.1.d0) pconsc=pconsc-1.d0
      if(pconsc.gt.1.d0) goto 16
   18 if(pconic.lt.0.d0) pconic=pconic+1.d0
      if(pconic.lt.0.d0) goto 18
   19 if(pconic.gt.1.d0) pconic=pconic-1.d0
      if(pconic.gt.1.d0) goto 19
      call jdph(hjd,phin,hjd0,period,dpdt,hjdout,faz)
   10 if(faz.gt.1.d0) faz=faz-1.d0
      if(faz.gt.1.d0) goto 10
   11 if(faz.lt.0.d0) faz=faz+1.d0
      if(faz.lt.0.d0) goto 11
      fazneg=faz-1.d0
      absdifs=dabs(faz-pconsc)
      absdifi=dabs(faz-pconic)
      absdifnegs=dabs(fazneg-pconsc)
      absdifnegi=dabs(fazneg-pconic)
      winsup=absdifs
      if(absdifnegs.lt.absdifs) winsup=absdifnegs
      wininf=absdifi
      if(absdifnegi.lt.absdifi) wininf=absdifnegi
      mtype=1
      if(wininf.lt.winsup) mtype=2
      return
      end
