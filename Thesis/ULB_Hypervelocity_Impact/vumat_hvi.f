C =====================================================================
C  vumat_hvi.f -- Abaqus/Explicit user material (VUMAT) for
C  hypervelocity impact on ceramic/aluminium shields:
C    Al2024 : Mie-Grueneisen EOS + Johnson-Cook plasticity and damage
C    NbC    : Johnson-Holmquist 2 ceramic model with bulking
C    one material for every body, FEM-to-SPH conversion flag (SDV 4)
C
C  Copyright 2026 Riccardo Savio
C  SPDX-License-Identifier: Apache-2.0
C
C  Licensed under the Apache License, Version 2.0 (see LICENSE).
C  Redistributions and derivative works must keep this notice and the
C  NOTICE file. Please cite this work if you use it (CITATION.cff).
C  https://github.com/riccardosavio02/academic-portfolio
C =====================================================================
      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, jInfoArray,
     2  stepTime, totalTime, dtArray, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), dtArray(2*(nblock)+1), 
     1  strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock), 
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock), jInfoArray(*)
C
      character*80 cmname
C
       parameter (zero = 0.d0, one = 1.d0, two = 2.d0)
       parameter (three = 3.d0, third = 1.d0/3.d0)
       parameter (half = 1.d0/2.d0)
       parameter (i_info_AnnealFlag = 1, i_info_Intpt = 2,
     1 i_info_layer = 3, i_info_kspt = 4, i_info_effModDefn = 5,
     2 i_info_ElemNumStartLoc = 6)
C      state variables (README): 1 EPSP equivalent plastic strain,
C      2 DMG JC damage, 3 EPSV accumulated volumetric log strain,
C      4 CONV conversion flag (1 = stays FEM, 0 = convert),
C      5 SIG1 largest principal stress, 6 PRESS pressure
       parameter (iEpsp = 1, iDmg = 2, iEpsV = 3)
       parameter (iConv = 4, iSig1 = 5, iPress = 6)
C
       pointer (ptrjElemNum, jElemNum)
       dimension jElemNum(nblock)
C      principal stresses of the block, from the utility vsprinc. The
C      manual declares eigVal(nblock,3): the first dimension is the REAL
C      block size, not maxblk. Stored as a vector, eigVal(km,j) is
C      eigVal(km + (j-1)*nblock). With (maxblk,3) the partial blocks read
C      uninitialised memory as sigma_1: spurious conversions
       dimension eigVal(3*maxblk)
C      trial and new deviator of one material point (6 components)
       dimension sTr(6), sNew(6)
C      state of every parent element, kept between calls in a common
C      block indexed by element label: internal and dissipated energy
C      per unit mass, EPSV, EPSP, DMG, and parCv = 1 once the element
C      has converted. A particle reads its parent's state at birth
C      (energy transfer, README). One process only (cpus=1, as every
C      deck with conversion); labels above maxPar are not stored.
       parameter (maxPar = 3000000)
       common /hvipar/ parEi(maxPar), parEd(maxPar), parEv(maxPar),
     1                 parEp(maxPar), parDm(maxPar), parCv(maxPar)
C      position and size of each parent at its last call: a particle
C      takes the state only if it lies within one element size of its
C      supposed parent (guards against a wrong label rule)
       common /hvipos/ parX(maxPar), parY(maxPar), parZ(maxPar),
     1                 parL(maxPar)
       integer km, i
C      1 = this point of the block is NbC (hvi_nbc), 0 = Al2024
       dimension iNbc(maxblk)
C
C      NbC (F4): a material whose name starts with NBC is Johnson-
C      Holmquist 2 with the polynomial EOS, in hvi_nbc (constants in the
C      order of the built-in ABQ_JH2, README "NbC branch")
       if (cmname(1:3) .eq. 'NBC') then
          do km = 1, nblock
             iNbc(km) = 1
          end do
          call hvi_nbc(nblock, ndir, nshr, nstatev, nprops, jInfoArray,
     1         stepTime, totalTime, dtArray, props, strainInc,
     2         stressOld, stateOld, enerInternOld, enerInelasOld,
     3         stressNew, stateNew, enerInternNew, enerInelasNew,
     4         iNbc, 0)
          return
       end if
C      F5: a material whose name starts with HVI holds Al2024 AND NbC
C      (one section for every body, README "One material"): props 1-19
C      Al2024 as above, 20-21 first and last label of the NbC elements,
C      22 energy transfer flag (1 = on), 25-56 the 32 NbC constants
C      (ABQ_JH2 order). The branch of each point is chosen by the label
C      of its element, for a particle by the label of its parent.
       iMulti = 0
       if (cmname(1:3) .eq. 'HVI') iMulti = 1
C
C      props (README): 1 G, 2 rho0, 3 c0, 4 s, 5 Gamma0 (Mie-Grueneisen
C      Us-Up), 6-10 JC A, B, n, C, epsdot0, 11-15 JC damage D1-D5,
C      16 sigma_max, 17 number of solid elements (a particle has a
C      larger element number)
       gmod = props(1)
       rho0 = props(2)
       c0 = props(3)
       s = props(4)
       gam0 = props(5)
C      Johnson-Cook yield (no thermal term). The exponent n is called
C      expn: a name starting with n would be an integer.
       a = props(6)
       b = props(7)
       expn = props(8)
       cjc = props(9)
       edot0 = props(10)
C      Johnson-Cook damage (Moore 2025 eq. 14, same form as Abaqus)
       d1 = props(11)
       d2 = props(12)
       d3 = props(13)
       d4 = props(14)
       d5 = props(15)
       sigMax = props(16)
       nSolid = nint(props(17))
C      props 24 = eps_c, conversion threshold on the
C      equivalent plastic strain (an element also converts when EPSP
C      reaches eps_c; 0 or absent = off, as vumat_hvi.f). The analogue
C      of Moore's calibrated eps_f of the NbC, for the Al2024 branch.
       epsC = zero
       if (nprops .ge. 24) epsC = props(24)
C      18 particles per direction PPD, 19 label of the first element of
C      this material. Abaqus numbers the particle of parent ordinal k as
C      Nsolid + (k - 1) PPD^3 + 1..PPD^3 (stage 0), so with contiguous
C      labels the parent of particle L is lPar0 + (L - Nsolid - 1)/PPD^3
C      (integer division). Energy transfer at birth only with 19 props
C      (iTransf = 1); with 17 the particles keep their own state.
       iTransf = 0
       nPpd3 = 1
       lPar0 = 1
       if (nprops .ge. 19) then
          iTransf = 1
          nPpd3 = nint(props(18))**3
          lPar0 = nint(props(19))
       end if
C      With 22 props or more the transfer flag is props(22) also
C      for a single material (a 24-props deck must not switch it on)
       if (nprops .ge. 22) iTransf = nint(props(22))
       if (iMulti .eq. 1) then
          iTransf = nint(props(22))
          lN0 = nint(props(20))
          lN1 = nint(props(21))
          lS0 = nint(props(23))
       end if
       twomu = two*gmod
C      time increment
       dt = dtArray(1)
C      bulk modulus of the EOS at rho = rho0 (packaging call)
       bulk = rho0*c0**2
C nblock : number of material points to be processed in this call to VUMAT
       iElemNumStartLoc = jInfoArray(i_info_ElemNumStartLoc)
       ptrjElemNum = loc(jInfoArray(iElemNumStartLoc))
C
C      F5: which points of the block are NbC, then hvi_nbc does them
       do km = 1, nblock
          iNbc(km) = 0
          if (iMulti .eq. 1) then
             lab = jElemNum(km)
             if (lab .gt. nSolid) then
                ip = lPar0 + (lab - nSolid - 1)/nPpd3
             else
                ip = lab
             end if
             if (ip .ge. lN0 .and. ip .le. lN1) iNbc(km) = 1
          end if
       end do
       if (iMulti .eq. 1) then
          call hvi_nbc(nblock, ndir, nshr, nstatev, 32, jInfoArray,
     1         stepTime, totalTime, dtArray, props(25), strainInc,
     2         stressOld, stateOld, enerInternOld, enerInelasOld,
     3         stressNew, stateNew, enerInternNew, enerInelasNew,
     4         iNbc, nSolid)
       end if
C
C      packaging call (first call, totalTime = stepTime = 0): the strain
C      increments are fictitious and only serve Abaqus to estimate the
C      wave speed. Linear elastic answer (K = rho0 c0^2, G), state
C      variables and energies left as they are.
       if (totalTime .eq. zero .and. stepTime .eq. zero) then
          alamda = bulk - twomu*third
          do km = 1, nblock
             if (iNbc(km) .eq. 1) cycle
             trace = strainInc(km,1) + strainInc(km,2) + strainInc(km,3)
             do i = 1, 3
                stressNew(km,i) = stressOld(km,i) + alamda*trace
     1                          + twomu*strainInc(km,i)
             end do
             do i = 4, 6
                stressNew(km,i) = stressOld(km,i)
     1                          + twomu*strainInc(km,i)
             end do
             do i = 1, nstatev
                stateNew(km,i) = stateOld(km,i)
             end do
             enerInternNew(km) = enerInternOld(km)
             enerInelasNew(km) = enerInelasOld(km)
          end do
          return
       end if
C
       do km = 1, nblock
            if (iNbc(km) .eq. 1) cycle
            trace = strainInc(km,1) + strainInc(km,2) + strainInc(km,3)
            lab = jElemNum(km)
C           state variables beyond the Al2024 six (NbC's DPBULK, EPRATE
C           in a multi-material deck) are carried unchanged
            do i = 7, nstatev
               stateNew(km,i) = stateOld(km,i)
            end do
C           energies at the start of the increment, per unit mass
            eOld = enerInternOld(km)
            edOld = enerInelasOld(km)
C           plastic strain and damage unchanged unless the element
C           yields
            stateNew(km,iEpsp) = stateOld(km,iEpsp)
            stateNew(km,iDmg) = stateOld(km,iDmg)
            epsVOld = stateOld(km,iEpsV)
C           birth of a particle: its SDV4 is still 1 (2 once born) and
C           its parent has converted. It takes the parent's energies,
C           EPSV, EPSP and DMG, so that the energy of the parent stays
C           in the books (Abaqus drops the parent's) and the pressure is
C           continuous. Before that, a dormant particle keeps its own.
            iBorn = 0
            if (iTransf .eq. 1 .and. lab .gt. nSolid .and.
     1          stateOld(km,iConv) .lt. 1.5) then
               ip = lPar0 + (lab - nSolid - 1)/nPpd3
C              the parent must be an element of this material
               if (ip .ge. max(1, lPar0) .and. ip .le. nSolid .and.
     1             ip .le. maxPar) then
                  if (parCv(ip) .gt. half) then
                     dist = sqrt((coordMp(km,1) - parX(ip))**2
     1                         + (coordMp(km,2) - parY(ip))**2
     2                         + (coordMp(km,3) - parZ(ip))**2)
                     if (dist .le. parL(ip)) then
                        iBorn = 1
                     else
C                       not its parent: keep the own state, and never
C                       try again (SDV4 = 2 below)
                        iBorn = -1
                     end if
                  end if
               end if
            end if
            if (iBorn .eq. 1) then
               eOld = parEi(ip)
               edOld = parEd(ip)
               epsVOld = parEv(ip)
               stateNew(km,iEpsp) = parEp(ip)
               stateNew(km,iDmg) = parDm(ip)
            end if
C           volume ratio J = V/V0 = rho0/rho from the accumulated
C           volumetric log strain (more accurate than density(km))
            epsV = epsVOld + trace
            stateNew(km,iEpsV) = epsV
            aj = exp(epsV)
            rho = rho0/aj
            dep = zero
            sy = zero
            pOld = -third*(stressOld(km,1) + stressOld(km,2)
     1                   + stressOld(km,3))
            if (jElemNum(km) .gt. nSolid) then
C              particle = failed material (He 2020): no deviatoric
C              stress, pressure only, no tension. Pressure from the
C              same Mie-Grueneisen EOS, with the particle's own J (its
C              EPSV, built from the strain increments it receives while
C              dormant) and its own internal energy; from its birth on,
C              the parent's state carried forward.
               call hvi_eos_mg(rho0, aj, eOld, c0, s, gam0, p)
               if (p .lt. zero) p = zero
C              SDV4 of a particle: 1 before birth, 2 from birth on
C              (3 if the parent found by the label rule was too far:
C              the transfer was refused, see SDV4 in the odb)
               if (iBorn .eq. -1) then
                  stateNew(km,iConv) = three
               else if (iBorn .eq. 1 .or.
     1                  stateOld(km,iConv) .gt. 1.5) then
                  stateNew(km,iConv) = max(two, stateOld(km,iConv))
               else
                  stateNew(km,iConv) = one
               end if
               do i = 1, 3
                  stressNew(km,i) = -p
               end do
               do i = 4, 6
                  stressNew(km,i) = zero
               end do
            else
C              element: pressure from the Mie-Grueneisen EOS with the
C              internal energy of the previous increment
               call hvi_eos_mg(rho0, aj, eOld, c0, s, gam0, p)
C              trial deviator s_tr = s_old + 2G de, de = dEps - tr/3
C              on the diagonal (s_old = sigma_old + p_old)
               do i = 1, 3
                  sTr(i) = stressOld(km,i) + pOld
     1                   + twomu*(strainInc(km,i) - third*trace)
               end do
               do i = 4, 6
                  sTr(i) = stressOld(km,i) + twomu*strainInc(km,i)
               end do
C              JC plasticity: radial return to the yield surface
               call hvi_radial_return(sTr, gmod, stateOld(km,iEpsp),
     1                  dt, a, b, expn, cjc, edot0, sNew, dep, sy)
               stateNew(km,iEpsp) = stateOld(km,iEpsp) + dep
C              JC damage, only when the element yields. Triaxiality
C              -p/q with q = sy, the Mises stress after the return.
               if (dep .gt. zero) then
                  triax = -p/sy
                  call hvi_jc_damage(stateOld(km,iDmg), dep, triax,
     1                 dep/dt, d1, d2, d3, d4, d5, edot0,
     2                 stateNew(km,iDmg))
               end if
               do i = 1, 3
                  stressNew(km,i) = sNew(i) - p
               end do
               do i = 4, 6
                  stressNew(km,i) = sNew(i)
               end do
            end if
            stateNew(km,iPress) = p
C           stress power per unit volume, 1/2 (sOld+sNew):dEps,
C           shear components counted twice (symmetric tensor)
            stressPower = half*(
     1        (stressOld(km,1)+stressNew(km,1))*strainInc(km,1)
     2      + (stressOld(km,2)+stressNew(km,2))*strainInc(km,2)
     3      + (stressOld(km,3)+stressNew(km,3))*strainInc(km,3)
     4      + two*(stressOld(km,4)+stressNew(km,4))*strainInc(km,4)
     5      + two*(stressOld(km,5)+stressNew(km,5))*strainInc(km,5)
     6      + two*(stressOld(km,6)+stressNew(km,6))*strainInc(km,6) )
            enerInternNew(km) = eOld + stressPower/rho
C           dissipated plastic work per unit mass
            enerInelasNew(km) = edOld + sy*dep/rho
C           an element leaves its state for its particles (only
C           needed with the energy transfer)
            if (iTransf .eq. 1 .and. lab .le. nSolid .and.
     1          lab .ge. 1 .and. lab .le. maxPar) then
               parEi(lab) = enerInternNew(km)
               parEd(lab) = enerInelasNew(km)
               parEv(lab) = epsV
               parEp(lab) = stateNew(km,iEpsp)
               parDm(lab) = stateNew(km,iDmg)
               parX(lab) = coordMp(km,1)
               parY(lab) = coordMp(km,2)
               parZ(lab) = coordMp(km,3)
               parL(lab) = charLength(km)
            end if
       end do
C
C      conversion flag (He 2020, Moore 2025): an element converts when
C      the JC damage reaches 1 or sigma_1 reaches sigma_max. Principal
C      values of the new stress for the whole block (unsorted), then
C      sigma_1 = the largest.
       call vsprinc(nblock, stressNew, eigVal, ndir, nshr)
       do km = 1, nblock
            if (iNbc(km) .eq. 1) cycle
            sig1 = max(eigVal(km), eigVal(km+nblock),
     1                 eigVal(km+2*nblock))
            stateNew(km,iSig1) = sig1
            lab = jElemNum(km)
            if (lab .gt. nSolid) then
C              a particle cannot convert (its SDV4 is set above)
               continue
            else if (stateOld(km,iConv) .eq. zero) then
C              conversion is one-way: once 0, stays 0
               stateNew(km,iConv) = zero
            else
               if (sig1 .ge. sigMax .or.
     1             stateNew(km,iDmg) .ge. one .or.
     2             (epsC .gt. zero .and.
     3              stateNew(km,iEpsp) .ge. epsC)) then
                  stateNew(km,iConv) = zero
               else
                  stateNew(km,iConv) = one
               end if
            end if
C           tell the particles that their parent has converted
            if (iTransf .eq. 1 .and. lab .le. nSolid .and.
     1          lab .ge. 1 .and. lab .le. maxPar) then
               if (stateNew(km,iConv) .eq. zero) parCv(lab) = one
            end if
       end do
C
C      multi-material mode: SDV9 BODY = 1 sphere, 2 NbC, 3 Al2024, for
C      every point, particles included (from the label of the element,
C      of the parent for a particle): the three bodies of one section
C      can be coloured apart in the viewer. props(23) = first label of
C      the sphere (0: no sphere distinction). Only with 9 SDVs.
       if (iMulti .eq. 1 .and. nstatev .ge. 9) then
          do km = 1, nblock
             lab = jElemNum(km)
             if (lab .gt. nSolid) then
                ip = lPar0 + (lab - nSolid - 1)/nPpd3
             else
                ip = lab
             end if
             body = three
             if (ip .ge. lN0 .and. ip .le. lN1) then
                body = two
             else if (lS0 .gt. 0 .and. ip .ge. lS0) then
                body = one
             end if
             stateNew(km,9) = body
          end do
       end if
       return
       end
C
C
      subroutine hvi_eos_mg(rho0, aj, e, c0, s, gam0, p)
C     Mie-Grueneisen EOS, linear Us-Up form (README; = Moore 2025
C     eq. 16 with S2 = S3 = a = 0). p > 0 in compression.
C     in:  rho0 reference density, aj = J = rho0/rho, e internal
C          energy per unit mass, c0, s (Us = c0 + s up), gam0 Gamma0
C     out: p pressure
      include 'vaba_param.inc'
      parameter (one = 1.d0, half = 0.5d0)
      eta = one - aj
      p = rho0*c0**2*eta*(one - half*gam0*eta)/(one - s*eta)**2
     1  + gam0*rho0*e
      return
      end
C
C
      subroutine hvi_jc_yield(ep, epdot, a, b, expn, cjc, edot0,
     1                        sy, dsde, dsdr)
C     Johnson-Cook yield stress without the thermal term (README):
C     sy = (A + B ep^n) (1 + C ln(epdot/edot0)), rate factor 1 when
C     epdot < edot0.
C     in:  ep equivalent plastic strain, epdot its rate, A, B, n, C,
C          edot0
C     out: sy yield stress, dsde = d sy/d ep, dsdr = d sy/d epdot
C          (the derivatives serve the Newton iteration of the return)
      include 'vaba_param.inc'
      parameter (zero = 0.d0, one = 1.d0)
C     strain hardening; ep^(n-1) is infinite at ep = 0, so the
C     derivative is left at 0 there
      if (ep .gt. zero) then
         hard = a + b*ep**expn
         dhard = expn*b*ep**(expn - one)
      else
         hard = a
         dhard = zero
      end if
C     strain-rate factor
      if (epdot .gt. edot0 .and. edot0 .gt. 0.d0) then
         rate = one + cjc*log(epdot/edot0)
         drate = cjc/epdot
      else
         rate = one
         drate = zero
      end if
      sy = hard*rate
      dsde = dhard*rate
      dsdr = hard*drate
      return
      end
C
C
      subroutine hvi_radial_return(sTr, gmod, ep, dt, a, b, expn,
     1                             cjc, edot0, sNew, dep, sy)
C     Radial return of the deviator onto the JC yield surface.
C     in:  sTr(6) trial deviator (shear components 4-6 = tensor
C          components), gmod shear modulus G, ep plastic strain at the
C          start of the increment, dt time increment, JC constants
C     out: sNew(6) new deviator, dep plastic strain increment,
C          sy yield stress at the end of the increment (0 if elastic)
C     The increment dep solves f(dep) = q_tr - 3G dep - sy = 0, with sy
C     taken at ep + dep and at the rate dep/dt (implicit in the rate,
C     as the native JC). Newton iteration kept inside the bracket
C     [xlo, xhi], where f changes sign; a step that leaves the bracket
C     is replaced by a bisection.
      include 'vaba_param.inc'
      dimension sTr(6), sNew(6)
      parameter (zero = 0.d0, one = 1.d0, half = 0.5d0, three = 3.d0,
     1           threeHalf = 1.5d0, tol = 1.d-6, maxIt = 50)
C     equivalent (Mises) trial stress, shear components counted twice
      q = sqrt(threeHalf*(sTr(1)**2 + sTr(2)**2 + sTr(3)**2
     1        + 2*(sTr(4)**2 + sTr(5)**2 + sTr(6)**2)))
C     static yield stress at the start of the increment
      call hvi_jc_yield(ep, zero, a, b, expn, cjc, edot0,
     1                  sy0, dsde, dsdr)
      if (q .le. sy0) then
C        elastic
         do i = 1, 6
            sNew(i) = sTr(i)
         end do
         dep = zero
         sy = zero
         return
      end if
C     plastic. Bracket: f(0) = q - sy0 > 0; at x = (q - sy0)/(3G)
C     f <= 0 because sy only grows with dep.
      g3 = three*gmod
      xlo = zero
      xhi = (q - sy0)/g3
      x = xhi
      do it = 1, maxIt
         call hvi_jc_yield(ep + x, x/dt, a, b, expn, cjc, edot0,
     1                     sy, dsde, dsdr)
         f = q - g3*x - sy
         if (abs(f) .le. tol*q) exit
         if (f .gt. zero) then
            xlo = x
         else
            xhi = x
         end if
         df = -g3 - dsde - dsdr/dt
         xNew = x - f/df
         if (xNew .le. xlo .or. xNew .ge. xhi) xNew = half*(xlo + xhi)
         x = xNew
      end do
      dep = x
C     scale the trial deviator back to the yield surface:
C     q_new = q - 3G dep (= sy at convergence)
      ratio = (q - g3*dep)/q
      do i = 1, 6
         sNew(i) = sTr(i)*ratio
      end do
      return
      end
C
C
      subroutine hvi_jc_damage(dOld, dep, triax, epdot, d1, d2, d3,
     1                         d4, d5, edot0, dNew)
C     Johnson-Cook damage (Moore 2025 eq. 14-15; same form as the
C     Abaqus JOHNSON COOK damage initiation, exp(-d3 eta)):
C     eps_f = [d1 + d2 exp(-d3 eta)] [1 + d4 ln(epdot/edot0)]
C             [1 + d5 T*],   D = sum dep/eps_f
C     in:  dOld damage at the start of the increment, dep plastic
C          strain increment, triax = eta = -p/q (> 0 in tension),
C          epdot plastic strain rate, d1-d5, edot0
C     out: dNew damage at the end of the increment
C     No thermal term: T* = 0 (the decks start below Troom), so d5 is
C     not used. Rate factor 1 below edot0, as in the yield stress.
      include 'vaba_param.inc'
      parameter (one = 1.d0)
      if (epdot .gt. edot0 .and. edot0 .gt. 0.d0) then
         rate = one + d4*log(epdot/edot0)
      else
         rate = one
      end if
C     the argument is capped: exp overflows in single precision above 88
      epsf = (d1 + d2*exp(min(-d3*triax, 50.d0)))*rate
      dNew = dOld + dep/epsf
      return
      end
C
C
      subroutine hvi_nbc(nblock, ndir, nshr, nstatev, nprops,
     1     jInfoArray, stepTime, totalTime, dtArray, props, strainInc,
     2     stressOld, stateOld, enerInternOld, enerInelasOld,
     3     stressNew, stateNew, enerInternNew, enerInelasNew,
     4     iMask, nSolidX)
C     NbC branch (F4): Johnson-Holmquist 2 with the polynomial EOS and
C     bulking, as the built-in ABQ_JH2 (manual example "High-velocity
C     impact of a ceramic target", JH-2 case). Constants in the order of
C     ABQ_JH2, 8 per line:
C       1 rho0  2 G  3 A  4 N  5 B  6 M  7 C  8 epsdot0
C       9 T  10 sigma_i,max*  11 sigma_f,max*  12 HEL  13 p_HEL  14 beta
C      17 D1  18 D2  19 eps_f,max  20 eps_f,min  21 FS  22 lDamage
C      25 K1  26 K2  27 K3
C      33 Nsolid (optional: label above it = particle)
C     Normalisation: P* = p/p_HEL, T* = T/p_HEL, sigma* = sigma/sigma_HEL
C     with sigma_HEL = 1.5 (HEL - p_HEL) (Johnson and Holmquist 1994).
C     State variables: 1 EPSP, 2 DMG (D), 3 EPSV, 4 CONV, 5 SIG1, 6 PRESS,
C     7 DPBULK (bulking pressure), 8 EPRATE (plastic strain rate).
C     Element: converts (CONV = 0) when EPSP >= FS, as ABQ_JH2 sets
C     MPSTATUS = 0. Particle: fractured
C     material, D = 1 (strength B (P*)^M, no tension).
      include 'vaba_param.inc'
      dimension jInfoArray(*), dtArray(*), props(nprops),
     1  strainInc(nblock,ndir+nshr), stressOld(nblock,ndir+nshr),
     2  stateOld(nblock,nstatev), enerInternOld(nblock),
     3  enerInelasOld(nblock), stressNew(nblock,ndir+nshr),
     4  stateNew(nblock,nstatev), enerInternNew(nblock),
     5  enerInelasNew(nblock)
      parameter (zero = 0.d0, one = 1.d0, two = 2.d0, half = 0.5d0,
     1  three = 3.d0, third = 1.d0/3.d0, six = 6.d0, tiny = 1.d-12)
      parameter (i_info_ElemNumStartLoc = 6)
      pointer (ptrjElemNum, jElemNum)
      dimension jElemNum(nblock)
      dimension eigVal(3*maxblk), sTr(6), sNew(6)
C     iMask(km) = 1: this point is NbC (all points with a NBC material,
C     the NbC ones of a multi-material HVI block); nSolidX > 0: number
C     of solid elements given by the caller (multi-material deck)
      dimension iMask(nblock)
C
      rho0 = props(1)
      gmod = props(2)
      aJh = props(3)
      expN = props(4)
      bJh = props(5)
      expM = props(6)
      cJh = props(7)
      edot0 = props(8)
      tMax = props(9)
      siMax = props(10)
      sfMax = props(11)
      hel = props(12)
      pHel = props(13)
      beta = props(14)
      d1 = props(17)
      d2 = props(18)
      efMax = props(19)
      efMin = props(20)
      fs = props(21)
      aLdam = props(22)
      ak1 = props(25)
      ak2 = props(26)
      ak3 = props(27)
C     no particle unless Nsolid is given
      nSolid = 2000000000
      if (nprops .ge. 33) nSolid = nint(props(33))
      if (nSolidX .gt. 0) nSolid = nSolidX
      sHel = 1.5d0*(hel - pHel)
      twomu = two*gmod
      dt = dtArray(1)
      iElemNumStartLoc = jInfoArray(i_info_ElemNumStartLoc)
      ptrjElemNum = loc(jInfoArray(iElemNumStartLoc))
C
C     packaging call: linear elastic with K1 and G
      if (totalTime .eq. zero .and. stepTime .eq. zero) then
         alamda = ak1 - twomu*third
         do km = 1, nblock
            if (iMask(km) .eq. 0) cycle
            trace = strainInc(km,1) + strainInc(km,2) + strainInc(km,3)
            do i = 1, 3
               stressNew(km,i) = stressOld(km,i) + alamda*trace
     1                         + twomu*strainInc(km,i)
            end do
            do i = 4, 6
               stressNew(km,i) = stressOld(km,i)
     1                         + twomu*strainInc(km,i)
            end do
            do i = 1, nstatev
               stateNew(km,i) = stateOld(km,i)
            end do
            enerInternNew(km) = enerInternOld(km)
            enerInelasNew(km) = enerInelasOld(km)
         end do
         return
      end if
C
      do km = 1, nblock
         if (iMask(km) .eq. 0) cycle
         trace = strainInc(km,1) + strainInc(km,2) + strainInc(km,3)
         lab = jElemNum(km)
         do i = 1, nstatev
            stateNew(km,i) = stateOld(km,i)
         end do
C        volume: J from the accumulated log strain, mu = rho/rho0 - 1
         epsV = stateOld(km,3) + trace
         stateNew(km,3) = epsV
         aj = exp(epsV)
         rho = rho0/aj
         amu = one/aj - one
C        damage at the start of the increment (a particle is fractured)
         dOld = stateOld(km,2)
         if (lab .gt. nSolid) dOld = one
         dpOld = zero
         if (nstatev .ge. 7) dpOld = stateOld(km,7)
C        plastic strain rate of the previous increment (SDV8): the rate
C        factor is explicit, as ABQ_JH2 (its SDV2 PEEQ_RATE; in the compression test
C        with an implicit rate the onset of yielding overshot)
         epdOld = zero
         if (nstatev .ge. 8) epdOld = stateOld(km,8)
C        pressure: polynomial in compression, linear in tension, plus
C        the bulking pressure; tension limited to T (1 - D)
         if (amu .ge. zero) then
            p = ak1*amu + ak2*amu**2 + ak3*amu**3 + dpOld
         else
            p = ak1*amu + dpOld
         end if
C        tensile limit: -T as ABQ_JH2 with lDamage = 0 (test nbc_tens of
C        p = -192 MPa at D = 1); -T (1 - D) with lDamage = 1 (the
C        original JH-2 and LS-DYNA); a particle carries no tension
         if (aLdam .gt. half) then
            pMin = -tMax*(one - dOld)
         else
            pMin = -tMax
         end if
         if (lab .gt. nSolid) pMin = zero
         if (p .lt. pMin) p = pMin
C        trial deviator and JH-2 return (strength at the start-of-
C        increment damage, rate implicit)
         pOld = -third*(stressOld(km,1) + stressOld(km,2)
     1                + stressOld(km,3))
         do i = 1, 3
            sTr(i) = stressOld(km,i) + pOld
     1             + twomu*(strainInc(km,i) - third*trace)
         end do
         do i = 4, 6
            sTr(i) = stressOld(km,i) + twomu*strainInc(km,i)
         end do
         call hvi_jh2_return(sTr, gmod, p, dOld, epdOld, aJh, expN,
     1        bJh, expM, cJh, edot0, tMax, siMax, sfMax, pHel, sHel,
     2        sNew, dep, sy)
         stateNew(km,1) = stateOld(km,1) + dep
         if (nstatev .ge. 8) stateNew(km,8) = dep/dt
C        damage D = sum dep / eps_f, eps_f = D1 (P* + T*)^D2 (elements)
         dNew = dOld
         if (lab .le. nSolid .and. dep .gt. zero) then
            pts = (p + tMax)/pHel
            if (pts .gt. zero) then
               ef = d1*pts**d2
            else
               ef = efMin
            end if
            if (ef .lt. efMin) ef = efMin
            if (efMax .gt. zero .and. ef .gt. efMax) ef = efMax
            if (ef .gt. tiny) then
               dNew = min(one, dOld + dep/ef)
            else
               dNew = one
            end if
            if (dNew .gt. dOld) then
               call hvi_jh2_yield(p, epdOld, dOld, aJh, expN, bJh,
     1              expM, cJh, edot0, tMax, siMax, sfMax, pHel, sHel,
     2              syO, dum)
               call hvi_jh2_yield(p, epdOld, dNew, aJh, expN, bJh,
     1              expM, cJh, edot0, tMax, siMax, sfMax, pHel, sHel,
     2              syN, dum)
C              softening: the deviator is brought onto the damaged
C              surface WITHOUT plastic strain. Counting the drop as
C              plastic strain fed the damage back (D 0.05 -> 1 in three
C              increments, EPSP +0.007); ABQ_JH2 damages gradually and
C              adds only 0.0015 of PEEQ in the same drop (compression test)
               qN = sqrt(1.5d0*(sNew(1)**2 + sNew(2)**2 + sNew(3)**2
     1            + two*(sNew(4)**2 + sNew(5)**2 + sNew(6)**2)))
               if (qN .gt. syN .and. qN .gt. tiny) then
                  do i = 1, 6
                     sNew(i) = sNew(i)*max(syN, zero)/qN
                  end do
               end if
C              bulking: the drop of elastic deviatoric energy
C              U = sy^2/(6G) becomes hydrostatic potential energy
C              (fraction beta)
               dU = (syO**2 - syN**2)/(six*gmod)
               if (beta .gt. zero .and. dU .gt. zero .and.
     1             nstatev .ge. 7) then
                  stateNew(km,7) = -ak1*amu + sqrt((ak1*amu + dpOld)**2
     1                           + two*beta*ak1*dU)
               end if
            end if
         end if
         if (lab .le. nSolid) stateNew(km,2) = dNew
         do i = 1, 3
            stressNew(km,i) = sNew(i) - p
         end do
         do i = 4, 6
            stressNew(km,i) = sNew(i)
         end do
         stateNew(km,6) = p
C        energies per unit mass
         stressPower = half*(
     1     (stressOld(km,1)+stressNew(km,1))*strainInc(km,1)
     2   + (stressOld(km,2)+stressNew(km,2))*strainInc(km,2)
     3   + (stressOld(km,3)+stressNew(km,3))*strainInc(km,3)
     4   + two*(stressOld(km,4)+stressNew(km,4))*strainInc(km,4)
     5   + two*(stressOld(km,5)+stressNew(km,5))*strainInc(km,5)
     6   + two*(stressOld(km,6)+stressNew(km,6))*strainInc(km,6) )
         enerInternNew(km) = enerInternOld(km) + stressPower/rho
         enerInelasNew(km) = enerInelasOld(km) + sy*dep/rho
      end do
C
C     conversion flag (elements only, one-way): EPSP >= FS, as the
C     MPSTATUS of ABQ_JH2; largest principal stress for the output
C     rows of points that are not NbC are written later by the Al
C     branch (multi-material mode): zero them so that vsprinc never
C     sees undefined numbers
      do km = 1, nblock
         if (iMask(km) .eq. 0) then
            do i = 1, 6
               stressNew(km,i) = zero
            end do
         end if
      end do
      call vsprinc(nblock, stressNew, eigVal, ndir, nshr)
      do km = 1, nblock
         if (iMask(km) .eq. 0) cycle
         stateNew(km,5) = max(eigVal(km), eigVal(km+nblock),
     1                        eigVal(km+2*nblock))
         lab = jElemNum(km)
         if (lab .gt. nSolid) then
            stateNew(km,4) = one
         else if (stateOld(km,4) .eq. zero) then
            stateNew(km,4) = zero
         else if (fs .gt. zero .and. stateNew(km,1) .ge. fs) then
            stateNew(km,4) = zero
         else if (fs .lt. zero .and. stateNew(km,2) .ge. one) then
C           FS < 0 (our convention, not ABQ_JH2's): convert at fracture,
C           D = 1, whatever the plastic strain (F5 sweep)
            stateNew(km,4) = zero
         else
            stateNew(km,4) = one
         end if
      end do
      return
      end
C
C
      subroutine hvi_jh2_yield(p, epdot, dmg, aJh, expN, bJh, expM,
     1     cJh, edot0, tMax, siMax, sfMax, pHel, sHel, sy, dsdr)
C     JH-2 strength: sigma* = sigma_i* - D (sigma_i* - sigma_f*),
C       sigma_i* = min(A (P* + T*)^N, sigma_i,max*) (1 + C ln epsdot*)
C       sigma_f* = min(B (P*)^M,      sigma_f,max*) (1 + C ln epsdot*)
C     (rate factor 1 below epsdot0; a cap <= 0 means no cap).
C     out: sy = sigma_HEL sigma*, dsdr = d sy / d epsdot
      include 'vaba_param.inc'
      parameter (zero = 0.d0, one = 1.d0)
      ps = p/pHel
      ts = tMax/pHel
      if (epdot .gt. edot0 .and. edot0 .gt. 0.d0) then
         rate = one + cJh*log(epdot/edot0)
         drate = cJh/epdot
      else
         rate = one
         drate = zero
      end if
C     The caps act on the quasi-static strength, the rate factor
C     multiplies the capped value: compression test against
C     ABQ_JH2, fractured plateau 759 = 704 (cap) x 1.078 (rate factor at
C     the plastic strain rate of the run, 7e4 1/s)
C     intact
      si = zero
      if (ps + ts .gt. zero) si = aJh*(ps + ts)**expN
      if (siMax .gt. zero .and. si .gt. siMax) si = siMax
      siR = si*rate
      dsi = si*drate
C     fractured
      sf = zero
      if (ps .gt. zero) sf = bJh*ps**expM
      if (sfMax .gt. zero .and. sf .gt. sfMax) sf = sfMax
      sfR = sf*rate
      dsf = sf*drate
      sy = sHel*(siR - dmg*(siR - sfR))
      dsdr = sHel*(dsi*(one - dmg) + dmg*dsf)
      return
      end
C
C
      subroutine hvi_jh2_return(sTr, gmod, p, dmg, epdot, aJh, expN,
     1     bJh, expM, cJh, edot0, tMax, siMax, sfMax, pHel, sHel, sNew,
     2     dep, sy)
C     Radial return of the trial deviator onto the JH-2 surface at the
C     pressure p, damage dmg and plastic strain rate epdot of the
C     previous increment (explicit, as ABQ_JH2): closed form,
C     dep = (q - sy)/(3G). No strength (sy = 0, e.g. full tension): the
C     deviator vanishes. sy is returned also when elastic (energy uses
C     it only with dep > 0).
      include 'vaba_param.inc'
      dimension sTr(6), sNew(6)
      parameter (zero = 0.d0, three = 3.d0, threeHalf = 1.5d0)
      q = sqrt(threeHalf*(sTr(1)**2 + sTr(2)**2 + sTr(3)**2
     1        + 2*(sTr(4)**2 + sTr(5)**2 + sTr(6)**2)))
      call hvi_jh2_yield(p, epdot, dmg, aJh, expN, bJh, expM, cJh,
     1     edot0, tMax, siMax, sfMax, pHel, sHel, sy, dsdr)
      if (q .le. sy) then
         do i = 1, 6
            sNew(i) = sTr(i)
         end do
         dep = zero
      else if (sy .le. zero) then
         do i = 1, 6
            sNew(i) = zero
         end do
         dep = q/(three*gmod)
         sy = zero
      else
         dep = (q - sy)/(three*gmod)
         do i = 1, 6
            sNew(i) = sTr(i)*(sy/q)
         end do
      end if
      return
      end
