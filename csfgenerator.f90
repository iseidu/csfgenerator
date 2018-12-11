!IS:
! Script for constructing csf list from DFT/MRCI output (out3 or cvsout3).
! The CI ordering of MOs are used and the MOs in mos are reordered accordingly.
! Files needed for running this script:
!   ---> out3 (or cvsout3)
!   ---> mos
! Outputs:
!   ---> csf list named as 'csf'//root number (e.g. csf1, csf2 ...)
!   ---> mos_new containing the MO ordering consistent with the CI ordering
!

      program csfgenerator

      implicit none

      integer  :: iext,nspace,nocc,iint,ierr,i,j,k,ncount,count1,ici,itm
      integer  :: nroots,int1,int2,int3,int4,int5,int6,int7,int8,iorb
      integer  :: nsaos,nchunk,nleft,nchk,ichk,ispace,socc,icount,l,&
                  itmp,msaos,norb,nel,docc,ifreeze,nfreeze
      logical  :: lexist,lxst
      integer, parameter :: dp = selected_real_kind(15, 307)
      character :: string1*20,string*40,char1*24,char2,char3*7,char4*6
      character :: char5*8,outfile*6 !outfile1*4,outfile2*5,outfile3*6
      character :: char6*22,char16*26,fmt1*80,fmt2*80,fmt3*80,fmt4*80
      character :: char7*7,char8*24,charr*55,str1*46,str2*46,str3,&
                   char0*26,string16*6,string17*30,string18*25,&
                   string19*16,fmt0*80,char00*7,char01*15
      character :: char9*5,char10*34,fmt5*80,char21*9,char22*2,label*2,&
                   char33*2,arg*32
      integer, allocatable :: iocc(:),ntm(:),nbocc(:,:),nci(:),&
                              tmpocc(:),isocc(:),init_ci(:),init_tm(:),&
                              ntm_f(:),iocc_f(:)
      real(dp)      :: coef
      real(dp), allocatable    :: enval(:),mo2as(:,:),cartas(:,:)
      logical                  :: lreadci,lfreeze

!c     Check the number of D + V in the dft/mrci output and 
!c     construct the base occupation number vector by searching
!c     'lowest external mo'.
!c     Occupied:   3
!c     Unoccupied: 0
!c     Singly occupied: 1 or 2

      lreadci = .false.
      lfreeze = .false.
      ifreeze = 0
      nfreeze = 0

      call get_command_argument(1,arg)
      if (len_trim(arg) == 0) then
         inquire(file='cvsout3',exist=lexist)
         if (lexist) then
            open(47,file='cvsout3',form='formatted',status='old')
         else
            inquire(file='out3',exist=lxst)
            if (lxst) then
               open(47,file='out3',form='formatted',status='old')
            else
               write(6,*) 'DFT/MRCI output file (cvsout3 or out3) not &
                           found!'
               stop
            endif
         endif
      else
         write(6,*)'Command line DFT/MRCI output file: ',trim(arg)
         open(47,file=arg,form='formatted',status='old')
2111     read(47,'(a30)',end=3111)string17
         if (string17.eq.' generating reference ci-space') then
            read(47,*)
            read(47,'(35x,i3)')nel
            read(47,'(35x,i3)')norb
            read(47,*)
            allocate(init_ci(norb),init_tm(norb),stat=ierr)
            if(ierr.ne.0) stop 'Memory allocation failed!'
            write(fmt0,'(a,i2,a)')'(a25,',norb,'(i4))'
            read(47,fmt0)string18,(init_ci(i),i=1,norb)
            read(47,fmt0)string18,(init_tm(i),i=1,norb)
4111        read(47,'(a16)')string19
            if (string19.eq.'# D internal mos') then
               backspace(47)
               read(47,'(18x,i3)')docc
            else
               goto 4111
            endif
            lreadci = .true.
            goto 3111
         else
            goto 2111
         endif
      endif

3111  continue
      rewind(47)

!     Check if calculation was run in c1, abort if it wasn't
      inquire(file='control',exist=lexist)
      if (lexist) then
         open(55,file='control',form='formatted',status='old')
989      read(55,'(a9)')char21
         if (char21.eq.'$symmetry') then
            backspace(55)
            read(55,'(a9,1x,a2)')char21,char22
         else
            goto 989
         endif
         rewind(55)
1009     read(55,'(a7)')char00
         if (char00.eq.'$freeze') then
            read(55,'(a15)')char01
            if (char01.ne.' implicit core=') then
               backspace(55)
               read(55,'(7x,i3,1x,i3)')ifreeze,nfreeze
               lfreeze = .true.
            else
               lfreeze = .false.
            endif
         else
            goto 1009
         endif
         goto 199
         close(55)
      else
         stop 'control file not found!'
      endif

199   continue
      if (char22.ne.'c1') stop 'Calculation has to be in c1 symmetry!'

!     Get CAS information
110   read(47,'(a20)',end=999)string1
      if (string1.ne.'  lowest external mo') goto 110

999   continue
!      write(6,*)string1

      backspace(47)
      read(47,'(a40,i3)')string,iext

      if (lreadci) then
         nspace = docc + norb
      else
         nspace = iext - 1
      endif

!     For allocating ntm
1110  read(47,'(a6)',end=1999)string16
      if (string16.ne.'frozen') goto 1110
1999  continue
      backspace(47)
      read(47,'(7x,i3)')ispace
      !write(6,*)ispace

!      write(6,*)nspace

!      do i = 1, 3
!         read(47,*)
!      enddo

!     Get the number of doubly occupied and singly occupied 
!     orbitals
      rewind(47)

      !allocate(iocc(nspace),isocc(nspace),stat=ierr)
      allocate(iocc(ispace),isocc(ispace),stat=ierr)
      if(ierr.ne.0) stop 'Memory allocation failed!'

      nocc = 0
      socc = 0
      icount = 0
      iocc = 0
      isocc = 0
910   read(47,'(a7)',end=111)char3
      if (char3.eq.' typ nr') then
         goto 111
      else
         goto 910
      endif

111   continue
!      write(6,*)char3

210   read(47,'(a24,i1)',end=777)char1,iint
!      write(6,*)iint
!     doubly occupied
      if (iint.eq.2) then
         nocc = nocc + 1
         icount = icount + 1
         iocc(icount) = 3
         goto 210
!     singly occupied
      else if (iint.eq.1) then
         socc = socc + 1
         icount = icount + 1
         iocc(icount) = 1
         isocc(socc) = icount
         goto 210
      else
         goto 777
      endif

777   continue


!c    Create the base occupation vector
!     Now incorporated into the previous
!      iocc = 0
!      do i = 1, nocc
!         iocc(i) = 3
!      enddo

!      write(6,*)nocc,iocc

!c     ci/tm orbital ordering for re-ordering mos
!c     itm contains the mo indices to be used
      rewind(47)
310   read(47,'(a7)',end=888)char3
      if (char3.eq.' typ nr') then
         goto 888
      else
         goto 310
      endif

888   continue

!      ispace = 3 * nspace
      allocate(ntm(ispace),stat=ierr)
      if(ierr.ne.0) stop 'Memory allocation failed!'

      ntm = 0
      ncount = 0
410   read(47,'(a7)',end=877)char3
!      write(6,*)char3
      if (char3.ne.'') then
         backspace(47)
         read(47,'(a6,2(1x,i3))')char4,ici,itm
         ncount = ncount + 1
         ntm(ncount) = itm
         goto 410
      else
         goto 877
      endif

877   continue
      !write(6,*)'Got here'
      !write(6,*)ispace,docc,norb
      !write(6,*)init_tm
      !write(6,*)lreadci
!     Correct reference space according CI read-in
      if (lreadci) then
         !write(6,*)ntm
         ntm(docc+1:docc+norb) = init_tm(1:norb)
      endif

      !write(6,*)'Got here'
      !write(6,*)lfreeze
      !write(6,*)nfreeze
!     If there are frozen core orbitals include them
      !if (nfreeze.ne.0) then
      if (lfreeze) then
         ispace = ispace + nfreeze
         allocate(ntm_f(ispace),iocc_f(ispace),stat=ierr)
         if (ierr.ne.0) stop 'Memory allocation failed!'
         ntm_f = 0
         do i = 1, nfreeze
            ntm_f(i) = i
            iocc_f(i) = 3
         enddo
         ncount = ncount + nfreeze
         nspace = nspace + nfreeze
         iocc_f(nfreeze+1:ispace) = iocc(1:ispace-nfreeze)
         do i = 1, ispace-nfreeze
            ntm_f(i+nfreeze) = nfreeze + ntm(i)
         enddo
         deallocate(iocc,ntm)
         allocate(ntm(ispace),iocc(ispace),stat=ierr)
         if (ierr.ne.0) stop 'Memory allocation failed!'
         iocc = iocc_f
         ntm = ntm_f
      endif
      !write(6,'(4i3)')ntm_f
      !write(6,*)
      !write(6,'(4i3)')ntm
      !write(6,*)
      !write(6,'(4i3)')iocc_f
      !write(6,*)
      !write(6,'(4i3)')iocc

!c     Search 'all roots converged' for the number of roots
!c     since it might differ from the number requested.
!510   read(47,'(a19)',end=866)char5
!      if (char5.ne.' loop partitioning:') then
!         goto 510
!      else
!         goto 866
!      endif
510   read(47,'(a8)',end=866)char5
      if (char5.ne.' nrprint') then
         goto 510
      else
         goto 866
      endif
866   continue
      backspace(47)
      read(47,'(17x,i3)')nroots
!      write(6,*)nroots

      !allocate(nbocc(nroots,nspace),nci(ispace),tmpocc(nspace),&
      allocate(nbocc(nroots,ispace),nci(ispace),tmpocc(ispace),&
               stat=ierr)
      if(ierr.ne.0) stop 'Memory allocation failed'

      nbocc = 0
      nci = 0
      tmpocc = 0

!c      do i = 1, nroots
!c         nbocc(i,:) = iocc(:)
!c      enddo
!      nbocc(1,:) = iocc(:)

      !write(6,*)ncount
      do i = 1, ncount
         nci(ntm(i)) = i
      enddo
!      write(6,*)init_tm
!      write(6,*)ntm
!      write(6,*)nci

!    Write format for output csfs
     write(fmt1,'(a,i2,a,a)')'(f10.7,',nspace,'(1x,i1),',&
                             '4(2x,a2,i2))'
     write(fmt2,'(a,i2,a,a,a)')'(f10.7,',nspace,'(1x,i1),',&
                             '1x,i2,a1,i2,','3(2x,a2,i2))'
     write(fmt3,'(a,i2,a,a,a)')'(f10.7,',nspace,'(1x,i1),',&
                             '2(1x,i2,a1,i2),','2(2x,a2,i2))'
     write(fmt4,'(a,i2,a,a,a)')'(f10.7,',nspace,'(1x,i1),',&
                             '3(1x,i2,a1,i2),','2x,a2,i2)'
     write(fmt5,'(a,i2,a,a,a)')'(f10.7,',nspace,'(1x,i1),',&
                             '4(1x,i2,a1,i2))'

     open(333,file='nspace',status='unknown')
     write(333,*)nspace
     close(333)
!     write(6,*)fmt2
!     write(6,*)fmt3
!     write(6,*)fmt4
!     write(6,*)fmt5

!     Contruct csfs from out3 (or cvsout3) for csf2det program
      do i = 1, nroots
         if (i.lt.10) then
            write(outfile,'(a3,i1)')'csf',i
!            open(102,outfile,status='unknown')
         else if (i.lt.100) then
            write(outfile,'(a3,i2)')'csf',i
         else
            write(outfile,'(a3,i3)')'csf',i
         endif

!         write(6,*)outfile

         open(102,file=outfile,status='unknown')

610      read(47,'(a22)',end=855)char6
         if (char6.ne.'  state         method') then
            goto 610
         else
            goto 855
         endif
855      continue
         do j = 1, 5
            read(47,*)
         enddo
         
         count1 = 0
710      read(47,'(2x,a7,a24,a2,a2,i1)')char7,char8,label,&
                                                char33,iorb
         tmpocc = iocc
         if (char7.ne.'') then
!           Debug write
            !write(6,*)label,iorb
            backspace(47)
            count1 = count1 + 1
            if ((label.eq.'GS'.and.iorb.eq.0).or.(label.eq.'GS'.and.&
                 iorb.eq.1)) then
               read(47,'(a5,f10.7)')char9,coef
               write(102,fmt1)coef,&
                       (tmpocc(k),k=1,nspace),'0:',0,'0:',0,'0:',0,&
                       '0:',0
!               write(6,*)iorb,coef
            else if (label.eq.'S'.and.iorb.eq.1) then
               read(47,'(a5,f10.7,a34,2(1x,i3))')char9,coef,char10,&
                                                int1,int2
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
               endif
               tmpocc(nci(int1)) = 1
               tmpocc(nci(int2)) = 3
               write(102,fmt1)coef,&
                    (tmpocc(k),k=1,nspace),'0:',0,'0:',0,'0:',0,&
                    '0:',0
            else if ((label.eq.'S'.and.iorb.eq.2).or.(label.eq.&
                       'S'.and.iorb.eq.3)) then
               read(47,'(a5,f10.7,a34,2(1x,i3))')char9,coef,char10,&
                                                 int1,int2
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
               endif
!c               ncoef(count1) = coef
               if (nci(int2).le.nspace) then
                  tmpocc(nci(int1)) = 1
                  tmpocc(nci(int2)) = 2
                  !write(102,'(f10.7,nspace(1x,i1),4(1x,i2))')coef,&
                  write(102,fmt1)coef,&
                       (tmpocc(k),k=1,nspace),'0:',0,'0:',0,'0:',0,&
                       '0:',0
               else
                  tmpocc(nci(int1)) = 1
                  write(102,fmt2)coef,&
                       (tmpocc(k),k=1,nspace),nci(int2),':',2,'0:',0,&
                       '0:',0,'0:',0
               endif
            else if (label.eq.'D'.and.iorb.eq.0) then
               read(47,'(a5,f10.7,a34,4(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
               endif
               if (nci(int3).le.nspace.and.nci(int4).le.nspace) then
                  if (nci(int3).eq.nci(int4)) then
                     if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 3
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 2
                       tmpocc(nci(int4)) = 3
                     endif
                  else
                     if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 2
                     endif
                  endif
                  write(102,fmt1)coef,&
                       (tmpocc(k),k=1,nspace),'0:',0,'0:',0,'0:',0,&
                       '0:',0
               !else if (nci(int3).le.nspace.and.nci(int4).gt.nspace) then
               !   if (nci(int1).eq.nci(int2)) then
               !     tmpocc(nci(int2)) = 0
               !     tmpocc(nci(int3)) = 1
               !   else
               !     tmpocc(nci(int1)) = 1
               !     tmpocc(nci(int2)) = 1
               !     tmpocc(nci(int3)) = 2
               !   endif
               !   write(102,fmt2)coef,&
               !        (tmpocc(k),k=1,nspace),nci(int4),':',2,'0:',0,&
               !        '0:',0,'0:',0
               else if (nci(int3).gt.nspace.and.nci(int4).gt.nspace) then
                  if (nci(int1).eq.nci(int2)) then
                    tmpocc(nci(int2)) = 0
                  else
                    tmpocc(nci(int1)) = 1
                    tmpocc(nci(int2)) = 2
                  endif
                  if (nci(int3).eq.nci(int4)) then
                     write(102,fmt2)coef,&
                       (tmpocc(k),k=1,nspace),nci(int3),':',3,&
                       '0:',0,'0:',0,'0:',0
                  else
                     write(102,fmt3)coef,&
                       (tmpocc(k),k=1,nspace),nci(int3),':',1,&
                       nci(int4),':',2,'0:',0,'0:',0
                  endif
               endif
            else if ((label.eq.'D'.and.iorb.eq.4).or.(label.eq.&
                      'D'.and.iorb.eq.2).or.(label.eq.'D'.and.&
                       iorb.eq.5)) then
               read(47,'(a5,f10.7,a34,4(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
               endif
               if (nci(int3).le.nspace.and.nci(int4).le.nspace) then
                  if (nci(int3).eq.nci(int4)) then
                     if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 3
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 2
                       tmpocc(nci(int4)) = 3
                     endif
                  else
                     if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 2
                     endif
                  endif
                  write(102,fmt1)coef,&
                       (tmpocc(k),k=1,nspace),'0:',0,'0:',0,'0:',0,&
                       '0:',0
               else if (nci(int3).le.nspace.and.nci(int4).gt.nspace) then
                  if (nci(int1).eq.nci(int2)) then
                    tmpocc(nci(int2)) = 0
                    tmpocc(nci(int3)) = 1
                  else
                    tmpocc(nci(int1)) = 1
                    tmpocc(nci(int2)) = 1
                    tmpocc(nci(int3)) = 2
                  endif
                  write(102,fmt2)coef,&
                       (tmpocc(k),k=1,nspace),nci(int4),':',2,'0:',0,&
                       '0:',0,'0:',0
               else if (nci(int3).gt.nspace.and.nci(int4).gt.nspace) then
                  if (nci(int1).eq.nci(int2)) then
                    tmpocc(nci(int2)) = 0
                  else
                    tmpocc(nci(int1)) = 1
                    tmpocc(nci(int2)) = 2
                  endif
                  if (nci(int3).eq.nci(int4)) then
                     write(102,fmt2)coef,&
                       (tmpocc(k),k=1,nspace),nci(int3),':',3,&
                       '0:',0,'0:',0,'0:',0
                  else
                     write(102,fmt3)coef,&
                       (tmpocc(k),k=1,nspace),nci(int3),':',1,&
                       nci(int4),':',2,'0:',0,'0:',0
                  endif
               endif
            else if (label.eq.'D'.and.iorb.eq.5) then
               read(47,'(a5,f10.7,a34,4(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
               endif
               tmpocc(nci(int1)) = 1
               tmpocc(nci(int2)) = 1
               tmpocc(nci(int3)) = 2
               tmpocc(nci(int4)) = 2
               write(102,fmt1)coef,&
                       (tmpocc(k),k=1,nspace),'0:',0,'0:',0,'0:',0,&
                       '0:',0
            else if (label.eq.'D'.and.iorb.eq.3) then
               read(47,'(a5,f10.7,a34,4(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
               endif
               do l = 1, socc
                  if (int1.eq.isocc(l)) then
                     tmpocc(nci(int1)) = 0
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 2
                     tmpocc(nci(int4)) = 2
                  else if (int2.eq.isocc(l)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 2
                     tmpocc(nci(int4)) = 2
                  else if (int3.eq.isocc(l)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 3
                     tmpocc(nci(int4)) = 2
                  else if (int4.eq.isocc(l)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 2
                     tmpocc(nci(int4)) = 3
                  endif
               enddo
               write(102,fmt1)coef,&
                       (tmpocc(k),k=1,nspace),'0:',0,'0:',0,'0:',0,&
                       '0:',0
            else if (label.eq.'D'.and.iorb.eq.1) then
               read(47,'(a5,f10.7,a34,4(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
               endif
               do l = 1, socc
                  if (int1.eq.isocc(l)) then
                     tmpocc(nci(int1)) = 0
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int4)) = 3
                  else if (int2.eq.isocc(l)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int4)) = 3
                  else if (int3.eq.isocc(l)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 3
                     tmpocc(nci(int4)) = 2
                  else if (int4.eq.isocc(l)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 2
                     tmpocc(nci(int4)) = 3
                  endif
               enddo
               write(102,fmt1)coef,&
                       (tmpocc(k),k=1,nspace),'0:',0,'0:',0,'0:',0,&
                       '0:',0
            else if (label.eq.'T'.and.iorb.eq.2) then
               read(47,'(a5,f10.7,a34,6(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4,&
                                                 int5,int6
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
                  int5 = int5 + nfreeze
                  int6 = int6 + nfreeze
               endif
               if (nci(int4).le.nspace.and.nci(int5).le.nspace.and.&
                            nci(int6).le.nspace) then
                  if (nci(int1).eq.nci(int2).and.nci(int4).eq.&
                          nci(int5)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int5)) = 3
                     tmpocc(nci(int6)) = 2
                  else if (nci(int1).eq.nci(int2).and.nci(int5).eq.&
                          nci(int6)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 2
                     tmpocc(nci(int6)) = 3
                  else if (nci(int1).eq.nci(int2).and.nci(int4).eq.&
                          nci(int6)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 3
                  else if (nci(int2).eq.nci(int3).and.nci(int4).eq.&
                          nci(int5)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int5)) = 3
                     tmpocc(nci(int6)) = 2
                  else if (nci(int2).eq.nci(int3).and.nci(int5).eq.&
                          nci(int6)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 2
                     tmpocc(nci(int6)) = 3
                  else if (nci(int2).eq.nci(int3).and.nci(int4).eq.&
                          nci(int6)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 3
                  endif
                  write(102,fmt1)coef,&
                       (tmpocc(k),k=1,nspace),'0:',0,'0:',0,'0:',0,&
                       '0:',0
               else if (nci(int4).le.nspace.and.nci(int5).le.&
                    nspace.and.nci(int6).gt.nspace) then
                  if (nci(int1).eq.nci(int2).and.nci(int4).eq.&
                         nci(int5)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int5)) = 3
                  else if (nci(int2).eq.nci(int3).and.nci(int4).eq.&
                         nci(int5)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int5)) = 3
                  endif
                  write(102,fmt2)coef,&
                       (tmpocc(k),k=1,nspace),nci(int6),':',2,'0:',0,&
                       '0:',0,'0:',0
               else if (nci(int4).le.nspace.and.nci(int5).gt.&
                         nspace.and.nci(int6).gt.nspace) then
                  if (nci(int1).eq.nci(int2)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 2
                  else if (nci(int2).eq.nci(int3)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 2
                  endif
                  write(102,fmt2)coef,&
                    (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                     '0:',0,'0:',0,&
                    '0:',0
               else if (nci(int4).gt.nspace.and.nci(int5).gt.&
                        nspace.and.nci(int6).gt.nspace) then
                  if (nci(int1).eq.nci(int2)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                  else if (nci(int2).eq.nci(int3)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                  endif
                  if (nci(int4).eq.nci(int5)) then
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                  else if (nci(int5).eq.nci(int6)) then
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                  else if (nci(int4).eq.nci(int6)) then
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                  endif
               endif
            else if (label.eq.'T'.and.iorb.eq.1) then
               read(47,'(a5,f10.7,a34,6(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4,&
                                                 int5,int6
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
                  int5 = int5 + nfreeze
                  int6 = int6 + nfreeze
               endif
               if (nci(int4).le.nspace.and.nci(int5).le.nspace.and.&
                            nci(int6).le.nspace) then
                  do l = 1, socc
                     if (int1.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int5)) = 3
                           tmpocc(nci(int6)) = 2
                        else if (int5.eq.int6) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 2
                           tmpocc(nci(int6)) = 3
                        endif
                     else if (int2.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int5)) = 3
                           tmpocc(nci(int6)) = 2
                        else if (int5.eq.int6) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 2
                           tmpocc(nci(int6)) = 3
                        endif
                     else if (int3.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int5)) = 3
                           tmpocc(nci(int6)) = 2
                        else if (int5.eq.int6) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 2
                           tmpocc(nci(int6)) = 3
                        endif
                     else if (int4.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 3
                           tmpocc(nci(int6)) = 3
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 3
                           tmpocc(nci(int6)) = 3
                        endif
                     else if (int5.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int5)) = 3
                           tmpocc(nci(int6)) = 3
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 3
                           tmpocc(nci(int6)) = 3
                        endif
                     else if (int6.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 3
                           tmpocc(nci(int6)) = 3
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 3
                           tmpocc(nci(int6)) = 3
                        endif
                     endif
                  enddo
                  write(102,fmt1)coef,&
                          (tmpocc(k),k=1,nspace),'0:',0,&
                           '0:',0,'0:',0,&
                          '0:',0
               else if (nci(int4).le.nspace.and.nci(int5).le.&
                        nspace.and.nci(int6).gt.nspace) then
                  do l = 1, socc
                     if (int1.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int5)) = 3
                        endif
                        write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int2.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int5)) = 3
                        endif
                        write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int3.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int5)) = 3
                        endif
                        write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int4.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 3
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 3
                        endif
                        write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int5.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int5)) = 3
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 3
                        endif
                        write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int6.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 3
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 3
                        endif
                        write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                           '0:',0,'0:',0,&
                          '0:',0
                     endif
                  enddo
               else if (nci(int4).le.nspace.and.nci(int5).gt.&
                        nspace.and.nci(int6).gt.nspace) then
                  do l = 1, socc
                     if (int1.eq.isocc(l)) then
                        if (int5.eq.int6) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 2
                        endif
                        write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int2.eq.isocc(l)) then
                        if (int5.eq.int6) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 2
                        endif
                        write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int3.eq.isocc(l)) then
                        if (int5.eq.int6) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 2
                        endif
                        write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int4.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 3
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 3
                        endif
                        write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int5.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int5)) = 3
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 3
                        endif
                        write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int6.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 3
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 3
                        endif
                        write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                           '0:',0,'0:',0,&
                          '0:',0
                     endif
                  enddo
               else if (nci(int4).gt.nspace.and.nci(int5).gt.&
                            nspace.and.&
                            nci(int6).gt.nspace) then
                  do l = 1, socc
                     if (int1.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                        else if (int5.eq.int6) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                        endif
                     else if (int2.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                        else if (int5.eq.int6) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                        endif
                     else if (int3.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                        else if (int5.eq.int6) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                        endif
                     else if (int4.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',3,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',3,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                        endif
                     else if (int5.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',3,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                        endif
                     else if (int6.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',3,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',3,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                        endif
                     endif
                  enddo
               endif
            else if (label.eq.'T'.and.iorb.eq.3) then
               read(47,'(a5,f10.7,a34,6(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4,&
                                                 int5,int6
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
                  int5 = int5 + nfreeze
                  int6 = int6 + nfreeze
               endif
               if (nci(int4).le.nspace.and.nci(int5).le.nspace.and.&
                            nci(int6).le.nspace) then
                  do l = 1, socc
                     if (int1.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int5)) = 3
                           tmpocc(nci(int6)) = 2
                        else if (int5.eq.int6) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 1
                           tmpocc(nci(int6)) = 3
                        else
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 1
                           tmpocc(nci(int5)) = 2
                           tmpocc(nci(int6)) = 2
                        endif
                     else if (int2.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int5)) = 3
                           tmpocc(nci(int6)) = 2
                        else if (int5.eq.int6) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 1
                           tmpocc(nci(int6)) = 3
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 1
                           tmpocc(nci(int5)) = 2
                           tmpocc(nci(int6)) = 2
                        endif
                     else if (int3.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int5)) = 3
                           tmpocc(nci(int6)) = 2
                        else if (int5.eq.int6) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 1
                           tmpocc(nci(int6)) = 3
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 1
                           tmpocc(nci(int5)) = 2
                           tmpocc(nci(int6)) = 2
                        endif
                     else if (int4.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 3
                           tmpocc(nci(int5)) = 2
                           tmpocc(nci(int6)) = 2
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 3
                           tmpocc(nci(int5)) = 2
                           tmpocc(nci(int6)) = 2
                        endif
                     else if (int5.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 2
                           tmpocc(nci(int5)) = 3
                           tmpocc(nci(int6)) = 2
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 2
                           tmpocc(nci(int5)) = 3
                           tmpocc(nci(int6)) = 2
                        endif
                     else if (int6.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 2
                           tmpocc(nci(int5)) = 2
                           tmpocc(nci(int6)) = 3
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 2
                           tmpocc(nci(int5)) = 2
                           tmpocc(nci(int6)) = 3
                        endif
                     endif
                  enddo
                  write(102,fmt1)coef,&
                          (tmpocc(k),k=1,nspace),'0:',0,&
                           '0:',0,'0:',0,&
                          '0:',0
               else if (nci(int4).le.nspace.and.nci(int5).le.&
                        nspace.and.nci(int6).gt.nspace) then
                  do l = 1, socc
                     if (int1.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int5)) = 3
                        else
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 1
                           tmpocc(nci(int5)) = 2
                        endif
                           write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int2.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int5)) = 3
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 1
                           tmpocc(nci(int5)) = 2
                        endif
                           write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int3.eq.isocc(l)) then
                        if (int4.eq.int5) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int5)) = 3
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 1
                           tmpocc(nci(int5)) = 2
                        endif
                           write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                     else if (int4.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 3
                           tmpocc(nci(int5)) = 2
                           write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 3
                           tmpocc(nci(int5)) = 2
                           write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                        endif
                     else if (int5.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 2
                           tmpocc(nci(int5)) = 3
                           write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 2
                           tmpocc(nci(int5)) = 3
                           write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                        endif
                     else if (int6.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 2
                           tmpocc(nci(int5)) = 2
                           write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                           '0:',0,'0:',0,&
                          '0:',0
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 2
                           tmpocc(nci(int5)) = 2
                           write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                           '0:',0,'0:',0,&
                          '0:',0
                        endif
                     endif
                  enddo
               else if (nci(int4).le.nspace.and.nci(int5).gt.&
                        nspace.and.nci(int6).gt.nspace) then
                  do l = 1, socc
                     if (int1.eq.isocc(l)) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 1
                         if (int5.eq.int6) then
                           write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                          '0:',0,'0:',0,&
                          '0:',0
                         else
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                         endif
                     else if (int2.eq.isocc(l)) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 1
                         if (int5.eq.int6) then
                           write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                          '0:',0,'0:',0,&
                          '0:',0
                         else
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                         endif
                     else if (int3.eq.isocc(l)) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 1
                         if (int5.eq.int6) then
                           write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                          '0:',0,'0:',0,&
                          '0:',0
                         else
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                         endif
                     else if (int4.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 3
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 3
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                        endif
                     else if (int5.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 2
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 2
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                        endif
                     else if (int6.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           tmpocc(nci(int4)) = 2
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           tmpocc(nci(int4)) = 2
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                        endif
                     endif
                  enddo
               else if (nci(int4).gt.nspace.and.nci(int5).gt.&
                        nspace.and.nci(int6).gt.nspace) then
                  do l = 1, socc
                     if (int1.eq.isocc(l)) then
                           tmpocc(nci(int1)) = 0
                           tmpocc(nci(int3)) = 0
                        if (int4.eq.int5) then
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                        else if (int5.eq.int6) then
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',1,&
                           nci(int6),':',3,'0:',2,&
                          '0:',0
                        else
                           write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',1,&
                           nci(int5),':',2,nci(int6),':',2,&
                          '0:',0
                        endif
                     else if (int2.eq.isocc(l)) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                        if (int4.eq.int5) then
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                        else if (int5.eq.int6) then
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',1,&
                           nci(int6),':',3,'0:',2,&
                          '0:',0
                        else
                           write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',1,&
                           nci(int5),':',2,nci(int6),':',2,&
                          '0:',0
                        endif
                     else if (int3.eq.isocc(l)) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 0
                        if (int4.eq.int5) then
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                        else if (int5.eq.int6) then
                           write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',1,&
                           nci(int6),':',3,'0:',2,&
                          '0:',0
                        else
                           write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',1,&
                           nci(int5),':',2,nci(int6),':',2,&
                          '0:',0
                        endif
                     else if (int4.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',3,&
                           nci(int5),':',2,nci(int6),':',2,&
                          '0:',0
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',3,&
                           nci(int5),':',2,nci(int6),':',2,&
                          '0:',0
                        endif
                     else if (int5.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int5),':',3,nci(int6),':',2,&
                          '0:',0
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int5),':',3,nci(int6),':',2,&
                          '0:',0
                        endif
                     else if (int6.eq.isocc(l)) then
                        if (int1.eq.int2) then
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int3)) = 1
                           write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int5),':',2,nci(int6),':',3,&
                          '0:',0
                        else if (int2.eq.int3) then
                           tmpocc(nci(int1)) = 1
                           tmpocc(nci(int3)) = 0
                           write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int5),':',2,nci(int6),':',3,&
                          '0:',0
                        endif
                     endif
                  enddo
               endif
            else if (label.eq.'T'.and.iorb.eq.4) then
               read(47,'(a5,f10.7,a34,6(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4,&
                                                 int5,int6
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
                  int5 = int5 + nfreeze
                  int6 = int6 + nfreeze
               endif
               if (nci(int4).le.nspace.and.nci(int5).le.nspace.and.&
                            nci(int6).le.nspace) then
                  if (int1.eq.int2) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 2
                  else if (int2.eq.int3) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 2
                  else if (int4.eq.int5) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 2
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int5)) = 3
                     tmpocc(nci(int6)) = 2
                  else if (int5.eq.int6) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 2
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 2
                     tmpocc(nci(int6)) = 3
                  endif
                  write(102,fmt1)coef,&
                          (tmpocc(k),k=1,nspace),'0:',0,&
                           '0:',0,'0:',0,&
                          '0:',0
               else if (nci(int4).le.nspace.and.nci(int5).le.&
                        nspace.and.nci(int6).gt.nspace) then
                  if (int1.eq.int2) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 2
                     tmpocc(nci(int5)) = 1
                     write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                  else if (int2.eq.int3) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 2
                     tmpocc(nci(int5)) = 1
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                          '0:',0,'0:',0,&
                          '0:',0
                  else if (int4.eq.int5) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 2
                     tmpocc(nci(int5)) = 3
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                          '0:',0,'0:',0,&
                          '0:',0
                  endif
               else if (nci(int4).le.nspace.and.nci(int5).gt.&
                        nspace.and.nci(int6).gt.nspace) then
                  if (int1.eq.int2) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 1
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                  else if (int2.eq.int3) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 1
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                  else if (int5.eq.int6) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 2
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 2
                     write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                          '0:',0,'0:',0,&
                          '0:',0
                  endif
               else if (nci(int4).gt.nspace.and.nci(int5).gt.&
                        nspace.and.nci(int6).gt.nspace) then
                  if (int1.eq.int2) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',1,&
                           nci(int5),':',2,nci(int6),':',2,&
                          '0:',0
                  else if (int2.eq.int3) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',1,&
                           nci(int5),':',2,nci(int6),':',2,&
                          '0:',0
                  else if (int4.eq.int5) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 2
                     tmpocc(nci(int3)) = 1
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                  else if (int5.eq.int6) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 2
                     tmpocc(nci(int3)) = 1
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                  endif
               endif
            else if (label.eq.'T'.and.iorb.eq.5) then
               read(47,'(a5,f10.7,a34,6(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4,&
                                                 int5,int6
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
                  int5 = int5 + nfreeze
                  int6 = int6 + nfreeze
               endif
               if (nci(int4).le.nspace.and.nci(int5).le.nspace.and.&
                            nci(int6).le.nspace) then
                  if (int1.eq.int2) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 2
                  else if (int2.eq.int3) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 2
                  else
                     do l = 1, socc
                        if (int4.eq.isocc(l)) then
                              tmpocc(nci(int1)) = 1
                              tmpocc(nci(int2)) = 1
                              tmpocc(nci(int3)) = 1
                              tmpocc(nci(int4)) = 3
                              tmpocc(nci(int5)) = 2
                              tmpocc(nci(int6)) = 2
                        else if (int5.eq.isocc(l)) then
                              tmpocc(nci(int1)) = 1
                              tmpocc(nci(int2)) = 1
                              tmpocc(nci(int3)) = 1
                              tmpocc(nci(int4)) = 2
                              tmpocc(nci(int5)) = 3
                              tmpocc(nci(int6)) = 2
                        else if (int6.eq.isocc(l)) then
                              tmpocc(nci(int1)) = 1
                              tmpocc(nci(int2)) = 1
                              tmpocc(nci(int3)) = 1
                              tmpocc(nci(int4)) = 2
                              tmpocc(nci(int5)) = 2
                              tmpocc(nci(int6)) = 3
                        endif
                     enddo
                  endif
                  write(102,fmt1)coef,&
                          (tmpocc(k),k=1,nspace),'0:',0,&
                           '0:',0,'0:',0,&
                          '0:',0
               else if (nci(int4).le.nspace.and.nci(int5).le.&
                        nspace.and.nci(int6).gt.nspace) then
                  if (int1.eq.int2) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 2
                     tmpocc(nci(int5)) = 2
                     write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                  else if (int2.eq.int3) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 2
                     tmpocc(nci(int5)) = 2
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                          '0:',0,'0:',0,&
                          '0:',0
                  else
                     do l = 1, socc
                        if (int4.eq.isocc(l)) then
                              tmpocc(nci(int1)) = 1
                              tmpocc(nci(int2)) = 1
                              tmpocc(nci(int3)) = 1
                              tmpocc(nci(int4)) = 3
                              tmpocc(nci(int5)) = 2
                              write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                        else if (int5.eq.isocc(l)) then
                              tmpocc(nci(int1)) = 1
                              tmpocc(nci(int2)) = 1
                              tmpocc(nci(int3)) = 1
                              tmpocc(nci(int4)) = 2
                              tmpocc(nci(int5)) = 3
                              write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
                        else if (int6.eq.isocc(l)) then
                              tmpocc(nci(int1)) = 1
                              tmpocc(nci(int2)) = 1
                              tmpocc(nci(int3)) = 1
                              tmpocc(nci(int4)) = 2
                              tmpocc(nci(int5)) = 2
                              write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                           '0:',0,'0:',0,&
                          '0:',0
                        endif
                     enddo
                  endif
               else if (nci(int4).le.nspace.and.nci(int5).gt.&
                        nspace.and.nci(int6).gt.nspace) then
                  if (int1.eq.int2) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 1
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                  else if (int2.eq.int3) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 1
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                  else
                     do l = 1, socc
                        if (int4.eq.isocc(l)) then
                              tmpocc(nci(int1)) = 1
                              tmpocc(nci(int2)) = 1
                              tmpocc(nci(int3)) = 1
                              tmpocc(nci(int4)) = 3
                              write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                        else if (int5.eq.isocc(l)) then
                              tmpocc(nci(int1)) = 1
                              tmpocc(nci(int2)) = 1
                              tmpocc(nci(int3)) = 1
                              tmpocc(nci(int4)) = 2
                              write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                        else if (int6.eq.isocc(l)) then
                              tmpocc(nci(int1)) = 1
                              tmpocc(nci(int2)) = 1
                              tmpocc(nci(int3)) = 1
                              tmpocc(nci(int4)) = 2
                              write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                        endif
                     enddo
                  endif
               else if (nci(int4).gt.nspace.and.nci(int5).gt.&
                        nspace.and.nci(int6).gt.nspace) then
                  if (int1.eq.int2) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',1,&
                           nci(int5),':',2,nci(int6),':',2,&
                          '0:',0
                  else if (int2.eq.int3) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',1,&
                           nci(int5),':',2,nci(int6),':',2,&
                          '0:',0
                  else
                     do l = 1, socc
                        if (int4.eq.isocc(l)) then
                              tmpocc(nci(int1)) = 1
                              tmpocc(nci(int2)) = 1
                              tmpocc(nci(int3)) = 1
                              write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',3,&
                           nci(int5),':',2,nci(int6),':',2,&
                          '0:',0
                        else if (int5.eq.isocc(l)) then
                              tmpocc(nci(int1)) = 1
                              tmpocc(nci(int2)) = 1
                              tmpocc(nci(int3)) = 1
                              write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int5),':',3,nci(int6),':',2,&
                          '0:',0
                        else if (int6.eq.isocc(l)) then
                              tmpocc(nci(int1)) = 1
                              tmpocc(nci(int2)) = 1
                              tmpocc(nci(int3)) = 1
                              write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int5),':',2,nci(int6),':',3,&
                          '0:',0
                        endif
                     enddo
                  endif
               endif
            else if (label.eq.'T'.and.iorb.eq.7) then
               read(47,'(a5,f10.7,a34,6(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4,&
                                                 int5,int6
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
                  int5 = int5 + nfreeze
                  int6 = int6 + nfreeze
               endif
               if (nci(int4).le.nspace.and.nci(int5).le.nspace.and.&
                            nci(int6).le.nspace) then
                  tmpocc(nci(int1)) = 1
                  tmpocc(nci(int2)) = 1
                  tmpocc(nci(int3)) = 1
                  tmpocc(nci(int4)) = 2
                  tmpocc(nci(int5)) = 2
                  tmpocc(nci(int6)) = 2
                  write(102,fmt1)coef,&
                          (tmpocc(k),k=1,nspace),'0:',0,&
                           '0:',0,'0:',0,&
                          '0:',0
               else if (nci(int4).le.nspace.and.nci(int5).le.&
                        nspace.and.nci(int6).gt.nspace) then
                  tmpocc(nci(int1)) = 1
                  tmpocc(nci(int2)) = 1
                  tmpocc(nci(int3)) = 1
                  tmpocc(nci(int4)) = 2
                  tmpocc(nci(int5)) = 2
                  write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                           '0:',0,'0:',0,&
                          '0:',0
               else if (nci(int4).le.nspace.and.nci(int5).gt.&
                        nspace.and.nci(int6).gt.nspace) then
                  tmpocc(nci(int1)) = 1
                  tmpocc(nci(int2)) = 1
                  tmpocc(nci(int3)) = 1
                  tmpocc(nci(int4)) = 2
                  write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',2,'0:',3,&
                          '0:',0
               else if (nci(int4).gt.nspace.and.nci(int5).gt.&
                        nspace.and.nci(int6).gt.nspace) then
                  tmpocc(nci(int1)) = 1
                  tmpocc(nci(int2)) = 1
                  tmpocc(nci(int3)) = 1
                  write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int5),':',2,nci(int6),':',2,&
                          '0:',0
               endif
            else if (label.eq.'T'.and.iorb.ne.2) then
               read(47,'(a5,f10.7,a34,6(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4,&
                                                 int5,int6
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
                  int5 = int5 + nfreeze
                  int6 = int6 + nfreeze
               endif
               if (nci(int4).le.nspace.and.nci(int5).le.nspace.and.&
                            nci(int6).le.nspace) then
                  if (nci(int4).eq.nci(int5)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 2
                     tmpocc(nci(int5)) = 3
                     tmpocc(nci(int6)) = 2
                  else if (nci(int5).eq.nci(int6)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 2
                     tmpocc(nci(int4)) = 2
                     tmpocc(nci(int6)) = 3
                  else if (nci(int4).eq.nci(int6)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 2
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 3
                  else if (nci(int1).eq.nci(int2)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 2
                  else if (nci(int1).eq.nci(int3)) then
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 2
                  else if (nci(int2).eq.nci(int3)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 2
                  else
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 2
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 2
                  endif
                  write(102,fmt1)coef,&
                       (tmpocc(k),k=1,nspace),'0:',0,'0:',0,'0:',0,&
                       '0:',0
               else if (nci(int4).le.nspace.and.nci(int5).le.&
                    nspace.and.nci(int6).gt.nspace) then
                  if (nci(int4).eq.nci(int5)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 2
                     tmpocc(nci(int5)) = 3
                  else if (nci(int1).eq.nci(int2)) then
                     tmpocc(nci(int1)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                  else if (nci(int1).eq.nci(int3)) then
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                  else if (nci(int2).eq.nci(int3)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                  else
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 2
                     tmpocc(nci(int5)) = 2
                  endif
                  write(102,fmt2)coef,&
                       (tmpocc(k),k=1,nspace),nci(int6),':',2,'0:',0,&
                       '0:',0,'0:',0
               else if (nci(int4).le.nspace.and.nci(int5).gt.&
                         nspace.and.nci(int6).gt.nspace) then
                  if (nci(int1).eq.nci(int2)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 2
                  else if (nci(int1).eq.nci(int3)) then
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 2
                  else if (nci(int2).eq.nci(int3)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 2
                  else
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 2
                     tmpocc(nci(int4)) = 2
                  endif
                  if (nci(int5).eq.nci(int6)) then
                     write(102,fmt2)coef,&
                       (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                        '0:',0,'0:',0,&
                       '0:',0
                  else
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',1,'0:',0,&
                          '0:',0
                  endif
               else if (nci(int4).gt.nspace.and.nci(int5).gt.&
                        nspace.and.nci(int6).gt.nspace) then
                  if (nci(int1).eq.nci(int2)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                  else if (nci(int1).eq.nci(int3)) then
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 0
                  else if (nci(int2).eq.nci(int3)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                  else
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 1
                  endif
                  if (nci(int4).eq.nci(int5)) then
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',3,&
                           nci(int6),':',2,'0:',0,&
                          '0:',0
                  else if (nci(int5).eq.nci(int6)) then
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                  else if (nci(int4).eq.nci(int6)) then
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                           nci(int6),':',3,'0:',0,&
                          '0:',0
                  else
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int4),':',2,&
                           nci(int5),':',2,nci(int6),':',2,&
                          '0:',0
                  endif
               endif
            else if (label.eq.'Q'.and.iorb.eq.3) then
               read(47,'(a5,f10.7,a34,2(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4,&
                                                 int5,int6,int7,int8
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
                  int5 = int5 + nfreeze
                  int6 = int6 + nfreeze
                  int7 = int7 + nfreeze
                  int8 = int8 + nfreeze
               endif
               if (nci(int5).le.nspace.and.nci(int6).le.nspace.and.&
                       nci(int7).le.nspace.and.nci(int8).le.nspace) then
                  if (nci(int5).eq.nci(int6)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 2
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 2
                              endif
                        else if (int2.eq.isocc(l)) then
                              if (int1.eq.int3) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 3
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 3
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 2
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 2
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 2
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 2
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int6)) = 3
                           tmpocc(nci(int7)) = 2
                           tmpocc(nci(int8)) = 2
                        endif
                     enddo
                  else if (nci(int5).eq.nci(int7)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 3
                                 tmpocc(nci(int8)) = 2
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 3
                                 tmpocc(nci(int8)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 3
                                 tmpocc(nci(int8)) = 2
                              endif
                        else if (int2.eq.isocc(l)) then
                              if (int1.eq.int3) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 3
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int8)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 3
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int8)) = 2
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 3
                                 tmpocc(nci(int8)) = 2
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 3
                                 tmpocc(nci(int8)) = 2
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int6)) = 2
                           tmpocc(nci(int7)) = 3
                           tmpocc(nci(int8)) = 2
                        endif
                     enddo
                  else if (nci(int5).eq.nci(int8)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 3
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 3
                              endif
                        else if (int2.eq.isocc(l)) then
                              if (int1.eq.int3) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 3
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 3
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 2
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 3
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 3
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int6)) = 2
                           tmpocc(nci(int7)) = 2
                           tmpocc(nci(int8)) = 3
                        endif
                     enddo
                  else if (nci(int6).eq.nci(int7)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 3
                                 tmpocc(nci(int8)) = 2
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 3
                                 tmpocc(nci(int8)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 3
                                 tmpocc(nci(int8)) = 2
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int8)) = 2
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 3
                                 tmpocc(nci(int8)) = 2
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 3
                                 tmpocc(nci(int8)) = 2
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int5)) = 2
                           tmpocc(nci(int7)) = 3
                           tmpocc(nci(int8)) = 2
                        endif
                     enddo
                  else if (nci(int6).eq.nci(int8)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 3
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 3
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 3
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 2
                                 tmpocc(nci(int8)) = 3
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int5)) = 2
                           tmpocc(nci(int7)) = 2
                           tmpocc(nci(int8)) = 3
                        endif
                     enddo
                  else if (nci(int7).eq.nci(int8)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int8)) = 3
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int8)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int8)) = 3
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int8)) = 3
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int8)) = 3
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int5)) = 2
                           tmpocc(nci(int6)) = 2
                           tmpocc(nci(int8)) = 3
                        endif
                     enddo
                  endif
                  write(102,fmt1)coef,&
                       (tmpocc(k),k=1,nspace),'0:',0,'0:',0,'0:',0,&
                       '0:',0
               else if (nci(int5).le.nspace.and.nci(int6).le.&
                        nspace.and.&
                       nci(int7).le.nspace.and.nci(int8).gt.nspace) then
                  if (nci(int5).eq.nci(int6)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                              endif
                        else if (int2.eq.isocc(l)) then
                              if (int1.eq.int3) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 3
                                 tmpocc(nci(int7)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 3
                                 tmpocc(nci(int7)) = 2
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int6)) = 3
                                 tmpocc(nci(int7)) = 2
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int6)) = 3
                           tmpocc(nci(int7)) = 2
                        endif
                     enddo
                  else if (nci(int5).eq.nci(int7)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 3
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 3
                              endif
                        else if (int2.eq.isocc(l)) then
                              if (int1.eq.int3) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 3
                                 tmpocc(nci(int6)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 3
                                 tmpocc(nci(int6)) = 2
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 3
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int6)) = 2
                                 tmpocc(nci(int7)) = 3
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int6)) = 2
                           tmpocc(nci(int7)) = 3
                        endif
                     enddo
                  else if (nci(int6).eq.nci(int7)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 3
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 3
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 3
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 3
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int7)) = 3
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int5)) = 2
                           tmpocc(nci(int7)) = 3
                        endif
                     enddo
                  endif
                  write(102,fmt2)coef,&
                       (tmpocc(k),k=1,nspace),nci(int8),':',2,'0:',0,&
                        '0:',0,'0:',0
               else if (nci(int5).le.nspace.and.nci(int6).le.&
                        nspace.and.&
                       nci(int7).gt.nspace.and.nci(int8).gt.nspace) then
                  if (nci(int5).eq.nci(int6)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int6)) = 3
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int6)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int6)) = 3
                              endif
                        else if (int2.eq.isocc(l)) then
                              if (int1.eq.int3) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 3
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int6)) = 3
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int6)) = 3
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int6)) = 3
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int6)) = 3
                        endif
                     enddo
                       write(102,fmt3)coef,&
                       (tmpocc(k),k=1,nspace),nci(int7),':',2,&
                        nci(int8),':',2,&
                       '0:',0,'0:',0
                  else if (nci(int7).eq.nci(int8)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 2
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 2
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 2
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                                 tmpocc(nci(int6)) = 2
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int5)) = 2
                           tmpocc(nci(int6)) = 2
                        endif
                     enddo
                       write(102,fmt2)coef,&
                       (tmpocc(k),k=1,nspace),nci(int8),':',3,'0:',0,&
                       '0:',0,'0:',0
                  endif
               else if (nci(int5).le.nspace.and.nci(int6).gt.&
                        nspace.and.&
                       nci(int7).gt.nspace.and.nci(int8).gt.nspace) then
                  if (nci(int6).eq.nci(int7)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int5)) = 2
                        endif
                     enddo
                       write(102,fmt3)coef,&
                       (tmpocc(k),k=1,nspace),nci(int7),':',3,&
                         nci(int8),':',2,&
                       '0:',0,'0:',0
                  else if (nci(int6).eq.nci(int8)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                              endif
                        else
                          tmpocc(nci(int2)) = 0
                          tmpocc(nci(int4)) = 0
                          tmpocc(nci(int5)) = 2
                        endif
                     enddo
                       write(102,fmt3)coef,&
                       (tmpocc(k),k=1,nspace),nci(int7),':',2,&
                        nci(int8),':',3,&
                       '0:',0,'0:',0
                  else if (nci(int7).eq.nci(int8)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int5)) = 2
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                                 tmpocc(nci(int5)) = 2
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                                 tmpocc(nci(int8)) = 2
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                                 tmpocc(nci(int5)) = 2
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                           tmpocc(nci(int5)) = 2
                        endif
                     enddo
                       write(102,fmt3)coef,&
                       (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                        nci(int8),':',3,&
                       '0:',0,'0:',0
                  endif
               else if (nci(int5).gt.nspace.and.nci(int6).gt.&
                       nspace.and.&
                       nci(int7).gt.nspace.and.nci(int8).gt.nspace) then
                  if (nci(int5).eq.nci(int6)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                              endif
                        else if (int2.eq.isocc(l)) then
                              if (int1.eq.int3) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int4)) = 0
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                        endif
                     enddo
                       write(102,fmt4)coef,&
                       (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                        nci(int7),':',2,&
                       nci(int8),':',2,'0:',0
                  else if (nci(int5).eq.nci(int7)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                              endif
                        else if (int2.eq.isocc(l)) then
                              if (int1.eq.int3) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int4)) = 0
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                        endif
                     enddo
                       write(102,fmt4)coef,&
                       (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                        nci(int7),':',3,&
                       nci(int8),':',2,'0:',0
                  else if (nci(int5).eq.nci(int8)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                              endif
                        else if (int2.eq.isocc(l)) then
                              if (int1.eq.int3) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int4)) = 0
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                        endif
                     enddo
                       write(102,fmt4)coef,&
                       (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                        nci(int7),':',2,&
                       nci(int8),':',3,'0:',0
                  else if (nci(int6).eq.nci(int7)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                              endif
                        else if (int2.eq.isocc(l)) then
                              if (int1.eq.int3) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int4)) = 0
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                        endif
                     enddo
                       write(102,fmt4)coef,&
                       (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                        nci(int7),':',3,&
                       nci(int8),':',2,'0:',0
                  else if (nci(int6).eq.nci(int8)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                              endif
                        else if (int2.eq.isocc(l)) then
                              if (int1.eq.int3) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int4)) = 0
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                        endif
                     enddo
                       write(102,fmt3)coef,&
                       (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                        nci(int7),':',2,&
                       nci(int8),':',3,'0:',0
                  else if (nci(int7).eq.nci(int8)) then
                     do l = 1, socc
                        if (int1.eq.isocc(l)) then
                              if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              else if (int2.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 3
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 0
                                 tmpocc(nci(int2)) = 1
                                 tmpocc(nci(int4)) = 3
                              endif
                        else if (int2.eq.isocc(l)) then
                              if (int1.eq.int3) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              else if (int3.eq.int4) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int4)) = 0
                              endif
                        else if (int3.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 1
                              endif
                        else if (int4.eq.isocc(l)) then
                              if (int1.eq.int2) then
                                 tmpocc(nci(int2)) = 0
                                 tmpocc(nci(int3)) = 1
                                 tmpocc(nci(int4)) = 0
                              else if (int2.eq.int3) then
                                 tmpocc(nci(int1)) = 1
                                 tmpocc(nci(int3)) = 0
                                 tmpocc(nci(int4)) = 0
                              endif
                        else
                           tmpocc(nci(int2)) = 0
                           tmpocc(nci(int4)) = 0
                        endif
                     enddo
                       write(102,fmt4)coef,&
                       (tmpocc(k),k=1,nspace),nci(int5),':',2,&
                        nci(int6),':',2,&
                       nci(int8),':',3,'0:',0
                  endif
               endif
            else if (iorb.eq.8.or.(label.eq.'Q'.and.iorb.eq.1)) then
               read(47,'(a5,f10.7,a34,2(1x,i3))')char9,coef,char10,&
                                                 int1,int2,int3,int4,&
                                                 int5,int6,int7,int8
               if (lfreeze) then
                  int1 = int1 + nfreeze
                  int2 = int2 + nfreeze
                  int3 = int3 + nfreeze
                  int4 = int4 + nfreeze
                  int5 = int5 + nfreeze
                  int6 = int6 + nfreeze
                  int7 = int7 + nfreeze
                  int8 = int8 + nfreeze
               endif
               if (nci(int5).le.nspace.and.nci(int6).le.nspace.and.&
                       nci(int7).le.nspace.and.nci(int8).le.nspace) then
                  if (nci(int5).eq.nci(int6).and.nci(int7).eq.&
                         nci(int8)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int8)) = 3
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int8)) = 3
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int8)) = 3
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 2
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int8)) = 3
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int8)) = 3
                     endif
                  else if ((nci(int5).eq.nci(int7).and.nci(int6).eq.&
                         nci(int8)).or.(nci(int5).eq.nci(int8).and.&
                         nci(int7).eq.nci(int6))) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 3
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 3
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 3
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 2
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 3
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 3
                     endif
                  else if (nci(int5).eq.nci(int6)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 1
                       tmpocc(nci(int8)) = 2
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 2
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 2
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 2
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 2
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 2
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 2
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 2
                     endif
                  else if (nci(int5).eq.nci(int7)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 1
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     endif
                  else if (nci(int5).eq.nci(int8)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 1
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     endif
                  else if (nci(int6).eq.nci(int7)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 1
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                       tmpocc(nci(int8)) = 2
                     endif
                  else if (nci(int6).eq.nci(int8)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 1
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 2
                       tmpocc(nci(int8)) = 3
                     endif
                  else if (nci(int7).eq.nci(int8)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 1
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int8)) = 3
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int8)) = 3
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int8)) = 3
                     endif
                  else
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 2
                     tmpocc(nci(int7)) = 2
                     tmpocc(nci(int8)) = 2
                  endif
                  write(102,fmt1)coef,&
                       (tmpocc(k),k=1,nspace),'0:',0,'0:',0,'0:',0,&
                       '0:',0
               else if (nci(int5).le.nspace.and.nci(int6).le.&
                     nspace.and.nci(int7).le.nspace.and.nci(int8).gt.&
                     nspace) then
                  if (nci(int5).eq.nci(int6)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 1
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                       tmpocc(nci(int7)) = 2
                     endif
                  else if (nci(int5).eq.nci(int7)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 1
                       tmpocc(nci(int7)) = 3
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 2
                       tmpocc(nci(int7)) = 3
                     endif
                  else if (nci(int6).eq.nci(int7)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 1
                       tmpocc(nci(int7)) = 3
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int7)) = 3
                     endif
                  else
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 2
                     tmpocc(nci(int7)) = 2
                  endif
                  write(102,fmt2)coef,&
                       (tmpocc(k),k=1,nspace),nci(int8),':',2,'0:',0,&
                       '0:',0,'0:',0
               else if (nci(int5).le.nspace.and.nci(int6).le.&
                     nspace.and.nci(int7).gt.nspace.and.nci(int8).gt.&
                     nspace) then
                  if (nci(int5).eq.nci(int6).and.nci(int7).eq.&
                      nci(int8)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 2
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                     endif
                     write(102,fmt2)coef,&
                          (tmpocc(k),k=1,nspace),nci(int8),':',3,&
                          '0:',0,'0:',0,'0:',0
                  else if (nci(int5).eq.nci(int6)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 2
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int6)) = 3
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 2
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int6)) = 3
                     endif
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int7),':',1,&
                          nci(int8),':',2,'0:',0,'0:',0
                  else if (nci(int7).eq.nci(int8)) then
                     if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 1
                       tmpocc(nci(int6)) = 2
                     else if (nci(int1).eq.nci(int2)) then
                       tmpocc(nci(int2)) = 0
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                     else if (nci(int1).eq.nci(int3)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                     else if (nci(int1).eq.nci(int4)) then
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                     else if (nci(int2).eq.nci(int3)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 0
                       tmpocc(nci(int4)) = 1
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                     else if (nci(int2).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                     else if (nci(int3).eq.nci(int4)) then
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int4)) = 0
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                     else
                       tmpocc(nci(int1)) = 1
                       tmpocc(nci(int2)) = 1
                       tmpocc(nci(int3)) = 1
                       tmpocc(nci(int4)) = 2
                       tmpocc(nci(int5)) = 2
                       tmpocc(nci(int6)) = 2
                     endif
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int8),':',3,&
                          '0:',0,'0:',0,'0:',0
                  else
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                     tmpocc(nci(int6)) = 2
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int7),':',2,&
                          nci(int8),':',2,'0:',0,'0:',0
                  endif
               else if (nci(int5).le.nspace.and.nci(int6).gt.&
                    nspace.and.nci(int7).gt.nspace.and.nci(int8).gt.&
                    nspace) then
                  if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int4)) = 0
                     tmpocc(nci(int5)) = 1
                  else if (nci(int1).eq.nci(int2)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                  else if (nci(int1).eq.nci(int3)) then
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                  else if (nci(int1).eq.nci(int4)) then
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 0
                     tmpocc(nci(int5)) = 2
                  else if (nci(int2).eq.nci(int3)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                  else if (nci(int2).eq.nci(int4)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 0
                     tmpocc(nci(int5)) = 2
                  else if (nci(int3).eq.nci(int4)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int4)) = 0
                     tmpocc(nci(int5)) = 2
                  else
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 1
                     tmpocc(nci(int5)) = 2
                  endif
                  if (nci(int6).eq.nci(int7)) then
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int7),':',3,&
                          nci(int8),':',2,'0:',0,'0:',0
                  else if (nci(int6).eq.nci(int8)) then
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int7),':',2,&
                          nci(int8),':',3,'0:',0,'0:',0
                  else if (nci(int7).eq.nci(int8)) then
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                          nci(int8),':',3,'0:',0,'0:',0
                  else
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',2,&
                          nci(int7),':',2,nci(int8),':',2,'0:',0
                  endif
               else if (nci(int5).gt.nspace.and.nci(int6).gt.&
                    nspace.and.nci(int7).gt.nspace.and.nci(int8).gt.&
                    nspace) then
                  if (nci(int1).eq.nci(int2).and.nci(int3).eq.&
                        nci(int4)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int4)) = 0
                  else if (nci(int1).eq.nci(int2)) then
                     tmpocc(nci(int2)) = 0
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 2
                  else if (nci(int1).eq.nci(int3)) then
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 2
                  else if (nci(int1).eq.nci(int4)) then
                     tmpocc(nci(int2)) = 1
                     tmpocc(nci(int3)) = 2
                     tmpocc(nci(int4)) = 0
                  else if (nci(int2).eq.nci(int3)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 0
                     tmpocc(nci(int4)) = 2
                  else if (nci(int2).eq.nci(int4)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int3)) = 2
                     tmpocc(nci(int4)) = 0
                  else if (nci(int3).eq.nci(int4)) then
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 2
                     tmpocc(nci(int4)) = 0
                  else
                     tmpocc(nci(int1)) = 1
                     tmpocc(nci(int2)) = 2
                     tmpocc(nci(int3)) = 1
                     tmpocc(nci(int4)) = 2
                  endif
                  if (nci(int5).eq.nci(int6).and.nci(int7).eq.&
                      nci(int8)) then
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                          nci(int8),':',3,'0:',0,&
                          '0:',0
                  else if ((nci(int5).eq.nci(int7).and.nci(int6).eq.&
                      nci(int8)).or.(nci(int5).eq.nci(int8).and.&
                      nci(int6).eq.nci(int7))) then
                     write(102,fmt3)coef,&
                          (tmpocc(k),k=1,nspace),nci(int7),':',3,&
                          nci(int8),':',3,'0:',0,&
                          '0:',0
                  else if (nci(int5).eq.nci(int6)) then
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',3,&
                          nci(int7),':',1,nci(int8),':',2,&
                          '0:',0
                  else if (nci(int5).eq.nci(int7)) then
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',1,&
                          nci(int7),':',3,nci(int8),':',2,&
                          '0:',0
                  else if (nci(int5).eq.nci(int8)) then
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int6),':',1,&
                          nci(int7),':',2,nci(int8),':',3,&
                          '0:',0
                  else if (nci(int6).eq.nci(int7)) then
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',1,&
                          nci(int7),':',3,nci(int8),':',2,&
                          '0:',0
                  else if (nci(int6).eq.nci(int8)) then
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',1,&
                          nci(int7),':',2,nci(int8),':',3,&
                          '0:',0
                  else if (nci(int7).eq.nci(int8)) then
                     write(102,fmt4)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',1,&
                          nci(int6),':',2,nci(int8),':',3,&
                          '0:',0
                  else
                     write(102,fmt5)coef,&
                          (tmpocc(k),k=1,nspace),nci(int5),':',1,&
                          nci(int6),':',2,nci(int7),':',1,&
                          nci(int8),':',2
                  endif
               endif
            endif
            goto 710
         else
            goto 844
         endif
844      continue
         close(102)
      enddo


      close(47)

!c     Read in the number of aos used in the construction of
!c     of the mos
      open(57,file='mos',status='old')
      do i = 1, 3
         read(57,*)
      enddo
      read(57,'(a55,i5)')charr,nsaos
      rewind(57)

!c     initialize enval, mo2as
      allocate(enval(nsaos),mo2as(nsaos,nsaos),stat=ierr)
      if (ierr.ne.0) stop 'Memory allocation failed!'

      read(57,'(a46)')str1
      read(57,'(a46)')str2
      read(57,'(a1)')str3

!      write(6,*)nsaos
!      write(6,*)str1
!      write(6,*)str2
!      write(6,*)str3
!c     Read in mos and energies
      nchunk = nsaos / 4
      nleft = nsaos - (4 * nchunk)
      if(nleft.ne.0) nchunk = nchunk + 1
!      write(6,*)nchunk

      do i = 1, nsaos
         read(57,'(a26,d20.14)')char16,enval(i)
         do j = 1, nchunk
            if(nleft.ne.0) then
               if(j.le.(nchunk-1)) then
                  nchk = 4
               else
                  nchk = nleft
               endif
            else
               nchk = 4
            endif
   
            ichk = (j - 1) * 4
            read(57,'(4d20.14)')(mo2as(i,ichk+k),k=1,nchk)
         enddo
      enddo

      close(57)
!      write(6,*)nchunk,nleft
!      write(6,*)enval(nsaos)
!      write(6,*)mo2as(nsaos,:)

!c     Write ordered mos as mos_new
      open(59,file='mos_new',form='formatted',status='unknown')
      write(59,'(a46)')str1
      write(59,'(a46)')str2
      write(59,'(a1)')str3
      do i = 1, nsaos
         if (i.le.ncount) then
            write(59,'(3x,i3,2x,a1,6x,a11,d20.14,3x,a6,i3,a22)')i,'a',&
               'eigenvalue=',enval(ntm(i)),'nsaos=',nsaos,''
            write(59,'(4d20.14)') (mo2as(ntm(i),j),j=1,nsaos)
         else
            write(59,'(3x,i3,2x,a1,6x,a11,d20.14,3x,a6,i3,a22)')i,'a',&
                'eigenvalue=',enval(i),'nsaos=',nsaos,''
            write(59,'(4d20.14)') (mo2as(i,j),j=1,nsaos)
         endif
      enddo
      write(59,'(a)')'$end'

      close(59)

      inquire(file='mos_cart',exist=lexist)
      if (lexist) then
         open(61,file='mos_cart',status='old')
         read(61,'(i4)') msaos

         allocate(cartas(nsaos,msaos),stat=ierr)
         if (ierr.ne.0) stop 'Memory allocation failed!'

         nchunk = msaos / 4
         nleft = msaos - (4 * nchunk)
         if(nleft.ne.0) nchunk = nchunk + 1
!         write(6,*)nchunk

         do i = 1, nsaos
            read(61,'(i4)') itmp
            do j = 1, nchunk
               if(nleft.ne.0) then
                  if(j.le.(nchunk-1)) then
                     nchk = 4
                  else
                     nchk = nleft
                  endif
               else
                  nchk = 4
               endif

               ichk = (j - 1) * 4
               read(61,'(4e20.14)')(cartas(i,ichk+k),k=1,nchk)
            enddo
         enddo

         close(61)
!         write(6,*)nchunk,nleft
!         write(6,*)enval(nsaos)
!         write(6,*)mo2as(nsaos,:)

!c        Write ordered mos as mos_new
         open(62,file='mos_new.cart',form='formatted',status='unknown')
         write(62,'(i4)')msaos
         do i = 1, nsaos
            if (i.le.ncount) then
               write(62,'(i4)')i
               write(62,'(4e20.14)') (cartas(ntm(i),j),j=1,msaos)
            else
               write(62,'(i4)')i
               write(62,'(4e20.14)') (cartas(i,j),j=1,msaos)
            endif
         enddo

         close(62)
         deallocate(cartas)
      endif

      deallocate(mo2as,ntm)

    end program csfgenerator
