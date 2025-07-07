!################## DMBFEPM.f90#############################
! FEPMに転位密度ベースモデル(Dislocation Mechanism-based model)を導入
! (笠井のプログラムを改変)
! 2021年8月開発開始 by T.Oya
! 2022年7月のバージョン
! 2024年8月小松改修(Taylor-DMB切り替え)
! 2025年4月小松改修(勉強用-複雑化の解消)
!----------------------------------------------------------------------
module prescribed
  implicit none
!----------------------------------------------------------------------
!要素数を変化させたい時は memx,memy,memz=()の数字を変化させる(1辺あたりの要素数)
integer,parameter::memx=10, memy=10, memz=10
double precision,parameter::elesize=1. !1辺の寸法
!mem1は１辺が含む節点の数、nodは多結晶モデル中の全ての節点数、melは多結晶モデル中の全ての要素数
integer,parameter::memx1=memx+1, memy1=memy+1, memz1=memz+1
integer,parameter::mel=memx*memy*memz, nod=memx1*memy1*memz1
integer,parameter::nod3=nod*3, nm3=3*memx1*memy1*memz1
integer,parameter::nband=3*(memx1*memy1+memx1+2), mm3=3*(memx1*memy1+memx1+2)
!----------------------------------------------------------------------
! 初期結晶方位aaa
! iso=1; isotropic 等方性(ランダムに生成)
! iso=0; anisotropic　異方性(ファイルから供給)
integer,parameter::iso=1

double precision,parameter::pi=3.141592653589793238462643 ! 円周率
!----------------------------------------------------------------------
! 変形条件
double precision::deex=0.0005!0.0005 !x方向単軸引張り時のひずみ増分
double precision::deey=0
double precision::deez=0
double precision::eexmax=0.1 ! 変形終了時のひずみ
double precision::eeymax=0
double precision::eezmax=0

!----------------------------------------------------------------------
! 材料パラメータ
double precision,parameter::yg=70000., pr=0.33 !ヤング率(MPa),ポアソン比
double precision,parameter::lambda=pr*yg/((1.+pr)*(1.-2.*pr)), rigid=yg/(2.*(1.+pr))
double precision,parameter::dt=0.5/(2*rigid)
double precision,parameter::taylor=3.07, hard=0., s0=160., e0=0., ykmin=5.
double precision,dimension(mel,12)::yield !分解せん断応力
double precision,parameter:: ss0=200, voS=400, voH=10
double precision,parameter:: sw0=200, sw1=0.0001, sw2=0.2


!----------------------------------------------------------------------
!計算過程で出てくる変数等の宣言
integer::istep, iter
double precision,dimension(6)::ss, eep
double precision::error
double precision,dimension(nod)::x, y, z
double precision,dimension(mel,6)::s, s1
double precision,dimension(mel)::phi, theta, psi
double precision,dimension(mel,6)::dep
double precision,dimension(6)::deep, deep0
double precision::yk
double precision::ex
double precision,dimension(mel,12)::ksi, act
integer::total_act
double precision::total_act_ratio
integer,dimension(mel,8)::nel
double precision,dimension(mel)::xc, yc, zc
double precision,dimension(mel,8)::xe, ye, ze, xr, yr, zr
double precision,dimension(3,3)::r
double precision,dimension(mel,12)::xv, yv, zv
double precision,dimension(6,6)::d
double precision::area
double precision,dimension(mel,6,12)::eschmid, wschmid
double precision,dimension(nod3,nband)::tk
double precision,dimension(24,24)::ak
double precision,dimension(4)::lsub
double precision,dimension(6,12)::bsub
double precision,dimension(3,12)::csub
double precision::volum
double precision,dimension(mel,6,24)::b
double precision,dimension(mel,24,6)::bd
double precision,dimension(mel,3,24)::c
double precision,dimension(mel)::tvol
double precision,dimension(mel,3)::dw, dwp
double precision,dimension(nod3)::du, f

double precision::gamaa
double precision,dimension(mel)::gamma_crystal, agamma
double precision,dimension(mel,12)::g, dg, dgg, tau0, tau
double precision,dimension(mel,12)::gammaa, gamb,gamc,gama,dggg

end module prescribed
!----------------------------------------------------------------------

!********************************************************************
!************************** start of main ***************************
!********************************************************************
program main
use prescribed
  implicit none
!———————————————————————————————————————————————————————————————————————————
double precision::rand
double precision::effstrain !,deep1, eep11, ss1
!double precision::dd, ee, ff, ddd, eee, fff
integer::i, l, m, n
double precision::t1,t2
double precision::eqs
double precision ,dimension(4*mel)::list
double precision ,dimension(12*mel)::slist
integer,dimension(8)::date_time
character*10::date_b(3)
character*20::numb
character*12::nameoffile
integer:: komai,komaj,komatt,komaty
double precision:: xmin,xmax,xminav,xmaxav,dxkom
character(len=6) :: path

!———————————————————————————————————————————————————————————————————————————
call cpu_time(t1) !計算時間計測開始
!———————————————————————————————————————————————————————————————————————————
!------Initial Euler angles----------------------------
if(iso.eq.1) then ! 等方性の場合はランダム生成
do m=1,mel
  call random_number(rand)
  phi(m)=acos(1.-2.*rand) !0 ≦ phi ≦ π
  theta(m)=2.*pi*rand !0 ≦ theta ≦ 2π
  psi(m)=2.*pi*rand   !0 ≦ psi ≦ 2π
end do
else
  open(1,file='anisotropy.csv') !異方性の場合の外部ファイル入力
do m=1,mel
  read(1,*) phi(m),theta(m),psi(m)
end do
  close(1)
endif


call initial
  effstrain=0
!---------start of deformation step(変形のステップのはじめ)--------------------
  istep=0
!---------header for csv file----------
  call date_and_time(date_b(1), date_b(2), date_b(3), date_time)
  write(*,*) '********************************************'
  write(*,*) '********Welcome to DMBFEPM simulation !********'
  write(*,*) '********************************************'
  write(*,*) '--フォルダ名を入力してください(５文字):--'
  read *, nameoffile
  call create_folder(nameoffile)
path =trim(nameoffile) // "/"
write(*,*) '********************************************'
open(3,file=path//'1sscurve.csv',status='replace')
  write(3,*) 'Date of Simulation:',',',date_b(1)
  write(3,*) 'Material parameters:'
  write(3,*) 'Young modulus (MPa)',',',yg
  write(3,*) 'Poisson ratio',',',pr
  write(3,*) 'Rigidity',',',rigid
  write(3,*) 'Latent hardening components:'
  write(3,*) 'Tensile strain increments:'
  write(3,*) 'deex',',',deex
  write(3,*) '-----',',','-----',',','-----',',','-----'
  write(3,*) 'true strain x',',','true strain y',',','true strain z',&
     &',','true stress x',',','true stress y',&
     &',','true stress z',',','ss(6)',',','eep(6)'
!-----end of header
!--------- start of outer iteration(反復の始まり)------------------------------
10000 continue
  istep=istep+1
!---transformation of coordinates with Euler angles（オイラー角を用いた座標変換）
 call trans
!----construction of rigidity matrix tk(剛性マトリックスの構築)-------------
 call tkmatrix
!---decomposition of TK matrix with Cholesky's method(コレスキー法による剛性行列の分解)---
 call matrix
!------make increments data zero(増加量の初期化)--------------------------
 call clear
!--------- start of inner iteration(逐次累積ループの始まり)----------------------
  iter= 0
20000 continue
  iter=iter+1
!-----determination of plastic strain(逐次累積法による塑性ひずみの予測)------
 call find
!--------- termination check of inner iteration(逐次累積ループの終了判定)--------
if(iter > 5000) then
  write(*,*) 'stop because of divergence'
  write(3,*) 'stop because of divergence'
  close(3)
 stop
endif
if (error > 1.e-5) goto 20000
!-----renewal of coordinates,stress and euler angles(座標、応力、オイラー角の更新)
 call renew
!--------- Total number of active slip systems ----------------------------
do m=1,mel
  do l=1,12
    act(m,l)=0
  enddo
enddo

do m=1,mel
  do l=1,12
    if (dg(m,l) /=0) then
      act(m,l)=1
    endif
  enddo
enddo

total_act=0 !活動すべり系の数(初期値0)
total_act_ratio=0. !活動すべり系率(初期値0.)
do m=1,mel
  do l=1,12
    total_act=total_act+act(m,l)
  enddo
enddo
total_act_ratio=total_act/(12.*memx*memy*memz)


!初期立ち上がりの検討のための、巨視的ひずみ計算
xmin=0
xmax=0
do komai=1,memy1
  do komaj=1,memz1
    komatt=1+(komai-1)*memx1+(komaj-1)*memx1*memy1
    komaty=memx1+(komai-1)*memx1+(komaj-1)*memx1*memy1
    xmin=xmin+x(komatt)
    xmax=xmax+x(komaty)
  end do    
end do
xminav=xmin/memy1/memz1
xmaxav=xmax/memy1/memz1
dxkom=xmaxav-xminav-1
!------------------------------------------------------------------------------
  write(3,*) eep(1),',',dxkom,',',ss(1)!,',',total_act,',',total_act_ratio!,'
  write(numb,*) istep
!------------------------------------------------------------------------------
! Effective plastic strain
effstrain=sqrt(2*(eep(1)**2+eep(2)**2+eep(3)**2)/3+(eep(4)**2+eep(5)**2+eep(6)**2)/3)




if(effstrain < eexmax ) goto 10000		!相当塑性ひずみがひずみ限界eexmaxに達していない場合再度計算を行う
!----------------------------------------------------------------------
!open(9,file='dg1.csv',status='replace') ! すべり系ごとのせん断ひずみ増分値の出力(活動/非活動の判定に使用)
!do m=1,mel
!    write(9,*) dg(m,1),dg(m,2),dg(m,3),dg(m,4),dg(m,5),dg(m,6),dg(m,7),dg(m,8),dg(m,9),dg(m,10),dg(m,11),dg(m,12)
!enddo
!do m=1,mel!交代率を求めるときの１次負荷の時の配列
!   do l=0,3
!       ddd=dg(m,l*3+1)
!       eee=dg(m,l*3+2)
!       fff=dg(m,l*3+3)
!       if(ddd /= 0 .or. eee/=0 .or. fff /= 0)then
!         list(4*(m-1)+(l+1))=1
!       else
!         list(4*(m-1)+(l+1))=0
!       endif
!   enddo
!   do l=1,12
!         slist(12*(m-1)+l)=act(m,l)
!   enddo
!enddo
!close(9)

call cpu_time(t2)
    print *,"cpu time:",t2-t1,"seconds"
stop
end program
!********************************************************************
!************************** end of main *****************************
!********************************************************************



!*************************************************************************
!************************* Subroutines ***********************************
!*************************************************************************

!********************************************************************
!****subroutine initial**********************************************
!********************************************************************
subroutine initial
use prescribed
 implicit none
!———————————————————————————————————————————————————————————————————————————
integer::i, j, k, l, m, n
double precision::dx, dy, dz
!———————————————————————————————————————————————————————————————————————————
!----- coordinates of nodal points (節点座標)----------------------------
    dx=elesize/memx         !x軸の要素間距離、1を要素数で割っている
    dy=elesize/memy         !y軸
    dz=elesize/memz         !z軸
    area=dx*memx*dy*memy	!xy断面積
!変形開始前の節点座標の出力
!open(21,file='beforenode',status='replace')	!
do k=1,memz1
 do j=1,memy1
  do i=1,memx1
    n=i+(j-1)*memx1+(k-1)*memx1*memy1   !結晶モデル中のn(=mem+1の三乗)個の節点について番号と座標の定義
    x(n)=dx*(i-1)
    y(n)=dy*(j-1)
    z(n)=dz*(k-1)
!write(21,'(I8,3f16.7)') n,x(n),y(n),z(n)
  end do
 end do
end do
!close(21)
!-----nodal points constructing an element　（要素を構成する節点）-----------
!open(22,file='elementsolid',status='replace')	!
do k=1,memz
 do j=1,memy
  do i=1,memx                        !一つの要素を構成する節点8個の番号を求めている
    m=i+(j-1)*memx+(k-1)*memx*memy
    nel(m,1)=i+(j-1)*memx1+(k-1)*memx1*memy1
    nel(m,2)=nel(m,1)+1
    nel(m,4)=i+j*memx1+(k-1)*memx1*memy1
    nel(m,3)=nel(m,4)+1
    nel(m,5)=nel(m,1)+memx1*memy1
    nel(m,6)=nel(m,5)+1
    nel(m,8)=nel(m,4)+memx1*memy1
    nel(m,7)=nel(m,8)+1
    xc(m)=(x(nel(m,1))+x(nel(m,2))+x(nel(m,3))+x(nel(m,4))+x(nel(m,5))+x(nel(m,6))+x(nel(m,7))+x(nel(m,8)))/8.
    yc(m)=(y(nel(m,1))+y(nel(m,2))+y(nel(m,3))+y(nel(m,4))+y(nel(m,5))+y(nel(m,6))+y(nel(m,7))+y(nel(m,8)))/8.
    zc(m)=(z(nel(m,1))+z(nel(m,2))+z(nel(m,3))+z(nel(m,4))+z(nel(m,5))+z(nel(m,6))+z(nel(m,7))+z(nel(m,8)))/8.
    do n=1,8
     xr(m,n)=x(nel(m,n))-xc(m)
     yr(m,n)=y(nel(m,n))-yc(m)
     zr(m,n)=z(nel(m,n))-zc(m)
    enddo

!write(22,'(10I8)') m,1,nel(m,1),nel(m,2),nel(m,3),nel(m,4),nel(m,5),nel(m,6),nel(m,7),nel(m,8)

  end do
 end do
end do
!close(22)
!---------------------------------------------------------------------
!--------D-matrix------（p34のD行列を定義）------------------------------
do i=1,6
 do j=1,6
      d(i,j)=0.!全ての要素を0とする
 end do
end do

do i=1,3
 do j=1,3
    d(i,j)=lambda!1-3にλを代入
 end do
    d(i,i)=lambda+2.*rigid!(1,1)(2,2)(3,3)にλ+2Gを代入
    d(i+3,i+3)=rigid!(4,4)(5,5)(6,6)にλを代入
end do
!---------------------------------------------------------------------
!-----make a fresh start----------------------------------------------
!(全ての要素(1<m<mel)に対して、応力s1(m,1)からs1(m,6)、gamma_crystal(m)、ひずみeep(1)からeep(6)の全てが0としている)
do m=1,mel
 do i=1,6
      s1(m,i)=0.
 end do
end do

do i=1,6
    eep(i)=0.
end do

return
end
!*****************************************************************************
!******* subroutine trans ****************************************************
!*****************************************************************************
subroutine trans
use prescribed
 implicit none
!———————————————————————————————————————————————————————————————————————————
integer::i, j, l, m, n, ii, jj
double precision,dimension(3,3,12)::ce, cw, cee, cww
double precision,dimension(3,12)::mat_a, mat_b
double precision,dimension(3)::as, bs
!----- fcc slip system (滑り系 表1.1のFCC金属のすべり系を再現)-------------------
 data mat_a/1, 1, 1, 1, 1, 1, 1, 1, 1,&
         & -1,-1,1, -1,-1,1, -1,-1,1,&
          & -1,1, 1, -1,1, 1, -1,1, 1,&
          & 1, -1, 1, 1, -1, 1, 1, -1, 1/
 data mat_b/ 0, 1, -1, -1, 0, 1, 1,-1, 0,&
          & 0, -1,-1, 1, 0, 1, -1, 1, 0,&
          & 0, 1, -1, 1, 0,1, -1, -1, 0,&
          & 0, -1,-1, -1, 0,1, 1,1, 0/
!———————————————————————————————————————————————————————————————————————————
!--- microscopic components of Schmid tensor(ミクロSchmidテンソル)------------
do l=1,12
 do i=1,3
      as(i)=mat_a(i,l)/sqrt(3.)!asはすべり面の単位法線ベクトル
      bs(i)=mat_b(i,l)/sqrt(2.)!bsはすべり方向の単位ベクトル
 end do
  do i=1,3
   do j=1,3
      ce(i,j,l)=(as(i)*bs(j)+as(j)*bs(i))/2.!ひずみについてのSchmidテンソル
      cw(i,j,l)=(as(i)*bs(j)-as(j)*bs(i))/2.!スピンについてのSchmidテンソル
   end do
  end do
end do
!---------------------------------------------
!----- transformation matrix------------------
do m=1,mel
  r(1,1)=cos(phi(m))*cos(theta(m))*cos(psi(m))-sin(theta(m))*sin(psi(m))
  r(1,2)=cos(phi(m))*sin(theta(m))*cos(psi(m))+cos(theta(m))*sin(psi(m))
  r(1,3)=-sin(phi(m))*cos(psi(m))
  r(2,1)=-cos(phi(m))*cos(theta(m))*sin(psi(m))-sin(theta(m))*cos(psi(m))
  r(2,2)=-cos(phi(m))*sin(theta(m))*sin(psi(m))+cos(theta(m))*cos(psi(m))
  r(2,3)=sin(phi(m))*sin(psi(m))
  r(3,1)=sin(phi(m))*cos(theta(m))
  r(3,2)=sin(phi(m))*sin(theta(m))
  r(3,3)=cos(phi(m))

  if(istep==1)then
    do n=1,8
      xe(m,n)=xc(m)+r(1,1)*xr(m,n)+r(1,2)*yr(m,n)+r(1,3)*zr(m,n)
      ye(m,n)=yc(m)+r(2,1)*xr(m,n)+r(2,2)*yr(m,n)+r(2,3)*zr(m,n)
      ze(m,n)=zc(m)+r(3,1)*xr(m,n)+r(3,2)*yr(m,n)+r(3,3)*zr(m,n)
    enddo
      xv(m,1)=xe(m,4)-xe(m,5)
      xv(m,2)=xe(m,5)-xe(m,2)
      xv(m,3)=xe(m,2)-xe(m,4)
      xv(m,4)=xe(m,2)-xe(m,7)
      xv(m,5)=xe(m,7)-xe(m,4)
      xv(m,6)=xe(m,4)-xe(m,2)
      xv(m,7)=xe(m,3)-xe(m,6)
      xv(m,8)=xe(m,6)-xe(m,1)
      xv(m,9)=xe(m,1)-xe(m,3)
      xv(m,10)=xe(m,1)-xe(m,8)
      xv(m,11)=xe(m,8)-xe(m,3)
      xv(m,12)=xe(m,3)-xe(m,1)
      yv(m,1)=ye(m,4)-ye(m,5)
      yv(m,2)=ye(m,5)-ye(m,2)
      yv(m,3)=ye(m,2)-ye(m,4)
      yv(m,4)=ye(m,2)-ye(m,7)
      yv(m,5)=ye(m,7)-ye(m,4)
      yv(m,6)=ye(m,4)-ye(m,2)
      yv(m,7)=ye(m,3)-ye(m,6)
      yv(m,8)=ye(m,6)-ye(m,1)
      yv(m,9)=ye(m,1)-ye(m,3)
      yv(m,10)=ye(m,1)-ye(m,8)
      yv(m,11)=ye(m,8)-ye(m,3)
      yv(m,12)=ye(m,3)-ye(m,1)
      zv(m,1)=ze(m,4)-ze(m,5)
      zv(m,2)=ze(m,5)-ze(m,2)
      zv(m,3)=ze(m,2)-ze(m,4)
      zv(m,4)=ze(m,2)-ze(m,7)
      zv(m,5)=ze(m,7)-ze(m,4)
      zv(m,6)=ze(m,4)-ze(m,2)
      zv(m,7)=ze(m,3)-ze(m,6)
      zv(m,8)=ze(m,6)-ze(m,1)
      zv(m,9)=ze(m,1)-ze(m,3)
      zv(m,10)=ze(m,1)-ze(m,8)
      zv(m,11)=ze(m,8)-ze(m,3)
      zv(m,12)=ze(m,3)-ze(m,1)

  endif
!-------------------------------------------------------
!---- macroscopic components of Schmid tensor-----------
 do l=1,12
  do i=1,3
   do j=1,3
      cee(i,j,l)=0
      cww(i,j,l)=0
    do ii=1,3
     do jj=1,3
      cee(i,j,l)=cee(i,j,l)+ce(ii,jj,l)*r(ii,i)*r(jj,j)
      cww(i,j,l)=cww(i,j,l)+cw(ii,jj,l)*r(ii,i)*r(jj,j)
     end do
    end do
   end do
  end do
    eschmid(m,1,l)=cee(1,1,l)
    eschmid(m,2,l)=cee(2,2,l)
    eschmid(m,3,l)=cee(3,3,l)
    eschmid(m,4,l)=cee(2,3,l)*2
    eschmid(m,5,l)=cee(3,1,l)*2
    eschmid(m,6,l)=cee(1,2,l)*2
    wschmid(m,1,l)=cww(2,3,l)
    wschmid(m,2,l)=cww(3,1,l)
    wschmid(m,3,l)=cww(1,2,l)
 end do
end do

return
end
!******************************************************************************
!***** subroutine tkmatrix*****************************************************
!******************************************************************************
subroutine tkmatrix
use prescribed
 implicit none
!———————————————————————————————————————————————————————————————————————————
integer::i, j, k, m, n, nb, np, jj, lli, llj
integer,dimension(24)::ll
!———————————————————————————————————————————————————————————————————————————
do i=1,nod3
 do j=1,nband
      tk(i,j)=0.
 end do
end do

do m=1,mel
!-----rigidity matrix of an element
  call akmatrix(m)
!---------------------------------------------
 do i=1,8
  do j=1,3
      ll(i*3+j-3)=3*nel(m,i)+j-3
  end do
 end do
 do i=1,24
  do j=1,24
      if (ll(i).le.ll(j)) then
       lli=ll(i)
       llj=ll(j)-ll(i)+1
       tk(lli,llj)=tk(lli,llj)+ak(i,j)
      endif
  end do
 end do
end do
!————————————————
!——boundary condition of TK-matrix——
do k=1,memz1
 do j=1,memy1
  do i=1,memx1
      nb=(i-1)*(j-1)*(k-1)*(i-memx1)*(j-memy1)*(k-memz1)
    if (nb.eq.0) then
      n=i+(j-1)*memx1+(k-1)*memx1*memy1
     do jj=1,3
      np=3*(n-1)+jj
      tk(np,1)=1.e20
     end do
    endif
  end do
 end do
end do

return
end
!****************************************************************************
!***** subroutine akmatrix**************************************************
!***************************************************************************
subroutine akmatrix(m)
use prescribed
 implicit none
!———————————————————————————————————————————————————————————————————————————
integer::i, j, ij, ii, jj, k, ik, jk, ik1, jk1, m
integer,dimension(12,12)::aksub
integer,dimension(12,6)::bdsub
integer,dimension(10,4)::lfour
!———————————————————————————————————————————————————————————————————————————
!—— tetrahedral elements in a cubic element——
data lfour/ 1,3,8,6,1,4,2,5,7,5,2,4,7,5,6,1,3,8,6,7,3,1,6,8,3,2,4,7,5,2,6,8,3,1,8,5,7,4,2,4 /
!————————————————
   tvol(m)=0.
do i=1,6
 do j=1,24
      b(m,i,j)=0.
 end do
end do

do k=1,3
 do j=1,24
      c(m,k,j)=0.
 end do
end do

do i=1,24
 do j=1,24
      ak(i,j)=0.
 end do
end do
!-----———————————————————————————————————————————————————————————————————————
do ij=1,10
 do j=1,4
      lsub(j)=nel(m,lfour(ij,j))
 end do
!————————————————
!—— B-matrix of a tetrahedral element——
 call bmatrix
!————————————————
do i=1,12
 do j=1,6
      bdsub(i,j)=0.
  do k=1,6
      bdsub(i,j)=bdsub(i,j)+bsub(k,i)*d(k,j)
  end do
 end do
end do

do i=1,12
 do j=1,12
      aksub(i,j)=0.
  do k=1,6
      aksub(i,j)=aksub(i,j)+bdsub(i,k)*bsub(k,j)*volum
  end do
 end do
end do
!————AKmatrix————————————
do i=1,4
 do j=1,4
      ik=(lfour(ij,i)-1)*3
      jk=(lfour(ij,j)-1)*3
      ik1=(i-1)*3
      jk1=(j-1)*3
  do ii=1,3
   do jj=1,3
      ak(ik+ii,jk+jj)=ak(ik+ii,jk+jj)+aksub(ik1+ii,jk1+jj)
   end do
  end do
 end do
end do
!———————Bmatrix of a cubic element------------------
do i=1,6
 do j=1,4
      jk=(lfour(ij,j)-1)*3
      jk1=(j-1)*3
  do jj=1,3
      b(m,i,jk+jj)=b(m,i,jk+jj)+bsub(i,jk1+jj)*volum
  end do
 end do
end do
!----- Cmatrix------
do k=1,3
 do j=1,4
      jk=(lfour(ij,j)-1)*3
      jk1=(j-1)*3
  do jj=1,3
      c(m,k,jk+jj)=c(m,k,jk+jj)+csub(k,jk1+jj)*volum
  end do
 end do
end do
    tvol(m)=tvol(m)+volum
end do
!-----———————————————————————————————————————————————————————————————————————
do i=1,6
 do j=1,24
      b(m,i,j)=b(m,i,j)/tvol(m)
 end do
end do

do k=1,3
 do j=1,24
      c(m,k,j)=c(m,k,j)/tvol(m)
 end do
end do

  tvol(m)=tvol(m)/2.

do i=1,24
 do j=1,24
      ak(i,j)=ak(i,j)/2.
 end do
end do
!————————————————----
!----- BD matrix-----
do i=1,24
 do j=1,6
      bd(m,i,j)=0.0
  do k=1,6
      bd(m,i,j) =bd(m,i,j)+b(m,k,i)*d(k,j)
  end do
 end do
end do

return
end
!***************************************************************************
!***** subroutine bmatrix***************************************************
!***************************************************************************
subroutine bmatrix
use prescribed
 implicit none
!———————————————————————————————————————————————————————————————————————————
integer::i, j, k, ik, i3
double precision::val
integer::ls
double precision,dimension(3,3)::ai, bi, ci, di
volum=0.
!———————————————————————————————————————————————————————————————————————————
do i=1,6
 do j=1,12
    bsub(i,j)=0.
 end do
end do
do k=1,3
 do j=1,12
    csub(k,j)=0.
 end do
end do
!-----———————————————————————————————————————————————————————————————————————
do i=1,4
 do k=1,3
    ik=i+k
  if(ik>4) ik=ik-4
    ls=lsub(ik)
    ai(k,1)=x(ls)
    ai(k,2)=y(ls)
    ai(k,3)=z(ls)
    bi(k,1)=1.
    bi(k,2)=y(ls)
    bi(k,3)=z(ls)
    ci(k,1)=x(ls)
    ci(k,2)=1.
    ci(k,3)=z(ls)
    di(k,1)=x(ls)
    di(k,2)=y(ls)
    di(k,3)=1.
 end do
!---------determinant----------
 call det(ai,val)
!————————————————
    volum=volum+val*(-1.)**(i-1)
    i3=(i-1)*3
!————————————————
 call det(bi,val)
!————————————————
    bsub(1,i3+1)=val*(-1.)**i
    bsub(6,i3+2)=val*(-1.)**i
    bsub(5,i3+3)=val*(-1.)**i
    csub(3,i3+2)=val*(-1.)**i
    csub(2,i3+3)=-val*(-1.)**i
!————————————————
 call det(ci,val)
!————————————————
    bsub(2,i3+2)=val*(-1.)**i
    bsub(6,i3+1)=val*(-1.)**i
    bsub(4,i3+3)=val*(-1.)**i
    csub(3,i3+1)=-val*(-1.)**i
    csub(1,i3+3)=val*(-1.)**i
!————————————————
 call det(di,val)
!————————————————
    bsub(3,i3+3)=val*(-1.)**i
    bsub(4,i3+2)=val*(-1.)**i
    bsub(5,i3+1)=val*(-1.)**i
    csub(1,i3+2)=-val*(-1.)**i
    csub(2,i3+1)=val*(-1.)**i
end do
!-----———————————————————————————————————————————————————————————————————————
do i=1,6
  do j=1,12
    bsub(i,j)=bsub(i,j)/volum
 end do
end do
do k=1,3
 do j=1,12
    csub(k,j)=csub(k,j)/(2.*volum)
 end do
end do
    volum=volum/6.
!-----———————————————————————————————————————————————————————————————————————
return
end
!****************************************************************************
!***** subroutine det（行列のdeterminant計算）*********************************
!****************************************************************************
subroutine det(a,val)
  implicit none
!———————————————————————————————————————————————————————————————————————————
double precision,dimension(3,3)::a
double precision::val
!———————————————————————————————————————————————————————————————————————————
  val=a(1,1)*a(2,2)*a(3,3)+a(2,1)*a(3,2)*a(1,3)+a(1,2)*a(2,3)*a(3,1)	&
  -a(1,3)*a(2,2)*a(3,1)-a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
return
end
!***************************************************************************
!***** subroutine clear (初期化)*********************************************
!***************************************************************************
subroutine clear
use prescribed
  implicit none
!———————————————————————————————————————————————————————————————————————————
integer::i,j,l,m
!———————————————————————————————————————————————————————————————————————————
do m=1,mel
 do i=1,6
      dep(m,i)=0.
 end do
 do j=1,3
      dw(m,j)=0.
 end do
 do l=1,12
      dg(m,l)=0.
 end do
end do

do i=1,6
    deep(i)=0.
    deep0(i)=0.
end do

return
end
!****************************************************************************
!***** subroutine find*******************************************************
!****************************************************************************
subroutine find
use prescribed
 implicit none
!-----—————————————————————————————————————————————————————————————————————
integer::i, j, k, l, ij, jj, nb, m, n, n3
double precision::dex, dey, dez, dgyz,dgzx, dgxy, xx
double precision,dimension(3)::duu
double precision,dimension(24)::fp
double precision,dimension(mel,6)::sw
double precision,dimension(6)::de
!-----—————————————————————————————————————————————————————————————————————
do n3=1,nod3
    f(n3)=0.
end do
!-----———————————————————————————————————————————————————————————————————————
!-----———————————————————————————————————————————————————————————————————————
do m=1,mel
!----- virtual stress due to spin スピンによって発生する仮想応力ω-————————————————
      sw(m,1)=2.*(s1(m,5)*dw(m,2)-s1(m,6)*dw(m,3))
      sw(m,2)=2.*(s1(m,6)*dw(m,3)-s1(m,4)*dw(m,1))
      sw(m,3)=2.*(s1(m,4)*dw(m,1)-s1(m,5)*dw(m,2))
      sw(m,4)=(s1(m,2)-s1(m,3))*dw(m,1)+s1(m,5)*dw(m,3)-s1(m,6)*dw(m,2)
      sw(m,5)=(s1(m,3)-s1(m,1))*dw(m,2)+s1(m,6)*dw(m,1)-s1(m,4)*dw(m,3)
      sw(m,6)=(s1(m,1)-s1(m,2))*dw(m,3)+s1(m,4)*dw(m,2)-s1(m,5)*dw(m,1)
!---- virtual force————————————————————————————————————————————————————————
do i=1,24
      fp(i)=0.
 do j=1,6
      fp(i)=fp(i)+(bd(m,i,j)*dep(m,j)-b(m,j,i)*(s1(m,j)+sw(m,j)))*tvol(m)
 end do
end do

do i=1,8
 do j=1,3
      n3=3*(nel(m,i)-1)+j
      ij=3*(i-1)+j
      f(n3)=f(n3)+fp(ij)
 end do
end do

end do
!-----———————————————————————————————————————————————————————————————————————
!-----——————————————————————————————————————————————————————————————————————
!-----boundary condition for nodal force————————————————————————————————————
!x方向単軸引張りの場合
  dex=deex
  dey=deep(2)-pr*(deex-deep(1))
  dez=deep(3)-pr*(deex-deep(1))
  dgyz=deep(4)
  dgzx=deep(5)
  dgxy=deep(6)
!y軸方向の場合
!    dex=deep(1)-pr*(deey-deep(2))
!    dey=deey
!    dez=deep(3)-pr*(deey-deep(2))
!    dgxy=deep(4)
!    dgyz=deep(5)
!    dgzx=deep(6)
!z軸方向の場合
!   dex=deep(1)-pr*(deez-deep(3))
!   dey=deep(2)-pr*(deez-deep(3))
!   dez=deez
!   dgyz=deep(4)
!   dgzx=deep(5)
!   dgxy=deep(6)
!——————————————————————————————————————————————————————————————————————————
!——————————————————————————————————————————————————————————————————————————
do k=1,memz1
  do j=1,memy1
    do i=1,memx1

    nb=(i-1)*(j-1)*(k-1)*(i-memx1)*(j-memy1)*(k-memz1)

!----境界条件設定-----------------------
if (nb.eq.0) then
n=i+(j-1)*memx1+(k-1)*memx1*memy1
!x軸方向 単軸引張の場合の境界条件
  duu(1)=dex*x(n)+1*dgxy*y(n)+1*dgzx*z(n)
  duu(2)=dey*y(n)+1*dgyz*z(n)/2.
  duu(3)=dez*z(n)+1*dgyz*y(n)/2.
!y軸方向 単軸引張の場合の境界
!   duu(1)=dex*x(n)+dgzx*z(n)/2.
!   duu(2)=dgxy*x(n)+dey*y(n)+dgyz*z(n)
!   duu(3)=dez*z(n)+dgzx*x(n)/2.
!z軸方向　単軸引張りの場合の境界
!   duu(1)=dex*x(n)+dgxy*y(n)/2.
!   duu(2)=dgxy*x(n)/2.+dey*y(n)
!   duu(3)=dgzx*x(n)+dgyz*y(n)+dez*z(n)
!z方向に無限長の平面ひずみ試験
!   duu(1)=dex*x(n)+dgxy*y(n)
!   duu(2)=dey*y(n)
!   duu(3)=0.
!-----おわり---------------------------
 do jj=1,3
    n3=3*(n-1)+jj
    f(n3)=duu(jj)*1.e20
 end do
endif
  end do
 end do
end do
!----- solution of the simultaneous equations——————————————————————————————
call solve
call yieeld
!-----—————————————————————————————————————————————————————————————————————
do m=1,mel
!-----strain and spin increments———————————————————————————————————————————
do k=1,6
      de(k)=0.
  do i=1,8
    do j=1,3
      de(k)=de(k)+b(m,k,3*(i-1)+j)*du(3*(nel(m,i)-1)+j)
    end do
  end do
end do
do k=1,3
      dw(m,k)=0.
  do i=1,8
    do j=1,3
      dw(m,k)=dw(m,k)+c(m,k,3*(i-1)+j)*du(3*(nel(m,i)-1)+j)
   end do
 end do
end do
!----- stress——————————————————————————————————————————————————————————————————
do i=1,6
      s(m,i)=s1(m,i)+sw(m,i)
 do j=1,6
      s(m,i)=s(m,i)+d(i,j)*(de(j)-dep(m,j))
 end do
end do
!------------------------------------------------——————————————————————————————
do l=1,12
    tau(m,l)=0.
  do i=1,6
    tau(m,l)=tau(m,l)+eschmid(m,i,l)*s(m,i)
  end do
!-----work-hardening function——————————————————————————————————————————————————
!    yk=yield(m,l)	!modified
!-----—————————————————————————————————————————————————————————————————————————
!----- successive accumulation method—————————————————————————————————————————
    xx=abs(dg(m,l))+(abs(tau(m,l))-yield(m,l))*dt
    if(xx<0.) xx=0.
    dg(m,l)=sign(xx,tau(m,l))
end do
!---- plastic strain and spin—————————————————————————————————————————————————
do i=1,6
    dep(m,i)=0.
 do l=1,12
    dep(m,i)=dep(m,i)+eschmid(m,i,l)*dg(m,l)
 end do
end do
do k=1,3
    dwp(m,k)=0.
  do l=1,12
    dwp(m,k)=dwp(m,k)+wschmid(m,k,l)*dg(m,l)
 end do
end do

end do

!----- macroscopic stress and strain——————————————————————————————————————————
do i=1,6
    deep(i)=0.
    ss(i)=0.
  do m=1,mel
    deep(i)=deep(i)+dep(m,i)
    ss(i)=ss(i)+s(m,i)
  end do
    deep(i)=deep(i)/mel
    ss(i)=ss(i)/mel
end do
!----- variation due to iteration—————————————————————————————————————————————
    error=0.
do i=1,6
    error=error+abs(deep(i)-deep0(i))
    deep0(i)=deep(i)
end do
    ex=ex+dex
return
end
!********************************************************************************
!***** subroutine renew**********************************************************
!********************************************************************************
subroutine renew
use prescribed
 implicit none
!-----——————————————————————————————————————————————————————————————————————
integer::i, k, l, m, n
double precision::dphi, dpsi, dtheta
double precision, dimension(3)::dwr
!-----——————————————————————————————————————————————————————————————————————
!---- renewal of coordinates————————————————————————————————————————————————
do n=1,nod
    x(n)=x(n)+du(3*n-2)
    y(n)=y(n)+du(3*n-1)
    z(n)=z(n)+du(3*n)
end do
!----- accumulation of plastic strain———————————————————————————————————————
do i=1,6
    eep(i)=eep(i)+deep(i)
end do
!-----——————————————————————————————————————————————————————————————————————
do m=1,mel
 do i=1,6
    s1(m,i)=s(m,i)
 end do
 do l=1,12
    gamma_crystal(m)=gamma_crystal(m)+abs(dg(m,l))	!absは絶対値
    dgg(m,l)=dgg(m,l)+dg(m,l)
    gammaa(m,l)=gammaa(m,l)+abs(dg(m,l))
 end do
!----- lattice rotation
 do i=1,3
    dwr(i)=dw(m,i)-dwp(m,i)
 end do
    dphi=-dwr(1)*sin(theta(m))+dwr(2)*cos(theta(m))
  if(abs(phi(m))<1.e-10) phi(m)=phi(m)+dphi
    dpsi=(dwr(1)*cos(theta(m))+dwr(2)*sin(theta(m)))/sin(phi(m))
    dtheta=-dpsi*cos(phi(m))+dwr(3)
    theta(m)=theta(m)+dtheta
    phi(m)=phi(m)+dphi
    psi(m)=psi(m)+dpsi

    xc(m)=(x(nel(m,1))+x(nel(m,2))+x(nel(m,3))+x(nel(m,4))+x(nel(m,5))+x(nel(m,6))+x(nel(m,7))+x(nel(m,8)))/8.
    yc(m)=(y(nel(m,1))+y(nel(m,2))+y(nel(m,3))+y(nel(m,4))+y(nel(m,5))+y(nel(m,6))+y(nel(m,7))+y(nel(m,8)))/8.
    zc(m)=(z(nel(m,1))+z(nel(m,2))+z(nel(m,3))+z(nel(m,4))+z(nel(m,5))+z(nel(m,6))+z(nel(m,7))+z(nel(m,8)))/8.
    r(1,1)=cos(phi(m))*cos(theta(m))*cos(psi(m))-sin(theta(m))*sin(psi(m))
    r(1,2)=cos(phi(m))*sin(theta(m))*cos(psi(m))+cos(theta(m))*sin(psi(m))
    r(1,3)=-sin(phi(m))*cos(psi(m))
    r(2,1)=-cos(phi(m))*cos(theta(m))*sin(psi(m))-sin(theta(m))*cos(psi(m))
    r(2,2)=-cos(phi(m))*sin(theta(m))*sin(psi(m))+cos(theta(m))*cos(psi(m))
    r(2,3)=sin(phi(m))*sin(psi(m))
    r(3,1)=sin(phi(m))*cos(theta(m))
    r(3,2)=sin(phi(m))*sin(theta(m))
    r(3,3)=cos(phi(m))
    do n=1,8
     xr(m,n)=x(nel(m,n))-xc(m)
     yr(m,n)=y(nel(m,n))-yc(m)
     zr(m,n)=z(nel(m,n))-zc(m)
     xe(m,n)=xc(m)+r(1,1)*xr(m,n)+r(1,2)*yr(m,n)+r(1,3)*zr(m,n)
     ye(m,n)=yc(m)+r(2,1)*xr(m,n)+r(2,2)*yr(m,n)+r(2,3)*zr(m,n)
     ze(m,n)=zc(m)+r(3,1)*xr(m,n)+r(3,2)*yr(m,n)+r(3,3)*zr(m,n)
    enddo
    xv(m,1)=xe(m,4)-xe(m,5)
    xv(m,2)=xe(m,5)-xe(m,2)
    xv(m,3)=xe(m,2)-xe(m,4)
    xv(m,4)=xe(m,2)-xe(m,7)
    xv(m,5)=xe(m,7)-xe(m,4)
    xv(m,6)=xe(m,4)-xe(m,2)
    xv(m,7)=xe(m,3)-xe(m,6)
    xv(m,8)=xe(m,6)-xe(m,1)
    xv(m,9)=xe(m,1)-xe(m,3)
    xv(m,10)=xe(m,1)-xe(m,8)
    xv(m,11)=xe(m,8)-xe(m,3)
    xv(m,12)=xe(m,3)-xe(m,1)
    yv(m,1)=ye(m,4)-ye(m,5)
    yv(m,2)=ye(m,5)-ye(m,2)
    yv(m,3)=ye(m,2)-ye(m,4)
    yv(m,4)=ye(m,2)-ye(m,7)
    yv(m,5)=ye(m,7)-ye(m,4)
    yv(m,6)=ye(m,4)-ye(m,2)
    yv(m,7)=ye(m,3)-ye(m,6)
    yv(m,8)=ye(m,6)-ye(m,1)
    yv(m,9)=ye(m,1)-ye(m,3)
    yv(m,10)=ye(m,1)-ye(m,8)
    yv(m,11)=ye(m,8)-ye(m,3)
    yv(m,12)=ye(m,3)-ye(m,1)
    zv(m,1)=ze(m,4)-ze(m,5)
    zv(m,2)=ze(m,5)-ze(m,2)
    zv(m,3)=ze(m,2)-ze(m,4)
    zv(m,4)=ze(m,2)-ze(m,7)
    zv(m,5)=ze(m,7)-ze(m,4)
    zv(m,6)=ze(m,4)-ze(m,2)
    zv(m,7)=ze(m,3)-ze(m,6)
    zv(m,8)=ze(m,6)-ze(m,1)
    zv(m,9)=ze(m,1)-ze(m,3)
    zv(m,10)=ze(m,1)-ze(m,8)
    zv(m,11)=ze(m,8)-ze(m,3)
    zv(m,12)=ze(m,3)-ze(m,1)
end do
 return
end
!********************************************************************************
!**** subroutine matrix**********************************************************
!********************************************************************************
subroutine matrix
use prescribed
 implicit none
!-----———————————————————————————————————————————————————————————————————————
integer::i, j, k, mp, mq, ki, kj
double precision::sum, temp
!-----———————————————————————————————————————————————————————————————————————
do i=1,nm3
    mp=nm3-i+1
    if(mm3.lt.mp) mp=mm3
  do 4 j=1,mp
     mq=mm3-j
    if(i-1.lt.mq) mq=i-1
     sum=tk(i,j)
    if (mq.lt.1) goto 2
  do k=1,mq
    ki=i-k
    kj=j+k
    sum=sum-tk(ki,k+1)*tk(ki,kj)
  end do
2 if(j.ne.1) goto 3
    if(sum.le.0.) goto 6
      temp=1./sqrt(sum)
      tk(i,j)=temp
    goto 4
3 tk(i,j)=sum*temp
4 continue
end do

 return
6 write(*,10) i,j,tk(i,j),sum
  write(1,10) i,j,tk(i,j),sum
10 format(6h0error,5x,2hu(,i3,1h,i2,2h)=,f12.8,3x,4hsum=,f12.8)

stop
end
!********************************************************************************
!***** subroutine solve *********************************************************
!********************************************************************************
subroutine solve
use prescribed
 implicit none
!-----——————————————————————————————————————————————————————————————————————
integer::i, j, k, ii, i1, ki
double precision::sum
!-----——————————————————————————————————————————————————————————————————————
do 3 i=1,nm3
    j=i-mm3+1
  if(i+1<=mm3) j=1
    sum=f(i)
    i1=i-1
  if(i1<j) goto 2
  do 1 k=j,i1
    ki=i-k+1
    sum=sum-tk(k,ki)*du(k)
  1 continue
 2 du(i)=sum*tk(i,1)
3 continue

do 6 ii=1,nm3
    i=nm3-ii+1
    j=i+mm3-1
  if(j.gt.nm3) j=nm3
    sum=du(i)
    i1=i+1
  if(i1.gt.j) goto 5
  do 4 k=i1,j
    ki=k-i+1
    sum=sum-tk(i,ki)*du(k)
  4 continue
 5 du(i)=sum*tk(i,1)
6 continue

return
end
!***end of program**********************************************************

subroutine yieeld
  use prescribed
  implicit none
  integer::m,n,l
  double precision::gam,a
  a=0
  do m=1,mel
    gam=gamma_crystal(m)
    do n=1,12
      gam=gam+abs(dg(m,n))
    enddo
      !yk=yk0*(gam)**hard
      yk=s0*(a+gam/taylor)**hard/taylor
      if(yk<ykmin) yk=ykmin
     do l=1,12 
      yield(m,l)=yk
     enddo 
  enddo
end

subroutine create_folder(folder)
  character(len=*), intent(in) :: folder
  character(len=100) :: command

  ! フォルダが存在しなければ作成
  command = 'mkdir -p ' // trim(folder)
  call system(command)

end subroutine create_folder
