      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
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
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
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
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
c     常量声明
      real*8,parameter:: x0=1.0d-30,PRECISION_GAMMA=5.0d-3,
     1  PRECISION_RHO=5.0d-3,DPEEQ_SHRINK=1.0d-5
      parameter(ZERO=0.0d0,ONE=1.0d0,TWO=2.0d0,THREE=3.0d0,HALF=0.5d0,
     1  SIX=6.0d0,MAXSTEP=50)
C
      real*8 E,mu,Rho_c0,Rho_w0,diameter,dpeeq,peeq,Rho,rho_w,rho_c,
     1  rhoc_old,rhow_old,pre_dpeeq,dpeeq_trial,volFrac,K_frac
      
c     状态变量内容
c     state(*,1) 屈服应力
c     state(*,2) 等效塑性应变
c     state(*,3) 位错壁密度
c     state(*,4) 位错胞密度
c     state(*,5) 位错胞直径
c     state(*,6) 位错壁体积分数
      
c     state(*,7) 总位错密度
      
c     用户输入参数赋予
      E = props(1) !弹性模量
      Mu = props(2) !泊松比
      
      Sig1 = props(3) !初始应力S1
      Alpha_star = props(4) !α*
      Beta_star = props(5) !β*
      K_c = props(6) !kc
      K_w = props(7) !kw
      N_c = props(8) !nc
      N_w = props(9) !nw
      M_star = props(10) !m*
      Burger = props(11) !b 伯氏矢量
      Taylor = props(12) !M 泰勒因子
      Eta = props(13) !η 常数
      volFrac_zero = props(14) !f_0 体积分数参数
      volFrac_infin = props(15) !f_infin 体积分数参数
      Gamma0 = props(16) !参考剪应变率
      K_zero = props(17) !K参数
      K_infin = props(18) !K参数
      Peeq_wave = props(19) !常数
      Beta = props(20) !常数
      Rho_c0 = props(21) !初始位错胞密度
      Rho_w0 = props(22) !初始位错壁密度
      
      

c     拉梅常数
      La = (Mu*E)/((ONE+Mu)*(ONE-TWO*Mu)) !lambda
      Lb = E/(TWO*(ONE+Mu)) !G
      
      if(stepTime .eq. ZERO)then !初始时间步
          do k=1,nblock
              trace = strainInc(k,1)+strainInc(k,2)+strainInc(k,3)
              
              stressNew(k,1) = TWO*Lb*strainInc(k,1) + trace*La
              stressNew(k,2) = TWO*Lb*strainInc(k,2) + trace*La
              stressNew(k,3) = TWO*Lb*strainInc(k,3) + trace*La
              stressNew(k,4) = TWO*Lb*strainInc(k,4)
              stressNew(k,5) = TWO*Lb*strainInc(k,5)
              stressNew(k,6) = TWO*Lb*strainInc(k,6)
          end do
      else !非初始时间步
          do k=1,nblock !计算试应力
              trace = strainInc(k,1)+strainInc(k,2)+strainInc(k,3)
              
              s1 = stressOld(k,1) + TWO*Lb*strainInc(k,1) + trace*La
              s2 = stressOld(k,2) + TWO*Lb*strainInc(k,2) + trace*La
              s3 = stressOld(k,3) + TWO*Lb*strainInc(k,3) + trace*La
              s4 = stressOld(k,4) + TWO*Lb*strainInc(k,4)
              s5 = stressOld(k,5) + TWO*Lb*strainInc(k,5)
              s6 = stressOld(k,6) + TWO*Lb*strainInc(k,6)
              
              !计算mise应力，判断是否屈服
              sm = (s1+s2+s3)/THREE
              s1 = s1 - sm
              s2 = s2 - sm
              s3 = s3 - sm
              
              mise = sqrt(1.5d0 * 
     1  ((s1*s1 + s2*s2 + s3*s3) + TWO*(s4*s4 + s5*s5 +s6*s6)))
              
              !传递量
              sig_y = stateOld(k,1)
              peeq = stateOld(k,2)
              rho_w = stateOld(k,3)
              rho_c = stateOld(k,4)
              diameter = stateOld(k,5)
              volFrac = stateOld(k,6)
              
              !不传递量
              Rho = StateOld(k,7)
              


              dpeeq = ZERO
              
              if(sig_y .EQ. ZERO)then !在未屈服时为相关变量赋值
                  peeq = ZERO
                  
                  rho_w = rho_w0
                  rho_c = rho_c0
                  
                  !相关变量初值计算
                  K_frac = K_zero
                  
                  volFrac = volFrac_zero
                  
                  Rho = volFrac*rho_w + (ONE-volFrac)*rho_c
                  
                  diameter = K_frac/sqrt(Rho)
                  
                  sig_y = sig1 + Taylor*Eta*Lb*Burger*
     1  (volFrac*sqrt(rho_w)+(ONE-volFrac)*sqrt(rho_c))

              end if
              
              if(mise .LE. sig_y)then !未屈服 试探应力就是实际应力
                  !应力计算
                  stressNew(k,1) = s1 + sm
                  stressNew(k,2) = s2 + sm
                  stressNew(k,3) = s3 + sm
                  stressNew(k,4) = s4
                  stressNew(k,5) = s5
                  stressNew(k,6) = s6
                  
                  !能量计算
                  strainEnergy = HALF * (
     1  (stressNew(k,1)+stressOld(k,1))*strainInc(k,1) + 
     2  (stressNew(k,2)+stressOld(k,2))*strainInc(k,2) + 
     3  (stressNew(k,3)+stressOld(k,3))*strainInc(k,3) ) + 
     4  (stressNew(k,4)+stressOld(k,4))*strainInc(k,4) + 
     5  (stressNew(k,5)+stressOld(k,5))*strainInc(k,5) + 
     6  (stressNew(k,6)+stressOld(k,6))*strainInc(k,6)
                  
                  enerInternNew(k) = enerInternOld(k) + 
     1  strainEnergy/density(k)
                  
                  enerInelasNew(k) = enerInelasOld(k)
              else !屈服 按位错密度本构模型计算
                  !计算试探等效塑性应变初值
                  
                  rhoc_old = rho_c
                  rhow_old = rho_w
                  d_old = diameter
                  
                  
          !迭代计算增量步开始时的等效应变增量(按牛顿迭代法求解)
          !该增量步作为为初始试探值

          !迭代常数声明
          jstep = ZERO
          dpeeq = (mise-Sig1)/(THREE*Lb) !由mise控制的假设迭代初值
          
          !计算迭代初值
          f = ONE
          do while(f .GT. ZERO)
              dpeeq = dpeeq*DPEEQ_SHRINK
              f = Sig1-mise+THREE*Lb*dpeeq+Taylor*Eta*Lb*Burger*
     1  (volFrac*dsqrt(rho_w)+(ONE-volFrac)*dsqrt(rho_c))*
     2  (((Taylor*dpeeq)/(Gamma0*dt))**(ONE/M_star))
          end do

          pre_dpeeq = -dpeeq
      !开始用牛顿迭代计算等效应变增量
      do while(abs((pre_dpeeq-dpeeq)/pre_dpeeq) .GT. PRECISION_GAMMA)
                  
          f = Sig1-mise+THREE*Lb*dpeeq+Taylor*Eta*Lb*Burger*
     1  (volFrac*dsqrt(rho_w)+(ONE-volFrac)*dsqrt(rho_c))*
     2  (((Taylor*dpeeq)/(Gamma0*dt))**(ONE/M_star))
          
          df = THREE*Lb + Taylor*Eta*Lb*Burger*
     1  (volFrac*dsqrt(rho_w)+(ONE-volFrac)*dsqrt(rho_c))*
     2  (((Taylor*dpeeq)/(Gamma0*dt))**((ONE/M_star)-ONE))*
     3  (Taylor/(Gamma0*dt))*(ONE/M_star)
      
          pre_dpeeq = dpeeq
          dpeeq = dpeeq - f/df
          
          !计数器
          jstep = jstep + ONE

          if(jstep .GT. MAXSTEP)then
              EXIT
          end if
      end do
      


          dpeeq_trial = -dpeeq

      !迭代计算试探位错密度与屈服应力
      !设置最多迭代步数，以防迭代不收敛时卡死程序
      kstep = ZERO      
      do while(kstep .LT. MAXSTEP)!迭代循环
          dpeeq_trail = dpeeq
          
          !体积分数
          volFrac = volFrac_infin + (volFrac_zero-volFrac_infin)*
     1  exp((-Taylor*(peeq+dpeeq))/Peeq_wave)
          
          !位错胞密度
          rho_c = rhoc_old + (
     1  dsqrt(rhow_old/THREE)*((Taylor*Alpha_star)/Burger) - 
     2  (SIX*Beta_star*Taylor)/
     3  (Burger*d_old*((ONE-volFrac)**(ONE/THREE))) - 
     4  K_c*Taylor*rhoc_old*(((Taylor*dpeeq)/(dt*Gamma0))**(-ONE/N_c))
     5  ) * dpeeq
          
          !位错壁密度
          rho_w = rhow_old + (
     1  (dsqrt(THREE*rhow_old)*Beta_star*Taylor*(ONE-volFrac))/
     2  (volFrac*Burger) + 
     3  (SIX*Beta_star*Taylor*((ONE-volFrac)**(TWO/THREE)))/
     4  (Burger*d_old*volFrac) - 
     5  K_w*Taylor*rhow_old*(((Taylor*dpeeq)/(dt*Gamma0))**(-ONE/N_w)) 
     6  ) * dpeeq
          
          !总位错密度
          Rho = volFrac*rho_w + (ONE-volFrac)*rho_c
          
          !晶胞直径表达式中的分子
          K_frac = K_infin + (K_zero-K_infin)*
     1  exp(-Beta*Taylor*(peeq+dpeeq))
          
          !晶胞直径
          diameter = K_frac/dsqrt(Rho)
          
          !计算位错更新后的等效塑性应变增量
          !迭代常数声明
          jstep = ZERO
          dpeeq = (mise-Sig1)/(THREE*Lb) !由mise控制的假设迭代初值
          
          !计算迭代初值
          f = ONE
          do while(f .GT. ZERO)
              dpeeq = dpeeq*DPEEQ_SHRINK
              f = Sig1-mise+THREE*Lb*dpeeq+Taylor*Eta*Lb*Burger*
     1  (volFrac*dsqrt(rho_w)+(ONE-volFrac)*dsqrt(rho_c))*
     2  (((Taylor*dpeeq)/(Gamma0*dt))**(ONE/M_star))
          end do
          
          pre_dpeeq = -dpeeq
      !开始用牛顿迭代计算等效应变增量
      do while(abs((pre_dpeeq-dpeeq)/pre_dpeeq) .GT. PRECISION_GAMMA)
                  
          f = Sig1-mise+THREE*Lb*dpeeq+Taylor*Eta*Lb*Burger*
     1  (volFrac*dsqrt(rho_w)+(ONE-volFrac)*dsqrt(rho_c))*
     2  (((Taylor*dpeeq)/(Gamma0*dt))**(ONE/M_star))
          
          df = THREE*Lb + Taylor*Eta*Lb*Burger*
     1  (volFrac*dsqrt(rho_w)+(ONE-volFrac)*dsqrt(rho_c))*
     2  (((Taylor*dpeeq)/(Gamma0*dt))**((ONE/M_star)-ONE))*
     3  (Taylor/(Gamma0*dt))*(ONE/M_star)
      
          pre_dpeeq = dpeeq
          dpeeq = dpeeq - f/df
          
          !计数器
          jstep = jstep + ONE

          if(jstep .GT. MAXSTEP)then
              EXIT
          end if
      end do
          
          
          
          if(abs((dpeeq_trial-dpeeq)/dpeeq_trial)
     1  .GT. PRECISION_GAMMA)then
              dpeeq_trial = dpeeq
          else
              EXIT
          end if
          
      !计数器
      kstep = kstep + ONE

      end do!迭代循环结束
      
      
      
      !屈服应力
      
      if(dpeeq .EQ. ZERO)then
          sig_y = mise
      else
          sig_y = Sig1 + (Taylor*Eta*Lb*Burger)*
     1  (volFrac*dsqrt(rho_w)+(ONE-volFrac)*dsqrt(rho_c))*
     2  (((Taylor*dpeeq)/(dt*Gamma0))**(ONE/M_star))
      end if

                  !根据迭代得出的等效塑性应变计算折减系数
                  factor = sig_y/mise
                  
                  !计算应力
                  stressNew(k,1) = s1*factor + sm
                  stressNew(k,2) = s2*factor + sm
                  stressNew(k,3) = s3*factor + sm
                  stressNew(k,4) = s4*factor
                  stressNew(k,5) = s5*factor
                  stressNew(k,6) = s6*factor
                  
                  !计算能量
                  strainEnergy = HALF * (
     1  (stressNew(k,1)+stressOld(k,1))*strainInc(k,1) + 
     2  (stressNew(k,2)+stressOld(k,2))*strainInc(k,2) + 
     3  (stressNew(k,3)+stressOld(k,3))*strainInc(k,3) ) + 
     4  (stressNew(k,4)+stressOld(k,4))*strainInc(k,4) + 
     5  (stressNew(k,5)+stressOld(k,5))*strainInc(k,5) + 
     6  (stressNew(k,6)+stressOld(k,6))*strainInc(k,6)
                  
                  plasEnergy = dpeeq*sig_y
                  
                  enerInternNew(k) = enerInternOld(k) + 
     1  strainEnergy/density(k)
                  
                  enerInelasNew(k) = enerInelasOld(k) +
     2  plasEnergy/density(k)
      
              end if
          stateNew(k,1) = sig_y
          stateNew(k,2) = peeq + dpeeq
          stateNew(k,3) = rho_w
          stateNew(k,4) = rho_c
          stateNew(k,5) = diameter
          stateNew(k,6) = volFrac
c
          stateNew(k,7) = Rho
          end do
      end if

      return
      end