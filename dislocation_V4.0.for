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
c     ��������
      real*8,parameter:: x0=1.0d-30,PRECISION_GAMMA=5.0d-3,
     1  PRECISION_RHO=5.0d-3,DPEEQ_SHRINK=1.0d-5
      parameter(ZERO=0.0d0,ONE=1.0d0,TWO=2.0d0,THREE=3.0d0,HALF=0.5d0,
     1  SIX=6.0d0,MAXSTEP=50)
C
      real*8 E,mu,Rho_c0,Rho_w0,diameter,dpeeq,peeq,Rho,rho_w,rho_c,
     1  rhoc_old,rhow_old,pre_dpeeq,dpeeq_trial,volFrac,K_frac
      
c     ״̬��������
c     state(*,1) ����Ӧ��
c     state(*,2) ��Ч����Ӧ��
c     state(*,3) λ����ܶ�
c     state(*,4) λ����ܶ�
c     state(*,5) λ���ֱ��
c     state(*,6) λ����������
      
c     state(*,7) ��λ���ܶ�
      
c     �û������������
      E = props(1) !����ģ��
      Mu = props(2) !���ɱ�
      
      Sig1 = props(3) !��ʼӦ��S1
      Alpha_star = props(4) !��*
      Beta_star = props(5) !��*
      K_c = props(6) !kc
      K_w = props(7) !kw
      N_c = props(8) !nc
      N_w = props(9) !nw
      M_star = props(10) !m*
      Burger = props(11) !b ����ʸ��
      Taylor = props(12) !M ̩������
      Eta = props(13) !�� ����
      volFrac_zero = props(14) !f_0 �����������
      volFrac_infin = props(15) !f_infin �����������
      Gamma0 = props(16) !�ο���Ӧ����
      K_zero = props(17) !K����
      K_infin = props(18) !K����
      Peeq_wave = props(19) !����
      Beta = props(20) !����
      Rho_c0 = props(21) !��ʼλ����ܶ�
      Rho_w0 = props(22) !��ʼλ����ܶ�
      
      

c     ��÷����
      La = (Mu*E)/((ONE+Mu)*(ONE-TWO*Mu)) !lambda
      Lb = E/(TWO*(ONE+Mu)) !G
      
      if(stepTime .eq. ZERO)then !��ʼʱ�䲽
          do k=1,nblock
              trace = strainInc(k,1)+strainInc(k,2)+strainInc(k,3)
              
              stressNew(k,1) = TWO*Lb*strainInc(k,1) + trace*La
              stressNew(k,2) = TWO*Lb*strainInc(k,2) + trace*La
              stressNew(k,3) = TWO*Lb*strainInc(k,3) + trace*La
              stressNew(k,4) = TWO*Lb*strainInc(k,4)
              stressNew(k,5) = TWO*Lb*strainInc(k,5)
              stressNew(k,6) = TWO*Lb*strainInc(k,6)
          end do
      else !�ǳ�ʼʱ�䲽
          do k=1,nblock !������Ӧ��
              trace = strainInc(k,1)+strainInc(k,2)+strainInc(k,3)
              
              s1 = stressOld(k,1) + TWO*Lb*strainInc(k,1) + trace*La
              s2 = stressOld(k,2) + TWO*Lb*strainInc(k,2) + trace*La
              s3 = stressOld(k,3) + TWO*Lb*strainInc(k,3) + trace*La
              s4 = stressOld(k,4) + TWO*Lb*strainInc(k,4)
              s5 = stressOld(k,5) + TWO*Lb*strainInc(k,5)
              s6 = stressOld(k,6) + TWO*Lb*strainInc(k,6)
              
              !����miseӦ�����ж��Ƿ�����
              sm = (s1+s2+s3)/THREE
              s1 = s1 - sm
              s2 = s2 - sm
              s3 = s3 - sm
              
              mise = sqrt(1.5d0 * 
     1  ((s1*s1 + s2*s2 + s3*s3) + TWO*(s4*s4 + s5*s5 +s6*s6)))
              
              !������
              sig_y = stateOld(k,1)
              peeq = stateOld(k,2)
              rho_w = stateOld(k,3)
              rho_c = stateOld(k,4)
              diameter = stateOld(k,5)
              volFrac = stateOld(k,6)
              
              !��������
              Rho = StateOld(k,7)
              


              dpeeq = ZERO
              
              if(sig_y .EQ. ZERO)then !��δ����ʱΪ��ر�����ֵ
                  peeq = ZERO
                  
                  rho_w = rho_w0
                  rho_c = rho_c0
                  
                  !��ر�����ֵ����
                  K_frac = K_zero
                  
                  volFrac = volFrac_zero
                  
                  Rho = volFrac*rho_w + (ONE-volFrac)*rho_c
                  
                  diameter = K_frac/sqrt(Rho)
                  
                  sig_y = sig1 + Taylor*Eta*Lb*Burger*
     1  (volFrac*sqrt(rho_w)+(ONE-volFrac)*sqrt(rho_c))

              end if
              
              if(mise .LE. sig_y)then !δ���� ��̽Ӧ������ʵ��Ӧ��
                  !Ӧ������
                  stressNew(k,1) = s1 + sm
                  stressNew(k,2) = s2 + sm
                  stressNew(k,3) = s3 + sm
                  stressNew(k,4) = s4
                  stressNew(k,5) = s5
                  stressNew(k,6) = s6
                  
                  !��������
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
              else !���� ��λ���ܶȱ���ģ�ͼ���
                  !������̽��Ч����Ӧ���ֵ
                  
                  rhoc_old = rho_c
                  rhow_old = rho_w
                  d_old = diameter
                  
                  
          !����������������ʼʱ�ĵ�ЧӦ������(��ţ�ٵ��������)
          !����������ΪΪ��ʼ��ֵ̽

          !������������
          jstep = ZERO
          dpeeq = (mise-Sig1)/(THREE*Lb) !��mise���Ƶļ��������ֵ
          
          !���������ֵ
          f = ONE
          do while(f .GT. ZERO)
              dpeeq = dpeeq*DPEEQ_SHRINK
              f = Sig1-mise+THREE*Lb*dpeeq+Taylor*Eta*Lb*Burger*
     1  (volFrac*dsqrt(rho_w)+(ONE-volFrac)*dsqrt(rho_c))*
     2  (((Taylor*dpeeq)/(Gamma0*dt))**(ONE/M_star))
          end do

          pre_dpeeq = -dpeeq
      !��ʼ��ţ�ٵ��������ЧӦ������
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
          
          !������
          jstep = jstep + ONE

          if(jstep .GT. MAXSTEP)then
              EXIT
          end if
      end do
      


          dpeeq_trial = -dpeeq

      !����������̽λ���ܶ�������Ӧ��
      !�����������������Է�����������ʱ��������
      kstep = ZERO      
      do while(kstep .LT. MAXSTEP)!����ѭ��
          dpeeq_trail = dpeeq
          
          !�������
          volFrac = volFrac_infin + (volFrac_zero-volFrac_infin)*
     1  exp((-Taylor*(peeq+dpeeq))/Peeq_wave)
          
          !λ����ܶ�
          rho_c = rhoc_old + (
     1  dsqrt(rhow_old/THREE)*((Taylor*Alpha_star)/Burger) - 
     2  (SIX*Beta_star*Taylor)/
     3  (Burger*d_old*((ONE-volFrac)**(ONE/THREE))) - 
     4  K_c*Taylor*rhoc_old*(((Taylor*dpeeq)/(dt*Gamma0))**(-ONE/N_c))
     5  ) * dpeeq
          
          !λ����ܶ�
          rho_w = rhow_old + (
     1  (dsqrt(THREE*rhow_old)*Beta_star*Taylor*(ONE-volFrac))/
     2  (volFrac*Burger) + 
     3  (SIX*Beta_star*Taylor*((ONE-volFrac)**(TWO/THREE)))/
     4  (Burger*d_old*volFrac) - 
     5  K_w*Taylor*rhow_old*(((Taylor*dpeeq)/(dt*Gamma0))**(-ONE/N_w)) 
     6  ) * dpeeq
          
          !��λ���ܶ�
          Rho = volFrac*rho_w + (ONE-volFrac)*rho_c
          
          !����ֱ�����ʽ�еķ���
          K_frac = K_infin + (K_zero-K_infin)*
     1  exp(-Beta*Taylor*(peeq+dpeeq))
          
          !����ֱ��
          diameter = K_frac/dsqrt(Rho)
          
          !����λ����º�ĵ�Ч����Ӧ������
          !������������
          jstep = ZERO
          dpeeq = (mise-Sig1)/(THREE*Lb) !��mise���Ƶļ��������ֵ
          
          !���������ֵ
          f = ONE
          do while(f .GT. ZERO)
              dpeeq = dpeeq*DPEEQ_SHRINK
              f = Sig1-mise+THREE*Lb*dpeeq+Taylor*Eta*Lb*Burger*
     1  (volFrac*dsqrt(rho_w)+(ONE-volFrac)*dsqrt(rho_c))*
     2  (((Taylor*dpeeq)/(Gamma0*dt))**(ONE/M_star))
          end do
          
          pre_dpeeq = -dpeeq
      !��ʼ��ţ�ٵ��������ЧӦ������
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
          
          !������
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
          
      !������
      kstep = kstep + ONE

      end do!����ѭ������
      
      
      
      !����Ӧ��
      
      if(dpeeq .EQ. ZERO)then
          sig_y = mise
      else
          sig_y = Sig1 + (Taylor*Eta*Lb*Burger)*
     1  (volFrac*dsqrt(rho_w)+(ONE-volFrac)*dsqrt(rho_c))*
     2  (((Taylor*dpeeq)/(dt*Gamma0))**(ONE/M_star))
      end if

                  !���ݵ����ó��ĵ�Ч����Ӧ������ۼ�ϵ��
                  factor = sig_y/mise
                  
                  !����Ӧ��
                  stressNew(k,1) = s1*factor + sm
                  stressNew(k,2) = s2*factor + sm
                  stressNew(k,3) = s3*factor + sm
                  stressNew(k,4) = s4*factor
                  stressNew(k,5) = s5*factor
                  stressNew(k,6) = s6*factor
                  
                  !��������
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