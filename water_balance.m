%%Clinton,MA run water balance model to validate with Eagleson
%
%% define meshgrid
[M, s_o]=meshgrid(0:0.001:1,0:0.001:1);
  
%Clinton site level parameters
%will need to 
 
 m_n=109; % mean number of storms per year % obtained from Eagles_on: Climate, s_oil, and Vegetation
  e_p=0.15; % average bare s_oil potential evaporation rate
  m_Pa=102.6; %mean annual precipitation in cm
  m_r=0.32; %mean storm depth--double check value
  m_b=2.98; %in units of days
  m_t=365; %mean rainy seas_on
  alp=1./m_b; %one over mean time between storms
  %Eagles_on (2002) has poiss_on pulse parameters for a number of stations in Appendix F
  %for Clinton, MA, I used parameters for Boston (Station ID 14 in eagles_on
  lambda=1.47.*2.54; %may be incorrect, changed inches to cm %scale parameter for gamma distribution of storm depths (1./cm)
  gamma_depth=0.50; %these may be incorrect %shape parameter for the gamma distribution of storm depths
   
  
  
 
%define s_oil paramters
%'clay')
  Ksat=0.72; %k saturation
  Matrix_pot=25; %s_oil matrix potential
  c =12; %for clay
  n=0.45; %s_oil porosity, for clay
  
  h_o=0.1; % small constant value for surface water retention
  k_v=1;% unstressed transpiration
  d=4.3;
  w=0; %assume no capillary rise
      
%ET model
 gamma_ratio = gammainc(gamma_depth, lambda.*h_o)./gamma(gamma_depth);
    
% this value gives the approximate for exfiltration capacity
ex_diffus_approx=0.3477-(0.0731.*c)+(0.0062.*c.^2)-0.0002.*c.^3;
 %ex diffus over 2 meters
%ex_diffus_approx= (1+1.85*(-1./(1.85+1))
 %calclate the indefinate integral for:
% syms s_o c x
% f= x.^((c+1)./2).*(s_o-x).^(0.85)
 %int(f,x,0, s_o)
 
 %d=c+1/2
 
 %result: x.^(c/2 + 1/2).*(s_o - x).^(17/20)
 
 %s_o.^(d).*(1./d+(5./3))+
 
 
 %ex_diffus_approx=(1./d+[(5./3)-1]*(s_o./1-s_o).^1+1./d+[(5./3)-2]*(s_o./1-s_o).^2+1./d+[(5./3)-3]*(s_o./1-s_o).^3+1./d+[(5./3)-4]*(s_o./1-s_o).^4);
 
 %ex_diffus_approx=ex_diffus_approx=(1-s_o)^d(1+1.85.*((-1)^1)./(1.85+1)+((-2)^1)./(1.85)+((-3)^1)./(1.85+3)+((-4)^1)./(1.85+4))
 B= ((1-M)./(1+M.*k_v-w./e_p))+((k_v.*(M.^2)+(1-M).*w./e_p)./((2.*(1+M.*k_v-w./e_p).^2)));
  
C= 1./2.*(M.*k_v-w./e_p).^(-2);
 
E=(alp.*n.*(c-3).*Ksat.*Matrix_pot.*ex_diffus_approx./(pi.*e_p.^2)).*s_o.^(c+5)./2;
  
%beta_s 
 beta=gamma_ratio-(1+(alp.*h_o./e_p)./lambda.*h_o).^(-1.*gamma_depth).*...
     (((gammainc(gamma_depth, (lambda.*h_o + (alp.*h_o)./e_p)))./gamma(gamma_depth)).*exp(-B.*E))+...
     (1 + gamma_ratio).*(1-exp(-B.*E-(alp.*h_o./e_p)).*(1 + M.*k_v + ((2.*B).^(1./2)).*E-(w./e_p))+...
     exp(-C.*E-alp.*h_o./e_p).*((M.*k_v +(2.*C).^(1./2)).*E-(w./e_p))+...
     (((2.*E).^(1./2).*exp(-alp.*h_o./e_p)).*((gammainc(3/2, C.*E)-gammainc(3./2,B.*E)))) + ...
     ((1+((alp.*h_o./e_p)./lambda.*h_o)).^(-1.*gamma_depth)).*...
     (gammainc(gamma_depth, lambda.*h_o+(alp.*h_o./e_p))./gamma(gamma_depth)).*...
     (sqrt(2.*E).*(gammainc(3./2, C.*E)-gammainc(3./2,B.*E))+...
     exp(-C.*E).*(M.*k_v +(((2.*C).^(1./2)).*E)-w./e_p)-...
     exp(-B.*E).*(M.*k_v+(((2.*B).^(1./2)).*E)-w./e_p)));
        
  %compute ET        
ET=((m_n.*e_p)./alp).*(1-M).*beta + M.*k_v ;

%in_diffus_approx=((5./3)+(1./2).*(c+1).*(1-s_o).^(1.425-0.0375.*(c+1)./2)).^(-1);

%used integrate to calculate integration of eqn 11 in Kerkhoff 
in_diffus_approx=(1-s_o).^(5./3).*((3.*(3.*s_o + 5)).*(1 - s_o).^(5./3)./40);
% deleted this part because s_o-s_o=0 ((3*(3*s_o + 5*s_o)*(s_o - s_o)^(5/3))/40))
%%%Runoff model

%mPA mean annual precipitation in mm, not sure if this needs to be changed to cm
  G=(Ksat.*((m_n.*m_r)./m_Pa)).*(((1+s_o.^c)./2)-(w./Ksat));
  sigma=((5.*n.*(c-3).*lambda.^2.*Ksat.*Matrix_pot.*in_diffus_approx.*m_r)./12.*pi.*(gamma_depth.^2)).^(1./3).*(1-s_o).^(2./3);
  Run =m_Pa.*(exp(-G-2.*sigma).*gamma(sigma+1).*(sigma.^(-sigma))+((m_t.*Ksat)./m_Pa).*(s_o.^c));
  scaled = Run./m_Pa;
%compute water balance
    bal=(ET + Run)-m_Pa; 
  %plot contour line
  figure(1)
 [C,h]=contour(M,s_o,bal);
clabel(C,h)
title('Clinton, clay')
xlabel('M')
ylabel('s_o')
%%
%for sandyloam
  Ksat=18 %k saturation
      Matrix_pot=200 %soil matrix potential
      c=4 %for sandy loam
      n=0.25 %soil porosity, for sandyy loam 
      
      %%these values set to constant for this purpose
      h_o=0.1 % small constant value for surface water retention
      k_v=1% unstressed transpiration
      %M=0.5 % ranges from 0 to 1 (but cant be one)
      w=0 %assume no capillary rise

      gamma_ratio = gammainc(gamma_depth, lambda.*h_o)./gamma(gamma_depth);
    
% this value gives the approximate for exfiltration capacity
ex_diffus_approx=0.3477-(0.0731.*c)+(0.0062.*c.^2)-0.0002.*c.^3;
 %ex diffus over 2 meters
%ex_diffus_approx= (1+1.85*(-1./(1.85+1))
 %calclate the indefinate integral for:
% syms s_o c x
% f= x.^((c+1)./2).*(s_o-x).^(0.85)
 %int(f,x,0, s_o)
 
 %d=c+1/2
 
 %result: x.^(c/2 + 1/2).*(s_o - x).^(17/20)
 
 %s_o.^(d).*(1./d+(5./3))+
 
 
 %ex_diffus_approx=(1./d+[(5./3)-1]*(s_o./1-s_o).^1+1./d+[(5./3)-2]*(s_o./1-s_o).^2+1./d+[(5./3)-3]*(s_o./1-s_o).^3+1./d+[(5./3)-4]*(s_o./1-s_o).^4);
 
 %ex_diffus_approx=ex_diffus_approx=(1-s_o)^d(1+1.85.*((-1)^1)./(1.85+1)+((-2)^1)./(1.85)+((-3)^1)./(1.85+3)+((-4)^1)./(1.85+4))
 B= ((1-M)./(1+M.*k_v-w./e_p))+((k_v.*(M.^2)+(1-M).*w./e_p)./((2.*(1+M.*k_v-w./e_p).^2)));
  
C= 1./2.*(M.*k_v-w./e_p).^(-2);
 
E=(alp.*n.*(c-3).*Ksat.*Matrix_pot.*ex_diffus_approx./(pi.*e_p.^2)).*s_o.^(c+5)./2;
  
%beta_s 
 beta=gamma_ratio-(1+(alp.*h_o./e_p)./lambda.*h_o).^(-1.*gamma_depth).*...
     (((gammainc(gamma_depth, (lambda.*h_o + (alp.*h_o)./e_p)))./gamma(gamma_depth)).*exp(-B.*E))+...
     (1 + gamma_ratio).*(1-exp(-B.*E-(alp.*h_o./e_p)).*(1 + M.*k_v + ((2.*B).^(1./2)).*E-(w./e_p))+...
     exp(-C.*E-alp.*h_o./e_p).*((M.*k_v +(2.*C).^(1./2)).*E-(w./e_p))+...
     (((2.*E).^(1./2).*exp(-alp.*h_o./e_p)).*((gammainc(3/2, C.*E)-gammainc(3./2,B.*E)))) + ...
     ((1+((alp.*h_o./e_p)./lambda.*h_o)).^(-1.*gamma_depth)).*...
     (gammainc(gamma_depth, lambda.*h_o+(alp.*h_o./e_p))./gamma(gamma_depth)).*...
     (sqrt(2.*E).*(gammainc(3./2, C.*E)-gammainc(3./2,B.*E))+...
     exp(-C.*E).*(M.*k_v +(((2.*C).^(1./2)).*E)-w./e_p)-...
     exp(-B.*E).*(M.*k_v+(((2.*B).^(1./2)).*E)-w./e_p)));
        
  %compute ET        
ET=((m_n.*e_p)./alp).*(1-M).*beta + M.*k_v ;
%in_diffus_approx=((5./3)+(1./2).*(c+1).*(1-s_o).^(1.425-0.0375.*(c+1)./2)).^(-1);

%used integrate to calculate integration of eqn 11 in Kerkhoff 
in_diffus_approx=(1-s_o).^(5./3).*((3.*(3.*s_o + 5)).*(1 - s_o).^(5./3)./40);
% deleted this part because s_o-s_o=0 ((3*(3*s_o + 5*s_o)*(s_o - s_o)^(5/3))/40))
%%%Runoff model

%mPA mean annual precipitation in mm, not sure if this needs to be changed to cm
  G=(Ksat.*((m_n.*m_r)./m_Pa)).*(((1+s_o.^c)./2)-(w./Ksat));
  sigma=((5.*n.*(c-3).*lambda.^2.*Ksat.*Matrix_pot.*in_diffus_approx.*m_r)./12.*pi.*(gamma_depth.^2)).^(1./3).*(1-s_o).^(2./3);
  Run =m_Pa.*(exp(-G-2.*sigma).*gamma(sigma+1).*(sigma.^(-sigma))+((m_t.*Ksat)./m_Pa).*(s_o.^c));
  scaled = Run./m_Pa;
%compute water balance
    bal=(ET + Run)-m_Pa; 
  %plot contour line
  figure(2)
 [C,h]=contour(M,s_o,bal);
clabel(C,h)
title('Clinton, sandyloam')
xlabel('M')
ylabel('s_o')


%%
% for silt
Ksat=8.64 %k saturation
        Matrix_pot=166 %soil matrix potential
        c=6 %for clayloam
        n=0.35 %soil porosity, for clay
        
        %%these values set to constant for this purpose
        h_o=0.1 % small constant value for surface water retention
        k_v=1% unstressed transpiration
        %M=0.5 % ranges from 0 to 1 (but cant be one)
        w=0 %assume no capillary rise
d=2.9
gamma_ratio = gammainc(gamma_depth, lambda.*h_o)./gamma(gamma_depth);
    
% this value gives the approximate for exfiltration capacity
ex_diffus_approx=0.3477-(0.0731.*c)+(0.0062.*c.^2)-0.0002.*c.^3;
 %ex diffus over 2 meters
%ex_diffus_approx= (1+1.85*(-1./(1.85+1))
 %calclate the indefinate integral for:
% syms s_o c x
% f= x.^((c+1)./2).*(s_o-x).^(0.85)
 %int(f,x,0, s_o)
 
 %d=c+1/2
 
 %result: x.^(c/2 + 1/2).*(s_o - x).^(17/20)
 
 %s_o.^(d).*(1./d+(5./3))+
 
 
 %ex_diffus_approx=(1./d+[(5./3)-1]*(s_o./1-s_o).^1+1./d+[(5./3)-2]*(s_o./1-s_o).^2+1./d+[(5./3)-3]*(s_o./1-s_o).^3+1./d+[(5./3)-4]*(s_o./1-s_o).^4);
 
 %ex_diffus_approx=ex_diffus_approx=(1-s_o)^d(1+1.85.*((-1)^1)./(1.85+1)+((-2)^1)./(1.85)+((-3)^1)./(1.85+3)+((-4)^1)./(1.85+4))
 B= ((1-M)./(1+M.*k_v-w./e_p))+((k_v.*(M.^2)+(1-M).*w./e_p)./((2.*(1+M.*k_v-w./e_p).^2)));
  
C= 1./2.*(M.*k_v-w./e_p).^(-2);
 
E=(alp.*n.*(c-3).*Ksat.*Matrix_pot.*ex_diffus_approx./(pi.*e_p.^2)).*s_o.^(c+5)./2;
  
%beta_s 
 beta=gamma_ratio-(1+(alp.*h_o./e_p)./lambda.*h_o).^(-1.*gamma_depth).*...
     (((gammainc(gamma_depth, (lambda.*h_o + (alp.*h_o)./e_p)))./gamma(gamma_depth)).*exp(-B.*E))+...
     (1 + gamma_ratio).*(1-exp(-B.*E-(alp.*h_o./e_p)).*(1 + M.*k_v + ((2.*B).^(1./2)).*E-(w./e_p))+...
     exp(-C.*E-alp.*h_o./e_p).*((M.*k_v +(2.*C).^(1./2)).*E-(w./e_p))+...
     (((2.*E).^(1./2).*exp(-alp.*h_o./e_p)).*((gammainc(3/2, C.*E)-gammainc(3./2,B.*E)))) + ...
     ((1+((alp.*h_o./e_p)./lambda.*h_o)).^(-1.*gamma_depth)).*...
     (gammainc(gamma_depth, lambda.*h_o+(alp.*h_o./e_p))./gamma(gamma_depth)).*...
     (sqrt(2.*E).*(gammainc(3./2, C.*E)-gammainc(3./2,B.*E))+...
     exp(-C.*E).*(M.*k_v +(((2.*C).^(1./2)).*E)-w./e_p)-...
     exp(-B.*E).*(M.*k_v+(((2.*B).^(1./2)).*E)-w./e_p)));
        
  %compute ET        
ET=((m_n.*e_p)./alp).*(1-M).*beta + M.*k_v ;

%in_diffus_approx=((5./3)+(1./2).*(c+1).*(1-s_o).^(1.425-0.0375.*(c+1)./2)).^(-1);

%used integrate to calculate integration of eqn 11 in Kerkhoff 
in_diffus_approx=(1-s_o).^(5./3).*((3.*(3.*s_o + 5)).*(1 - s_o).^(5./3)./40);
% deleted this part because s_o-s_o=0 ((3*(3*s_o + 5*s_o)*(s_o - s_o)^(5/3))/40))
%%%Runoff model

%mPA mean annual precipitation in mm, not sure if this needs to be changed to cm
  G=(Ksat.*((m_n.*m_r)./m_Pa)).*(((1+s_o.^c)./2)-(w./Ksat));
  sigma=((5.*n.*(c-3).*lambda.^2.*Ksat.*Matrix_pot.*in_diffus_approx.*m_r)./12.*pi.*(gamma_depth.^2)).^(1./3).*(1-s_o).^(2./3);
  Run =m_Pa.*(exp(-G-2.*sigma).*gamma(sigma+1).*(sigma.^(-sigma))+((m_t.*Ksat)./m_Pa).*(s_o.^c));
  scaled = Run./m_Pa;
  
%compute water balance
    bal=(ET + Run)-m_Pa; 
  %plot contour line
  figure(3)
 [C,h]=contour(M,s_o,bal);
clabel(C,h)
title('Clinton, silt')
xlabel('M')
ylabel('s_o')
%%

%for clayloam
%%
% 
Ksat=2.0 %k saturation
        Matrix_pot=19 %soil matrix potential
        c=10 %for clayloam
        n=0.35 %soil porosity, for clay
        d=4.3
        %%these values set to constant for this purpose
        h_o=0.1 % small constant value for surface water retention
        k_v=1% unstressed transpiration
        %M=0.5 % ranges from 0 to 1 (but cant be one)
        w=0 %assume no capillary rise
d=2.9
gamma_ratio = gammainc(gamma_depth, lambda.*h_o)./gamma(gamma_depth);
    
% this value gives the approximate for exfiltration capacity
ex_diffus_approx=0.3477-(0.0731.*c)+(0.0062.*c.^2)-0.0002.*c.^3;
 %ex diffus over 2 meters
%ex_diffus_approx= (1+1.85*(-1./(1.85+1))
 %calclate the indefinate integral for:
% syms s_o c x
% f= x.^((c+1)./2).*(s_o-x).^(0.85)
 %int(f,x,0, s_o)
 
 %d=c+1/2
 
 %result: x.^(c/2 + 1/2).*(s_o - x).^(17/20)
 
 %s_o.^(d).*(1./d+(5./3))+
 
 
 %ex_diffus_approx=(1./d+[(5./3)-1]*(s_o./1-s_o).^1+1./d+[(5./3)-2]*(s_o./1-s_o).^2+1./d+[(5./3)-3]*(s_o./1-s_o).^3+1./d+[(5./3)-4]*(s_o./1-s_o).^4);
 
 %ex_diffus_approx=ex_diffus_approx=(1-s_o)^d(1+1.85.*((-1)^1)./(1.85+1)+((-2)^1)./(1.85)+((-3)^1)./(1.85+3)+((-4)^1)./(1.85+4))
 B= ((1-M)./(1+M.*k_v-w./e_p))+((k_v.*(M.^2)+(1-M).*w./e_p)./((2.*(1+M.*k_v-w./e_p).^2)));
  
C= 1./2.*(M.*k_v-w./e_p).^(-2);
 
E=(alp.*n.*(c-3).*Ksat.*Matrix_pot.*ex_diffus_approx./(pi.*e_p.^2)).*s_o.^(c+5)./2;
  
%beta_s 
 beta=gamma_ratio-(1+(alp.*h_o./e_p)./lambda.*h_o).^(-1.*gamma_depth).*...
     (((gammainc(gamma_depth, (lambda.*h_o + (alp.*h_o)./e_p)))./gamma(gamma_depth)).*exp(-B.*E))+...
     (1 + gamma_ratio).*(1-exp(-B.*E-(alp.*h_o./e_p)).*(1 + M.*k_v + ((2.*B).^(1./2)).*E-(w./e_p))+...
     exp(-C.*E-alp.*h_o./e_p).*((M.*k_v +(2.*C).^(1./2)).*E-(w./e_p))+...
     (((2.*E).^(1./2).*exp(-alp.*h_o./e_p)).*((gammainc(3/2, C.*E)-gammainc(3./2,B.*E)))) + ...
     ((1+((alp.*h_o./e_p)./lambda.*h_o)).^(-1.*gamma_depth)).*...
     (gammainc(gamma_depth, lambda.*h_o+(alp.*h_o./e_p))./gamma(gamma_depth)).*...
     (sqrt(2.*E).*(gammainc(3./2, C.*E)-gammainc(3./2,B.*E))+...
     exp(-C.*E).*(M.*k_v +(((2.*C).^(1./2)).*E)-w./e_p)-...
     exp(-B.*E).*(M.*k_v+(((2.*B).^(1./2)).*E)-w./e_p)));
        
  %compute ET        
ET=((m_n.*e_p)./alp).*(1-M).*beta + M.*k_v ;

%in_diffus_approx=((5./3)+(1./2).*(c+1).*(1-s_o).^(1.425-0.0375.*(c+1)./2)).^(-1);

%used integrate to calculate integration of eqn 11 in Kerkhoff 
in_diffus_approx=(1-s_o).^(5./3).*((3.*(3.*s_o + 5)).*(1 - s_o).^(5./3)./40);
% deleted this part because s_o-s_o=0 ((3*(3*s_o + 5*s_o)*(s_o - s_o)^(5/3))/40))
%%%Runoff model

%mPA mean annual precipitation in mm, not sure if this needs to be changed to cm
  G=(Ksat.*((m_n.*m_r)./m_Pa)).*(((1+s_o.^c)./2)-(w./Ksat));
  sigma=((5.*n.*(c-3).*lambda.^2.*Ksat.*Matrix_pot.*in_diffus_approx.*m_r)./12.*pi.*(gamma_depth.^2)).^(1./3).*(1-s_o).^(2./3);
  Run =m_Pa.*(exp(-G-2.*sigma).*gamma(sigma+1).*(sigma.^(-sigma))+((m_t.*Ksat)./m_Pa).*(s_o.^c));
  scaled = Run./m_Pa;
  
%compute water balance
    bal=(ET + Run)-m_Pa; 
  %plot contour line
  figure(4)
 [C,h]=contour(M,s_o,bal);
clabel(C,h)
title('Clinton, clayloam')
xlabel('M')
ylabel('s_o')
%% for SantaPaula
 %SantaPaula
  m_n=15.7; % mean number of storms per year % obtained from Eagleson: Climate, Soil, and Vegetation
      e_p=0.27; % average bare soil potential evaporation rate
      m_Pa=54.4; %mean annual precipitation in cm
      m_r=1.4; %mean storm depth--double check value
      m_b=10.42; %in units of days
      m_t=212; %mean rainy season
     
      alp=1/m_b; %one over mean time between storms
      %Eagleson (2002) has poisson pulse parameters for a number of stations in Appendix F
      %for Clinton, MA, I used parameters for Boston (Station ID 14 in eagleson
      gamma_depth=0.37; 
      lambda=gamma_depth/(m_Pa/m_n); %
      %lambda= 1.86*2.54;
        %c=4.95
        %define s_oil paramters
%'clay')
  Ksat=0.72; %k saturation
  Matrix_pot=25; %s_oil matrix potential
  c =12; %for clay
  n=0.45; %s_oil porosity, for clay
  
  h_o=0.1; % small constant value for surface water retention
  k_v=1;% unstressed transpiration
  d=4.3;
  w=0; %assume no capillary rise
      
%ET model
 gamma_ratio = gammainc(gamma_depth, lambda.*h_o)./gamma(gamma_depth);
    
 
% this value gives the approximate for exfiltration capacity
%fun = @(x,c,s_o) (x.^(c/2+1/2)).*(s_o-x).^(0.85);


%for i=1:1:s_o;
%q(i)=integral(@(x)fun(x,c,s_o(i)),0,s_o(i));
%end

%need to get q to match s_o
 %q2=repmat([q'],1,1001);
%ex_diffus_approx=(s_o.^(c./2+1./2)).*1.85.*s_o.^(-1.85).*q2; 
% this value gives the approximate for exfiltration capacity
ex_diffus_approx=0.3477-(0.0731.*c)+(0.0062.*c.^2)-0.0002.*c.^3;
 %ex diffus over 2 meters
%ex_diffus_approx= (1+1.85*(-1./(1.85+1))
 %calclate the indefinate integral for:
% syms s_o c x
% f= x.^((c+1)./2).*(s_o-x).^(0.85)
 %int(f,x,0, s_o)
 
 %d=c+1/2
 
 %result: x.^(c/2 + 1/2).*(s_o - x).^(17/20)
 
 %s_o.^(d).*(1./d+(5./3))+
 
 
 %ex_diffus_approx=(1./d+[(5./3)-1]*(s_o./1-s_o).^1+1./d+[(5./3)-2]*(s_o./1-s_o).^2+1./d+[(5./3)-3]*(s_o./1-s_o).^3+1./d+[(5./3)-4]*(s_o./1-s_o).^4);
 
 %ex_diffus_approx=ex_diffus_approx=(1-s_o)^d(1+1.85.*((-1)^1)./(1.85+1)+((-2)^1)./(1.85)+((-3)^1)./(1.85+3)+((-4)^1)./(1.85+4))
 B= ((1-M)./(1+M.*k_v-w./e_p))+((k_v.*(M.^2)+(1-M).*w./e_p)./((2.*(1+M.*k_v-w./e_p).^2)));
  
C= 1./2.*(M.*k_v-w./e_p).^(-2);
 
E=(alp.*n.*(c-3).*Ksat.*Matrix_pot.*ex_diffus_approx./(pi.*e_p.^2)).*s_o.^(c+5)./2;
  
%beta_s 
 beta=gamma_ratio-(1+(alp.*h_o./e_p)./lambda.*h_o).^(-1.*gamma_depth).*...
     (((gammainc(gamma_depth, (lambda.*h_o + (alp.*h_o)./e_p)))./gamma(gamma_depth)).*exp(-B.*E))+...
     (1 + gamma_ratio).*(1-exp(-B.*E-(alp.*h_o./e_p)).*(1 + M.*k_v + ((2.*B).^(1./2)).*E-(w./e_p))+...
     exp(-C.*E-alp.*h_o./e_p).*((M.*k_v +(2.*C).^(1./2)).*E-(w./e_p))+...
     (((2.*E).^(1./2).*exp(-alp.*h_o./e_p)).*((gammainc(3/2, C.*E)-gammainc(3./2,B.*E)))) + ...
     ((1+((alp.*h_o./e_p)./lambda.*h_o)).^(-1.*gamma_depth)).*...
     (gammainc(gamma_depth, lambda.*h_o+(alp.*h_o./e_p))./gamma(gamma_depth)).*...
     (sqrt(2.*E).*(gammainc(3./2, C.*E)-gammainc(3./2,B.*E))+...
     exp(-C.*E).*(M.*k_v +(((2.*C).^(1./2)).*E)-w./e_p)-...
     exp(-B.*E).*(M.*k_v+(((2.*B).^(1./2)).*E)-w./e_p)));
        
  %compute ET        
ET=((m_n.*e_p)./alp).*(1-M).*beta + M.*k_v ;

%in_diffus_approx=((5./3)+(1./2).*(c+1).*(1-s_o).^(1.425-0.0375.*(c+1)./2)).^(-1);

%used integrate to calculate integration of eqn 11 in Kerkhoff 
in_diffus_approx=(1-s_o).^(5./3).*((3.*(3.*s_o + 5)).*(1 - s_o).^(5./3)./40);
% deleted this part because s_o-s_o=0 ((3*(3*s_o + 5*s_o)*(s_o - s_o)^(5/3))/40))
%%%Runoff model

%mPA mean annual precipitation in mm, not sure if this needs to be changed to cm
  G=(Ksat.*((m_n.*m_r)./m_Pa)).*(((1+s_o.^c)./2)-(w./Ksat));
  sigma=((5.*n.*(c-3).*lambda.^2.*Ksat.*Matrix_pot.*in_diffus_approx.*m_r)./12.*pi.*(gamma_depth.^2)).^(1./3).*(1-s_o).^(2./3);
  Run =m_Pa.*(exp(-G-2.*sigma).*gamma(sigma+1).*(sigma.^(-sigma))+((m_t.*Ksat)./m_Pa).*(s_o.^c));
  scaled = Run./m_Pa;
%compute water balance
    bal=(ET + Run)-m_Pa; 
  %plot contour line
  figure(5)
 [C,h]=contour(M,s_o,bal);
clabel(C,h)
title('Santa Paula, clay')
xlabel('M')
ylabel('s_o')
%%

%for sandyloam
  Ksat=18; %k saturation
      Matrix_pot=200; %soil matrix potential
      c=4; %for sandy loam
      n=0.25; %soil porosity, for sandyy loam 
      
      %%these values set to constant for this purpose
      h_o=0.1; % small constant value for surface water retention
      k_v=1;% unstressed transpiration
      %M=0.5 % ranges from 0 to 1 (but cant be one)
      w=0; %assume no capillary rise

      gamma_ratio = gammainc(gamma_depth, lambda.*h_o)./gamma(gamma_depth);
    
% this value gives the approximate for exfiltration capacity
%ex_diffus_approx=0.3477-(0.0731.*c)+(0.0062.*c.^2)-0.0002.*c.^3;
 %ex diffus over 2 meters
%ex_diffus_approx= (1+1.85*(-1./(1.85+1))
 %calclate the indefinate integral for:
% syms s_o c x
% f= x.^((c+1)./2).*(s_o-x).^(0.85)
 %int(f,x,0, s_o)
 
 %d=c+1/2
 
 %result: x.^(c/2 + 1/2).*(s_o - x).^(17/20)
 
 %s_o.^(d).*(1./d+(5./3))+
 
 % this value gives the approximate for exfiltration capacity
fun = @(x,c,s_o) (x.^(c/2+1/2)).*(s_o-x).^(0.85);


for i=1:1:s_o;
q(i)=integral(@(x)fun(x,c,s_o(i)),0,s_o(i));
end

%need to get q to match s_o
 q2=repmat([q'],1,1001);
ex_diffus_approx=(s_o.^(c./2+1./2)).*1.85.*s_o.^(-1.85).*q2; 
 %ex_diffus_approx=(1./d+[(5./3)-1]*(s_o./1-s_o).^1+1./d+[(5./3)-2]*(s_o./1-s_o).^2+1./d+[(5./3)-3]*(s_o./1-s_o).^3+1./d+[(5./3)-4]*(s_o./1-s_o).^4);
 
 %ex_diffus_approx=ex_diffus_approx=(1-s_o)^d(1+1.85.*((-1)^1)./(1.85+1)+((-2)^1)./(1.85)+((-3)^1)./(1.85+3)+((-4)^1)./(1.85+4))
 B= ((1-M)./(1+M.*k_v-w./e_p))+((k_v.*(M.^2)+(1-M).*w./e_p)./((2.*(1+M.*k_v-w./e_p).^2)));
  
C= 1./2.*(M.*k_v-w./e_p).^(-2);
 
E=(alp.*n.*(c-3).*Ksat.*Matrix_pot.*ex_diffus_approx./(pi.*e_p.^2)).*s_o.^(c+5)./2;
  
%beta_s 
 beta=gamma_ratio-(1+(alp.*h_o./e_p)./lambda.*h_o).^(-1.*gamma_depth).*...
     (((gammainc(gamma_depth, (lambda.*h_o + (alp.*h_o)./e_p)))./gamma(gamma_depth)).*exp(-B.*E))+...
     (1 + gamma_ratio).*(1-exp(-B.*E-(alp.*h_o./e_p)).*(1 + M.*k_v + ((2.*B).^(1./2)).*E-(w./e_p))+...
     exp(-C.*E-alp.*h_o./e_p).*((M.*k_v +(2.*C).^(1./2)).*E-(w./e_p))+...
     (((2.*E).^(1./2).*exp(-alp.*h_o./e_p)).*((gammainc(3/2, C.*E)-gammainc(3./2,B.*E)))) + ...
     ((1+((alp.*h_o./e_p)./lambda.*h_o)).^(-1.*gamma_depth)).*...
     (gammainc(gamma_depth, lambda.*h_o+(alp.*h_o./e_p))./gamma(gamma_depth)).*...
     (sqrt(2.*E).*(gammainc(3./2, C.*E)-gammainc(3./2,B.*E))+...
     exp(-C.*E).*(M.*k_v +(((2.*C).^(1./2)).*E)-w./e_p)-...
     exp(-B.*E).*(M.*k_v+(((2.*B).^(1./2)).*E)-w./e_p)));
        
  %compute ET        
ET=((m_n.*e_p)./alp).*(1-M).*beta + M.*k_v ;
%in_diffus_approx=((5./3)+(1./2).*(c+1).*(1-s_o).^(1.425-0.0375.*(c+1)./2)).^(-1);

%used integrate to calculate integration of eqn 11 in Kerkhoff 
in_diffus_approx=(1-s_o).^(5./3).*((3.*(3.*s_o + 5)).*(1 - s_o).^(5./3)./40);
% deleted this part because s_o-s_o=0 ((3*(3*s_o + 5*s_o)*(s_o - s_o)^(5/3))/40))
%%%Runoff model

%mPA mean annual precipitation in mm, not sure if this needs to be changed to cm
  G=(Ksat.*((m_n.*m_r)./m_Pa)).*(((1+s_o.^c)./2)-(w./Ksat));
  sigma=((5.*n.*(c-3).*lambda.^2.*Ksat.*Matrix_pot.*in_diffus_approx.*m_r)./12.*pi.*(gamma_depth.^2)).^(1./3).*(1-s_o).^(2./3);
  Run =m_Pa.*(exp(-G-2.*sigma).*gamma(sigma+1).*(sigma.^(-sigma))+((m_t.*Ksat)./m_Pa).*(s_o.^c));
  scaled = Run./m_Pa;
%compute water balance
    bal=(ET + Run)-m_Pa; 
  %plot contour line
  figure(6)
 [C,h]=contour(M,s_o,bal);
clabel(C,h)
title('Santa Paula, Sandyloam soil')
xlabel('M')
ylabel('s_o')


%%
% for silt
Ksat=8.64; %k saturation
        Matrix_pot=166; %soil matrix potential
        c=6; %for clayloam
        n=0.35; %soil porosity, for clay
        
        %%these values set to constant for this purpose
        h_o=0.1; % small constant value for surface water retention
        k_v=1;% unstressed transpiration
        %M=0.5 % ranges from 0 to 1 (but cant be one)
        w=0; %assume no capillary rise
d=2.9
gamma_ratio = gammainc(gamma_depth, lambda.*h_o)./gamma(gamma_depth);
    
% this value gives the approximate for exfiltration capacity
%ex_diffus_approx=0.3477-(0.0731.*c)+(0.0062.*c.^2)-0.0002.*c.^3;
 %ex diffus over 2 meters
%ex_diffus_approx= (1+1.85*(-1./(1.85+1))
 %calclate the indefinate integral for:
% syms s_o c x
% f= x.^((c+1)./2).*(s_o-x).^(0.85)
 %int(f,x,0, s_o)
 
 %d=c+1/2
 
 %result: x.^(c/2 + 1/2).*(s_o - x).^(17/20)
 
 %s_o.^(d).*(1./d+(5./3))+
 % this value gives the approximate for exfiltration capacity
fun = @(x,c,s_o) (x.^(c/2+1/2)).*(s_o-x).^(0.85);


for i=1:1:s_o;
q(i)=integral(@(x)fun(x,c,s_o(i)),0,s_o(i));
end

%need to get q to match s_o
 q2=repmat([q'],1,1001);
ex_diffus_approx=(s_o.^(c./2+1./2)).*1.85.*s_o.^(-1.85).*q2; 
 
 %ex_diffus_approx=(1./d+[(5./3)-1]*(s_o./1-s_o).^1+1./d+[(5./3)-2]*(s_o./1-s_o).^2+1./d+[(5./3)-3]*(s_o./1-s_o).^3+1./d+[(5./3)-4]*(s_o./1-s_o).^4);
 
 %ex_diffus_approx=ex_diffus_approx=(1-s_o)^d(1+1.85.*((-1)^1)./(1.85+1)+((-2)^1)./(1.85)+((-3)^1)./(1.85+3)+((-4)^1)./(1.85+4))
 B= ((1-M)./(1+M.*k_v-w./e_p))+((k_v.*(M.^2)+(1-M).*w./e_p)./((2.*(1+M.*k_v-w./e_p).^2)));
  
C= 1./2.*(M.*k_v-w./e_p).^(-2);
 
E=(alp.*n.*(c-3).*Ksat.*Matrix_pot.*ex_diffus_approx./(pi.*e_p.^2)).*s_o.^(c+5)./2;
  
%beta_s 
 beta=gamma_ratio-(1+(alp.*h_o./e_p)./lambda.*h_o).^(-1.*gamma_depth).*...
     (((gammainc(gamma_depth, (lambda.*h_o + (alp.*h_o)./e_p)))./gamma(gamma_depth)).*exp(-B.*E))+...
     (1 + gamma_ratio).*(1-exp(-B.*E-(alp.*h_o./e_p)).*(1 + M.*k_v + ((2.*B).^(1./2)).*E-(w./e_p))+...
     exp(-C.*E-alp.*h_o./e_p).*((M.*k_v +(2.*C).^(1./2)).*E-(w./e_p))+...
     (((2.*E).^(1./2).*exp(-alp.*h_o./e_p)).*((gammainc(3/2, C.*E)-gammainc(3./2,B.*E)))) + ...
     ((1+((alp.*h_o./e_p)./lambda.*h_o)).^(-1.*gamma_depth)).*...
     (gammainc(gamma_depth, lambda.*h_o+(alp.*h_o./e_p))./gamma(gamma_depth)).*...
     (sqrt(2.*E).*(gammainc(3./2, C.*E)-gammainc(3./2,B.*E))+...
     exp(-C.*E).*(M.*k_v +(((2.*C).^(1./2)).*E)-w./e_p)-...
     exp(-B.*E).*(M.*k_v+(((2.*B).^(1./2)).*E)-w./e_p)));
        
  %compute ET        
ET=((m_n.*e_p)./alp).*(1-M).*beta + M.*k_v ;

%in_diffus_approx=((5./3)+(1./2).*(c+1).*(1-s_o).^(1.425-0.0375.*(c+1)./2)).^(-1);

%used integrate to calculate integration of eqn 11 in Kerkhoff 
in_diffus_approx=(1-s_o).^(5./3).*((3.*(3.*s_o + 5)).*(1 - s_o).^(5./3)./40);
% deleted this part because s_o-s_o=0 ((3*(3*s_o + 5*s_o)*(s_o - s_o)^(5/3))/40))
%%%Runoff model

%mPA mean annual precipitation in mm, not sure if this needs to be changed to cm
  G=(Ksat.*((m_n.*m_r)./m_Pa)).*(((1+s_o.^c)./2)-(w./Ksat));
  sigma=((5.*n.*(c-3).*lambda.^2.*Ksat.*Matrix_pot.*in_diffus_approx.*m_r)./12.*pi.*(gamma_depth.^2)).^(1./3).*(1-s_o).^(2./3);
  Run =m_Pa.*(exp(-G-2.*sigma).*gamma(sigma+1).*(sigma.^(-sigma))+((m_t.*Ksat)./m_Pa).*(s_o.^c));
  scaled = Run./m_Pa;
  
%compute water balance
    bal=(ET + Run)-m_Pa; 
  %plot contour line
  figure(7)
 [C,h]=contour(M,s_o,bal);
clabel(C,h)
title('Santa Paula, silt loam')
xlabel('M')
ylabel('s_o')
%%

%for clayloam
%%

Ksat=2.0; %k saturation
        Matrix_pot=19; %soil matrix potential
        c=10; %for clayloam
        n=0.35; %soil porosity, for clay
        d=4.3;
        %%these values set to constant for this purpose
        h_o=0.1; % small constant value for surface water retention
        k_v=1;% unstressed transpiration
        %M=0.5 % ranges from 0 to 1 (but cant be one)
        w=0; %assume no capillary rise

gamma_ratio = gammainc(gamma_depth, lambda.*h_o)./gamma(gamma_depth);
    
% this value gives the approximate for exfiltration capacity
% this value gives the approximate for exfiltration capacity
fun = @(x,c,s_o) (x.^(c/2+1/2)).*(s_o-x).^(0.85);


%for i=1:1:s_o;
%q(i)=integral(@(x)fun(x,c,s_o(i)),0,s_o(i));
%end

%need to get q to match s_o
 %q2=repmat([q'],1,1001);
%ex_diffus_approx=(s_o.^(c./2+1./2)).*1.85.*s_o.^(-1.85).*q2; 
ex_diffus_approx=0.3477-(0.0731.*c)+(0.0062.*c.^2)-0.0002.*c.^3;
 %ex diffus over 2 meters
%ex_diffus_approx= (1+1.85*(-1./(1.85+1))
 %calclate the indefinate integral for:
% syms s_o c x
% f= x.^((c+1)./2).*(s_o-x).^(0.85)
 %int(f,x,0, s_o)
 
 %d=c+1/2
 
 %result: x.^(c/2 + 1/2).*(s_o - x).^(17/20)
 
 %s_o.^(d).*(1./d+(5./3))+
 
 
 %ex_diffus_approx=(1./d+[(5./3)-1]*(s_o./1-s_o).^1+1./d+[(5./3)-2]*(s_o./1-s_o).^2+1./d+[(5./3)-3]*(s_o./1-s_o).^3+1./d+[(5./3)-4]*(s_o./1-s_o).^4);
 
 %ex_diffus_approx=ex_diffus_approx=(1-s_o)^d(1+1.85.*((-1)^1)./(1.85+1)+((-2)^1)./(1.85)+((-3)^1)./(1.85+3)+((-4)^1)./(1.85+4))
 B= ((1-M)./(1+M.*k_v-w./e_p))+((k_v.*(M.^2)+(1-M).*w./e_p)./((2.*(1+M.*k_v-w./e_p).^2)));
  
C= 1./2.*(M.*k_v-w./e_p).^(-2);
 
E=(alp.*n.*(c-3).*Ksat.*Matrix_pot.*ex_diffus_approx./(pi.*e_p.^2)).*s_o.^(c+5)./2;
  
%beta_s 
 beta=gamma_ratio-(1+(alp.*h_o./e_p)./lambda.*h_o).^(-1.*gamma_depth).*...
     (((gammainc(gamma_depth, (lambda.*h_o + (alp.*h_o)./e_p)))./gamma(gamma_depth)).*exp(-B.*E))+...
     (1 + gamma_ratio).*(1-exp(-B.*E-(alp.*h_o./e_p)).*(1 + M.*k_v + ((2.*B).^(1./2)).*E-(w./e_p))+...
     exp(-C.*E-alp.*h_o./e_p).*((M.*k_v +(2.*C).^(1./2)).*E-(w./e_p))+...
     (((2.*E).^(1./2).*exp(-alp.*h_o./e_p)).*((gammainc(3/2, C.*E)-gammainc(3./2,B.*E)))) + ...
     ((1+((alp.*h_o./e_p)./lambda.*h_o)).^(-1.*gamma_depth)).*...
     (gammainc(gamma_depth, lambda.*h_o+(alp.*h_o./e_p))./gamma(gamma_depth)).*...
     (sqrt(2.*E).*(gammainc(3./2, C.*E)-gammainc(3./2,B.*E))+...
     exp(-C.*E).*(M.*k_v +(((2.*C).^(1./2)).*E)-w./e_p)-...
     exp(-B.*E).*(M.*k_v+(((2.*B).^(1./2)).*E)-w./e_p)));
        
  %compute ET        
ET=((m_n.*e_p)./alp).*(1-M).*beta + M.*k_v ;

%in_diffus_approx=((5./3)+(1./2).*(c+1).*(1-s_o).^(1.425-0.0375.*(c+1)./2)).^(-1);

%used integrate to calculate integration of eqn 11 in Kerkhoff 
in_diffus_approx=(1-s_o).^(5./3).*((3.*(3.*s_o + 5)).*(1 - s_o).^(5./3)./40);
% deleted this part because s_o-s_o=0 ((3*(3*s_o + 5*s_o)*(s_o - s_o)^(5/3))/40))
%%%Runoff model

%mPA mean annual precipitation in mm, not sure if this needs to be changed to cm
  G=(Ksat.*((m_n.*m_r)./m_Pa)).*(((1+s_o.^c)./2)-(w./Ksat));
  sigma=((5.*n.*(c-3).*lambda.^2.*Ksat.*Matrix_pot.*in_diffus_approx.*m_r)./12.*pi.*(gamma_depth.^2)).^(1./3).*(1-s_o).^(2./3);
  Run =m_Pa.*(exp(-G-2.*sigma).*gamma(sigma+1).*(sigma.^(-sigma))+((m_t.*Ksat)./m_Pa).*(s_o.^c));
  scaled = Run./m_Pa;
  
%compute water balance
    bal=(ET + Run)-m_Pa; 
  %plot contour line
  figure(8)
 [C,h]=contour(M,s_o,bal);
clabel(C,h)
title('Santa Paula, clay loam')
xlabel('M')
ylabel('s_o')