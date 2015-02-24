%% define meshgrid
[M, s_o]=meshgrid(0:0.001:1,0:0.001:1);

%% for los pinos minimum
m_h=3.92 %mean storm depth in mm
 m_n=10; % mean number of storms per year % obtained from Eagles_on: Climate, s_oil, and Vegetation
  e_p=0.32; % average bare s_oil potential evaporation rate
  m_Pa=27.5; %mean annual precipitation in cm
  m_r=0.09; %mean storm depth--double check value
  m_b=3.17; %in units of days
  m_t=365; %mean rainy seas_on
  alp=1./m_b; %one over mean time between storms
  %Eagles_on (2002) has poiss_on pulse parameters for a number of stations in Appendix F
  %for Clinton, MA, I used parameters for Boston (Station ID 14 in eagles_on
  gamma_depth=0.50; %these may be incorrect %shape parameter for the gamma distribution of storm depths
  lambda=gamma_depth/m_h; %may be incorrect, changed inches to cm %scale parameter for gamma distribution of storm depths (1./cm)
  
  
%% for soil 1
Ksat=0.0028; %k saturation
  Matrix_pot=343.31; %s_oil matrix potential
  c =13.83; %for clay
  n=0.52; %s_oil porosity, for clay
 %k_one=0.000000000141
  h_o=0.1; % small constant value for surface water retention
  k_v=1;% unstressed transpiration
  w=0; %assume no capillary rise
  %%
  %% 
%ET model
% this value gives the approximate for exfiltration capacity
%ex_diffus_approx=0.3477-(0.0731.*c)+(0.0062.*c.^2)-0.0002.*c.^3;
 
 %calclate the indefinate integral for:
 fun = @(x,c,s_o) (x.^(c/2+1/2)).*(s_o-x).^(0.85);
%s_o_col=s_o(1:1001,1);
q= zeros(1001);
for i=1:1001;
q(i)=integral(@(x)fun(x,c,s_o(i)),0,s_o(i));
end

%need to get q to match s_o
 q2=repmat(q(:,1),1,1001);
ex_diffus_approx=(s_o.^(c./2+1./2)).*1.85.*(s_o.^(-1.85)).*q2; 
 
 
 
B= ((1-M)./(1+M.*k_v-w./e_p))+((k_v.*(M.^2)+(1-M).*w./e_p)./((2.*(1+M.*k_v-w./e_p).^2)));
  
C= 1./2.*(M.*k_v-w./e_p).^(-2);
 
E=(alp.*n.*(c-3).*Ksat.*Matrix_pot.*ex_diffus_approx./(pi.*e_p.^2)).*s_o.^(c+5)./2;

gamma_ratio = gammainc(gamma_depth, lambda.*h_o)./gamma(gamma_depth);
 gamma_bar_evap=gammainc(gamma_depth, (lambda.*h_o + alp.*h_o./e_p))./gamma(gamma_depth);
 time_evap=(1+ ((alp.*h_o./e_p)./(lambda.*h_o))).^gamma_depth;
 gamma_diff=(gammainc(3/2, C.*E)-gammainc(3/2,B.*E));
 C_mkv=M.*k_v+((2.*C).^(1/2)).*E-(w./e_p);
 B_mkv=M.*k_v+((2.*B).^(1/2)).*E-(w./e_p);


beta= gamma_ratio - time_evap.*gamma_bar_evap.*exp(-B.*E)+ ((1 - gamma_ratio).*((1-exp(-B.*E-(alp.*h_o./e_p)).*(1+B_mkv)+...
     exp(-C.*E-(alp.*h_o./e_p)).*C_mkv)+((2.*E).^(0.5)).*exp(-alp.*h_o./e_p).*gamma_diff))+...
     time_evap.*gamma_bar_evap.*(((2.*E).^(0.5).*gamma_diff) + exp(-C.*E).*C_mkv - exp(-B.*E).*B_mkv);

%compute ET        
ET=((m_n.*e_p)./alp).*(1-M).*beta + M.*k_v ;

%in_diffus_approx=((5./3)+(1./2).*(c+1).*(1-s_o).^(1.425-0.0375.*(c+1)./2)).^(-1);
 fun = @(x,c,s_o) (x.^(c/2+1/2)).*(x-s_o).^(2/3);
%s_o_col=s_o(1001);
q3= zeros(1001);
for i=1:1001;
q4(i)=integral(@(x)fun(x,c,s_o(i)),s_o(i),1);
end
q5=repmat(q4(1,:),1001,1);

%used integrate to calculate integration of eqn 11 in Kerkhoff 
%create this into actual integral as well
in_diffus_approx=((1-s_o).^(5./3)).*q5;

%%%Runoff model

%mPA mean annual precipitation in mm, not sure if this needs to be changed to cm
  G=(Ksat.*((m_n.*m_r)./m_Pa)).*(((1+s_o.^c)./2)-(w./Ksat));
  
  sigma=((5.*n.*(c-3).*lambda.^2.*Ksat.*Matrix_pot.*in_diffus_approx.*m_r)./12.*pi.*(gamma_depth.^2)).^(1./3).*(1-s_o).^(2./3);
  
  Run =m_Pa.*(exp(-G-2.*sigma).*gamma(sigma+1).*(sigma.^(-sigma))+((m_t.*Ksat)./m_Pa).*(s_o.^c));
 
%compute water balance
    bal=(ET + Run)-m_Pa; 
 %plot contour line
  figure(1)
 [C,h]=contour(M,s_o,bal);
clabel(C,h)
title('Los Pinos, soil1')
xlabel('M')
ylabel('s_o')
%%
%%for soil 2
Ksat=0.65; %k saturation
  Matrix_pot=52.53; %s_oil matrix potential
  c =17.07; %for clay
  n=0.47; %s_oil porosity, for clay
 %k_one=0.000000000141
  h_o=0.1; % small constant value for surface water retention
  k_v=1;% unstressed transpiration
  w=0; %assume no capillary rise
%% 
%% 
%ET model
% this value gives the approximate for exfiltration capacity
%ex_diffus_approx=0.3477-(0.0731.*c)+(0.0062.*c.^2)-0.0002.*c.^3;
 
 %calclate the indefinate integral for:
 fun = @(x,c,s_o) (x.^(c/2+1/2)).*(s_o-x).^(0.85);
%s_o_col=s_o(1:1001,1);
q= zeros(1001);
for i=1:1001;
q(i)=integral(@(x)fun(x,c,s_o(i)),0,s_o(i));
end

%need to get q to match s_o
 q2=repmat(q(:,1),1,1001);
ex_diffus_approx=(s_o.^(c./2+1./2)).*1.85.*(s_o.^(-1.85)).*q2; 
 
 
 
B= ((1-M)./(1+M.*k_v-w./e_p))+((k_v.*(M.^2)+(1-M).*w./e_p)./((2.*(1+M.*k_v-w./e_p).^2)));
  
C= 1./2.*(M.*k_v-w./e_p).^(-2);
 
E=(alp.*n.*(c-3).*Ksat.*Matrix_pot.*ex_diffus_approx./(pi.*e_p.^2)).*s_o.^(c+5)./2;

gamma_ratio = gammainc(gamma_depth, lambda.*h_o)./gamma(gamma_depth);
 gamma_bar_evap=gammainc(gamma_depth, (lambda.*h_o + alp.*h_o./e_p))./gamma(gamma_depth);
 time_evap=(1+ ((alp.*h_o./e_p)./(lambda.*h_o))).^gamma_depth;
 gamma_diff=(gammainc(3/2, C.*E)-gammainc(3/2,B.*E));
 C_mkv=M.*k_v+((2.*C).^(1/2)).*E-(w./e_p);
 B_mkv=M.*k_v+((2.*B).^(1/2)).*E-(w./e_p);


beta= gamma_ratio - time_evap.*gamma_bar_evap.*exp(-B.*E)+ ((1 - gamma_ratio).*((1-exp(-B.*E-(alp.*h_o./e_p)).*(1+B_mkv)+...
     exp(-C.*E-(alp.*h_o./e_p)).*C_mkv)+((2.*E).^(0.5)).*exp(-alp.*h_o./e_p).*gamma_diff))+...
     time_evap.*gamma_bar_evap.*(((2.*E).^(0.5).*gamma_diff) + exp(-C.*E).*C_mkv - exp(-B.*E).*B_mkv);

%compute ET        
ET=((m_n.*e_p)./alp).*(1-M).*beta + M.*k_v ;

%in_diffus_approx=((5./3)+(1./2).*(c+1).*(1-s_o).^(1.425-0.0375.*(c+1)./2)).^(-1);
 fun = @(x,c,s_o) (x.^(c/2+1/2)).*(x-s_o).^(2/3);
%s_o_col=s_o(1001);
q3= zeros(1001);
for i=1:1001;
q4(i)=integral(@(x)fun(x,c,s_o(i)),s_o(i),1);
end
q5=repmat(q4(1,:),1001,1);

%used integrate to calculate integration of eqn 11 in Kerkhoff 
%create this into actual integral as well
in_diffus_approx=((1-s_o).^(5./3)).*q5;

%%%Runoff model

%mPA mean annual precipitation in mm, not sure if this needs to be changed to cm
  G=(Ksat.*((m_n.*m_r)./m_Pa)).*(((1+s_o.^c)./2)-(w./Ksat));
  
  sigma=((5.*n.*(c-3).*lambda.^2.*Ksat.*Matrix_pot.*in_diffus_approx.*m_r)./12.*pi.*(gamma_depth.^2)).^(1./3).*(1-s_o).^(2./3);
  
  Run =m_Pa.*(exp(-G-2.*sigma).*gamma(sigma+1).*(sigma.^(-sigma))+((m_t.*Ksat)./m_Pa).*(s_o.^c));
 
%compute water balance
    bal=(ET + Run)-m_Pa; 
 %plot contour line
  figure(2)
 [C,h]=contour(M,s_o,bal);
clabel(C,h)
title('Los Pinos, soil2')
xlabel('M')
ylabel('s_o')
%%
%%for soil 3
Ksat=0.10; %k saturation
  Matrix_pot=93.33; %s_oil matrix potential
  c =9.32; %for clay
  n=0.43; %s_oil porosity, for clay
 %k_one=0.000000000141
  h_o=0.1; % small constant value for surface water retention
  k_v=1;% unstressed transpiration
  w=0; %assume no capillary rise
%%
%% 
%ET model
% this value gives the approximate for exfiltration capacity
%ex_diffus_approx=0.3477-(0.0731.*c)+(0.0062.*c.^2)-0.0002.*c.^3;
 
 %calclate the indefinate integral for:
 fun = @(x,c,s_o) (x.^(c/2+1/2)).*(s_o-x).^(0.85);
%s_o_col=s_o(1:1001,1);
q= zeros(1001);
for i=1:1001;
q(i)=integral(@(x)fun(x,c,s_o(i)),0,s_o(i));
end

%need to get q to match s_o
 q2=repmat(q(:,1),1,1001);
ex_diffus_approx=(s_o.^(c./2+1./2)).*1.85.*(s_o.^(-1.85)).*q2; 
 
 
 
B= ((1-M)./(1+M.*k_v-w./e_p))+((k_v.*(M.^2)+(1-M).*w./e_p)./((2.*(1+M.*k_v-w./e_p).^2)));
  
C= 1./2.*(M.*k_v-w./e_p).^(-2);
 
E=(alp.*n.*(c-3).*Ksat.*Matrix_pot.*ex_diffus_approx./(pi.*e_p.^2)).*s_o.^(c+5)./2;

gamma_ratio = gammainc(gamma_depth, lambda.*h_o)./gamma(gamma_depth);
 gamma_bar_evap=gammainc(gamma_depth, (lambda.*h_o + alp.*h_o./e_p))./gamma(gamma_depth);
 time_evap=(1+ ((alp.*h_o./e_p)./(lambda.*h_o))).^gamma_depth;
 gamma_diff=(gammainc(3/2, C.*E)-gammainc(3/2,B.*E));
 C_mkv=M.*k_v+((2.*C).^(1/2)).*E-(w./e_p);
 B_mkv=M.*k_v+((2.*B).^(1/2)).*E-(w./e_p);


beta= gamma_ratio - time_evap.*gamma_bar_evap.*exp(-B.*E)+ ((1 - gamma_ratio).*((1-exp(-B.*E-(alp.*h_o./e_p)).*(1+B_mkv)+...
     exp(-C.*E-(alp.*h_o./e_p)).*C_mkv)+((2.*E).^(0.5)).*exp(-alp.*h_o./e_p).*gamma_diff))+...
     time_evap.*gamma_bar_evap.*(((2.*E).^(0.5).*gamma_diff) + exp(-C.*E).*C_mkv - exp(-B.*E).*B_mkv);

%compute ET        
ET=((m_n.*e_p)./alp).*(1-M).*beta + M.*k_v ;

%in_diffus_approx=((5./3)+(1./2).*(c+1).*(1-s_o).^(1.425-0.0375.*(c+1)./2)).^(-1);
 fun = @(x,c,s_o) (x.^(c/2+1/2)).*(x-s_o).^(2/3);
%s_o_col=s_o(1001);
q3= zeros(1001);
for i=1:1001;
q4(i)=integral(@(x)fun(x,c,s_o(i)),s_o(i),1);
end
q5=repmat(q4(1,:),1001,1);

%used integrate to calculate integration of eqn 11 in Kerkhoff 
%create this into actual integral as well
in_diffus_approx=((1-s_o).^(5./3)).*q5;

%%%Runoff model

%mPA mean annual precipitation in mm, not sure if this needs to be changed to cm
  G=(Ksat.*((m_n.*m_r)./m_Pa)).*(((1+s_o.^c)./2)-(w./Ksat));
  
  sigma=((5.*n.*(c-3).*lambda.^2.*Ksat.*Matrix_pot.*in_diffus_approx.*m_r)./12.*pi.*(gamma_depth.^2)).^(1./3).*(1-s_o).^(2./3);
  
  Run =m_Pa.*(exp(-G-2.*sigma).*gamma(sigma+1).*(sigma.^(-sigma))+((m_t.*Ksat)./m_Pa).*(s_o.^c));
 
%compute water balance
    bal=(ET + Run)-m_Pa; 
 %plot contour line
  figure(3)
 [C,h]=contour(M,s_o,bal);
clabel(C,h)
title('Los Pinos, soil3')
xlabel('M')
ylabel('s_o')

%% for soil 4 

Ksat=0.083; %k saturation
  Matrix_pot=75.08; %s_oil matrix potential
  c =6.44; %for clay
  n=0.42; %s_oil porosity, for clay
 %k_one=0.000000000141
  h_o=0.1; % small constant value for surface water retention
  k_v=1;% unstressed transpiration
  w=0; %assume no capillary rise
 %%
%ET model
% this value gives the approximate for exfiltration capacity
%ex_diffus_approx=0.3477-(0.0731.*c)+(0.0062.*c.^2)-0.0002.*c.^3;
 
 %calclate the indefinate integral for:
 fun = @(x,c,s_o) (x.^(c/2+1/2)).*(s_o-x).^(0.85);
%s_o_col=s_o(1:1001,1);
q= zeros(1001);
for i=1:1001;
q(i)=integral(@(x)fun(x,c,s_o(i)),0,s_o(i));
end

%need to get q to match s_o
 q2=repmat(q(:,1),1,1001);
ex_diffus_approx=(s_o.^(c./2+1./2)).*1.85.*(s_o.^(-1.85)).*q2; 
 
 
 
B= ((1-M)./(1+M.*k_v-w./e_p))+((k_v.*(M.^2)+(1-M).*w./e_p)./((2.*(1+M.*k_v-w./e_p).^2)));
  
C= 1./2.*(M.*k_v-w./e_p).^(-2);
 
E=(alp.*n.*(c-3).*Ksat.*Matrix_pot.*ex_diffus_approx./(pi.*e_p.^2)).*s_o.^(c+5)./2;

gamma_ratio = gammainc(gamma_depth, lambda.*h_o)./gamma(gamma_depth);
 gamma_bar_evap=gammainc(gamma_depth, (lambda.*h_o + alp.*h_o./e_p))./gamma(gamma_depth);
 time_evap=(1+ ((alp.*h_o./e_p)./(lambda.*h_o))).^gamma_depth;
 gamma_diff=(gammainc(3/2, C.*E)-gammainc(3/2,B.*E));
 C_mkv=M.*k_v+((2.*C).^(1/2)).*E-(w./e_p);
 B_mkv=M.*k_v+((2.*B).^(1/2)).*E-(w./e_p);


beta= gamma_ratio - time_evap.*gamma_bar_evap.*exp(-B.*E)+ ((1 - gamma_ratio).*((1-exp(-B.*E-(alp.*h_o./e_p)).*(1+B_mkv)+...
     exp(-C.*E-(alp.*h_o./e_p)).*C_mkv)+((2.*E).^(0.5)).*exp(-alp.*h_o./e_p).*gamma_diff))+...
     time_evap.*gamma_bar_evap.*(((2.*E).^(0.5).*gamma_diff) + exp(-C.*E).*C_mkv - exp(-B.*E).*B_mkv);

%compute ET        
ET=((m_n.*e_p)./alp).*(1-M).*beta + M.*k_v ;

%in_diffus_approx=((5./3)+(1./2).*(c+1).*(1-s_o).^(1.425-0.0375.*(c+1)./2)).^(-1);
 fun = @(x,c,s_o) (x.^(c/2+1/2)).*(x-s_o).^(2/3);
%s_o_col=s_o(1001);
q3= zeros(1001);
for i=1:1001;
q4(i)=integral(@(x)fun(x,c,s_o(i)),s_o(i),1);
end
q5=repmat(q4(1,:),1001,1);

%used integrate to calculate integration of eqn 11 in Kerkhoff 
%create this into actual integral as well
in_diffus_approx=((1-s_o).^(5./3)).*q5;

%%%Runoff model

%mPA mean annual precipitation in mm, not sure if this needs to be changed to cm
  G=(Ksat.*((m_n.*m_r)./m_Pa)).*(((1+s_o.^c)./2)-(w./Ksat));
  
  sigma=((5.*n.*(c-3).*lambda.^2.*Ksat.*Matrix_pot.*in_diffus_approx.*m_r)./12.*pi.*(gamma_depth.^2)).^(1./3).*(1-s_o).^(2./3);
  
  Run =m_Pa.*(exp(-G-2.*sigma).*gamma(sigma+1).*(sigma.^(-sigma))+((m_t.*Ksat)./m_Pa).*(s_o.^c));
 
%compute water balance
    bal=(ET + Run)-m_Pa; 
 %plot contour line
  figure(4)
 [C,h]=contour(M,s_o,bal);
clabel(C,h)
title('los pinos, soil4')
xlabel('M')
ylabel('s_o')
%%
