Pr.a<-ET.model.a+Runoff.a
ET.model.a<-((m.n*e.p)/alpha)*(1-M)*beta.s +M*k.v
#m.n=mean number rainstorms per year
#Beta.s=avg bare soil evap efficiency
#M=vegetation canopy denisity
#ET.model.a=evaporation component of water balance

#treat M as a state variable

#parameters
m.n<-15 # mean number of storms per year
e.p<-0.15
m.Pa<-111.3
m.r<-0.32
m.t<-365
alpha<-5 # made this up
lambda<-0.5
gamma.depth<-0.5

Ksat<-0.72
Matrix.pot<-25
c<-12
n<-0.45
h.o<-0.1
k.v<-1
M<-0.5 # ranges from 0 to 1 (but cant be one)
w<-0

###Beta.s model
library(gsl)
gamma.ratio<-gamma_inc(gamma.depth, lambda*h.o)/gamma(gamma.depth)
B<- ((1-M)/(1+M*k.v-w/e.p))+(k.v*(M^2)+(1-M)*w/e.p)/((2*(1+M*k.v-w/e.p)))

C<- 1/2*(M*k.v-w/e.p)^(-2)

E<-((alpha*n*(c-3)*Ksat*Matrix.pot*ex.diffus)/(pi*e.p^2))*s.o^((C+5)/2)

beta.s<-(gamma.ratio)-((1+(alpha*h.o/e.p)/lambda*h.o)^(-1*gamma.depth))*((gamma_inc(gamma.depth, lambda*h.o +alpha*lambda*h.o/e.p))/gamma(gamma.depth))*exp(-B*E)+(1+gamma.ratio)*
  (1-exp(-B*E-alpha*h.o/e.p)*(1+M*k.v+(2*B)^(1/2)*E-w/e.p)+
     exp(-C*E-alpha*h.o/e.p)*(M*k.v +(2*C)^(1/2)*E-w/e.p)+
     (2*E)^(1/2)*exp(-alpha*h.o/e.p)*(gamma_inc(3/2, C*E)-gamma_inc(3/2,B*E) + 
     ((1+(alpha*h.o/e.p)/lambda*h.o)^(-1*gamma.depth))*(gamma_inc(gamma.depth, lambda*h.o+alpha*h.o/e.p )/gamma(gamma.depth))*(sqrt(2*E)*(gamma_inc(3/2, C*E)-gamma_inc(3/2,B*E))+
           exp(-C*E)*(M*k.v +(2*C)^(1/2)*E-w/e.p) -
             exp(-B*E)*(M*k.v+(2*B)^(1/2)*E-w/e.p)



#need to define exfiltration diffusivity
s.o<- 0.8 #made this up , but should be between 0 and 1 
ex.diffus<-s.o^((C+1)/2)*(1.85*s.o^(-1.85)*s.o^(((C+1)/2)+1)*(s.o-s.o^2)^0.85)

#in.diffus and out.diffus are incorrectly specified

###Runoff model
mPa<-1100 # mean annual precipitation
G<-(Ksat*(m.n*m.r)/mPa)*((1+s.o^C)/2)-w/Ksat
in.diffus<-(1-s.o)^(5/3)*((1^(((C+1)/2)+1)*(1^2-s.o)^(2/3) - (s.o^(((C+1)/2)+1)*((s.o^2-s.o)^(2/3))))
sigma<-((5*n*(c-3)*lambda^2*Ksat*Matrix.pot*in.diffus*m.r)/12*pi*gamma.depth^2)^(1/3)*(1-s.o)^(2/3)

Runoff.a<-mPa*(exp(-G-2*sigma)*gamma(sigma+1)*sigma^(-sigma)+((m.t*Ksat)/mPa)*s.o^C)
