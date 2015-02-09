Pr.a<-ET.model.a+Runoff.a
ET.exp<-mPa-Runoff.a

ET.model.a<-matrix(0,10,1)
for(i in 1:length(M)){
    ET.model.a[i]<-((m.n*e.p)/alpha)*(1-M[i])*beta.s(M[i],k.v) + M[i]*k.v
  
}

Runoff.model<-matrix(0,1, length(s.o))
for(i in 1:length(s.o)){
Runoff.model[i]<-Runoff(s.o[i])
}

#rownames(ET.model.a)<-M

dif<-ET.model.a - ET.exp
#ET.model.a<-((m.n*e.p)/alpha)*(1-M)*beta.s(M,k.v) + M*k.v

#m.n=mean number rainstorms per year
#Beta.s=avg bare soil evap efficiency
#M=vegetation canopy denisity
#ET.model.a=evaporation component of water balance

#treat M and k.v as state variables
M<-c(0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, .9999)
k.v<-0.097/e.p #unstressed composite transpiration rate/e.p
#parameters
m.n<-15 # mean number of storms per year
e.p<-0.15 # average bare soil potential evaporation rate
m.Pa<-111.3
m.r<-0.32
m.t<-365
alpha<-1/m.t 
lambda<-0.5
gamma.depth<-0.5

Ksat<-0.72
Matrix.pot<-25
c<-12 #for clay
n<-0.45 #soil porosity, for clay
h.o<-0.1 # small constant value for surface water retention
k.v<-1# unstressed transpiration
M<-0.5 # ranges from 0 to 1 (but cant be one)
w<-0 #assume no capillary rise

###Beta.s model

require(gsl)
gamma.ratio<-gamma_inc(gamma.depth, lambda*h.o)/gamma(gamma.depth)

#define beta.s function
beta.s<-function(M, k.v){
  B<- ((1-M)/(1+M*k.v-w/e.p))+(k.v*(M^2)+(1-M)*w/e.p)/((2*(1+M*k.v-w/e.p)))
  
  C<- 1/2*(M*k.v-w/e.p)^(-2)
  
  E<-((alpha*n*(c-3)*Ksat*Matrix.pot*ex.diffus)/(pi*e.p^2))*s.o^((C+5)/2)
  
  beta<-(gamma.ratio)-((1+(alpha*h.o/e.p)/lambda*h.o)^(-1*gamma.depth))*((gamma_inc(gamma.depth, lambda*h.o +alpha*lambda*h.o/e.p))/gamma(gamma.depth))*exp(-B*E)+(1+gamma.ratio)*
  (1-exp(-B*E-alpha*h.o/e.p)*(1+M*k.v+(2*B)^(1/2)*E-w/e.p)+
     exp(-C*E-alpha*h.o/e.p)*(M*k.v +(2*C)^(1/2)*E-w/e.p)+
     (2*E)^(1/2)*exp(-alpha*h.o/e.p)*(gamma_inc(3/2, C*E)-gamma_inc(3/2,B*E) + 
     ((1+(alpha*h.o/e.p)/lambda*h.o)^(-1*gamma.depth))*(gamma_inc(gamma.depth, lambda*h.o+alpha*h.o/e.p )/gamma(gamma.depth))*(sqrt(2*E)*(gamma_inc(3/2, C*E)-gamma_inc(3/2,B*E))+
           exp(-C*E)*(M*k.v +(2*C)^(1/2)*E-w/e.p) -
             exp(-B*E)*(M*k.v+(2*B)^(1/2)*E-w/e.p))))

}




#define exfiltration diffusivity
s.o<-c(0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, .9999) #made this up , but should be between 0 and 1,and need to solve for s.o that satisfies water balance closure 


Runoff<-function(s.o){
## define the integrated function
  integrand.ex <- function(x) {x^((c+1)/2)*(s.o-x)^(0.85)}
## integrate the function from 0 to s.o
  b<-integrate(integrand.ex, lower = 0, upper = s.o)
#b$value provides the value of the intergral
  ex.diffus<-s.o^((c+1)/2)*(1.85*s.o^(-1.85)*b$value) #in.diffus and out.diffus are incorrectly specified

###Runoff model
  mPa<-1100 # mean annual precipitation
  G<-(Ksat*(m.n*m.r)/mPa)*((1+s.o^c)/2)-w/Ksat

  integrand.in<-function(x){x^((c+1)/2)*(x-s.o)^(2/3)}
  b.in<-integrate(integrand.in, lower = s.o, upper = 1)
  in.diffus<-(1-s.o)^(5/3)*b.in$value

  sigma<-((5*n*(c-3)*lambda^2*Ksat*Matrix.pot*in.diffus*m.r)/12*pi*gamma.depth^2)^(1/3)*(1-s.o)^(2/3)

  Runoff.a<-mPa*(exp(-G-2*sigma)*gamma(sigma+1)*sigma^(-sigma)+((m.t*Ksat)/mPa)*s.o^C)
}
