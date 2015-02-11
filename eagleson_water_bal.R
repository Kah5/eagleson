#Eagleson Water Balance model with s.o and M as state variables
#This is coded up with parameters from Clinton, MA clay parameters
#Need to check lambda parameters and gamma.depth


###Goal: 
###find s.o and M such that Pr.modeled = annual Precip
###in other words, such that ET.model.a(M,kv)-Runoff.model(s.o)-Annual Precip=0
M.1<-as.list(runif(100)) #samples from uniform distribution to get values from 0 to 1
s.o<-M.1 #take the same values for s.o (also from 0 to 1)
#set.seed(10)

#create matrix of all possible M and s.o combinations to run in the water balance function
m.so<-expand.grid(M.1,s.o)


#basic water balance
#water.bal<-function(M,s.o){
#abs(ET.model.a(M)+Runoff(s.o)-mPa)
#}


#ET.model.a<-((m.n*e.p)/alpha)*(1-M)*beta.s(M,k.v) + M*k.v

#m.n=mean number rainstorms per year
#Beta.s=avg bare soil evap efficiency
#M=vegetation canopy denisity
#ET.model.a=evaporation component of water balance

#treat M and k.v as state variables

water.bal<-function(M,s.o){
#k.v<-0.097/e.p #unstressed composite transpiration rate/e.p

k.v<-1


#####################################
#define site based model parameters##
#####################################
m.n<-15 # mean number of storms per year
e.p<-0.15 # average bare soil potential evaporation rate
m.Pa<-111.3 #mean annual precipitation in cm
m.r<-0.32 #mean storm depth--double check value
m.b<-3.0 #in units of days
m.t<-365 #mean rainy season
alpha<-1/m.b #one over mean time between storms
lambda<-0.5 #may be incorrect
gamma.depth<-0.5 #these may be incorrect

Ksat<-0.72 #k saturation
Matrix.pot<-25 #soil matrix potential
c<-12 #for clay
n<-0.45 #soil porosity, for clay
h.o<-0.1 # small constant value for surface water retention
k.v<-1# unstressed transpiration
M<-0.5 # ranges from 0 to 1 (but cant be one)
w<-0 #assume no capillary rise

#####################
#####Beta.s model####
#####################
#bets.s model needs k.v and s.o
require(gsl)
gamma.ratio<-gamma_inc(gamma.depth, lambda*h.o)/gamma(gamma.depth)

#define beta.s function
beta.s<-function(M, s.o){
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



##########################################
#define exfiltration diffusivity & runoff#
##########################################


#Runoff Function:

Runoff<-function(s.o){
## define the integrated function
  integrand.ex <- function(x) {x^((c+1)/2)*(s.o-x)^(0.85)}
## integrate the function from 0 to s.o
#b$value provides the value of the intergral
  ex.diffus<-s.o^((c+1)/2)*(1.85*s.o^(-1.85)*integrate(integrand.ex, lower = 0, upper = s.o)$value) #in.diffus and out.diffus are incorrectly specified
###Runoff model
  mPa<-1113 # mean annual precipitation in mm, not sure if this needs to be changed to cm
  G<-(Ksat*(m.n*m.r)/mPa)*((1+s.o^c)/2)-w/Ksat

  integrand.in<-function(x){x^((c+1)/2)*(x-s.o)^(2/3)}
  b.in<-integrate(integrand.in, lower = s.o, upper = 1)
  in.diffus<-(1-s.o)^(5/3)*b.in$value

  sigma<-((5*n*(c-3)*lambda^2*Ksat*Matrix.pot*in.diffus*m.r)/12*pi*gamma.depth^2)^(1/3)*(1-s.o)^(2/3)

  Runoff.a<-mPa*(exp(-G-2*sigma)*gamma(sigma+1)*sigma^(-sigma)+((m.t*Ksat)/mPa)*s.o^c)
}


#basic ET.model.a functions repeated from above
ET.model.a<-function(M,s.o){
  ET.model.a<-((m.n*e.p)/alpha)*(1-M)*beta.s(M,s.o) + M*k.v
}


  abs(ET.model.a(M,s.o)+Runoff(s.o)-mPa)
  
}


test3<-mapply(water.bal, m.so$Var1, m.so$Var2) # this is to test if R can handle the large m.so

#test3 now contains the output of water.bal for all combinations of 100x100 values between 0 and 1
#to determine which values of M and s.o close the water balance, we need to keep values where abs(ET.model.a(M,s.o)+Runoff(s.o)-mPa) is ~ 0

#lets say we are willing to accept M and s.o values that fall within +/-5% of the Mean annual precidp
wb<-cbind(m.so, test3)
wb.small<-wb[wb$test3 < (mPa*0.005),]

plot(wb.small$Var1,wb.small$Var2,)

plot(wb.small$Var1,wb.small$Var2, main="Isoclines of parameter combinations (M,s.o) that satisfiy water balance closre for Clinton, MA", xlab="M, Canopy density", ylab="Equlibrium soil moisture, s.o")
