#this is another form of eaglesons water balance equations where I 
#use the optimize function in R to optimize the parameters
M<-0.4
s.o<-seq(0.1,0.9,0.001)

p<-seq(0.1,0.99,0.1)
y<-function(p=p,x=x.r) -x^2+p*1/exp(x^3)
optimx(s.o,water)
water<-function(M=M,s.o=s.o,type="clayloam"){
m.n<-15.7 # mean number of storms per year # obtained from Eagleson: Climate, Soil, and Vegetation
e.p<-0.27 # average bare soil potential evaporation rate
m.Pa<-54.4 #mean annual precipitation in cm
m.r<-1.4 #mean storm depth--double check value
m.b<-10.42 #in units of days
m.t<-212 #mean rainy season
alpha<-1/m.b #one over mean time between storms
#Eagleson (2002) has poisson pulse parameters for a number of stations in Appendix F
#for Clinton, MA, I used parameters for Boston (Station ID 14 in eagleson
lambda<-1.86*2.54 #may be incorrect, changed inches to cm #scale parameter for gamma distribution of storm depths (1/cm)
gamma.depth<-0.25 

if(type =="clay"){
  Ksat<-0.72 #k saturation
  Matrix.pot<-25 #soil matrix potential
  c<-12 #for clay
  n<-0.45 #soil porosity, for clay
  
  h.o<-0.1 # small constant value for surface water retention
  k.v<-0.937/e.p# unstressed transpiration
  #M<-0.5 # ranges from 0 to 1 (but cant be one)
  w<-0 #assume no capillary rise
}else{
  if(type =="clayloam"){
    Ksat<-2.0 #k saturation
    Matrix.pot<-19 #soil matrix potential
    c<-10 #for clayloam
    n<-0.35 #soil porosity, for clay
    
    ##these values set to constant for this purpose
    h.o<-0.1 # small constant value for surface water retention
    k.v<-0.937/e.p# unstressed transpiration
    #M<-0.5 # ranges from 0 to 1 (but cant be one)
    w<-0 #assume no capillary rise
  }else{
    if(type =="siltloam"){
      Ksat<-8.64 #k saturation
      Matrix.pot<-166 #soil matrix potential
      c<-6 #for clayloam
      n<-0.35 #soil porosity, for clay
      
      ##these values set to constant for this purpose
      h.o<-0.1 # small constant value for surface water retention
      k.v<-0.937/e.p# unstressed transpiration
      #M<-0.5 # ranges from 0 to 1 (but cant be one)
      w<-0 #assume no capillary rise
    }else{
      if(type =="sandyloam"){
        Ksat<-18 #k saturation
        Matrix.pot<-200 #soil matrix potential
        c<-4 #for clayloam
        n<-0.25 #soil porosity, for clay
        
        ##these values set to constant for this purpose
        h.o<-0.1 # small constant value for surface water retention
        k.v<-0.937/e.p# unstressed transpiration
        #M<-0.5 # ranges from 0 to 1 (but cant be one)
        w<-0 #assume no capillary rise
        
      }}}
  
  #k.v<-0.097/e.p #unstressed composite transpiration rate/e.p
  #k.v<-1
  #####################
  #####Beta.s model####
  #####################
  #bets.s model needs k.v and s.o
  require(gsl)
  gamma.ratio<-gamma_inc(gamma.depth, lambda*h.o)/gamma(gamma.depth)
  
  #define beta.s function
  
}
beta.s<-function(M, s.o){
  B<- ((1-M)/(1+M*k.v-w/e.p))+((k.v*(M^2)+(1-M)*w/e.p)/((2*(1+M*k.v-w/e.p)^2)))
  
  C<- 1/2*(M*k.v-w/e.p)^(-2)
  
  E<-((alpha*n*(c-3)*Ksat*Matrix.pot*ex.diffus)/(pi*e.p^2))*s.o^((c+5)/2)
  
  beta<-(gamma.ratio)-((1+(alpha*h.o/e.p)/lambda*h.o)^(-1*gamma.depth))*((gamma_inc(gamma.depth, (lambda*h.o+alpha*h.o/e.p)))/gamma(gamma.depth))*exp(-B*E)+
    (1+gamma.ratio)*(1-exp(-B*E-alpha*h.o/e.p)*(1+M*k.v+(2*B)^(1/2)*E-w/e.p)+
                       exp(-C*E-alpha*h.o/e.p)*(M*k.v +(2*C)^(1/2)*E-w/e.p)+
                       ((2*E)^(1/2))*exp(-alpha*h.o/e.p)*(gamma_inc(3/2, C*E)-gamma_inc(3/2,B*E)) + 
                       ((1+(1+(alpha*h.o/e.p)/lambda*h.o))^(-1*gamma.depth))*(gamma_inc(gamma.depth, lambda*h.o+(alpha*h.o/e.p))/gamma(gamma.depth))*(sqrt(2*E)*(gamma_inc(3/2, C*E)-gamma_inc(3/2,B*E))+
                                                                                                                                                        exp(-C*E)*(M*k.v +(2*C)^(1/2)*E-w/e.p) -
                                                                                                                                                        exp(-B*E)*(M*k.v+(2*B)^(1/2)*E-w/e.p)))
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
  mPa<-111.3 # mean annual precipitation in mm, not sure if this needs to be changed to cm
  G<-(Ksat*(m.n*m.r)/mPa)*(((1+s.o^c)/2)-(w/Ksat))
  
  integrand.in<-function(x){x^((c+1)/2)*(x-s.o)^(2/3)}
  b.in<-integrate(integrand.in, lower = s.o, upper = 1)
  in.diffus<- (1-s.o)^(5/3)*b.in$value
  
  sigma<-((5*n*(c-3)*lambda^2*Ksat*Matrix.pot*in.diffus*m.r)/12*pi*gamma.depth^2)^(1/3)*(1-s.o)^(2/3)
  
  Runoff.a<-mPa*(exp(-G-2*sigma)*gamma(sigma+1)*sigma^(-sigma)+((m.t*Ksat)/m.Pa)*s.o^c)
}


#basic ET.model.a functions repeated from above
ET.model.a<-function(M, s.o){
  ET.model.a<-((m.n*e.p)/alpha)*(1-M)*beta.s(M,s.o) + M*k.v
}
(ET.model.a(M,s.o)+Runoff(s.o)-m.Pa)^2

}

#this works
for(v in s.o){
  print(optimize(water,upper=range(s.o)[2],lower=range(s.o)[1],s.o=v))}

#you could also apply the optimize function over


ll <- lapply(s.o,function(i)
  (optimize(water,upper=range(s.o)[2],lower=range(s.o)[1], s.o=i)))
ll[[which.min(sapply(ll,'[[','minimum'))]]
#but this is taking the minmum soil moisture value given, not sure

#psoptim
lapply(s.o, function(i)
(psoptim(par=rep(0.5,2), water, gr = NULL, lower = -1, upper = 1, control = list()))
