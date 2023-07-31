# -*- coding: utf-8 -*-

#functions for Bellman equation in the case of two dimensions (capital,stock of emissions)


#function to maximise
#x is control 2D-array of length lx
#u is state 2D-array of size lx*dim
#t is time 1D-array of size 1 or lx
#eps is stochastic choc of mean 0
#v is maximand value of size lx*1
def maximand(t,x,u,eps):
    x=maximum(x,zeros(shape(x)))
    x=minimum(x,ones(shape(x)))
    c=x[:,[0]]#c is relative consumption
    a=x[:,[1]]#a is relative abatement
    k=(A[t,:]**(1/(1-alpha)))*L[t,:]*u[:,[0]]#u[:,[0]] is capital per efficient capita
    s=u[:,[1]]#total emissions
    conso=damage_factor(temperature(s),eps)*(1-abatement_cost(a,theta1[t,:]))*production(A[t,:],k,L[t,:])/L[t,:]*c
    v=1/L[[0],:]*L[t,:]*utility(conso,s)
    return v;

#partial derivative with respect to total consumption of maximand
#u is state 2D-array of size lx*dim
#eps is stochastic choc of mean 0
#t is time 1D-array of size 1 or lx
#v is maximand value of size lx*1
def partialc_maximand(t,x,u,eps):
    x=maximum(x,zeros(shape(x)))
    x=minimum(x,ones(shape(x)))
    c=x[:,[0]]#c is relative consumption
    a=x[:,[1]]#a is relative abatement
    k=(A[t,:]**(1/(1-alpha)))*L[t,:]*u[:,[0]]#u[:,[0]] is capital per efficient capita
    s=u[:,[1]]#total emissions
    conso=damage_factor(temperature(s),eps)*(1-abatement_cost(a,theta1[t,:]))*production(A[t,:],k,L[t,:])/L[t,:]*c
    v=1/L[[0],:]*marginal_utility(conso,s)
    return v;


#law of motion of state variable
#t is time (scalar)
#x is control 2D-array of length lt
#u is state 2D-array of size lt*dim
#eps is stochastic choc
#un is next state variable of size lt*dim
def law_motion(t,x,u,eps):
    x=maximum(x,zeros(shape(x)))
    x=minimum(x,ones(shape(x)))
    c=x[:,[0]]#c is relative consumption
    a=x[:,[1]]#a is relative abatement
    k=(A[[t-1],:]**(1/(1-alpha)))*L[[t-1],:]*u[:,[0]]#u[:,[0]] is capital per efficient capita
    s=u[:,[1]]#total emissions
    y=production(A[[t-1],:],k,L[[t-1],:])*damage_factor(temperature(s),eps)
    kn=((1-delta))*k*damage_factor(temperature(s),eps)**(damage_capital/alpha)+y*(1-abatement_cost(a,theta1[[t-1],:]))*(1-c)
    sn=s+emissions(a,y,sigm[[t-1],:])
    un=concatenate((kn/((A[[t],:]**(1/(1-alpha)))*L[[t],:]),sn),axis=1)
    return un;

#partial derivative of law_motion with respect to s
def partials_law_motion(t,x,u,eps):
    x=maximum(x,zeros(shape(x)))
    x=minimum(x,ones(shape(x)))
    c=x[:,[0]]#c is relative consumption
    a=x[:,[1]]#a is relative abatement
    k=(A[[t-1],:]**(1/(1-alpha)))*L[[t-1],:]*u[:,[0]]#u[:,[0]] is capital per efficient capita
    s=u[:,[1]]#total emissions
    y=damage_factor(temperature(s),eps)*production(A[[t-1],:],k,L[[t-1],:])
    dy=partialt_damage_factor(temperature(s),eps)*partials_temperature(s)*production(A[[t-1],:],k,L[[t-1],:])
    dkn=((1-delta))*k*(damage_capital/alpha)*damage_factor(temperature(s),eps)**(damage_capital/alpha-1)*partialt_damage_factor(temperature(s),eps)*partials_temperature(s)+(1-abatement_cost(a,theta1[[t-1],:]))*dy*(1-c)
    dsn=1+partialy_emissions(a,y,sigm[[t-1],:])*dy
    dun=concatenate((dkn/((A[[t],:]**(1/(1-alpha)))*L[[t],:]),dsn),axis=1)
    return dun;


#terminal scrap-value
#sum of the utilities from T+1 to infty
def terminal_value(t,u):
    if not(t==T+1):
        print('Are you sure terminal_value function is correctly called?')
    kpec=u[:,[0]]
    s=u[:,[1]]
    k=(A[[t-1],:]**(1/(1-alpha)))*L[[t-1],:]*kpec
    kn=(A[[t],:]**(1/(1-alpha)))*L[[t],:]*kpec#constant capital per efficient capita
    #consumption for constant capital per efficient capita and total abatment
    eps=0
    a=1
    #a=0
    conso=(damage_factor(temperature(s),eps)*(1-abatement_cost(a,theta1[[t-1],:]))*production(A[[t-1],:],k,L[[t-1],:])+((1-delta))*k*damage_factor(temperature(s),eps)**(damage_capital/alpha)-kn)/L[[t-1],:]
    #adjustement in terminal constraint adapted from Barr Manne (GA is growth rate per year, not per period) 
    v=1/L[[0],:]*L[[t-1],:]*1/(1-bet*((1+GA[[t-1],:])**(deltaT*(1-theta)/(1-alpha))))*utility(conso,s)
    return v;


# function that defines general risk-sensitive preferences
def Bommier_function(u1,u2,pr):
    return (-(bet/epsilon_b)*log((1-pr)*exp(-epsilon_b*u1)+pr*exp(-epsilon_b*u2)))

