# -*- coding: utf-8 -*-


#functions for Bellman equation in the case of one dimension (capital)

#function to maximise
#x is control 2D-array of length lx
#u is state 2D-array of size lx*dim
#t is time 1D-array of size 1 or lx
#eps is stochastic choc
#v is maximand value of size lx*1
def maximand(t,x,u,eps):
    #x is proportion of consumed production
    #u is capital per efficient capita
    x=maximum(x,zeros(shape(x)))
    x=minimum(x,ones(shape(x)))
    k=(A[t,:]**(1/(1-alpha)))*L[t,:]*u
    conso=production(A[t,:],k,L[t,:])/(L[t,:])*x
    s=s0
    v=1/L[[0],:]*L[t,:]*utility(conso,s)
    return v;

#partial derivative with respect to total consumption of maximand
#x is control 2D-array of length lx
#u is state 2D-array of size lx*dim
#eps is stochastic choc of mean 0
#t is time 1D-array of size 1 or lx
#v is maximand value of size lx*1
def partialc_maximand(t,x,u,eps):
    x=maximum(x,zeros(shape(x)))
    x=minimum(x,ones(shape(x)))
    k=(A[t,:]**(1/(1-alpha)))*L[t,:]*u
    conso=production(A[t,:],k,L[t,:])/(L[t,:])*x
    s=s0
    v=1/L[[0],:]*marginal_utility(conso,s)
    return v;


#law of motion of state variable
#t is time (scalar)
#x is control 2D-array of size lx*dim
#u is state 2D-array of size lx*dim
#eps is stochastic choc
#un is next state variable of size lx*dim
def law_motion(t,x,u,eps):
    #x is proportion of consumed production
    #u is capital per efficient capita
    x=maximum(x,zeros(shape(x)))
    x=minimum(x,ones(shape(x)))
    k=(A[[t-1],:]**(1/(1-alpha)))*L[[t-1],:]*u
    kn=(1-delta)*k+production(A[[t-1],:],k,L[[t-1],:])*(1-x)
    un=kn/((A[[t],:]**(1/(1-alpha)))*L[[t],:])
    return un;

#terminal scrap-value
#sum of the utilities from T+1 to infty
def terminal_value(t,u):
    if not(t==T+1):
        print('Are you sure terminal_value function is correctly called?')
    k=(A[[t-1],:]**(1/(1-alpha)))*L[[t-1],:]*u
    kn=(A[[t],:]**(1/(1-alpha)))*L[[t],:]*u
    #consumption for constant capital per efficient capita
    conso=(production(A[[t-1],:],k,L[[t-1],:])+(1-delta)*k-kn)/(L[[t-1],:])
    #adjustement in terminal constraint adapted from Barr Manne (GA is growth rate per year, not per period) 
    s=s0
    v=1/L[[0],:]*L[[t-1],:]*1/(1-(bet)*((1+GA[[t-1],:])**(deltaT*(1-theta)/(1-alpha))))*utility(conso,s)
    return v;
    
# function that defines general risk-sensitive preferences
def Bommier_function2(u1,u2,pr):
    return (-(1/epsilon_b)*log((1-pr)*exp(-epsilon_b*u1)+pr*exp(-epsilon_b*u2)))

# function that defines general risk-sensitive preferences
def Bommier_function(u1,u2,pr):
    return (-(bet/epsilon_b)*log((1-pr)*exp(-epsilon_b*u1)+pr*exp(-epsilon_b*u2)))

