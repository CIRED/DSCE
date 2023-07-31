# -*- coding: utf-8 -*-

#definition of CES function
if subs==1:
    def constant_elasticity_substitution(c1,c2):
        CES=(c1**(1-gamma))*(c2**gamma)
        return CES;
else:
    def constant_elasticity_substitution(c1,c2):
        CES=((1-gamma)*c1**((subs-1)/subs)+gamma*(c2**((subs-1)/subs)))**(subs/(subs-1))
        return CES;

#definition of utility function
if theta==1:
    def utility(consumption,carbon_stock):
        c=maximum(consumption,10**(-8)*ones(shape(consumption)))
        E=environment_from_temperature(temperature(carbon_stock))
        CES=constant_elasticity_substitution(c,E)
        u=log(CES)
        return u;
else:
    def utility(consumption,carbon_stock):
        c=maximum(consumption,10**(-8)*ones(shape(consumption)))
        E=environment_from_temperature(temperature(carbon_stock))
        CES=constant_elasticity_substitution(c,E)
        #change utility f
        #u=(CES**(1-theta)-1)/(1-theta)
        u=CES**(1-theta)/(1-theta)
        return u;

def marginal_utility(consumption,carbon_stock):
    c=maximum(consumption,10**(-8)*ones(shape(consumption)))
    E=environment_from_temperature(temperature(carbon_stock))
    CES=constant_elasticity_substitution(c,E)
    du=CES**(-theta)*(1-gamma)*(c/CES)**(-1/subs)
    return du;

#total production without damages per period
def production(TFP,capital,labor):
    y=deltaT*TFP*(labor**(1-alpha))*capital**alpha
#    y=TFP*(labor**(1-alpha))*capital**alpha
    return y;

#abatement cost in percentage GDP
def abatement_cost(ab_level,total_cost):
    ab=maximum(ab_level,0)
    ab=minimum(ab,1)
    b=total_cost*ab**theta2    
    return b;

#marginal abatement cost in $/tC02
def marginal_abatement_cost(ab_level,total_cost,carbon_intensity):
    ab=maximum(ab_level,0)
    ab=minimum(ab,1)
    b=total_cost*theta2*ab**(theta2-1)/carbon_intensity    
    return b;
    
def emissions(ab_level,production,carbon_intensity):#emissions per period
    ab=maximum(ab_level,0)
    ab=minimum(ab,1)
    e=carbon_intensity*(1-ab)*production
    return e;

def temperature(cum_e):
    v=beta_temp*cum_e
    return v;

#state of the environment as a function of temperature
def environment_from_temperature(temp):
    E=E0/(1+ae*temp**2)
    return E;

def damage_factor(temp,eps):#damages factor, correcting TFP, from temperature and shock (caution: positive shock=more productivity=less damages)
    #omega=maximum(10**(-8), 1 /(1+pi1*temp+pi2*temp**pi3)-eps)
    #New damage factor from Nordhaus (2016)
    omega=maximum(10**(-8), 1-(pi1*temp+pi2*temp**pi3)-eps)
    return omega;
