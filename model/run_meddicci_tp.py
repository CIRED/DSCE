# -*- coding: utf-8 -*-

#Reading parameters and calculating exogenous trends
exec(open(model_folder+'read_parameters.py').read())

#Calling functions
if tp==0:
    exec(open(model_folder+'f_model.py').read()) #IAM functions

elif tp==1:
    exec(open(model_folder+'f_model_tp.py').read()) #IAM functions

# Call functions for approximation with Chebyshev polynomials
# For the discounted expected utility model
if epsilon_b==0:
    if tp==0:
        exec(open(model_folder+'f_approx_chebyshev.py').read())
    elif tp==1:
        exec(open(model_folder+'f_approx_chebyshev_tp.py').read())

# For risk-sensitive preferences
else:
    exec(open(model_folder + 'f_approx_chebyshev_tp_Bommier.py').read())

#Functions for Bellman equation for various dimensions of the problem
if dim==1:
    exec(open(model_folder+'f_bellman_1.py').read())
elif dim==2:
    if tp==0:
        exec(open(model_folder+'f_bellman_2.py').read())
    elif tp==1:
        exec(open(model_folder+'f_bellman_2_tp.py').read())
elif dim==3:
    exec(open(model_folder+'f_bellman_3.py').read())

#Defining approximation interval
exec(open(model_folder+'approx_interval.py').read())

#Saving the boundaries of interpolation
savetxt(run_folder+'u_min.csv', u_min, delimiter=";", header='inferior boundary of interpolation dim '+str(dim)+' stochastic '+str(stochastic))
savetxt(run_folder+'u_max.csv', u_max, delimiter=";", header='superior boundary of interpolation dim '+str(dim)+' stochastic '+str(stochastic))

#Initializing some value (first guess of controls, u0 for states)
if dim==1:
    guess0=0.5*ones((1,1))
    u0=u0trend[[0],:]
if dim==2:
    guess0=array((0.5,0.1)).reshape(1,2)
    u0=concatenate((K0/((A0**(1/(1-alpha)))*POP0)*ones((1,1)),s0*ones((1,1))),axis=1)
if dim==3:
    guess0=array((0.5,0.1)).reshape(1,2)
    u0=concatenate((K0/((A0**(1/(1-alpha)))*POP0)*ones((1,1)),s0*ones((1,1)),A0*ones((1,1))),axis=1)
    
#Interpolation to compute coefficients of Chebyshev polynomials
if tp==0:
    exec(open(model_folder+'interpolate.py').read()) 
    savetxt(run_folder+'coef_V.csv', coef_V, delimiter=";", header='coef_V dim '+str(dim)+' stochastic '+str(stochastic))     

elif tp==1:
    u_threshold_max=threshold_max/beta_temp
    u_threshold_min=threshold_min/beta_temp
    u_max_tp=genfromtxt(run_folder+'/u_max.csv', dtype=float, delimiter=';')
    u_min_tp=genfromtxt(run_folder+'/u_min.csv', dtype=float, delimiter=';')
    u_upper_pretp=genfromtxt(run_folder+'/u_max.csv', dtype=float, delimiter=';')

    for i in list(range(0,T)):
        u_min_tp[i,1]=max(u_threshold_min,u_min_tp[i,1])
        u_max_tp[i,1]=min(u_threshold_max,u_max_tp[i,1])
        u_upper_pretp[i,1]=min(u_threshold_min,u_upper_pretp[i,1])

    u_max_posttipping=u_max.copy()
    u_min_posttipping=u_min_tp.copy()
    u_max_uncertain=u_max_tp.copy()
    u_min_uncertain=u_min_tp.copy()
    u_max_pretipping=u_upper_pretp.copy()
    u_min_pretipping=u_min.copy()

    savetxt(run_folder+'u_min_pretipping.csv', u_min_pretipping, delimiter=";", header='min value of tipping poing in cumu emissions '+str(dim)+' stochastic '+str(stochastic))
    savetxt(run_folder+'u_max_pretipping.csv', u_max_pretipping, delimiter=";", header='max value of tipping poing in cumu emissions '+str(dim)+' stochastic '+str(stochastic))
    savetxt(run_folder+'u_min_posttipping.csv', u_min_posttipping, delimiter=";", header='min value of tipping poing in cumu emissions '+str(dim)+' stochastic '+str(stochastic))
    savetxt(run_folder+'u_max_posttipping.csv', u_max_posttipping, delimiter=";", header='max value of tipping poing in cumu emissions '+str(dim)+' stochastic '+str(stochastic))
    savetxt(run_folder+'u_min_uncertain.csv', u_min_uncertain, delimiter=";", header='min value of tipping poing in cumu emissions '+str(dim)+' stochastic '+str(stochastic))
    savetxt(run_folder+'u_max_uncertain.csv', u_max_uncertain, delimiter=";", header='max value of tipping poing in cumu emissions '+str(dim)+' stochastic '+str(stochastic))

    exec(open(model_folder+'interpolate_tp_posttipping.py').read())
    savetxt(run_folder+'coef_V_posttipping.csv', coef_V_posttipping, delimiter=";", header='coef_V dim '+str(dim)+' stochastic '+str(stochastic))         
 
#Two different possibilities: either a certain tipping point or an incertain tipping point
    #First interpolation in uncertain zone
    exec(open(model_folder+'interpolate_tp_uncertain.py').read())
    savetxt(run_folder+'coef_V_uncertain.csv', coef_V_uncertain, delimiter=";", header='coef_V dim '+str(dim)+' stochastic '+str(stochastic))

    #Second interpolation in pre-tipping zone
    exec(open(model_folder+'interpolate_tp_pretipping.py').read())
    savetxt(run_folder+'coef_V_pretipping.csv', coef_V_pretipping, delimiter=";", header='coef_V dim '+str(dim)+' stochastic '+str(stochastic))

#We will do stochastic draws later
nbofdraws=draws
draws=1

#Simulate (1 draw)
if tp==0:
    exec(open(model_folder+'simulate.py').read())
elif tp==1:
    exec(open(model_folder+'simulate_tp.py').read())

savetxt(run_folder+'state_V.csv', u_m, delimiter=";", header='state variables dim '+str(dim)+' stochastic '+str(stochastic))
savetxt(run_folder+'control.csv', x_m, delimiter=";", header='control variables dim '+str(dim)+' stochastic '+str(stochastic))

