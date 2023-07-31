# -*- coding: utf-8 -*-

#importing modules
from numpy import *
import scipy.integrate as si
import csv
from scipy.optimize import fmin
from scipy.optimize import fmin_cg
from scipy.optimize import fmin_bfgs
from operator import itemgetter
import time
import math
import pyparsing as pypars
import matplotlib.pyplot as plt
import os
import functools
import copy
import numpy as np
import sys
import pandas as pd

#defining the paths for folders                         
pathsep=os.sep
preprod_folder='..'+pathsep+'run'+pathsep
model_folder='..'+pathsep+'model'+pathsep
param_folder='..'+pathsep+'parameters'+pathsep
outputs_folder='..'+pathsep+'outputs'+pathsep
postprod_folder='..'+pathsep+'postproduction'+pathsep
figures_folder='..'+pathsep+'figures'+pathsep

#defining the study and the set of runs to compute or analyze
exec(open(preprod_folder+'study.py').read())

#reading the study
table_study=loadtxt(preprod_folder+study+'.csv',dtype=str,delimiter=';')
table_param_names=table_study[0,2:]

#run identification numbers available in name_of_study.csv
run_id=loadtxt(preprod_folder+study+'.csv',dtype=int,delimiter=';',skiprows=1,usecols=(0,))
if table_study[0,1]!='prev_run':
	print("You need to define the previous runs in the second column of the study file, please correct file "+study+".csv")
	exit()
prev_run_id=loadtxt(preprod_folder+study+'.csv',dtype=int,delimiter=';',skiprows=1,usecols=(1,))
table_param_values_temp=loadtxt(preprod_folder+study+'.csv',dtype=None,delimiter=';',skiprows=1)
table_param_values=table_param_values_temp[:,2:]
#nb_digits_run_id=len(str(max(run_id)))
nb_digits_run_id=4

#checking if the study is well-defined
exec(open(preprod_folder+'check_study.py').read())

if well_defined==0:
    print("The study is not well defined, please correct file "+study+".csv according to the errors messages above.")
    exit()

#restricting the operations to run whose id number is available (otherwise the run is not defined)
to_run=intersect1d(to_run,run_id)
to_analyze=intersect1d(to_analyze,run_id)

if to_run.any!=[0]: 
    #checking the runs that have already been computed
    to_run_still=[]
    already_run=[]
    runs_outputs=array([None] * shape(run_id)[0])
    list_outputs=os.listdir(outputs_folder)
    for run_no in to_run:
        already_done= [s for s in list_outputs if study+'_run'+str(run_no).zfill(nb_digits_run_id) in s]
        if not already_done:
            to_run_still=to_run_still+[run_no]
        elif shape(already_done)[0]!=1:
            print("The outputs folder contains more than one output folder for run number "+str(run_no)+". Please clean the outputs folder.")
            exit()
        else:
            already_run=already_run+[run_no]
            runs_outputs[run_id==run_no]=already_done[0]
    if already_run!=[]:
        if size(already_run)==1:
            print("Run "+str(already_run)+" has already been computed")
        else:
            print("Runs "+str(already_run)+" have already been computed")
    if to_run_still!=[]:
        if size(already_run)==1:
            print("Only run "+str(to_run_still)+" will be computed now")
        else:
            print("Only runs "+str(to_run_still)+" will be computed now")
    
    #running
    for run_no in to_run_still:	
       
        #finding the previous run results to use for approximation interval definition if dim is not 1
        dim=int(table_param_values[run_id==run_no,table_param_names=='dim'])
        stochastic=int(table_param_values[run_id==run_no,table_param_names=='stochastic'])
        if dim!=1:
            prev_run=int(prev_run_id[run_id==run_no])
            #checking if the output folder for previous run exists
            list_outputs=os.listdir(outputs_folder)
            already_done= [s for s in list_outputs if study+'_run'+str(prev_run).zfill(nb_digits_run_id) in s]
            if not already_done:
                print("Cannot compute run number "+str(run_no))
                print("Previous run outputs not available, please compute before the previous run, run number "+str(prev_run))
                continue
            elif shape(already_done)[0]!=1:
                print("There are more than one output folder for previous run, run number "+str(prev_run)+". Please clean the outputs folder.")
                continue
            else:
                prev_outputs= already_done[0]
                if os.path.isfile(outputs_folder+prev_outputs+pathsep+'simulate_break.csv'):
                    print("Cannot compute run number "+str(run_no))
                    print("Simulation of previous run (run number "+str(prev_run)+") has broken")
                    continue
                else:
                    prev_run_folder=outputs_folder+prev_outputs+pathsep
                    
#get the corresponding no damage scenario (idealine) and check if computed
            ideal_run=int(ideal_run_id[run_id==run_no])
            #checking if the output folder for idealine run exists
            list_outputs=os.listdir(outputs_folder)
            already_done= [s for s in list_outputs if study+'_run'+str(ideal_run).zfill(nb_digits_run_id) in s]
            if not already_done:
                print("Baseline run outputs not available, please compute before the idealine run, run number "+str(ideal_run))
                continue
            elif shape(already_done)[0]!=1:
                print("There are more than one output folder for idealine run, run number "+str(ideal_run)+". Please clean the outputs folder.")
                continue
            else:
                ideal_outputs= already_done[0]
                if os.path.isfile(outputs_folder+ideal_outputs+pathsep+'simulate_break.csv'):
                    print("Simulation of idealine run (run number "+str(ideal_run)+") has broken")
                    continue
                else:
                    ideal_run_folder=outputs_folder+ideal_outputs+pathsep                    
         
        #defining the run folder (where parameters are stored and results will be writen)
        current_time=time.strftime("%Y-%m-%d")+'-'+time.strftime("%Hh%M")
        run_folder=outputs_folder+study+'_run'+str(run_no).zfill(nb_digits_run_id)+'-'+current_time+pathsep
        os.makedirs(run_folder)
        
        #storing the parameters in the run folder
        exec(open(preprod_folder+'create_parameters_files.py').read())
                    
        #running the model
        print('Executing run number '+str(run_no))    
        #execfile(model_folder+'run_meddicci.py')
        
        exec(open(model_folder+'run_meddicci_tp.py').read())
        #execfile(model_folder+'run_meddicci_marginal.py')
        
        print('Run number '+str(run_no)+' finished')


#rechecking the runs that have already been computed
to_run_still=[]
already_run=[]
runs_outputs=array([None] * shape(run_id)[0])
list_outputs=os.listdir(outputs_folder)
for run_no in to_analyze:
    already_done= [s for s in list_outputs if study+'_run'+str(run_no).zfill(nb_digits_run_id) in s]
    if not already_done:
        to_run_still=to_run_still+[run_no]
    elif shape(already_done)[0]!=1:
        print("The outputs folder contains more than one output folder for run number "+str(run_no)+". Please clean the outputs folder.")
        exit()
    else:
        already_run=already_run+[run_no]
        runs_outputs[run_id==run_no]=already_done[0]
if already_run!=[]:
    if size(already_run)==1:
        print("Run "+str(already_run)+" has already been computed and will be analyzed")
    else:
        print("Runs "+str(already_run)+" have already been computed and will be analyzed")
if to_run_still!=[]:
    if size(to_run_still)==1:
        print("Run "+str(to_run_still)+" has not been computed yet, so its analysis will not be possible") 
    else:
        print("Runs "+str(to_run_still)+" have not been computed yet, so their analysis will not be possible")

#analyzing results
for run_no in to_analyze:
    #getting the run outputs folder
    run_outputs=runs_outputs[run_id==run_no][0]
    if not run_outputs:
        print("Cannot analyze run number "+str(run_no)+", it has not been run")
    else:
        run_folder=outputs_folder+run_outputs+pathsep
        if os.path.isfile(run_folder+'simulate_break.csv'):
            print("Run number "+str(run_no)+" has broken")
        exec(open(model_folder+'read_parameters.py').read())
        
        
        if tp==0:
            exec(open(model_folder+'f_model.py').read()) #IAM functions    
            exec(open(model_folder+'f_approx_chebyshev.py').read()) #functions for approximation with chebyshev polynomials'''

        
        elif tp==1:
            if threshold_min==threshold_max:
                exec(open(model_folder+'f_approx_chebyshev_tp_certain.py').read())         
            else:
                if theta==gamma_EZ:
                    exec(open(model_folder+'f_approx_chebyshev_tp.py').read())
                else:                    
                    exec(open(model_folder+'f_approx_chebyshev_tp_EZ.py').read())

            
            if deterministic==0:
                exec(open(model_folder+'f_model_tp.py').read()) #IAM functions  
                #execfile(model_folder+'f_model_tp_marginal.py') #IAM functions
            elif deterministic==1:    
                exec(open(model_folder+'f_model_deterministic.py').read()) #IAM functions    
                #execfile(model_folder+'f_model_deterministic_marginal.py')
            
#saving trends up to T
        savetxt(run_folder+'time.csv', deltaT*time_horizon[0:T,:], delimiter=";", header='time in years dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'population.csv', L[0:T,:], delimiter=";", header='population  dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'growth.csv', GA[0:T,:], delimiter=";", header='annual productivity growth at the beginning of period dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'carbon_intensity.csv', sigm[0:T,:], delimiter=";", header='CO2 intensity of GDP dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'abat_cost_factor.csv', theta1[0:T,:], delimiter=";", header='total abatment cost dim '+str(dim)+' stochastic '+str(stochastic))

#getting the results
        u=genfromtxt(run_folder+'state_V.csv', dtype=float, delimiter=';').reshape(T,dim)
        kpec=u[:,[0]]
        x=genfromtxt(run_folder+'control.csv', dtype=float, delimiter=';')
        if dim==1:
            productivity=A[0:T,:]
            x=x.reshape(T,1)
            abat=zeros((T,1))
        if dim==2:
            cumulatedemissions=u[:,[1]]
            productivity=A[0:T,:]
            abat=x[:,[1]]
        if dim==3:
            cumulatedemissions=u[:,[1]]
            productivity=u[:,[2]]
            abat=x[:,[1]]

#computing all economic variables
        rel_conso=x[:,[0]]
        #economic variables
        capitalpercapita=kpec*(productivity**(1/(1-alpha)))
        capital=capitalpercapita*L[0:T,:]
        product=production(productivity,capital,L[0:T,:])
        if dim>=2:
            temp=temperature(cumulatedemissions)
            #product=damage_factor(temp,zeros((T,1)))*product
        frac_abat_cost=abatement_cost(abat,theta1[0:T,:])
        abat_cost=frac_abat_cost*product
        consumption=product*rel_conso*(1-frac_abat_cost)
        frac_consumption=rel_conso*(1-frac_abat_cost)
        consumptionpercapita=consumption/L[0:T,:]
        grossinvestment=product*(1-rel_conso)*(1-frac_abat_cost)
        frac_grossinvestment=(1-rel_conso)*(1-frac_abat_cost)
        netinvestment=grossinvestment-delta*capital
        frac_netinvestment=netinvestment/product
        frac_depreciation=delta*capital/product
        growth_rate=(productivity[1:T,:]-productivity[0:T-1,:])**(1/deltaT)-ones((T-1,1))

#saving all economic variables
        savetxt(run_folder+'mean_growth.csv', growth_rate, delimiter=";", header='mean annual growth between periods dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'productivity.csv', productivity, delimiter=";", header='total factor productivity dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'capital.csv', capital, delimiter=";", header='capital dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'capital_per_capita.csv', capitalpercapita, delimiter=";", header='capital per capita dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'consumption.csv', consumption, delimiter=";", header='consumption dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'production.csv', product, delimiter=";", header='production after damages dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'consumption_per_capita.csv', consumptionpercapita, delimiter=";", header='consumption per capita dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'gross_investment.csv', grossinvestment, delimiter=";", header='gross investment dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'net_investment.csv', netinvestment, delimiter=";", header='net investment dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'frac_consumption.csv', frac_consumption, delimiter=";", header='consumption in fraction of current GDP dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'frac_gross_investment.csv', frac_grossinvestment, delimiter=";", header='gross investment in fraction of current GDP dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'frac_net_investment.csv', frac_netinvestment, delimiter=";", header='net investment in fraction of current GDP dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'frac_depreciation.csv', frac_depreciation, delimiter=";", header='depreciation in fraction of current GDP dim '+str(dim)+' stochastic '+str(stochastic))
        #load maximand in f_bellman for direct computation of intertemporal social welfare (caution: this does not work when rra\neq theta)
        
        if tp==0:        
            exec(open(model_folder+'f_bellman_'+str(dim)+'.py').read())
        elif tp==1:
            exec(open(model_folder+'f_bellman_'+str(dim)+'_tp.py').read())
            m=deg+1

        sw_time_horizon=time_horizon[0:T,:].reshape((T))
        sw_x=vstack((x[0:T-1,:],x[[T-2],:]))
        sw_eps=zeros((T,1))

#save climatic variables
        emissionsperperiod=emissions(abat,product,sigm[0:T,:])
        savetxt(run_folder+'emissions.csv', emissionsperperiod, delimiter=";", header='emissions per period dim '+str(dim)+' stochastic '+str(stochastic))
        #special computing when dim 1 (for idealine)
        if dim==1:
            cumulatedemissions=s0+emissionsperperiod.cumsum()
            temp=temperature(cumulatedemissions)

        savetxt(run_folder+'cumulated_emissions.csv', cumulatedemissions, delimiter=";", header='cumulated emissions dim '+str(dim)+' stochastic '+str(stochastic))
        savetxt(run_folder+'temperature.csv', temp, delimiter=";", header='temperature dim '+str(dim)+' stochastic '+str(stochastic))
#computing and saving climate damages and costs-related variables
        if dim>=2:
            savetxt(run_folder+'abatement.csv', abat, delimiter=";", header='abatement dim '+str(dim)+' stochastic '+str(stochastic))
            savetxt(run_folder+'abatement_cost.csv', abat_cost, delimiter=";", header='abatement cost in dollars dim '+str(dim)+' stochastic '+str(stochastic))
            marg_abat_cost=marginal_abatement_cost(abat,theta1[0:T,:],sigm[0:T,:])
            savetxt(run_folder+'marginal_abatement_cost.csv', marg_abat_cost, delimiter=";", header='marginal abatement cost in $/tCO2 dim '+str(dim)+' stochastic '+str(stochastic))
            savetxt(run_folder+'frac_abatement_cost.csv', frac_abat_cost, delimiter=";", header='abatement cost in fraction of current GDP dim '+str(dim)+' stochastic '+str(stochastic))


#computing SCC along the optimal line without crossing the tipping point
            weight_cheb=fill_cheb(deg,dim)
            
            if tp==0:
                coef_V=genfromtxt(run_folder+'/coef_V.csv', dtype=float, delimiter=';')
                
            elif tp==1:
                coef_V_pretipping=genfromtxt(run_folder+'/coef_V_pretipping.csv', dtype=float, delimiter=';')
                coef_V_uncertain=genfromtxt(run_folder+'/coef_V_uncertain.csv', dtype=float, delimiter=';')
                coef_V_posttipping=genfromtxt(run_folder+'/coef_V_posttipping.csv', dtype=float, delimiter=';')
                #info=genfromtxt(run_folder+'/info_tp.csv', dtype=float, delimiter=';')

            # How small must dS be?
            dS=u[:,[1]]*0.001
            
            if dim==2:
                Delta_u=hstack((zeros((T,1)),dS))
            if dim==3:
                Delta_u=hstack((zeros((T,1)),dS,zeros((T,1))))
            uplus=u+Delta_u
            
            if tp==0:
                
                u_max=genfromtxt(run_folder+'/u_max.csv', dtype=float, delimiter=';')
                u_min=genfromtxt(run_folder+'/u_min.csv', dtype=float, delimiter=';')
                u_adj=minimum(maximum(u,u_min),u_max)
                uplus_adj=minimum(maximum(uplus,u_min),u_max)
                u_min_adj=adjust_min(u_min,u_max,deg+1)
                u_max_adj=adjust_max(u_min,u_max,deg+1)

                # For all preferences
                SCC=-bet*(approx_Chebyshev(coef_V[1:T,:],weight_cheb,uplus_adj[1:T,:],u_min_adj[1:T,:],u_max_adj[1:T,:])-approx_Chebyshev(coef_V[1:T,:],weight_cheb,u_adj[1:T,:],u_min_adj[1:T,:],u_max_adj[1:T,:]))/(dS[1:T,:]*partialc_maximand(sw_time_horizon,sw_x,u,sw_eps)[:T-1,:])

                savetxt(run_folder+'social_cost_carbon.csv', SCC, delimiter=";", header='social cost of carbon in $/tCO2 dim '+str(dim)+' stochastic '+str(stochastic))
                              
            elif tp==1:
                u_max_pretipping=genfromtxt(run_folder+'/u_max_pretipping.csv', dtype=float, delimiter=';')
                u_min_pretipping=genfromtxt(run_folder+'/u_min_pretipping.csv', dtype=float, delimiter=';')
                u_adj=minimum(maximum(u,u_min_pretipping),u_max_pretipping)
                uplus_adj=minimum(maximum(uplus,u_min_pretipping),u_max_pretipping)
                u_min_adj=adjust_min(u_min_pretipping,u_max_pretipping,deg+1)
                u_max_adj=adjust_max(u_min_pretipping,u_max_pretipping,deg+1)
                u_max_uncertain=genfromtxt(run_folder+'/u_max_uncertain.csv', dtype=float, delimiter=';')
                u_min_uncertain=genfromtxt(run_folder+'/u_min_uncertain.csv', dtype=float, delimiter=';')
                u_adj3=minimum(maximum(u,u_min_uncertain),u_max_uncertain)
                uplus_adj3=minimum(maximum(uplus,u_min_uncertain),u_max_uncertain)
                u_min_adj3=adjust_min(u_min_uncertain,u_max_uncertain,deg+1)
                u_max_adj3=adjust_max(u_min_uncertain,u_max_uncertain,deg+1)
                u_max_posttipping=genfromtxt(run_folder+'/u_max_posttipping.csv', dtype=float, delimiter=';')
                u_min_posttipping=genfromtxt(run_folder+'/u_min_posttipping.csv', dtype=float, delimiter=';')
                u_adj2=minimum(maximum(u,u_min_posttipping),u_max_posttipping)
                uplus_adj2=minimum(maximum(uplus,u_min_posttipping),u_max_posttipping)
                u_min_adj2=adjust_min(u_min_posttipping,u_max_posttipping,deg+1)
                u_max_adj2=adjust_max(u_min_posttipping,u_max_posttipping,deg+1)
                #h is the probability of tipping
                S_max=(threshold_max/beta_temp)*ones((T-1,1))
                S_min=(threshold_min/beta_temp)*ones((T-1,1))
                hproba=h(u[0:T-1,[1]],u[1:T,[1]],S_min,S_max)
                hplus=h(u[0:T-1,[1]],u[1:T,[1]]+Delta_u[1:T,[1]],S_min,S_max)
                hplus2=h(uplus[0:T-1,[1]],u[1:T,[1]]+Delta_u[0:T-1,[1]],S_min,S_max)

                #increased probability if A_{t+1} increases
                dh=(hplus-hproba)/dS[1:T,:]

                #increased probability if both A_t and A_{t+1} increase
                dh2=(h(uplus[0:T-1,[1]],u[1:T,[1]]+Delta_u[0:T-1,[1]],S_min,S_max)-h(u[0:T-1,[1]],u[1:T,[1]],S_min,S_max))/dS[1:T,:]
                
                #We define indicators to know which bellman function applies.
                def indic_pre(S):
                    return maximum(0,sign(S_min-S))
                
                def indic_unc(S):
                    return maximum(0,sign(sign(S-S_min)+0.5))*maximum(0,sign(S_max-S))

                def indic_post(S):
                    return maximum(0,sign(S-S_max))

                #we will need partialCu
                partialCu=(1-indic_post(u[:T-1,[1]]))*partialc_maximand_pretipping(sw_time_horizon,sw_x,u,sw_eps)[:T-1,:]+indic_post(u[:T-1,[1]])*partialc_maximand_posttipping(sw_time_horizon,sw_x,u,sw_eps)[:T-1,:]
                partialCun=(1-indic_post(u[1:T,[1]]))*partialc_maximand_pretipping(sw_time_horizon,sw_x,u,sw_eps)[1:T,:]+indic_post(u[1:T,[1]])*partialc_maximand_posttipping(sw_time_horizon,sw_x,u,sw_eps)[1:T,:]                

                #Direct method to calculate SCC / EU preferences
                if epsilon_b==0:
                    E_V=indic_pre(u[1:T,[1]])*approx_Chebyshev(coef_V_pretipping[1:T,:],weight_cheb,u_adj[1:T,:],u_min_adj[1:T,:],u_max_adj[1:T,:])+(1-indic_pre(u[1:T,[1]]))*((1-hproba)*approx_Chebyshev(coef_V_uncertain[1:T,:],weight_cheb,u_adj3[1:T,:],u_min_adj3[1:T,:],u_max_adj3[1:T,:])+hproba*approx_Chebyshev(coef_V_posttipping[1:T,:],weight_cheb,u_adj2[1:T,:],u_min_adj2[1:T,:],u_max_adj2[1:T,:]))
                    E_V_plus=indic_pre(uplus[1:T,[1]])*approx_Chebyshev(coef_V_pretipping[1:T,:],weight_cheb,uplus_adj[1:T,:],u_min_adj[1:T,:],u_max_adj[1:T,:])+(1-indic_pre(uplus[1:T,[1]]))*((1-hplus)*approx_Chebyshev(coef_V_uncertain[1:T,:],weight_cheb,uplus_adj3[1:T,:],u_min_adj3[1:T,:],u_max_adj3[1:T,:])+hplus*approx_Chebyshev(coef_V_posttipping[1:T,:],weight_cheb,uplus_adj2[1:T,:],u_min_adj2[1:T,:],u_max_adj2[1:T,:]))
                    d_E_V1=(E_V_plus-E_V)/dS[1:T,:]
                    SCC=-bet*d_E_V1/partialCu

                #risk-sensitive preferences
                else:
                    E_V_Bommier=indic_pre(u[1:T,[1]])*approx_Chebyshev(coef_V_pretipping[1:T,:],weight_cheb,u_adj[1:T,:],u_min_adj[1:T,:],u_max_adj[1:T,:])+(1-indic_pre(u[1:T,[1]]))*Bommier_function2(approx_Chebyshev(coef_V_uncertain[1:T,:],weight_cheb,u_adj3[1:T,:],u_min_adj3[1:T,:],u_max_adj3[1:T,:]),approx_Chebyshev(coef_V_posttipping[1:T,:],weight_cheb,u_adj2[1:T,:],u_min_adj2[1:T,:],u_max_adj2[1:T,:]),hproba)
                    E_V_plus_Bommier=indic_pre(uplus[1:T,[1]])*approx_Chebyshev(coef_V_pretipping[1:T,:],weight_cheb,uplus_adj[1:T,:],u_min_adj[1:T,:],u_max_adj[1:T,:])+(1-indic_pre(uplus[1:T,[1]]))*Bommier_function2(approx_Chebyshev(coef_V_uncertain[1:T,:],weight_cheb,uplus_adj3[1:T,:],u_min_adj3[1:T,:],u_max_adj3[1:T,:]),approx_Chebyshev(coef_V_posttipping[1:T,:],weight_cheb,uplus_adj2[1:T,:],u_min_adj2[1:T,:],u_max_adj2[1:T,:]),hplus)
                    d_E_V1_Bommier=(E_V_plus_Bommier-E_V_Bommier)/dS[1:T,:]
                    SCC=-bet*d_E_V1_Bommier/partialCu

                print('SCC run'+str(run_no))
                print(SCC[0])
                savetxt(run_folder+'social_cost_carbon.csv', SCC, delimiter=";", header='social cost of carbon in $/tCO2 dim '+str(dim)+' stochastic '+str(stochastic))

                #Calculation for stochastic draws.
                #define nb of draws
                nbofdraws=10

                #define run where we want the draws
                LL=[4,7,11,13,60,63,65,67,69,74,77,79,81,83]

                if run_no in LL:
                        if not threshold_max==threshold_min:
                                exec(open(preprod_folder+'stochasticdraws3.py').read())

                #exec(open(preprod_folder+'stochasticdraws.py').read())                                
                #average_SCC=[-1]
                #probability=[-1]

# Create a matrix that contains all info we want
listofrho1y=[0.015]
listoftheta=[0.8]
listofgamma=[0.5,1.5,2.,2.5,3.,3.5,5.,7.,10.,14.,17.,20.]
listofpi2=[0.0028,0.004,0.005,0.006,0.008,0.010,0.020,0.050,0.1]
SCC_to_analyze=to_analyze
MatrixofSCC=np.zeros((size(SCC_to_analyze),12))
MatrixofDecomp = np.zeros((size(SCC_to_analyze), 11))
compteur=0
nbofbroken=0
Notrun=[]

for run_no in SCC_to_analyze:
#for run_no in array([5083,5084]):
    #getting the run outputs folder
    run_outputs=runs_outputs[run_id==run_no][0]
    
    if not run_outputs:
        print("Cannot analyze run number "+str(run_no)+", it has not been run")
        Notrun.append(run_no)
    else:
        run_folder=outputs_folder+run_outputs+pathsep
        
        broken=0
        SCC0=0
        finalT=0

        exec(open(model_folder+'read_parameters.py').read())
        
        if os.path.isfile(run_folder+'simulate_break.csv'):
            #print "Run number "+str(run_no)+" has broken"
            broken=1
            nbofbroken=nbofbroken+1
            
                
        if pi2!=0:
            
            ValueOfSCC=genfromtxt(run_folder+'/social_cost_carbon.csv', dtype=float, delimiter=';')
            SCC0=ValueOfSCC[0]
            temperaturepath=genfromtxt(run_folder+'/temperature.csv', dtype=float, delimiter=';')
            finalT=temperaturepath[T-1]
        MatrixofSCC[compteur,]=[run_no,rho1y,theta,epsilon_b,J,broken,tp,deterministic,SCC0,finalT,threshold_max,pi2]

    compteur=compteur+1

dataf=pd.DataFrame(MatrixofSCC)
dataf2=pd.DataFrame(MatrixofDecomp)
dataf.to_excel(preprod_folder+'SCC_matrix_'+study+'.xlsx')
dataf2.to_excel(preprod_folder+'Decomp_matrix_'+study+'.xlsx')



        
