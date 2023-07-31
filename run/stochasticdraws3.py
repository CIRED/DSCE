# -*- coding: utf-8 -*-

#EMISSIONSSAVE=emissionsperperiod
import pandas as pd

guess0=array((0.5,0.1)).reshape(1,2)

u0=concatenate((K0/((A0**(1/(1-alpha)))*POP0)*ones((1,1)),s0*ones((1,1))),axis=1)
expectedproba=0.
averagetime=0
MatrixofDecomp = np.zeros((20, 7))

average_SCC=zeros((T-1,1))

TIME=deltaT*time_horizon[0:T,:]


#plt.plot(TIME[:50],SCC_firstmethod[:50],label='without crossing')
#plt.legend()

#plt.figure()
for wuwu in list(range(nbofdraws)):
    print('draw nb '+str(wuwu))
    # if the threshold is certain, strong discontinuity : use 3 functions to make sure we're not missing the optimal point
    #exec(open(model_folder+'simulate_tp_3interpolations_draws.py').read())
    exec(open(model_folder+'simulate_tp.py').read())

    if info_tp[0]==1:
        print(info_tp[1])
        averagetime=averagetime+info_tp[1]
        expectedproba=expectedproba+info_tp[0]/float(nbofdraws)
        print(averagetime)
    tp_yes = np.array([info_tp[0]])

    if wuwu==0:
        tipping=np.array([info_tp[0]])
    else:
        tipping=concatenate((tipping,tp_yes),axis=0)
    kpec=u[:,[0]]                  
    cumulatedemissions=u[:,[1]]
    productivity=A[0:T,:]
    abat=x[:,[1]]
    #computing all economic variables
    rel_conso=x[:,[0]]
    #economic variables
    capitalpercapita=kpec*(productivity**(1/(1-alpha)))

    if wuwu==0:
        capitalperc=capitalpercapita
    else:
        capitalperc=concatenate((capitalperc,capitalpercapita),axis=1)
    capital=capitalpercapita*L[0:T,:]
    
    temp=temperature(cumulatedemissions)
    #product=production(productivity,capital,L[0:T,:])
    #We take net produc only for the case when no tipping ; if the analysis is made (traj for instance), drop runs for which tipping point crossed.
    product=damage_factor_pretipping(temp,eps)*production(productivity,capital,L[0:T,:])

    consumption=product*rel_conso*(1-frac_abat_cost)
    consumptionpercapita=consumption/L[0:T,:]

    if wuwu==0:
        consumptionpercapitat=consumptionpercapita
    else:
        consumptionpercapitat=concatenate((consumptionpercapitat,consumptionpercapita),axis=1)
        
    if wuwu==0:
        tempt=temp
    else:
        tempt=concatenate((tempt,temp),axis=1)

    if info_tp[0]==1:
        product=damage_factor_posttipping(temp,zeros((T,1)))*product
    else:
        product=damage_factor_pretipping(temp,zeros((T,1)))*product
        
    if wuwu==0:
        productt=product
    else:
        productt=concatenate((productt,product),axis=1)
        
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
   
   
    exec(open(model_folder+'f_bellman_'+str(dim)+'_tp.py').read())
    m=deg+1

    sw_time_horizon=time_horizon[0:T,:].reshape((T))
    sw_x=vstack((x[0:T-1,:],x[[T-2],:]))
    sw_eps=zeros((T,1))

        
#SCC
    weight_cheb=fill_cheb(deg,dim)                                          
    emissionsperperiod=emissions(abat,product,sigm[0:T,:])
    dS=u[:,[1]]*0.001
    Delta_u=hstack((zeros((T,1)),dS))
    uplus=u+Delta_u        
    S_max=(threshold_max/beta_temp)*ones((T-1,1))
    S_min=(threshold_min/beta_temp)*ones((T-1,1))
    hproba=h(u[0:T-1,[1]],u[1:T,[1]],S_min,S_max)                
    hplus=h(u[0:T-1,[1]],uplus[1:T,[1]],S_min,S_max)
    
# for CRRA  preferences
    if epsilon_b==0:
        E_V=indic_pre(u[1:T,[1]])*approx_Chebyshev(coef_V_pretipping[1:T,:],weight_cheb,u_adj[1:T,:],u_min_adj[1:T,:],u_max_adj[1:T,:])+(1-indic_pre(u[1:T,[1]]))*((1-hproba)*approx_Chebyshev(coef_V_uncertain[1:T,:],weight_cheb,u_adj3[1:T,:],u_min_adj3[1:T,:],u_max_adj3[1:T,:])+hproba*approx_Chebyshev(coef_V_posttipping[1:T,:],weight_cheb,u_adj2[1:T,:],u_min_adj2[1:T,:],u_max_adj2[1:T,:]))
        E_V_plus=indic_pre(uplus[1:T,[1]])*approx_Chebyshev(coef_V_pretipping[1:T,:],weight_cheb,uplus_adj[1:T,:],u_min_adj[1:T,:],u_max_adj[1:T,:])+(1-indic_pre(uplus[1:T,[1]]))*((1-hplus)*approx_Chebyshev(coef_V_uncertain[1:T,:],weight_cheb,uplus_adj3[1:T,:],u_min_adj3[1:T,:],u_max_adj3[1:T,:])+hplus*approx_Chebyshev(coef_V_posttipping[1:T,:],weight_cheb,uplus_adj2[1:T,:],u_min_adj2[1:T,:],u_max_adj2[1:T,:]))
        d_E_V1=(E_V_plus-E_V)/dS[1:T,:]
        SCC_firstmethod=-bet*d_E_V1/partialc_maximand_pretipping(sw_time_horizon,sw_x,u,sw_eps)[:T-1,:]

#for general risk-sensitive preferences (Hansen, Sargent)
    else:
        E_V_Bommier=indic_pre(u[1:T,[1]])*(-(1/epsilon_b)*log(exp(-epsilon_b*approx_Chebyshev(coef_V_pretipping[1:T,:],weight_cheb,u_adj[1:T,:],u_min_adj[1:T,:],u_max_adj[1:T,:]))))+(1-indic_pre(u[1:T,[1]]))*Bommier_function2(approx_Chebyshev(coef_V_uncertain[1:T,:],weight_cheb,u_adj3[1:T,:],u_min_adj3[1:T,:],u_max_adj3[1:T,:]),approx_Chebyshev(coef_V_posttipping[1:T,:],weight_cheb,u_adj2[1:T,:],u_min_adj2[1:T,:],u_max_adj2[1:T,:]),hproba)
        E_V_plus_Bommier=indic_pre(uplus[1:T,[1]])*(-(1/epsilon_b)*log(exp(-epsilon_b*approx_Chebyshev(coef_V_pretipping[1:T,:],weight_cheb,uplus_adj[1:T,:],u_min_adj[1:T,:],u_max_adj[1:T,:]))))+(1-indic_pre(uplus[1:T,[1]]))*Bommier_function2(approx_Chebyshev(coef_V_uncertain[1:T,:],weight_cheb,uplus_adj3[1:T,:],u_min_adj3[1:T,:],u_max_adj3[1:T,:]),approx_Chebyshev(coef_V_posttipping[1:T,:],weight_cheb,uplus_adj2[1:T,:],u_min_adj2[1:T,:],u_max_adj2[1:T,:]),hplus)
        d_E_V1=(E_V_plus_Bommier-E_V_Bommier)/dS[1:T,:]
        SCC_firstmethod=-bet*d_E_V1/partialc_maximand_pretipping(sw_time_horizon,sw_x,u,sw_eps)[:T-1,:]

#for risk-sensitive preferences
    if info_tp[0]==1:
        dV_post=(-(1/epsilon_b)*log(exp(-epsilon_b*(approx_Chebyshev(coef_V_posttipping[1:T,:],weight_cheb,uplus_adj2[1:T,:],u_min_adj2[1:T,:],u_max_adj2[1:T,:]))))-(-(1/epsilon_b)*log(exp(-epsilon_b*(approx_Chebyshev(coef_V_posttipping[1:T,:],weight_cheb,u_adj2[1:T,:],u_min_adj2[1:T,:],u_max_adj2[1:T,:]))))))/(dS[1:T,:])
        SCC_aftertipping=-bet*dV_post/partialc_maximand_posttipping(sw_time_horizon,sw_x,u,sw_eps)[:T-1,:]
        SCC_firstmethod[info_tp[1]:]=SCC_aftertipping[info_tp[1]:]  
    
    #average_SCC=average_SCC+SCC_firstmethod/nbofdraws
    if wuwu==0:
        SCC_save=SCC_firstmethod
    else:
        SCC_save=concatenate((SCC_save,SCC_firstmethod),axis=1)
    if wuwu==0:
        emissions_save=u[:,[1]]
    else:
        emissions_save=concatenate((emissions_save,u[:,[1]]),axis=1)

for wuwu in list(range(nbofdraws)):
#Value function decomposition
#From Lemoine and Traeger (2016)
    dS = u[:, [1]] * 0.001
    dmu = x[:, [1]] * 0.001
    Delta_u = hstack((zeros((T, 1)), dS))
    uplus = u + Delta_u
    S_max = (threshold_max / beta_temp) * ones((T - 1, 1))
    S_min = (threshold_min / beta_temp) * ones((T - 1, 1))
    hproba = h(u[0:T - 1, [1]], u[1:T, [1]], S_min, S_max)
    hplus = h(u[0:T - 1, [1]], uplus[1:T, [1]], S_min, S_max)

 ###########################
    #Calculation of MHE neutral
    ###########################
    MHE1 = (hplus - hproba) / dS[1:T, :]
    MHE2 = (uplus[1:T, [1]] - u[1:T, [1]]) / dmu[0:T - 1, :]

    # If additive preferences
    VPRE = approx_Chebyshev(coef_V_pretipping[1:T, :], weight_cheb, u_adj[1:T, :], u_min_adj[1:T, :],u_max_adj[1:T, :])
    VPOST = approx_Chebyshev(coef_V_posttipping[1:T, :], weight_cheb, u_adj2[1:T, :], u_min_adj2[1:T, :],u_max_adj2[1:T, :])

    MHE3 = VPRE - VPOST
    MHE_neutral = - MHE1 * MHE2 * MHE3

    #Calculation of DWI neutral
    DWI1 = hproba
    DWI3 = MHE2

    # Certain pre-tipping VF, uncertain boundaries
    u_max_uncertain = genfromtxt(run_folder + '/u_max_uncertain.csv', dtype=float, delimiter=';')
    u_min_uncertain = genfromtxt(run_folder + '/u_min_uncertain.csv', dtype=float, delimiter=';')
    u_adj = minimum(maximum(u, u_min_uncertain), u_max_uncertain)
    uplus_adj = minimum(maximum(uplus, u_min_uncertain), u_max_uncertain)
    u_min_adj = adjust_min(u_min_uncertain, u_max_uncertain, deg + 1)
    u_max_adj = adjust_max(u_min_uncertain, u_max_uncertain, deg + 1)

    u_max_posttipping = genfromtxt(run_folder + '/u_max_posttipping.csv', dtype=float, delimiter=';')
    u_min_posttipping = genfromtxt(run_folder + '/u_min_posttipping.csv', dtype=float, delimiter=';')
    u_adj2 = minimum(maximum(u, u_min_posttipping), u_max_posttipping)
    uplus_adj2 = minimum(maximum(uplus, u_min_posttipping), u_max_posttipping)
    u_min_adj2 = adjust_min(u_min_posttipping, u_max_posttipping, deg + 1)
    u_max_adj2 = adjust_max(u_min_posttipping, u_max_posttipping, deg + 1)

    E_V_pre = approx_Chebyshev(coef_V_pretipping[1:T, :], weight_cheb, u_adj[1:T, :], u_min_adj[1:T, :],u_max_adj[1:T, :])
    E_V_plus_pre = approx_Chebyshev(coef_V_pretipping[1:T, :], weight_cheb, uplus_adj[1:T, :], u_min_adj[1:T, :],u_max_adj[1:T, :])
    E_V_post = approx_Chebyshev(coef_V_posttipping[1:T, :], weight_cheb, u_adj2[1:T, :], u_min_adj2[1:T, :],u_max_adj2[1:T, :])
    E_V_plus_post = approx_Chebyshev(coef_V_posttipping[1:T, :], weight_cheb, uplus_adj2[1:T, :],u_min_adj2[1:T, :], u_max_adj2[1:T, :])

    DWI2a = (E_V_pre - E_V_pre) / dS[1:T, :]
    DWI2b = (E_V_plus_post - E_V_post) / dS[1:T, :]

    DWI2 = DWI2a - DWI2b
    DWI_neutral = - DWI1 * DWI2 * DWI3

    #####Let us add temporal risk aversion
    A = 1/((1-hproba) * exp(-epsilon_b * VPRE) + hproba * exp(-epsilon_b * VPOST))

    #Calculation of MHE risk-averse
    MHE4 = exp(-epsilon_b * VPRE) - exp(- epsilon_b * VPOST)
    MHE_sensitive = - A / epsilon_b * MHE1 * MHE2 * MHE4

    #Calculation of DWI risk-averse
    DWI_sensitive = A * DWI1 * MHE2 * (DWI2a * exp(- epsilon_b * VPRE) - DWI2b * (- epsilon_b * VPOST))

    #Calculation of MPre risk-averse
    MPRE_sensitive = A * DWI2a * MHE2 * exp( - epsilon_b * VPRE)

    #Calculation of the complete effect
 #   complete1 = (uplus[1:T, [1]] - u[1:T, [1]])/dS[0:T-1, :]
#    complete1 = (uplus[2:T, [1]] - u[2:T, [1]])/dS[1:T-1, :]
    complete2 = exp(- epsilon_b * VPRE) * A

    if wuwu==0:
        DWI_neutral_save=DWI_neutral
        MHE_neutral_save=MHE_neutral
        MHE_sensitive_save=MHE_sensitive
        DWI_sensitive_save=DWI_sensitive
        MPRE_sensitive_save=MPRE_sensitive
  #      complete1_save=complete1
        complete2_save=complete2
        MHE_immediate_save=MHE_sensitive[0]
        DWI_immediate_save=DWI_sensitive[0]
        MHE_immediate_save_neutral=MHE_neutral[0]
        DWI_immediate_save_neutral=DWI_neutral[0]

        hproba_save=hproba
    else:
        DWI_neutral_save=concatenate((DWI_neutral_save,DWI_neutral),axis=1)
        MHE_neutral_save=concatenate((MHE_neutral_save,MHE_neutral),axis=1)
        MHE_sensitive_save=concatenate((MHE_sensitive_save,MHE_sensitive),axis=1)
        DWI_sensitive_save=concatenate((DWI_sensitive_save,DWI_sensitive),axis=1)
        MPRE_sensitive_save=concatenate((MPRE_sensitive_save,MPRE_sensitive),axis=1)
   #     complete1_save=concatenate((complete1_save,complete1),axis=1)
        complete2_save=concatenate((complete2_save,complete2),axis=1)
        MHE_immediate_save=concatenate((MHE_immediate_save, MHE_sensitive[0]))
        DWI_immediate_save=concatenate((DWI_immediate_save, DWI_sensitive[0]))
        MHE_immediate_save_neutral=concatenate((MHE_immediate_save_neutral, MHE_neutral[0]))
        DWI_immediate_save_neutral=concatenate((DWI_immediate_save_neutral, DWI_neutral[0]))

        hproba_save=concatenate((hproba_save, hproba))

    # print average_SCC
    averagetime = 0
    # averagetime=averagetime/(nbofdraws*expectedproba)
    probaoftipping = array([[expectedproba], [averagetime]])
    print(probaoftipping)

    savetxt(run_folder + 'tipping.csv', tipping, delimiter=";",header='probability of tipping + average period of tipping ')
    savetxt(run_folder + 'proba.csv', probaoftipping, delimiter=";",header='probability of tipping + average period of tipping ')
    savetxt(run_folder + 'hproba_tot.csv', hproba_save, delimiter=";", header='expected probability of tipping')
    savetxt(run_folder + 'productt.csv', productt, delimiter=";",header='probability of tipping + average period of tipping ')
    savetxt(run_folder + 'cumem.csv', tempt, delimiter=";",header='probability of tipping + average period of tipping ')
    savetxt(run_folder + 'tempt.csv', tempt, delimiter=";",header='probability of tipping + average period of tipping ')
    savetxt(run_folder + 'conspercapita.csv', consumptionpercapitat, delimiter=";",header='probability of tipping + average period of tipping')
    savetxt(run_folder + 'traj_SCC.csv', SCC_save, delimiter=";",header='average value of SCC based on ' + str(nbofdraws) + ' runs')
    savetxt(run_folder + 'emissions_save.csv', emissions_save, delimiter=";",header='average value of SCC based on ' + str(nbofdraws) + ' runs')

    # decomposition
    savetxt(run_folder + 'dwineutral.csv', DWI_neutral_save, delimiter=";",header='average value of SCC based on ' + str(nbofdraws) + ' runs')
    savetxt(run_folder + 'mheneutral.csv', MHE_neutral_save, delimiter=";",header='average value of SCC based on ' + str(nbofdraws) + ' runs')
    savetxt(run_folder + 'dwi.csv', DWI_sensitive_save, delimiter=";",header='average value of SCC based on ' + str(nbofdraws) + ' runs')
    savetxt(run_folder + 'mhe.csv', MHE_sensitive_save, delimiter=";",header='average value of SCC based on ' + str(nbofdraws) + ' runs')
    savetxt(run_folder + 'mhe_immediate.csv', MHE_immediate_save, delimiter=";",header='average value of SCC based on ' + str(nbofdraws) + ' runs')
    savetxt(run_folder + 'dwi_immediate.csv', DWI_immediate_save, delimiter=";",header='average value of SCC based on ' + str(nbofdraws) + ' runs')
    savetxt(run_folder + 'complete2.csv', complete2_save, delimiter=";", header='average value of SCC based on ' + str(nbofdraws) + ' runs')

MHE_immediate_exp=mean(MHE_immediate_save)
DWI_immediate_exp=mean(DWI_immediate_save)
MHE_immediate_exp_neutral=mean(MHE_immediate_save_neutral)
DWI_immediate_exp_neutral=mean(DWI_immediate_save_neutral)

time = list(range(1, T))
beta=[bet]*(T-1)

for t in range(0,T-1):
    beta[t]=beta[t]**time[t]

complete_MHE=np.empty((T,nbofdraws))
complete_MHE[:]=np.NaN
complete_DWI=np.empty((T,nbofdraws))
complete_DWI[:]=np.NaN

for t in list(range(0,nbofdraws)):
    for i in list(range(1, T-1)):
        complete_MHE[i-1][t] = beta[i-1]
        complete_MHE[i-1][t] = complete_MHE[i-1][t] * complete2_save[i-1][t]
        complete_MHE[i-1][t] = complete_MHE[i-1][t] * MHE_sensitive_save[i-1][t]
        complete_DWI[i-1][t] = beta[i-1] * complete2_save[i-1][t]*DWI_sensitive_save[i-1][t]

complete_MHE_rs=nansum(complete_MHE, axis=0)
complete_DWI_rs=nansum(complete_DWI, axis=0)
complete_MHE_rs=nanmean(complete_MHE_rs, axis=0)
complete_DWI_rs=nanmean(complete_DWI_rs, axis=0)
complete_MHE_rs=[complete_MHE_rs]
complete_DWI_rs=[complete_DWI_rs]
MHE_immediate_exp_neutral=[MHE_immediate_exp_neutral]
DWI_immediate_exp_neutral=[DWI_immediate_exp_neutral]
DWI_immediate_exp=[DWI_immediate_exp]
MHE_immediate_exp=[MHE_immediate_exp]
savetxt(run_folder + 'complete_MHE_rs.csv', complete_MHE_rs, delimiter=";", header='average value of SCC based on ' + str(nbofdraws) + ' runs')
savetxt(run_folder + 'complete_DWI_rs.csv', complete_DWI_rs, delimiter=";", header='average value of SCC based on ' + str(nbofdraws) + ' runs')
savetxt(run_folder + 'MHE_immediate_exp_neutral.csv', MHE_immediate_exp_neutral, delimiter=";", header='average value of SCC based on ' + str(nbofdraws) + ' runs')
savetxt(run_folder + 'DWI_immediate_exp_neutral.csv', DWI_immediate_exp_neutral, delimiter=";", header='average value of SCC based on ' + str(nbofdraws) + ' runs')
savetxt(run_folder + 'DWI_immediate_exp.csv', DWI_immediate_exp, delimiter=";", header='average value of SCC based on ' + str(nbofdraws) + ' runs')
savetxt(run_folder + 'MHE_immediate_exp.csv', MHE_immediate_exp, delimiter=";", header='average value of SCC based on ' + str(nbofdraws) + ' runs')


