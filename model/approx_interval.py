# -*- coding: utf-8 -*-

#defining approximation interval 
if dim==1:
    strip_width=0.3

    #stationary capital per efficient capita from ramsey formula
    ustat=(alpha*deltaT/(delta*ones((T+2,1))+rho*ones((T+2,1))+theta*deltaT*(GA/(1-alpha))))**(1/(1-alpha))

    #trend of initial capital
    u0trend=K0/((A0**(1/(1-alpha)))*POP0)*ones((T+2,1))

    #mixing trends
    mix_param_min=0.95
    ucentered_min=ustat*(ones((T+2,1))-(mix_param_min**(deltaT*time_horizon)))+u0trend*(mix_param_min**(deltaT*time_horizon))
    mix_param_max=0.95
    ucentered_max=ustat*(ones((T+2,1))-(mix_param_max**(deltaT*time_horizon)))+u0trend*(mix_param_max**(deltaT*time_horizon))

    #boundary of interpolation
    u_min=(1-strip_width)*ucentered_min[0:T,:]
    u_max=(1+strip_width)*ucentered_max[0:T,:]


if dim==2:
    u_m=genfromtxt(prev_run_folder+'state_V.csv', dtype=float, delimiter=';')
    x_prev_run=genfromtxt(prev_run_folder+'control.csv', dtype=float, delimiter=';')
    u_ideal=genfromtxt(ideal_run_folder+'state_V.csv', dtype=float, delimiter=';')
    k_ideal=k_ramsey=(A[0:T,:]**(1/(1-alpha)))*L[0:T,:]*u_ideal.reshape((T,1))

    if shape(u_m)==(T,):
        k_ramsey=k_ramsey=(A[0:T,:]**(1/(1-alpha)))*L[0:T,:]*u_m.reshape((T,1))
        k_min=0.5*k_ramsey
        kpec_min=0.5*u_m.reshape((T,1))
        k_max=1.2*k_ramsey
        kpec_max=1.2*u_m.reshape((T,1))

    else:
        k_ramsey=(A[0:T,:]**(1/(1-alpha)))*L[0:T,:]*u_m[:,0].reshape((T,1))
        k_min=0.5*k_ramsey
        kpec_min=0.5*u_m[:,0].reshape((T,1))
        k_max=1.2*k_ramsey
        kpec_max=1.2*u_m[:,0].reshape((T,1))

    emissions_max=zeros((T,1))
    emissions_min=zeros((T,1))
    emissions_ideal=zeros((T,1))
    ab_pess=zeros((T,1))

    for j in list(range(1,T+1)):
        emissions_max[j-1,:]=emissions(ab_pess[j-1,:],production(A[j-1,:],k_max[j-1,:],L[[j-1],:]),sigm[[j-1],:])
        emissions_min[j-1,:]=emissions(0.8,production(A[j-1,:],k_min[j-1,:],L[[j-1],:]),sigm[[j-1],:])
        emissions_ideal[j-1,:]=emissions(0,production(A[j-1,:],k_ideal[j-1,:],L[[j-1],:]),sigm[[j-1],:])

    s_max=1.1*s0*ones((T,1))
    s_min=0.9*s0*ones((T,1))
    s_max_ideal=1.1*s0*ones((T,1))

    for j in list(range(1,T)):
        s_max[j,:]=s_max[j-1,:]+1.1*emissions_max[j-1,:]
        s_max_ideal[j,:]=s_max_ideal[j-1,:]+1.1*emissions_ideal[j-1,:]

    if shape(u_m)==(T,2):
       s_max=1.1*u_m[:,1].reshape((T,1))

       if x_prev_run[0,1]<0.01:
           s_min=0.9*s0*ones((T,1))
           if pi3==6:
               for j in list(range(1,T)):
                   s_max[j,:]=minimum(s_max[j,:],6/beta_temp*ones(shape(s_max[j,:])))
       else:
           s_min=0.9*s0*ones((T,1))

    for j in list(range(1,T)):
        if j<int(50/deltaT)+2:
            s_max[j,:]=maximum(s_max[j,:],s_max_ideal[j,:])
        else:
            s_max[j,:]=maximum(s_max[j,:],s_max_ideal[int(50/deltaT)+2,:])

    u_min=concatenate((kpec_min,s_min),axis=1)
    u_max=concatenate((kpec_max,s_max),axis=1)

            
if dim==3:
    x_det=genfromtxt(prev_run_folder+'control.csv', dtype=float, delimiter=';')
    u_det=genfromtxt(prev_run_folder+'state_V.csv', dtype=float, delimiter=';')
    u_ideal=genfromtxt(ideal_run_folder+'state_V.csv', dtype=float, delimiter=';')
    k_ideal=k_ramsey=(A[0:T,:]**(1/(1-alpha)))*L[0:T,:]*u_ideal.reshape((T,1))
    temp_exo=temperature(u_det[:,[1]])
    k_min=0.7*u_det[:,0].reshape((T,1))
    k_max=1.2*u_det[:,0].reshape((T,1))
    s_max=1.01*u_det[:,1].reshape((T,1))

    if x_det[0,1]<0.01:
        s_min=0.99*s0*ones((T,1))
        if kappa3==6:
            for j in list(range(1,T)):
                s_max[j,:]=minimum(s_max[j,:],6/beta_temp*ones(shape(s_max[j,:])))
    else:
        s_min=0.4*u_det[:,1].reshape((T,1))
    emissions_ideal=zeros((T,1))

    for j in list(range(1,T+1)):
        emissions_ideal[j-1,:]=emissions(0,production(A[j-1,:],k_ideal[j-1,:],L[[j-1],:]),sigm[[j-1],:])
    s_max_ideal=1.01*s0*ones((T,1))

    for j in list(range(1,T)):
        s_max_ideal[j,:]=s_max_ideal[j-1,:]+emissions_ideal[j-1,:]

    for j in list(range(1,T)):
        if j<int(50/deltaT)+2:
            s_max[j,:]=maximum(s_max[j,:],s_max_ideal[j,:])
        else:
            s_max[j,:]=maximum(s_max[j,:],s_max_ideal[int(50/deltaT)+2,:])

    a_max=1.001*A0*ones((T,1))
    a_min=0.999*A0*ones((T,1))

    for j in list(range(1,T)):
        a_min[j,:]=a_min[j-1,:]*(1+growth(j,a_min[j-1,:],temperature(s_max[j-1,:]),0))
        a_max[j,:]=a_max[j-1,:]*(1+growth(j,a_max[j-1,:],temperature(s_min[j-1,:]),0))

    u_max=concatenate((k_max,s_max,a_max),axis=1)
    u_min=concatenate((k_min,s_min,a_min),axis=1)

