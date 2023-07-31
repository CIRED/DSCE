# -*- coding: utf-8 -*-

# Simulation of optimal path
# forward iteration using the approximation of the value functions 

#Declaration of matrixes to store results
x=zeros((T,shape(guess0)[1]))
u=zeros((T,dim))
x_s=zeros((T,shape(guess0)[1],draws))
u_s=zeros((T,dim,draws))
shocks=zeros((T,draws))
out_of_bounds=0

for j in list(range(draws)):
    #Initialization of state variables
    u[[0],:]=u0
    #Iteration
    guess=guess0
    for t_step in list(range(1,T)):
        print('simulation in draw '+str(j+1)+' and step '+str(t_step))
        if j!=0:
            guess=x_s[t_step,:,j-1].reshape(shape(guess0))
        #Compute optimal control variable (consumption and abattement) at time t by the minimization of Bellman function
        (xopt,fopt, nb_iter, nb_funcalls, exitflag)=fmin(Bell_max,guess,xtol=tolmin, ftol=tolmin,full_output=1,disp=0)#xopt is the minimum and fopt the value of the function at its minimum
        #print('error: '+str(guess-xopt))        
        guess=xopt.copy()
        xopt=maximum(zeros(shape(xopt)),xopt)#transpose because x is vector
        xopt=minimum(ones(shape(xopt)),xopt)
        x[[t_step-1],:]=xopt
        if stochastic:
            min_sh=min_shock(u[[t_step-1],:])
            max_sh=max_shock(u[[t_step-1],:])
            mode_sh=mode_shock(u[[t_step-1],:])
            unif_stand=random.uniform(0,1)
            if (law_shock=="triang") or (law_shock=="rectang"):
                exec("eps=inv_repart_"+law_shock+"(unif_stand,min_sh,max_sh,mode_sh)")
            if law_shock=="betalaw":
                eps=fsolve(randombeta,0.5*(min_sh+max_sh))
                eps=max(min(eps,max_sh),min_sh)
            u[[t_step],:]=law_motion(t_step,x[[t_step-1],:],u[[t_step-1],:],eps)
        else:
            eps=0
            u[[t_step],:]=law_motion(t_step,x[[t_step-1],:],u[[t_step-1],:],eps)	
        shocks[t_step-1,j]=eps
        if (u[[t_step],:]>u_max[[t_step],:]).any() or (u[[t_step],:]<u_min[[t_step],:]).any():
            print("state variable out of bounds in simulate, draw "+str(j)+", time "+str(t_step))
            savetxt(run_folder+'simulate_break.csv', [j,t_step], delimiter=";")
            out_of_bounds=1
            break
	
    if j==0:
        x_s=x.reshape((T,shape(guess0)[1],1))
        u_s=u.reshape((T,dim,1))
    else:
        x_s=concatenate((x_s,x.reshape((T,shape(guess0)[1],1))),axis=2)
        u_s=concatenate((u_s,u.reshape((T,dim,1))),axis=2)
    
    if out_of_bounds==1:
        break
        
x_m=mean(x_s,axis=2)
u_m=mean(u_s,axis=2)

