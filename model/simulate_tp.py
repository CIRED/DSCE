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
info_tp=array([0,0])
import random
uncertainzone=0
havewetipped=0
table_xopt=zeros((T,4,2))
draws=1 #draw = 1 for now, stochastic draw later

for j in list(range(draws)):
    #Initialization of state variables
    u[[0],:]=u0
    
    if u0[0][1]>u_threshold_min:
        print('exploring uncertain zone')
        uncertainzone=1

        if u0[0][1]>u_threshold_max:
            print('tipping point has been crossed')
            havewetipped=1

    #Iteration
    guess=guess0
    guess_out=guess0
        
    for t_step in list(range(1,T)):
        print('simulation in draw '+str(j+1)+' and step '+str(t_step))
        if j!=0:
            guess=x_s[t_step,:,j-1].reshape(shape(guess0))
            
        #the code is : first variable indicates whether we have
        choice=array([0,0])
        
        #Compute optimal control variable (consumption and abattement) at time t by the minimization of Bellman function
        if uncertainzone==0:
            (xopt,fopt, nb_iter, nb_funcalls, exitflag)=fmin(Bell_max_pretipping_3interpolations,guess,xtol=tolmin, ftol=tolmin,full_output=1,disp=0)#xopt is the minimum and fopt the value of the function at its minimum
            guess=xopt.copy()
            xopt=maximum(zeros(shape(xopt)),xopt)#transpose because x is vector
            xopt=minimum(ones(shape(xopt)),xopt)
            value=-Bell_max_pretipping_3interpolations(xopt)
            print(xopt,value)

        elif havewetipped==0:
            (xopt,fopt, nb_iter, nb_funcalls, exitflag)=fmin(Bell_max_uncertain_2,guess,xtol=tolmin, ftol=tolmin,full_output=1,disp=0)#xopt is the minimum and fopt the value of the function at its minimum
            (yopt,fopt, nb_iter, nb_funcalls, exitflag)=fmin(Bell_max_uncertain_boundary,guess_bound,xtol=tolmin, ftol=tolmin,full_output=1,disp=0)#xopt is the minimum and fopt the value of the function at its minimum
            guess_bound=yopt.copy()
            yopt=maximum(zeros(shape(yopt)),yopt)#transpose because x is vector
            yopt=minimum(ones(shape(yopt)),yopt)
            xopt_zero=array([[yopt[0],1]])
            value_zero=-Bell_max_uncertain_boundary(yopt)
            guess=xopt.copy()
            xopt=maximum(zeros(shape(xopt)),xopt)#transpose because x is vector
            xopt=minimum(ones(shape(xopt)),xopt)
            value=-Bell_max_uncertain_2(xopt)

            if value_zero>value:
                xopt=xopt_zero.copy()
                print('zero is optimal')

        elif havewetipped==1:
            (xopt,fopt, nb_iter, nb_funcalls, exitflag)=fmin(Bell_max_posttipping,guess,xtol=tolmin, ftol=tolmin,full_output=1,disp=0)#xopt is the minimum and fopt the value of the function at its minimum
        x[[t_step-1],:]=xopt

        eps=0
        if uncertainzone==0:
            u[[t_step],:]=law_motion_pretipping(t_step,x[[t_step-1],:],u[[t_step-1],:],eps)
            u_step=u[[t_step],:]

            if u_step[0][1]>u_threshold_min:
                print('exploring uncertain zone')
                uncertainzone=1

        elif havewetipped==0:

            u[[t_step],:]=law_motion_pretipping(t_step,x[[t_step-1],:],u[[t_step-1],:],eps)
            u_step=u[[t_step],:]

            if u_step[0][1]>u_threshold_max:
                print('tipped point crossed!')
                havewetipped=1
                info_tp=array([1,t_step])

            prob=h(u[t_step-1,1],u_step[0][1],u_threshold_min,u_threshold_max)
            rnumber=random.random()

            if rnumber<prob:
                havewetipped=1
                info_tp=array([1,t_step])
                print('tipping point crossed')

        elif havewetipped==1:
            u[[t_step],:]=law_motion_posttipping(t_step,x[[t_step-1],:],u[[t_step-1],:],eps)
                        
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

savetxt(run_folder+'info_tp.csv', info_tp, delimiter=";", header='info on tipping point'+str(dim)+' stochastic '+str(stochastic))

x_m=mean(x_s,axis=2)
u_m=mean(u_s,axis=2)

