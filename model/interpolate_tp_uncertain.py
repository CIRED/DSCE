# -*- coding: utf-8 -*-

m=deg+1
#number of terms in the complete Chebyshev approximation
nb_terms_cheb=n_terms_cheb(deg,dim)
nb_terms_cheb=int(nb_terms_cheb)

#matrix of weight (size: nb_terms_cheb *d)
#each line correspond to a multi-indice
weight_cheb=fill_cheb(deg,dim)

#indicatrice of the multi-indice (size nb_terms_cheb *1)
#= number of positive indices
d_cheb=ind_cheb(weight_cheb)
coef_V_uncertain=zeros((T,nb_terms_cheb))

#setting the boundary of approximation
u=zeros((T,dim))
z=nodes_Chebyshev(dim,m)
z_adjust=-cos(pi/(2*m))*ones((1,dim))
guess=guess0
guess_out=guess0

guess_store=zeros((T,shape(guess0)[1],m**dim))
error_guess=zeros((T,shape(guess0)[1],m**dim))
guess_store_out=zeros((T,shape(guess0)[1],m**dim))
error_guess_out=zeros((T,shape(guess0)[1],m**dim))
guess_bound=array([0.5])
print('interpolation uncertain')

for t_step in list(range(T,0,-1)):
    #print('interpolation in uncertain zone in step '+str(t_step))
    u_adjust=(z_adjust+1)/(-2*z_adjust)*(u_max_uncertain[[t_step-1],:]-u_min_uncertain[[t_step-1],:])

    for i in list(range(m**dim)):
        if t_step!=T:
            if i%m==0:
                guess=guess_store[t_step,:,i]
                guess_out=guess_store_out[t_step,:,i]
                
            else:
                guess=guess_store[t_step-1,:,i-1]+guess_store[t_step,:,i]-guess_store[t_step,:,i-1]
                guess_out=guess_store_out[t_step-1,:,i-1]+guess_store_out[t_step,:,i]-guess_store_out[t_step,:,i-1]

        #these are the point of chebyshev
        u[[t_step-1],:]=(z[[i],:]+1)/2 *(u_max_uncertain[[t_step-1],:]-u_min_uncertain[[t_step-1],:]+2*u_adjust)+u_min_uncertain[[t_step-1],:]-u_adjust
        (xopt,fopt, nb_iter, nb_funcalls, exitflag)=fmin(Bell_max_uncertain_2,guess,xtol=tolmin, ftol=tolmin,full_output=1,disp=0)#xopt is the minimum and fopt the value of the function at its minimum
        xopt=maximum(zeros(shape(xopt)),xopt)
        xopt=minimum(ones(shape(xopt)),xopt)
        approx_value=-Bell_max_uncertain_2(xopt)
        guess_store[t_step-1,:,i]=xopt.copy()
        error_guess[t_step-1,:,i]=guess-xopt
        guess=xopt.copy()
        (yopt,fopt, nb_iter, nb_funcalls, exitflag)=fmin(Bell_max_uncertain_boundary,guess_bound,xtol=tolmin, ftol=tolmin,full_output=1,disp=0)#xopt is the minimum and fopt the value of the function at its minimum
        xopt_zero=array([[yopt[0],1.]])
        yopt=maximum(zeros(shape(yopt)),yopt)
        yopt=minimum(ones(shape(yopt)),yopt)
        guess_bound=yopt.copy()
        approx_value_zero=-Bell_max_uncertain_boundary(yopt)

        if approx_value_zero>approx_value:
            approx_value=approx_value_zero
            
        for j in list(range(nb_terms_cheb)):
            coef_V_uncertain[t_step-1,j]=coef_V_uncertain[t_step-1,j]+2**(d_cheb[[j],:])/(m**dim)*approx_value*Chebyshev(weight_cheb[[j],:],z[[i],:])
