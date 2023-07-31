# -*- coding: utf-8 -*-

#Functions for approximation with chebyshev polynomials

#number of Chebyshev polynomials of total degree less or equal to deg in dimension dim
#z, deg, dim are positive integers
def n_terms_cheb(deg,dim):
    z=prod(list(range(1+deg,1+deg+dim)))/prod(list(range(1,1+dim)))
    return z;

#indicator of a weight (number of positive partial degree)
#weight is lw*dim matrix
#d is lw*1 column
def ind_cheb(weight):
    d=minimum(weight,ones(shape(weight)))
    d=sum(d,axis=1).reshape((shape(d)[0],1))
    return d;

#matrix of weight of Chebyshev polynomials of total degree less or egal to deg in dimension d
#deg is integer
#dim is positive integer
#weight_cheb is n_terms_cheb(deg,dim)*dim matrix 
def fill_cheb(deg,dim):
    if dim==1:
        weight_cheb=arange(deg+1).reshape((deg+1,1))
    elif deg==0:
        weight_cheb=zeros((1,dim))
    else:
        weight_cheb1=fill_cheb(deg,dim-1)
        s1=shape(weight_cheb1)[0]
        weight_cheb1=concatenate((zeros((s1,1)),weight_cheb1),axis=1)
        weight_cheb2=fill_cheb(deg-1,dim)
        s2=shape(weight_cheb2)[0]
        weight_cheb2=concatenate((ones((s2,1)),zeros((s2,dim-1))),axis=1)+weight_cheb2
        weight_cheb=concatenate((weight_cheb1,weight_cheb2),axis=0)
    return weight_cheb;

#alternative methods to compute Chebyshev polynomials in two steps
#first in one dimension
#weight is integer
#x is lx*1 column
def Chebyshev_one_dim(weight,x):
    if weight==0:
        z=1
    elif weight==1:
        z=x
    elif weight==2:
        z=2*x**2-1
    elif weight==3:
        z=4*x**3-3*x
    elif weight==4:
        z=8*x**4-8*x**2+1
    elif weight==5:
        z=16*x**5-20*x**3+5*x
    else: 
        z=2*x*Chebyshev_one_dim(weight-1,x)-Chebyshev_one_dim(weight-2,x)
    return z;

#then in any dimension dim
#weight is 1*dim row
#x is lx*dim matrix
#z is lx*1 column
def Chebyshev(weight,x):
    z=0
    d=shape(weight)[1]
    if d==shape(x)[1]:
        z=1
        for i in list(range(d)):
            z=z*Chebyshev_one_dim(weight[:,[i]],x[:,[i]])	
    else:
        print('Incompatible dimensions in Chebyshev')
    return z;

#computes the value at x of the approximation by coefficients coef of Chebyshev polynomials of degree deg in dimension dim on [x_min,x_max]
#x, x_min, x_max are lx *dim matrixes
#coef is matrix of coefficients of approximation. Each column corresponds to a weight. If coef has several rows, each row represents the approximation of a function.
#it is either of size 1*nb_terms_cheb, in which case the approximate function is the same for each row of x
#or it is of size lx*nb_terms_cheb, in which case the approximate function is different for each row of x
#z is of size lx*1
def approx_Chebyshev(coef,weight,x,x_min,x_max):
    z=0
    nb_terms=shape(coef)[1]
    if  shape(weight)[0]==nb_terms:
        for i in list(range(nb_terms)):
            z=z+coef[:,[i]]*Chebyshev(weight[[i],:],(2*x-x_min-x_max)/(x_max-x_min))
    else:
        print('check the Chebyshev approximation, something wrong in dimensions')
    return z;


#calculates the m Chebyshev nodes in dimension d
#z is a m^d*d matrix
def nodes_Chebyshev(d,m):
    if d==1:
        ind=arange(1,m+1).reshape((m,1))
        z=-cos((2*ind-1)*pi/(2*m)) 
    else:   
        z2=nodes_Chebyshev(d-1,m)
        #boustrophedon way through the nodes 
        #at each step, only one coordinate changes
        if m%2==0:
            z2=kron(ones((m//2,1)),concatenate((z2,z2[::-1])))
        elif m%2==1:
            z2=concatenate((kron(ones(((m-1)//2,1)),concatenate((z2,z2[::-1]))),z2))
        else:
            print('Disp problems with m: non integer number of nodes')
        z1=nodes_Chebyshev(1,m)
        z1=kron(z1,ones((m**(d-1),1)))
        z=concatenate((z1,z2),axis=1)
    return z;

#adjusts the lower bound of a domain [x_min,x_max] with m nodes
#m is positive integer
#x_min,x_max are 1*dim rows
#v is 1*dim row
def adjust_min(x_min,x_max,m):
    z_adjust=-cos(pi/(2*m))
    x_adjust=(z_adjust+1)/(-2*z_adjust)*(x_max-x_min)
    v=x_min-x_adjust
    return v;

#adjusts the upper bound of a domain [x_min,x_max] with m nodes
#m is positive integer
#x_min,x_max are 1*dim rows
#v is 1*dim row
def adjust_max(x_min,x_max,m):
    z_adjust=-cos(pi/(2*m))
    x_adjust=(z_adjust+1)/(-2*z_adjust)*(x_max-x_min)
    v=x_max+x_adjust
    return v;

#penality of x outside [x_min,x_max]
#x,x_min,x_max of size lx*dim
#e of size lx*1
def out_interval(x,x_min,x_max):
    if shape(x)[1]==3:
        x=x[:,0:2]
        x_min=x_min[:,0:2]
        x_max=x_max[:,0:2]
    x_save=x.copy()
    x=maximum(x,x_min)
    x_err=x-x_save
    x_save=x.copy()
    x=minimum(x,x_max)
    x_err=x_save-x+x_err
    x_err=x_err/(x_max-x_min)
    e=sum(x_err,axis=1).reshape((shape(x_err)[0],1))
    e=e+e**2   
    return e;

def Bell_max(x):
    #caution: in help, it is written that x should be a column vector. If this is true then x should be transpose for x is a line in the model functions;
    x_step=x.reshape(shape(guess0))
    if t_step>T:
        print('Something wrong: t_step is beyond time horizon')
    else:
        #the state variable
        u_step=u[[t_step-1],:]

    # checking if t_step is final 
    if t_step==T:
        eps=0
        un=law_motion(t_step,x_step,u_step,eps)

        #value with penalities for x_step outside bounds
        v1=maximand(array([t_step-1]),x_step,u_step,eps)+bet*terminal_value(t_step+1,un)-out_penality*out_interval(x_step,zeros(shape(x_step)),ones(shape(x_step)))

    #expected value
        v=float(-v1)#minus because we will minimize this function (so as to maximize the Bellman function)

    elif t_step<T:

        #coefficient approximating the next step value function
        coefn=coef_V[[t_step],:]
        eps=0
        un=law_motion(t_step,x_step,u_step,eps)
        un_min=u_min[[t_step],:]
        un_max=u_max[[t_step],:]

        #next-stage state variable should be in the domain of the approximation of next-stage value function
        un_adj=minimum(maximum(un,un_min),un_max)
        un_min_adj=adjust_min(un_min,un_max,m)
        un_max_adj=adjust_max(un_min,un_max,m)
        v1=maximand(array([t_step-1]),x_step,u_step,eps)+bet*approx_Chebyshev(coefn,weight_cheb,un_adj,un_min_adj,un_max_adj)-out_penality*out_interval(x_step,zeros(shape(x_step)),ones(shape(x_step)))-out_penality*out_interval(un,un_min,un_max)

    #expected value
        v=float(-v1)#minus because we will minimize this function (so as to maximize the Bellman function)
    return v;
