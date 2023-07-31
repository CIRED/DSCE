# -*- coding: utf-8 -*-


#reading general parameters
gamma=0
#subs=3
E0=222.1112799
ae=0.1
#T0=0.8

with open(run_folder+'param_gen.csv', 'r') as f:
    reader = csv.reader(f,  delimiter=';')
    for row in reader:
        if row[0]=='model':
            exec(row[0]+'=row[1]')
        elif row[0]=='stochastic' or row[0]=='draws' or row[0]=='deg':
            exec(row[0]+'=int(row[1])')
        else:
            exec(row[0]+'=float(row[1])')

#reading model parameters
with open(run_folder+'param_model.csv', 'r') as f:
    reader = csv.reader(f,  delimiter=';')
    for row in reader:
        if row[0]=='dim' or row[0]=='T' or row[0]=='deltaT':
            exec(row[0]+'=int(row[1])')
        else:
            exec(row[0]+'=float(row[1])')

#adjust exponential rates from 1 year to period (deltaT years)
rho=(1+rho1y)**deltaT-1
delta=1-(1-delta1y)**deltaT
bet=1/(1+rho)

#s is cumulated emissions in TtC
s0=T0/beta_temp#ratio between observed temperature increase in 2005 and beta_temp

#calculating trends
time_horizon=arange(T+2).reshape((T+2,1))#a little further than horizon for terminal value

# Another way to compute population (Nordhaus 2016)
L=POP0*ones((T+2,1))
for t in list(range(T+1)):
    L[[t+1],:]=(L[[t],:]*(POPASYM/L[[t],:])**GPOP0)

# New dynamics for productivity (Nordhaus 2016) / change production function accordingly
GA=GA0*exp(-DELA*deltaT*time_horizon)
A=A0*ones((T+2,1))
for t in list(range(T+1)):
    A[[t+1],:]=A[[t],:]/(1-GA[[t],:])
GA=GA/deltaT

#New emissions factor (Nordhaus 2016)
GSIGMA=-GSIGMA*ones((T+2,1))
for t in list(range(T+1)):
   GSIGMA[[t+1],:]=GSIGMA[[t],:]*((1-DSIG)**deltaT)
sigm=SIG0*ones((T+2,1))
for t in list(range(T+1)):
    sigm[[t+1],:]=sigm[[t],:]*exp(GSIGMA[[t],:]*deltaT)

#New abatement cost Nordhaus (2016)
PBACKTIME=PBACK*(1-GBACK)**(time_horizon)
theta1=PBACKTIME*sigm/theta2

    


