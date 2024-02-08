from casadi import *

T = 24.0 #changing from 10 to 24 hours to be "realistic"#10. # Time horizon
N = 20 # number of control intervals

# Declare model variables
x1 = MX.sym('x1')
x2 = MX.sym('x2')
x3 = MX.sym('x3')
x = vertcat(x1, x2,x3)
u = MX.sym('u')
#how to get in parameters in the best way? 
#and how to make the timeseries of weather data?
#Parameters in my system
rho=1.0
a_loss=3.0
T_outside=-10 #keeping this constant for now...
c_house=80 #double check if this should be changed..
beta=4
E_bmin=20
Ebmax=80
V=1000
w_tes=2 #denne må endres til å kunne variere med tanke på pumpen, dette er jo totale gjennomstrømningen i systemet...
pd_max=200 #kWh
TES_max=80 #degrees celcius
TES_min= 50 #degrees celcius
#in the last system several of the power-funtions were functions of t, how to do this descrete?
def wind():
    return 20 
def solar():
    return 200
def elecload():
    return 230
def elkjel():
    return 116
def q_l():
    return 0.1

pw=wind()
ppv=solar()
pl=elecload()
pel=elkjel()
q_loss=q_l()

#does this now work or does it need to be in a function?? 
#def make_xdot():
    #if (x2< E_bmin) or (x1<19.0):
       # pd=250
    #elif (x2>Ebmax-2.0) and (x1>22.0):
        #pd=0
    #else:
       # pd=200
        #Power balance for battery
    
#defining the pd value so based on u being from 0 to 1 we have our pd value
pd= pd_max*u

pb= pd+pw+ppv-pl-pel #power balance for the battery: important-> pd=u now since we control everything throug pd!!

#defining the temperature of the water entering the TES based on the power produced in the DG
def temp_in(u):
    t_in= TES_min + u*(TES_max-TES_min)
    return t_in


t_in_tes=temp_in(u)

#Model equations!!!!!!!!!
xdot= vertcat((w_tes*(t_in_tes-x1)-q_loss)/(rho*V), (w_tes*(x1-x2)-a_loss*(x2-T_outside))/c_house , pb/beta )
    
#T_inn=T_in(pd) #need to create a pd-function possibly...

#xdot = make_xdot()
r_house=22.0 #reference temperature for ideal housetemperature...
c_h=20.0 #weighing of the different components of the objective function...
c_el=2.0
c_co2=10.0
# Objective term -> uttrykk for cost-funksjon
L= u**2*c_co2 + (pd + pw+ ppv + pb - pl - pel)**2*c_el + (x2 -r_house)**2*c_h
#L = x1**2 + x2**2 + x3**2 + u**2 # for minst cost må x1,x2 og u lik null! Her må jeg implementere en faktisk relevant cost-funksjon!!!!

# Formulate discrete time dynamics
if False:
   # CVODES from the SUNDIALS suite
   dae = {'x':x, 'p':u, 'ode':xdot, 'quad':L}
   F = integrator('F', 'cvodes', dae, 0, T/N)
else:
   # Fixed step Runge-Kutta 4 integrator
   M = 4 # RK4 steps per interval
   DT = T/N/M #dette er time-step
   f = Function('f', [x, u], [xdot, L]) #f(xk,uk), dette er funksjoin man vil integrere, ender opp med x_dot og L
   X0 = MX.sym('X0', 3) #init state, changed from 2 to 3 since we have  x-es!!!!
   U = MX.sym('U')
   X = X0
   Q = 0
   for j in range(M): #mer nøyaktig versjon, er runge kutta 4, de ulike k-likningene osv, her finner vi neste state
       k1, k1_q = f(X, U)
       k2, k2_q = f(X + DT/2 * k1, U)
       k3, k3_q = f(X + DT/2 * k2, U)
       k4, k4_q = f(X + DT * k3, U)
       X=X+DT/6*(k1 +2*k2 +2*k3 +k4)
       Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q)
   F = Function('F', [X0, U], [X, Q],['x0','p'],['xf','qf']) #xf: x final, Qf : final cost, F er integralet av f

# Evaluate at a test point
Fk = F(x0=[50.0,20.0,40.0],p=0.4) #tror dette er startpunkt men må sjekke ut?
#Fk = F(x0=[0.2,0.3,],p=0.4)
print(Fk['xf'])
print(Fk['qf'])

# Start with an empty NLP, initialiserer et tom nlp for å bruke multiple shooting
w=[]
w0 = []
lbw = []
ubw = []
J = 0
g=[]
lbg = []
ubg = []

# "Lift" initial conditions, vet hva init tilstander er og lager equality constraint, casadi har ikke dette, men sier at 
#lb og ub er det samme
Xk = MX.sym('X0', 3) #changed to three since we now have three x-es!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
w += [Xk]
lbw += [50.0, 20.0 , 40.0] #init conditions!!!!!!!!!!!!!!
ubw += [50.0, 20.0 , 40.0]
w0 += [50.0, 20.0 , 40.0] #her begynner casaadi å søke, må være feasible!!!!, ikke del på null

# Formulate the NLP
for k in range(N):
    # New NLP variable for the control
    Uk = MX.sym('U_' + str(k))
    w   += [Uk]
    lbw += [-1] #dette er grensene for u 
    ubw += [1] #w er decision variable, xuuuxuxuxu #trying to see if u gets bigger now
    w0  += [0]

    # Integrate till the end of the interval
    Fk = F(x0=Xk, p=Uk) #x-en på slutt av første intervalll
    Xk_end = Fk['xf']
    J=J+Fk['qf'] #inkrementerer cost

    # New NLP variable for state at end of interval
    Xk = MX.sym('X_' + str(k+1), 3) #changed this one from 2 to 3 because we now have 3 x-es!!!!!!!!!!!!!
    w   += [Xk]
    
    #changed now to have more realistic limits, let's see!!
    lbw += [50.0, 18.0, 20.0] #må adde en tredje her siden tre states!!
    ubw += [ 80.0, 25.0,90.0] 
    w0+= [50.0, 20.0 , 40.0] #[65.0, 20.0 , 70.0] #this does not affect, what happens here???????????????????????????????????????????
   # w0  += [0.0,0.0,0.0] ### Is this where you want it to end??
    # Add equality constraint
    g   += [Xk_end-Xk] #blå minus rød fra video, multiple shoot constrainten!!! bruker g for vanlige constraints også
    #både equality og inequality constraints havner her, om ubegrenset: upper bound uendelig for eksempel
    lbg += [0, 0, 0]
    ubg += [0, 0, 0] #changed this from [0,0,0] and it now evaluates objective function more times...

# Create an NLP solver
prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)} #kom tilbake til parametere som varierer, om de inngår i difflikningene
solver = nlpsol('solver', 'ipopt', prob);

# Solve the NLP, den initialiseres
sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
w_opt = sol['x'].full().flatten()

# Plot the solution, sånn henter man ut variablene
x1_opt = w_opt[0::4]
x2_opt = w_opt[1::4]
x3_opt = w_opt[2::4] # adder tilstand nr 3 her, satser p åat det funker...
u_opt = w_opt[3::4]

tgrid = [T/N*k for k in range(N+1)]

import matplotlib.pyplot as plt
plt.figure(1)
plt.clf()
plt.plot(tgrid, x1_opt, '--')
plt.plot(tgrid, x2_opt, '-')
plt.plot(tgrid,x3_opt, '.')
plt.step(tgrid, vertcat(DM.nan(1), u_opt), '-.')
plt.xlabel('t')
plt.legend(['x1','x2', 'x3','u'])
plt.grid()
plt.show()