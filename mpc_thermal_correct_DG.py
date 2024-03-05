from casadi import *
#To do: how to get the electric boilers implemented in just the thermal system now??? ->
# You could for instance use the electric boiler and assume a cost for electricity? 
#Then it is very likely (depending on the cost) that either the DG or the EB is preferred based on cost

T = 24.0 #changing from 10 to 24 hours to be "realistic"#10. # Time horizon
N = 20 # number of control intervals

# Declare model variables
x1 = MX.sym('x1')
x2 = MX.sym('x2')
x = vertcat(x1, x2)
u1 = MX.sym('u1') #capacity percentage DG
u2 = MX.sym('u2') #capacity percentage boiler
u3 = MX.sym('u3') #percentage of total mass flow flowing into DG
u4 = MX.sym('u4') #percentage of total mass flow flowing into boiler
u= vertcat(u1,u2,u3,u4)
#how to get in parameters in the best way? 
#and how to make the timeseries of weather data?
#Parameters in my system
rho=1.0
#a_loss=4.0
#c_house= 3000#250 #80 #double check if this should be changed..
V_rad=3000
V=12000 #1000
w_tes=200 #denne må endres til å kunne variere med tanke på pumpen, dette er jo totale gjennomstrømningen i systemet...
#pd_max=200 #kWh
DG_max=85 #degrees celcius
TES_min= 50 #degrees celcius
boiler_max=85
boiler_min=50
#in the last system several of the power-funtions were functions of t, how to do this descrete?

def q_l():
    return 0.1

q_loss=q_l()
def q_sorad(): 
    return 300.0
q_rad=q_sorad() #defining heat being used in the house...

#temperature of water passing through the DG
def temp_dg():
    t_in= (x2+u1*(DG_max-x2))*u3
    return t_in
t_dg=temp_dg()
#temperature of water passing by the DG (not entering and being heated)
t_not_dg=x2*(1-u3)
#total temp og mixed water heated and not heated by the DG
t_mix = t_dg + t_not_dg
#temperature of the mixed water after being heated by boiler
def t_boile():
    t_in=u4*(t_mix+u2*(boiler_max-t_mix))
    return t_in
t_boil= t_boile()
#water bypassing the boiler
t_not_boil= t_mix*(1-u4)

#all the water going into the TES
t_in_tes=t_boil+t_not_boil

xdot= vertcat((w_tes*(t_in_tes-x1)-q_loss)/(rho*V), (w_tes*(x1-x2)-q_rad)/V_rad) #now state nr two is the temperature of water coming out of the house
#this will depend on wht the temperature of the house and what heat the user decides to put on.....


#when c_h=0: seeing if the
c_X1=0.5 #weighing of the different components of the objective function...
C_X2=0.5
c_boiler=3.5
c_co2=2.0 #seeing what the temperatures end up with now
#reference temperatures to ensure high enough temperature in the "house", still don't know what these bounds should be...
x1_ref=71.0
x2_ref=64.0

# Objective term -> uttrykk for cost-funksjon
L= u1**2*c_co2 + C_X2*(x2 - x2_ref)**2 + c_X1*(x1-x1_ref)**2 + c_boiler*u2**2 
#here in the objective function, the usage of the DG and boiler is punished, but not the water flowing through
#I think this makes sense, but might need to be looked at...


# Formulate discrete time dynamics|
if False:
   # CVODES from the SUNDIALS suite
   dae = {'x':x, 'p':u, 'ode':xdot, 'quad':L}
   F = integrator('F', 'cvodes', dae, 0, T/N)
else:
   # Fixed step Runge-Kutta 4 integrator
   M = 4 # RK4 steps per interval
   DT = T/N/M #dette er time-step
   f = Function('f', [x, u], [xdot, L]) #f(xk,uk), dette er funksjoin man vil integrere, ender opp med x_dot og L
   X0 = MX.sym('X0', 2) #init state,
   U = MX.sym('U', 4) #sier her at det er fire u-er!!!
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
Fk = F(x0=[67.0,55.0],p=0.4) #tror dette er startpunkt men må sjekke ut?
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
Xk = MX.sym('X0', 2) #changed to three since we now have three x-es!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
w += [Xk]
lbw += [67.0, 55.0 ] #init conditions!!!!!!!!!!!!!!
ubw += [67.0, 55.0 ]
w0 += [67.0, 55.0 ] #her begynner casaadi å søke, må være feasible!!!!, ikke del på null

# Formulate the NLP
for k in range(N):
    # New NLP variable for the control
    Uk = MX.sym('U_' + str(k), 4) # have to put ,4 in this formulation????????????????????????????????????????????????????
    w   += [Uk]
    lbw += [0,0,0,0] #dette er grensene for u (here it is taken into account that there ar 4 u's)
    ubw += [1,1,1,1] #w er decision variable, xuuuxuxuxu #trying to see if u gets bigger now
    w0  += [0,0,0,0]

    # Integrate till the end of the interval
    Fk = F(x0=Xk, p=Uk) #x-en på slutt av første intervalll
    Xk_end = Fk['xf']
    J=J+Fk['qf'] #inkrementerer cost

    # New NLP variable for state at end of interval
    Xk = MX.sym('X_' + str(k+1), 2) 
    w   += [Xk]
    
    #changed now to have more realistic limits, let's see!!
    lbw += [40.0, 20.0] #må adde en tredje her siden tre states!!
    ubw += [ 80.0, 70.0] 
    w0+= [67.0, 55.0 ] #
  
    # Add equality constraint
    g   += [Xk_end-Xk] #blå minus rød fra video, multiple shoot constrainten!!! bruker g for vanlige constraints også
    #både equality og inequality constraints havner her, om ubegrenset: upper bound uendelig for eksempel
    lbg += [0, 0]
    ubg += [0, 0] #changed this from [0,0,0] and it now evaluates objective function more times...

# Create an NLP solver
prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)} #kom tilbake til parametere som varierer, om de inngår i difflikningene
solver = nlpsol('solver', 'ipopt', prob);

# Solve the NLP, den initialiseres
sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
w_opt = sol['x'].full().flatten()

# Plot the solution, sånn henter man ut variablene
x1_opt = w_opt[0::6]
x2_opt = w_opt[1::6]
#x3_opt = w_opt[2::4] # adder tilstand nr 3 her, satser p åat det funker...
u1_opt = w_opt[2::6]
u2_opt = w_opt[3::6]
u3_opt = w_opt[4::6]
u4_opt = w_opt[5::6]


tgrid = [T/N*k for k in range(N+1)]

import matplotlib.pyplot as plt
plt.figure(1)
plt.clf()
plt.plot(tgrid, x1_opt, '--')
plt.plot(tgrid, x2_opt, '-')
#plt.plot(tgrid,x3_opt, '.')
plt.step(tgrid, vertcat(DM.nan(1), u1_opt), '-.') #her i plottingen kan det være vanskelig å få det riktig hmmmm....
plt.step(tgrid, vertcat(DM.nan(1), u2_opt), '-.') #prøver å få plotta alle u-ene, får se hva som skjer...
plt.step(tgrid, vertcat(DM.nan(1), u3_opt), '-.')
plt.step(tgrid, vertcat(DM.nan(1), u4_opt), '-.')
plt.xlabel('t')
plt.legend(['x1:Temp_TES','x2:Temp_after_LOAD','u1:power_%_DG','u2:power_%_boiler','u3:%_mass_flow_DG','u4:%_mass_flow_boiler'])
plt.grid()
plt.show()