from casadi import *
#To do: how to get the electric boilers implemented in just the thermal system now??? ->
# You could for instance use the electric boiler and assume a cost for electricity? 
#Then it is very likely (depending on the cost) that either the DG or the EB is preferred based on cost

T = 24.0*60*60 #changing from 10 to 24 hours to be "realistic"#10. # Time horizon
N = 144 # number of control intervals ten minutes in hourform


# Declare model variables
x1 = MX.sym('x1')
x2 = MX.sym('x2')
x3 = MX.sym('x3')
x4 = MX.sym('x4')
x = vertcat(x1, x2, x3, x4)
u1 = MX.sym('u1') #capacity percentage DG
u2 = MX.sym('u2') #capacity percentage boiler
u3 = MX.sym('u3') #percentage of total mass flow flowing into DG
u4 = MX.sym('u4') #percentage of total mass flow flowing into boiler
u= vertcat(u1,u2,u3,u4)

#Parameters in my system
cp=4.186 # Joule/(gram*degree celcius) = kJoule/kg*degree celcius
rho=1.0  #1000 kg/m3 = 1 kg/dm3 = 1 kg/l
c_tes=3.0
c_dg_heatcap = 2000 #000 
c_boiler_heatcap = 2000 #000                                                                                                                                                                                                                                     
c_rad_heatcap=3 #000
V_tes=12000 #liter 
V_dg= 1
V_boiler= 1
V_rad = 1000
#the unit is kilo/sec (liter/sec)
w_tot= 2 # 2l/s
#in the last system several of the power-funtions were functions of t, how to do this descrete?
q_loss= 0.1



def q_sorad(): 
    return 14.0 ## to get it from minutes to seconds!!!
q_rad=q_sorad() #defining heat being used in the house...

#defining the heat that the DG returns when on (being used at 70%): 
DG_heat = 22.0 # kW heat kiloJoule/s ##########BECAUSE SAMPLING EVERY TEN MINUTES!
DG_el = 140.0 # kW power
B_MAX_HEAT = 128.0#kW heat?

q_DG= u1*DG_heat
q_BOILER= u2*B_MAX_HEAT

t_mix_ratio = u3*x1 + x4*(1-u3)
t_in_tes_ratio = u4*x2 + t_mix_ratio*(1-u4) #viktig at det ikke bare er x1(1-u4)
xdot= vertcat((cp*u3*w_tot*(x4-x1)*90 + q_DG-q_loss)/(rho*c_dg_heatcap*V_dg*90),
               (cp*u4*w_tot*(t_mix_ratio-x2)*90 -q_loss + q_BOILER)/(rho*c_dg_heatcap*V_boiler*90),
                  (cp*w_tot*(t_in_tes_ratio - x3)*90 -q_loss)/(rho*V_tes*c_tes*90),
                   (cp*w_tot*(x3-x4)*90 - q_rad - q_loss)/(rho*c_dg_heatcap*V_rad*90)) #now state nr two is the temperature of water coming out of the house
#this will depend on wht the temperature of the house and what heat the user decides to put on.....

#weighing of the different components of the objective function...
c_X3=20.5 
c_x1=0.0
c_x2=0.0

c_boiler=0.01
c_co2=0.01 #seeing what the temperatures end up with now
#reference temperatures to ensure high enough temperature in the "house", still don't know what these bounds should be...
x3_ref=65.0/90
x1_ref=75.0/90
x2_ref=75.0/90


# Objective term -> uttrykk for cost-funksjon
L= u1**2*c_co2  + c_X3*(x3-x3_ref)**2 + c_boiler*u2**2  + c_x1*(x1-x1_ref)**2 + c_x2*(x2-x2_ref)**2  #+ u3**2*c_co2 + u4**2*c_boiler
#-u3**2 -u4**2# u3 og u4 er for å se om det går noenting gjennom der nå...
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
   X0 = MX.sym('X0', 4) #init state,
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
Fk = F(x0=[62.0/90, 62.0/90, 62.0/90, 62.0/90 ],p=0.4) #tror dette er startpunkt men må sjekke ut?
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
Xk = MX.sym('X0', 4) #changed to three since we now have three x-es!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
w += [Xk]
#lbw += [67.0, 66.5, 66.0, 53.0 ]  #init conditions!!!!!!!!!!!!!!
lbw +=[62.0/90, 62/90, 62.0/90, 62.0/90 ]
ubw += [62.0/90, 62/90, 62.0/90, 62.0/90 ]
w0 += [62.0/90, 62.0/90, 62.0/90, 62.0/90 ] #her begynner casaadi å søke, må være feasible!!!!, ikke del på null

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
    Xk = MX.sym('X_' + str(k+1), 4) 
    w   += [Xk]
    
    #changed now to have more realistic limits, let's see!!
    lbw += [40.0/90, 40.0/90, 40.0/90, 30.0/90 ] # temperatur av TES skal eeeeeeeeegt ikke gå lavere enn 65 men tester dette....
    ubw += [90.0/90, 90.0/90, 90.0/90, 90.0/90 ]
    w0 += [62.0/90, 62.0/90, 62.0/90, 62.0/90 ]  #endret her til 78 på dg
  
    # Add equality constraint
    g   += [Xk_end-Xk] #blå minus rød fra video, multiple shoot constrainten!!! bruker g for vanlige constraints også
    #både equality og inequality constraints havner her, om ubegrenset: upper bound uendelig for eksempel
    lbg += [0, 0, 0, 0]
    ubg += [0, 0, 0, 0] #changed this from [0,0,0] and it now evaluates objective function more times...

# Create an NLP solver
prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)} #kom tilbake til parametere som varierer, om de inngår i difflikningene
solver = nlpsol('solver', 'ipopt', prob);

# Solve the NLP, den initialiseres
sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
w_opt = sol['x'].full().flatten()

# Plot the solution, sånn henter man ut variablene
x1_opt = w_opt[0::8]
x2_opt = w_opt[1::8]
x3_opt = w_opt[2::8]
x4_opt = w_opt[3::8]
u1_opt = w_opt[4::8]
u2_opt = w_opt[5::8]
u3_opt = w_opt[6::8]
u4_opt = w_opt[7::8]


tgrid = [T/N*k for k in range(N+1)]
t_values = np.linspace(0, T, N+1)  # Adjusted initialization

import matplotlib.pyplot as plt
plt.figure(1)
plt.clf()
plt.plot(t_values/(60*60), x1_opt*90, '--')
plt.plot(t_values/(60*60), x2_opt*90, '-')
plt.plot(t_values/(60*60),x3_opt*90, '.')
plt.plot(t_values/(60*60),x4_opt*90, '.-')
plt.ylabel('Temperature degrees celcius')
plt.xlabel('T: hours')
plt.legend(['x1:temp water dg','x2:temp water boiler','x3:temp water tes','x4:temp water ahouse'])
#prøver å få det til to plots, let's see ...
plt.figure(2)
plt.step(t_values/(60*60), vertcat(DM.nan(1), u1_opt), '-.') #her i plottingen kan det være vanskelig å få det riktig hmmmm....
plt.step(t_values/(60*60), vertcat(DM.nan(1), u2_opt), '-.') #prøver å få plotta alle u-ene, får se hva som skjer...
plt.step(t_values/(60*60), vertcat(DM.nan(1), u3_opt), '-.')
plt.step(t_values/(60*60), vertcat(DM.nan(1), u4_opt), '-.')
plt.ylabel('Percentage')
plt.xlabel('T: hours')
plt.legend(['u1:power_%_DG','u2:power_%_boiler','u3:%_mass_flow_DG','u4:%_mass_flow_boiler'])
plt.grid()
plt.show()
