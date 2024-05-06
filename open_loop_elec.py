from casadi import *
#To do: how to get the electric boilers implemented in just the thermal system now??? ->
# You could for instance use the electric boiler and assume a cost for electricity? 
#Then it is very likely (depending on the cost) that either the DG or the EB is preferred based on cost

T = 24.0*60*60 #changing from 10 to 24 hours to be "realistic"#10. # Time horizon
N = 144 # number of control intervals ten minutes in hourform
def get_ppv():
    for k in range(N):
        if k<70:
            ppv=110
        else:
            ppv=10
        return ppv
    

# Declare model variables
x1 = MX.sym('x1')
x2 = MX.sym('x2')
x3 = MX.sym('x3')
x4 = MX.sym('x4')
x5 = MX.sym('x5')
x = vertcat(x1, x2, x3, x4, x5) #Adding a new state here, x5 
u1 = MX.sym('u1') #capacity percentage DG
u2 = MX.sym('u2') #capacity percentage boiler
u3 = MX.sym('u3') #percentage of total mass flow flowing into DG
u4 = MX.sym('u4') #percentage of total mass flow flowing into boiler
u5 = MX.sym('u5') #Decisioni variable for the curtailed pv, "removing" excess pv power so battery does not get too much powah
u6 = MX.sym('u6') #decision variable to control pb
u= vertcat(u1,u2,u3,u4,u5,u6)

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
BOILER_MAX_HEAT = 128.0#kW heat?
PPV_MAX=275 #kw
Battery_max=240 #kW charging and discharging
q_DG= u1*DG_heat
q_BOILER= u2*BOILER_MAX_HEAT

#making a constraint that does not lest water being let in to the DG unless it is turned on: ##############
#u3<=u1
#Let's see if this changes things
t_mix_ratio = u3*x1 + x4*(1-u3)
t_in_tes_ratio = u4*x2 + t_mix_ratio*(1-u4) #viktig at det ikke bare er x1(1-u4)

#power balance for the electrical systems:  ############################### NEW ###########################
beta=  111.94 # 403*1000/(60*60), should be correct #4 #this can also be chaaaanged 
#pw=55 #kw just guessing for now
pl=90 #kw, this is what the hotel/house needs
ppv=get_ppv()
pcurt =u5*ppv
#pel=q_BOILER/0.98
pel=(u2*BOILER_MAX_HEAT)/0.98
pd=u1*DG_el
pb=u6*Battery_max
##################Hvordan definere dette annerledes?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?
Powerbalace = pb+ ppv+pd-pel-pl -pcurt
###################################### NEW """"""###############################################
xdot= vertcat((cp*u3*w_tot*(x4-x1)*90 + q_DG-q_loss)/(rho*c_dg_heatcap*V_dg*90),
               (cp*u4*w_tot*(t_mix_ratio-x2)*90 -q_loss + q_BOILER)/(rho*c_dg_heatcap*V_boiler*90),
                  (cp*w_tot*(t_in_tes_ratio - x3)*90 -q_loss)/(rho*V_tes*c_tes*90),
                   (cp*w_tot*(x3-x4)*90 - q_rad - q_loss)/(rho*c_dg_heatcap*V_rad*90),
                   ##### NNNNEEEEEEWWWW###############################
                   (pb/beta)) #this most likely needs to be scaled mehhhhhhhhhhhhhhhhhhhhhh
###################### NEW #################################################################



#weighing of the different components of the objective function...
c_X3=25.5 
c_x1=0.0
c_x2=0.0
c_powerbalance= 500.0
c_boiler=0.1
c_curt=0.1
c_co2=0.11 #seeing what the temperatures end up with now
#reference temperatures to ensure high enough temperature in the "house", still don't know what these bounds should be...
x3_ref=65.0/90
x1_ref=75.0/90
x2_ref=75.0/90

#added objective term to punish the powerbalance!!!!!!
# Objective term -> uttrykk for cost-funksjon
L= u1**2*c_co2  + c_X3*(x3-x3_ref)**2 + c_boiler*u2**2  + c_x1*(x1-x1_ref)**2 + c_x2*(x2-x2_ref)**2  +c_powerbalance*Powerbalace**2 +c_curt*u5**2#+ u3**2*c_co2 + u4**2*c_boiler
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
   X0 = MX.sym('X0', 5) #init state,
   U = MX.sym('U', 6) #sier her at det er fire u-er!!!
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
Fk = F(x0=[62.0/90, 62.0/90, 62.0/90, 62.0/90,62.0/90 ],p=0.4) #tror dette er startpunkt men må sjekke ut?
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
Xk = MX.sym('X0', 5) #changed to three since we now have three x-es!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
w += [Xk]
#lbw += [67.0, 66.5, 66.0, 53.0 ]  #init conditions!!!!!!!!!!!!!!
lbw +=[62.0/90, 62/90, 62.0/90, 62.0/90, 62.0/100 ]
ubw += [62.0/90, 62/90, 62.0/90, 62.0/90, 62.0/100 ]
w0 += [62.0/90, 62.0/90, 62.0/90, 62.0/90, 62.0/100 ] #her begynner casaadi å søke, må være feasible!!!!, ikke del på null

pb_values = []
pl_values = []
pd_values = []
pel_values = []
ppv_values = []
pcurt_values = []
powerbal_values = []


# Formulate the NLP
for k in range(N):
    # New NLP variable for the control
    Uk = MX.sym('U_' + str(k), 6) # have to put ,4 in this formulation????????????????????????????????????????????????????
    w   += [Uk]
    lbw += [0,0,0,0,0,-1] #-1 the limit for pb because it can be negative! (batteriet ledes ut...)
    ubw += [1,1,1,1,1,1] #w er decision variable, xuuuxuxuxu #trying to see if u gets bigger now
    
    w0  += [0,0,0,0,0,0]

    # Integrate till the end of the interval
    Fk = F(x0=Xk, p=Uk) #x-en på slutt av første intervalll
    Xk_end = Fk['xf']
    J=J+Fk['qf'] #inkrementerer cost

    # New NLP variable for state at end of interval
    Xk = MX.sym('X_' + str(k+1), 5) 
    w   += [Xk]
    
    #changed now to have more realistic limits, let's see!!
    lbw += [40.0/90, 40.0/90, 40.0/90, 30.0/90, 20.0/100 ] # temperatur av TES skal eeeeeeeeegt ikke gå lavere enn 65 men tester dette....
    ubw += [90.0/90, 90.0/90, 90.0/90, 90.0/90, 90.0/100 ]
    w0 += [62.0/90, 62.0/90, 62.0/90, 62.0/90, 62.0/100 ]  #endret her til 78 på dg
  
    # Add equality constraint
    g   += [Xk_end-Xk] #blå minus rød fra video, multiple shoot constrainten!!! bruker g for vanlige constraints også
    #både equality og inequality constraints havner her, om ubegrenset: upper bound uendelig for eksempel
    lbg += [0, 0, 0, 0, 0]
    ubg += [0, 0, 0, 0, 0] #changed this from [0,0,0] and it now evaluates objective function more times...
    if k<70:
        ppv=110
    else:ppv=10
    
     

#pb_values.append(f(X[0::5]))
#print("here is pb_values: ", pb_values, "here is the type: ",type(pb_values))


#print(" Anothe one, np_pb: ",np_pb," the type of it: ", type(np_pb), "length: ",len(np_pb))

    

# Create an NLP solver
prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)} #kom tilbake til parametere som varierer, om de inngår i difflikningene
solver = nlpsol('solver', 'ipopt', prob);
#print("here is length w0: ",len(w0),"here is length lbw: ",len(lbw), "here is length ubw: ", len(ubw), "here is length lbg: ",len(lbg),"here is length ubg: ",len(ubg))

# Solve the NLP, den initialiseres
sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
w_opt = sol['x'].full().flatten()

# Plot the solution, sånn henter man ut variablene
x1_opt = w_opt[0::11]
x2_opt = w_opt[1::11]
x3_opt = w_opt[2::11]
x4_opt = w_opt[3::11]
x5_opt = w_opt[4::11] #here is the battery state mehhhhh
u1_opt = w_opt[5::11]
u2_opt = w_opt[6::11]
u3_opt = w_opt[7::11]
u4_opt = w_opt[8::11]
u5_opt = w_opt[9::11]
u6_opt = w_opt[10::11]

#pel=(u2*B_MAX_HEAT)/0.98
#pd=u1*DG_el
#ppv=get_ppv()
for i in range(N):
    if i<70:
        ppv=110
    else: 
        ppv=10
    pb_values.append(u6_opt[i]*Battery_max)
    pel_values.append((BOILER_MAX_HEAT*u2_opt[i])/0.98)
    pd_values.append(u1_opt[i]*DG_el)
    ppv_values.append(ppv)
    pcurt_values.append(ppv*u5_opt[i])
    pl_values.append(pl)
    powerbal_values.append(u6_opt[i]*Battery_max+ ppv+u1_opt[i]*DG_el-(BOILER_MAX_HEAT*u2_opt[i])/0.98 - pl -ppv*u5_opt[i])

np_pb=np.array(pb_values)
#print("Here the final result of pb_values is: ", np_pb)
np_pl=np.array(pl_values)
np_pd=np.array(pd_values)
np_pel=np.array(pel_values)
np_ppv=np.array(ppv_values)
np_pcurt=np.array(pcurt_values)
np_powerbal=np.array(powerbal_values)
print('Here are the values of the powerbalance: ', np_powerbal)


tgrid = [T/N*k for k in range(N+1)]
t_values = np.linspace(0, T, N+1)  # Adjusted initialization
#print("LENGTH OF TVALUSS: ", len(t_values))
t_vals_pbosv = np.linspace(0, T, N)

import matplotlib.pyplot as plt
plt.figure(1)
plt.clf()
plt.plot(t_values/(60*60), x1_opt*90, '--')
plt.plot(t_values/(60*60), x2_opt*90, '-')
plt.plot(t_values/(60*60),x3_opt*90, '.')
plt.plot(t_values/(60*60),x4_opt*90, '.-')
plt.plot(t_values/(60*60), x5_opt*90,'.-.') #plotting the new state and scaling mehhhh
plt.ylabel('Temperature degrees celcius/battery percentage')
plt.xlabel('T: hours')
plt.legend(['x1:temp water dg','x2:temp water boiler','x3:temp water tes','x4:temp water ahouse','x5:percentage battery'])
#prøver å få det til to plots, let's see ...
plt.figure(2)
plt.step(t_values/(60*60), vertcat(DM.nan(1), u1_opt), '-.') #her i plottingen kan det være vanskelig å få det riktig hmmmm....
plt.step(t_values/(60*60), vertcat(DM.nan(1), u2_opt), '-.') #prøver å få plotta alle u-ene, får se hva som skjer...
plt.step(t_values/(60*60), vertcat(DM.nan(1), u3_opt), '-.')
plt.step(t_values/(60*60), vertcat(DM.nan(1), u4_opt), '-.')
plt.step(t_values/(60*60), vertcat(DM.nan(1), u5_opt), '--')
plt.step(t_values/(60*60), vertcat(DM.nan(1), u6_opt), '--')
#print('here is u6_opt lezzgo: ',u6_opt)

plt.ylabel('Percentage')
plt.xlabel('T: hours')
plt.legend(['u1:power_%_DG','u2:power_%_boiler','u3:%_mass_flow_DG','u4:%_mass_flow_boiler','u5: ppv curtailed','u6: pb'])
#plt.grid()
#plt.show()

plt.figure(3)
#plotting the new state and scaling mehhhh
plt.plot(t_vals_pbosv/(60*60), np_pb, '.' )
plt.plot(t_vals_pbosv/(60*60),np_pcurt, '--')
plt.plot(t_vals_pbosv/(60*60),np_pel, '-')
plt.plot(t_vals_pbosv/(60*60),np_pd,'.-')
plt.plot(t_vals_pbosv/(60*60),np_pl,'.-' )
plt.plot(t_vals_pbosv/(60*60), np_ppv, '-.')
plt.plot(t_vals_pbosv/(60*60), np_powerbal, '-.')
plt.ylabel('Power in kW')
plt.xlabel('T: hours')
plt.legend(['pb: power in/out of batttery','pcurt: curtailed PV power','pel: power used electrical boiler'
            , 'pd:power prod Diesel', 'pl: power used Isfjorden','ppv: power prod PV', 'powerbalance,0 ideally'])
plt.grid()
plt.show()
