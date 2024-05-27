from casadi import *
#To do: how to get the electric boilers implemented in just the thermal system now??? ->
# You could for instance use the electric boiler and assume a cost for electricity? 
#Then it is very likely (depending on the cost) that either the DG or the EB is preferred based on cost
import casadi as ca
import numpy as np
T = 24.0*60*60 #changing from 10 to 24 hours to be "realistic"#10. # Time horizon
N = 144 # number of control intervals ten minutes in hourform
def get_ppv():
    for k in range(N):
        if k<70:
            ppv=100
        else:
            ppv=50
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

#A symbolic vector of the data of ppv power production
ppv = MX.sym('ppv')


Data_ppv_values = []
mu, sigma = 2, 0.5
s = np.random.normal(mu, sigma, 288)
for i in range(2*N):
   # Data_ppv_values.append(80+s[i])
    if i<70:
        Data_ppv_values.append(100)
    else: Data_ppv_values.append(50)

#print(np.array(v_list)) 
np_D_ppv=np.array(Data_ppv_values)
Data_ppv_values_plott = []
for i in range(N):
    Data_ppv_values_plott.append(80+s[i])
    #if i<70:
        #Data_ppv_values_plott.append(100)
    #else: Data_ppv_values_plott.append(50)

#print(np.array(v_list)) 
np_D_ppv_plott=np.array(Data_ppv_values_plott)

pb_values = []
pl_values = []
pd_values = []
pel_values = []
ppv_values = []
pcurt_values = []
powerbal_values = []


def params_xdot():
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
    q_loss= 0.1
    q_rad=14.0 #defining heat being used in the house...
    #defining the heat that the DG returns when on (being used at 70%): 
    DG_heat = 22.0 # kW heat kiloJoule/s ##########BECAUSE SAMPLING EVERY TEN MINUTES!
    DG_el = 140.0 # kW power
    BOILER_MAX_HEAT = 128.0#kW heat?
    Battery_max=240 #kW charging and discharging
    q_DG= u1*DG_heat
    q_BOILER= u2*BOILER_MAX_HEAT
    t_mix_ratio = u3*x1 + x4*(1-u3)
    t_in_tes_ratio = u4*x2 + t_mix_ratio*(1-u4) #viktig at det ikke bare er x1(1-u4)

    beta=  111.94 # 403*1000/(60*60), should be correct #4 #this can also be chaaaanged 
    #pw=55 #kw just guessing for now
    pl=90 #kw, this is what the hotel/house needs
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
                (pb/beta))
    
    #weighing of the different components of the objective function...
    c_X3=1.3

    c_co2=0.0119 #seeing what the temperatures end up with now
    c_boiler=0.001
    c_u3=0.0011
    c_u4=0.0011
    c_curt=0.012
    c_pb = 0.6
    c_powerbalance= 0.055
    #reference temperatures to ensure high enough temperature in the "house", still don't know what these bounds should be...
    x3_ref=66.0/90
    #added objective term to punish the powerbalance!!!!!!
    # Objective term -> uttrykk for cost-funksjon
    L= u1**2*c_co2  + c_X3*(x3-x3_ref)**2 + c_boiler*u2**2 + c_powerbalance*Powerbalace**2 + c_pb*u6**2 + c_curt*u5**2 + c_u3*u3**2 + c_u4*u4**2
    return xdot,L,Powerbalace

def rk4():
     # Fixed step Runge-Kutta 4 integrator
   x_dot=params_xdot()[0]
   L=params_xdot()[1]
   Powerbalace=params_xdot()[2]
   M = 4 # RK4 steps per interval
   DT = T/N/M #dette er time-step 
   f = Function('f', [x, u, ppv], [x_dot,L ]) #f(xk,uk), dette er funksjoin man vil integrere, ender opp med x_dot og L
   X0 = MX.sym('X0', 5) #init state,
   U = MX.sym('U', 6) #sier her at det er fire u-er!!!
   PPV = MX.sym('PPV',1)
   X = X0
   Q = 0
   for j in range(M): #mer nøyaktig versjon, er runge kutta 4, de ulike k-likningene osv, her finner vi neste state
       k1, k1_q = f(X, U, PPV)
       k2, k2_q = f(X + DT/2 * k1, U, PPV)
       k3, k3_q = f(X + DT/2 * k2, U, PPV)
       k4, k4_q = f(X + DT * k3, U, PPV)
       X=X+DT/6*(k1 +2*k2 +2*k3 +k4)
       Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q)
   F = Function('F', [X0, U, PPV], [X, Q],['x0','p', 'ppv'],['xf','qf']) #xf: x final, Qf : final cost, F er integralet av f
   return F

def plant(x0,u_current, PPV,F): #take in F intead of calling on the rk4() all the time
   F_new=F(x0=x0,p=u_current,ppv=PPV) #should not have to do the brackets for p
   return F_new["xf"] #this is to not struggle with the wrong format and only get the state values

# Formulate discrete time dynamics|
if False:
   # CVODES from the SUNDIALS suite
   dae = {'x':x, 'p':u, 'ode':xdot, 'quad':L}
   F = integrator('F', 'cvodes', dae, 0, T/N)
else:
  F=rk4()

# Evaluate at a test point
#Fk = F(x0=[62.0/90, 62.0/90, 62.0/90, 62.0/90 ],p=[0.2,0.2,0.2,0.2]) #tror dette er startpunkt men må sjekke ut?
#Fk = F(x0=[0.2,0.3,],p=0.4)
#print(Fk['xf'])
#print(Fk['qf'])
w=[]
wPPV = []
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

#init conditions!!!!!!!!!!!!!!
x0_init=[62.0/90,62.0/90,62.0/90,62.0/90, 62.0/90]

lbw +=x0_init
ubw += x0_init
w0 += x0_init #her begynner casaadi å søke, må være feasible!!!!, ikke del på null

# Formulate the NLP
def formulating(w, wPPV,w0,lbw,ubw,J,g,lbg,ubg,Xk):
    # Start with an empty NLP, initialiserer et tom nlp for å bruke multiple shooting
    
    for k in range(N):
        # New NLP variable for the control
        Uk = MX.sym('U_' + str(k), 6) # have to put ,4 in this formulation????????????????????????????????????????????????????
        w   += [Uk]
        lbw += [0,0,0,0,0,-1] #dette er grensene for u (here it is taken into account that there ar 4 u's)
        ubw += [1,1,1,1,1,1] #w er decision variable, xuuuxuxuxu #trying to see if u gets bigger now
    
        w0  += [0,0,0,0,0,0]
        PPVk= MX.sym('PPV_' + str(k),1)
        wPPV += [PPVk]
        # Integrate till the end of the interval
        Fk = F(x0=Xk, p=Uk,ppv=PPVk ) #x-en på slutt av første intervalll
        Xk_end = Fk['xf']
        J=J+Fk['qf'] #inkrementerer cost

        # New NLP variable for state at end of interval
        Xk = MX.sym('X_' + str(k+1), 5) 
        w   += [Xk]
    
        #changed now to have more realistic limits, let's see!!
        lbw += [40.0/90, 40.0/90, 40.0/90, 30.0/90, 20.0/90 ] # temperatur av TES skal eeeeeeeeegt ikke gå lavere enn 65 men tester dette....
        ubw += [90.0/90, 90.0/90, 90.0/90, 90.0/90, 90.0/90 ]
        w0 += [62.0/90,62.0/90,62.0/90,62.0/90, 62.0/90]  
  
        # Add equality constraint
        g   += [Xk_end-Xk] #blå minus rød fra video, multiple shoot constrainten!!! bruker g for vanlige constraints også
        #både equality og inequality constraints havner her, om ubegrenset: upper bound uendelig for eksempel
        lbg += [0, 0, 0, 0, 0]
        ubg += [0, 0, 0, 0, 0] #changed this from [0,0,0] and it now evaluates objective function more times...
    return w,wPPV,g,w0,lbw,ubw,lbg,ubg,J,Xk_end,Fk

# Create an NLP solver
w,wPPV,g,w0,lbw,ubw,lbg,ubg,J,Xk_end,Fk=formulating(w,wPPV, w0,lbw,ubw,J,g,lbg,ubg,Xk)
#print("here is w0 yasss: ",w0, "here is the length: ",len(w0))
def update_wo(u_guess_0, x0_init):
    w0=[]
    w0 += x0_init.tolist()
    for k in range(N):
        w0 += u_guess_0.tolist() 
        w0 += x0_init.tolist()  # CONVERTING BOTH TO LIST FOR IT TO BECOME ONEBIG ARRAY...
    return w0
#print("Here is w0 before the for-loop: ", len(w0), w0)
def update_lbw_ubw(x0_init):
    lbw=[]
    ubw=[]
    lbw+=x0_init.tolist()
    ubw+=x0_init.tolist()
    for k in range(N):
        lbw += [0,0,0,0,0,-1] #dette er grensene for u (here it is taken into account that there ar 4 u's)
        ubw += [1,1,1,1,1,1]
        lbw += [40.0/90, 40.0/90, 40.0/90, 30.0/90, 20.0/90 ] # temperatur av TES skal eeeeeeeeegt ikke gå lavere enn 65 men tester dette....
        ubw += [90.0/90, 90.0/90, 90.0/90, 90.0/90, 90.0/90 ]
    return lbw,ubw


prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g), 'p': vertcat(*wPPV)} #kom tilbake til parametere som varierer, om de inngår i difflikningene
solver = nlpsol('solver', 'ipopt', prob);

#initializing the final output we want to plot...
control_results = []
final_state_results = []

#initializing the first u_guesses and x_guesses
num_controls = 6
num_states = 5
u_guess = np.array(w0[: num_controls * N]).reshape(N, num_controls).T
x_guess = np.array(w0[num_controls * N :]).reshape(N + 1, num_states).T

# Solve the NLP
for i in range(N):
    ppv=np_D_ppv[i] #does this even help???
    mpc_ppv = np_D_ppv[i:i+N]
    final_state_results.append(x0_init) 
    #solver skal ha lengdte på 1156 inn som w0,lbw og ubw, mens lbg og ubg skal ha leNgde 576
    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg, p=mpc_ppv) ###################
    w_opt = sol['x'].full().flatten()
    #print("w_opt før u_guess: ", len(w_opt),w_opt[0:29])
    u_guess = w_opt[: num_controls * N].reshape(N, num_controls).T 
    x_guess = w_opt[num_controls * N :].reshape(N + 1, num_states).T
    #print(w_opt[4:8])
    #control_results.append(w_opt[4:8]) 
    control_results.append(w_opt[5:11] )# Need to check this out: 6 u's so think this is the first control input to implement!
    #THIS THE MISSING PART OF PPV?
    plantd_next=plant(x0_init,w_opt[5:11], np_D_ppv[i],F)
    #THIS THE MISSING PART OF PPV?
    x0_init=np.array(plantd_next).flatten() #to get x0_init on the right format as a [ 1 2 3]
    u_guess_0=w_opt[5:11]
    w0=update_wo(u_guess_0,x0_init)
    lbw,ubw= update_lbw_ubw(x0_init)
    
  
########################################################################################

final_state_results_done= [item for sublist in final_state_results for item in sublist]


x1_opt = final_state_results_done[0::5]
x2_opt = final_state_results_done[1::5]
x3_opt=final_state_results_done[2::5]
x4_opt=final_state_results_done[3::5]
x5_opt=final_state_results_done[4::5]


control_results_done = [item for sublist in control_results for item in sublist]
u1_opt = control_results_done[0::6]
u2_opt = control_results_done[1::6]
u3_opt= control_results_done[2::6]
u4_opt = control_results_done[3::6]
u5_opt = control_results_done[4::6]
u6_opt = control_results_done[5::6]

for i in range(N):
    BOILER_MAX_HEAT= 128.0#kW heat
    DG_el = 140.0 #kW
    pl=90 #kW
    ppv=np_D_ppv[i]
    Battery_max=240 #kW
    pb_values.append(u6_opt[i]*Battery_max)
    pel_values.append((BOILER_MAX_HEAT*u2_opt[i])/0.98)
    pd_values.append(u1_opt[i]*DG_el)
    #ppv_values.append(ppv)
    pcurt_values.append(ppv*u5_opt[i])
    pl_values.append(pl)
    powerbal_values.append(u6_opt[i]*Battery_max+ ppv+u1_opt[i]*DG_el-(BOILER_MAX_HEAT*u2_opt[i])/0.98 - pl -ppv*u5_opt[i])

np_pb=np.array(pb_values)
#print("Here the final result of pb_values is: ", np_pb)
np_pl=np.array(pl_values)
np_pd=np.array(pd_values)
np_pel=np.array(pel_values)
#np_ppv=np.array(ppv_values)
np_pcurt=np.array(pcurt_values)
np_powerbal=np.array(powerbal_values)
#print('Here are the values of the powerbalance: ', np_powerbal)



x1_done_np = np.array(x1_opt)
x2_done_np = np.array(x2_opt)
x3_done_np = np.array(x3_opt)
x4_done_np = np.array(x4_opt)
x5_done_np = np.array(x5_opt)

u1_done_np = np.array(u1_opt)
u2_done_np = np.array(u2_opt)
u3_done_np = np.array(u3_opt)
u4_done_np = np.array(u4_opt)
u5_done_np = np.array(u5_opt)
u6_done_np = np.array(u6_opt)
#calculating the KPI::
DG_heat = 22.0 # kW
KPI=0
for i in range(N):
    KPI+= u1_done_np[i]

KPI=KPI*DG_heat*0.57
print('Here is the KPI for this test: ', KPI)

#tgrid = [T/N*k for k in range(N+1)]
tgrid = [(T/N*k)/(60*60) for k in range(N)]
#tgrid = [(T/N*k)/(60*60) for k in range(N+1)]

#t_values = np.linspace(0, T, N+1)  # Adjusted initialization
t_values = np.linspace(0, T, N)  # Adjusted initialization
t_vals_pbosv = np.linspace(0, T, N)
#print("here is t_values: ",t_values)

import matplotlib.pyplot as plt
plt.figure(1)
plt.clf()
plt.plot(tgrid, x1_done_np*90, '--')
plt.plot(tgrid, x2_done_np*90, '-')
plt.plot(tgrid,x3_done_np*90, '.')
plt.plot(tgrid,x4_done_np*90, '-')
plt.plot(tgrid,x5_done_np*90, '-')
plt.xlabel('T: hours')
plt.ylabel('Water Temperature degrees celcius/battery percentage')
plt.legend(['x1:temp water dg','x2:temp water boiler','x3:temp water tes','x4:temp water ahouse', 'x5: Battery percentage'])

# Plotting
plt.figure(2)
plt.plot(tgrid, u1_done_np, '--') 
plt.plot(tgrid, u2_done_np, '-.') 
plt.plot(tgrid, u3_done_np, '--') 
plt.plot(tgrid, u4_done_np, '-.') 
plt.plot(tgrid, u5_done_np*1000, '-.')
plt.plot(tgrid, u6_done_np*300, '-.')
plt.xlabel('T: hours')
plt.ylabel('Percentage')
plt.legend(['u1:power_%_DG','u2:power_%_boiler','u3:%_mass_flow_DG','u4:%_mass_flow_boiler', 'u5: ppv curtailed','u6: pb'])

plt.figure(3)
#plotting the new state and scaling mehhhh
plt.plot(t_vals_pbosv/(60*60), np_pb*300, '--' )
plt.plot(t_vals_pbosv/(60*60),np_pcurt*1000, '--')
plt.plot(t_vals_pbosv/(60*60),np_pel, '-')
plt.plot(t_vals_pbosv/(60*60),np_pd,'.-')
plt.plot(t_vals_pbosv/(60*60),np_pl,'.-' )
plt.plot(t_vals_pbosv/(60*60), np_D_ppv_plott, '-.')
plt.plot(t_vals_pbosv/(60*60), np_powerbal, '-.')
plt.ylabel('Power in kW')
plt.xlabel('T: hours')
plt.legend(['pb: power in/out of batttery','pcurt: curtailed PV power','pel: power used electrical boiler'
            , 'pd:power prod Diesel', 'pl: power used Isfjorden','ppv: power prod PV', 'powerbalance'])


plt.grid()
plt.show()


