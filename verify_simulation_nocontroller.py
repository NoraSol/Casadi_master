#Parameters in my system
cp=4.186
rho=1.0
V_dg = 3000
V_boiler = 3000                                                                                                                                                                                                                                     
V_rad=3000
V_tes=12000 #liter
w_tot=500 #denne må endres til å kunne variere med tanke på pumpen, dette er jo totale gjennomstrømningen i systemet...
#in the last system several of the power-funtions were functions of t, how to do this descrete?
q_loss= 0.1

def q_sorad(): 
    return 17000.0
q_rad=q_sorad() #defining heat being used in the house...

#defining the heat that the DG returns when on (being used at 70%): 
DG_heat = 22000.0 # kW heat
DG_el = 140000.0 # kW power
B_MAX_HEAT = 128000.0 #kW heat?