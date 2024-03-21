import numpy as np
from casadi import *
import matplotlib.pyplot as plt

x1 = MX.sym('x1')
x2 = MX.sym('x2')
x3 = MX.sym('x3')
x4 = MX.sym('x4')
x = vertcat(x1, x2, x3, x4)
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
    return 14.0

q_rad = q_sorad()

# defining the heat that the DG returns when on (being used at 70%):
DG_heat = 22.0 # kW heat kiloJoule/s ##########BECAUSE SAMPLING EVERY TEN MINUTES!
DG_el = 140.0 # kW power
B_MAX_HEAT = 128.0 #kW heat?
#now the u's are constant!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
u1 = 0.09
u2 = 0.097
u3 = 0.4
u4 = 0.4
q_DG = u1 * DG_heat
q_BOILER = u2 * B_MAX_HEAT
t_mix_ratio = u3 * x1 + x4 * (1 - u3)
t_in_tes_ratio = u4 * x2 + t_mix_ratio * (1 - u4)

# Defining xdot
xdot = vertcat((cp * u3 * w_tot * (x4 - x1)*90 + q_DG - q_loss) / (rho * V_dg * c_dg_heatcap*90),
               (cp * u4 * w_tot * (t_mix_ratio - x2)*90 - q_loss + q_BOILER) / (rho * V_boiler * c_boiler_heatcap*90),
               (cp * w_tot * (t_in_tes_ratio - x3)*90 - q_loss) / (rho * V_tes * c_tes*90),
               (cp * w_tot * (x3 - x4)*90 - q_rad - q_loss) / (rho * V_rad * c_rad_heatcap*90))

# Define initial conditions
x0_init = [70.0/90, 68.0/90, 62.0/90, 62.0/90]

# Time parameters
T = 24.0*60*60 # Total time
N = 144  # Number of time steps
dt = T / N  # Time step size

# Initializing lists to store states and time
x_values = []
t_values = np.linspace(0, T, N+1)  # Adjusted initialization

# Initializing the initial state
x_current = x0_init

# Defining the system dynamics as a CasADi function
xdot_function = Function('xdot_function', [x], [xdot])


# Simulating the system dynamics over time
#using euler, this is without a controller, just to get out the states....
for i in range(N+1):  # Adjusted range to ensure N+1 time steps
    # Storing the current state
    x_values.append(x_current)

    # Computing the next state using the system dynamics function
    xdot_val = xdot_function(x_current)
    x_next = np.array(x_current) + dt * np.array(xdot_val.full().flatten())  # Flatten the xdot_val

    # Updating the current state for the next iteration
    x_current = x_next

# Converting the state values to arrays for plotting
x_values = np.array(x_values)*90

# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(t_values, x_values[:, 0], label='x1: Temperature in DG Water')
plt.plot(t_values, x_values[:, 1], label='x2: Temperature in Boiler Water')
plt.plot(t_values, x_values[:, 2], label='x3: Temperature in TES Water')
plt.plot(t_values, x_values[:, 3], label='x4: Temperature in Water after house')
plt.xlabel('Time')
plt.ylabel('Temperature (°C)')
plt.title('Development of Temperatures Over Time')
plt.legend()
plt.grid(True)
plt.show()

