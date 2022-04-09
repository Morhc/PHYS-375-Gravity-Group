# Main code for generating the main sequence stars
# NO GRAVITY MODIFICATION! Let's copy this working code into a new file solely for the gravity modification

import os
import numpy as np
import pandas as pd
from standards import Standards as s
import tools as tools
import scipy.integrate as int
import matplotlib.pyplot as plt
#import star_info as star_info

np.seterr(all='ignore')

def solveODEs(r_init, r_fin, y0, steps):
    '''Function to solve the ODEs (equation 2 in the project description)'''

    r_values = np.linspace(r_init, r_fin, steps)
    t_span = (r_init, r_fin)

    solutions = int.solve_ivp(tools.dydr, t_span, y0, 'RK45', t_eval=r_values) # outputs an array 'y' which contain the solutions to the ODEs which are described dydr from tools.py

    rho = solutions.y[0]
    T = solutions.y[1]
    M = solutions.y[2]
    L = solutions.y[3]
    tau = solutions.y[4]
    r = solutions.t # t here mean the 'time points' which in this case are the radius values or rvals

    return r, rho, T, M, L, tau

def bisection(rhoc_min, rhoc_max, r_initial, r_final, y, steps):
    '''This will use the bisection method to vary and determine the value for
    for rho_c such that the lumnosity boundary condition is  met
    '''

    rhoc_mid = (rhoc_min + rhoc_max)/2 # determine mid point for rhoc given the min and max inputs

    y_min = np.copy(y) # new intial value array for rhoc_min
    y_min[0] = rhoc_min

    y_mid = np.copy(y) # new intial value array for rhoc_mid
    y_mid[0] = rhoc_mid

    y_max = np.copy(y) # new intial value array for rhoc_max
    y_max[0] = rhoc_max

    # Now we want to use rhoc min,mid, and max to solve the ODEs individually
    solutions_min = solveODEs(r_initial, r_final, y_min, steps)
    solutions_mid = solveODEs(r_initial, r_final, y_mid, steps)
    solutions_max = solveODEs(r_initial, r_final, y_max, steps)

    # Need to find R_star, L_star, and T_star for the three solutions
    tau_min = solutions_min[5]
    tau_mid = solutions_mid[5]
    tau_max = solutions_max[5]

    i_min = np.argmin(tau_min[-1] - tau_min - 2/3)
    i_mid = np.argmin(tau_mid[-1] - tau_mid - 2/3)
    i_max = np.argmin(tau_max[-1] - tau_max - 2/3)

    R_Star_min = solutions_min[0][i_min]
    R_Star_mid = solutions_mid[0][i_mid]
    R_Star_max = solutions_max[0][i_max]

    L_Star_min = solutions_min[4][i_min]
    L_Star_mid = solutions_mid[4][i_mid]
    L_Star_max = solutions_max[4][i_max]

    T_Star_min = solutions_min[2][i_min]
    T_Star_mid = solutions_mid[2][i_mid]
    T_Star_max = solutions_max[2][i_max]

    # Using the luminosity_check (equation 17) for the valoues found using rhoc min,max, mid
    fmin = tools.luminosity_check(R_Star_min, L_Star_min, T_Star_min)
    fmid = tools.luminosity_check(R_Star_mid, L_Star_mid, T_Star_mid)
    fmax = tools.luminosity_check(R_Star_max, L_Star_max, T_Star_max)

    # Note: if 'f' is > 0 then our L(R_star) is greater than the theoretical

    if fmin*fmax > 0:
        print ("Change initial conditions, no rho_c possible.")
        quit()
    else:
        #initial condition for rhoc
        current = rhoc_mid
        #arbitrary initializing
        previous, diff = 0, 2

        while abs(diff) > 1:
            if fmin*fmid < 0: # Either fmin or fmid overshot while the other undershot
                rhoc_max = rhoc_mid

                rhoc_mid = (rhoc_min+rhoc_max)/2

                # Now we re-solve ODEs using the new rhoc_mid
                y_mid = np.copy(y) # new intial value array for rhoc_mid
                y_mid[0] = rhoc_mid

                soln = solveODEs(r_initial, r_final, y_mid, steps)
                tau = soln[5]
                i = np.argmin(tau[-1] - tau - 2/3)


                R_Star_mid = soln[0][i]
                T_Star_mid = soln[2][i]
                L_Star_mid = soln[4][i]

                fmid = tools.luminosity_check(R_Star_mid, L_Star_mid, T_Star_mid)

            elif fmid*fmax < 0: # Either fmid or fmax overshot while the other undershot
                rhoc_min = rhoc_mid
                rhoc_mid = (rhoc_min+rhoc_max)/2

                # Now we re-solve ODEs using the new rhoc_mid
                y_mid =  np.copy(y) # new intial value array for rhoc_mid
                y_mid[0] = rhoc_mid

                soln = solveODEs(r_initial, r_final, y_mid, steps)
                tau = soln[5]
                i = np.argmin(tau[-1] - tau - 2/3)


                R_Star_mid = soln[0][i]
                T_Star_mid = soln[2][i]
                L_Star_mid = soln[4][i]

                fmid = tools.luminosity_check(R_Star_mid, L_Star_mid, T_Star_mid)

            elif fmid == 0:
                return rhoc_mid

            previous = current
            current = rhoc_mid
            diff = current - previous

    return rhoc_mid

###############################################################################
######################### Setting the Initial conditions ######################
###############################################################################
r_initial = 1
r_final = 1e9
steps = 10000 #number of steps

Tc = 8.23544e+06  # (core temperature in K) we will choose Tc to be the MS paramter

# intializing an array that contains the intial values for density, temperature, mass, luminosity and tau
y = np.array([0, Tc, 0, 0, 0])

rhoc_min = 300 # min value for rho_c in g/cm^3
rhoc_max = 500000 # max value for rho_c in g/cm^3
#our solution to f(rho_c)
rho_c = bisection(rhoc_min, rhoc_max, r_initial, r_final, y, steps)

# re-defining the previos 'y' array to now contain the newly determine rho_c value
y_new =  np.array([rho_c, Tc, 0, 0, 0])

# Solving the ODEs
solutions_final = solveODEs(r_initial, r_final, y_new, steps)

# Storing the solutions to the ODEs appropriately
r_values = solutions_final[0]
rho_values = solutions_final[1]
T_values = solutions_final[2]
M_values = solutions_final[3]
L_values = solutions_final[4]
tau_values = solutions_final[5]

i = np.argmin(tau_values[-1] - tau_values - 2/3)


# Declaring R_Star, Rho_Star, T_Star, M_Star, and L_Star
R_Star = r_values[i]
Rho_Star = rho_values[i]
M_Star = M_values[i]
L_Star = L_values[i]
T_Star = (L_Star/(4*np.pi*s.sb*(R_Star**2)))**(1/4) #T_values[index_at_surface]

# Saving all values related to pressure
P_values = tools.pressure(rho_values, T_values)
P_degen_values = tools.P_degen(rho_values)
P_ideal_values = tools.P_ideal(rho_values, T_values)
P_photongas_values = tools.P_photongas(T_values)
P_c = tools.pressure(rho_c, Tc) # calcualting the central pressure

# Saving all values related to kappa
# Note all of these kappa values are in units of cm**2/g
# Can convert from m**2/kg to cm**2/g by multiplying by 10
kappa_values = []
for i in range(0, len(rho_values)):
    kappa_values.append( tools.kappa(rho_values[i], T_values[i]) * 10)

kappa_ff_values = tools.kappa_ff(rho_values, T_values) * 10
kappa_es_values = np.full(len(kappa_ff_values),tools.kappa_es() * 10)
kappa_Hminus_values = tools.kappa_Hminus(rho_values, T_values) * 10

# Saving all values related to dL/dr
dL_dr_values = tools.dL_dr(rho_values , T_values, r_values)
dL_dr_pp_values = tools.dL_dr_PP(rho_values , T_values, r_values)
dL_dr_cno_values = tools.dL_dr_CNO(rho_values , T_values, r_values)

# Saving all values for dln(P)/dln(T) partial derivative
dlnP_dlnT_values = tools.dlnP_dlnT(rho_values , T_values)

print('Saving the data')
### SAVING OUT THE DATA ###
summary_file = os.path.join(s.data_folder, 'generated_stars.csv')

with open(summary_file, 'a+') as f:
    f.write(f'{rho_c},{Tc},{R_Star},{M_Star},{L_Star},{T_Star}\n')

df = pd.DataFrame(data=zip(r_values/R_Star, rho_values, L_values, M_values, T_values, P_ideal_values,
                           P_degen_values, P_photongas_values, P_values, dlnP_dlnT_values,
                           kappa_values, kappa_ff_values, kappa_Hminus_values, kappa_es_values,
                           dL_dr_values, dL_dr_pp_values, dL_dr_cno_values),
                  columns=['r/R', 'rho(r)', 'L(r)', 'M(r)', 'T(r)', 'Pgas', 'Pdeg', 'Ppho', 'P', 'dlogP/dlogT',
                           'kappa', 'kappaff', 'kappaHm', 'kappaes', 'dL/dr', 'dLpp/dr', 'dLcno/dr'])

df.to_csv(rf"{os.path.join(s.data_folder, 'star_1.csv')}", index=False)

quit()

print('R_star is', R_Star/s.Rsun)
print('Rho_star is', Rho_Star)
print('T_star is', T_Star)
print('M_star/M_sun is', M_Star/s.Msun)
print('L_star/L_sun is', L_Star/s.Lsun)
print('L_star is', L_Star)
print('rho_c is', rho_c)

# Plotting normalized density, L,M, and T
plt.plot(r_values/R_Star, rho_values/rho_c, 'k' ,label='rho/rho_c')
plt.plot(r_values/R_Star, L_values/L_Star, '-.b',label='L/L_*')
plt.plot(r_values/R_Star, M_values/M_Star, '-g',label='M/M_*')
plt.plot(r_values/R_Star, T_values/Tc, ':r', label='T/Tc')
plt.xlim(0,1)
plt.legend(loc='best')
plt.show()

# Plotting pressure results
plt.plot(r_values/R_Star, P_values/P_c, 'k' ,label='P')
plt.plot(r_values/R_Star, P_degen_values/P_c, '--r' ,label='P_degen')
plt.plot(r_values/R_Star, P_ideal_values/P_c, '-.b' ,label='P_ideal')
plt.plot(r_values/R_Star, P_photongas_values/P_c, ':g' ,label='P_photongas')
plt.xlim(0,1)
plt.legend(loc='best')
plt.show()

# Plotting kappa results
# In cm**2/g
plt.plot(r_values/R_Star, np.log10( kappa_values), 'k' ,label='kappa')
plt.plot(r_values/R_Star, np.log10( kappa_es_values), '--b' ,label='kappa_es')
plt.plot(r_values/R_Star, np.log10( kappa_ff_values), ':g' ,label='kappa_ff')
plt.plot(r_values/R_Star, np.log10( kappa_Hminus_values), '-.r' ,label='kappa_Hminus')
plt.xlim(0,1)
plt.legend(loc='best')
plt.show()

# Plotting all dL/dr results
plt.plot(r_values/R_Star, dL_dr_values*(R_Star/L_Star), 'k' ,label='dL/dr')
plt.plot(r_values/R_Star, dL_dr_pp_values*(R_Star/L_Star), '--r' ,label='dL/dr PP')
plt.plot(r_values/R_Star, dL_dr_cno_values*(R_Star/L_Star), '-.b' ,label='dL/dr CNO')
plt.xlim(0,1)
plt.legend(loc='best')
plt.show()

# Plotting ln(dP/dT) # NOT SURE ABOUT THIS LAST ONE!
plt.plot(r_values/R_Star, dlnP_dlnT_values, 'k' ,label='dln(P)/dln(T)')
plt.xlim(0,1)
plt.legend(loc='best')
plt.show()

'''
The arrays that hold the values of interest are:
r_values
rho_values
T_values
M_values
L_values
tau_values (we never actually plot this so probably don't need this)

kappa_values
kappa_ff_values
kappa_es_values
kappa_Hminus_values
dL_dr_values
dL_dr_pp_values
dL_dr_cno_values
dlnP_dlnT_values

The constants:
rho_c
Tc
R_Star
Rho_Star
T_Star
M_Star
L_Star
'''
# Below is my attempt at using the star_info.py and pandas
#data = pd.DataFrame(data=zip(M_Star, r_values, rho_values, T_values, M_values, L_values, dL_dr_values, dL_dr_pp_values, dL_dr_cno_values, dlnP_dlnT_values, r_values/R_Star, rho_values/rho_c), columns=['M (kg)', 'r (m)', 'rho(r) (kg/cm^3)', 'T(r) (K)', 'M(r) (kg)', 'L(r) (J/s)', 'dL/dr (J/s/m)', 'dLpp/dr (J/s/m)', 'dLcno/dr (J/s/m)', 'dlogP/dlogT (J/s/m)', 'r/R', 'rho(r)/rhoc'])
#star_info.plot_star(data, info = None, savepath="")

################################################################################
############################     DEPREACATED CODE    ###########################
################################################################################
# ------------------------
# Section 2.2.1 (Step 2)
# ------------------------
# Now we account for the two main complications which arises from the behaviour near r = 0 and the practical determination of R_star

# We want to integrate the ODEs until del_tau (defined in eqution 16) is << 1 and we introdue a mass limit. We stop when M > 10^3 solar masses




# while (np.any(del_tau) >  del_tau_threshold) and (M < MassLimit):
#
#     y5,hnew = tools.RKF45Method_adaptive(tools.dydr,r,y,h)
#
#     r_values.append(r + h)
#     y_values.append(y5)
#
#     # Channging variales around for next iteration of the 'while' loop
#     # I.e., hnew becomes h for the next iteration and y5 becomes y for the next iteration
#     h = hnew
#     y = y5
#     r += h
#
#     rho = y[0]
#     T = y[1]
#     M = y[2]
#     L = y[3]
#     tau = y[4]
#
#     del_tau = tools.del_tau(rho,T, r, M, L )
