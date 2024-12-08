import scipy as sc
import matplotlib . pyplot as plt
import numpy as np

#gamma function found in dρ/dr
def gamma(y):   
    x = y**(2/3)
    return x/(3*(1+x)**0.5)

def f(t,y):#function t=radius y[0]=mass y[1]=density

    #dm/dr=3*r^2*ρ
    f1 = 4*np.pi*y[1]*(t**2)
    
    #dρ/dr=-Cc*ρ*m/r^2*gamma(ρ/ρ0)
    f2 = -(constant*y[0]*y[1])/(gamma(y[1]/rho_0)*(t**2))
    return [f1,f2]

def hit_edge(t,y): return y[1]-0.01 #if density hits 0.01 then trigger event
hit_edge.terminal = True #makes the event terminates solve_ivp
hit_edge.direction = -1 #terminates if it goes from positive to negative


#t0 is the starting radius for interation
t0 = 0.00001

#tf is the final radius set to 1
tf = 1

#t is an array used for t_eval
t = np.linspace(t0,tf,20000)

#Ye is number of electrons per nucleon
#ye is an array filled with ye value for carbon then iron
Ye = [0.5,0.46428]

#for labelling and colouring the graphs
name=['carbon','iron']
colour=['blue','green']
info=[]

#=============================mass plots=================================

#to loop over for carbon and for iron
for j in range(2):
    
    #rho_0 = Mp*Me^3*c^3/3*pi^2*hbar^3*Ye
    #rho_0 = 9.79e8/Ye kg/m^3
    #rho_0 = 165986/Ye Mo/Ro^3
    rho_0 = 165986/Ye[j] 
    
    #our differerntial equation dρ/dr has a constant in it.
    #constant = G*Mp/Ye*Me*c^2
    constant = ((6.67e-11)*(1.67e-27))/((Ye[j])*(9.11e-31)*(9e16)) #m/kg
    constant = constant*(1.98e30/6.95e8) #solar units

    #some arrays to place the ending mass and radii for different densities
    X = []
    Y = []
    
    #a range of small densities
    y0_set = np.linspace(1e-4,2e5,300)
    
    #loop over the differernt densities
    for i in range(1,len(y0_set)):
        
        #solves the differerntial equation for the given density
        sol = sc.integrate.solve_ivp(f,[t0,tf],[0,y0_set[i]],t_eval=t,events = hit_edge)

        #appends the final radius component to a list
        X.append(sol.t[-1])
        
        #appends the final mass component to a list
        Y.append(sol.y[0][-1])
        
        #this checks if it is carbon
        if j == 0:
            
            #if so this will plot the individual lines for each density
            plt.plot(sol.t,sol.y[0],c='r')

    #repeat same process for medium densities    
    y0_set = np.linspace(1e5,2e8,300)
    for i in range(1,len(y0_set)):
        sol = sc.integrate.solve_ivp(f,[t0,tf],[0,y0_set[i]],t_eval=t,events = hit_edge)
        X.append(sol.t[-1])
        Y.append(sol.y[0][-1])
        if j == 0:
            plt.plot(sol.t,sol.y[0],c='b')
    
    #repeat same process for high densities
    y0_set = np.linspace(1e8,2e10,300)
    for i in range(1,len(y0_set)):
        sol = sc.integrate.solve_ivp(f,[t0,tf],[0,y0_set[i]],t_eval=t,events = hit_edge)
        X.append(sol.t[-1])
        Y.append(sol.y[0][-1])
        if j == 0:
            plt.plot(sol.t,sol.y[0],c='g')
            
    #this just shows that first plot
    if j == 0:
        plt.title('varying density from 1e-5 - 2e10')
        plt.ylabel('Mass / Mo') 
        plt.xlabel('radius / Ro')
        plt.axhline(0,c='k')
        plt.axvline(0,c='k')
        # plt.grid()
        plt.show()

    #plots the mass radius lines for both carbon and iron
    plt.plot(Y, X, label = name[j], c = colour[j], zorder = 1)
    info.append(f'for {name[j]} core white dwarfs the mass limit is {sol.y[0][-1]} Mo')
    

#plot on real world data to compare.
mass = [1.053,0.48,0.50]
mass_error = [0.028,0.02,0.05]
radius = [0.0074,0.0124,0.0115]
radius_error = [0.0006,0.0005,0.0012]
c = ['red', 'darkviolet', 'dodgerblue']
name = ['Sirus B', '40 Eri B', 'Stein 2051']
for k in range(len(mass)): 
       plt.errorbar(mass[k],radius[k],radius_error[k],mass_error[k],linestyle = "", capsize = 4, marker = 'o' ,markerEdgeColor = 'k', c = c[k], label = name[k], zorder = 2)

#SDSS DATA
SDSS_mass = [0.48,0.84,0.54,0.54,0.80,0.56,0.54,0.52,0.98,0.49]
SDSS_radius= [0.0154, 0.0099, 0.0143, 0.0142, 0.0105, 0.0139, 0.0143, 0.0146, 0.0083, 0.0152]
SDSS_mass_error = [0.06,0.14,0.06,0.04,0.03,0.06,0.04,0.04,0.07,0.03]
SDSS_radius_error = [0.0012, 0.0018, 0.0011, 0.0008, 0.0003, 0.001, 0.0008, 0.0007, 0.0008, 0.0007]
plt.errorbar(SDSS_mass,SDSS_radius,SDSS_radius_error,SDSS_mass_error,linestyle = "", capsize = 3, marker = 'o',markerEdgeColor = 'navy', c = 'k', label = 'SDSS set', zorder = 2)
#sets the plot to look nice
plt.title('Mass Profile Plot')
plt.xlabel('Mass (Mo)')
plt.ylabel('radius (Ro)')
plt.axhline(0,c='k')
plt.axvline(0,c='k')
# plt.grid()
plt.legend()
plt.show()

#=============================density plots=========================
#now for the density radius plot

name=['carbon','iron']
#again loop of both Ye for carbon and iron
for j in range(2):
 
    #same constants as before
    rho_0 = 165986/Ye[j] 
    constant = ((6.67e-11)*(1.67e-27))/((Ye[j])*(9.11e-31)*(9e16))#m/kg
    constant = constant*(1.98e30/6.95e8)

    #solve and plot the equations
    sol = sc.integrate.solve_ivp(f,[t0,tf],[0,rho_0],t_eval=t,events = hit_edge)
    info.append(f'for {name[j]} core white dwarfs the maximum radius size is {sol.t[-1]} Ro')
    plt.plot(sol.t,sol.y[1]/rho_0,label = name[j])

#plot details
plt.title('density profile plot')
plt.xlabel('radius (Ro)')
plt.ylabel('ρ / ρo')
plt.axhline(0,c='k')
plt.axvline(0,c='k')
plt.legend()
plt.show()

print()
print('mass')
print(info[0])
print(info[1])
print()
print('density')
print(info[2])
print(info[3])
'''
quantifying errors:
curve_fit = chi fit minimization
chi squared over degrees of freedom,good fit = 1
confidence level ~ sigma / standard-deviation
degrees if freedom is number of points you have - number of parameter you are fitting
'''

'''
#SDSS DATA
SDSS_mass = [0.48,0.84,0.54,0.54,0.80,0.56,0.54,0.52,0.98,0.49]
SDSS_radius = [1.54,0.99,1.43,1.42,1.05,1.39,1.43,1.46,0.83,1.52]
SDSS_mass_error = [0.06,0.14,0.06,0.04,0.03,0.06,0.04,0.04,0.07,0.03]
SDSS_radius_error = [0.12,0.18,0.11,0.08,0.03,0.10,0.08,0.07,0.08,0.07]
'''