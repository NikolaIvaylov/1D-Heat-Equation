import numpy as np
from matplotlib import pyplot as plt, animation
import math

"--- 1D Heat Equation Simunlation using the Explicit Method (both 2D and 3D plots) ---"

'''--- Variables you can play with ---'''
#Rod Length:
L = 10 #total lenght of the rod (meters)
length_step = 0.5 #length step (meters), determines the time step. Should be less than 1
#Other Variables:
k = 0.2 #thermal diffusivity constant, should be less than 1
totally_insulated = True #can be False. If it's True, the border temperature won't be taken into account
left_border_temp = 10 #the temperature at the left border (°C)
right_border_temp = 10 #the temperature at the right border (°C)
#Time:
T = 120 #for how many seconds the simulation should run?
#Initial Temperature:
'''--- The initial temperature is a sinusoid, with:
   >>> Amplitude 'A'
   >>> Frequency 'w' ---'''
A = 50 #°C
w = 1 #Hz

'''--- Variables you cannot play with ---'''
#Time:
time_step = ((length_step) ** 2) / (2 * k) #the time step is calculated here so the scheme is stable
time = np.arange(0, T, time_step) #list of all values of t
#Rod Lenght:
length = np.arange(0, L, length_step) #list of all values of l
#Initial Temperature Values:
initial_temperatures = []
for n in range(0,len(length)):
    temp_at_l = A * math.sin(w*length[n])
    initial_temperatures.append(temp_at_l)
#list of lists to store all values of all temperatures in time:
all_temperatures = [initial_temperatures] #initial temperatures are included by default

#case of the ends NOT being isnulated:
if totally_insulated == False:
    #nested loop to apply Explicit Method formulas on all values of the length for each value of time:
    for t in time:
        next_temperatures = [] #list to save newly calculated values
        for n in range(0, len(length)):
            prev_temperatures = all_temperatures[-1] #get all values of the previous temperatures
            #series of if() statements to check for Border Conditions (BC):
            if n == 0:
                new_temp = prev_temperatures[n] + (k * ((left_border_temp + prev_temperatures[n+1]) / 2 - prev_temperatures[n]) * time_step)
            elif n == (len(length)-1):
                new_temp = prev_temperatures[n] + (k * ((prev_temperatures[n-1] + right_border_temp) / 2 - prev_temperatures[n]) * time_step)
            else:
                new_temp = prev_temperatures[n] + (k * ((prev_temperatures[n-1] + prev_temperatures[n+1]) / 2 - prev_temperatures[n]) * time_step)
            next_temperatures.append(new_temp) #add the temperature calculated for length[n] for time[t]
        all_temperatures.append(next_temperatures) #add all temperatures calculated for the whole length for time[t]
#case of the ends being isnulated:
else:
    #nested loop to apply Explicit Method formulas on all values of the length for each value of time:
    for t in time:
        next_temperatures = [] #list to save newly calculated values
        for n in range(0, len(length)):
            prev_temperatures = all_temperatures[-1] #get all values of the previous temperatures
            #series of if() statements to check for Border Conditions (BC):
            if n == 0:
                new_temp = prev_temperatures[n] + (k * ((prev_temperatures[n+1]) - prev_temperatures[n]) * time_step)
            elif n == (len(length)-1):
                new_temp = prev_temperatures[n] + (k * ((prev_temperatures[n-1]) - prev_temperatures[n]) * time_step)
            else:
                new_temp = prev_temperatures[n] + (k * ((prev_temperatures[n-1] + prev_temperatures[n+1]) / 2 - prev_temperatures[n]) * time_step)
            next_temperatures.append(new_temp) #add the temperature calculated for length[n] for time[t]
        all_temperatures.append(next_temperatures) #add all temperatures calculated for the whole length for time[t]
    
#output logic:
def output():
    #options of the 2D plot:
    plt.rcParams["figure.figsize"] = [L, (abs(A/10))]
    plt.rcParams["figure.autolayout"] = True
    plt.rcParams['figure.dpi'] = 300
    fig, ax = plt.subplots()
    ax.set(xlim=(0, L), ylim=(-A, A))
    line, = ax.plot(length, all_temperatures[0], color='k', lw=2)
    #display all temperatures in a 2D plot:
    def animate(i):
        line.set_ydata(all_temperatures[i])
    anim = animation.FuncAnimation(fig, animate, interval=100, frames=len(time) - 1)
    anim.save(f'Insulated_{totally_insulated}, (Δx = {length_step}, k = {k}).gif') #save animation as a GIF file
    plt.show()
    
    #options of the 3D plot:
    plt.figure(dpi=350)
    plt.axes(projection='3d')
    plt.title(f"Heat Equation Simulation (Δx = {length_step}, k = {k})")
    plt.xlabel("Time (s)")
    plt.ylabel("Lenght (m)")
    #display all temperatures in a 3D plot:
    for i in range(0, len(time)):
        plt.plot([time[i]]*len(length), length, all_temperatures[i])
    plt.show()
    
output()