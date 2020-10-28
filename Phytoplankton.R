############################################
# Phytoplankton mathematical modeling 2020 #
############################################

# Chosen Values
u0  = 1000                # Initial amount of phytoplankton
vp  = 0.9                 # Velocity of phytoplancton downwards
k   = 4                   # Diffusion constant
dt  = 0.5                 # Something we chose
dz  = 0.5                 # Something we chose                                        

# Dichlet Boundaries 
L = 100                  # Since phytoplankton is not found in depths deeper then 100m
t = seq(0, L, by = dt)   # Here L is used since 100 is a nice number for the time
N = L/dz                 # Number of special intervals

# initial calculation of n as Gauss distributed phytoplankton
# the density of phytoplankton is decreasing the deeper they are found.
# The calculation should look like this:

z = seq(0, L, by = dz)   # Here L is used as max due to death in lower sea levels
u = seq(u0,0,,N+1)       # initialize the Population u as 10000 phytoplankton to 0 over time
u_next = u               # Save for later use  

#--------------#
#Equation parts:
#--------------#

# The differential equation has multiple parts:
# 1. The velocitiy of the phytoplanktopn
A    =  (vp*dt) / dz
# 2. The difusion rate of the Phytoplankton
B    =  (k*dt) / (dz^2)
# 3. Finally the f the function for the which takes in account the amount of phytoplankton
#    here only the values of +1 or -1 is used
S    =  seq(100,10,, N+1)
ffac =  ifelse(cumsum(S)*dz > 100, -1, 1)
f    =  ffac / 100

#----------------#
# Plot of the PDE:
#----------------#
#Plot of initial conditions 

plot (t,u)

#change over time 
for (i in 2:N){ 
  
  if (i == 1) {
    u_next[i] = (f[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u[i]))*dt
    plot(t,u_next,type = "l", xlab = 'Time', ylab = 'Population of Phytoplankton')
  }
  Sys.sleep(0.05)
  u_next[i] = (f[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u[i-1]))*dt
  
  Sys.sleep(0.05)
  plot(t,u_next,type = "l", xlab = 'Time', ylab = 'Population of Phytoplankton')
  
  
}

u_next            


# Additional bondnary condition needs to be added from the second function

#Questions to ask
# we expect exponentail decay
# the implementation does not work as intended -> why does u_next just calculate 2 different values?


# Next Things to do:
# Speak with Alish and inform him about what matthias told me
# try to get this project done soonish since we do not have much time left :-S