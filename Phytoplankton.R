############################################
# Phytoplankton mathematical modeling 2020 #
############################################

# Chosen Values
u0  = 1000                # Inital amount of Phytoplankton
vp  = 0.1                 # Velocity of phytoplancton downwards
k   = 0.1                 # Diffusion constant
dt  = 0.5                 # Something we chose
dz  = 0.5                 # Something we chose                                        

# Dichlet Boundaries 
L = 100                  # Since phytoplankton is not found in depths deeper then 100m
t = seq(0, L, by = dt)   # Here L is used since 100 is a nice number for the time
N = L/dz                 # Number of special intervals

# initial calculation of n as gaus distributed Phytoplankton
# the density of Phytoplankton is decreasing the deeper they are found.
# The calculation should look like this:

z = seq(0, L, by = dz)   # Here L is used as max due to death in lower sea levels
u = seq(u0,0,,N+1)       # initialize the Population u as 10000 Phytoplankton to 0 over time
u_next = u               # Save for later use  

#--------------#
#Equation parts:
#--------------#

# The differential equation has multiple parts:
# 1. The velocitiy of the Phytoplanktopn
A    =  (vp*dt) / dz
# 2. The difusion rate of the Phytoplankton
B    =  (k*dt) / (dz^2)
# 3. Finally the f the function for the which takes in acount the amount of Phytoplankton
#    here only the values of +1 or -1 is used
S    =  seq(100,10,, N+1)
ffac =  ifelse(cumsum(S)*dz > 100, -1, 1)
f    =  ffac / 100
#----------------#
# Plot of the PDE:
#----------------#
#Plot of initlai conditions  
plot (t,u)

#change over time 
for (i in 2:N){ 
  if (i == 1) {
    u_next[i] = (f[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+0))*dt
    plot(t,u_next,type = "l", xlab = 'Time', ylab = 'Population of Phytoplankton')
  }
  Sys.sleep(0.1)
  u_next[i] = (f[i]*S-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u[i-1]))*dt
  # Initial example of the function
  Sys.sleep(0.0005)
  plot(t,u_next,type = "l", xlab = 'Time', ylab = 'Population of Phytoplankton')
}
u_next            