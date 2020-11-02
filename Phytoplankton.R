############################################
# Phytoplankton mathematical modeling 2020 #
############################################

# Chosen Values
u0  = 1000                # Initial amount of phytoplankton
vp  = 0.9                 # Velocity of phytoplankton downwards
k   = 4                   # Diffusion constant
dt  = 0.5                 # Something we chose
dz  = 0.5                 # Something we chose                                        

# Dichlet Boundaries 
L  = 100                  # Since phytoplankton is not found in depths deeper then 100m
LT = 100                  # Final time point
t  = seq(0, L, by = dt)   # Here L is used since 100 is a nice number for the time
N  = L/dz                 # Number of special intervals
NT = LT/dt                # Number of time intervals

# initial calculation of n as Gauss distributed phytoplankton
# the density of phytoplankton is decreasing the deeper they are found.
# The calculation should look like this:


z = seq(0, L, by = dz)   # Here L is used as max due to death in lower sea levels
t = seq(0, LT, by = dz)  # Here LT is used as max due to death in lower sea levels
u = seq(u0,0,,N+1)       # initialize the Population u as 10000 phytoplankton to
# 0 over time
u_next = u               # Save for later use  

#--------------#
#Equation parts:
#--------------#

# The differential equation has multiple parts:
# 1. The velocitiy of the phytoplanktopn
A    =  (vp*dt) / dz
# 2. The difusion rate of the Phytoplankton
B    =  (k*dt) / (dz^2)

# 3. part will be coosen alter in the for loop since it depends on the 
# amount of phytoplankton

#----------------#
# Plot of the PDE:
#----------------#

#Plot of initial conditions 
plot (t,u)

# 3. Finally the f the function for the which takes in account the amount of 
# phytoplankton here only the values of +1 or -1 is used.



for (j in 1:NT) {
  #change over time 
  S =  cumsum(u*dz)           # Determin new amount of accumulated Phytoplankton
  f =  ifelse(S > 5000, -1, 1)# Determine if the amount of Phytoplankton grows 
  # or decays
  for (i in 1:N){
    if (i == 1) {
      Sys.sleep(0.1)
      u_gv= u_next[2]-((vp*2*dz)/k)*u_next[1]
      u_next[1]= (f[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u_gv))*dt+u[i]
      plot(t,u_next,type = "l", xlab = 'Time', ylab = 'Population of Phytoplankton')
    } else {
      Sys.sleep(0.1)
      u_next[i] = (f[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u[i]))*dt+u[i]
      plot(t,u_next,type = "l", xlab = 'Time', ylab = 'Population of Phytoplankton')
    }
  }
}

u_next