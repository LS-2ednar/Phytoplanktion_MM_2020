############################################
# Phytoplankton mathematical modeling 2020 #
############################################

# change this to your working directory
#setwd("C:/Users/wiese/OneDrive - ZHAW/0_HS20/Mathematical Modelling/MM_Phytoplankton")

# Chosen Values
u0  = 1000                # Initial amount of phytoplankton
vp  = 0.1                 # Velocity of phytoplankton downwards
k   = 4                   # Diffusion constant
dt  = 0.125               # Something we chose
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

z = seq(0, L,  by = dz)   # Here L is used as max due to death in lower sea levels
t = seq(0, LT, by = dt)   # Here LT is used as max due to death in lower sea levels
u = seq(u0,0,,N+1)        # initialize the Population u as 10000 phytoplankton to
                          # 0 over time
u_next = u                # Save for later use  

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
plot(z,u,type="l", ylab = 'Phytoplankton', xlab = 'z')

# 3. Finally the f the function for the which takes in account the amount of 
# phytoplankton here only the values of +1 or -1 is used.

for (j in 1:NT) {
#change over time 
  S =  cumsum(u*dz)                    # Determin new amount of accumulated Phytoplankton
  f =  ifelse(S > 5000, -0.005, 1.1)   # Determine if the amount of Phytoplankton grows 
                                       # or decays
  #print(j)
  for (i in 1:N){
    if (i == 1) {
      u_gv= u[2]+((vp*2*dz)/k)*u[1]
      u_next[1]= (f[i]*u[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u_gv))*dt+u[i]
     
    } else {
      u_next[i] = (f[i]*u[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u[i-1]))*dt+u[i]
      
    }
  }
  u=u_next
  plot(z,u_next,type = "l", xlab = 'z', ylab = 'Phytoplankton', main = 'Population of Phytoplankton', col = 'blue')
  Sys.sleep(0.1)
}

u_next
  
###                                     ###
#Initial conditions Conditions Plot saving# 
###                                     ###

jpeg(file="chosen_conditions.jpeg")
plot(z,u_next,type = "l", xlab = 'z', ylab = 'Phytoplankton', main = 'Population of Phytoplankton', col = 'blue')
dev.off()


################################################################################
#                               Stability tests                                #                   
################################################################################

#------------------------------------------------------------------------------#
#                              changing dt                                     #
#------------------------------------------------------------------------------#

################################
#Changeing dt from 0.125 to 0.5#
################################

# Chosen Values
dt  = 0.5                 # Value which was changed  
dz  = 0.5                                                         

# Dichlet Boundaries 
L  = 100                  # Since phytoplankton is not found in depths deeper then 100m
LT = 100                  # Final time point
t  = seq(0, L, by = dt)   # Here L is used since 100 is a nice number for the time
N  = L/dz                 # Number of special intervals
NT = LT/dt                # Number of time intervals

z = seq(0, L,  by = dz)   # Here L is used as max due to death in lower sea levels
t = seq(0, LT, by = dt)   # Here LT is used as max due to death in lower sea levels
u = seq(u0,0,,N+1)        # initialize the Population u as 10000 phytoplankton to
                          # 0 over time
u_next = u                # Save for later use  

#Recalculating Function Parts
A    =  (vp*dt) / dz                        
B    =  (k*dt) / (dz^2)

#Plot of initial conditions 
plot (t,u,type="l", ylab = 'Phytoplankton', xlab = 'z')

for (j in 1:NT) {
  #change over time 
  S =  cumsum(u*dz)                    # Determin new amount of accumulated Phytoplankton
  f =  ifelse(S > 5000, -0.005, 1.1)   # Determine if the amount of Phytoplankton grows 
  # or decays
  #print(j)
  for (i in 1:N){
    if (i == 1) {
      u_gv= u[2]+((vp*2*dz)/k)*u[1]
      u_next[1]= (f[i]*u[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u_gv))*dt+u[i]
      
    } else {
      u_next[i] = (f[i]*u[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u[i-1]))*dt+u[i]
      
    }
  }
  u=u_next
  plot(z,u_next,type = "l", xlab = 'z', ylab = 'Phytoplankton', main = 'Population of Phytoplankton (dt = 0.5)', col = 'red')
  Sys.sleep(0.1)
}

#Save this Plot
jpeg(file="dt_0_5.jpeg")
plot(z,u_next,type = "l", xlab = 'z', ylab = 'Phytoplankton', main = 'dt = 0.5', col = 'red')
dev.off()
u_next
####################################
#Changeing dt from 0.125 to 0.00125#
####################################

# Chosen Values
dt  = 0.00125             # Value which was changed  
dz  = 0.5                                                         

# Dichlet Boundaries 
L  = 100                  # Since phytoplankton is not found in depths deeper then 100m
LT = 100                  # Final time point
t  = seq(0, L, by = dt)   # Here L is used since 100 is a nice number for the time
N  = L/dz                 # Number of special intervals
NT = LT/dt                # Number of time intervals

z = seq(0, L,  by = dz)   # Here L is used as max due to death in lower sea levels
t = seq(0, LT, by = dt)   # Here LT is used as max due to death in lower sea levels
u = seq(u0,0,,N+1)        # initialize the Population u as 10000 phytoplankton to
                          # 0 over time
u_next = u                # Save for later use  

#Recalculating Function Parts
A    =  (vp*dt) / dz                        
B    =  (k*dt) / (dz^2)

#Plot of initial conditions 
plot (z,u,type="l", ylab = 'Phytoplankton', xlab = 'z')

for (j in 1:NT) {
  #change over time 
  S =  cumsum(u*dz)                    # Determin new amount of accumulated Phytoplankton
  f =  ifelse(S > 5000, -0.005, 1.1)   # Determine if the amount of Phytoplankton grows 
  # or decays
  #print(j)
  for (i in 1:N){
    if (i == 1) {
      u_gv= u[2]+((vp*2*dz)/k)*u[1]
      u_next[1]= (f[i]*u[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u_gv))*dt+u[i]
      
    } else {
      u_next[i] = (f[i]*u[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u[i-1]))*dt+u[i]
      
    }
  }
  u=u_next
  plot(z,u_next,type = "l", xlab = 'z', ylab = 'Phytoplankton', main = 'dt = 0.0125', col = 'darkgreen')
  Sys.sleep(0.1)
}

jpeg(file="dt_0_00125.jpeg")
plot(z,u_next,type = "l", xlab = 'z', ylab = 'Phytoplankton', main = 'dt = 0.0125', col = 'darkgreen')
dev.off()
u_next

#------------------------------------------------------------------------------#
#                              changeling dz                                   #
#------------------------------------------------------------------------------#

#################################
#Changeing dz from 0.5 to 0.125#
#################################

# Chosen Values
dt  = 0.125                 
dz  = 1.25                # Value which was changed                                    

# Dichlet Boundaries 
L  = 100                  # Since phytoplankton is not found in depths deeper then 100m
LT = 100                  # Final time point
t  = seq(0, L, by = dt)   # Here L is used since 100 is a nice number for the time
N  = L/dz                 # Number of special intervals
NT = LT/dt                # Number of time intervals

z = seq(0, L,  by = dz)   # Here L is used as max due to death in lower sea levels
t = seq(0, LT, by = dt)   # Here LT is used as max due to death in lower sea levels
u = seq(u0,0,,N+1)        # initialize the Population u as 10000 phytoplankton to
                          # 0 over time
u_next = u                # Save for later use  

#Recalculating Function Parts
A    =  (vp*dt) / dz                        
B    =  (k*dt) / (dz^2)

#Plot of initial conditions 
plot (z,u,type="l", ylab = 'Phytoplankton', xlab = 'z')

for (j in 1:NT) {
  #change over time 
  S =  cumsum(u*dz)                    # Determin new amount of accumulated Phytoplankton
  f =  ifelse(S > 5000, -0.005, 1.1)   # Determine if the amount of Phytoplankton grows 
  # or decays
  #print(j)
  for (i in 1:N){
    if (i == 1) {
      u_gv= u[2]+((vp*2*dz)/k)*u[1]
      u_next[1]= (f[i]*u[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u_gv))*dt+u[i]
      
    } else {
      u_next[i] = (f[i]*u[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u[i-1]))*dt+u[i]
      
    }
  }
  u=u_next
  plot(z,u_next,type = "l", xlab = 'z', ylab = 'Phytoplankton [-]', main = 'Population of Phyoplankton')
  Sys.sleep(0.1)
}

plot(z,u_next,type = "l", xlab = 'z', ylab = 'Phytoplankton', main = 'dz = 0.0125', col = 'darkgreen')


jpeg(file="dz_0_0125.jpeg")
plot(z,u_next,type = "l", xlab = 'z', ylab = 'Phytoplankton', main = 'dz = 0.0125', col = 'darkgreen')
dev.off()
u_next
####################################
#Changeing dz from 0.125 to 12.5#
####################################

# Chosen Values
dt  = 0.125              
dz  = 12.5                # Value which was changed                                            

# Dichlet Boundaries 
L  = 100                  # Since phytoplankton is not found in depths deeper then 100m
LT = 100                  # Final time point
t  = seq(0, L, by = dt)   # Here L is used since 100 is a nice number for the time
N  = L/dz                 # Number of special intervals
NT = LT/dt                # Number of time intervals

z = seq(0, L,  by = dz)   # Here L is used as max due to death in lower sea levels
t = seq(0, LT, by = dt)   # Here LT is used as max due to death in lower sea levels
u = seq(u0,0,,N+1)        # initialize the Population u as 10000 phytoplankton to
                          # 0 over time
u_next = u                # Save for later use  

#Recalculating Function Parts
A    =  (vp*dt) / dz                        
B    =  (k*dt) / (dz^2)

#Plot of initial conditions 
plot (z,u,type="l", ylab = 'Phytoplankton', xlab = 'z')

for (j in 1:NT) {
  #change over time 
  S =  cumsum(u*dz)                    # Determin new amount of accumulated Phytoplankton
  f =  ifelse(S > 5000, -0.005, 1.1)   # Determine if the amount of Phytoplankton grows 
  # or decays
  #print(j)
  for (i in 1:N){
    if (i == 1) {
      u_gv= u[2]+((vp*2*dz)/k)*u[1]
      u_next[1]= (f[i]*u[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u_gv))*dt+u[i]
      
    } else {
      u_next[i] = (f[i]*u[i]-A*(u[i+1]-u[i])+B*(u[i+1]-2*u[i]+u[i-1]))*dt+u[i]
      
    }
  }
  u=u_next
  plot(z,u_next,type = "l", xlab = 'z', ylab = 'Phytoplankton', main = 'Population of Phyoplankton (dz = 12.5)', col = 'red')
  Sys.sleep(0.1)
}

jpeg(file="dz_12_5.jpeg")
plot(z,u_next,type = "l", xlab = 'z', ylab = 'Phytoplankton', main = 'dz = 12.5', col = 'red')
dev.off()
u_next