############################################
# Phytoplankton mathematical modeling 2020 #
############################################

# Chosen Values
n0 = 100000              # Phytoplankton initial amount
vp = 0.1                 # Needs to be bigger then 0 but which number
k  = 0.1                 # Needs to be bigger then 0 but which number
dt = 0.5                 # Something we chose
dz = 0.5                 # Something we chose                                        
S  = 0.3                 # Something we chose

N = L/dz                 # Number of special intervals

# Dichlet Boundaries
L   = 100                # Since phytoplankton is not found depths deeper then 100m
z = seq(0, L , by = dz)  # Here L is used as max due to death in lower sea levels
t = seq(0, L , by = dt)  # Here L is used since 100 is a nice number for the time

# Initial example of the function                                               # Here work is needed
n = seq(n0,0,,2)
t = seq(0,L,,2)
plot(t,n,type = "l", xlab = 'Time', ylab = 'Population of Phytoplankton')

# Calculation of phytoplankton over time                                        # Here work is needed
S = cumsum(N)*dz

