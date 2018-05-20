#=
Compute the effective medium stiffness tensor using the Voigt, Reuss, and
Hill average
=#

using Distributions, Homogenization

# We'll do an z-x-z euler rotation where both z-rotations are uniformly
# distributed and the x-distribution is normal. We'll do this for
nfab = 5
m = 5000 # number of crystals in our sample
rotation_sequence = ["z" "x" "z"]


# the incidence will be along the z-axis for simplicity
theta = 0
rho = 910 # density in kg/m^3

# Define the stiffness tensor for a single crystal at a specific temperature
T = -10
P = 0.01 #1MPa
C = ice_stiffness(T, P)
S = inv(C)


# we're only going to go to 40 percent porosity
npor = 1 #Int(round(0.4*m) )
pore_replace = sample(collect(1:41), npor-1 )

# define the standard deviations for each fabric
stds = [sqrt(pi/2)]# sqrt(pi/4) sqrt(pi/8) sqrt(pi/16) sqrt(pi/32)]
nfab = length(stds)



# Allocate space
vpv = zeros(npor, nfab)
vpr = zeros(npor, nfab)
vph = zeros(npor, nfab)

vsvv = zeros(npor, nfab)
vsvr = zeros(npor, nfab)
vsvh = zeros(npor, nfab)

vshv = zeros(npor, nfab)
vshr = zeros(npor, nfab)
vshh = zeros(npor, nfab)

# We'll save the VRH estimates for zero porosity and the maximum porosity
Cv = zeros(6, 6, 2)
Cr = zeros(6, 6, 2)
Ch = zeros(6, 6, 2)

# Initialize values that we need from the for loop
Cvoigt = []
Creuss = []
Chill = []

for k = 1:nfab
    euler_angles = [rand(m) rand(Normal(0, stds[k] ), m) rand(m) ].*180./pi #

    # reset the voigt matrices
    Cvoigt =  zeros(6,6)
    Creuss = zeros(6,6)

    # Let's assume that each crystal is uniform in size so the fractional volume is
    # the same so we can scale it at the end
    #F = ones(m)./m


    # Compute the homogenized stiffness tensor
    for i = 1:m
        R = euler_rot( euler_angles[i,:], rotation_sequence)
        M,N = bond_matrix(R)
        Cvoigt = Cvoigt + M*C*M' # + F[i].*( M*C*M' )
        Creuss = Creuss + N*S*N' #F[i].*( N*S*N' )
    end

    Creuss = inv(Creuss)
    Chill = (Cvoigt + Creuss)/(2)

    Cv[:, :, 1] = Cvoigt
    Cr[:, :, 1] = Creuss
    Ch[:, :, 1] = Chill

    vpv[1, k], vsvv[1,k], vshv[1,k] = seismic_velocity(Cvoigt, rho, theta)
    vpr[1, k], vsvr[1,k], vshr[1,k] = seismic_velocity(Creuss, rho, theta)
    vph[1, k], vsvh[1,k], vshh[1,k] = seismic_velocity(Chill, rho, theta)

    # get a value for each porosity
    for j = 2:npor

        # reset the voigt matrices
        Cvoigt =  zeros(6,6)
        Creuss = zeros(6,6)
        #F[pore_replace[j-1]] = 0

        for i = 1:m
            R = euler_rot( euler_angles[i,:], rotation_sequence)
            M,N = bond_matrix(R)
            Cvoigt = Cvoigt + M*C*M' #Cvoigt + F[i].*( M*C*M' )
            Creuss = Creuss + N*S*N' #F[i].*( N*S*N' )
        end

        Cvoigt = Cvoigt./m
        Creuss = inv(Creuss)./m
        Chill = (Cvoigt + Creuss)/(2)

        vpv[j, k], vsvv[j,k], vshv[j,k] = seismic_velocity(Cvoigt, rho, theta)
        vpr[j, k], vsvr[j,k], vshr[j,k] = seismic_velocity(Creuss, rho, theta)
        vph[j, k], vsvh[j,k], vshh[j,k] = seismic_velocity(Chill, rho, theta)

    end

end


Cv[:, :, 2] = Cvoigt
Cr[:, :, 2] = Creuss
Ch[:, :, 2] = Chill

# Plot the results
using RCall
x = collect(0:npor-1)

# change velocities to km/s
vpv = vpv.*10000
vpr = vpr.*10000
vph = vph.*10000
vsvv = vsvv.*10000
vsvr = vsvr.*10000
vsvh = vsvh.*10000
vshv = vshv.*10000
vshr = vshr.*10000
vshh = vshh.*10000


writecsv("voigt_pvel_neg2C_01MPa.csv", vpv)
writecsv("reuss_pvel_neg2C_01MPa.csv", vpr)
writecsv("hill_pvel_neg2C_01MPa.csv", vph)

writecsv("voigt_svvel_neg2C_01MPa.csv", vsvv)
writecsv("reuss_svvel_neg2C_01MPa.csv", vsvr)
writecsv("hill_svvel_neg2C_01MPa.csv", vsvh)

ylimits = [ minimum(vpr[:,1]) maximum(vpv[:,1]) ]
#=
R"""
#png("VRH_Estimates_neg2C.png")
plot($x, $vph[,1], col = "black", type = "l", ylim = $ylimits, xlab = "Porosity (%)", ylab = "P-wave Velocity (m/s)", main = paste("Temperature =", $T)  )
lines($x, $vpv[,1], col = "black")
lines($x, $vpr[,1], col = "blue")

lines($x, $vph[,$nfab], col = "red", lty = 2)
lines($x, $vpv[,$nfab], col = "red", lty = 2)
lines($x, $vpr[,$nfab], col = "red", lty = 2)
#dev.off()
"""
=#
