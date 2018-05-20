#=
Compute the effective medium stiffness tensor using the Voigt, Reuss, and
Hill average
=#

using Distributions, Homogenization

# We'll do an z-x-z euler rotation where both z-rotations are uniformly
# distributed and the x-distribution is normal. We'll do this for
nfab = 5
m = 500 # number of crystals in our sample
rotation_sequence = ["z" "x" "z"]


# the incidence will be along the z-axis for simplicity
theta = 0
rho = 910 # density in kg/m^3

# Define the stiffness tensor for a single crystal at a specific temperature
T = -2
P = 0.01 #1MPa



# we're only going to go to 40 percent porosity
npor = 1 #Int(round(0.4*m) )
pore_replace = sample(collect(1:41), npor-1 )

# define the standard deviations for each fabric
stds = [sqrt(pi/2)]# sqrt(pi/4) sqrt(pi/8) sqrt(pi/16) sqrt(pi/32)]
nfab = length(stds)
temps = collect(-40:0.5:-2)
ntemp = length(temps)


# Allocate space
vpv = zeros(ntemp, nfab)
vpr = zeros(ntemp, nfab)
vph = zeros(ntemp, nfab)

vsvv = zeros(ntemp, nfab)
vsvr = zeros(ntemp, nfab)
vsvh = zeros(ntemp, nfab)

vshv = zeros(ntemp, nfab)
vshr = zeros(ntemp, nfab)
vshh = zeros(ntemp, nfab)

# We'll save the VRH estimates for zero porosity and the maximum porosity
Cv = zeros(6, 6, 2)
Cr = zeros(6, 6, 2)
Ch = zeros(6, 6, 2)

# Initialize values that we need from the for loop
Cvoigt = []
Creuss = []
Chill = []

for j = 1:ntemp
    C = ice_stiffness(temps[j], P)
    S = inv(C)
    for k = 1:nfab
        euler_angles = [rand(m) rand(m) rand(m) ].*180./pi #

        # reset the voigt matrices
        Cvoigt =  zeros(6,6)
        Creuss = zeros(6,6)

        # Compute the homogenized stiffness tensor
        for i = 1:m
            R = euler_rot( euler_angles[i,:], rotation_sequence)
            M,N = bond_matrix(R)
            Cvoigt = Cvoigt + ( M*C*M' )
            Creuss = Creuss + N*S*N'
        end

        Cvoigt = Cvoigt./m
        Creuss = m.*inv(Creuss)
        Chill = (Cvoigt + Creuss)/(2)

        Cv[:, :, 1] = Cvoigt
        Cr[:, :, 1] = Creuss
        Ch[:, :, 1] = Chill

        vpv[j, k], vsvv[j,k], vshv[j,k] = seismic_velocity(Cvoigt, rho, theta)
        vpr[j, k], vsvr[j,k], vshr[j,k] = seismic_velocity(Creuss, rho, theta)
        vph[j, k], vsvh[j,k], vshh[j,k] = seismic_velocity(Chill, rho, theta)


    end

end


vpv = vpv.*10000
vpr = vpr.*10000
vph = vph.*10000
vsvv = vsvv.*10000
vsvr = vsvr.*10000
vsvh = vsvh.*10000
vshv = vshv.*10000
vshr = vshr.*10000
vshh = vshh.*10000


writecsv("voigt_pvel_temp", vpv)
writecsv("reuss_pvel_temp.csv", vpr)
writecsv("hill_pvel_temp.csv", vph)

writecsv("voigt_svvel_temp.csv", vsvv)
writecsv("reuss_svvel_temp.csv", vsvr)
writecsv("hill_svvel_temp.csv", vsvh)
