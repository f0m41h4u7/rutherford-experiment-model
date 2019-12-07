###########################################################
# Build particle trajectory depending on a target parameter
###########################################################

using Plots; gr()

k = 1/(4*pi*8.85*10^(-12)) #constant: k=1/(4pi*e0)
e = 1.6*10^(-19) #charge of electron in coulombs
alpha_q = 2*e #charge of alpha-particle in coulombs
alpha_m = 6.64456*10^(-27) #mass of alpha-particle in kilos
Au_q = 79*e #charge of Au particle in coulombs
Au_d = 10^(-14) #diameter of Au atom in meters

v0 = 10^7 #initial velocity of alpha-particle

wholeTime = 10^(-13)
dt = 10^(-15)

h = 0.00001
tmp = k * Au_q * alpha_q * dt / alpha_m

params = [1e-14, 1e-13, 2e-13, 5e-13, -1e-13, -1e-14, 4e-13]

tx = [0, 1337]
ty = [0, 1337]
p = scatter(tx, ty, 
	lw=kostyl, 
	xlims = (-3e-13, 5e-13), 
	ylims = (-3e-13, 7e-13),
	label = "Au",
	formatter = :scientific
)

for b in params	
	i = 2
	x = Float64[]
	y = Float64[]
	v_x = Float64[]
	v_y = Float64[]

	push!(v_x, v0)
	push!(v_y, 0)

	push!(x, -(10*Au_d))
	push!(y, b)

	for t in dt:dt:wholeTime
		push!(v_x, v_x[i-1] + h * tmp * x[i-1] / ( ( x[i-1]^2 + y[i-1]^2 ) ^ (3/2) ) )
		push!(v_y, v_y[i-1] + h * tmp * y[i-1] / ( ( x[i-1]^2 + y[i-1]^2 ) ^ (3/2) ) )

		push!(x, x[i-1] + (h/2) * (v_x[i-1] + v_x[i]) * dt)
		push!(y, y[i-1] + (h/2) * (v_y[i-1] + v_y[i]) * dt)

		i += 1
	end
	plot!(p, x, y, label = b)
end

savefig(p, "plot1.png")
