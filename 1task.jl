###########################################################
# Particle trajectory depending on a target parameter
###########################################################

using Plots; gr()

k = 1/(4*pi*8.85e-12) 
e = 1.6e-19 
alpha_q = 2*e 
alpha_m = 6.64456e-27 
Au_q = 79*e 
Au_d = 1e-14

v0 = 1e7 

wholeTime = 1e-13
dt = 1e-15

h = 0.00001
tmp = k * Au_q * alpha_q * dt / alpha_m

params = [-1e-14, -1e-13, 1e-14, 1e-13, 2e-13, 4e-13, 5e-13]

tx = [0, 1337]
ty = [0, 1337]
p = scatter(tx, ty, 
	xlims = (-3e-13, 5e-13), 
	ylims = (-3e-13, 7e-13),
	label = "Au",
	formatter = :scientific,
	title = "Depending on b"
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
