##########################################################
# Particle trajectory with different initial conditions
##########################################################

using Plots; gr()

k = 1/(4*pi*8.85e-12) 
e = 1.6e-19 
wholeTime = 1e-13
dt = 1e-20
b = 5e-14
Au_d = 10^(-14) #diameter of Au atom in meters
h = 0.00001 #accuracy parameter

#################### Depending on q #######################
v0 = 1e+7
m = 6.64456e-27
Q = 79*e

tx = [0, 1337]
ty = [0, 1337]
p_q = scatter(tx, ty,
        xlims = (-3e-13, 5e-13),
        ylims = (-3e-13, 7e-13),
        label = "Au",
        formatter = :scientific,
	title="Depending on q"
)

charges = [-2.85*e, -1.03*e, 0.23e, 1e-23, 3*e, 5*e]

for q in charges
	tmp = k * Q * q * dt / m
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
        plot!(p_q, x, y, label = q)
end

#################### Depending on Q #######################
v0 = 1e+7
m = 6.64456e-27
q = 2*e

tx = [0, 1337]
ty = [0, 1337]
p_Q = scatter(tx, ty,
        xlims = (-3e-13, 5e-13),
        ylims = (-3e-13, 7e-13),
        label = "Au",
        formatter = :scientific,
	title="Depending on Q"
)

charges = [-110*e, -50*e, -10*e, 5*e, 26*e, 40*e, 99*e]

for Q in charges
        tmp = k * Q * q * dt / m
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
        plot!(p_Q, x, y, label = Q)
end

#################### Depending on m #######################
v0 = 1e+7
Q = 79*e
q = 2*e

tx = [0, 1337]
ty = [0, 1337]
p_m = scatter(tx, ty,
        xlims = (-3e-13, 5e-13),
        ylims = (-3e-13, 7e-13),
        label = "Au",
        formatter = :scientific,
	title="Depending on m"
)

masses = [3e-29, 3e-27, 6e-27, 8e-27, 1e-26, 2.5e-26, 4e-25]

for m in masses
        tmp = k * Q * q * dt / m
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
        plot!(p_m, x, y, label = m)
end

#################### Depending on v0  #######################
m = 6e-27
Q = 79*e
q = 2*e

tx = [0, 1337]
ty = [0, 1337]
p_v0 = scatter(tx, ty,
        xlims = (-3e-13, 5e-13),
        ylims = (-3e-13, 7e-13),
        label = "Au",
        formatter = :scientific,
        title="Depending on v0"
)

speeds = [1e6, 5e6, 8.9e6, 9.43e6, 1.23e7, 1.43e7, 3e7]

for v0 in speeds
        tmp = k * Q * q * dt / m
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
        plot!(p_v0, x, y, label = v0)
end

plot(p_q, p_Q, p_m, p_v0, layout = 4)
plot!(size=(1000,1000))
savefig("plot2.png")
