# depaint.jl
# Integrates atomic trajectories in the presence of a rastering laser dipole force.  
module depaint

using Plots, DifferentialEquations

function getGaussianStylus(width,height)
	function stylus(x::Real,y::Real)
		return height*exp(-(x.^2 + y.^2)/width^2)
	end
	return stylus
end

function getPath1(R,fradial,ftheta)
	function spiralPath(t)
		return R*(mod(fradial*t,2)-1)*[cos(ftheta*t),sin(ftheta*t)]
	end
	return spiralPath
end

function evolve(path,Rstylus,amp,init,t0,tf)
	# Integrates motion of an atom under laser dipole force
		# path:		<function> path of laser raster; function of time
		# stylus:	<function> laser beam shape; function of (x,y)
		# Rcutoff:	<positive number> maximum extent of laser beam; used for optimization.
		# init:		<length 4 vector> atom initial condition, in format <x,y,px,py>
		# t0:		<number> initial time
		# tf:		<number> final time to integrate to
		# vel:		<positive number> approximate velocity of the laser path; used for optimization.
	
	# Initialize state variables
	x = init[1:2]	# Current position
	p = init[3:4]	# Current momentum
	t = t0			# Current time
	L = path(t0)	# Laser position
	
	function dudt(u,p,t)
		c = path(t)
		a = 2*amp*exp(-((u[1]-c[1])^2+(u[2]-c[2])^2)/(Rstylus^2))/Rstylus^2
		return [u[3],u[4],(u[1]-c[1])*a,(u[2]-c[2])*a]
	end
	println(dudt(init,0,0))
	prob = ODEProblem(dudt,init,(t0,tf))
	println("hi")
	sol = solve(prob)
	return sol
end

end