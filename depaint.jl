# depaint.jl
# Integrates atomic trajectories in the presence of a rastering laser dipole force.  

using Plots, DifferentialEquations

function getGaussianStylus(width,height)
	function stylus(x::Real,y::Real)
		return height*exp(-(x.^2 + y.^2)/width^2)
	end
	return stylus
end

function getPath

function evolve(path,stylus,Rcutoff,init,t0,tf,vel)
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
	states = []		# will store state data over time, as (t,x,y,px,py) vectors
	
"""	while t < tf
		if norm(x-L) > Rcutoff
			
		else
			
		end
		
		
	end"""
end