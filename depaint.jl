# depaint.jl
# Integrates atomic trajectories in the presence of a rastering laser dipole force.  
module depaint

using Plots, DifferentialEquations, ApproxFun

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

function evolve(path,Rstylus,amp,init,tf)
	# Integrates motion of an atom under laser dipole force
		# path:		<function> path of laser raster; function of time
		# stylus:	<function> laser beam shape; function of (x,y)
		# Rcutoff:	<positive number> maximum extent of laser beam; used for optimization.
		# init:		<length 4 vector> atom initial condition, in format <x,y,px,py>
		# t0:		<number> initial time
		# tf:		<number> final time to integrate to
		# vel:		<positive number> approximate velocity of the laser path; used for optimization.
	
	t0=0		# Initial time
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

function see(times,path,Rstylus,sol)
	# Plots laser and particle trajectories
	p = [path(t)[i] for t in times, i=1:2]
	s = [sol(t)[i] for t in times, i=1:2]
	xmax = findmax(cat(p[:,1],s[:,1];dims=1))[1]
	xmin = findmin(cat(p[:,1],s[:,1];dims=1))[1]
	ymax = findmax(cat(p[:,2],s[:,2];dims=1))[1]
	ymin = findmin(cat(p[:,2],s[:,2];dims=1))[1]
	
	mradius = 500/2*Rstylus/(xmax-xmin)
	@gif for i in 1:length(times)
		println("Progress: ",i," of ",length(times))
		scatter([p[i,1]],[p[i,2]],legend=false,xlim=(xmin-Rstylus,xmax+Rstylus),ylim=(ymin-Rstylus,ymax+Rstylus),size=(500,500),marker=(:circle,mradius,.4,:red))
		scatter!([s[i,1]],[s[i,2]],marker=(:circle,10,.7,:blue))
	end
end

function unpaint(path,potential,tf,maxvel=1.0)
	# Determines the speed required to yield a path which traces out a desired potential.
	#dudt(u,p,t) = 
end

function ensemble(path,Rstylus,amp,temp,tf)
	# Determines evolution of ensemble of atoms.
end

end