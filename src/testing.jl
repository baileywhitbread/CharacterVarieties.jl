## Testing E-polynomials
function palindrome_X(G::FiniteCoxeterGroup,g,puncture_min,puncture_max)
	# Checks palindromicity of X with g=genus and n=puncture_min,...,puncture_max
	d=group_type_data(G)
	for n in puncture_min:puncture_max
		try 
			print("Checking E(X;q) palindromic when (G,g,n)=("*xrepr(rio(),G)*",$g,$n): ")
			if ispalindromic(fast_EX(G,g,n,d))
				println("Yes")
			else
				println("No")
			end
		catch err 
			println("$err")
		end
	end
end



function euler_X(G::FiniteCoxeterGroup,g,puncture_min,puncture_max)
	# Checks Euler characteristic of X with g=genus and n=puncture_min,...,puncture_max
	d=group_type_data(G)
	for n in puncture_min:puncture_max
		try 
			print("If (G,g,n)=("*xrepr(rio(),G)*",$g,$n) then ")
			if fast_EX(G,g,n,d)(1)==0
				println("χ(X)=0")
			else
				println("χ(X)=",fast_EX(G,g,n,d)(1))
			end
		catch err 
			println("$err")
		end
	end
end

function nonnegative_Y(G::FiniteCoxeterGroup,g,puncture_min,puncture_max)
	# Checks negativity of coefficients of E(Y;q) with g=genus and n=puncture_min,...,puncture_max
	d=algebra_type_data(G)
	for n in puncture_min:puncture_max
		try 
			print("Checking coefficients of E(Y;q) when (G,g,n)=("*xrepr(rio(),G)*",$g,$n): ")
			if isnonnegative(fast_EY(G,g,n,d))
				println("All non-negative")
			else
				println("Negative coefficients found")
			end
		catch err 
			println("$err")
		end
	end
end




function nonnegative_X(G::FiniteCoxeterGroup,g,puncture_min,puncture_max)
	# Checks negativity of coefficients of E(X;q) with g=genus and n=puncture_min,...,puncture_max
	d=group_type_data(G)
	for n in puncture_min:puncture_max
		try 
			print("Checking coefficients of E(X;q) when (G,g,n)=("*xrepr(rio(),G)*",$g,$n): ")
			if isnonnegative(fast_EX(G,g,n,d))
				println("All non-negative")
			else
				println("Negative coefficients found")
			end
		catch err 
			println("$err")
		end
	end
end



function log_nonnegative_Y(G::FiniteCoxeterGroup,g,puncture_min,puncture_max)
	# Checks negativity of coefficients of E(Y;q) with g=genus and n=puncture_min,...,puncture_max
	log_name = "EY_nonnegative_coeff_"*xrepr(rio(),G)*"_g="*string(g)*"_n="*string(puncture_min)*"..."*string(puncture_max)*"_"*randstring(12)*".txt"
	io = open(log_name, "w+")
	logger = SimpleLogger(io)
	global_logger(logger)
	d=algebra_type_data(G)
	for n in puncture_min:puncture_max
		try 
			if isnonnegative(fast_EY(G,g,n,d))
				@info("Coefficients all non-negative when (G,g,n)=("*xrepr(rio(),G)*",$g,$n): ")
			else
				@info("Negative coefficients found when (G,g,n)=("*xrepr(rio(),G)*",$g,$n): ")
			end
		catch err 
			println("$err")
		end
	end
	close(io)
end

function check_dim_X(G::FiniteCoxeterGroup,genus_min,genus_max,puncture_min,puncture_max)
	for n in puncture_min:puncture_max
		for g in genus_min:genus_max
			if dimension_XY(G,g,n) == degree(fast_EX(G,g,n,d))
				println("Dimension correct when (G,g,n)=("*xrepr(rio(),G)*",$g,$n): ")
			else
				println("Dimension wrong when (G,g,n)=("*xrepr(rio(),G)*",$g,$n): ")
			end
		end
	end
end

function check_dim_Y(G::FiniteCoxeterGroup,genus_min,genus_max,puncture_min,puncture_max)
	d=algebra_type_data(G)
	for n in puncture_min:puncture_max
		for g in genus_min:genus_max
			if dimension_XY(G,g,n) == degree(fast_EY(G,g,n,d))
				println("Dimension correct when (G,g,n)=("*xrepr(rio(),G)*",$g,$n): ")
			else
				println("Dimension wrong when (G,g,n)=("*xrepr(rio(),G)*",$g,$n): ")
			end
		end
	end
end