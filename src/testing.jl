## Testing E-polynomials
function palindrome_X(G::FiniteCoxeterGroup,genus::Int64,puncture_min::Int64,puncture_max::Int64)
	# Checks palindromicity of X with g=genus and n=puncture_min,...,puncture_max
	d=group_type_data(G)
	for n in puncture_min:puncture_max
		try 
			print("Checking E(X;q) palindromic when G=",G,", g=",genus," and n=",n,": ")
			if ispalindromic(fast_EX(G,genus,n,d))
				println("Yes")
			else
				println("No")
			end
		catch err
			if isa(err,OverflowError)
				println("Overflow error")
				#break
			elseif isa(err,ErrorException)
				println("Error exception")
				#break
			else
				println(err)
				#break
			end
		end
	end
end



function euler_X(G::FiniteCoxeterGroup,genus::Int64,puncture_min::Int64,puncture_max::Int64)
	# Checks Euler characteristic of X with g=genus and n=puncture_min,...,puncture_max
	d=group_type_data(G)
	for n in puncture_min:puncture_max
		try 
			print("Computing χ(X) when G=",G,", g=",genus," and n=",n,": ")
			if fast_EX(G,genus,n,d)(1)==0
				println("χ(X)=0")
			else
				println("χ(X)=",fast_EX(G,genus,n,d)(1))
			end
		catch err
			if isa(err,OverflowError)
				println("Overflow error")
				#break
			elseif isa(err,ErrorException)
				println("Error exception")
				#break
			else
				println(err)
				#break
			end
		end
	end
end

function nonnegative_Y(G::FiniteCoxeterGroup,genus::Int64,puncture_min::Int64,puncture_max::Int64)
	# Checks negativity of coefficients of E(Y;q) with g=genus and n=puncture_min,...,puncture_max
	d=algebra_type_data(G)
	for n in puncture_min:puncture_max
		try 
			print("Checking coefficients of E(Y;q) when G=",G,", g=",genus," and n=",n,": ")
			if isnonnegative(fast_EY(G,genus,n,d))
				println("All non-negative")
			else
				println("Negative coefficients found")
			end
		catch err
			if isa(err,OverflowError)
				println("Overflow error")
				#break
			elseif isa(err,ErrorException)
				println("Error exception")
				#break
			else
				println(err)
				#break
			end
		end
	end
end




function nonnegative_X(G::FiniteCoxeterGroup,genus::Int64,puncture_min::Int64,puncture_max::Int64)
	# Checks negativity of coefficients of E(X;q) with g=genus and n=puncture_min,...,puncture_max
	d=group_type_data(G)
	for n in puncture_min:puncture_max
		try 
			print("Checking coefficients of E(X;q) when G=",G,", g=",genus," and n=",n,": ")
			if isnonnegative(fast_EX(G,genus,n,d))
				println("All non-negative")
			else
				println("Negative coefficients found")
			end
		catch err
			if isa(err,OverflowError)
				println("Overflow error")
				#break
			elseif isa(err,ErrorException)
				println("Error exception")
				#break
			else
				println(err)
				#break
			end
		end
	end
end



function log_nonnegative_Y(G::FiniteCoxeterGroup,genus::Int64,puncture_min::Int64,puncture_max::Int64)
	# Checks negativity of coefficients of E(Y;q) with g=genus and n=puncture_min,...,puncture_max
	log_name = "EY_nonnegative_coeff_"*xrepr(rio(),G)*"_g="*string(genus)*"_n="*string(puncture_min)*"..."*string(puncture_max)*"_"*randstring(12)*".txt"
	io = open(log_name, "w+")
	logger = SimpleLogger(io)
	global_logger(logger)
	d=algebra_type_data(G)
	for n in puncture_min:puncture_max
		try 
			if isnonnegative(fast_EY(G,genus,n,d))
				@info("Coefficients all non-negative when",n)
			else
				@info("Negative coefficients found when",n)
			end
		catch err
			if isa(err,OverflowError)
				println("Overflow error")
				#break
			elseif isa(err,ErrorException)
				println("Error exception")
				#break
			else
				println(err)
				#break
			end
		end
	end
	close(io)
end

function check_dim_X(G::FiniteCoxeterGroup,genus_min::Int64,genus_max::Int64,puncture_min::Int64,puncture_max::Int64)
	d=group_type_data(G)
	for n in puncture_min:puncture_max
		for g in genus_min:genus_max
			if dimension_XY(G,g,n) == fast_EX(G,g,n,d)
				println("Dimension correct when (g,n) is ","(",g,",",n,")")
			end
		end
	end
end

function check_dim_Y(G::FiniteCoxeterGroup,genus_min::Int64,genus_max::Int64,puncture_min::Int64,puncture_max::Int64)
	d=group_type_data(G)
	for n in puncture_min:puncture_max
		for g in genus_min:genus_max
			if dimension_XY(G,g,n) == fast_EY(G,g,n,d)
				println("Dimension correct when (g,n) is ","(",g,",",n,")")
			end
		end
	end
end