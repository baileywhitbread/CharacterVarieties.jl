## Testing E-polynomials
function palindrome_X(G::FiniteCoxeterGroup,genus::Int64,puncture_min::Int64,puncture_max::Int64)
	# Checks palindromicity of X with g=genus and n=puncture_min,...,puncture_max
	d=group_type_data(G)
	for n in puncture_min:puncture_max
		try 
			print("Checking E(X;q) palindromic when G=",G,", g=",genus," and n=",n,": ")
			if ispalindromic(bigint_EX(G,genus,n,d))
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
			if fast_bigint_EX(G,genus,n,d)(1)==0
				println("χ(X)=0")
			else
				println("χ(X)=",fast_bigint_EX(G,genus,n,d)(1))
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
			if isnonnegative(fast_bigint_EY(G,genus,n,d))
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
			if isnonnegative(fast_bigint_EX(G,genus,n,d))
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
	random_string = randstring(12)
	log_name = "EY_coefficients_G_"*string(G)*"_"*random_string*".txt"
	io = open(log_name, "w+")
	global_logger(SimpleLogger(io))
	d=algebra_type_data(G)
	for n in puncture_min:puncture_max
		try 
			@info("Checking coefficients of E(Y;q) when G=",G,", g=",genus," and n=",n,": ")
			if isnonnegative(fast_bigint_EY(G,genus,n,d))
				@info("All non-negative")
			else
				@info("Negative coefficients found")
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