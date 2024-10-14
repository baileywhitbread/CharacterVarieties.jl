function log_nonnegative_X(G::FiniteCoxeterGroup,g::Union{BigInt,Integer},puncture_min::Union{BigInt,Integer},puncture_max::Union{BigInt,Integer})
	# Checks negativity of coefficients of E(X;q) with g=genus and n=puncture_min,...,puncture_max
	log_name = "EX_nonnegative_coeff_"*xrepr(rio(),G)*"_g="*string(g)*"_n="*string(puncture_min)*"..."*string(puncture_max)*"_"*randstring(12)*".txt"
	io = open(log_name, "w+")
	logger = SimpleLogger(io)
	global_logger(logger)
	d = group_type_data(G)
	log_text = ""
	for n in puncture_min:puncture_max
		try 
			if isnonnegative(fast_EX(G,g,n,d))
				log_text *= "\nCoeffs of EX all non-neg when (G,g,n)=("*xrepr(rio(),G)*",$g,$n)"
			else
				log_text *= "\nNegative coeffs found when (G,g,n)=("*xrepr(rio(),G)*",$g,$n)"
			end
		catch err 
			println("$err")
		end
	end
	@info(log_text)
	close(io)
end

function log_nonnegative_Y(G::FiniteCoxeterGroup,g::Union{BigInt,Integer},puncture_min::Union{BigInt,Integer},puncture_max::Union{BigInt,Integer})
	# Checks negativity of coefficients of E(Y;q) with g=genus and n=puncture_min,...,puncture_max
	log_name = "EY_nonnegative_coeff_"*xrepr(rio(),G)*"_g="*string(g)*"_n="*string(puncture_min)*"..."*string(puncture_max)*"_"*randstring(12)*".txt"
	io = open(log_name, "w+")
	logger = SimpleLogger(io)
	global_logger(logger)
	d = algebra_type_data(G)
	log_text = ""
	for n in puncture_min:puncture_max
		try 
			if isnonnegative(fast_EY(G,g,n,d))
				log_text *= "\nCoeffs of EY all non-neg when (G,g,n)=("*xrepr(rio(),G)*",$g,$n)"
			else
				log_text *= "\nNegative coeffs found when (G,g,n)=("*xrepr(rio(),G)*",$g,$n)"
			end
		catch err 
			println("$err")
		end
	end
	@info(log_text)
	close(io)
end

function fast_log_nonnegative_Y(G::FiniteCoxeterGroup,g::Union{BigInt,Integer},puncture_min::Union{BigInt,Integer},puncture_max::Union{BigInt,Integer},type_data::Any)
	# Checks negativity of coefficients of E(Y;q) with g=genus and n=puncture_min,...,puncture_max
	log_name = "EY_nonnegative_coeff_"*xrepr(rio(),G)*"_g="*string(g)*"_n="*string(puncture_min)*"..."*string(puncture_max)*"_"*randstring(12)*".txt"
	io = open(log_name, "w+")
	logger = SimpleLogger(io)
	global_logger(logger)
	log_text = ""
	for n in puncture_min:puncture_max
		try 
			if isnonnegative(fast_EY(G,g,n,type_data))
				log_text *= "\nCoeffs of EY all non-neg when (G,g,n)=("*xrepr(rio(),G)*",$g,$n)"
			else
				log_text *= "\nNegative coeffs found when (G,g,n)=("*xrepr(rio(),G)*",$g,$n)"
			end
		catch err 
			println("$err")
		end
	end
	@info(log_text)
	close(io)
end