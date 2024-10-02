## Display human-readable tables
function group_type_table(G::FiniteCoxeterGroup;summands=false,n=1::BigInt)
	if summands == false
		d = group_type_data(G)
		num_of_types = size(d)[1]
		# The next lines make the |L(Fq)| and ρ(1) columns readable
		# By converting the entries from Pol to CycPol
		for i in [2,3]
			d[i*num_of_types+1:(i+1)*num_of_types] = CycPol.(d[i*num_of_types+1:(i+1)*num_of_types])
		end
		clabels = ["|Φ(L)+|","|L(Fq)|","ρ(1)","χᵨ(1)","|W(L)|","|[L]|","ν(L)"];
		rlabels = xrepr.(Ref(rio()),d[:,1]); 
		# xrepr(rio(), __ ) is a string of __ when printed on the REPR
		repr_d = xrepr.(Ref(rio()),d[:,2:size(d)[2]]);
		println("A G-type is a W-orbit [L,ρ] where ")
		println("L is an endoscopy group of G containing T")
		println("ρ is a principal unipotent character of L(Fq)")
		println("Φ(L)+ is the set of positive roots of L")
		println("|L(Fq)| is the size of L(Fq)")
		println("ρ(1) is the degree of the unipotent character ρ")
		println("χᵨ(1) is the degree of the Weyl group character associated to ρ")
		println("W(L) is the Weyl group of L")
		println("[L] is the orbit of L under the W-action")
		println("ν(L) is an integer only depending on L")
		println("")
		return showtable(repr_d;col_labels=clabels,rows_label="Types [L,ρ]",row_labels=rlabels)
	elseif summands == true
		d = group_type_data(G)
		num_of_types = size(d)[1]
		# Add Mτ and Sτ to this data
		Mtau_column = Array{Any}(nothing,num_of_types,1)
		Stau_column = Array{Any}(nothing,num_of_types,1)
		for i in 1:num_of_types
			Mtau_column[i] = Pol{Rational{BigInt}}(fast_Mtau(G,i,d))
			Stau_column[i] = Pol{Rational{BigInt}}(fast_Stau(G,n,i,d))
		end
		d = hcat(d,Mtau_column)
		d = hcat(d,Stau_column)
		# The next lines make the |L(Fq)|, ρ(1), Mτ and Sτ columns readable
		# By converting the entries from Pol to CycPol
		for i in [2,3,8,9]
			d[i*num_of_types+1:(i+1)*num_of_types] = CycPol.(d[i*num_of_types+1:(i+1)*num_of_types])
		end
		clabels = ["|Φ(L)+|","|L(Fq)|","ρ(1)","χᵨ(1)","|W(L)|","|[L]|","ν(L)","Mτ","Sτ"];
		rlabels = xrepr.(Ref(rio()),d[:,1]); 
		# xrepr(rio(), __ ) is a string of __ when printed on the REPR
		repr_d = xrepr.(Ref(rio()),d[:,2:size(d)[2]]);
		println("A G-type is a W-orbit [L,ρ] where ")
		println("L is an endoscopy group of G containing T")
		println("ρ is a principal unipotent character of L(Fq)")
		println("Φ(L)+ is the set of positive roots of L")
		println("|L(Fq)| is the size of L(Fq)")
		println("ρ(1) is the degree of the unipotent character ρ")
		println("χᵨ(1) is the degree of the Weyl group character associated to ρ")
		println("W(L) is the Weyl group of L")
		println("[L] is the orbit of L under the W-action")
		println("ν(L) is an integer only depending on L")
		println("Mτ(q) is the q-mass of τ")
		println("Sτ(q) is the character sum of τ")
		println("")
		return showtable(repr_d;col_labels=clabels,rows_label="Types [L,ρ]",row_labels=rlabels)
	else 
		ArgumentError("The optional summands argument must be true or false; default is false")
	end
end





function algebra_type_table(G::FiniteCoxeterGroup;summands=false,g=1::BigInt,n=1::BigInt)
	if summands == false
		d = algebra_type_data(G)
		num_of_types = size(d)[1]
		# The next lines make the |N| and Q_T^L(N) columns readable
		# By converting the entries from Pol to CycPol
		for i in [2,3]
			d[i*num_of_types+1:(i+1)*num_of_types] = CycPol.(d[i*num_of_types+1:(i+1)*num_of_types])
		end
		clabels = ["d(τ)","|N|","Q_T^L(N)","|[L]|","µ(L,G)"];
		rlabels = xrepr.(Ref(rio()),d[:,1]); 
		# xrepr(rio(), __ ) is a string of __ when printed on the REPR
		repr_d = xrepr.(Ref(rio()),d[:,2:size(d)[2]]);
		println("A g-type is a W-orbit [L,N] where ")
		println("L is a Levi subgroup of G containing T")
		println("N is an L(Fq)-orbit of a nilpotent element of Lie(L)(Fq)")
		println("d(τ) is dim(L)-degree(|N|)")
		println("Q_T^L(N) is the Greens function evaluated at the unipotent L(Fq)-conjugacy class associated to N")
		println("[L] is the orbit of L under the W-action")
		println("µ(L,G) is the Mobius function of the poset of Levi subgroups of G containing T")
		println("")
		return showtable(repr_d;col_labels=clabels,rows_label="Types [L,N]",row_labels=rlabels)
	elseif summands == true
		d = algebra_type_data(G)
		num_of_types = size(d)[1]
		# Add qdτ and Hτ to this data
		qdtau_column = Array{Any}(nothing,num_of_types,1)
		Htau_column = Array{Any}(nothing,num_of_types,1)
		for i in 1:num_of_types
			qdtau_column[i] = Pol{BigInt}(fast_qdtau(G,i,d))
			Htau_column[i] = Pol{Rational{BigInt}}(fast_Htau(G,n,i,d))
		end
		d = hcat(d,qdtau_column)
		d = hcat(d,Htau_column)
		# The next lines make the |N| and Q_T^L(N) columns readable
		# By converting the entries from Pol to CycPol
		for i in [2,3]
			d[i*num_of_types+1:(i+1)*num_of_types] = CycPol.(d[i*num_of_types+1:(i+1)*num_of_types])
		end
		clabels = ["d(τ)","|N|","Q_T^L(N)","|[L]|","µ(L,G)","qdτ","Hτ"];
		rlabels = xrepr.(Ref(rio()),d[:,1]); 
		# xrepr(rio(), __ ) is a string of __ when printed on the REPR
		repr_d = xrepr.(Ref(rio()),d[:,2:size(d)[2]]);
		println("A g-type is a W-orbit [L,N] where ")
		println("L is a Levi subgroup of G containing T")
		println("N is an L(Fq)-orbit of a nilpotent element of Lie(L)(Fq)")
		println("d(τ) is dim(L)-degree(|N|)")
		println("Q_T^L(N) is the Greens function evaluated at the unipotent L(Fq)-conjugacy class associated to N")
		println("[L] is the orbit of L under the W-action")
		println("µ(L,G) is the Mobius function of the poset of Levi subgroups of G containing T")
		println("qdτ is q^d(τ)")
		println("Hτ(q) is the character sum of τ")
		println("")
		return showtable(repr_d;col_labels=clabels,rows_label="Types [L,N]",row_labels=rlabels)
	else 
		ArgumentError("The optional summands argument must be true or false; default is false")
	end
end