function L(member::Dict{Integer,Array{Integer}}, node::Dict{Integer,Array{Float64}})
	l = Dict{Integer, Float64}()
	for e = 1:length(member)
		l[e] = hypot(node[member[e][2]][1]-node[member[e][1]][1], node[member[e][2]][2]-node[member[e][1]][2])
	end
	return l
end

function LocalK(member::Dict{Integer,Array{Integer}}, node::Dict{Integer,Array{Float64}}, E::Dict{Integer, Float64}, ν::Dict{Integer, Float64}, A::Dict{Integer, Float64}, I:: Dict{Integer, Array{Float64}})
	K = Dict{Integer, Array{Float64}}()
	l = L(member,node)
	for e = 1:length(member)
		K[e] = zeros(Float64,12,12)
		K[e][1,1] = E[e]*A[e]/l[e]
		K[e][1,7] = -E[e]*A[e]/l[e]
		K[e][2,2] = E[e]*12*I[e][3]/l[e]^3
		K[e][2,6] = E[e]*6*I[e][3]/l[e]^2
		K[e][2,8] = -E[e]*12*I[e][3]/l[e]^3
		K[e][2,12] = E[e]*6*I[e][3]/l[e]^2
		K[e][3,3] = E[e]*12*I[e][2]/l[e]^3
		K[e][3,5] = -E[e]*6*I[e][2]/l[e]^2
		K[e][3,9] = -E[e]*12*I[e][2]/l[e]^3
		K[e][3,11] = -E[e]*6*I[e][2]/l[e]^2
		K[e][4,4] = E[e]*I[e][1]/2/(1+ν[e])/l[e]
		K[e][4,10] = -E[e]*I[e][1]/2/(1+ν[e])/l[e]
		K[e][5,3] = -E[e]*6*I[e][2]/l[e]^2
		K[e][5,5] = E[e]*4*I[e][2]/l[e]
		K[e][5,9] = E[e]*6*I[e][2]/l[e]^2
		K[e][5,11] = E[e]*2*I[e][2]/l[e]
		K[e][6,2] = E[e]*6*I[e][3]/l[e]^2
		K[e][6,6] = E[e]*4*I[e][3]/l[e]
		K[e][6,8] = -E[e]*6*I[e][3]/l[e]^2
		K[e][6,12] = E[e]*2*I[e][3]/l[e]
		K[e][7,1] = -E[e]*A[e]/l[e]
		K[e][7,7] = E[e]*A[e]/l[e]
		K[e][8,2] = -E[e]*12*I[e][3]/l[e]^3
		K[e][8,6] = -E[e]*6*I[e][3]/l[e]^2
		K[e][8,8] = E[e]*12*I[e][3]/l[e]^3
		K[e][8,12] = -E[e]*6*I[e][3]/l[e]^2
		K[e][9,3] = -E[e]*12*I[e][2]/l[e]^3
		K[e][9,5] = E[e]*6*I[e][2]/l[e]^2
		K[e][9,9] = E[e]*12*I[e][2]/l[e]^3
		K[e][9,11] = E[e]*6*I[e][2]/l[e]^2
		K[e][10,4] = -E[e]*I[e][1]/2/(1+ν[e])/l[e]
		K[e][10,10] = E[e]*I[e][1]/2/(1+ν[e])/l[e]
		K[e][11,3] = -E[e]*6*I[e][2]/l[e]^2
		K[e][11,5] = E[e]*2*I[e][2]/l[e]
		K[e][11,9] = E[e]*6*I[e][2]/l[e]^2
		K[e][11,11] = E[e]*4*I[e][2]/l[e]
		K[e][12,2] = E[e]*6*I[e][3]/l[e]^2
		K[e][12,6] = E[e]*2*I[e][3]/l[e]
		K[e][12,8] = -E[e]*6*I[e][3]/l[e]^2
		K[e][12,12] = E[e]*4*I[e][3]/l[e]
	end
	return K
end