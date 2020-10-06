### Nodes
begin
	node = Dict{Integer,Array{Float64}}()
	node[1] = [4000/tan(30*π/180),4000]
	node[2] = [10928,0]
	node[3] = [0,0]
	node
end

### Members
begin
	member = Dict{Integer,Array{Integer}}()
	member[1] = [1,2]
	member[2] = [1,3]
	member[3] = [2,3]
	member
end

node[member[1][2]][1]

### Modulus of Elasticity
begin
	E = Dict{Integer, Float64}()
	E[1] = 200
	E[2] = 200
	E[3] = 200
	E
end

### Cross-Sectional Area
begin
	A = Dict{Integer, Float64}()
	A[1] = 20000
	A[2] = 15000
	A[3] = 18000
	A
end



function L(member::Dict{Integer,Array{Integer}}, node::Dict{Integer,Array{Float64}})
	l = Dict{Integer, Float64}()
	for e = 1:length(member)
		l[e] = hypot(node[member[e][2]][1]-node[member[e][1]][1], node[member[e][2]][2]-node[member[e][1]][2])
	end
	return l
end

function LocalK(E::Dict, A::Dict, I::Dict, ν::Dict, Member::Dict{Integer,Array{Integer}}, Node::Dict{Integer,Array{Float64}})
	K = Dict{}()
	Length = L(Member,Node)
	for e = 1:length(Member)
		K[e] = E[e]*[[A[e]/Length[e],0,0,0,0,0,-A[e]/Length[e],0,0,0,0,0],
		[0,12*I[e][2]/Length[e]^3,0,0,0,6*I[e][2]/Length[e]^2,0,-12*I[e][2]/Length[e]^3,0,0,0,6*I[e][2]/Length[e]^2],
		[0,0,12*I[e][1]/Length[e]^3,0,-6*I[e][1]/Length[e]^2,0,0,0,-12*I[e][1]/Length[e]^3,0,-6*I[e][1]/Length[e]^2,0],
		[0,0,0,(I[e][1]+I[e][2])/(2*(1+ν[e])*Length[e]^3),0,0,0,0,0,-(I[e][1]+I[e][2])/(2*(1+ν[e])*Length[e]^3),0,0],
		[0,0,-6*I[e][1]/Length[e]^2,0,4*I[e][1]/Length[e],0,0,0,6*I[e][1]/Length[e]^2,0,2*I[e][1]/Length[e],0],
		[0,6*I[e][2]/Length[e]^2,0,0,0,4*I[e][2]/Length[e],0,0,-6*I[e][2]/Length[e]^2,0,0,0,2*I[e][2]/Length[e],0],
		[-A[e]/Length[e],0,0,0,0,0,A[e]/Length[e],0,0,0,0,0],
		[0,-12*I[e][2]/Length[e]^3,0,0,0,-6*I[e][2]/Length[e]^2,0,12*I[e][2]/Length[e]^3,0,0,0,-6*I[e][2]/Length[e]^2],
		[0,0,-12*I[e][1]/Length[e]^3,0,6*I[e][1]/Length[e]^2,0,0,0,12*I[e][1]/Length[e]^3,0,6*I[e][1]/Length[e]^2,0],
		[0,0,0,-(I[e][1]+I[e][2])/(2*(1+ν[e])*Length[e]^3),0,0,0,0,0,(I[e][1]+I[e][2])/(2*(1+ν[e])*Length[e]^3),0,0],
		[0,0,-6*I[e][1]/Length[e]^2,0,2*I[e][1]/Length[e],0,0,0,6*I[e][1]/Length[e]^2,0,4*I[e][1]/Length[e],0],
		[0,6*I[e][2]/Length[e]^2,0,0,0,2*I[e][2]/Length[e],0,0,-6*I[e][2]/Length[e]^2,0,0,0,4*I[e][2]/Length[e],0]]
	end
	return K
end