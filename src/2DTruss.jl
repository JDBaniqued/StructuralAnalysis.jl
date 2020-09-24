# Truss Finite Element Analysis

using LinearAlgebra

## Local Stiffness Matrix"

function L(member::Dict{Integer,Array{Integer}}, node::Dict{Integer,Array{Float64}})
	l = Dict{Integer, Float64}()
	for e = 1:length(member)
		l[e] = hypot(node[member[e][2]][1]-node[member[e][1]][1], node[member[e][2]][2]-node[member[e][1]][2])
	end
	return l
end

function cosϕ(member::Dict{Integer,Array{Integer}}, node::Dict{Integer,Array{Float64}})
	cosphi = Dict{Integer, Float64}()
	for e = 1:length(member)
		hor = node[member[e][2]][1]-node[member[e][1]][1]
		hyp = hypot(node[member[e][2]][1]-node[member[e][1]][1], node[member[e][2]][2]-node[member[e][1]][2])
		cosphi[e] = hor/hyp
	end
	return cosphi
end

function sinϕ(member::Dict{Integer,Array{Integer}}, node::Dict{Integer,Array{Float64}})
	sinphi = Dict{Integer, Float64}()
	for e = 1:length(member)
		hor = node[member[e][2]][2]-node[member[e][1]][2]
		hyp = hypot(node[member[e][2]][1]-node[member[e][1]][1], node[member[e][2]][2]-node[member[e][1]][2])
		sinphi[e] = hor/hyp
	end
	return sinphi
end

function AxialLocalK(E::Dict, A::Dict, Member::Dict{Integer,Array{Integer}}, Node::Dict{Integer,Array{Float64}})
	k = Dict{}()
	Length = L(Member,Node)
	Cosϕ = cosϕ(Member,Node)
	Sinϕ = sinϕ(Member,Node)
	for e = 1:length(Member)
		k[e] = E[e]*A[e]/Length[e]*
	[[Cosϕ[e]^2, Sinϕ[e]*Cosϕ[e], -Cosϕ[e]^2, -Sinϕ[e]*Cosϕ[e]],
	[Sinϕ[e]*Cosϕ[e], Sinϕ[e]^2, -Sinϕ[e]*Cosϕ[e], -Sinϕ[e]^2],
	[-Cosϕ[e]^2, -Sinϕ[e]*Cosϕ[e], Cosϕ[e]^2, Sinϕ[e]*Cosϕ[e]],
	[-Sinϕ[e]*Cosϕ[e], -Sinϕ[e]^2, Sinϕ[e]*Cosϕ[e], Sinϕ[e]^2]]
	end
	return k
end

function AxialGlobalK(k::Dict,member::Dict, node::Dict)
	K = zeros(length(node)*2, length(node)*2)
	for e = 1:length(member)
		a = member[e][1] * 2 - 1
		b = member[e][1] * 2
		c = member[e][2] * 2 - 1
		d = member[e][2] * 2
		K[a,a] += k[e][1][1]
		K[a,b] += k[e][1][2]
		K[a,c] += k[e][1][3]
		K[a,d] += k[e][1][4]
		K[b,a] += k[e][2][1]
		K[b,b] += k[e][2][2]
		K[b,c] += k[e][2][3]
		K[b,d] += k[e][2][4]
		K[c,a] += k[e][3][1]
		K[c,b] += k[e][3][2]
		K[c,c] += k[e][3][3]
		K[c,d] += k[e][3][4]
		K[d,a] += k[e][4][1]
		K[d,b] += k[e][4][2]
		K[d,c] += k[e][4][3]
		K[d,d] += k[e][4][4]
	end
	return K
end

### Determine free nodal directions
function FreeNodalDir(u::Dict)
	Free = Integer[]
	for i = 1:length(u)
		if ismissing(u[i][1])
			push!(Free, i*2-1)
		end
		if ismissing(u[i][2])
			push!(Free, i*2)
		end
	end
	return Free
end

### Augmented Global Stiffness Matrix
function AugmentGlobal(Free::Array,K::Array)
	K_Aug = zeros(length(Free),length(Free))
	for m = 1:length(Free)
		for n = 1:length(Free)
			global K_aug
			K_Aug[m,n] = K[Free[m],Free[n]]
		end
	end
	K_Aug
end

## Global Force Matrix"
function GlobalForce(Free::Array,f::Dict)
	Force = zeros(length(free))
	for m = 1:length(free)
		if free[m]%2 != 0
			Force[m]= f[0.5*free[m]+0.5][1]
		end
		if free[m]%2 == 0
			Force[m]= f[0.5*free[m]][2]
		end
	end
	return Force
end