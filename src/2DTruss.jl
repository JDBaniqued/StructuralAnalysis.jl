### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ e15dd650-fb20-11ea-0430-db6f140f8052
md"# Truss Finite Element Analysis"

# ╔═╡ 25074f30-fb21-11ea-2287-4d9af37e41fa
md"## Model Description"

# ╔═╡ 33d91a20-fb21-11ea-0090-d32238c11685
md"### Nodes"

# ╔═╡ 52b9eaf0-fb21-11ea-1621-470a7d53e6f6
begin
	node = Dict{Integer,Array{Float64}}()
	node[1] = [4000/tan(30*π/180),4000]
	node[2] = [10928,0]
	node[3] = [0,0]
	node
end

# ╔═╡ ccf5e840-fbbd-11ea-033e-059d3630806d
md"### Members"

# ╔═╡ d2451810-fbb4-11ea-2620-9982c7135dc7
begin
	member = Dict{Integer,Array{Integer}}()
	member[1] = [1,2]
	member[2] = [1,3]
	member[3] = [2,3]
	member
end

# ╔═╡ db18fca0-fbbd-11ea-1dd3-0f6de312c90c
md"### Modulus of Elasticity"

# ╔═╡ ee8adc90-fbbd-11ea-1717-839320c5edbb
begin
	E = Dict{Integer, Float64}()
	E[1] = 200
	E[2] = 200
	E[3] = 200
	E
end

# ╔═╡ 183efa80-fbbe-11ea-120f-c91e6a8f883d
md"### Cross Sectional Area"

# ╔═╡ 248adde0-fbbe-11ea-2190-87ab49de85be
begin
	A = Dict{Integer, Float64}()
	A[1] = 20000
	A[2] = 15000
	A[3] = 18000
	A
end

# ╔═╡ 16f2b450-fbd2-11ea-30bd-7751c90f730d
md"### Forces"

# ╔═╡ 297419c0-fbd2-11ea-3fc7-1f200584578c
begin
	f = Dict{Integer, Array}()
	f[1] = [500*cos(40*π/180),500*sin(40*π/180)]
	f[2] = [0,missing]
	f[3] = [missing,missing]
	f
end

# ╔═╡ b6d185f0-fbd2-11ea-2d92-192bed8ee971
md"### Displacements"

# ╔═╡ c6346d50-fbd2-11ea-0029-cfd7141e59f1
begin
	u = Dict{Integer, Array}()
	u[1] = [missing,missing]
	u[2] = [missing,0]
	u[3] = [0,0]
	u
end

# ╔═╡ f2b33ede-fb20-11ea-21d4-5792f5aeb11b
md"## Local Stiffness Matrix"

# ╔═╡ 1685c720-fb21-11ea-1ce1-310f5e226283
function L(a::Array{Integer})
	return hypot(node[a[1]][1]-node[a[2]][1],node[a[1]][2]-node[a[2]][2])
end

# ╔═╡ df170212-fbb9-11ea-248a-8152e27dfa4a
function cosϕ(a::Array{Integer})
	hor = node[a[2]][1]-node[a[1]][1]
	hyp = hypot(node[a[1]][1]-node[a[2]][1],node[a[1]][2]-node[a[2]][2])
	return hor/hyp
end

# ╔═╡ e424c5a0-fbbc-11ea-2a89-8f910bc843d3
function sinϕ(a::Array{Integer})
	ver = node[a[2]][2]-node[a[1]][2]
	hyp = hypot(node[a[1]][1]-node[a[2]][1],node[a[1]][2]-node[a[2]][2])
	ver/hyp
end

# ╔═╡ a55ba700-fc28-11ea-0f52-93d14e401382
function AxialLocalK(E::Dict, A::Dict, Member::Dict)
	k = Dict{}()
	for e = 1:length(member)
		k[e] = E[e]*A[e]/L(Member[e])*
	[[cosϕ(Member[e])^2, sinϕ(Member[e])*cosϕ(Member[e]), -cosϕ(Member[e])^2, -sinϕ(Member[e])*cosϕ(Member[e])],
	[sinϕ(Member[e])*cosϕ(Member[e]), sinϕ(Member[e])^2, -sinϕ(Member[e])*cosϕ(Member[e]), -sinϕ(Member[e])^2],
	[-cosϕ(Member[e])^2, -sinϕ(Member[e])*cosϕ(Member[e]), cosϕ(Member[e])^2, sinϕ(Member[e])*cosϕ(Member[e])],
	[-sinϕ(Member[e])*cosϕ(Member[e]), -sinϕ(Member[e])^2, sinϕ(Member[e])*cosϕ(Member[e]),sinϕ(Member[e])^2]]
	end
	return k
end

# ╔═╡ ee27e020-fc2d-11ea-0a1b-6704ccef1be5
k = AxialLocalK(E,A,member)

# ╔═╡ e5c5c640-fbce-11ea-1bda-b553d206003c
md"## Global Stiffness Matrix"

# ╔═╡ 95c1bf30-fc2a-11ea-371d-7d891001a9ff
function AxialGlobalK(k::Dict,member::Dict)
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

# ╔═╡ 3f9b1f40-fc2d-11ea-3ea5-5de05b43505f
K = AxialGlobalK(k,member)

# ╔═╡ 66d04890-fc2f-11ea-3f7d-af3e8fb0c568
typeof(u)

# ╔═╡ c10995c0-fbd8-11ea-1079-d379ab6eb740
md"### Determine free nodal directions"

# ╔═╡ a2200a62-fc30-11ea-063e-a5ef1c6260d0
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

# ╔═╡ c73cecf0-fc30-11ea-2d5a-cf996ba3d2fe
free = FreeNodalDir(u)

# ╔═╡ 8b517ffe-fc32-11ea-1d4e-071aec9b2fcb
typeof(free)

# ╔═╡ e6862480-fbd8-11ea-2ccb-b75584838f33
md"### Augmented Global Stiffness Matrix"

# ╔═╡ 512b7bb0-fc32-11ea-1ebd-61241a27a867
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

# ╔═╡ ab200c30-fc32-11ea-3b2b-71c328c69150
KAug = AugmentGlobal(free,K)

# ╔═╡ 447e5c90-fbd1-11ea-2352-0717d77ffd3a
md"## Global Force Matrix"

# ╔═╡ 4c7e51e0-fc33-11ea-1628-cb4a385cdb1f
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

# ╔═╡ 8c587340-fc33-11ea-0a0a-edb42cb25639
F = GlobalForce(free,f)

# ╔═╡ a4d21630-fbdc-11ea-398d-5d59079dbe16
begin
	using LinearAlgebra
	U = inv(KAug) * F
end

# ╔═╡ 3d4ffa90-fbdc-11ea-17b5-f91d8ad196f2
md"## Determining the displacements"

# ╔═╡ Cell order:
# ╟─e15dd650-fb20-11ea-0430-db6f140f8052
# ╟─25074f30-fb21-11ea-2287-4d9af37e41fa
# ╟─33d91a20-fb21-11ea-0090-d32238c11685
# ╠═52b9eaf0-fb21-11ea-1621-470a7d53e6f6
# ╟─ccf5e840-fbbd-11ea-033e-059d3630806d
# ╠═d2451810-fbb4-11ea-2620-9982c7135dc7
# ╟─db18fca0-fbbd-11ea-1dd3-0f6de312c90c
# ╠═ee8adc90-fbbd-11ea-1717-839320c5edbb
# ╟─183efa80-fbbe-11ea-120f-c91e6a8f883d
# ╠═248adde0-fbbe-11ea-2190-87ab49de85be
# ╟─16f2b450-fbd2-11ea-30bd-7751c90f730d
# ╠═297419c0-fbd2-11ea-3fc7-1f200584578c
# ╟─b6d185f0-fbd2-11ea-2d92-192bed8ee971
# ╠═c6346d50-fbd2-11ea-0029-cfd7141e59f1
# ╟─f2b33ede-fb20-11ea-21d4-5792f5aeb11b
# ╟─1685c720-fb21-11ea-1ce1-310f5e226283
# ╟─df170212-fbb9-11ea-248a-8152e27dfa4a
# ╟─e424c5a0-fbbc-11ea-2a89-8f910bc843d3
# ╟─a55ba700-fc28-11ea-0f52-93d14e401382
# ╠═ee27e020-fc2d-11ea-0a1b-6704ccef1be5
# ╟─e5c5c640-fbce-11ea-1bda-b553d206003c
# ╟─95c1bf30-fc2a-11ea-371d-7d891001a9ff
# ╠═3f9b1f40-fc2d-11ea-3ea5-5de05b43505f
# ╠═66d04890-fc2f-11ea-3f7d-af3e8fb0c568
# ╟─c10995c0-fbd8-11ea-1079-d379ab6eb740
# ╟─a2200a62-fc30-11ea-063e-a5ef1c6260d0
# ╠═c73cecf0-fc30-11ea-2d5a-cf996ba3d2fe
# ╠═8b517ffe-fc32-11ea-1d4e-071aec9b2fcb
# ╟─e6862480-fbd8-11ea-2ccb-b75584838f33
# ╟─512b7bb0-fc32-11ea-1ebd-61241a27a867
# ╠═ab200c30-fc32-11ea-3b2b-71c328c69150
# ╟─447e5c90-fbd1-11ea-2352-0717d77ffd3a
# ╟─4c7e51e0-fc33-11ea-1628-cb4a385cdb1f
# ╠═8c587340-fc33-11ea-0a0a-edb42cb25639
# ╟─3d4ffa90-fbdc-11ea-17b5-f91d8ad196f2
# ╠═a4d21630-fbdc-11ea-398d-5d59079dbe16
