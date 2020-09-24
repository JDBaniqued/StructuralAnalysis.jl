using StructuralAnalysis
sa = StructuralAnalysis

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

### Forces
begin
	f = Dict{Integer, Array}()
	f[1] = [500*cos(40*π/180),500*sin(40*π/180)]
	f[2] = [0,missing]
	f[3] = [missing,missing]
	f
end

### Displacements
begin
	u = Dict{Integer, Array}()
	u[1] = [missing,missing]
	u[2] = [missing,0]
	u[3] = [0,0]
	u
end

k = sa.AxialLocalK(E,A,member,node)
KG = sa.AxialGlobalK(k,member,node)
Free = sa.FreeNodalDir(u)
KGAug = sa.AugmentGlobal(Free,KG)
Force = sa.GlobalForce(Free,f)
freedisp = inv(KGAug) * Force