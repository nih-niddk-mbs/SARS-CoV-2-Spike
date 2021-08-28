using DelimitedFiles
using Statistics
# using PyPlot

"""
Code to compute the Coefficient of Variation of selected Minima in coarse-grained
hydropathy using Moret and Zebende (MZ) Solvent Susceptible Surface Area exponents of a
target sequence compared to mutations of a reference sequence


"""

"""
scan_window(S)
Find the local minima and compute measures over a range of windows

S: Vector of the sequence, each element is a single letter string representing an amino acid

"""

function scan_window(S::Vector)
    sds = Array{Float64,1}(undef,0)
    range = 11:55
    for w in range
        psi = score(S,MZ,w,0)
        out = findminima(psi,40,[437,552,677,786,917,1155])
        push!(sds,out[3])
    end
    return range,sds,findmin(sds)
end

"""
findminima(S,w,centers)

Find the minima located within window w of specified center locations
SARS-CoV-2 MZ35 minima are located at [441,547,661,787,917];
"""
function findminima(psi,w=20,centers = [437,552,677,786,917,1155])
    v = Vector{Float64}(undef,0)
    locations = Vector{Int}(undef,0)
    for center in centers
        a = findmin(psi[center-w:center+w,2])
        push!(v,a[1])
        push!(locations,psi[center-w-1+a[2],1])
    end
    w = sort(v)[1:4]
    return std(v), mean(v), std(v)/mean(v),v,locations,std(w)/mean(w)
end

"""
findmaxima(S,w,centers )

Find the local maxima located within w of specified center locations
SARS-CoV-2 MZ35 maxima are located at [441,547,661,787,917];
"""
function findmaxima(psi,w=5,centers = [131,232,378,491,610,755,892,1115,1225])
    centers .-= 18
    v = Vector{Float64}(undef,0)
    locations = Vector{Int}(undef,0)
    for center in centers
        a = findmax(psi[center-w:center+w,2])
        push!(v,a[1])
        push!(locations,psi[center-w-1+a[2],1])
    end
    return std(v), mean(v), std(v)/mean(v),v,locations
end


"""
score(v,scale,w,window=0)

Compute the coarse-grained score given a scoring scale (e.g. MZ scale)
"""
function score(v,scale,w,window=0)
    h = score(v,scale)
    w2 = div(w,2)
    hcat(Vector(w2+1:length(v)-w2),[ avg(h,i-w2,i+w2,window) for i in w2+1:length(v)-w2 ])
end

function score(v,scale)
    [ scale[v] for v in v]
end

"""
avg(v,first,last,window)

Compute the average over a window between size and first
with window type Hanning or Rectangular
"""
function avg(v,first,last,window)
    delta = last - div(last+first,2)
    if window == 1
        hann = 2*cos.(pi*Vector(-delta:delta)/(2*delta+1)).^2
    # gauss = exp.(-(Vector(-delta:delta)/(2*delta)).^2)
        return mean(hann .* v[first:last])
    else
        return mean(v[first:last])
    end
end

"""
Random mutation control:

1) Reference sequence is mutated with random insertions, deletions, and subsitutions
2) Coarse-grained MZ score of mutated sequence is computed
3) SM is computed
4) Multiple samples are taken to obtain a distribution of the SM
5) Probabilty that target SM is obtained randomly is computed
"""

"""
significance(samples,S,window,w=35)

Compute significance of target SM given SM samples of mutated reference sequence
"""
function significance(samples,S,window,w=35)
    psi = score(S,MZ,w,window)
    standard = findedges(psi[:,2])[1]
    a = hist(samples,100)
    cdf = cumsum(a[1])/sum(a[1])
    ind =findfirst(a[2] .>= standard)
    cdf[ind]
end

"""
sample(S,total,inserts,deletes,subs,window=0,w=35,scale=MZ)

Collect sample SMs of mutated reference distribution S
"""
function sample(S,total,inserts,deletes,subs,window=0,w=35,scale=MZ)
    samples = Array{Float64,1}(undef,0)
    for i in 1:total
        S1 = mutation(S,inserts,deletes,subs)
        psi = score(S1,scale,w,window)
        push!(samples,findedges(psi[:,2])[1])
    end
    return samples
end


"""
mutation(S,inserts::Vector,deletes::Vector,subs::Int)

Mutate sequence S
inserts: Vector of tuples of the ranges of insertions  e.g. (372,374)
deletes: Vector of tuples of the ranges of deletions
subs: Number of substitutions
"""
function mutation(S,inserts::Vector,deletes::Vector,subs::Int)
    S1 = copy(S)
    for s in 1:subs
        substitution!(S1)
    end
    for n in inserts
        S1 = insertion(S1,n)
    end
    for d in deletes
        S1 = deletion(S1,d)
    end
    return S1
end

"""
mutation(S,subs,deletes)

deletes: Vector of tuples of the ranges of deletions
subs: Number of substitutions
"""

function mutation(S,subs,deletes)
    S1 = copy(S)
    for s in subs
        substitution!(S1,s)
    end
    for d in deletes
        S1 = deletion(S1,d)
    end
    return S1
end


"""
substitution!(S)

mutate a random amino acid
"""
function substitution!(S,mutation::Tuple)
    ind = mutation[1]
    S[ind] = mutation[2]
end
function substitution!(S)
    ind = rand(1:length(S))
    S[ind] = rand(aminoacids)
end

"""
insertion(S,n)

Insert n consecutive random amino acids
"""

function insertion(S::Vector,n::Int)
    N = length(S)
    S1 = Vector(undef,N+n)
    i = rand(Vector(1:N))
    S1[1:i] = S[1:i]
    S1[i+n+1:end] = S[i+1:end]
    for j = 1:n
        S1[i+j] = rand(aminoacids)
    end
    return S1
end

"""
deletion(S,n::Tuple)

Delete n specific amino acids
"""
function deletion(S,n::Tuple)
    N = length(S)
    S1 = Vector(undef,N-n[2]+n[1]-1)
    i = n[1] - 1
    S1[1:i] = S[1:i]
    S1[i+1:end] = S[n[2]+1:end]
    return S1
end
"""
deletion(S,n::Int)

Delete n random amino acids
"""

function deletion(S,n::Int)
    N = length(S)
    S1 = Vector(undef,N-n)
    i = rand(Vector(1:N-n))
    S1[1:i] = S[1:i]
    S1[i+1:end] = S[i+n+1:end]
    return S1
end

"""
Dictionary for MZ hydropathy exponents (shifted and rescaled)

"""
MZ =
Dict(
  "Q" => 105,
  "W" => 174,
  "T" => 135,
  "C" => 246,
  "P" => 121,
  "V" => 238,
  "L" => 197,
  "M" => 221,
  "N" => 113,
  "H" => 152,
  "A" => 157,
  "D" => 87,
  "G" => 156,
  "E" => 94,
  "Y" => 222,
  "I" => 222,
  "S" => 100,
  "K" => 69,
  "R" => 78,
  "F" => 218
)



"""
functions for reading sequences and
making dictionaries

"""
function makeMZdictionary(file::String)
    scale = load_MZscale(file)
    Dict( scale[i,1] => scale[i,2] for i in 1:size(scale,1))
end

function makeMZdictionary(scale::Vector)
    Dict( scale[i,1] => scale[i,2] for i in 1:size(scale,1))
end

function load_MZscale(file::String)
    readdlm(file,',')
end

function load_sequence(file)
    readdlm(file,',')[:,2]
end




"""
Codon transitions
"""

Sr(scale)=rand(scale[:,1],1255)

"""
subsingle!(S)

Mutations due to a single nucleotide substitution
"""
function subsingle!(S)
    ind = rand(1:length(S))
    S[ind] = rand(transitions[S[ind]])
end

transitions =
Dict(
"F"=> ["L","S"],
"L"=> ["F","S","P","I"],
"I"=> ["L","T","M"],
"M"=> ["I","T","V"],
"V"=> ["M","A"],
"S"=> ["F","L","P","Y","R","N"],
"P"=> ["S","L","T","H","Q"],
"T"=> ["P","I","M","A","K","N"],
"A"=> ["T","V","D","E"],
"Y"=> ["S","C"],
"H"=> ["P","Q","R"],
"Q"=> ["H","P","N","R"],
"N"=> ["Q","T","K","S"],
"K"=> ["N","T","D","R"],
"D"=> ["K","A","E","G"],
"E"=> ["D","A","G"],
"C"=> ["Y"],
"W"=> ["R"],
"R"=> ["W","H","Q","S","K","G"],
"G"=> ["R","D","E"]
)

aminoacids = [
"F",
"L",
"I",
"M",
"V",
"S",
"P",
"T",
"A",
"Y",
"H",
"Q",
"N",
"K",
"D",
"E",
"C",
"W",
"R",
"G"
]
