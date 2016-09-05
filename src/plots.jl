include("mycolors.jl")

function showmatrix(matrix::AbstractMatrix; minmax=[0,1], cmap="binary")
    fig = plt[:figure]()
    ax = fig[:add_axes]([0,0,1,1])
    ax[:imshow](matrix,interpolation="nearest",cmap=cmap, vmin=minmax[1], vmax=minmax[2])
    ax[:axis]("off")
    fig
end

function zhist{T<:Real}(zscores::Vector{T}; lim=3.0)
    x = collect(linspace(-3.5,3.5,40))
    fig = plt[:figure](figsize=(5,3))
    ax = fig[:add_subplot](111)
    ax[:hist](zscores, bins=x, color=BLUE, ec="white",normed=1)
    ax[:plot](x, pdf(Normal(),x),"--",color=RED, lw=2)
    ax[:set_ylim]([0,0.5])
    ax[:set_xlim]([-lim,lim])
    fig
end

function manhattanplot{T<:Real}(zscores::Vector{T}, highlight::AbstractVector{Bool}=Bool[])
    N = length(zscores)
    if length(highlight) != N
        highlight=BitArray(N)
        highlight[:] = false
    end

    neglogp = -log10(z2p(zscores))
    snpnumber = collect(1:N)
    range = [0,N+1]
    Bonf1 = -log10(0.05/N)
    Bonf2 = -log10(0.01/N)
    
    fig = plt[:figure](figsize=(5,3))
    ax = fig[:add_subplot](111)
    ax[:plot](snpnumber[~highlight], neglogp[~highlight], "xk", alpha=0.5)
    ax[:plot](snpnumber[highlight], neglogp[highlight], "o", color=RED)

    ax[:plot](range, [Bonf2, Bonf2], "--", color=GREEN)
    ax[:plot](range, [Bonf1, Bonf1], ":", color=GREEN)
    ax[:set_xlim](range)
    fig
end

function probhist(sample::Vector)
    @assert 0.0 <= minimum(sample) < maximum(sample) <= 1.0
    beta = fit(Beta, sample)
    mm = mean(sample)
    p = linspace(0,1,50)
    
    fig = plt[:figure](figsize=(5,3))
    ax = fig[:add_subplot](111)
    ax[:hist](sample, p, color=BLUE, alpha=1, normed=1)
    ax[:plot](p, pdf(beta,p), "--k", alpha=1, lw=1)
end

function uniform_quantile(N::Int)
    q = linspace(0,1,N+1)[1:end-1]
    q += 0.5*q[2]
    q 
end

function stratify_mask(scores::Vector{Float64}, N::Integer)
    @assert N >= 2
    lim = linspace(minimum(scores), maximum(scores), N+1)
    res = BitArray(length(scores), N)
    res[:,1] = scores .<= lim[2]
    for i in 2:N
        res[:,i] = lim[i] .< scores .<= lim[i+1]
    end
    res
end

function mask_cummulative_lim(scores::Vector{Float64}, lim::AbstractVector)
    mask = BitArray(length(scores), length(lim)+1)
    mask[1,:] = true
    for i in 2:size(mask)[2]
        mask[:,i] = scores .<= lim[i]
    end
    mask
end

function zz_scratterplot(zscores1::Vector{Float64}, zscores2::Vector{Float64};
                         cscores::Vector{Float64}=Float64[],
                         cscoreN::Int=10,
                         highlight::AbstractVector{Bool}=BitArray(0),
                         limthresh::Float64=12.0,
                         BF::Bool=false,
                         )

    N = length(zscores1)
    highlight = length(highlight)==0 ?  BitArray(N)*false : highlight

    @assert N == length(zscores2) == length(highlight)
    @assert N == length(cscores) || length(cscores) == 0

    fig = plt[:figure](figsize=(7,7))
    ax = fig[:add_subplot](111)
    ax[:set_aspect]("equal")
    maxlim = min(1.03*max( 1, maximum(zscores1), maximum(zscores2)),limthresh)
    minlim = max(1.03*min(-1, minimum(zscores1), minimum(zscores2)),-limthresh)
    ax[:set_xlim]([minlim,maxlim])
    ax[:set_ylim]([minlim,maxlim])
    ax[:plot]([0,0], [minlim,maxlim], "-k", alpha=0.5)
    ax[:plot]([minlim,maxlim], [0,0], "-k", alpha=0.5)

    if length(cscores) == 0
        ax[:plot](zscores1[~highlight], zscores2[~highlight],".k",color=GWASTools.BLUE,alpha=0.3)
        ax[:plot](zscores1[highlight], zscores2[highlight],".",color=GWASTools.RED)
    else
        cmap = ColorMap("viridis")
        display(cmap)
        mask = stratify_mask(cscores,cscoreN)
        mask[highlight,:] = false
        val = linspace(0, 1, cscoreN)
        for i in 1:cscoreN
            ax[:plot](zscores1[mask[:,i]],zscores2[mask[:,i]], ".", color=cmap(val[i]), alpha=0.3)
        end
        ax[:plot](zscores1[highlight], zscores2[highlight], ".", color=GWASTools.RED)
    end

    if BF
        theta = linspace(0,2π)
        r = p2z(5*10^-8.0)
        ax[:plot](r*sin(theta), r*cos(theta),"--k",lw=2)
    end
end
