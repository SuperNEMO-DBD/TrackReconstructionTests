### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 61e83e95-d6ee-4f7f-95b0-41f746dce9e0
import Pkg

# ╔═╡ 09741a91-6f9f-451b-a64d-cd7e91f5e76e
Pkg.add(["PythonPlot", "PlutoUI", "UnROOT", "DataFrames", "StatsPlots", "FHist", "StatsBase", "Distributions", "Turing", "LinearAlgebra", "Optim", "LaTeXStrings", "Measurements"]);

# ╔═╡ 1f9e1c97-8f54-43a3-b5dd-cc4f6c315abd
begin
	using PlutoUI, UnROOT, DataFrames, StatsPlots, FHist, StatsBase, Distributions, Turing, LinearAlgebra, Optim, LaTeXStrings, Measurements
end

# ╔═╡ 476fec1b-f99f-429c-8ea7-75bd08d307fe
html"<button onclick='present()'>present</button>"

# ╔═╡ 3d3e035d-201c-4752-a05c-daa335e9beb8
md"""
# This notebook is part of the series of track reconstruction tests defined for comparison of two track reconstruction algorithms in SuperNEMO. The full test is described [here](https://github.com/SuperNEMO-DBD/TrackReconstructionTests/issues/2)
"""

# ╔═╡ 223fd264-2a10-43b6-a46c-da06589387b0
md"""
# In this notebook in particular, we compare the CAT vs Simulated vertices!
"""

# ╔═╡ 08a73ae4-7859-4e24-9f2b-2c12b05f3eee
md"""
# First, we initiate packages and load data.
"""

# ╔═╡ 972c1676-2016-44de-857c-f95a6133a373
pythonplot()

# ╔═╡ 4ae5cfbf-4ca7-451f-90b4-9d3df59be2c7
theme(:dao)

# ╔═╡ 5df18baa-0294-4486-9e17-bb49165e85f7
default(size = (600, 400), fontfamily="serif", legend=:outerright)

# ╔═╡ 14961967-bdd4-4661-9838-6033fd1d2c04
# load root file
f2nu = ROOTFile("../data/CAT/bb2nu/bb2nu_1M_events.root")

# ╔═╡ 48b77fff-6ce8-44d5-a792-2e0bfbab2c1f
md"""
# There are $(length(keys(f2nu["tree"]))) variables saved in the root file, we don't need all of them though. Just the ones concerning the vertex positions: reconstructed and simulated.
"""

# ╔═╡ 8d52909c-2b68-4d65-ba9a-e4edd8573c52
begin
	# extract reconstructed vertex positions into dataframe
	df_raw = LazyTree(f2nu, "tree", [
		"y1Reconstructed","z1Reconstructed","y2Reconstructed","z2Reconstructed",
		"y1Simulated","z1Simulated","y2Simulated","z2Simulated",
	]) |> DataFrame
	describe(df_raw)
end

# ╔═╡ 9d13a802-a591-4dd7-aecd-d6ac19f18ad3
md"""
### Notice that for each event there are 4 y-components and 4 z-components of the vertices. For example for the y-component, for a single event there are: y1Reco, y2Reco, y1Sim, y2Sim. The simulated vertices are both the same, however the reconstructed are not! We do not know which is which. We combine the events so that we have only 2 variables per axis: yReco, ySim. We do this by concatenating the columns y1, and y2, essentialy doubling the number of rows:
"""

# ╔═╡ 76707366-cb2d-4e8e-b1be-8cb1e91da18e
begin
	df_concat = DataFrame( 
		yReco = vcat(df_raw.y1Reconstructed, df_raw.y1Reconstructed),
		ySim = vcat(df_raw.y1Simulated, df_raw.y1Simulated),
		zReco = vcat(df_raw.z1Reconstructed, df_raw.z1Reconstructed),
		zSim = vcat(df_raw.z1Simulated, df_raw.z1Simulated),
	)
	describe(df_concat)
end

# ╔═╡ 48dde7ab-2782-4ef7-a43a-2d2f777d1532
md"""
# The overall efficiency of the data-cuts is: eff = $(nrow(df_raw) / 1e6 *100 |> round)%
"""

# ╔═╡ 28a60d88-78ec-443f-97e6-41116b3812f7
md"""
# Analysis:

### The goal here is to create 2 histograms:
1. 1D histogram of $\Delta y$
2. 1D histogram of $\Delta z$

### For each histogram we extract [68%, 95%, 99%] credible intervals (interpreted as simplified resolution):
"""

# ╔═╡ 9a1738be-25bb-483b-ace6-c19201c53db1
md"""
# First we create 2 new columns: 
* ``dy \equiv y_{rec} - y_{sim}``
* ``dz \equiv z_{rec} - z_{sim}``

"""

# ╔═╡ 7ef237dc-67a4-4daa-8419-ba58630c11d8
# transform df to get Δy, Δz, d data
begin
	df_concat.dy = df_concat.yReco .- df_concat.ySim
	df_concat.dz = df_concat.zReco .- df_concat.zSim
	describe(df_concat)
end

# ╔═╡ ee9d7340-667b-4d39-b397-87e2014c7f27
md"""
## Now that the important data is generated, let's take a look at what the dy, dz, dyz look like:

  1. We can see that for each variable there are outlier values (take a look at the max and compare to mean a median); i.e. `maximum(dy)` = **$( round(Int, maximum(df_concat.dy)) ) mm!!!**
  2. Comparing the mean and median values for dy and dz shows that a larger error is to be expected in dz compared to dy
  
"""

# ╔═╡ 8b05885c-51d8-4d26-b503-ab4a1d23594a
md"""
## Concerning the outlier values: 
- let's take a look at how many outlier there are by studying the quantiles of dy and dz
- [TO DO] we should take a look in flvisualize at these events and figure out what's happening! (Should figure out a way of getting the event ID of the outliers)
"""

# ╔═╡ 1bf770ac-4e73-4480-994f-116ee2a2ebc5
begin
	alpha = 0.995
	lq, uq = (1.0 - alpha)/2.0, 1.0-(1.0-alpha)/2.0
	@show dy_995q = quantile(df_concat.dy, [lq, uq])
	@show dz_995q = quantile(df_concat.dz, [lq, uq])
end

# ╔═╡ 90674156-6647-441f-a5fd-de961ebf7f61
md"""
It can be seen that 99.5% of events have vertices reconstructed that fall well within (a somewhat arbitrary) range of say ~ $\Delta \in (-150, 150) mm$ - it should be safe to ignore the crazy outliers and use this the binning of -150:150.
We filter these values.
"""

# ╔═╡ cc5d29fe-24ae-45b8-aeee-7e4e89d46575
# drop outliers:
df = filter(row -> -150.0 < row.dy < 150.0 &&  -150.0 < row.dz < 150.0, df_concat);

# ╔═╡ 158e7525-7348-4b00-a906-4ce91bee83ce
describe(df)

# ╔═╡ 1cf788c8-3d5b-4614-9824-2de60fb1fa51
md"""
# Now we can visualize the data and plot the credible intervals:
"""

# ╔═╡ cb8697f4-4f30-4fc4-8d6a-35e908798acc
slider = @bind bw Slider(vcat(1, collect(5:5:20)), default = 5)

# ╔═╡ b18679eb-6247-4c49-8f73-7b435a2bef2a
md"""
binWidth = $bw
"""

# ╔═╡ 885843e1-7f0a-4d6e-8321-5ac5dd16c58e
begin
	h1dy = StatsBase.fit(Histogram{Float32}, df.dy, -150.0:bw:150.0)
	h1dz = StatsBase.fit(Histogram{Float32}, df.dz, -150.0:bw:150.0)
	nothing
end

# ╔═╡ 1130b08a-df77-4e67-bf1d-f9b49154f023
let
	p1 = plot(h1dy, xlabel = "dy / mm", ylabel ="counts / $(bw) mm", label ="", c=1, legend=:outertop, lw = ifelse(bw==1, 0, 1), title="CAT vs Simulation dy")
	p2 = plot(h1dz, xlabel = "dz / mm", ylabel ="counts / $(bw) mm", label ="", c=2,legend=:outertop, lw = ifelse(bw==1, 0, 1), title="CAT vs Simulation dz")
	p3 = histogram2d(df.dy, df.dz, bins = (-240:bw:240,-240:bw:240), aspect =1, xlabel ="dy", ylabel="dz", title="2D Histogram", colorbar_scale=:log10)

	plot(p1, p2, p3, layout=grid(2,2), size = (700,600), thickness_scaling= 1.1, dpi =200)
	
	savefig("../plots/CAT_vs_SIM_2nubb.png")
	current()
end

# ╔═╡ 35d8c337-9017-4b72-9961-c21b792abdb6
md"""
## Next we fit the distributions to extract 68% credible interval - to be used as sigma uncertainty (?)
"""

# ╔═╡ c6274998-e9f4-4624-84e0-0bdefa67344c
md"""
- We use the following distributions to find which is best: Normal, Cauchy, Laplace distribution 

### We can create the Turing model as:
"""

# ╔═╡ db8739ea-de36-4e0a-8624-790edc39c606
md"""
To create models in Turing can be as easy as specifying the priors and their distribution (in this case we have 2 parameters: ``\mu`` and ``\sigma``) and the likelihood.
"""

# ╔═╡ ae57ded6-e27f-429a-b484-91d3605ef9bd
@model function model_Laplace(d)
	# priors
	μ ~ Normal(0, 0.1)
	σ ~ Uniform(1e-5, 1e3)

	# likelihood
	d .~ Laplace(μ, σ)
end

# ╔═╡ 0136ddcf-9d79-44de-be73-3ff66ee19728
md"""
Same for dz.
"""

# ╔═╡ c61f17a2-e5a7-49fc-8045-0a9376896aff
@model function model_Normal(d)
	# priors
	μ ~ Normal(0, 0.1)
	σ ~ Uniform(1e-5, 1e3)

	# likelihood
	d .~ Normal(μ, σ)
end

# ╔═╡ 0a984258-1fb9-4159-9bb6-0293175d4744
@model function model_Cauchy(d)
	# priors
	μ ~ Normal(0, 0.1)
	σ ~ Uniform(1e-5, 1e3)

	# likelihood
	d .~ Cauchy(μ, σ)
end

# ╔═╡ b1d94a2f-84ce-4de2-b8e2-e188c5c9019d
md"""
## Create chain object used in MCMC

> For each model we use the Metropolis-Hastings algorithm with 1000 samples per chain, 1 chain. This is very simple task, so such *loose* conditions should be fine.
"""

# ╔═╡ 32bff16a-ad4c-4d3f-a064-60b5d9086c32
begin
	chain_dy_L = Turing.sample(model_Laplace(df.dy), MH(), 1000);
	chain_dy_C = Turing.sample(model_Cauchy(df.dy), MH(), 1000);
	chain_dy_N = Turing.sample(model_Normal(df.dy), MH(), 1000);
	nothing # to show no output
end

# ╔═╡ fb25d005-a68f-426a-b0bd-2269e8c5ae31
begin
	chain_dz_L = Turing.sample(model_Laplace(df.dz), MH(), 1000);
	chain_dz_C = Turing.sample(model_Cauchy(df.dz), MH(), 1000);
	chain_dz_N = Turing.sample(model_Normal(df.dz), MH(), 1000);
	nothing # to show no output
end

# ╔═╡ e021b46f-edd3-444b-bef8-c1593e9a6963
md"""
# We can investigate the fitted distributions: 
"""

# ╔═╡ 14b117a1-3a3f-4df7-b47a-99f055e2c282
describe(chain_dy_L)

# ╔═╡ 7fc26e8b-4db0-4d05-b194-00a21186f1e3
describe(chain_dy_C)

# ╔═╡ e916d332-68d5-4fde-bfc8-6e65f78827e4
describe(chain_dy_N)

# ╔═╡ e2256956-1eef-4f1d-b46b-0a46e49ba19c
md"""
# Now we extract the mean fitted parameters and create the fitting function
"""

# ╔═╡ ba5a743d-2be9-41cb-866b-39abeb5feb38
begin
	# now we extract the fitted parameters:
	(mu_dy_L, sigma_dy_L) = mean(chain_dy_L[:,1,:]), mean(chain_dy_L[:,2,:])
	(mu_dy_C, sigma_dy_C) = mean(chain_dy_C[:,1,:]), mean(chain_dy_C[:,2,:])
	(mu_dy_N, sigma_dy_N) = mean(chain_dy_N[:,1,:]), mean(chain_dy_N[:,2,:])
	(mu_dz_L, sigma_dz_L) = mean(chain_dz_L[:,1,:]), mean(chain_dz_L[:,2,:])
	(mu_dz_C, sigma_dz_C) = mean(chain_dz_C[:,1,:]), mean(chain_dz_C[:,2,:])
	(mu_dz_N, sigma_dz_N) = mean(chain_dz_N[:,1,:]), mean(chain_dz_N[:,2,:])
end

# ╔═╡ 38ea7de3-c3dd-484b-9e38-8fd277c9a1c1
begin
	# fit functions for plotting:
	fit_function_L( mu, sigma, x ) = pdf( Laplace(mu, sigma), x )
	fit_function_C( mu, sigma, x ) = pdf( Cauchy(mu, sigma), x )
	fit_function_N( mu, sigma, x ) = pdf( Normal(mu, sigma), x )
end

# ╔═╡ eba9e296-ad1f-4884-acc3-65acc477247c
md"""
# Finally we can visialize the results:
First dy
"""

# ╔═╡ 0221a814-5e62-4b57-a2e6-a7b3ae266aa4
savefig( "../plots/CAT_vs_SIM_dy_fit.png" )

# ╔═╡ c469af66-eee7-4802-9174-a767adacbcb3
md"""
#
Now dz
"""

# ╔═╡ ddc4094b-dd23-4db5-9512-d932baf587a6
savefig( "../plots/CAT_vs_SIM_dz_fit.png" )

# ╔═╡ e78dccc3-0dc3-477c-9542-a40c8dbcb2db
md"""
# Finally we find that the uncertainty on the CAT vertex reconstruction precision is:
* ``68\% CI \approx (-30.0, 29.9)`` mm with rms_y = 36.2 mm 
* ``68\% CI \approx (-28.7, 28.7)`` mm with rms_z = 36.7 mm 
"""

# ╔═╡ c94dd369-5440-4e3c-a2c5-6c4ac8311ae8
md"""
# Further studies:

- It would be nice to investigate the uncertainty as a function of energy, escape direction w.r.t. foil (see Miro's thesis)
- Maybe it's not the greatest idea in the world to fit the distribution, probably not much useful information is obtained this way, should it be enough to compare data quantiles? 

"""

# ╔═╡ 4be27a94-2d31-4229-9d5d-317f275b7f9f
"""
    r-squared coefficient of fit as defined in (https://en.wikipedia.org/wiki/Coefficient_of_determination)
"""
function r2(fit_function, histogram)
    yᵢ = histogram.weights
    xᵢ = collect( midpoints(histogram.edges[1]) )
    fᵢ = fit_function.(xᵢ)
    ybar = mean(yᵢ)

    SStot = sum( (yᵢ .- ybar).^2 )
    SSres = sum( (yᵢ .- fᵢ  ).^2 )

    return 1-( SSres / SStot)
end

# ╔═╡ 594e62b0-c302-45d7-852e-d4c6c7674cad
with(gr()) do
	theme(:dao)
	fdyL(x) = step(h1dy.edges[1])*nrow(df)*fit_function_L( mu_dy_L, sigma_dy_L, x)
	fdyC(x) = step(h1dy.edges[1])*nrow(df)*fit_function_C( mu_dy_C, sigma_dy_C, x)
	fdyN(x) = step(h1dy.edges[1])*nrow(df)*fit_function_N( mu_dy_N, sigma_dy_N, x)
	
	pmSigmadyL = quantile(Laplace(mu_dy_L, sigma_dy_L), [0.15865, 0.84135]) # [0.15865, 0.84135] is equal to 68% central quantiles
	CILabeldyL = "68% CI " *L"\in~" * "($(round(pmSigmadyL[1], sigdigits=3)), $(round(pmSigmadyL[2], sigdigits=3))) mm" 

	pmSigmadyC = quantile(Cauchy(mu_dy_C, sigma_dy_C), [0.15865, 0.84135]) # [0.15865, 0.84135] is equal to 68% central quantiles
	CILabeldyC = "68% CI " *L"\in~" * "($(round(pmSigmadyC[1], sigdigits=3)), $(round(pmSigmadyC[2], sigdigits=3))) mm" 
	
	plot(1, title = "dy distribution", dpi =200, label="", legend =:outerbottom, legend_column=2, size = (700, 700))
	plot!( # plot Laplace fit
		-200:0.1:200, x-> fdyL(x), 
		label ="Laplace " *L"r^2"* " = $(round(r2(fdyL, h1dy), sigdigits=4))\n" * CILabeldyL, labeltitle="fit options", c=3, lw=2
	)
	plot!( # plot Cauchy fit
		-200:0.1:200, x-> fdyC(x), 
		label ="Cauchy " *L"r^2"* " = $(round(r2(fdyC, h1dy), sigdigits=4))\n" * CILabeldyC,  
		c=4, lw=2
	)
	plot!( # plot Normal fit
		-200:0.1:200, x-> fdyN(x), 
		label ="Normal " *L"r^2"* " = $(round(r2(fdyN, h1dy), sigdigits=4))",  
		c=5, lw=2
	)
	plot!( # plot data
		h1dy, 
		label ="data", 
		c=2, lw = 2, st=:step, widen =:false
	)
	ylims!(0, 1.3*maximum(h1dy.weights))
	xlabel!("dy [mm]")
	ylabel!("counts / $(bw) mm")
end

# ╔═╡ 0b549061-2199-4339-86c7-1bf1dc4a508b
with(gr()) do
	theme(:dao)
	fdzL(x) = step(h1dz.edges[1])*nrow(df)*fit_function_L( mu_dz_L, sigma_dz_L, x)
	fdzC(x) = step(h1dz.edges[1])*nrow(df)*fit_function_C( mu_dz_C, sigma_dz_C, x)
	fdzN(x) = step(h1dz.edges[1])*nrow(df)*fit_function_N( mu_dz_N, sigma_dz_N, x)
	
	pmSigmadzL = quantile(Laplace(mu_dz_L, sigma_dz_L), [0.15865, 0.84135]) # [0.15865, 0.84135] is equal to 68% central quantiles
	CILabeldzL = "68% CI " *L"\in~" * "($(round(pmSigmadzL[1], sigdigits=3)), $(round(pmSigmadzL[2], sigdigits=3))) mm" 

	pmSigmadzC = quantile(Cauchy(mu_dz_C, sigma_dz_C), [0.15865, 0.84135]) # [0.15865, 0.84135] is equal to 68% central quantiles
	CILabeldzC = "68% CI " *L"\in~" * "($(round(pmSigmadzC[1], sigdigits=3)), $(round(pmSigmadzC[2], sigdigits=3))) mm" 
	
	plot(1, title = "dz distribution", dpi =200, label="", legend =:outerbottom, legend_column=2, size = (700, 700))
	plot!( # plot Laplace fit
		-200:0.1:200, x-> fdzL(x), 
		label ="Laplace " *L"r^2"* " = $(round(r2(fdzL, h1dz), sigdigits=4))\n" * CILabeldzL, labeltitle="fit options", c=3, lw=2
	)
	plot!( # plot Cauchy fit
		-200:0.1:200, x-> fdzC(x), 
		label ="Cauchy " *L"r^2"* " = $(round(r2(fdzC, h1dz), sigdigits=4))\n" * CILabeldzC,  
		c=4, lw=2
	)
	plot!( # plot Normal fit
		-200:0.1:200, x-> fdzN(x), 
		label ="Normal " *L"r^2"* " = $(round(r2(fdzN, h1dz), sigdigits=4))",  
		c=5, lw=2
	)
	plot!( # plot data
		h1dz, 
		label ="data", 
		c=2, lw = 2, st=:step, widen =:false
	)
	ylims!(0, 1.3*maximum(h1dz.weights))
	xlabel!("dy [mm]")
	ylabel!("counts / $(bw) mm")
end

# ╔═╡ ac2dfaef-1de4-428e-93d3-6d5c7a89d526
let fh1y = Hist1D( df.dy; binedges= -150:bw:150 )
	theme(:dao)
	plot(bincenters(fh1y), bincounts(fh1y), st=:stepmid, legend=:topright, label="data")
	q1 = quantile(df.dy, [0.15865, 0.84135])
	s1 = restrict(fh1y, floor(Int,q1[1]), ceil(Int,q1[2]))
	plot!(bincenters(s1), bincounts(s1), st=:stepmid, c= 1, f=0, fa=0.4, lw =0, label="68% of data in ($(q1[1]|>round), $(q1[2]|>round)) mm")
	xlabel!("dy [mm]")
	ylabel!("counts / $bw mm")
	ylims!(0, 5e4)
end

# ╔═╡ b317e8ea-a13f-4206-872a-ed90da5dbdc8
let fh1z = Hist1D( df.dz; binedges= -150:bw:150 )
	theme(:dao)
	plot(bincenters(fh1z), bincounts(fh1z), st=:stepmid, legend=:topright, label="data", c=2)
	q1 = quantile(df.dz, [0.15865, 0.84135])
	s1 = restrict(fh1z, floor(Int,q1[1]), ceil(Int,q1[2]))
	plot!(bincenters(s1), bincounts(s1), st=:stepmid, c= 2, f=0, fa=0.4, lw =0, label="68% of data in ($(q1[1]|>round), $(q1[2]|>round)) mm")
	xlabel!("dz [mm]")
	ylabel!("counts / $bw mm")
	ylims!(0, 5e4)
end

# ╔═╡ Cell order:
# ╟─476fec1b-f99f-429c-8ea7-75bd08d307fe
# ╟─3d3e035d-201c-4752-a05c-daa335e9beb8
# ╟─223fd264-2a10-43b6-a46c-da06589387b0
# ╠═08a73ae4-7859-4e24-9f2b-2c12b05f3eee
# ╠═61e83e95-d6ee-4f7f-95b0-41f746dce9e0
# ╠═09741a91-6f9f-451b-a64d-cd7e91f5e76e
# ╠═1f9e1c97-8f54-43a3-b5dd-cc4f6c315abd
# ╠═972c1676-2016-44de-857c-f95a6133a373
# ╠═4ae5cfbf-4ca7-451f-90b4-9d3df59be2c7
# ╠═5df18baa-0294-4486-9e17-bb49165e85f7
# ╠═14961967-bdd4-4661-9838-6033fd1d2c04
# ╟─48b77fff-6ce8-44d5-a792-2e0bfbab2c1f
# ╠═8d52909c-2b68-4d65-ba9a-e4edd8573c52
# ╟─9d13a802-a591-4dd7-aecd-d6ac19f18ad3
# ╠═76707366-cb2d-4e8e-b1be-8cb1e91da18e
# ╟─48dde7ab-2782-4ef7-a43a-2d2f777d1532
# ╟─28a60d88-78ec-443f-97e6-41116b3812f7
# ╟─9a1738be-25bb-483b-ace6-c19201c53db1
# ╠═7ef237dc-67a4-4daa-8419-ba58630c11d8
# ╟─ee9d7340-667b-4d39-b397-87e2014c7f27
# ╠═8b05885c-51d8-4d26-b503-ab4a1d23594a
# ╠═1bf770ac-4e73-4480-994f-116ee2a2ebc5
# ╟─90674156-6647-441f-a5fd-de961ebf7f61
# ╠═cc5d29fe-24ae-45b8-aeee-7e4e89d46575
# ╠═158e7525-7348-4b00-a906-4ce91bee83ce
# ╠═1cf788c8-3d5b-4614-9824-2de60fb1fa51
# ╠═cb8697f4-4f30-4fc4-8d6a-35e908798acc
# ╟─b18679eb-6247-4c49-8f73-7b435a2bef2a
# ╠═885843e1-7f0a-4d6e-8321-5ac5dd16c58e
# ╠═1130b08a-df77-4e67-bf1d-f9b49154f023
# ╠═35d8c337-9017-4b72-9961-c21b792abdb6
# ╠═c6274998-e9f4-4624-84e0-0bdefa67344c
# ╠═db8739ea-de36-4e0a-8624-790edc39c606
# ╠═ae57ded6-e27f-429a-b484-91d3605ef9bd
# ╠═0136ddcf-9d79-44de-be73-3ff66ee19728
# ╠═c61f17a2-e5a7-49fc-8045-0a9376896aff
# ╠═0a984258-1fb9-4159-9bb6-0293175d4744
# ╠═b1d94a2f-84ce-4de2-b8e2-e188c5c9019d
# ╠═32bff16a-ad4c-4d3f-a064-60b5d9086c32
# ╠═fb25d005-a68f-426a-b0bd-2269e8c5ae31
# ╠═e021b46f-edd3-444b-bef8-c1593e9a6963
# ╠═14b117a1-3a3f-4df7-b47a-99f055e2c282
# ╠═7fc26e8b-4db0-4d05-b194-00a21186f1e3
# ╠═e916d332-68d5-4fde-bfc8-6e65f78827e4
# ╠═e2256956-1eef-4f1d-b46b-0a46e49ba19c
# ╠═ba5a743d-2be9-41cb-866b-39abeb5feb38
# ╠═38ea7de3-c3dd-484b-9e38-8fd277c9a1c1
# ╠═eba9e296-ad1f-4884-acc3-65acc477247c
# ╠═594e62b0-c302-45d7-852e-d4c6c7674cad
# ╠═0221a814-5e62-4b57-a2e6-a7b3ae266aa4
# ╠═c469af66-eee7-4802-9174-a767adacbcb3
# ╠═0b549061-2199-4339-86c7-1bf1dc4a508b
# ╠═ddc4094b-dd23-4db5-9512-d932baf587a6
# ╠═e78dccc3-0dc3-477c-9542-a40c8dbcb2db
# ╠═c94dd369-5440-4e3c-a2c5-6c4ac8311ae8
# ╠═4be27a94-2d31-4229-9d5d-317f275b7f9f
# ╠═ac2dfaef-1de4-428e-93d3-6d5c7a89d526
# ╠═b317e8ea-a13f-4206-872a-ed90da5dbdc8
