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

# ╔═╡ 6b9feb78-371e-446d-83ca-cc0da8c65526
import Pkg

# ╔═╡ c9e89e94-1185-42c9-bab2-02628fb7b481
Pkg.add(["PythonPlot", "PlutoUI", "UnROOT", "DataFrames", "StatsPlots", "FHist", "StatsBase", "Distributions", "Turing", "LinearAlgebra", "Optim", "LaTeXStrings", "Measurements"]);

# ╔═╡ d3370197-b6c6-436c-8a54-30118a741278
begin
	using PlutoUI, UnROOT, DataFrames, StatsPlots, FHist, StatsBase, Distributions, Turing, LinearAlgebra, Optim, LaTeXStrings, Measurements
end

# ╔═╡ 8d855d70-1557-49ae-b995-bd8d4cf7a455
html"<button onclick='present()'>present</button>"

# ╔═╡ c6572c84-de77-40e1-9ca4-d9e72f73ac76
html"""
<style>
  main {
    max-width: 950px;
  }
</style>
"""

# ╔═╡ e4e23a9c-cfbd-467b-a3a5-f2477fc66531
md"""
# Initiate packages and load data.
"""

# ╔═╡ 9ef3374f-bd99-4a41-893a-e231329397a1
pythonplot()

# ╔═╡ 6d3a8414-b746-4a74-abe2-5b400aae8628
theme(:dao)

# ╔═╡ f17ca753-1ffc-4ba2-a178-02741c9711ba
default(size = (600, 400), fontfamily="serif", legend=:outerright)

# ╔═╡ aed436eb-08ba-47af-9379-f0417b2636d6
# load root file
f2nu = ROOTFile("../data/TrackIT/bb2nu/TIT_bb2nu_1M_events.root")

# ╔═╡ 384fedb4-db7e-401d-9b1c-0f83ce6b6070
md"""
# There are 20 variables saved in the root file, we don't need all of them though.
"""

# ╔═╡ 4416fbcd-9b42-482d-9c50-3ce262bcb4c9
begin
	# extract reconstructed vertex positions into dataframe
	df_raw = LazyTree(f2nu, "tree", [
		"y1Reconstructed","z1Reconstructed","y2Reconstructed","z2Reconstructed",
	]) |> DataFrame
	describe(df_raw)
end

# ╔═╡ e7c17231-ab6f-406b-8603-95b6908d921d
md"""
# Analysis:

### The goal here is to create 3 histograms:
1. 1D histogram of $\Delta y$
2. 1D histogram of $\Delta z$
3. 2D histogram of $\Delta d_{zy}$

### For each histogram we fit and extract 68% credible interval (interpreted as sigma), FWHM and rms:
"""

# ╔═╡ 9aa7fb83-e2d4-404c-aea8-a0c8b3740c71
md"""
# First we create 3 new columns: 
* ``dy \equiv y_1 - y_2``
* ``dz \equiv z_1 - z_2``
* ``dyz \equiv \sqrt( dy^2 + dz^2 )``
"""

# ╔═╡ 01de8d24-b00b-4586-af68-c774012fbc94
# transform df to get Δy, Δz, d data
begin
	df_raw.dy = df_raw.y1Reconstructed .- df_raw.y2Reconstructed
	df_raw.dz = df_raw.z1Reconstructed .- df_raw.z2Reconstructed
	df_raw.dyz = @. sqrt( df_raw.dy^2 + df_raw.dz^2)
	df_raw
end

# ╔═╡ cc1e04d1-89bb-47b6-b518-2e785e8d9c61
md"""
## Now that the important data is generated, let's take a look at what the dy, dz, dyz look like:

  1. We can see that for each variable there are outlier values (take a look at the max and compare to mean a median); i.e. `maximum(dy)` = **$( round(Int, maximum(df_raw.dy)) ) mm!!!**
  2. Comparing the mean and median values for dy and dz shows that a larger error is to be expected in dz compared to dy
  
"""

# ╔═╡ c1c26732-7dd3-478a-9ded-8baf49526e3b
describe(df_raw)

# ╔═╡ cd289464-6616-4c05-ae84-42cee12bd4f5
md"""
## Concerning the outlier values: 
- let's take a look at how many outlier there are by studying the quantiles of dy and dz
- [TO DO] we should take a look in flvisualize at these events and figure out what's happening! (Should figure out a way of getting the event ID of the outliers)
"""

# ╔═╡ 6469d868-b94e-473e-8f40-7d17c4edfda9
begin
	alpha = 0.995
	lq, uq = (1.0 - alpha)/2.0, 1.0-(1.0-alpha)/2.0
	@show dy_995q = quantile(df_raw.dy, [lq, uq])
	@show dz_995q = quantile(df_raw.dz, [lq, uq])
end

# ╔═╡ aa42f3ba-fef3-44f4-a346-7a85d5148169
md"""
It can be seen that 99.5% of events have vertices reconstructed that fall well within (a somewhat arbitrary) range of say ~ $\Delta y \in (-240, 240) mm$; $\Delta z \in (-220, 220) mm$ - it should be safe to ignore the crazy outliers and use this the binning of -240:240.
We filter these values.
"""

# ╔═╡ 9dadfd80-f852-4a6b-a0a7-0c4009781d26
# drop outliers:
df = filter(row -> -240.0 < row.dy < 240.0 &&  -220.0 < row.dz < 220.0, df_raw);

# ╔═╡ 07d5bdea-ac7e-4d53-8915-48981535954f
describe(df)

# ╔═╡ 920de5f0-d74a-4eda-ab35-3b4b15b16014
md"""
# Now we can visualize raw data:
"""

# ╔═╡ 1759b08f-ce37-4e87-b0dd-cfe6f05a3989
slider = @bind bw Slider(vcat(1, collect(5:5:20)), default = 5)

# ╔═╡ 52b1a9aa-dfc5-4e25-aa22-8771e0a2ba0c
md"""
binWidth = $bw
"""

# ╔═╡ d57b7942-d584-4c55-9752-63ab0cfaeca6
begin
	h1dy = StatsBase.fit(Histogram{Float32}, df.dy, -240.0:bw:240.0)
	h1dz = StatsBase.fit(Histogram{Float32}, df.dz, -220.0:bw:220.0)
	nothing
end

# ╔═╡ 2852fc0e-d5b4-48d2-bfc7-f315edd84b59
begin
	rms(v::Vector{<:Real}) = std(v; corrected=:false, mean=0)
	nothing
end

# ╔═╡ 26b0aaf0-a9b3-4b65-894f-07d54d19f0a2
with(pythonplot()) do 
	theme(:dao; fontfamily="serif")
	p1 = plot(h1dy, xlabel = "dy / mm", ylabel ="counts / $(bw) mm", label ="", c=1, legend=:outertop, lw = ifelse(bw==1, 0, 0.05*bw), title="TIT dy")
	p2 = plot(h1dz, xlabel = "dz / mm", ylabel ="counts / $(bw) mm", label ="", c=2,legend=:outertop, lw = ifelse(bw==1, 0, 0.05*bw), title="TIT dz")
	p3 = histogram2d(df.dy, df.dz, bins = (-240:bw:240,-240:bw:240), aspect =1, xlabel ="dy", ylabel="dz", title="2D Histogram", colorbar_scale=:log10)

	p = plot(p1, p2, p3, layout=grid(2,2), size = (700,600), thickness_scaling= 1.1, dpi =200)
	
	savefig(p, "../plots/TIT_2nubb.png")
	current()
end

# ╔═╡ f1335a10-8219-46f8-84ee-61c923efa0bd
md"""
## Next we fit the distributions to extract 68% credible interval - to be used as sigma uncertainty (?)
"""

# ╔═╡ 7fc6f77c-c4cd-4eab-817e-078aca48b091
md"""
- We use the following distributions to find which is best: Normal, Cauchy, Laplace distribution 

### We can create the Turing model as:
"""

# ╔═╡ 042caa7a-b290-455e-85bd-0235bab47b80
md"""
To create models in Turing can be as easy as specifying the priors and their distribution (in this case we have 2 parameters: ``\mu`` and ``\sigma``) and the likelihood.
"""

# ╔═╡ bcc2cc51-5b09-4771-9845-39a61244744e
@model function model_Laplace(d)
	# priors
	μ ~ Normal(0, 0.1)
	σ ~ Uniform(1e-5, 1e3)

	# likelihood
	d .~ Laplace(μ, σ)
end

# ╔═╡ 103406a2-97f6-4329-9bc6-b9cc4e333146
md"""
Same for dz.
"""

# ╔═╡ d5d9b097-0abc-4540-b0ed-4202958fecb9
@model function model_Normal(d)
	# priors
	μ ~ Normal(0, 0.1)
	σ ~ Uniform(1e-5, 1e3)

	# likelihood
	d .~ Normal(μ, σ)
end

# ╔═╡ 79ff14f4-a0d6-49c1-b7a5-b62ac1f103d7
@model function model_Cauchy(d)
	# priors
	μ ~ Normal(0, 0.1)
	σ ~ Uniform(1e-5, 1e3)

	# likelihood
	d .~ Cauchy(μ, σ)
end

# ╔═╡ afbad1c3-78a6-4fac-a0bf-616c3e71cdaa
md"""
## Create chain object used in MCMC

> For each model we use the Metropolis-Hastings algorithm with 1000 samples per chain, 1 chain. This is very simple task, so such *loose* conditions should be fine.
"""

# ╔═╡ f8636e78-b7b8-460e-9fac-21640d1293cf
begin
	# sample with metropolis-hastings, we create 10^4 samples per chain, 4 chains:
	chain_dy_L = Turing.sample(model_Laplace(df.dy), MH(), 1000);
	chain_dy_C = Turing.sample(model_Cauchy(df.dy), MH(), 1000);
	chain_dy_N = Turing.sample(model_Normal(df.dy), MH(), 1000);
	nothing # to show no output
end

# ╔═╡ f9d85343-f5f0-4ed1-bfe5-5747e39d677e
begin
	# sample with metropolis-hastings, we create 10^4 samples per chain, 4 chains:
	chain_dz_L = Turing.sample(model_Laplace(df.dz), MH(), 1000);
	chain_dz_C = Turing.sample(model_Cauchy(df.dz), MH(), 1000);
	chain_dz_N = Turing.sample(model_Normal(df.dz), MH(), 1000);
	nothing # to show no output
end

# ╔═╡ 652f3745-349b-4022-a507-12d75ab34cb9
md"""
# We can investigate the fitted distributions: 
"""

# ╔═╡ bbcc79d2-891a-47d3-9c9f-7f5635ac3289
describe(chain_dy_L)

# ╔═╡ 58abb385-6b8f-4b1b-9595-455fb2d6a8fb
describe(chain_dy_C)

# ╔═╡ 9b5a7d31-1008-4270-ba15-d6aca369b1ca
describe(chain_dy_N)

# ╔═╡ 6ac84871-0590-45c2-8593-736b119b7c32
md"""
# Now we extract the mean fitted parameters and create the fitting function
"""

# ╔═╡ eef3674e-efbf-4e23-88cc-6838c5628966
begin
	# now we extract the fitted parameters:
	(mu_dy_L, sigma_dy_L) = mean(chain_dy_L[:,1,:]), mean(chain_dy_L[:,2,:])
	(mu_dy_C, sigma_dy_C) = mean(chain_dy_C[:,1,:]), mean(chain_dy_C[:,2,:])
	(mu_dy_N, sigma_dy_N) = mean(chain_dy_N[:,1,:]), mean(chain_dy_N[:,2,:])
	(mu_dz_L, sigma_dz_L) = mean(chain_dz_L[:,1,:]), mean(chain_dz_L[:,2,:])
	(mu_dz_C, sigma_dz_C) = mean(chain_dz_C[:,1,:]), mean(chain_dz_C[:,2,:])
	(mu_dz_N, sigma_dz_N) = mean(chain_dz_N[:,1,:]), mean(chain_dz_N[:,2,:])
end

# ╔═╡ f978aae0-43c3-43c6-9bf3-f3f5e3aaee6f
begin
	# fit functions for plotting:
	fit_function_L( mu, sigma, x ) = pdf( Laplace(mu, sigma), x )
	fit_function_C( mu, sigma, x ) = pdf( Cauchy(mu, sigma), x )
	fit_function_N( mu, sigma, x ) = pdf( Normal(mu, sigma), x )
end

# ╔═╡ 69964897-74bf-48b7-b88e-604e2607dca1
md"""
# Finally we can visialize the results:
First dy
"""

# ╔═╡ a2615c9b-1c40-4089-8eab-607b5272d999
md"""
#
Now dz
"""

# ╔═╡ 9fb1af4d-f75b-4349-9c0c-e2ed59a26c7c
md"""
# Finally we find that the uncertainty on the CAT vertex reconstruction precision is:
* ``68\% CI \approx (-28, 28)`` mm
* ``68\% CI \approx (-29, 29)`` mm 
"""

# ╔═╡ 34eb944c-93d5-4503-8d93-8c80e244a906
md"""
# Further studies:

- It would be nice to investigate the uncertainty as a function of energy, escape direction w.r.t. foil (see Miro's thesis)
- Maybe it's not the greatest idea in the world to fit the distribution, probably not much useful information is obtained this way, should it be enough to compare data quantiles? 

"""

# ╔═╡ 76988721-9948-4bc3-819a-de1e62827010
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

# ╔═╡ 7581177f-6990-43fd-b734-e5425045b31d
with(gr()) do
	theme(:dao)
	fdyL(x) = step(h1dy.edges[1])*nrow(df)*fit_function_L( mu_dy_L, sigma_dy_L, x)
	fdyC(x) = step(h1dy.edges[1])*nrow(df)*fit_function_C( mu_dy_C, sigma_dy_C, x)
	fdyN(x) = step(h1dy.edges[1])*nrow(df)*fit_function_N( mu_dy_N, sigma_dy_N, x)
	
	pmSigmadyL = quantile(Laplace(mu_dy_L, sigma_dy_L), [0.15865, 0.84135]) # [0.15865, 0.84135] is equal to 68% central quantiles
	CILabeldyL = "68% CI " *L"\in~" * "($(round(pmSigmadyL[1], sigdigits=3)), $(round(pmSigmadyL[2], sigdigits=3))) mm" 

	pmSigmadyC = quantile(Cauchy(mu_dy_C, sigma_dy_C), [0.15865, 0.84135]) # [0.15865, 0.84135] is equal to 68% central quantiles
	CILabeldyC = "68% CI " *L"\in~" * "($(round(pmSigmadyC[1], sigdigits=3)), $(round(pmSigmadyC[2], sigdigits=3))) mm" 
	
	p=plot(1, title = "dy distribution", dpi =200, label="", legend =:outerbottom, legend_column=2, size = (700, 700))
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
		label ="data \nrms = $(round(rms(df.dy), digits=2))", 
		c=2, lw = 2, st=:step, widen =:false
	)
	ylims!(0, 1.3*maximum(h1dy.weights))
	xlabel!("dy [mm]")
	ylabel!("counts / $(bw) mm")
	savefig(p, "../plots/TIT_dy_fit.png" )
	p
end

# ╔═╡ 4d217f3a-340f-4edf-8f25-536ea0a6559e
with(gr()) do
	theme(:dao)
	fdzL(x) = step(h1dz.edges[1])*nrow(df)*fit_function_L( mu_dz_L, sigma_dz_L, x)
	fdzC(x) = step(h1dz.edges[1])*nrow(df)*fit_function_C( mu_dz_C, sigma_dz_C, x)
	fdzN(x) = step(h1dz.edges[1])*nrow(df)*fit_function_N( mu_dz_N, sigma_dz_N, x)
	
	pmSigmadzL = quantile(Laplace(mu_dz_L, sigma_dz_L), [0.15865, 0.84135]) # [0.15865, 0.84135] is equal to 68% central quantiles
	CILabeldzL = "68% CI " *L"\in~" * "($(round(pmSigmadzL[1], sigdigits=3)), $(round(pmSigmadzL[2], sigdigits=3))) mm" 

	pmSigmadzC = quantile(Cauchy(mu_dz_C, sigma_dz_C), [0.15865, 0.84135]) # [0.15865, 0.84135] is equal to 68% central quantiles
	CILabeldzC = "68% CI " *L"\in~" * "($(round(pmSigmadzC[1], sigdigits=3)), $(round(pmSigmadzC[2], sigdigits=3))) mm" 
	
	p=plot(1, title = "dz distribution", dpi =200, label="", legend =:outerbottom, legend_column=2, size = (700, 700))
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
		label ="data \nrms = $(round(rms(df.dz), digits=2))", 
		c=2, lw = 2, st=:step, widen =:false
	)
	ylims!(0, 1.3*maximum(h1dz.weights))
	xlabel!("dy [mm]")
	ylabel!("counts / $(bw) mm")
	savefig(p, "../plots/TIT_dz_fit.png" )
	p
end

# ╔═╡ 0af1a061-05f8-4a70-a8c5-abe54c0e7170
begin
	df.dyAbs = abs.( df.y1Reconstructed .- df.y2Reconstructed )
	df.dzAbs = abs.( df.z1Reconstructed .- df.z2Reconstructed )
end

# ╔═╡ 7b094836-c9e4-4a51-a5f3-a41762942a0e
begin
	h1dyAbs = StatsBase.fit(Histogram, df.dyAbs,0:bw:200)
	h1dzAbs = StatsBase.fit(Histogram, df.dzAbs,0:bw:200)
end

# ╔═╡ 3d2bebc8-1b2a-4a82-bcff-68d921f7db9d
@model function modelExp(d)
	lambda ~ Uniform(1e-3, 1e3)

	d .~ Exponential(lambda)
end

# ╔═╡ bcdf7c56-a0f9-43f4-b3ff-6d48393ff124
begin
	chain_dy_E = Turing.sample(modelExp(df.dyAbs), NUTS(0.65), 1000)
	chain_dz_E = Turing.sample(modelExp(df.dyAbs), NUTS(0.65), 1000)
end

# ╔═╡ 5e53595d-4c89-4ce7-b1f4-40d72723f5b1
plot(chain_dz_E)

# ╔═╡ 9f09aa68-39fc-4e5f-9c5d-23745bab288c
begin
	l_dyAbs = mean(chain_dy_E[:lambda])
	l_dzAbs = mean(chain_dz_E[:lambda])
end

# ╔═╡ c0b85bdc-e7d4-4197-ad34-8a3806582db4
pdyAbs=let
	gr()
	theme(:dao, fontfamily="serif")
	p=plot(widen=:false,label="",title=L"dy \equiv |y1 - y2|", legend=:outerright, legend_column=1, thickness_scaling=1.2)
	fExp(x) = step(h1dyAbs.edges[1])*nrow(df)*pdf(Exponential(l_dyAbs),x)
	q1 = round(quantile(Exponential(l_dyAbs), 0.6827), sigdigits=3)
	r2Exp = r2(fExp, h1dyAbs)
	plot!(h1dyAbs, label="data", legendtitle="Exponential fit: \n"* L"p(x)=\lambda e^{-\lambda x}",st=:step, fa=0.3, f=0, lw=3)
	plot!(0:1:200, x->fExp(x), label="λ = $(round(l_dyAbs,sigdigits=3))\nr2 = $(round(r2Exp, sigdigits=3) )\n68% CI: (0,$q1)mm", lw=3,size(1200,600),)
	xlabel!("dy [mm]")
	ylabel!("counts / $bw mm")
	savefig(p, "../plots/TIT_dyAbs_fit.png" )
	
	p
end

# ╔═╡ 0397b7e7-8fe9-49a7-b1f5-2473f7e9d9b2
pdzAbs=let
	gr()
	theme(:dao, fontfamily="serif")
	p=plot(widen=:false,label="",title=L"dz \equiv |z1 - z2|", legend=:outerright, legend_column=1, thickness_scaling=1.2)
	fExp(x) = step(h1dzAbs.edges[1])*nrow(df)*pdf(Exponential(l_dzAbs),x)
	q1 = round(quantile(Exponential(l_dzAbs), 0.6827), sigdigits=3)
	
	r2Exp = r2(fExp, h1dzAbs)
	plot!(h1dzAbs, label="data", legendtitle="Exponential fit: \n"* L"p(x)=\lambda e^{-\lambda x}",st=:step, fa=0.3, f=0, lw=3)
	plot!(0:1:200, x->fExp(x), label="λ = $(round(l_dzAbs,sigdigits=3))\nr2 = $(round(r2Exp, sigdigits=3) )\n68% CI: (0,$q1)mm", lw=3)
	xlabel!("dz [mm]")
	ylabel!("counts / $bw mm")
	savefig(p, "../plots/TIT_dzAbs_fit.png" )
	p
end

# ╔═╡ Cell order:
# ╟─8d855d70-1557-49ae-b995-bd8d4cf7a455
# ╟─c6572c84-de77-40e1-9ca4-d9e72f73ac76
# ╟─e4e23a9c-cfbd-467b-a3a5-f2477fc66531
# ╠═6b9feb78-371e-446d-83ca-cc0da8c65526
# ╠═c9e89e94-1185-42c9-bab2-02628fb7b481
# ╠═d3370197-b6c6-436c-8a54-30118a741278
# ╠═9ef3374f-bd99-4a41-893a-e231329397a1
# ╠═6d3a8414-b746-4a74-abe2-5b400aae8628
# ╠═f17ca753-1ffc-4ba2-a178-02741c9711ba
# ╠═aed436eb-08ba-47af-9379-f0417b2636d6
# ╟─384fedb4-db7e-401d-9b1c-0f83ce6b6070
# ╠═4416fbcd-9b42-482d-9c50-3ce262bcb4c9
# ╠═e7c17231-ab6f-406b-8603-95b6908d921d
# ╠═9aa7fb83-e2d4-404c-aea8-a0c8b3740c71
# ╠═01de8d24-b00b-4586-af68-c774012fbc94
# ╠═cc1e04d1-89bb-47b6-b518-2e785e8d9c61
# ╠═c1c26732-7dd3-478a-9ded-8baf49526e3b
# ╠═cd289464-6616-4c05-ae84-42cee12bd4f5
# ╠═6469d868-b94e-473e-8f40-7d17c4edfda9
# ╠═aa42f3ba-fef3-44f4-a346-7a85d5148169
# ╠═9dadfd80-f852-4a6b-a0a7-0c4009781d26
# ╠═07d5bdea-ac7e-4d53-8915-48981535954f
# ╟─920de5f0-d74a-4eda-ab35-3b4b15b16014
# ╠═1759b08f-ce37-4e87-b0dd-cfe6f05a3989
# ╟─52b1a9aa-dfc5-4e25-aa22-8771e0a2ba0c
# ╠═d57b7942-d584-4c55-9752-63ab0cfaeca6
# ╟─2852fc0e-d5b4-48d2-bfc7-f315edd84b59
# ╠═26b0aaf0-a9b3-4b65-894f-07d54d19f0a2
# ╟─f1335a10-8219-46f8-84ee-61c923efa0bd
# ╠═7fc6f77c-c4cd-4eab-817e-078aca48b091
# ╠═042caa7a-b290-455e-85bd-0235bab47b80
# ╠═bcc2cc51-5b09-4771-9845-39a61244744e
# ╟─103406a2-97f6-4329-9bc6-b9cc4e333146
# ╠═d5d9b097-0abc-4540-b0ed-4202958fecb9
# ╠═79ff14f4-a0d6-49c1-b7a5-b62ac1f103d7
# ╟─afbad1c3-78a6-4fac-a0bf-616c3e71cdaa
# ╠═f8636e78-b7b8-460e-9fac-21640d1293cf
# ╠═f9d85343-f5f0-4ed1-bfe5-5747e39d677e
# ╟─652f3745-349b-4022-a507-12d75ab34cb9
# ╠═bbcc79d2-891a-47d3-9c9f-7f5635ac3289
# ╠═58abb385-6b8f-4b1b-9595-455fb2d6a8fb
# ╠═9b5a7d31-1008-4270-ba15-d6aca369b1ca
# ╟─6ac84871-0590-45c2-8593-736b119b7c32
# ╠═eef3674e-efbf-4e23-88cc-6838c5628966
# ╠═f978aae0-43c3-43c6-9bf3-f3f5e3aaee6f
# ╟─69964897-74bf-48b7-b88e-604e2607dca1
# ╠═7581177f-6990-43fd-b734-e5425045b31d
# ╟─a2615c9b-1c40-4089-8eab-607b5272d999
# ╠═4d217f3a-340f-4edf-8f25-536ea0a6559e
# ╠═9fb1af4d-f75b-4349-9c0c-e2ed59a26c7c
# ╟─34eb944c-93d5-4503-8d93-8c80e244a906
# ╠═76988721-9948-4bc3-819a-de1e62827010
# ╠═0af1a061-05f8-4a70-a8c5-abe54c0e7170
# ╠═7b094836-c9e4-4a51-a5f3-a41762942a0e
# ╠═3d2bebc8-1b2a-4a82-bcff-68d921f7db9d
# ╠═bcdf7c56-a0f9-43f4-b3ff-6d48393ff124
# ╠═5e53595d-4c89-4ce7-b1f4-40d72723f5b1
# ╠═9f09aa68-39fc-4e5f-9c5d-23745bab288c
# ╠═c0b85bdc-e7d4-4197-ad34-8a3806582db4
# ╠═0397b7e7-8fe9-49a7-b1f5-2473f7e9d9b2
