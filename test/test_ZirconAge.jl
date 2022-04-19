using Test
using GeoParams, LinearAlgebra, DelimitedFiles
@testset "ZirconAge.jl" begin

# declare constant variables
s2y 			= 365.0*24.0*3600.0 							# second to year
ZirconData  	= ZirconAgeData();		# default data

# Parameters for Zircon statistical analysis
n_analyses		= 300											# number of synthetic zircon analyses

# read input file (this step should be skipped, when using tracers)		
filename 		= 	"Data/Tt15_st_Bl_rad5_0.0126.txt"				# path to input file 
Tt_paths		=	readdlm(filename,',')
n_paths			= 	size(Tt_paths,2)-1							# -1 ro remove first column that corresponds to the time
time_years 		= 	Tt_paths[:,1]./s2y							# time fron seconds to years

# Here a matrix is constructed that only contains the temperature information:
Tt_paths_Temp 	= 	Tt_paths[:,2:end]	
	
# first: test case in which we provide the Tt-path as matrix 
prob, ages_eruptible, number_zircons, T_av_time, T_sd_time = compute_zircons_Ttpath(time_years, Tt_paths_Temp, ZirconData=ZirconData)
time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average  = zircon_age_PDF(ages_eruptible, number_zircons, bandwidth=1e5, n_analyses=n_analyses)

# add tests to check that results are consistent
@test sum(number_zircons[:,200]) == 40479.0
@test sum(number_zircons)==5.411985e6
@test  prob[100] ≈ 5.5432526143365145e-6

# A second way to generate the input is having Vector{Vector} for both time and Tt-path. 
time_years_vecs = Vector{Vector}(undef,size(Tt_paths_Temp,2));
Tt_paths_Temp_vecs		= Vector{Vector}(undef,size(Tt_paths_Temp,2));
for i=1:size(Tt_paths_Temp,2)
	if i<400
		time_years_vecs[i] 		= time_years
		Tt_paths_Temp_vecs[i] 	= Tt_paths_Temp[:,i]
	else
		# slightly adjust the size of the vectors
		time_years_vecs[i] 		= time_years[3:end]
		Tt_paths_Temp_vecs[i] 	= Tt_paths_Temp[3:end,i]
	end

end

# the calling routine is the same, but we get one additional output vector:
time_years1, prob1, ages_eruptible1, number_zircons1, T_av_time1, T_sd_time1  = compute_zircons_Ttpath(time_years_vecs, Tt_paths_Temp_vecs)

# add tests to check that results are consistent
@test sum(number_zircons1[:,200]) == 40479.0
@test sum(number_zircons1)==5.411985e6
@test  prob1[100] ≈ 5.5432526143365145e-6


# Do the same but with a single routine that also returns the PDF's 
# Note that given the randomness, you'll always get different results
time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average, time_years, 
	prob2, ages_eruptible, number_zircons2, T_av_time, T_sd_time = compute_zircon_age_PDF(time_years_vecs, Tt_paths_Temp_vecs)

@test sum(number_zircons2[:,200]) == 40479.0
@test sum(number_zircons2)==5.411985e6
@test  prob2[100] ≈ 5.5432526143365145e-6

# plotting routine (commented, but this is how you'd call it):
# Plot_ZirconAge_PDF(time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average)


# Julia translation of the of the R-script provided by Gregor Weber & Luca Caricchi 
# Used in the publication "Estimating the current size and state of subvolcanic magma reservoirs"

# 15/04/2022


using DelimitedFiles
using Loess
using Statistics
using StatsBase, KernelDensity, Loess
using Plots
using Test
using Parameters

"""

	Declare function to calculate zircon fraction as function of a temperature profile
"""
function zircon_fraction(T::AbstractArray{_T}, max_x_zr::_T) 	where _T
	A = (1.62.-1.8*(10^4)*exp.((-10^4)./(T .+ 273.15))).*max_x_zr 
	A[A .<= 0.0] .= 0.0

	return A
end

@with_kw_noshow struct ZirconAgeData
	Tsat::Float64 			= 825.0		# Maximum zircon saturation temperature [C]
	Tmin::Float64 			= 690.0 	# [C] Minimum zircon saturation Temperature [C]
	Tsol::Float64 			= 690.0		# [C] Solidus temperature [C]
	Tcal_max::Float64		= 800.0		# max temperature to calculate zircon fraction
	Tcal_step::Float64 		= 1.0		# temperature step to caclulate zircon fraction (resolution of Zircon saturation curve discretization)
	max_x_zr::Float64 		= 0.001		# max fraction zircons at solidus
	zircon_number::Int64    = 100.0		# number of required zircons 
	time_zr_growth::Float64 = 0.7e6		# Minimum time within T saturation range (This is what the method used in the R script, boils down too)
	# -> remain in the Zr saturation zone more than 1/3 of the time the Tt path with the longest time in the saturation zone
	
end

"""
	Loess fit of zircon number over saturation range
"""
function loess_fit_zircon_sat(ZirconData::ZirconAgeData)
	@unpack Tsol, Tcal_step, Tcal_max, Tsat, max_x_zr, zircon_number, Tmin = ZirconData								

	# Get cumulative Zircon fraction (after Tierney et al., 2016; Geology)
	T 				= range(Tsol, stop = Tcal_max, step = Tcal_step  )
	Tfit 			= range(Tsol, stop = Tsat, 	   length = length(T))

	x_zircon 		= zircon_fraction(T,max_x_zr)			# compute zircon fraction 
	n_zircon 		= zero(x_zircon)						# number of zircons

	n_zircon[2:end] = -diff(x_zircon)						# differentiate zircon fraction
	n_zircon[1] 	= n_zircon[2]							# set first value

	n_zircon		= ceil.(( (n_zircon*zircon_number)/maximum(n_zircon)).-minimum(floor.((n_zircon*zircon_number)/maximum(n_zircon))))

	# fit n_zircon_N with loess regression
	n_zircon_N_fit 	= loess(Tfit, n_zircon, span=1.0)

	return n_zircon_N_fit
end

"""

This employs a loess fit function to compute the number of zircons for each of the Tt paths
"""
function  compute_number_zircons!(n_zr::AbstractArray{_T,N}, Tt_paths_Temp::AbstractArray{_T,N}, ZirconData::ZirconAgeData) where {_T,N}
	@unpack Tmin, Tsat, Tsol = ZirconData
	
	n_zircon_N_fit  = 	loess_fit_zircon_sat(ZirconData)			# loess fit through zircon saturation
	
	for i in 1:size(Tt_paths_Temp,2)
		for j in 1:size(Tt_paths_Temp,1)
			T = Tt_paths_Temp[j,i]
			if (T > Tsol) & (T < Tsat)
				# The following line takes a lot of time, because Loess is not type-stable. 
				# There is an open PR in the package that may fix it
				dat::_T      = Loess.predict(n_zircon_N_fit, T)		
				n_zr[j,i]    = floor(dat)
			else
				n_zr[j,i]    = 0.0
			end
		end
	end

	nothing
end

"""
	prob, ages_eruptible, number_zircons, T_av_time, T_sd_time  compute_zircons_Ttpath(time_years::AbstractArray{Float64,1}, Tt_paths_Temp::AbstractArray{Float64,2}; ZirconData::ZirconAgeData)

This computes the number of zircons produced from a series of temperature-time path's. 
The Tt-paths are stored in a 2D matrix `Tt_paths_Temp` with rows being the temperature at time `time_years`.

Input:
====
- `time_years` : vector of length `nt` with the time in years (since the beginning of the simulation) of the points provided
- `Tt_paths_Temp` : array of size `(nt,npaths)`` with the temperature of every path.

Output:
- `prob` : a vector that gives the relative probability that a zircon with a given age exists
- `ages_eruptible` : age of eruptble magma
- `number_zircons` : 2D array of size `(nt,)`
- `T_av_time`: vector of size `nt` that contains the average T of the paths
- `T_sd_time`: vector of size `nt` that contains the standard deviation of the T of the paths

This routine was developed based on an R-routine provided as electronic supplement in the paper:
- Weber, G., Caricchi, L., Arce, J.L., Schmitt, A.K., 2020. Determining the current size and state of subvolcanic magma reservoirs. Nat Commun 11, 5477. https://doi.org/10.1038/s41467-020-19084-2

"""
function compute_zircons_Ttpath(time_years::AbstractArray{Float64,1}, Tt_paths_Temp::AbstractArray{Float64,2}; ZirconData::ZirconAgeData = ZirconAgeData())

	@unpack Tmin, Tsat, Tsol, time_zr_growth = ZirconData
	
	Tt_paths_Temp1 = copy(Tt_paths_Temp)
	

	Δt 				=	diff(time_years)[1]							# timestep [yrs]
	time_er_min 	= 	maximum(time_years)							# backward count

	# find all the Tt paths that go through the zircon saturation range & are at the end of the path still below Tsat (otherwise Zr are not yet crystallized)
	ID_col_er 		= findall( (maximum(Tt_paths_Temp1,dims=1).>Tmin) .& (Tt_paths_Temp1[end,:]' .< Tsat))
	n_zr			= zero(Tt_paths_Temp1)
	compute_number_zircons!(n_zr, Tt_paths_Temp1, ZirconData)		# computes the number of zircons for every path

	# find the number of timesteps for every path, during which the temperature is > Tmin and < Tsat
	length_trace 	= Vector{Float64}(undef,length(ID_col_er))
	for i in 1:length(ID_col_er)
		length_trace[i]    = length( findall( (Tt_paths_Temp1[:,ID_col_er[i][2]] .> Tmin) .& (Tt_paths_Temp1[:,ID_col_er[i][2]] .< Tsat)) )
	end
	
	# the next several lines can likely be achieved in a more elegant way...
	id				= findall( length_trace .== maximum(length_trace)) 
	ID_col_lgst_tr 	= ID_col_er[ id[1] ][2]
	
	id 				= findall( Tt_paths_Temp1[:,ID_col_lgst_tr] .< Tsat) 
	VALID_min_time 	= findmin( Tt_paths_Temp1[id,ID_col_lgst_tr]) 
	ID_min_time		= VALID_min_time[2]
	
	max_age_spread	= maximum(length_trace)*Δt					# This is defined among all selected paths
	
	T_av_time_1 	= replace!(Tt_paths_Temp1, 0.0 => NaN)
	T_av_time 		= zeros(size(T_av_time_1,1))
	T_sd_time 		= zeros(size(T_av_time_1,1))
	
	# get the average temperature of the Tt paths and the standard deviation
	for i in 1:size(Tt_paths_Temp1,1)
		T_av_time[i]	= mean(filter(!isnan, T_av_time_1[i,:]))
		T_sd_time[i]	= std(filter( !isnan, T_av_time_1[i,:]))
	end

	# I clarified the R function because the minimum step length to grow a zircon is simply a ratio of the maximum trace between Tmin and Tsol
	# This makes sense as we only deal with fractions here. Because no mass is provided the real zircon size cannot possibly be determined
	min_step_n		= floor( (time_zr_growth/max_age_spread)*(max_age_spread/Δt) )
	
	# find the Tt paths that have a number of timesteps in the saturation range greater than the defined min_step_n ()
	# this is to mimic that it takes some time to grow zircons
	id				= findall( length_trace .> min_step_n) 
	ID_col_er_1		= getindex.(ID_col_er[id], [2])
	
	int_zr_sat		= collect(Float64,  ID_min_time:1.0:(time_er_min/Δt)-min_step_n)
	int_zr_sat		= floor.(Int64,int_zr_sat)
	
	T_av_time_slct 	= Vector{Float64}(undef,length(int_zr_sat)-1) .= 0.0
	
	for i in 1:(size(int_zr_sat,1)-1)
		id2				= ID_col_er_1[ findall( (Tt_paths_Temp1[int_zr_sat[i],ID_col_er_1[:]] .> Tmin) .& (Tt_paths_Temp1[int_zr_sat[i],ID_col_er_1[:]] .< Tsat)) ]
		if isempty(id2) == true
			T_av_time_slct[i] = NaN
		else
			T_av_time_slct[i] = median(filter(!isnan, Tt_paths_Temp1[int_zr_sat[i],id2]))
		end
	end
	
	replace!(Tt_paths_Temp1, NaN => 0.0)
	ID_col_er			= getindex.(ID_col_er, [2])
	for i in 1:length(ID_col_er)
		k 				= maximum(findall( (Tt_paths_Temp1[:,ID_col_er[i]]) .== 0.0 ))
		Tt_paths_Temp1[1:k,ID_col_er[i]] .= 0.0
	end

	zr_select			= zero(Tt_paths_Temp1)
	zr_select[Tt_paths_Temp1 .> 0.0] .= 1.0	
	n_zrc2_0			= zr_select.*n_zr						# filters out those Tt path that are still >Tsat @ the end 
	number_zircons      = n_zrc2_0[:,ID_col_er_1];
	n_measurable_ages 	= sum(number_zircons, dims=2)	
	sz 					= size(number_zircons,1)
	ages_eruptible		= collect(Float64,  1.0:Δt:sz*Δt)
	
	# probability that a certain zircon is sampled, dependens on how many of a given age ara available:
	prob 				= n_measurable_ages/sum(n_measurable_ages)
	prob 				= prob[:,1]

	return prob, ages_eruptible, number_zircons, T_av_time, T_sd_time
end


"""
	time_years, prob, ages_eruptible, number_zircons, T_av_time, T_sd_time  = compute_zircons_Ttpath(time_years::Vector{Vector{Float64}}, Tt_paths_Temp::Vector{Vector{Float64}}; ZirconData::ZirconAgeData = ZirconAgeData())

This accepts Vector{Vector} as input for time and temperature of each Tt-path. Here, the length of the vector can be variable between different points.

Internally, we interpolate this into a 2D matrix and a longer vector that includes all paths and a single vector with times 

"""
function compute_zircons_Ttpath(time_years_vecs::Vector{Vector}, Tt_paths_Temp_vecs::Vector{Vector}; ZirconData::ZirconAgeData = ZirconAgeData())

	# Create a single vector with the time values:
	time_years = unique(reduce(vcat,unique(time_years_vecs))); 

	# Add the vectors to an array with T values
	Tt_paths_Temp = zeros(Float64,length(time_years), length(Tt_paths_Temp_vecs))
	for i in eachindex(Tt_paths_Temp_vecs)
		istart = findall(time_years .== time_years_vecs[i][1])
		iend   = findall(time_years .== time_years_vecs[i][end])
		Tt_paths_Temp[istart[1]:iend[1], i] .= Tt_paths_Temp_vecs[i]
	end

	# call main routine
	prob, ages_eruptible, number_zircons, T_av_time, T_sd_time = compute_zircons_Ttpath(time_years, Tt_paths_Temp, ZirconData=ZirconData)

	# return, including time_years
	return time_years, prob, ages_eruptible, number_zircons, T_av_time, T_sd_time 

end

"""
	zircon_age_PDF(ages_eruptible::AbstractArray{Float64,1}, number_zircons::AbstractArray{Float64,2}, bandwidth=1e5, n_analyses=300)

Compute probability density functions for zircon age path's describes in `number_zircons` with age `ages_eruptible` (both computed ).
`bandwidth` is the smoothening window of the resulting curves (in years), whereas `n_analyses` are the number of analyses done. 	

"""
function zircon_age_PDF(ages_eruptible::AbstractArray{Float64,1}, number_zircons::AbstractArray{Float64,2}; bandwidth=1e5, n_analyses=300)

	# compute PDF for each of the zircon Tt-paths:
	PDF_zircons = []
	time_Ma = [];
	for i in 1:size(number_zircons,2)
		n_meas 			 = number_zircons[:,i]
		px	  			 = n_meas/sum(n_meas)		# probability to have a certain age
		
		# random numbers selected according to the probability
		smp				 = sample( (maximum(ages_eruptible) .- ages_eruptible)/1e6, Weights(px), n_analyses, replace=true)
		y 				 = kde(smp, bandwidth=bandwidth/1e6)

		# store data
		push!(PDF_zircons, 	y.density)
		push!(time_Ma, 		y.x)
	end

	n_measurable_ages   = sum(number_zircons, dims=2)
	pxAv	  			= n_measurable_ages[:,1]./sum(n_measurable_ages[:,1])
	smpAv				= sample( (maximum(ages_eruptible) .- ages_eruptible)/1e6, Weights(pxAv), n_analyses, replace=true)
	yAv 				= kde(smpAv, bandwidth=bandwidth/1e6)
	time_Ma_average     = Vector(yAv.x);
	PDF_zircon_average  = Vector(yAv.density);
	
	return time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average
end


"""
	time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average, time_years, prob, ages_eruptible, number_zircons, T_av_time, T_sd_time = compute_zircon_age_PDF(time_years_vecs::Vector{Vector}, Tt_paths_Temp_vecs::Vector{Vector}; ZirconData::ZirconAgeData = ZirconAgeData(), bandwidth=bandwidth, n_analyses=300)

This computes the PDF (probability density function) with zircon age data from Vectors with Tt-paths	

"""
function compute_zircon_age_PDF(time_years_vecs::Vector{Vector}, Tt_paths_Temp_vecs::Vector{Vector}; ZirconData::ZirconAgeData = ZirconAgeData(), bandwidth=1e5, n_analyses=300)

	# Compute the probability that a zircon of certain age is sampled:
	time_years, prob, ages_eruptible, number_zircons, T_av_time, T_sd_time = compute_zircons_Ttpath(time_years_vecs, Tt_paths_Temp_vecs, ZirconData=ZirconData);

	# Use this to compute PDF curves: 
	time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average  = zircon_age_PDF(ages_eruptible, number_zircons, bandwidth=bandwidth, n_analyses=n_analyses)

	return time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average, time_years, prob, ages_eruptible, number_zircons, T_av_time, T_sd_time

end


"""
	plt = Plot_ZirconAge_PDF(time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average)

Creates a plot of the Zircon Age probability density function from the parameters in a simulation
"""
function Plot_ZirconAge_PDF(time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average)

	plt = Plots.plot(time_Ma[1], PDF_zircons[1], color=:lightgray,linewidth=0.1, 
				xlabel="Time [Ma]", ylabel="probability []", title = "Zircon age probability distribution", legend=:none)
	for i in 2:length(PDF_zircons)
		plt = Plots.plot!(time_Ma[i], PDF_zircons[i], color=:lightgray,linewidth=0.1)
	end
	Plots.plot!(time_Ma_average, PDF_zircon_average, color=:black,linewidth=2.)
	
	display(plt)

	return plt
end

# declare constant variables
s2y 			= 365.0*24.0*3600.0 							# second to year
ZirconData  	= ZirconAgeData();		# default data

# Clarified way to define the zircon growth as a minimum time within T saturation range (This is what the method used in the R script, boils down too)
# -> remain in the Zr saturation zone more than 1/3 of the time the Tt path with the longest time in the saturation zone
time_zr_growth 	= 0.7e6											# ref should be 0.7

# Parameters for Zircon statistical analysis
n_analyses		= 300											# number of synthetic zircon analyses

# read input file (this step should be skipped, when using tracers)		
filename 		= 	"Tt15_st_Bl_rad5_0.0126.txt"				# path to input file 
Tt_paths		=	readdlm(filename,',')
n_paths			= 	size(Tt_paths,2)-1							# -1 ro remove first column that corresponds to the time
time_years 		= 	Tt_paths[:,1]./s2y							# time fron seconds to years

# Here a matrix is constructed that only contains the temperature information:
Tt_paths_Temp 	= 	Tt_paths[:,2:end]	
	
# first: test case in which we provide the Tt-path as matrix 
prob, ages_eruptible, number_zircons, T_av_time, T_sd_time = compute_zircons_Ttpath(time_years, Tt_paths_Temp, ZirconData=ZirconData)
time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average  = zircon_age_PDF(ages_eruptible, number_zircons, bandwidth=1e5, n_analyses=n_analyses)

# add tests to check that results are consistent
@test sum(number_zircons[:,200]) == 40479.0
@test sum(number_zircons)==5.411985e6
@test  prob[100] ≈ 5.5432526143365145e-6

# A second way to generate the input is having Vector{Vector} for both time and Tt-path. 
time_years_vecs = Vector{Vector}(undef,size(Tt_paths_Temp,2));
Tt_paths_Temp_vecs		= Vector{Vector}(undef,size(Tt_paths_Temp,2));
for i=1:size(Tt_paths_Temp,2)
	if i<400
		time_years_vecs[i] 		= time_years
		Tt_paths_Temp_vecs[i] 	= Tt_paths_Temp[:,i]
	else
		# slightly adjust the size of the vectors
		time_years_vecs[i] 		= time_years[3:end]
		Tt_paths_Temp_vecs[i] 	= Tt_paths_Temp[3:end,i]
	end

end

# the calling routine is the same, but we get one additional output vector:
time_years1, prob1, ages_eruptible1, number_zircons1, T_av_time1, T_sd_time1  = compute_zircons_Ttpath(time_years_vecs, Tt_paths_Temp_vecs)

# add tests to check that results are consistent
@test sum(number_zircons1[:,200]) == 40479.0
@test sum(number_zircons1)==5.411985e6
@test  prob1[100] ≈ 5.5432526143365145e-6


# Do the same but with a single routine that also returns the PDF's 
# Note that given the randomness, you'll always get different results
time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average, time_years, 
	prob2, ages_eruptible, number_zircons2, T_av_time, T_sd_time = compute_zircon_age_PDF(time_years_vecs, Tt_paths_Temp_vecs)

@test sum(number_zircons2[:,200]) == 40479.0
@test sum(number_zircons2)==5.411985e6
@test  prob2[100] ≈ 5.5432526143365145e-6

	
# plotting routine:
Plot_ZirconAge_PDF(time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average)


#=	
    # Plot Zircon age probability distribution
	# these are the ploting routines using Makie, which is currently not a dependency of GeoParams (but may become one @ some stage)

	f = Figure()
	Axis(f[1, 1], xlabel = "Age [Myr]", ylabel = "Kernel density [ ]", title = "Zircon age probability distribution")
	for i in 1:length(PDF_zircons)
		CairoMakie.lines!(time_Ma[i]/1e6, PDF_zircons[i], color="gray66",linewidth=0.25)
	end
	CairoMakie.lines!(time_Ma_average/1e6, PDF_zircon_average, color="grey0",linewidth=2.)
	CairoMakie.xlims!(-1e5/1e6,1.5e6/1e6)
	save("Zircon_probability_plot_800kyrs.png",f)

	# plot evolution of the average temperature since the onset of magma injection
	f = Figure()
	Axis(f[1, 1], ylabel = "Temperature [°C]", xlabel = "Time [Myr]", title = "Evolution of the average temperature of the magmatic system")
	CairoMakie.lines!(time_years./1e6, T_av_time, color="grey0",linewidth=0.5)
	f
	save("Average_T.png",f)
=#




end

