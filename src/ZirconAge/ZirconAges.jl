module ZirconAges

# Julia translation of a R-script provided by Gregor Weber & Luca Caricchi to compute zircon ages from Tt paths
# Used in the Nat. Comm. publication "Estimating the current size and state of subvolcanic magma reservoirs"
# 15/04/2022, Nico Riel & Boris Kaus

import Base.Threads
using Parameters
using Loess, Statistics, StatsBase, KernelDensity, Loess

export  ZirconAgeData, 
        compute_zircon_age_PDF,  compute_zircons_Ttpath,  # calculation routines
        zircon_age_PDF, 
        compute_zircons_convert_vecs2mat
    
"""
    Declare function to calculate zircon fraction as function of a temperature profile
"""
function zircon_fraction(T::AbstractArray{_T}, max_x_zr::_T) 	where _T
    A = (1.62.-1.8*(10^4)*exp.((-10^4)./(T .+ 273.15))).*max_x_zr 
    A[A .<= 0.0] .= 0.0

    return A
end

"""
    ZirconAgeData

Struct that holds default parameters for the calculations
"""
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
    Threads.@threads for i in 1:size(Tt_paths_Temp,2)
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
function compute_zircons_Ttpath(time_years::AbstractArray{_T,1}, Tt_paths_Temp::AbstractArray{_T,2}; ZirconData::ZirconAgeData = ZirconAgeData()) where _T

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
    if isempty(id)
        max_Ptpath = maximum(length_trace)*Δt 
        error("I don't have a single Pt-path that is sufficiently long within time_zr_growth (=$(time_zr_growth) yrs). 
            The Longest Pt-path I have is $(max_Ptpath) years. 
            Decrease this value within the ZirconDataAge struct with ZirconData=ZirconAgeData(time_zr_growth=0.1e6) & rerun.")
    end
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
        ind = findall( (Tt_paths_Temp1[:,ID_col_er[i]]) .== 0.0 );
        if ~isempty(ind)
            k 				= maximum(ind)
            Tt_paths_Temp1[1:k,ID_col_er[i]] .= 0.0
        end
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
function compute_zircons_Ttpath(time_years_vecs::Vector{Vector{_T}}, Tt_paths_Temp_vecs::Vector{Vector{_T}}; ZirconData::ZirconAgeData = ZirconAgeData()) where _T

    # convert to a vector with time & matrix with T values at every timestep
    time_years, Tt_paths_Temp = compute_zircons_convert_vecs2mat(time_years_vecs, Tt_paths_Temp_vecs)

    # call main routine
    prob, ages_eruptible, number_zircons, T_av_time, T_sd_time = compute_zircons_Ttpath(time_years, Tt_paths_Temp, ZirconData=ZirconData)

    # return, including time_years
    return time_years, prob, ages_eruptible, number_zircons, T_av_time, T_sd_time 

end


"""
    time_years, Ttpaths_mat = compute_zircons_convert_vecs2mat(time_years_vecs::Vector{Vector{Float64}}, Tt_paths_Temp_vecs::Vector{Vector{Float64}})     

This converts a vector with Vectors contain time and temperature path's to a single time vector and a matrix that combines all vectors
"""
function compute_zircons_convert_vecs2mat(time_years_vecs::Vector{Vector{_T}}, Tt_paths_Temp_vecs::Vector{Vector{_T}}) where _T
    
    # Create a single vector with the time values:
    time_years = unique(reduce(vcat,unique(time_years_vecs))); 

    # Add the vectors to an array with T values
    Tt_paths_mat = zeros(_T,length(time_years), length(Tt_paths_Temp_vecs))
    for i in eachindex(Tt_paths_Temp_vecs)
        istart = findall(time_years .== time_years_vecs[i][1])
        iend   = findall(time_years .== time_years_vecs[i][end])
        Tt_paths_mat[istart[1]:iend[1], i] .= Tt_paths_Temp_vecs[i]
    end

    return time_years, Tt_paths_mat

end




"""
    zircon_age_PDF(ages_eruptible::AbstractArray{Float64,1}, number_zircons::AbstractArray{Float64,2}, bandwidth=1e5, n_analyses=300)

Compute probability density functions for zircon age path's describes in `number_zircons` with age `ages_eruptible` (both computed ).
`bandwidth` is the smoothening window of the resulting curves (in years), whereas `n_analyses` are the number of analyses done. 	

"""
function zircon_age_PDF(ages_eruptible::AbstractArray{_T,1}, number_zircons::AbstractArray{_T,2}; bandwidth=1e5, n_analyses=300) where _T

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
function compute_zircon_age_PDF(time_years_vecs::Vector{Vector{_T}}, Tt_paths_Temp_vecs::Vector{Vector{_T}}; ZirconData::ZirconAgeData = ZirconAgeData(), bandwidth=1e5, n_analyses=300) where _T
    
    # Compute the probability that a zircon of certain age is sampled:
    time_years, prob, ages_eruptible, number_zircons, T_av_time, T_sd_time = compute_zircons_Ttpath(time_years_vecs, Tt_paths_Temp_vecs, ZirconData=ZirconData);
    
    # Use this to compute PDF curves: 
    time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average  = zircon_age_PDF(ages_eruptible, number_zircons, bandwidth=bandwidth, n_analyses=n_analyses)

    return time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average, time_years, prob, ages_eruptible, number_zircons, T_av_time, T_sd_time

end

end