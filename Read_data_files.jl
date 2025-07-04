# import libraries
using CSV
using DataFrames
using UnPack

function str_convert(rho, tau, gamma, alpha)
    """Function to convert dispersal parameters to a string format without decimal points."""
    rho_str = replace(string(rho), "." => "")      # rho
    tau_str = replace(string(tau), "." => "")      # tau
    gamma_str = replace(string(gamma), "." => "")  # gamma
    alpha_str = replace(string(alpha), "." => "")  # alpha

    return rho_str, tau_str, gamma_str, alpha_str  # return converted strings
end

function read_data(params,JOB_ID)
    """Reads in data files  for the total wildtype and Wolbachia-infected populations as well as free mosquito 
    populations and the tracked households based on the provided simulation parameters and job ID."""
    @unpack rho, tau, gamma, alpha, H, rel_t, rel_size, rel_type = params # unpack parameters
    # Convert parameters to string format for file naming
    rho_str, tau_str, gamma_str, alpha_str = str_convert(rho, tau, gamma, alpha)
    ## Read in the total wildtype and Wolbachia-infected mosquito populations data, 
    ## as well as free population only versions
    # read in the data files for the total wildtype and Wolbachia-infected populations
    df_m = CSV.read(joinpath("Initial_paper_results", "$(JOB_ID)", "$(JOB_ID)", "m_mosqs_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"), DataFrame)
    df_w = CSV.read(joinpath("Initial_paper_results", "$(JOB_ID)", "$(JOB_ID)", "w_mosqs_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"), DataFrame)
    m_array = Array(df_m) # convert to array
    w_array = Array(df_w)
    # read in the data files for the free mosquito populations
    df_free_m = CSV.read(joinpath("Initial_paper_results", "$(JOB_ID)", "$(JOB_ID)", "free_m_mosqs_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"), DataFrame)
    df_free_w = CSV.read(joinpath("Initial_paper_results", "$(JOB_ID)", "$(JOB_ID)", "free_w_mosqs_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"), DataFrame)
    free_m_array = Array(df_free_m) # convert to array
    free_w_array = Array(df_free_w)
    # read in the data files for the tracked households (Wolbachia-infected males and females)
    df_track_male = CSV.read(joinpath("Initial_paper_results", "$(JOB_ID)", "$(JOB_ID)", "track_male_w_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"), DataFrame)
    df_track_fem = CSV.read(joinpath("Initial_paper_results", "$(JOB_ID)", "$(JOB_ID)", "track_fem_w_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"), DataFrame)
    track_male = Array(df_track_male)  # convert to array
    track_fem = Array(df_track_fem)

    return m_array, w_array, free_m_array, free_w_array, track_male, track_fem  # return all data arrays
end
