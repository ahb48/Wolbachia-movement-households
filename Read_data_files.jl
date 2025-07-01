using CSV
using DataFrames
using UnPack

function str_convert(rho, tau, gamma, alpha)
    rho_str = replace(string(rho), "." => "")
    tau_str = replace(string(tau), "." => "")
    gamma_str = replace(string(gamma), "." => "")
    alpha_str = replace(string(alpha), "." => "")

    return rho_str, tau_str, gamma_str, alpha_str
end

function read_data(params,JOB_ID)
    @unpack rho, tau, gamma, alpha, H, rel_t, rel_size, rel_type = params
    ### Convert dispersal parameters to a string and remove the decimal point
    ## We will use this to call the correct data files
    
    rho_str, tau_str, gamma_str, alpha_str = str_convert(rho, tau, gamma, alpha)
    ## Read in the total wildtype and Wolbachia-infected mosquito populations data, 
    ## as well as free population only versions

    df_m = CSV.read(joinpath("sensitivity-anal", "$(JOB_ID)", "$(JOB_ID)", "m_mosqs_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"), DataFrame)
    df_w = CSV.read(joinpath("sensitivity-anal", "$(JOB_ID)", "$(JOB_ID)", "w_mosqs_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"), DataFrame)
    m_array = Array(df_m) # convert to array
    w_array = Array(df_w)
    df_free_m = CSV.read(joinpath("sensitivity-anal", "$(JOB_ID)", "$(JOB_ID)", "free_m_mosqs_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"), DataFrame)
    df_free_w = CSV.read(joinpath("sensitivity-anal", "$(JOB_ID)", "$(JOB_ID)", "free_w_mosqs_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"), DataFrame)
    free_m_array = Array(df_free_m)
    free_w_array = Array(df_free_w)
    df_track_male = CSV.read(joinpath("sensitivity-anal", "$(JOB_ID)", "$(JOB_ID)", "track_male_w_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"), DataFrame)
    df_track_fem = CSV.read(joinpath("sensitivity-anal", "$(JOB_ID)", "$(JOB_ID)", "track_fem_w_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"), DataFrame)
    track_male = Array(df_track_male)
    track_fem = Array(df_track_fem)

    return m_array, w_array, free_m_array, free_w_array, track_male, track_fem
end
