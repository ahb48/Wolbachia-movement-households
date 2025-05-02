############### PACKAGES ##################

using Statistics   # loading Julia packages
using Distributed
using CSV
using DataFrames

############### PARAMETERS ################
JOB_ID = ARGS[1]  # get the job ID from the command line

# Define parameters dictionary
fem_w0 = 0  # how many wildtypes and infected mosquitoes starting off with
male_w0 = 0
fem_m0 = 8
male_m0 = 8
m_free0 = 0 # no. of free wild-type
w_free0 = 0 # no. of free Wolbachia

# model parameters
## vary these for different dispersal scenarios
gamma = parse(Float64, ARGS[2])    # dispersal rate from release site to household
rho = parse(Float64, ARGS[3])      # rate visit the households
tau = parse(Float64, ARGS[4])      # rate leave the households
alpha = parse(Float64, ARGS[5])     # dispersal rate of females moving from household directly to another household

phi = 0.85      # wolbachia fitness effect on birth rate
K = 30          # carrying capacity of mosquitoes per household
d = 12/100      # death rate of mosquitoes
k = 0.3         # larval density parameter
h = 0.19*100^k  # larval density parameter
b = 0.54        # per capita birth rate of mosquitoes (wildtypes)
H = 1000         # number of households
u = 1           # vertical transmission probability
v = 1           # CI effect, proportion of non-viable offspring from infected male and wildtype female birth

rel_t = 150   # initial release time
delta_rel = 0 # time in days between releases, currently set to a single release inside the SSA code
## remember this is per household for the household releases so should use smaller release numbers
## currently set for community wide release
rel_size = 8*H # total number of Wolbachia-infected mosquitoes released into the community

t_start = 0    
t_end = 500   # start time and end time (days) of simulation

num_simulations = 100  # number of realisations of SSA to run, 10,000 is a good number for the paper

parameters = Dict(           # dictionary of parameters
    :fem_m0 => fem_m0,
    :male_m0 => male_m0,
    :fem_w0 => fem_w0,
    :male_w0 => male_w0,
    :m_free0 => m_free0,
    :w_free0 => w_free0,
    :rho => rho,
    :phi => phi,
    :b => b,
    :K => K,
    :d => d,
    :h => h,
    :k => k,
    :u => u,
    :v => v,
    :tau => tau,
    :H => H,
    :t_start => t_start,
    :t_end => t_end,
    :seed => 1234,
    :delta_rel => delta_rel,
    :rel_size => rel_size,
    :gamma => gamma,
    :alpha => alpha,
    :rel_t => rel_t
)

############### SIMULATION ################

# Add worker processes
# remove any workers hanging around e.g. after code crashes
if length(workers()) > 1
    rmprocs(workers())
end

addprocs(6)  # Adjust the number based on your system's cores
@everywhere using SharedArrays   # Load the SharedArrays package on all worker processes
@everywhere using CSV            # Load CSV package on all worker processes
@everywhere using DataFrames     # Load DataFrames package on all worker processes
# Define your module with the function to be run in parallel on all worker processes
@everywhere module Wolbachia_coupled_para ## for household release need to change script
    include("SSA_rel_hhold_fem_switch.jl")  
end

result_length = length(t_start:t_end) # number of time points to store results
weeks = round(Int,result_length/7)    # number of weeks to store results

## creating shared arrays to store results
## stores for all realisations together
m_results = SharedArray{Float64}(result_length,num_simulations)  # stores total number of wildtype mosquitoes each day
w_results = SharedArray{Float64}(result_length,num_simulations)  # stores total number of Wolbachia-infected mosquitoes each day

free_m_results = SharedArray{Float64}(result_length,num_simulations) # stores total number of free wildtype mosquitoes each day
free_w_results = SharedArray{Float64}(result_length,num_simulations) # stores total number of free Wolbachia-infected mosquitoes each day

# stores number of mosquitoes of each type in the tracked household
track_fem_m_results = SharedArray{Int}(result_length,num_simulations) # no. wildtype females
track_male_m_results = SharedArray{Int}(result_length,num_simulations) # no. wildtype males
track_fem_w_results = SharedArray{Int}(result_length,num_simulations)  # no. Wolbachia-infected females
track_male_w_results = SharedArray{Int}(result_length,num_simulations) # no. Wolbachia-infected males

# Convert dispersal parameters to a string and remove the decimal point
rho_str = replace(string(rho), "." => "")
tau_str = replace(string(tau), "." => "")
gamma_str = replace(string(gamma), "." => "")
alpha_str = replace(string(alpha), "." => "")

# Run the simulations in parallel; @time just times the whole block, @sync ensures all 
# processes are done before moving on
@time begin
    @sync @distributed for i in 1:num_simulations
        # Initialize shared arrays for the number of mosquitoes of each type in each household
        # we need to store this data separately for each simulation as too large to store in a single array
        # store per week indtead of per day
        fem_m_hholds = SharedArray{Int}(weeks, H)  # stores number of wildtype females in each household each week
        male_m_hholds = SharedArray{Int}(weeks, H) # stores number of wildtype males in each household each week
        fem_w_hholds = SharedArray{Int}(weeks, H)  # stores number of Wolbachia-infected females in each household each week
        male_w_hholds = SharedArray{Int}(weeks, H) # stores number of Wolbachia-infected males in each household each week

        # Run the simulation
        # collecting the results from the simulation
        m_vec, w_vec, free_m_vec, free_w_vec, num_house_vec, key_house_vec, track_fem_m, track_male_m, 
        track_fem_w, track_male_w = Wolbachia_coupled_para.gillespie(parameters)
        
        # Store results
        m_results[:, i] = m_vec               # wildtypes each day, for sim i
        w_results[:, i] = w_vec               # Wolbachia-infected each day, for sim i
        free_m_results[:, i] = free_m_vec     # free wildtypes each day, for sim i
        free_w_results[:, i] = free_w_vec     # free Wolbachia each day, for sim i
        track_fem_m_results[:, i] = track_fem_m    # tracked household no. female wildtypes each day, for sim i
        track_male_m_results[:, i] = track_male_m  # tracked household no. male wildtypes each day, for sim i
        track_fem_w_results[:, i] = track_fem_w    # tracked household no. female W-infected each day, for sim i
        track_male_w_results[:, i] = track_male_w  # tracked household no. male W-infected each day, for sim i

        ### Process household states data
        # removing the empty entries  from the vectors, i.e. where no mosquitoes left
        num_house_vec = filter(x -> !isempty(x), num_house_vec)  # for data on number of households in each household state each day
        key_house_vec = filter(x -> !isempty(x), key_house_vec)  # for data on household states occupied each day
        vec_length = length(num_house_vec)  # number of days with non-empty data

        for t in 1:vec_length  # loop over each day
            num1 = 1 # set the start index for the household data
            for j in 1:length(key_house_vec[t])  # loop over number of household states occupied at time point
                num = num_house_vec[t][j]  # number of households in this state
                # check household has valid key (so is a household) and there is atleast one household in this state
                if key_house_vec[t][j] != [999, 999, 999, 999] && num > 0
                    num2 = num1 + num - 1  # calculate where the data for this household state
                    fem_m_hholds[t, num1:num2] .= key_house_vec[t][j][1]  # record no. of female wildtypes for all the household in this state
                    male_m_hholds[t, num1:num2] .= key_house_vec[t][j][2] # male wildtypes
                    fem_w_hholds[t, num1:num2] .= key_house_vec[t][j][3]  # female Wolbachia-infected
                    male_w_hholds[t, num1:num2] .= key_house_vec[t][j][4] # male Wolbachia-infected
                    num1 = num2 + 1 # update the start index for the next household state
                end
            end
        end

        # Convert to DataFrames
        df_fem_m_hholds = DataFrame(fem_m_hholds, :auto) # auto sets the column names
        df_male_m_hholds = DataFrame(male_m_hholds, :auto)
        df_fem_w_hholds = DataFrame(fem_w_hholds, :auto)
        df_male_w_hholds = DataFrame(male_w_hholds, :auto)

        # Write to CSV files
        # the strings of values distinguich the different parameter sets
        CSV.write(joinpath(JOB_ID,"fem_m_hholds_$(i)_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_relhhold.csv"), df_fem_m_hholds)
        CSV.write(joinpath(JOB_ID,"male_m_hholds_$(i)_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_relhhold.csv"), df_male_m_hholds)
        CSV.write(joinpath(JOB_ID,"fem_w_hholds_$(i)_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_relhhold.csv"), df_fem_w_hholds)
        CSV.write(joinpath(JOB_ID,"male_w_hholds_$(i)_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_relhhold.csv"), df_male_w_hholds)
    end
end
rmprocs(workers())  # Remove the worker processes

############### SAVE RESULTS ################
## save the rest of the results which are stored simulations all together 
## make the data frames
df_m = DataFrame(m_results, :auto)  # total wildtypes each day each simulation
df_w = DataFrame(w_results, :auto)  # total Wolbachia-infected each day each simulation

df_free_m = DataFrame(free_m_results, :auto)  # total free wildtypes each day each simulation
df_free_w = DataFrame(free_w_results, :auto)  # total free Wolbachia each day each simulation

df_track_fem_m = DataFrame(track_fem_m_results, :auto) # tracked household female wildtypes each day each simulation
df_track_male_m = DataFrame(track_male_m_results, :auto) # tracked household male wildtypes each day each simulation
df_track_fem_w = DataFrame(track_fem_w_results, :auto) # tracked household female Wolbachia each day each simulation
df_track_male_w = DataFrame(track_male_w_results, :auto) # tracked household male Wolbachia each day each simulation

### write the CSV files, again using the parameter values in the file name
# for total wildtypes and Wolbachia-infected mosquitoes
CSV.write(joinpath(JOB_ID,"m_mosqs_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_relhhold.csv"), df_m)
CSV.write(joinpath(JOB_ID,"w_mosqs_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_relhhold.csv"), df_w)
# for total free wildtypes and Wolbachia-infected mosquitoes
CSV.write(joinpath(JOB_ID,"free_m_mosqs_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_relhhold.csv"), df_free_m)
CSV.write(joinpath(JOB_ID,"free_w_mosqs_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_relhhold.csv"), df_free_w)
# for tracked household mosquitoes of each type
CSV.write(joinpath(JOB_ID,"track_fem_m_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_relhhold.csv"), df_track_fem_m)
CSV.write(joinpath(JOB_ID,"track_male_m_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_relhhold.csv"), df_track_male_m)
CSV.write(joinpath(JOB_ID,"track_fem_w_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_relhhold.csv"), df_track_fem_w)
CSV.write(joinpath(JOB_ID,"track_male_w_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_relhhold.csv"), df_track_male_w)


