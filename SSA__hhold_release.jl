#import Pkg
#Pkg.add("StatsBase")
#Pkg.add("Distributions")
#Pkg.add("Random")
#Pkg.add("UnPack")

# load required packages
using Random
using Distributions
using UnPack
using StatsBase

## define variables that don't change here
# households events array
const events_array = ["birth_fem_m", "birth_male_m", "birth_fem_w", "birth_male_w", "death_fem_m", 
                        "death_male_m", "death_fem_w", "death_male_w", "leave_m", "leave_w", "switch_m", "switch_w"]
# free events array
const events_array_free = ["leave_fm", "leave_fw", "death_fm", "death_fw"]
# list of release events i.e. infected male or female from release enters a household
const events_array_release = ["release_male", "release_fem", "death_rel_male", "death_rel_fem"] 
# dictionary of changes to household states for each event that occurs
const event_changes = Dict(         
            "birth_fem_m" => (1, 0, 0, 0),
            "birth_male_m" => (0, 1, 0, 0),
            "birth_fem_w" => (0, 0, 1, 0),
            "birth_male_w" => (0, 0, 0, 1),
            "death_fem_m" => (-1, 0, 0, 0),
            "death_male_m" => (0, -1, 0, 0),
            "death_fem_w" => (0, 0, -1, 0),
            "death_male_w" => (0, 0, 0, -1),
            "leave_m" => (0, -1, 0, 0),
            "leave_w" => (0, 0, 0, -1),
            "switch_m" => (-1, 0, 0, 0),
            "switch_w" => (0, 0, -1, 0)
        )

mutable struct Household_mut
    """Define the Household struct. We allow our structs to be mutable so that we can updates 
    household attributes. A single instance of this struct represennts the households in a certain state."""
    fem_m::Int      # no. of wildtype female mosquitoes in the household
    male_m::Int     # no. of wildtype male mosquitoes in the household
    fem_w::Int      # no. of Wolbachia-infected female mosquitoes in the household
    male_w::Int     # no. of Wolbachia-infected male mosquitoes in the household
    num::Int        # no. of households of this type
    event_weights::Vector{Float64}  # event weights for this household type
    prop_ind::Float64  # propensity of an individual household of this type (sum of event weights)
    prop::Float64  # total propensity of this household type
    event_dist::Categorical{Float64}  # distribution of events for this household type
end


function birth_m(params::Dict, fem_m::Int, male_m::Int, fem_w::Int, male_w::Int, zm::Float64, zw::Float64)
    """Calculates birth rate of wild-type mosquitoes in household. zm and zw are precalculated."""
    @unpack b,K,h,k,u,v,phi = params  # extracting required parameter values from dictionary
    fem_sum = fem_m + fem_w  # precalculate sum of wildtype and Wolbachia-infected female mosquitoes
    male_sum = male_m + male_w # precalculate the male sum
    # need both males and females to reproduce
    if (fem_sum == 0) || (fem_sum >= K) || (male_sum == 0) # if the household is empty or full, return 0
        return 0
    else # otherwise, calculate the birth rate
        return b*zm*exp(-h*(zm + zw)^k)
    end
end

function birth_w(params::Dict, fem_m::Int, male_m::Int, fem_w::Int, male_w::Int, zm::Float64, zw::Float64)
    """Calculates birth rate of Wolbachia-infected mosquitoes in household. zm and zw are precalculated."""
    @unpack b,K,h,k,v,phi = params  # extracting required parameter values from dictionary
    fem_sum = fem_m + fem_w  # precalculate sum of wildtype and Wolbachia-infected mosquitoes
    male_sum = male_m + male_w # precalculate the male sum
    # need both males and females to reproduce
    if (fem_w == 0) || (fem_sum >= K) || (male_sum == 0) # if the household is empty, full, or has no Wolbachia-infected mosquitoes,
        return 0
    else  # otherwise, calculate the birth rate
        return b*zw*exp(-h*(zm + zw)^k)
    end
end

# defines functions for calculating death, leave and switch rates from the household
death(x, d) = d * x
leave(x, tau) = tau * x
switch(x, alpha) = alpha * x

function prop_ind(params::Dict, fem_m::Int, male_m::Int, fem_w::Int, male_w::Int)
    """Function to initialise the total propensity of a household. Used in initialise_gillespie()."""
    @unpack b,K,h,k,d,tau,phi,v,u,alpha = params  # extracting required parameter values from dictionary
    
    m_plus_w = fem_m + male_m + fem_w + male_w  # precalculate sum of wildtype and Wolbachia-infected mosquitoes
    death_fem_m_val = death(fem_m, d) # calculate death rate of wildtype female mosquitoes in household
    death_male_m_val = death(male_m, d) # calculate death rate of wildtype male mpsquitoes in household
    death_fem_w_val = death(fem_w, d) # calculate death rate of Wolbachia-infected female mosquitoes in household
    death_male_w_val = death(male_w, d) # calculate death rate of Wolbachia-infected male mosquitoes in household
    leave_m_val = leave(male_m, tau) # calculate leave rate of wildtype male mosquitoes in household
    leave_w_val = leave(male_w, tau) # calculate leave rate of Wolbachia-infected male osquitoes in household
    switch_m_val = switch(fem_m, alpha) # calculate the rate at which wildtype female mosquitoes switch households
    switch_w_val = switch(fem_w, alpha) # calculate the rate at which Wolbachia-infected females switch households

    if (m_plus_w == 0) # if the household is empty,
        # return event weights
        return [0, 0, 0, 0, death_fem_m_val, death_male_m_val, death_fem_w_val, 
                death_male_w_val, leave_m_val, leave_w_val, switch_m_val, switch_w_val] # birth rates are 0
    else
        # else find and include the birth rates
        male_sum = male_m + male_w  # precalculate male sum
        term1 = (1-v)*phi
        # precalculate zm and zw terms
        zm = (male_m*(fem_m + fem_w*term1) + male_w*(fem_m*(1-u) + fem_w*term1))/male_sum
        zw = v*phi*fem_w
        # calculate birth rates
        birth_m_val = birth_m(params, fem_m, male_m, fem_w, male_w, zm, zw)
        birth_w_val = birth_w(params, fem_m, male_m, fem_w, male_w, zm, zw)
        # return event weights
        return [birth_m_val, birth_m_val, birth_w_val, birth_w_val, 
                death_fem_m_val, death_male_m_val, death_fem_w_val, death_male_w_val,
                leave_m_val, leave_w_val, switch_m_val, switch_w_val]
    end
end

function update_number!(increment::Int, household::Household_mut, params::Dict)
    """Function that updates the total number of households in the given household state and total propensity.
    Call this function whenever the number of households with that state changes."""
    household.num += increment  # updates attribute for number of households in this state
    # updates attribute for total propensity of this household type
    household.prop = household.prop_ind * household.num  
end

function get_household_event(household::Household_mut)
    """Function that determines which event occurs. The SSA returns that a broad event has occurred 
    for a household of a given state, this function splits that broad event into a more specific event."""
    # The list of events and their weights are a fixed attribute of the object, just need to select one at random
    # Normalize the weights to sum to 1
    event_index = rand(household.event_dist)
    event = events_array[event_index]
    # Return the event to be processed
    return event
end

###################################

mutable struct Track_mut
    """Define the Track struct. We allow our structs to be mutable so that we can update 
    track attributes. We use a single instance of this struct to track a single household. This records 
    the numbers of each type of mosquito in the tracked household over time."""
    fem_m::Int   # no. of wildtype females in the tracked household
    male_m::Int  # no. of wildtype males in the tracked household
    fem_w::Int   # no. of Wolbachia-infected females in the tracked household
    male_w::Int  # no. of Wolbachia-infected males in the tracked household
    num::Int    # no. of households of this type, always only one tracked household
    event_weights::Vector{Float64}  # event weights for this household type
    prop_ind::Float64  # propensity of an individual household of this type (sum of event weights)
    prop::Float64  # total propensity of this household type
    event_dist::Categorical{Float64}  # distribution of events for this household type
end

### can use most of the functions defined for household for track as well

function get_household_event(track::Track_mut)
    """Function that determines which event occurs. The SSA returns that a broad event has occurred 
    for a household of a given configuration, this function splits that broad event into a more specific event.
    This is the same as the function defined for the household struct, but takes the an instance track 
    struct as input."""
    # The list of events and their weights are a fixed attribute of the object, just need to select one at random
    # Normalize the weights to sum to 1
    event_index = rand(track.event_dist)
    event = events_array[event_index]
    # Return the event to be processed
    return event
end

###################################

mutable struct Release_mut
    """Define the Release struct. We allow our structs to be mutable so that we can update 
    release attributes. We use a single instance of this struct to record the number of released mosquitoes
    of each type."""
    fem_m::Int   # released female wildtype mosquitoes, set to zero
    male_m::Int  # released male wildtype mosquitoes, set to zero
    fem_w::Int   # free female Wolbachia-infected mosquitoes
    male_w::Int  # free male Wolbachia-infected mosquitoes
    num::Int  # number of households in this state - always equal to 0 since not a household state
    prop::Float64  # sum of all event rates for the release mosquitoes
    prop_ind::Vector{Float64}  # propensities of the individual release mosquito types
    event_dist::Categorical{Float64}  # distribution of release events
end

function update_rates!(params_dict::Dict, release::Release_mut)
    """Updates the total propensity of the free mosquitoes."""
    @unpack gamma, d = params_dict # unpack the parameter values from the dictionary
    # array of propensities for the individual free mosquito types
    release.prop_ind = [gamma * release.male_w, gamma * release.fem_w, d * release.male_w, d * release.fem_w] 
    release.prop = sum(release.prop_ind) # total propensity of the released mosquitoes
    # set the event distribution
    # we use a ternary operator to set the event distribution to some accepted placeholder values if the total propensity is 0
    release.event_dist = release.prop == 0 ? Categorical([0.25, 0.25, 0.25, 0.25]) : Categorical(release.prop_ind / release.prop) 
end

function update_number!(increment::Int, mtype::AbstractString, release::Release_mut, params_dict::Dict)
    """Updates the number of released mosquitoes recorded in states_dict."""
    if mtype == "male_w"    # if it's the wildtype mosquitoes that are changing
        release.male_w += increment   # add required increment to free wildtypes male mosquitoes
    else
        release.fem_w += increment   # else add increment to free Wolbachia-infected male mosquitoes
    end

    # update free mosquitoes total propensity
    update_rates!(params_dict,release)
end

function get_release_event(release::Release_mut)
    """Chooses an event at random with weights corresponding to propensities 
    (which depend on population size)."""
    event_indx = rand(release.event_dist)  # sampling event index
    event = events_array_release[event_indx]  # get the event corresponding to the index
    return event  # return event
end

###################################

mutable struct Free_mut
    """Define the free struct. We allow our structs to be mutable so that we can update 
    free attributes. We use a single instance of this struct to record the number of free mosquitoes."""
    fem_m::Int   # free female wildtype mosquitoes, set to zero
    male_m::Int  # free male wildtype mosquitoes
    fem_w::Int   # free female Wolbachia-infected mosquitoes, set to zero
    male_w::Int  # free male Wolbachia-infected mosquitoes
    num::Int   # number of households in this state - always equal to 0 since not a household state
    prop::Float64  # sum of all event rates for the free mosquitoes
    prop_ind::Vector{Float64}  # propensities of the individual free mosquitoe types
    event_dist::Categorical{Float64} # distribution of free events
end

function update_rates!(params_dict::Dict, free::Free_mut)
    """Updates the total propensity of the free mosquitoes."""
    @unpack rho, d = params_dict # unpack the parameter values from the dictionary
    free.prop_ind = [rho * free.male_m, rho * free.male_w, d * free.male_m, d * free.male_w] # array of propensities for the individual free mosquito types
    free.prop = sum(free.prop_ind) # total propensity of the free mosquitoes
    # set the event distribution
    # we use a ternary operator to set the event distribution to some accepted placeholder values if the total propensity is 0
    free.event_dist = free.prop == 0 ? Categorical([0.25, 0.25, 0.25, 0.25]) : Categorical(free.prop_ind / free.prop) # set the event distribution
end

function update_number!(increment::Int, mtype::AbstractString, free::Free_mut, params_dict::Dict)
    """Updates the number of free mosquitoes recorded in states_dict."""
    if mtype == "male_m"    # if it's the wildtype mosquitoes that are changing
        free.male_m += increment   # add required increment to free wildtypes male mosquitoes
    else
        free.male_w += increment   # else add increment to free Wolbachia-infected male mosquitoes
    end

    # update free mosquitoes total propensity
    update_rates!(params_dict,free)
end

function get_free_event(free::Free_mut)
    """Chooses an event at random with weights corresponding to propensities 
    (which depend on population size)."""
    event_indx = rand(free.event_dist)  # sampling event index
    event = events_array_free[event_indx]  # get the event corresponding to the index
    return event  # return event
end

###################################
# saves a dictionary to contain the household keys already collected       
const string_dict = Dict{Tuple{Int64, Int64, Int64, Int64}, String}()  

function make_key_house(fem_m::Int, male_m::Int, fem_w::Int, male_w::Int)
    """Function to create dictionary keys that summarise household attributes. If key already exists,
    takes it from the dictionary, otherwise creates a new key and adds it to the dictionary."""
    get!(string_dict, (fem_m, male_m, fem_w, male_w)) do
        string(fem_m, "-", male_m, "-", fem_w, "-", male_w)
    end
end

function do_event(states_dict::Dict, event_key::AbstractString, event::AbstractString, params_dict::Dict)
    """Gillespie produces an event_key - the state in which an event occurs, and an event 
    - the name of the event that occurs. This function determines which other states are affected by the 
    event, and then updates the state in which the event occurs, and any other states affected."""
    @unpack K = params_dict  # unpack the parameter values from the dictionary

    no_key = "999-999-999-999"  # placeholder key
    # initialise key for household/space that event occurs in and key of new household state of initially affected household (if there is one)
    house_key_out, house_key_in = event_key, no_key 

    function sample_household(states_dict, K)
        """Samples a household for the free/ released mosquito to enter. Uses the number of households
        in each state as weights for the sampling. Non household states are automatically excluded since num 
        is always 0 for these states."""
        key_list = collect(keys(states_dict))  # get the keys of the states dictionary
        # using ternary operator to set the weights to 0 if the household is full
        weights_list = [(value.fem_m + value.male_m + value.fem_w + value.male_w < 2*K) ? value.num : 0 for value in values(states_dict)]
        weights_list /= sum(weights_list) # normalise the weights to sum to 1
        dist = Categorical(weights_list) # create a categorical distribution from the weights
        return key_list[rand(dist)] # return the key of the sampled household
    end

    function switch_household(states_dict, K, house_key_out)
        """Samples a new household for the female mosquito to switch to. This function behaves the same way 
        as sample_household, but excludes the household that the mosquito is leaving from the sampling."""
        key_list = collect(keys(states_dict)) # get the keys of the states dictionary
        weights_list = zeros(Float64, length(key_list)) # initialise the weights list
        
        total_weight = 0.0 # initialise the total weight
        for (i, key) in enumerate(key_list)  # loop through the keys and keep index i
            value = states_dict[key]  # get the household state
            # set the weights to 0 if household is full
            if value.fem_m + value.male_m + value.fem_w + value.male_w >= 2*K
                weights_list[i] = 0.0
            # take one household away from the weight of the household the mosquito is leaving
            elseif key == house_key_out
                weights_list[i] = value.num - 1.0
            else  # else can just take the weight of the household
                weights_list[i] = value.num
            end
            total_weight += weights_list[i] # add the weight to the total
        end
        
        if total_weight == 0 # if the sum of the weights is 0, throw an error
            throw(DomainError(weights_list, "Sum of weights is zero"))
        end
        
        weights_list /= total_weight  # normalise the weights to sum to 1
        dist = Categorical(weights_list)  # create a categorical distribution from the weights
        return key_list[rand(dist)]  # return the key of the sampled household
    end

    function extract_household_state(states_dict, house_key)
        """Function to extract the household state correspinding to a given key."""
        state = states_dict[house_key]  # get the household state from the dictionary
        # return the number of mosquitoes of each type in the household
        return state.fem_m, state.male_m, state.fem_w, state.male_w 
    end

    function update_track(states_dict, fem_m, male_m, fem_w, male_w, params_dict)
        """Updates the tracked household instance."""
        if fem_m + male_m + fem_w + male_w > 0   # if the household is not empty, update the tracked household instance
            event_ws = prop_ind(params_dict, fem_m, male_m, fem_w, male_w)  # get the event weights for this household
            prop_i = sum(event_ws)  # total propensity for the individual household
            dist = event_ws / prop_i  # normalize the event weights to get the event distribution
            # update the tracked household instance in the dictionary
            states_dict["track"] = Track_mut(fem_m, male_m, fem_w, male_w, 1, dist, prop_i, prop_i * 1, Categorical(dist))
        else  # if the household is empty, update the tracked household instance with all events set to 0,
            # but have to have a distribution that sums to 1 for the Categorical function or it throws an error
            states_dict["track"] = Track_mut(fem_m, male_m, fem_w, male_w, 1, zeros(Int, 12), 0, 0, Categorical([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
        end
    end

    function handle_event(event_key, event, states_dict, params_dict)
        """Function to handle events that occur within a household. Updates the numbers of each mosquito
        type in the current household, depending on the event that occurs."""
        fem_m, male_m, fem_w, male_w = extract_household_state(states_dict, event_key) # get the household state
        d_fem_m, d_male_m, d_fem_w, d_male_w = event_changes[event] # get the changes to the household state
        fem_m += d_fem_m # update the numbers of each mosquito type in the household
        male_m += d_male_m
        fem_w += d_fem_w
        male_w += d_male_w
        return fem_m, male_m, fem_w, male_w # return the updated numbers of each mosquito type
    end

    function process_key_in(house_key_out, house_key_in, states_dict, params_dict, fem_m, male_m, fem_w, male_w)
        """Function to make the key of the new household state after a household event occurs. If the event
        occurs in the tracked household it updates the tracked household instance instead."""
        house_key_in = no_key # initialise the key of the new household state
        if house_key_out == "track" # if the event occurred in the tracked household
            update_track(states_dict, fem_m, male_m, fem_w, male_w, params_dict) # update the tracked household instance
        else # otherwise make the new house key
            house_key_in = make_key_house(fem_m, male_m, fem_w, male_w)
        end
        
        return house_key_in  # return the key of the new household state (will return the placeholder for no key if the event occurred in the tracked household)
    end

    if event_key == "release" ## if event occurred in the release state
        if event == "death_rel_male"  # if a released male dies
            update_number!(-1, "male_w", states_dict["release"], params_dict) # update number of males in release pool
            house_key_out = no_key # no household to update
        elseif event == "death_rel_fem" # if a released female dies
            update_number!(-1, "fem_w", states_dict["release"], params_dict) # update number of females in release pool
            house_key_out = no_key # no household to update
        else # else a released mosquito moves into a household
            house_key_out = sample_household(states_dict, K) # sample household for mosquito to enter
            # get the household state of the household the mosquito moves into
            fem_m, male_m, fem_w, male_w = extract_household_state(states_dict, house_key_out)
            if event == "release_male" #  if released male enters household
                male_w += 1 # add one to W-infected males in household
                update_number!(-1, "male_w", states_dict["release"], params_dict) # update number of males in released pool
            else  # else released female enters household
                fem_w += 1  # add one to W-infected females in the household
                update_number!(-1, "fem_w", states_dict["release"], params_dict)  # update number of females in released pool
            end
            # creates the key of the updated household (or updates the tracked household)
            house_key_in = process_key_in(house_key_out, house_key_in, states_dict, params_dict, fem_m, male_m, fem_w, male_w)
        end

    elseif event_key == "free" ## if event occurred in the free state
        if event == "death_fm"  # if a free wildtype dies (male)
            # update the number of free wildtypes in the free pool
            update_number!(-1, "male_m", states_dict["free"], params_dict)
            house_key_out = no_key  # no household to update
        elseif event == "death_fw" # if a free Wolbachia-infected dies
            # update the number of free Wolbachia-infected in the free pool
            update_number!(-1, "male_w", states_dict["free"], params_dict)
            house_key_out = no_key # no household to update
        else # free mosquito moves into a household
            house_key_out = sample_household(states_dict, K)  # sample household for mosquito to enter
            # get the household state of the household the mosquito moves into
            fem_m, male_m, fem_w, male_w = extract_household_state(states_dict, house_key_out)
            if event == "leave_fm" # if a free wildtype enters a household
                # update the number of free wildtypes in the free pool
                update_number!(-1, "male_m", states_dict["free"], params_dict)
                male_m += 1 # add one to the wildtype males in the household entered
            else # if a free Wolbachia-infected enters a household
                # update the number of free Wolbachia-infected in the free pool
                update_number!(-1, "male_w", states_dict["free"], params_dict)
                male_w += 1 # add one to the Wolbachia-infected males in the household entered
            end
            # creates the key of the updated household (or updates the tracked household)
            house_key_in = process_key_in(house_key_out, house_key_in, states_dict, params_dict, fem_m, male_m, fem_w, male_w)
        end
    else # event occurred within the household
        # get the new household state of household event occurred in
        fem_m, male_m, fem_w, male_w = handle_event(event_key, event, states_dict, params_dict)
        # creates the key of the updated household (or updates the tracked household)
        house_key_in = process_key_in(event_key, house_key_in, states_dict, params_dict, fem_m, male_m, fem_w, male_w)
        
        if event == "leave_m" # if a wildtype (male) leave event
            update_number!(1, "male_m", states_dict["free"], params_dict) # update the number of free wildtype males
        elseif event == "leave_w" # if a Wolbachia-infected (male) leave event
            update_number!(1, "male_w", states_dict["free"], params_dict) # update the number of free Wolbachia males
        
        elseif event in ["switch_m", "switch_w"] # otherwise female switching households event
            house_key_switch_to = switch_household(states_dict, K, house_key_out) # sample household switches to
            # get the household state of the household the female switches to
            fem_m_, male_m_, fem_w_, male_w_ = extract_household_state(states_dict, house_key_switch_to)
            # updates the number of households in the household state just switched to (household state of that household changes)
            if house_key_switch_to != "track" # assuming household switch to is not the tracked household
                update_household_number!(states_dict, house_key_switch_to, -1, params_dict) 
            end
            if event == "switch_m" # if event is a wildtype-female switching households
                fem_m_ += 1 # add one to the number of wildtype females in the household
            elseif event == "switch_w" # if event is a Wolbachia-infected female switching households
                fem_w_ += 1 # add one to the number of Wolbachia-infected females in the household
            end
        
            # make a house key for the new household state, after switching into household
            house_key_switch_in = make_key_house(fem_m_, male_m_, fem_w_, male_w_)
            if house_key_switch_to == "track" # if the household the female switches to is the tracked household
                update_track(states_dict, fem_m_, male_m_, fem_w_, male_w_, params_dict) # update the tracked household instance
            else # if the household the female switches to is any other household
                # update the number of households in the new state of the household switched to
                update_household_number!(states_dict, house_key_switch_in, 1, params_dict)
            end
        end

    end
    ## updates the number of households in the household state the event initially occured in
    update_household_number!(states_dict, house_key_out, -1, params_dict)
    ## updates the number of households in the household state the event changes the initial household to
    update_household_number!(states_dict, house_key_in, 1, params_dict)
end

function update_household_number!(states_dict::Dict, house_key::AbstractString, inc::Int, params_dict::Dict)
    """Updates the number of households in state house_key by amount inc.
    If the state does not exist, create it; if the number in state drops to 0, remove it. But if key is
    not a valid household key, do nothing."""

    no_key = "999-999-999-999"  # placeholder key
    # If house_key is the placeholder value or track,
    # no household to update, event must have affected free/released population or the tracked household only
    if (house_key === no_key) || (house_key == "track")
        return
    end
    # If this household type already exists, its key will be in the dictionary
    # So call its update_number function, which will also recalculate its propensity
    household = get(states_dict, house_key, nothing)
    if household !== nothing
        update_number!(inc, household, params_dict)
        # If the number of households in the state is 0, remove it
        if household.num == 0
            delete!(states_dict, house_key)
        end
    else   # if the household type does not exist, create it
        fem_m, male_m, fem_w, male_w = get_value(house_key)  # get the number of mosquitoes of each type in the household
        if fem_m + male_m + fem_w + male_w > 0   # if the household is not empty, add it to the dictionary
            event_ws = prop_ind(params_dict, fem_m, male_m, fem_w, male_w)  # get the event weights for this household
            prop_i = sum(event_ws)  # total propensity for the individual household
            dist = event_ws/prop_i  # normalise the event weights to get the event distribution
            # add the household to the dictionary, have to create a new instance of the Household_mut struct
            states_dict[house_key] = Household_mut(fem_m, male_m, fem_w, male_w, inc, dist, prop_i, prop_i*inc, Categorical(dist))
        else 
            # if the household is empty, all events are 0
            # but have to have a distribution that sums to 1 for the Categorical function or it throws an error
            # this is just a placeholder, as prop is 0 so the household will never be chosen,
            # but may be entered by free mosquitoes
            states_dict[house_key] = Household_mut(fem_m, male_m, fem_w, male_w, inc, zeros(Int, 12), 0, 0, Categorical([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]))
        end
    end
end

function get_value(house_key::String)
    """Converts string from household key to integer array."""
    # get the number of mosquitoes of type component in a household from its key
    # key is of the form m-w so need to cut string at the -
    return map(x -> parse(Int, x), split(house_key, '-'))  # gives a list of integers
end

function get_population_size(states_dict::Dict)
    """Sums up the total population size from the state dictionary of the wildtypes, Wolbachia-infected
    free wildtypes, free W-infected and male and female wildtypes and W-infected seperately."""
    # initialise the total population sizes as empty
    m_pop, w_pop, free_m_pop, free_w_pop, male_m_pop, fem_m_pop, male_w_pop, fem_w_pop = 0, 0, 0, 0, 0, 0, 0, 0 
    release_fem_pop, release_male_pop = 0, 0 # initialise the total population sizes

    if haskey(states_dict, "free")       # if the free state exists, add its population to the corresponding totals
        free_state = states_dict["free"]
        m_pop += free_state.male_m
        w_pop += free_state.male_w
        free_m_pop += free_state.male_m
        free_w_pop += free_state.male_w
        male_m_pop += free_state.male_m
        male_w_pop += free_state.male_w
    end

    if haskey(states_dict, "release")  # if the release state exists, add its population to the corresponding totals
        release_state = states_dict["release"]
        w_pop += release_state.male_w + release_state.fem_w
        fem_w_pop += release_state.fem_w
        male_w_pop += release_state.male_w
        release_fem_pop += release_state.fem_w
        release_male_pop += release_state.male_w
    end

    if haskey(states_dict, "track")  # if the track state exists, add its population to the corresponding totals
        track_state = states_dict["track"]
        m_pop += track_state.fem_m + track_state.male_m
        w_pop += track_state.fem_w + track_state.male_w
        male_m_pop += track_state.male_m
        fem_m_pop += track_state.fem_m
        male_w_pop += track_state.male_w
        fem_w_pop += track_state.fem_w
    end

    for (key, value) in states_dict    # looping through the states dictionary
        if (key != "free") && (key != "release") && (key != "track")   # if a household state
            m_pop += (value.fem_m + value.male_m) * value.num  # need to multiply by the number of households in this state
            w_pop += (value.fem_w + value.male_w) * value.num
            male_m_pop += value.male_m * value.num
            fem_m_pop += value.fem_m * value.num
            male_w_pop += value.male_w * value.num
            fem_w_pop += value.fem_w * value.num
        end
    end
    # returns the total numbers of mosquito types
    return m_pop, w_pop, free_m_pop, free_w_pop, male_m_pop, fem_m_pop, male_w_pop, fem_w_pop, release_fem_pop, release_male_pop 
end

function do_release(states_dict::Dict, params_dict::Dict)
    """Function to release a number of Wolbachia-infected mosquitoes into each household."""
    @unpack K, rel_male, rel_fem = params_dict  # unpack the parameter values from the dictionary
    rel_size = rel_fem + rel_male  # total number of released mosquitoes
    keys_list = collect(keys(states_dict))  # get the keys of the states dictionary
    # create a new dictionary to store the updated household states
    states_dict_ = Dict{String, Union{Household_mut, Free_mut, Release_mut, Track_mut}}()
    states_dict_["free"] = states_dict["free"]        # copy the free state from the original dictionary
    states_dict_["release"] = states_dict["release"]  # copy the release state from the original dictionary
    
    for key in keys_list    # loop through the keys of the states dictionary
        if key in ["free", "release"]  # if the key is free or release, skip it
            continue
        end
        
        household = states_dict[key]  # get the household state from the dictionary
        H_ = household.num   # number of households in this state
        # get the number of mosquitoes of each type in the household
        fem_m, male_m, fem_w, male_w = household.fem_m, household.male_m, household.fem_w, household.male_w
        m_plus_w = fem_m + male_m + fem_w + male_w   # total number of mosquitoes in the household
        
        if m_plus_w + rel_size > 2 * K   # if the household is full, we need to sample mosquitoes
            # Calculate the number of mosquitoes to sample
            num_to_sample = round(Int, max(2 * K - m_plus_w, 0))
            
            # Combine the weights of female and male mosquitoes
            weights = vcat(ones(rel_fem), ones(rel_male))
            
            # Create an array of indices representing the mosquitoes
            mosquitoes = 1:rel_size
            
            # Sample the mosquitoes
            sampled_mosquitoes = sample(mosquitoes, Weights(weights), num_to_sample, replace=false)
            
            # Determine the number of female and male mosquitoes sampled
            num_females_sampled = count(x -> x <= rel_fem, sampled_mosquitoes)
            num_males_sampled = count(x -> x > rel_fem, sampled_mosquitoes)
            
            fem_w += num_females_sampled   # add sampled females
            male_w += num_males_sampled    # add sampled males
        else    # if household not full
            fem_w += rel_fem    # add females
            male_w += rel_male   # add females
        end
        # create household key
        key_ = key == "track" ? "track" : make_key_house(fem_m, male_m, fem_w, male_w)
        event_ws = prop_ind(params_dict, fem_m, male_m, fem_w, male_w)  # event weights for the initial household state
        prop_i = sum(event_ws)   # propensity of an individual household of this type
        dist = event_ws / prop_i   # normalize the event weights to get the event distribution
        states_dict_[key_] = Household_mut(fem_m, male_m, fem_w, male_w, H_, dist, prop_i, prop_i * H_, Categorical(dist)) 
    end
    return states_dict_   # return new state dictionary
end



function initialise_gillespie(params_dict::Dict)
    """Initialise dictionary for initial household configurations, free mosquitoes, 
    released mosquitoes and tracked household. Initialises arrays to record various information 
    and set time steps."""
    # Initially, all households are in state (m0, w0); make a key corresponding to this household state
    @unpack fem_m0, male_m0, fem_w0, male_w0, H, m_free0, w_free0, rho, t_start, t_end, rel_t = params_dict  # unpack the parameter values
    
    # Define a union type, so that dictionary can take entries of either struct
    UnionType = Union{Household_mut, Free_mut, Release_mut, Track_mut}

    # Initialise states dictionary as empty
    states_dict = Dict{String, UnionType}() 
    # add all households (total number H-1, because we track one specific household) in the initial state
    event_ws = prop_ind(params_dict, fem_m0, male_m0, fem_w0, male_w0)  # event weights for the initial household state
    prop_i = sum(event_ws)   # propensity of an individual household of this type
    dist = event_ws / prop_i   # normalise the event weights to get the event distribution
    # add the household to the dictionary, have to create a new instance of the Household_mut struct
    states_dict[make_key_house(fem_m0, male_m0, fem_w0, male_w0)] = Household_mut(fem_m0, male_m0, fem_w0, male_w0, H-1, dist, prop_i, prop_i*(H-1), Categorical(dist))
    # define separate dictionary entry for the household we are tracking, using an instance of the Track_mut struct
    states_dict["track"] = Track_mut(fem_m0, male_m0, fem_w0, male_w0, 1, dist, prop_i, prop_i*1, Categorical(dist))

    # Add the free state to the states dictionary, initially m_free0, w_free0 in the free state
    if m_free0 + w_free0 == 0  # if the free state is empty, set prop to 0
        # here the event_dist entry is just a place holder and its values shouldn't matter
        # as prop is zero
        states_dict["free"] = Free_mut(0, m_free0*H, 0, w_free0*H, 0, 0, [0.0,0.0,0.0,0.0], Categorical([0.25, 0.25, 0.25, 0.25]))
    else # if the free state is not empty, calculate the propensities
        propf_mw = [rho * m_free0 * H, rho * w_free0 * H, d * m_free0 * H, d * w_free0 * H]
        propf = sum(propf_mw)
        states_dict["free"] = Free_mut(0, m_free0 * H, 0, w_free0 * H, 0, propf, propf_mw, Categorical(propf_mw / propf))
    end

    # Add the release state to the states dictionary, initially 0 in the release state
    states_dict["release"] = Release_mut(0, 0, 0, 0, 0, 0, [0.0,0.0,0.0,0.0], Categorical([0.25, 0.25, 0.25, 0.25]))

    ## Set up the time points
    t = t_start                # Set current time to start time
    t_vec = t_start:t_end      # Time points for Gillespie output, daily
    len_t_vec = length(t_vec)  # Precalculate length of t_vec

    i_out = 2                  # Set time point counter to first point after start point
    t_out = t_vec[i_out]       # Set time for next output

    # To store no. of wild-type mosquitoes at each time point, can fill in first time point
    m_vec = zeros(Float64, len_t_vec)
    m_vec[1] = (fem_m0 + male_m0) * H + m_free0 * H
    # To store no. of Wolbachia-infected mosquitoes at each time point, can fill in first time point
    w_vec = zeros(Float64, len_t_vec)
    w_vec[1] = (fem_w0 + male_w0) * H + w_free0 * H
    # Stores number of free wiltype and Wolbachia-infected mosquitoes at each time point
    # Have initialised by filling in all time points with the initial value, this is more efficient
    free_m_vec = fill(m_free0 * H, len_t_vec)
    free_w_vec = fill(w_free0 * H, len_t_vec)
    # initialise male wildtype vector
    male_m_vec = zeros(Float64, len_t_vec)
    male_m_vec[1] = male_m0 * H
    # initialise female wildtype vector
    fem_m_vec = zeros(Float64, len_t_vec)
    fem_m_vec[1] = fem_m0 * H
    # initialise female Wolbachia-infected vector
    male_w_vec = zeros(Float64, len_t_vec)
    male_w_vec[1] = male_w0 * H
    # initialise male Wolbachia-infected vector
    fem_w_vec = zeros(Float64, len_t_vec)
    fem_w_vec[1] = fem_w0 * H

    # To store time exact time points used
    tOut = zeros(Float64, len_t_vec)

    # To store number of households of each type at each time point
    # only store every week
    num_house_vec = [Vector{Int64}() for _ in 1:round(len_t_vec/7)]
    # To store list of household keys once a week
    key_house_vec = [Vector{Any}() for _ in 1:round(len_t_vec/7)]  
    # can fill in first time point
    # returns a list of the number of households of each type at the first time point
    num_house_vec[1] = [value.num for value in values(states_dict)]
    # returns a list of the household keys at the first time point, uses ternary operator to set key to [999,999,999,999] if the key is free
    key_house_vec[1] = [(key == "free" || key == "release") ? [999,999,999,999] : [value.fem_m, value.male_m, value.fem_w, value.male_w] for (key, value) in states_dict]
    
    # To store the household state of the tracked household at each time point
    track_fem_m = zeros(Int, len_t_vec)  # initialise array for female wildtypes
    track_fem_m[1] = fem_m0  # fill in the first time point
    track_male_m = zeros(Int, len_t_vec) # initialise array for male wildtypes
    track_male_m[1] = male_m0 # fill in the first time point
    track_fem_w = zeros(Int, len_t_vec)  # initialise array for female Wolbachia-infected
    track_fem_w[1] = fem_w0 # fill in the first time point
    track_male_w = zeros(Int, len_t_vec) # initialise array for male Wolbachia-infected
    track_male_w[1] = male_w0 # fill in the first time point

    # To store the number of released
    release_fem_vec = zeros(Int, len_t_vec)  # initialise array for number of released female
    release_fem_vec[1] = 0  # fill in the first time point
    release_male_vec = zeros(Int, len_t_vec)  # initialise array for number of released males
    release_male_vec[1] = 0  # fill in the first time point

    # return all the initialised values/ arrays
    return states_dict, t, i_out, t_out, t_vec, m_vec, w_vec, free_m_vec, free_w_vec, male_m_vec, fem_m_vec, male_w_vec, fem_w_vec, tOut, rel_t, 
    num_house_vec, key_house_vec, track_fem_m, track_male_m, track_fem_w, track_male_w, release_fem_vec, release_male_vec 
end

function gillespie(params_dict::Dict)
    """Function to carry out the Gillespie SSA and return results. 
    Many of the functions we have created above feed into this."""
    # Initialise states dictionary, time counters for simulation, arrays to store results and some parameter values etc...
    states_dict, t, i_out, t_out, t_vec, m_vec, w_vec, free_m_vec, free_w_vec, male_m_vec, fem_m_vec, male_w_vec, fem_w_vec, tOut, rel_t, 
    num_house_vec, key_house_vec, track_fem_m, track_male_m, track_fem_w, track_male_w, release_fem_vec, release_male_vec = initialise_gillespie(params_dict)
    rel_time = true  # boolean to check if release event has occurred, we only want one release event to occur
    t_week = 7.0  # week interval to record number of households
    w_out = 2 # index week counter to record number of households
    weeks = round(length(t_vec)/7)

    while true  # Keeps iterating until break when passes final output time, or population goes extinct

        # Extract the keys and propensities from the states dictionary
        keys_list = collect(keys(states_dict))
        # Julia should be able to recognise which function/struct referring to depending on whether 
        # value is household or free etc
        # vector of total propensities for each household state
        prop_list = [value.prop for value in values(states_dict)]
        # sum of all propensities
        prop_sum = sum(prop_list)
        # Find the time until the next reaction
        delta_t = -log(rand()) / prop_sum

        # Find the next reaction
        # Determine the state involved in the next reaction by
        # selecting state key from list randomly with weights given by the event propensities
        weights = prop_list / prop_sum     # normalise the propensities
        dist = Categorical(weights)        # create a discrete distribution with the weights
        event_key = keys_list[rand(dist)]  # sample the key of the next state to be involved in an event

        # Event_key is the household state, or free/release/tracked state, where an event has occurred
        # Process that event
        if event_key == "free"  # If event occurs in free population
            do_event(states_dict, event_key, get_free_event(states_dict[event_key]), params_dict)
        elseif event_key == "release"  # If event occurs in released population
            do_event(states_dict, event_key, get_release_event(states_dict[event_key]), params_dict)
        elseif event_key == "track"  # If event occurs in tracked household
            do_event(states_dict, event_key, get_household_event(states_dict[event_key]), params_dict)
        else      # If event occurs in household population
            do_event(states_dict, event_key, get_household_event(states_dict[event_key]), params_dict)
        end

        # Update time
        t += delta_t

        # check if time for release
        if t >= rel_t && rel_time # update the the household states after a release 
            states_dict = do_release(states_dict, params_dict)
            rel_time = false # this ensures only a single release event occurs
        end

        # record the total population sizes of all the types at the current time point
        m_pop, w_pop, free_m_pop, free_w_pop, male_m_pop, fem_m_pop, male_w_pop, fem_w_pop, release_fem_pop, release_male_pop = get_population_size(states_dict) 
        
        if t >= t_out                       # Check if past next time point            
            tOut[i_out] = t                 # If so, record time event occurred (should be close to time point)
            m_vec[i_out] = m_pop            # Record no. wild-type mosquitoes at time step
            w_vec[i_out] = w_pop            # Record no. Wolbachia-infected mosquitoes at time step
            free_m_vec[i_out] = free_m_pop  # Record no. free wild-type mosquitoes at time step
            free_w_vec[i_out] = free_w_pop  # Record no. free Wolbachia-infected mosquitoes at time step
            male_m_vec[i_out] = male_m_pop  # male wildtypes
            fem_m_vec[i_out] = fem_m_pop    # female wildtypes
            male_w_vec[i_out] = male_w_pop  # male Wolbachia-infected
            fem_w_vec[i_out] = fem_w_pop    # female Wolbachia-infected

            # record release population sizes
            release_fem_vec = release_fem_pop
            release_male_vec = release_male_pop

            # record the state of the tracked household
            track_fem_m[i_out] = states_dict["track"].fem_m    # female wildtypes
            track_male_m[i_out] = states_dict["track"].male_m  # male wildtypes
            track_fem_w[i_out] = states_dict["track"].fem_w    # female Wolbachia-infected
            track_male_w[i_out] = states_dict["track"].male_w  # male Wolbachia-infected

            ## want to record number of households of each type at each time point
            # record number of households of each type present
            if t >= t_week && w_out <= weeks  # only record every week
                # get the number of households in each state present at the current time point
                num_house_vec[w_out] = [value.num for value in values(states_dict)]
                # record all the household states present at the current time point, non household states are returned as [999,999,999,999]
                key_house_vec[w_out] = [(key == "free" || key == "release") ? [999,999,999,999] : [value.fem_m, value.male_m, value.fem_w, value.male_w] for (key, value) in states_dict]
                t_week += 7  # update the week counter, next time will record results
                w_out += 1   # update the index week counter
            end

            # Stop if we've reached or passed the final time point
            if t >= t_vec[end]
                break
            end

            # Update the time point for the next (daily) output
            i_out += 1
            t_out = t_vec[i_out]
        end
        
        if m_pop + w_pop == 0   # If population has become extinct, end simulation
            break
        end
    end
    # return all the results we have recorded
    return m_vec, w_vec, free_m_vec, free_w_vec, num_house_vec, key_house_vec, track_fem_m, track_male_m, track_fem_w, track_male_w, release_fem_pop, release_male_pop
end
