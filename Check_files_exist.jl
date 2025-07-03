# import libraries
using CSV
using DataFrames
using Statistics

function generate_filename(type, i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
    """Function to generate the file name based on the parameters."""
    return "$(type)_$(i)_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"
end

function check_and_read_file(type, i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
    """Function to check if the file exists and read it into a DataFrame."""
    filename = generate_filename(type, i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
    # if the file exists
    if isfile(joinpath("Initial_paper_results", "$(JOB_ID)", "$(JOB_ID)",filename))
        #println("File $filename exists. Reading the file.")
        # read the file into a DataFrame
        df = CSV.read(joinpath("Initial_paper_results", "$(JOB_ID)", "$(JOB_ID)",filename), DataFrame)
        return df   # return the DataFrame
    else  # if the file does not exist
        #println("File $filename does not exist.")
        return nothing    # don't return anything
    end
end

function check_file_exist(alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, num_simulations, JOB_ID)
    """Function to check if the files exist for all simulations and sum the number of wildtypes 
    in the household and the number of simulations."""
    # Initialize an array to accumulate the sum of arrays
    sum_array = nothing  # Initialize array sum as nothing
    ct = 0 # initialize count as 0
    
    for i in 1:num_simulations # loop through the number of simulations
        # Check if the files exist for male and female wildtypes in the households
        df_f = check_and_read_file("fem_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_m = check_and_read_file("male_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        if df_f !== nothing && df_m !== nothing  # if both DataFrames are not nothing
            # Convert DataFrame to array
            arr_f = Matrix(df_f)
            arr_m = Matrix(df_m)
            sum = arr_f .+ arr_m  # sum the male and female arrays
        
            # Accumulate the sum of arrays
            if sum_array === nothing   # if sum_array is nothing, initialize it
                sum_array = sum
            else
                sum_array .+= sum      # otherwise, add the sum to the existing sum_array
            end
            ct += 1   # increment the count of simulations with data
        else
            println("DataFrame for simulation $i does not exist. Skipping.")
        end
    end
    return sum_array, ct  # return array sum and count of simulations
end

function hhold_count(alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, num_simulations, JOB_ID)
    """Function to count the number of households containing no mosquitoes at each time step. 
    First checks if the files exist for all simulations."""
    # Initialize an array to accumulate the sum of arrays
    no_count = nothing
    ct = 0 # initialize count as 0

    for i in 1:num_simulations # loop through the number of simulations
        # Check if the files exist
        df_f = check_and_read_file("fem_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_m = check_and_read_file("male_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_fw = check_and_read_file("fem_w_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_mw = check_and_read_file("male_w_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        # if all DataFrames are not nothing
        if df_f !== nothing && df_m !== nothing && df_fw !== nothing && df_mw !== nothing
            # Convert DataFrame to array
            arr_f = Matrix(df_f)
            arr_m = Matrix(df_m)
            arr_fw = Matrix(df_fw)
            arr_mw = Matrix(df_mw)
            # return the number of households containing no wildtypes at each time step
            if no_count === nothing # if no_count is nothing, initialize it
                #print("Initializing no_wild_count")
                # count the number of households containing no wildtypes at each time step
                no_count = [count(==(0), row) for row in eachrow(arr_f.+arr_m+arr_fw+arr_mw)]
            else # otherwise, add the count to the existing no_count
                no_count .+= [count(==(0), row) for row in eachrow(arr_f.+arr_m+arr_fw+arr_mw)]
                print(no_count)
            end
            ct += 1  # increment the count of simulations with data
        else  # if any DataFrame is nothing
            println("DataFrame for simulation $i does not exist. Skipping.")
        end
    end
    return no_count, ct  # return the count of households with no wildtypes and the count of simulations with data
end

function hhold_wild_count(alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, num_simulations, JOB_ID)
    """Function to count the number of households containing no wildtypes at each time step.
    First checks if the files exist for all simulations."""
    # Initialize an array to accumulate the sum of arrays
    no_wild_count = nothing
    ct = 0  # initialize count as 0

    for i in 1:num_simulations  # loop through the number of simulations
        # Check if the files exist
        df_f = check_and_read_file("fem_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_m = check_and_read_file("male_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        if df_f !== nothing && df_m !== nothing  # if both DataFrames are not nothing
            # Convert DataFrame to array
            arr_f = Matrix(df_f)
            arr_m = Matrix(df_m)
            # return the number of households containing no wildtypes at each time step
            if no_wild_count === nothing   # if no_wild_count is nothing, initialize it
                #print("Initializing no_wild_count")
                # count the number of households containing no wildtypes at each time step
                no_wild_count = [count(==(0), row) for row in eachrow(arr_f.+arr_m)]
            else
                no_wild_count .+= [count(==(0), row) for row in eachrow(arr_f.+arr_m)]
                #print(no_wild_count)
            end
            ct += 1  # increment the count of simulations with data
        else  # if any DataFrame is nothing
            println("DataFrame for simulation $i does not exist. Skipping.")
        end
    end
    return no_wild_count, ct  # return the count of households with no wildtypes and the count of simulations with data
end

function var_hholds(alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, num_simulations, JOB_ID)
    """Function to calculate the variance in the number of households containing wildtypes and 
    Wolbachia-infected mosquitoes at each time step. First checks if the files exist for all simulations."""
    # Initialize arrays to accumulate the arrays
    arr_fem_m = nothing
    arr_male_m = nothing
    arr_fem_w = nothing
    arr_male_w = nothing
    arr_w_tot = nothing
    arr_m_tot = nothing

    for i in 1:num_simulations  # loop through the number of simulations
        # Check if the files exist
        df_fem_m = Check.check_and_read_file("fem_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_male_m = Check.check_and_read_file("male_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_fem_w = Check.check_and_read_file("fem_w_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_male_w = Check.check_and_read_file("male_w_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        # Convert DataFrame to array
        if arr_fem_m === nothing  # if arr_fem_m is nothing, initialize it
            arr_fem_m = Matrix(df_fem_m)
        else  # horizontally concatenates the arrays so can calculate the variance over all the sims for each time step
            arr_fem_m = hcat(arr_fem_m, Matrix(df_fem_m))
        end

        if arr_male_m === nothing  # if arr_male_m is nothing, initialize it
            arr_male_m = Matrix(df_male_m)
        else
            arr_male_m = hcat(arr_male_m, Matrix(df_male_m)) 
        end

        if arr_fem_w === nothing # if arr_fem_w is nothing, initialize it
            arr_fem_w = Matrix(df_fem_w)
        else
            arr_fem_w = hcat(arr_fem_w, Matrix(df_fem_w))
        end

        if arr_male_w === nothing # if arr_male_w is nothing, initialize it
            arr_male_w = Matrix(df_male_w)
        else
            arr_male_w = hcat(arr_male_w, Matrix(df_male_w))
        end

        df_w_tot = df_male_w .+ df_fem_w  # Total number of Wolbachia-infected mosquitoes in the household
        if arr_w_tot === nothing
            arr_w_tot = Matrix(df_w_tot)
        else
            arr_w_tot = hcat(arr_w_tot, Matrix(df_w_tot))
        end

        df_m_tot = df_male_m .+ df_fem_m  # Total number of wildtype mosquitoes in the household
        if arr_m_tot === nothing
            arr_m_tot = Matrix(df_m_tot)
        else
            arr_m_tot = hcat(arr_m_tot, Matrix(df_m_tot))
        end
    end

    # Calculate the variance for each row
    row_var_fem_m = [var(row) for row in eachrow(arr_fem_m)]
    row_var_male_m = [var(row) for row in eachrow(arr_male_m)]
    row_var_fem_w = [var(row) for row in eachrow(arr_fem_w)]
    row_var_male_w = [var(row) for row in eachrow(arr_male_w)]
    row_var_w_tot = [var(row) for row in eachrow(arr_w_tot)]
    row_var_m_tot = [var(row) for row in eachrow(arr_m_tot)]
    # Return the variances for each type of mosquito
    return row_var_fem_m, row_var_male_m, row_var_fem_w, row_var_male_w, row_var_w_tot, row_var_m_tot
end
