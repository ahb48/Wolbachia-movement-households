using CSV
using DataFrames
using Statistics

# Function to generate the file name
function generate_filename(type, i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
    return "$(type)_$(i)_$(alpha_str)_$(gamma_str)_$(rho_str)_$(tau_str)_$(H)__$(rel_t)_$(rel_size)_$(rel_type).csv"
end

# Function to check if the file exists and read it
function check_and_read_file(type, i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
    filename = generate_filename(type, i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
    if isfile(joinpath("Initial_paper_results", "$(JOB_ID)", "$(JOB_ID)",filename))
        #println("File $filename exists. Reading the file.")
        df = CSV.read(joinpath("Initial_paper_results", "$(JOB_ID)", "$(JOB_ID)",filename), DataFrame)
        return df
    else
        #println("File $filename does not exist.")
        return nothing
    end
end

function check_file_exist(alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, num_simulations, JOB_ID)
    # Initialize an array to accumulate the sum of arrays
    sum_array = nothing
    ct = 0
    
    for i in 1:num_simulations
        df_f = check_and_read_file("fem_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_m = check_and_read_file("male_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        if df_f !== nothing && df_m !== nothing
            # Convert DataFrame to array
            arr_f = Matrix(df_f)
            arr_m = Matrix(df_m)
            sum = arr_f .+ arr_m
        
            # Accumulate the sum of arrays
            if sum_array === nothing
                sum_array = sum
            else
                sum_array .+= sum
            end
            ct += 1
        else
            println("DataFrame for simulation $i does not exist. Skipping.")
        end
    end
    return sum_array, ct
end

function hhold_count(alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, num_simulations, JOB_ID)
    # Initialize an array to accumulate the sum of arrays
    no_count = nothing
    ct = 0

    for i in 1:num_simulations
        df_f = check_and_read_file("fem_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_m = check_and_read_file("male_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_fw = check_and_read_file("fem_w_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_mw = check_and_read_file("male_w_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        if df_f !== nothing && df_m !== nothing && df_fw !== nothing && df_mw !== nothing
            # Convert DataFrame to array
            arr_f = Matrix(df_f)
            arr_m = Matrix(df_m)
            arr_fw = Matrix(df_fw)
            arr_mw = Matrix(df_mw)
            # return the number of households containing no wildtypes at each time step
            if no_count === nothing 
                print("Initializing no_wild_count")
                no_count = [count(==(0), row) for row in eachrow(arr_f.+arr_m+arr_fw+arr_mw)]
            else
                no_count .+= [count(==(0), row) for row in eachrow(arr_f.+arr_m+arr_fw+arr_mw)]
                print(no_count)
            end
            ct += 1
        else
            println("DataFrame for simulation $i does not exist. Skipping.")
        end
    end
    return no_count, ct
end

function hhold_wild_count(alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, num_simulations, JOB_ID)
    # Initialize an array to accumulate the sum of arrays
    no_wild_count = nothing
    ct = 0

    for i in 1:num_simulations
        df_f = check_and_read_file("fem_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_m = check_and_read_file("male_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        if df_f !== nothing && df_m !== nothing
            # Convert DataFrame to array
            arr_f = Matrix(df_f)
            arr_m = Matrix(df_m)
            # return the number of households containing no wildtypes at each time step
            if no_wild_count === nothing 
                print("Initializing no_wild_count")
                no_wild_count = [count(==(0), row) for row in eachrow(arr_f.+arr_m)]
            else
                no_wild_count .+= [count(==(0), row) for row in eachrow(arr_f.+arr_m)]
                print(no_wild_count)
            end
            ct += 1
        else
            println("DataFrame for simulation $i does not exist. Skipping.")
        end
    end
    return no_wild_count, ct
end

function var_hholds(alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, num_simulations, JOB_ID)
    # Initialize arrays to accumulate the arrays
    arr_fem_m = nothing
    arr_male_m = nothing
    arr_fem_w = nothing
    arr_male_w = nothing
    arr_w_tot = nothing
    arr_m_tot = nothing

    for i in 1:num_simulations
        df_fem_m = Check.check_and_read_file("fem_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_male_m = Check.check_and_read_file("male_m_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_fem_w = Check.check_and_read_file("fem_w_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        df_male_w = Check.check_and_read_file("male_w_hholds", i, alpha_str, gamma_str, rho_str, tau_str, H, rel_t, rel_size, rel_type, JOB_ID)
        # Convert DataFrame to array
        if arr_fem_m === nothing
            arr_fem_m = Matrix(df_fem_m)
        else  # horizontally concatenates the arrays so can calculate the variance over all the sims for each time step
            arr_fem_m = hcat(arr_fem_m, Matrix(df_fem_m))
        end

        if arr_male_m === nothing
            arr_male_m = Matrix(df_male_m)
        else
            arr_male_m = hcat(arr_male_m, Matrix(df_male_m))
        end

        if arr_fem_w === nothing
            arr_fem_w = Matrix(df_fem_w)
        else
            arr_fem_w = hcat(arr_fem_w, Matrix(df_fem_w))
        end

        if arr_male_w === nothing
            arr_male_w = Matrix(df_male_w)
        else
            arr_male_w = hcat(arr_male_w, Matrix(df_male_w))
        end

        df_w_tot = df_male_w .+ df_fem_w
        if arr_w_tot === nothing
            arr_w_tot = Matrix(df_w_tot)
        else
            arr_w_tot = hcat(arr_w_tot, Matrix(df_w_tot))
        end

        df_m_tot = df_male_m .+ df_fem_m
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

    return row_var_fem_m, row_var_male_m, row_var_fem_w, row_var_male_w, row_var_w_tot, row_var_m_tot
end
