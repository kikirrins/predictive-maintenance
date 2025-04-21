import os
import scipy.io
import numpy as np

# Function to load mat and extract the primary data structure
def load_mat_struct(filename):
    """Loads a .mat file, expecting a single main struct variable.
       Uses squeeze_me=True to simplify array shapes and
       struct_as_record=False to load structs as objects with attributes.
    """
    try:
        # Load the .mat file
        mat_data = scipy.io.loadmat(filename, squeeze_me=True, struct_as_record=False)

        # Filter out MATLAB's metadata keys
        data_keys = [k for k in mat_data if not k.startswith('__')]

        if not data_keys:
            raise ValueError(f"No data variable found in {filename}")
        if len(data_keys) > 1:
            print(f"Warning: Multiple variables found in {filename}. Using '{data_keys[0]}'.")

        # Assume the first non-metadata key contains the main data struct
        main_key = data_keys[0]
        return mat_data[main_key]
    except FileNotFoundError:
        print(f"Error: File {filename} not found.")
        return None
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

# --- Configuration ---
# Number of experiments per condition to process
# (Use 2 for the reduced dataset, up to 10 for the full dataset)
numExperiments = 2

# Original data file names from the unzipped archive
# These correspond to different health conditions (number of broken bars)
files = [
    "struct_rs_R1.mat",    # Healthy (0 broken bars)
    "struct_r1b_R1.mat",   # 1 broken bar
    "struct_r2b_R1.mat",   # 2 broken bars
    "struct_r3b_R1.mat",   # 3 broken bars
    "struct_r4b_R1.mat",   # 4 broken bars
]

# Corresponding health condition labels
health = [
    "healthy",
    "broken_bar_1",
    "broken_bar_2",
    "broken_bar_3",
    "broken_bar_4",
]

# Sampling frequencies (Hz)
Fs_vib = 7600   # Vibration signals
Fs_elec = 50000 # Electrical signals

# Target directory for the new ensemble member .mat files
output_folder = 'data_files'

# --- Data Extraction Process ---

# Create the output directory if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    print(f"Created directory: {output_folder}")

# Iterate over each original data file (representing a health condition)
for i, (filename, health_condition) in enumerate(zip(files, health)):
    print(f'Processing data file {filename} ({health_condition})')

    # Load the main data structure from the .mat file
    dataset = load_mat_struct(filename)

    # If loading failed, skip to the next file
    if dataset is None:
        continue

    try:
        # Get the names of the load levels (e.g., 'torque05', 'torque10', ...)
        # These are fields within the loaded 'dataset' structure object
        if not hasattr(dataset, '_fieldnames'):
             print(f"Warning: Could not find load level fields in {filename}. Skipping.")
             continue
        load_levels = dataset._fieldnames

        # Iterate over each load level (torque) found in the file
        for load_level in load_levels:
            # Get the array of experiments for the current load level
            # Accessing the field using getattr()
            experiments_array = getattr(dataset, load_level)

            # Determine how many experiments are actually available for this condition
            # Handle case where experiments_array might be a single struct if only 1 experiment exists
            if isinstance(experiments_array, np.ndarray):
                num_available_experiments = len(experiments_array)
            elif hasattr(experiments_array, '_fieldnames'): # Check if it's a single struct-like object
                 num_available_experiments = 1
                 experiments_array = [experiments_array] # Wrap in a list for consistent iteration
            else:
                 print(f"Warning: Unexpected data format for experiments under {load_level} in {filename}. Skipping load level.")
                 continue


            # Limit processing to the desired number of experiments or available ones
            experiments_to_process = min(numExperiments, num_available_experiments)
            if experiments_to_process < numExperiments:
                 print(f"\tNote: Only {num_available_experiments} experiments available for {load_level}, processing {experiments_to_process}.")


            # Iterate through each selected experiment
            for k in range(experiments_to_process):
                # Get the data for the current experiment (k-th experiment)
                current_experiment = experiments_array[k]

                # Prepare a dictionary to hold the data for the new .mat file
                data_to_save = {}

                # Get the names of the signals (fields within the experiment structure)
                if not hasattr(current_experiment, '_fieldnames'):
                     print(f"Warning: Could not find signal fields in experiment {k+1} for {load_level} in {filename}. Skipping experiment.")
                     continue
                signal_names = current_experiment._fieldnames

                # Copy each signal from the experiment data to the save dictionary
                for signal_name in signal_names:
                    signal_value = getattr(current_experiment, signal_name)
                    data_to_save[signal_name] = signal_value

                # Add metadata (operating conditions, constants) to the dictionary
                data_to_save['Health'] = health_condition
                data_to_save['Load'] = str(load_level) # Ensure load level is saved as a string
                data_to_save['Fs_vib'] = Fs_vib
                data_to_save['Fs_elec'] = Fs_elec

                # Construct the output filename according to the specified format
                # Example: rotor0b_torque05_experiment01.mat
                output_filename_base = f'rotor{i}b_{load_level}_experiment{k+1:02d}'
                output_filename = os.path.join(output_folder, f'{output_filename_base}.mat')

                print(f'\tCreating the member data file {os.path.basename(output_filename)}')

                # Save the dictionary to a new .mat file
                # `do_compression=True` is generally good practice.
                # `scipy.io.savemat` defaults to MATLAB v5 format unless data requires newer.
                # For strict v7.3 compatibility (e.g., large variables), consider the 'hdf5storage' library.
                scipy.io.savemat(output_filename, data_to_save, do_compression=True)

    except AttributeError as ae:
         print(f"Error accessing attributes in {filename} (potentially unexpected structure): {ae}")
    except Exception as e:
        print(f"An unexpected error occurred while processing {filename}: {e}")

print("\nFinished extracting ensemble member data.")
