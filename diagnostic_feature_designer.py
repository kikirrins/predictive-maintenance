import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np

# --- Import functions/classes from construct_datastore.py ---
try:
    # We need the datastore class and the read function
    from construct_datastore import FileEnsembleDatastore, read_member_data
except ImportError:
    print("Error: Could not import from 'construct_datastore.py'.")
    print("Make sure 'construct_datastore.py' is in the same directory or your PYTHONPATH.")
    exit()

# --- Configuration ---
data_folder = 'data_files'       # Folder with the .mat files

# --- Helper Function for Color Mapping ---
def get_color_map(unique_conditions):
    """Creates a mapping from condition labels to sequential colors (Green to Red)."""

    def get_broken_bars(condition_str):
        """Extracts the number of broken bars from the condition string."""
        if condition_str == "healthy":
            return 0
        try:
            # Assumes format like "broken_bar_X"
            return int(condition_str.split('_')[-1])
        except (ValueError, IndexError):
            # Return a high number for unknown formats to place them at the 'bad' end
            print(f"Warning: Could not parse broken bars from '{condition_str}'. Assigning default order.")
            return float('inf')

    # Sort conditions based on the number of broken bars
    sorted_conditions = sorted(unique_conditions, key=get_broken_bars)

    # Choose a sequential colormap (Green=Low/Healthy, Red=High/Most Broken)
    # RdYlGn_r goes Red -> Yellow -> Green, so we reverse it implicitly by mapping
    # healthy (0 bars) to the green end (near 1.0) and max bars to the red end (near 0.0).
    # Let's use RdYlGn_r directly. Green is near 1, Red is near 0.
    cmap = plt.cm.RdYlGn_r
    num_conditions = len(sorted_conditions)

    # Create normalized values (0 to 1) corresponding to the sorted order
    norm_values = np.linspace(0, 1, num_conditions)

    # Create the color map dictionary
    color_map = {condition: cmap(norm_values[i])
                 for i, condition in enumerate(sorted_conditions)}

    return color_map

# --- Main Script ---
if __name__ == "__main__":

    print("Initializing datastore...")
    location = os.path.abspath(data_folder)
    ens = FileEnsembleDatastore(location, '.mat')

    if ens.NumMembers == 0:
        print(f"Error: No .mat files found in {location}. Cannot plot.")
        exit()

    # Configure the datastore to read necessary variables
    ens.ReadFcn = read_member_data
    ens.SelectedVariables = ["Vib_acpi_env", "Health"] # Only need these for this plot
    ens.DataVariables = ens.SelectedVariables # Keep it simple

    # --- Read All Member Data ---
    print(f"Reading data for all {ens.NumMembers} members...")
    all_member_data = []
    ens.reset()
    while True:
        member_data = ens.read() # Reads based on ens.SelectedVariables
        if member_data is None: # End of ensemble reached or error
            # Check if end of ensemble was reached cleanly
            if ens._current_member_index >= ens.NumMembers:
                print("Finished reading all members.")
            else:
                print(f"Warning: Stopped reading prematurely after member index {ens._current_member_index-1}.")
            break
        all_member_data.append(member_data)

    if not all_member_data:
        print("Error: No data was successfully read from any members.")
        exit()

    # --- Prepare for Plotting ---
    print("Preparing data for plotting...")
    # Get unique health conditions for coloring and legend
    unique_health_conditions = sorted(list(set(md['Health'] for md in all_member_data if md is not None and 'Health' in md)))
    color_map = get_color_map(unique_health_conditions)

    # --- Create Plot ---
    fig, ax = plt.subplots(figsize=(12, 7))
    plotted_labels = set() # Keep track of labels already added to the legend

    print("Plotting signal traces...")
    for member_data in all_member_data:
        if member_data is None: continue # Skip if read failed for this member

        health = member_data.get('Health')
        signal = member_data.get('Vib_acpi_env')

        if health is None or not isinstance(signal, pd.Series) or signal.empty:
            print(f"Warning: Skipping a member due to missing Health ({health}) or invalid signal data.")
            continue

        color = color_map.get(health, 'gray') # Default to gray if condition somehow not in map
        label = health

        # Plot signal vs time index
        time_index_sec = signal.index.total_seconds()

        # Add label only once per health condition for a clean legend
        if label not in plotted_labels:
            ax.plot(time_index_sec, signal.values, color=color, label=label, alpha=0.8)
            plotted_labels.add(label)
        else:
            ax.plot(time_index_sec, signal.values, color=color, alpha=0.8)

    # --- Final Touches ---
    ax.set_title('Signal Trace: Vib_acpi_env/Data grouped by Health')
    ax.set_xlabel('Time (seconds)')
    ax.set_ylabel('Vib_acpi_env Amplitude')

    # Set x-axis limits to match the typical 10-11s window from get_signal
    # Adjust if your get_signal function uses a different window
    ax.set_xlim(10.100, 10.150) # Zoomed limits as requested

    # Add legend - place it nicely
    ax.legend(title="Health Condition", bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.grid(True)

    plt.tight_layout(rect=[0, 0, 0.85, 1]) # Adjust layout to make space for the legend outside
    plt.show()

    print("\nVisualization complete.") 