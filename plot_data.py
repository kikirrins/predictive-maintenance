import os
import matplotlib.pyplot as plt
import scipy.signal
import pandas as pd
import numpy as np # Added numpy for sample index array

# --- Import functions/classes from construct_datastore.py ---
try:
    from construct_datastore import FileEnsembleDatastore, read_member_data, calculate_envelope # Import calculate_envelope if needed directly
except ImportError:
    print("Error: Could not import from 'construct_datastore.py'.")
    print("Make sure 'construct_datastore.py' is in the same directory or your PYTHONPATH.")
    exit()

# --- Configuration ---
data_folder = 'data_files'       # Folder with the .mat files
Fs_vib = 7600                    # Sampling frequency for vibration signals (Hz)

# --- Main Script ---
if __name__ == "__main__":

    print("Initializing datastore for plotting...")
    location = os.path.abspath(data_folder)
    ens = FileEnsembleDatastore(location, '.mat')

    if ens.NumMembers == 0:
        print(f"Error: No .mat files found in {location}. Cannot plot.")
        exit()

    # Configure the datastore to read necessary variables
    ens.ReadFcn = read_member_data
    # Ensure Vib_acpi and Vib_acpi_env are selected
    ens.SelectedVariables = ["Vib_acpi", "Vib_acpi_env", "Health", "Load"]
    # Define DataVariables (can be the same as SelectedVariables or a superset)
    ens.DataVariables = ens.SelectedVariables # Keep it simple for now

    print("Reading data for the first member...")
    ens.reset() # Ensure reading starts from the first member
    member_data = ens.read()

    if member_data is None:
        print("Failed to read the first member. Cannot plot.")
        exit()

    # Extract data for plotting
    vib_acpi_signal = member_data.get('Vib_acpi')
    vib_acpi_env_signal = member_data.get('Vib_acpi_env')
    member_filename = os.path.basename(ens.files[ens._current_member_index-1]) # Get filename of member just read

    # --- Create Figure with Two Subplots ---
    fig, axes = plt.subplots(2, 1, figsize=(10, 10)) # 2 rows, 1 column

    # --- Plot 1: Power Spectrum of Vib_acpi (on axes[0]) ---
    print("Plotting Power Spectrum (kHz)...")
    if isinstance(vib_acpi_signal, pd.Series) and not vib_acpi_signal.empty:
        # Use axes[0].psd to plot on the first subplot
        axes[0].psd(vib_acpi_signal.values, NFFT=1024, Fs=Fs_vib)
        axes[0].set_title(f'Power Spectrum of Vib_acpi (Member: {member_filename})')

        # Convert MATLAB annotation coordinates (relative figure units) to data coordinates.
        arrow_start_x_hz = 1500
        arrow_start_y_db = -20
        arrow_end_x_hz = 1200
        arrow_end_y_db = -30

        axes[0].annotate('Fault frequency region of interest',
                     xy=(arrow_end_x_hz, arrow_end_y_db), # Point to this location (in Hz)
                     xytext=(arrow_start_x_hz, arrow_start_y_db), # Text starts here
                     arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8))

        # Set x-axis to kHz
        current_ticks_hz = axes[0].get_xticks()
        # Filter out negative ticks if they appear due to autoscaling
        current_ticks_hz = [tick for tick in current_ticks_hz if tick >= 0]
        axes[0].set_xticks(current_ticks_hz) # Keep the tick positions
        axes[0].set_xticklabels([f'{x/1000:.1f}' for x in current_ticks_hz]) # Format as kHz
        axes[0].set_xlabel('Frequency (kHz)') # Update label
        axes[0].set_ylabel('Power Spectral Density (dB/Hz)') # Ensure y label is descriptive
        axes[0].grid(True)

    else:
        print("Could not plot Vib_acpi spectrum: Signal is missing, empty, or not a Pandas Series.")
        axes[0].set_title(f'Vib_acpi Power Spectrum Data Unavailable (Member: {member_filename})')

    # --- Plot 2: Envelope of Bandpassed Vib_acpi (on axes[1]) ---
    print("Plotting Envelope and Bandpassed Signal...")
    # Check if both the original signal (to be bandpassed) and the pre-calculated envelope are available
    if (isinstance(vib_acpi_signal, pd.Series) and not vib_acpi_signal.empty and
        isinstance(vib_acpi_env_signal, pd.Series) and not vib_acpi_env_signal.empty):

         # Re-calculate the bandpassed signal 'y' for plotting purposes
         # (This matches the first step within calculate_envelope)
         sos = scipy.signal.butter(5, [900, 1300], btype='bandpass', fs=Fs_vib, output='sos')
         y_bandpassed = scipy.signal.sosfiltfilt(sos, vib_acpi_signal.values)

         # Determine number of samples to plot (matching MATLAB example)
         num_samples_to_plot = 800
         plot_indices = np.arange(num_samples_to_plot)

         if len(y_bandpassed) < num_samples_to_plot:
              # Adjust if the signal segment is shorter
              num_samples_to_plot = len(y_bandpassed)
              plot_indices = np.arange(num_samples_to_plot)
              plot_data_bandpassed = y_bandpassed
              plot_data_envelope = vib_acpi_env_signal.values[:num_samples_to_plot] # Envelope should have same length
         else:
              plot_data_bandpassed = y_bandpassed[:num_samples_to_plot]
              plot_data_envelope = vib_acpi_env_signal.values[:num_samples_to_plot]


         # Plot the bandpassed signal (lighter color/alpha)
         axes[1].plot(plot_indices, plot_data_bandpassed, label='Bandpassed Signal', alpha=0.7)
         # Plot the envelope (usually more prominent)
         axes[1].plot(plot_indices, plot_data_envelope, label='Envelope')


         axes[1].set_title(f'Bandpassed Vib_acpi Signal and Envelope (First {num_samples_to_plot} Samples)')
         axes[1].set_xlabel('Sample Index')
         axes[1].set_ylabel('Amplitude')
         # Set axis limits similar to MATLAB example
         axes[1].axis([0, num_samples_to_plot, -0.5, 0.5])
         axes[1].grid(True)
         axes[1].legend() # Add legend

    else:
        print("Could not plot Envelope/Bandpassed Signal: Required data (Vib_acpi or Vib_acpi_env) is missing or invalid.")
        axes[1].set_title(f'Vib_acpi Bandpassed/Envelope Data Unavailable\n(Member: {member_filename})')
        axes[1].grid(True)

    # --- Final Display ---
    plt.tight_layout() # Adjust layout to prevent overlapping titles/labels
    plt.show()

    print("Plotting complete.") 