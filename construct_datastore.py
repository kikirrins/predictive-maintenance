import os
import glob
import scipy.io
import scipy.signal
import numpy as np
import pandas as pd
# Optional: import matplotlib.pyplot as plt

# --- Helper Functions ---

def write_member_data(filename, data_dict):
    """
    Appends data to an existing .mat file or creates a new one.

    Loads the existing file, updates it with new data, and saves it back.
    Note: This is less efficient than MATLAB's direct append.

    Args:
        filename (str): Path to the .mat file.
        data_dict (dict): Dictionary containing data to add/update.
    """
    try:
        # Load existing data if the file exists
        try:
            existing_data = scipy.io.loadmat(filename)
            # Remove MATLAB metadata keys before merging
            existing_data = {k: v for k, v in existing_data.items() if not k.startswith('__')}
        except FileNotFoundError:
            existing_data = {}
        except Exception as e:
            print(f"Warning: Could not load existing data from {filename}. Overwriting if keys conflict. Error: {e}")
            existing_data = {}

        # Update with new data
        existing_data.update(data_dict)

        # Save the combined data
        scipy.io.savemat(filename, existing_data, do_compression=True)

    except Exception as e:
        print(f"Error writing to {filename}: {e}")


def get_signal(mat_data, signame, fs):
    """
    Extracts a 1.0 second portion of a signal from loaded .mat data,
    starting 10 seconds into the measurement.

    Args:
        mat_data (dict): Data loaded from a .mat file.
        signame (str): Name of the signal variable.
        fs (float): Sampling frequency of the signal.

    Returns:
        pd.Series: A pandas Series containing the signal segment,
                   indexed by time in seconds. Returns None on error.
    """
    try:
        signal_data = mat_data[signame].flatten() # Ensure 1D array
        n = len(signal_data)
        t_full = np.arange(n) / fs

        # Find indices for the time window [10.0, 11.0]
        indices = np.where((t_full >= 10.0) & (t_full <= 11.0))[0]

        if len(indices) == 0:
             print(f"Warning: No data found in the 10s-11s window for {signame}.")
             # Return an empty Series with appropriate index type if needed, or None
             return pd.Series([], dtype=float, index=pd.to_timedelta([], unit='s'), name=signame)


        time_segment = t_full[indices]
        data_segment = signal_data[indices]

        # Create a pandas Series with a TimedeltaIndex
        time_index = pd.to_timedelta(time_segment, unit='s')
        series = pd.Series(data_segment, index=time_index, name=signame)
        return series

    except KeyError:
        print(f"Error: Signal '{signame}' not found in the data.")
        return None
    except Exception as e:
        print(f"Error processing signal {signame}: {e}")
        return None

def calculate_envelope(signal_data, fs, band=None):
    """Calculates the envelope of a signal, optionally after bandpass filtering."""
    if band:
        # Design Butterworth bandpass filter
        sos = scipy.signal.butter(5, band, btype='bandpass', fs=fs, output='sos')
        # Apply filter
        filtered_signal = scipy.signal.sosfiltfilt(sos, signal_data)
    else:
        filtered_signal = signal_data

    # Calculate envelope using Hilbert transform
    analytic_signal = scipy.signal.hilbert(filtered_signal)
    envelope = np.abs(analytic_signal)
    return envelope

def calculate_envelope_spectrum(signal_series, fs, band):
    """
    Calculates the envelope spectrum of a signal within a specific band.
    Approximates MATLAB's envspectrum(..., 'Method', 'hilbert', 'Band', band).
    """
    signal_data = signal_series.values
    # 1. Bandpass filter
    sos = scipy.signal.butter(5, band, btype='bandpass', fs=fs, output='sos')
    filtered_signal = scipy.signal.sosfiltfilt(sos, signal_data)

    # 2. Calculate envelope
    analytic_signal = scipy.signal.hilbert(filtered_signal)
    envelope = np.abs(analytic_signal)

    # 3. Calculate spectrum (PSD) of the envelope using Welch
    # Adjust nperseg for desired frequency resolution vs. variance trade-off
    nperseg = min(len(envelope), 2048) # Example segment length
    freqs, psd = scipy.signal.welch(envelope, fs=fs, nperseg=nperseg, scaling='density')

    # Return as a DataFrame similar to MATLAB's [F, ES] output
    spectrum_df = pd.DataFrame({'Frequency': freqs, 'PSD': psd})
    return spectrum_df


def read_member_data(filename, variables):
    """
    Reads specified variables from a single member .mat file.
    Mimics the behavior of the MATLAB readMemberData function, including
    generating synthetic signals on the fly.

    Args:
        filename (str): Path to the .mat file.
        variables (list): List of variable names to read/generate.

    Returns:
        pd.Series: A pandas Series containing the requested data for one member.
                   Values that are arrays/dataframes are kept as objects.
    """
    try:
        # Load the entire .mat file first (less efficient than matfile but simpler)
        mat_data = scipy.io.loadmat(filename, squeeze_me=True, struct_as_record=False)
        # Extract sampling frequencies, assuming they exist
        fs_elec = mat_data.get('Fs_elec', None)
        fs_vib = mat_data.get('Fs_vib', None)
        if fs_elec is None or fs_vib is None:
             print(f"Warning: Sampling frequencies not found in {filename}. Synthetic signals may fail.")


        member_data = {}
        for var in variables:
            val = None
            try:
                if var in ['Health', 'Load']:
                    # Condition variables
                    val = mat_data[var]
                elif var in ['Va', 'Vb', 'Vc', 'Ia', 'Ib', 'Ic']:
                    # Electrical signals (extract 1s segment)
                    val = get_signal(mat_data, var, fs_elec)
                elif var in ['Vib_acpi', 'Vib_carc', 'Vib_acpe', 'Vib_axial', 'Vib_base', 'Trigger']:
                    # Vibration signals (extract 1s segment)
                    val = get_signal(mat_data, var, fs_vib)
                elif var == 'Vib_acpi_env':
                    # Synthetic envelope for Vib_acpi
                    base_sig_name = 'Vib_acpi'
                    base_signal = get_signal(mat_data, base_sig_name, fs_vib)
                    if base_signal is not None:
                         # Calculate envelope after bandpass filtering [900, 1300] Hz
                         envelope = calculate_envelope(base_signal.values, fs_vib, band=[900, 1300])
                         # Create Series with the same time index as the base signal segment
                         val = pd.Series(envelope, index=base_signal.index, name=var)

                elif var == 'Ia_env_ps':
                     # Synthetic envelope spectrum for Ia
                     base_sig_name = 'Ia'
                     base_signal = get_signal(mat_data, base_sig_name, fs_elec)
                     if base_signal is not None:
                          # Calculate envelope spectrum for band [900, 1300] Hz
                          val = calculate_envelope_spectrum(base_signal, fs_elec, band=[900, 1300])

                # Add cases for other synthetic signals like 'Ia_env' if needed
                # elif var == 'Ia_env':
                #     # ... similar logic using fs_elec ...

                else:
                    # Try reading other variables directly (e.g., features added later)
                    if var in mat_data:
                         val = mat_data[var]
                    else:
                         print(f"Warning: Variable '{var}' not found directly in {filename} and not a known synthetic signal.")

            except KeyError:
                 print(f"Warning: Variable '{var}' or its dependencies not found in {filename}.")
            except Exception as e:
                 print(f"Error processing variable '{var}' in {filename}: {e}")


            # Store the value (even if None)
            member_data[var] = val

        # Return as a pandas Series (allows different data types, including objects)
        return pd.Series(member_data)

    except FileNotFoundError:
        print(f"Error: File not found: {filename}")
        return None
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        return None


# --- FileEnsembleDatastore Class ---

class FileEnsembleDatastore:
    """
    Mimics MATLAB's fileEnsembleDatastore for managing collections of .mat files.
    """
    def __init__(self, location, extension='.mat'):
        self.location = os.path.abspath(location)
        self.extension = extension
        if not os.path.isdir(self.location):
            raise ValueError(f"Location does not exist or is not a directory: {self.location}")

        self.files = sorted(glob.glob(os.path.join(self.location, f'*{self.extension}')))
        if not self.files:
            print(f"Warning: No files with extension '{self.extension}' found in {self.location}")

        self.ReadFcn = None
        self.WriteToMemberFcn = None
        self.DataVariables = []
        self.ConditionVariables = []
        self.SelectedVariables = []
        self._current_member_index = 0 # For iterating with read()

    @property
    def NumMembers(self):
        return len(self.files)

    def reset(self):
        """Resets the iteration index for the read() method."""
        self._current_member_index = 0

    def read(self):
        """
        Reads data from the next member file using the assigned ReadFcn
        and SelectedVariables. Returns None if all members have been read
        or if ReadFcn is not set.
        """
        if self.ReadFcn is None:
            print("Error: ReadFcn is not set.")
            return None
        if not self.SelectedVariables:
             print("Warning: SelectedVariables is empty. Nothing to read.")
             return pd.Series(dtype=object) # Return empty Series consistent with read_member_data

        if self._current_member_index >= self.NumMembers:
            print("End of ensemble reached.")
            return None # Or raise StopIteration? Let's return None for now.

        current_file = self.files[self._current_member_index]
        print(f"Reading member: {os.path.basename(current_file)}")
        try:
            data = self.ReadFcn(current_file, self.SelectedVariables)
            self._current_member_index += 1
            return data
        except Exception as e:
            print(f"Error reading member {current_file} using ReadFcn: {e}")
            # Decide whether to advance index or not. Let's advance to avoid getting stuck.
            self._current_member_index += 1
            return None # Indicate error for this member

    def write(self, data_dict):
         """
         Writes data to the 'current' member file using WriteToMemberFcn.
         Note: The concept of 'current' member for writing might need refinement
         depending on the workflow. This implementation writes to the member
         that *would* be read next.
         """
         if self.WriteToMemberFcn is None:
             print("Error: WriteToMemberFcn is not set.")
             return
         if self._current_member_index >= self.NumMembers:
             print("Error: Cannot write, end of ensemble reached or ensemble is empty.")
             return

         current_file = self.files[self._current_member_index]
         print(f"Writing to member: {os.path.basename(current_file)}")
         try:
             self.WriteToMemberFcn(current_file, data_dict)
         except Exception as e:
             print(f"Error writing to member {current_file} using WriteToMemberFcn: {e}")


    def __repr__(self):
        repr_str = f"{self.__class__.__name__} with properties:\n\n"
        repr_str += f"               Location: {self.location}\n"
        repr_str += f"              Extension: {self.extension}\n"
        repr_str += f"             NumMembers: {self.NumMembers}\n"
        repr_str += f"                ReadFcn: {'Set' if self.ReadFcn else 'Not set'}\n"
        repr_str += f"       WriteToMemberFcn: {'Set' if self.WriteToMemberFcn else 'Not set'}\n"
        repr_str += f"          DataVariables: {len(self.DataVariables)} variables\n"
        repr_str += f"     ConditionVariables: {len(self.ConditionVariables)} variables\n"
        repr_str += f"      SelectedVariables: {len(self.SelectedVariables)} variables\n"
        # Simple display for lists if short, otherwise just count
        selected_vars_display = self.SelectedVariables if len(self.SelectedVariables) < 8 else f"[{len(self.SelectedVariables)} variables]"
        repr_str += f"    SelectedVariables: {selected_vars_display}\n"
        repr_str += f"     CurrentFileIndex: {self._current_member_index}\n"

        # Add more details if needed, e.g., list first few files
        if self.files:
             repr_str += f"                Files: [ '{os.path.basename(self.files[0])}', ... ]\n"
        else:
             repr_str += f"                Files: []\n"

        return repr_str

# --- Main Execution Example --- (if run as script)
if __name__ == "__main__":

    # Define the folder containing the extracted .mat files
    data_folder = 'data_files' # Assumes this folder is in the same directory

    # 1. Construct the File Ensemble Datastore
    print("Constructing FileEnsembleDatastore...")
    location = os.path.abspath(data_folder)
    ens = FileEnsembleDatastore(location, '.mat')

    if ens.NumMembers == 0:
        print(f"Error: No .mat files found in {location}. Please run prepare_data.py first.")
    else:
        # 2. Assign Read and Write Functions
        print("Assigning Read/Write functions...")
        ens.ReadFcn = read_member_data
        ens.WriteToMemberFcn = write_member_data # Using the helper function defined above

        # 3. Define Variables
        print("Defining variables...")
        ens.DataVariables = [
            "Va", "Vb", "Vc", "Ia", "Ib", "Ic",
            "Vib_acpi", "Vib_carc", "Vib_acpe", "Vib_axial", "Vib_base", "Trigger"
        ]
        ens.ConditionVariables = ["Health", "Load"]

        # Initially select a subset of variables
        ens.SelectedVariables = ["Ia", "Vib_acpi", "Health", "Load"]

        # Add synthetic signals to DataVariables and SelectedVariables
        synthetic_data_vars = ["Vib_acpi_env", "Ia_env_ps"] # Add 'Ia_env' here if implemented
        ens.DataVariables.extend(synthetic_data_vars)
        ens.SelectedVariables.extend(synthetic_data_vars)


        # 4. Examine the ensemble configuration
        print("\nEnsemble Configuration:")
        print(ens)

        # 5. Read the first member
        print("\nReading the first member:")
        ens.reset() # Ensure we start from the beginning
        member_data = ens.read()

        if member_data is not None:
            print("\nData read from the first member:")
            # Displaying the Series - might be verbose for large data
            # print(member_data)
            # More concise display: print keys and types/shapes
            for k, v in member_data.items():
                 if isinstance(v, (pd.Series, pd.DataFrame, np.ndarray)):
                      print(f"  {k}: {type(v).__name__} with shape {getattr(v, 'shape', 'N/A')}")
                 else:
                      print(f"  {k}: {v} ({type(v).__name__})")


            # --- Optional: Plotting for Verification ---
            # Requires matplotlib: pip install matplotlib
            # Uncomment the import at the top and this section to enable plotting

            # print("\nGenerating plots for verification...")
            # import matplotlib.pyplot as plt

            # # Plot Power Spectrum of Vib_acpi
            # vib_acpi_signal = member_data.get('Vib_acpi')
            # if vib_acpi_signal is not None and not vib_acpi_signal.empty:
            #     fs_vib_read = 7600 # Assume known or read from file if stored differently
            #     plt.figure(figsize=(10, 5))
            #     plt.psd(vib_acpi_signal.values, NFFT=1024, Fs=fs_vib_read) # Basic PSD plot
            #     # Or use Welch for better estimate:
            #     # freqs, psd_welch = scipy.signal.welch(vib_acpi_signal.values, fs=fs_vib_read, nperseg=1024)
            #     # plt.semilogy(freqs, psd_welch)
            #     plt.title(f'Power Spectrum of Vib_acpi (Member: {os.path.basename(ens.files[0])})')
            #     plt.xlabel('Frequency (Hz)')
            #     plt.ylabel('Power/Frequency (dB/Hz)')
            #     plt.grid(True)
            #     # Add annotation similar to MATLAB example
            #     plt.annotate('Fault frequency region of interest', xy=(1100, -80), xytext=(1500, -60), # Adjust coordinates as needed
            #                  arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8))
            #     plt.xlim(0, fs_vib_read / 2)
            #     plt.show()
            # else:
            #      print("Could not plot Vib_acpi spectrum (data missing or empty).")


            # # Plot Envelope of Vib_acpi_env
            # vib_acpi_env_signal = member_data.get('Vib_acpi_env')
            # if vib_acpi_env_signal is not None and not vib_acpi_env_signal.empty:
            #     plt.figure(figsize=(10, 5))
            #     # Plot using the time index directly
            #     plt.plot(vib_acpi_env_signal.index.total_seconds(), vib_acpi_env_signal.values)
            #     plt.title(f'Envelope of Bandpassed Vib_acpi (Member: {os.path.basename(ens.files[0])})')
            #     plt.xlabel('Time (s)')
            #     plt.ylabel('Amplitude')
            #     plt.ylim(-0.5, 0.5) # Match MATLAB example axis limits
            #     plt.grid(True)
            #     plt.show()
            # else:
            #      print("Could not plot Vib_acpi_env (data missing or empty).")

        else:
            print("Failed to read the first member.")
