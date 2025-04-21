import os
import zipfile
import urllib.request

# Set the filename and URL for the reduced dataset
filename = "experimental_database_short.zip"
url = "https://ssd.mathworks.com/supportfiles/predmaint/broken-rotor-bar-fault-data/" + filename

# Download the file if it doesn't exist
if not os.path.exists(filename):
    print(f"Downloading {filename}...")
    urllib.request.urlretrieve(url, filename)
    print("Download complete!")
else:
    print(f"File {filename} already exists.")

# Extract the compressed data files
print("Extracting files...")
with zipfile.ZipFile(filename, 'r') as zip_ref:
    zip_ref.extractall(".")
print("Extraction complete!")
