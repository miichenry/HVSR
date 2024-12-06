import os

# Directory containing the files and subdirectories
input_directory = "../output/output_DFA_filtered_v11Nov/" #"/srv/beegfs/scratch/shares/cdff/DPM/Processing/STACK_phase_only_rma"

# Output text file
output_file = "DFA_filtered_v11Nov.txt"

# Ensure the directory for the output file exists
output_dir = os.path.dirname(output_file)
if output_dir:  # Only attempt to create directory if it's non-empty
    os.makedirs(output_dir, exist_ok=True)

# Initialize an empty list to hold all file paths
file_list = []

# Walk through the directory tree
for root, dirs, files in os.walk(input_directory):
    for file in files:
        #if file.endswith('.h5'):  # Filter to include only .h5 files
        # Construct the full file path
        file_path = os.path.join(root, file)
        file_list.append(file_path)

# Write the file paths to the output file
try:
    with open(output_file, 'w') as f:
        for file_path in file_list:
            f.write(file_path + '\n')
    print(f"File list created: {output_file}")
except Exception as e:
    print(f"An error occurred while writing to the file: {e}")

