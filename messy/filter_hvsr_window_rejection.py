import pandas as pd
import os
import logging

logging.basicConfig(level=logging.INFO)

def clean_csv(input_csv, output_txt):
    try:
        # Step 1: Extract the header line (last commented line)
        with open(input_csv) as f:
            header_line = ""
            for line in f:
                if line.startswith("#"):
                    header_line = line.strip("# ").strip()  # Remove `#` and any extra whitespace

        # Step 2: Define column names from the header
        columns = header_line.split(",")

        # Step 3: Read the CSV, skipping commented lines, and set the header
        df = pd.read_csv(input_csv, comment='#', names=columns, skiprows=1)

        # Step 4: Select specific columns
        df = df[['frequency (Hz)', 'mean curve (lognormal)', 'mean curve std (lognormal)']]
        
        # Step 5: Save the cleaned data to a new text file without header and with whitespace delimiter
        df.to_csv(output_txt, index=False, header=False, sep=" ")
        logging.info(f"Successfully cleaned the CSV and saved it as {output_txt}.")
    
    except Exception as e:
        logging.error(f"An error occurred: {e}")

def process_all_csv_in_directory(directory):
    # Loop over all files in the given directory
    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            input_csv = os.path.join(directory, filename)
            output_txt = os.path.join(outdirectory, f"{os.path.splitext(filename)[0]}.txt")
            clean_csv(input_csv, output_txt)

# Example usage
directory = 'output/output_window_rejection_v13Nov'  # Replace with the path to your directory containing CSV files
outdirectory = 'output/output_window_rejection_v13Nov_filtered'
os.makedirs(outdirectory, exist_ok=True)  # Create the directory if it doesn't exist
process_all_csv_in_directory(directory)
