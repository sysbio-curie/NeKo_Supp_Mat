import csv
import sys
import os

def convert_txt_to_csv(file_name):
    # Check if the file has a .txt extension
    if not file_name.endswith('.txt'):
        print("Error: The input file does not have a .txt extension.")
        return
    
    # Set the output file name with a .csv extension
    output_file = file_name.replace('.txt', '.csv')
    
    # Open the TXT file and the CSV file
    with open(file_name, 'r', newline='') as txt_file:
        with open(output_file, 'w', newline='') as csv_file:
            tsv_reader = csv.reader(txt_file, delimiter='\t')
            csv_writer = csv.writer(csv_file, delimiter=',')
            
            # Write all rows from the TXT file to the CSV file
            for row in tsv_reader:
                csv_writer.writerow(row)

    print(f"Converted {file_name} to {output_file} successfully.")

if __name__ == "__main__":
    # Ensure a file name is provided as an argument
    if len(sys.argv) != 2:
        print("Usage: python script.py file_name.txt")
    else:
        file_name = sys.argv[1]
        convert_txt_to_csv(file_name)
