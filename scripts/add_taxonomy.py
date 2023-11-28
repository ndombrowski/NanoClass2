from Bio import SeqIO
import argparse
import subprocess

def add_taxonomy(input_file, taxonomy_file, output_file):
    # Create a dictionary to store taxonomy information
    taxonomy_dict = {}

    # Read taxonomy file and store information in the dictionary
    with open(taxonomy_file, "r") as tax_file:
        for line in tax_file:
            parts = line.strip().split()
            seq_header = parts[0]
            taxonomy = " ".join(parts[1:])
            taxonomy_dict[seq_header] = taxonomy

    # Update the headers in the fasta file and write to a temporary file
    temp_output_file = "temp_output.fna"
    with open(temp_output_file, "w") as output:
        for record in SeqIO.parse(input_file, "fasta"):
            seq_header = record.id
            taxonomy = taxonomy_dict.get(seq_header, "No_taxonomy_info")
            new_header = f">{seq_header} {taxonomy}"
            output.write(f"{new_header}\n{record.seq}\n")

    # Use sed to remove "D_0__", "D_1__", etc. prefixes from taxonomy strings
    subprocess.run(["sed", "-i", 's/D_[0-9]__//g', temp_output_file])

    # Rename the temporary file to the final output file
    subprocess.run(["mv", temp_output_file, output_file])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add taxonomy information to fasta file headers.")
    parser.add_argument("-i", "--input", help="Input fasta file", required=True)
    parser.add_argument("-t", "--taxonomy", help="Taxonomy file", required=True)
    parser.add_argument("-o", "--output", help="Output fasta file with taxonomy", required=True)
    
    args = parser.parse_args()

    add_taxonomy(args.input, args.taxonomy, args.output)
