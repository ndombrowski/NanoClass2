import argparse
import sys 

def filter_rows(input_file, output_file, query_coverage_threshold):
    try:
        if input_file == '-':
            data = sys.stdin.read()
        else:
            with open(input_file, 'r') as file:
                data = file.read()

        # Split the input data into rows
        rows = [line.split() for line in data.split('\n')]

        # Create a dictionary to store the rows based on column 1
        row_dict = {}
        for row in rows:
            if len(row) > 0:
                key = row[0]
                if key not in row_dict:
                    row_dict[key] = []
                row_dict[key].append(row)

        # Filter the rows based on conditions
        filtered_rows = []
        for key, values in row_dict.items():
            # Sort the values based on column 4, then column 3
            sorted_values = sorted(values, key=lambda x: (int(x[3]), int(x[2])), reverse=True)
            
            # Check if values in column 3 and 4 are different for the first and last elements
            if (sorted_values[0][2], sorted_values[0][3]) != (sorted_values[-1][2], sorted_values[-1][3]):
                # Add only the top value
                filtered_rows.append(sorted_values[0])
            else:
                # Add all values if both column 3 and 4 are the same
                filtered_rows.extend(sorted_values)

        #only keep rows above the query_coverage_threshold
        filtered_rows = [row for row in filtered_rows if float(row[9]) >= float(query_coverage_threshold)]
        
        # Write the result to the output file or stdout
        if output_file == '-':
            sys.stdout.write('\n'.join('\t'.join(row) for row in filtered_rows))
        else:
            with open(output_file, 'w') as file:
                for row in filtered_rows:
                    file.write('\t'.join(row) + '\n')

    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter rows based on specified conditions.")
    parser.add_argument("-i", "--input_file", help="Input file path", required=True)
    parser.add_argument("-o", "--output_file", help="Output file path", required=True)
    parser.add_argument("-c", "--query_coverage_threshold", help="Output file path", required=True)

    args = parser.parse_args()
    filter_rows(args.input_file, args.output_file, args.query_coverage_threshold)

