import os
import sys
import random
import math
import time
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

def apply_quantum_seed():
    print("Initializing Quantum Circuit...")
    
    # Create a quantum circuit with 16 qubits and 16 classical bits
    qc = QuantumCircuit(16, 16)
    
    # Apply a Hadamard (H) gate to every qubit to put them in superposition
    for i in range(16):
        qc.h(i)
        
    # Measure the qubits, collapsing their quantum state into classical bits
    qc.measure(range(16), range(16))
    
    # Run the circuit
    simulator = AerSimulator()
    job = simulator.run(qc, shots=1)
    result = job.result()
    
    # Extract the resulting binary string and convert it to a decimal integer
    counts = result.get_counts()
    binary_string = list(counts.keys())[0]
    quantum_seed = int(binary_string, 2)
    
    print(f"Quantum collapse resulted in binary: {binary_string}")
    print(f"Success! Python is now seeded with quantum randomness: {quantum_seed}")
    
    # Seed the standard random module with our quantum number
    random.seed(quantum_seed)

def read_file_into_string(input_file, ord_range):
    the_file = open(input_file, 'r')
    current_char = the_file.read(1)
    file_string = ""
    length = len(ord_range)
    while current_char != "":
        i = 0
        while i < length:
            if ord(current_char) >= ord_range[i][0] and ord(current_char) <= ord_range[i][1]:
                file_string = file_string + current_char
                i = length
            else:
                i = i + 1
        current_char = the_file.read(1)
    the_file.close()
    return file_string


def remove_all_spaces(the_string):
    length = len(the_string)
    new_string = ""
    for i in range(length):
        if the_string[i] != " ":
            new_string = new_string + the_string[i]
    return new_string


def integerize(the_string):
    length = len(the_string)
    stripped_string = "0"
    for i in range(0, length):
        if ord(the_string[i]) >= 48 and ord(the_string[i]) <= 57:
            stripped_string = stripped_string + the_string[i]
    resulting_int = int(stripped_string)
    return resulting_int


def convert_to_list_of_int(the_string):
    list_of_integers = []
    location = 0
    finished = False
    while finished == False:
        found_comma = the_string.find(',', location)
        if found_comma == -1:
            finished = True
        else:
            list_of_integers.append(integerize(
                the_string[location:found_comma]))
            location = found_comma + 1
            if the_string[location:location + 5] == "NOTE=":
                finished = True
    return list_of_integers


def build_distance_matrix(num_cities, distances, city_format):
    dist_matrix = []
    i = 0
    if city_format == "full":
        for j in range(num_cities):
            row = []
            for k in range(0, num_cities):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    elif city_format == "upper_tri":
        for j in range(0, num_cities):
            row = []
            for k in range(j):
                row.append(0)
            for k in range(num_cities - j):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    else:
        for j in range(0, num_cities):
            row = []
            for k in range(j + 1):
                row.append(0)
            for k in range(0, num_cities - (j + 1)):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    if city_format == "upper_tri" or city_format == "strict_upper_tri":
        for i in range(0, num_cities):
            for j in range(0, num_cities):
                if i > j:
                    dist_matrix[i][j] = dist_matrix[j][i]
    return dist_matrix


input_file = "AISearchfile012.txt"
if len(sys.argv) > 1:
    input_file = sys.argv[1]
if len(sys.argv) > 2:
    max_time = int(sys.argv[2])

path_for_city_files = "city-files"
path_to_input_file = os.path.join(path_for_city_files, input_file)

if os.path.isfile(path_to_input_file):
    ord_range = [[32, 126]]
    file_string = read_file_into_string(path_to_input_file, ord_range)
    file_string = remove_all_spaces(file_string)
    print(f"Reading input file: {input_file}")
else:
    print(f"*** error: The city file {input_file} does not exist.")
    sys.exit()

location = file_string.find("SIZE=")
if location == -1:
    print(f"*** error: {input_file} is incorrectly formatted (Missing SIZE).")
    sys.exit()

comma = file_string.find(",", location)
if comma == -1:
    sys.exit()

num_cities_as_string = file_string[location + 5:comma]
num_cities = integerize(num_cities_as_string)

comma = comma + 1
stripped_file_string = file_string[comma:]
distances = convert_to_list_of_int(stripped_file_string)

counted_distances = len(distances)
if counted_distances == num_cities * num_cities:
    city_format = "full"
elif counted_distances == (num_cities * (num_cities + 1))/2:
    city_format = "upper_tri"
elif counted_distances == (num_cities * (num_cities - 1))/2:
    city_format = "strict_upper_tri"
else:
    print(f"*** error: Distance count mismatch in {input_file}.")
    sys.exit()

dist_matrix = build_distance_matrix(num_cities, distances, city_format)

start_time = time.time()

apply_quantum_seed()

def calculate_tour_len(tour):
    tour_length = 0
    for i in range(0, num_cities - 1):
        tour_length = tour_length + dist_matrix[tour[i]][tour[i + 1]]
    tour_length = tour_length + dist_matrix[tour[num_cities - 1]][tour[0]]
    return tour_length


# Generate random tour to start at
current_tour = list(range(num_cities))
random.shuffle(current_tour)
current_tour_len = calculate_tour_len(current_tour)
tour_length = current_tour_len
tour = list(current_tour)

# Calculate the range of all the edges in the graph
largest_edge = 0
smallest_edge = math.inf
for row in dist_matrix:
    for dist in row:
        largest_edge = max(dist, largest_edge)
        smallest_edge = min(dist, smallest_edge)
edge_range = largest_edge - smallest_edge

# Set the temperature parameters
min_temp = 0.000000000001
temp = 10 * edge_range / (-math.log(0.5))
decay_constant = 0.99999

while temp > min_temp:
    if time.time() - start_time > max_time:
        break
    # Select 2 random indexes for the current tour (which are not the same)
    i = random.randint(0, num_cities-1)
    j = (i + random.randint(1, num_cities-1)) % num_cities

    # Make sure that i is smaller than j
    if i > j:
        i, j = j, i

    # When i and j are on the edge, make sure i is the last index and j is the first index
    if (i == 0 and j == num_cities - 1):
        i, j = j, i

    # Find the cities affected by switching positions i and j
    city1 = current_tour[i-1]
    city2 = current_tour[i]
    city3 = current_tour[(i+1) % num_cities]
    city4 = current_tour[j-1]
    city5 = current_tour[j]
    city6 = current_tour[(j+1) % num_cities]

    # If the indexes are adjacent, you only have to remove and add 2 edges
    if (j == i + 1) or (j == 0 and i == num_cities - 1):
        removed_cost = dist_matrix[city1][city2] + dist_matrix[city5][city6]
        added_cost = dist_matrix[city1][city5] + dist_matrix[city2][city6]
    # Else add and remove all of the edges affected
    else:
        removed_cost = dist_matrix[city1][city2] + dist_matrix[city2][city3] + \
            dist_matrix[city4][city5] + dist_matrix[city5][city6]
        added_cost = dist_matrix[city1][city5] + dist_matrix[city5][city3] + \
            dist_matrix[city4][city2] + dist_matrix[city2][city6]

    length_change = added_cost - removed_cost

    # If the tour has shrunk, always move to this new tour
    if length_change < 0:
        current_tour_len += length_change
        current_tour[i], current_tour[j] = current_tour[j], current_tour[i]
        # If the tour is a new global best, update it
        if current_tour_len < tour_length:
            tour = list(current_tour)
            tour_length = current_tour_len

    # If the tour length has increased, decide with probability based on length change and temp
    else:
        probability = random.random()
        if probability < math.e ** (-length_change / temp):
            current_tour_len += length_change
            current_tour[i], current_tour[j] = current_tour[j], current_tour[i]

    # Decrease the temp
    temp = temp * decay_constant

print(f"Optimal tour length found: {tour_length}")
print(f"Best tour path: {tour}")


# Calculate time taken
end_time = time.time()
time_taken = end_time - start_time
print(f"Time taken: {time_taken}")

# Ensure the 'tours' directory exists
tours_dir = "quantum_tours"
if not os.path.exists(tours_dir):
    os.makedirs(tours_dir)

# Create a clean filename based on the input file
base_name = os.path.basename(input_file).split('.')[0]
output_file_name = os.path.join(tours_dir, f"{base_name}_result.txt")

should_save = False

# Check if a previous run exists and compare it
if os.path.exists(output_file_name):
    prev_tour_length = math.inf
    prev_time_taken = math.inf
    
    # Read the existing file to find previous stats
    with open(output_file_name, 'r') as f:
        for line in f:
            if line.startswith("Tour Length:"):
                prev_tour_length = int(line.split(":")[1].strip())
            elif line.startswith("Time Taken:"):
                prev_time_taken = float(line.split(":")[1].replace("seconds", "").strip())
                
    # Compare the new run against the old run
    if tour_length < prev_tour_length:
        print(f"New best tour found! Length {tour_length} beats previous {prev_tour_length}.")
        should_save = True
    elif tour_length == prev_tour_length and time_taken < prev_time_taken:
        print(f"Tied tour length of {tour_length}, but faster! {time_taken:.4f}s beats {prev_time_taken:.4f}s.")
        should_save = True
    else:
        print(f"Existing tour is better or equal. Not saving. (Previous: {prev_tour_length} len in {prev_time_taken:.4f}s)")
else:
    # No previous file exists, so we should save it
    print("No previous tour found. Saving new results.")
    should_save = True

# Write the data to the file if the conditions are met
if should_save:
    with open(output_file_name, 'w') as f:
        f.write(f"Input File: {input_file}\n")
        f.write(f"Time Taken: {time_taken:.4f} seconds\n")
        f.write(f"Tour Length: {tour_length}\n")
        f.write(f"Tour Path: {tour}\n")
    print(f"Results successfully saved to {output_file_name}")