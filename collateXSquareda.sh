# If the data files that contain the collated data do not already exists create them.
rm LatticeSpacingsExpectationXSquared.dat

touch LatticeSpacingsExpectationXSquared.dat


# Iterate through each set of results.
for results in "$@"*
do
	# Get the Tempering parameter from the input file and append it to each collated data file.s
	awk '/^Lattice-Spacing: /{printf $NF}' $results/input.txt >> LatticeSpacingsExpectationXSquared.dat
	
	# Put an empty column in each file.
	printf " " >> LatticeSpacingsExpectationXSquared.dat
 
	# Get the value of <X^2> and its error from the results file.
	awk '/^(X\^2): /{printf $(NF-2)}' $results/results.txt >> LatticeSpacingsExpectationXSquared.dat

	# Put an empty column in each file.
	printf " " >> LatticeSpacingsExpectationXSquared.dat
	awk '/^(X\^2): /{print $(NF)}' $results/results.txt >> LatticeSpacingsExpectationXSquared.dat

done

# Remove duplicate lines from files.
sort -u LatticeSpacingsExpectationXSquared.dat -o LatticeSpacingsExpectationXSquared.dat
