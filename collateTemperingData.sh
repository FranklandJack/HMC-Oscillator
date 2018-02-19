# If the data files that contain the collated data do not already exists create them.
rm TemperingVsX.dat


touch TemperingVsX.dat


# Iterate through each set of results.
for results in "$@"*
do
	# Get the Tempering parameter from the input file and append it to each collated data file.s
	awk '/^Tempering-Parameter: /{printf $NF}' $results/input.txt >> TemperingVsX.dat
	
	# Put an empty column in each file.
	printf " " >> TemperingVsX.dat
 
	# Get the value of <X> and its error from the results file.
	awk '/^X: /{printf $(NF-2)}' $results/results.txt >> TemperingVsX.dat
	# Put an empty column in each file.
	printf " " >> TemperingVsX.dat
	awk '/^X: /{print $(NF)}' $results/results.txt >> TemperingVsX.dat

done

# Remove duplicate lines from files.
sort -u TemperingVsX.dat -o TemperingVsX.dat
