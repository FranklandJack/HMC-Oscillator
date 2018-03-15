for i in `seq 0 10`;
do
	deltaT=$(python -c "print($i/1000.0+1.0)")
	./hmc --anharmonic -l 1 -f 4 -d 0.1 -N 10 -a 0.1 -L 1000 -b 100000 -T $deltaT 
done