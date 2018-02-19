for i in `seq 0 20`;
do
	deltaT=$(python -c "print($i/10000.0+1.0)")
	./hmc --anharmonic -l 1 -f 4 -d 0.02 -N 50 -T $deltaT 
done