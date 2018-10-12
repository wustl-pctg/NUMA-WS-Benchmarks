CILK_NWORKERS=1 taskset -c 0 ./hull  -r 5   ../geometryData/data/2DonSphere_10M > 1.log
CILK_NWORKERS=4 taskset -c 0-4 ./hull  -r 5   ../geometryData/data/2DonSphere_10M > 4.log
CILK_NWORKERS=8 taskset -c 0-7 ./hull  -r 5   ../geometryData/data/2DonSphere_10M > 8.log
CILK_NWORKERS=16 taskset -c 0-15 ./hull  -r 5   ../geometryData/data/2DonSphere_10M > 16.log
CILK_NWORKERS=24 taskset -c 0-23 ./hull  -r 5  ../geometryData/data/2DonSphere_10M > 24.log
CILK_NWORKERS=32 taskset -c 0-31 ./hull  -r 5   ../geometryData/data/2DonSphere_10M > 32.log
