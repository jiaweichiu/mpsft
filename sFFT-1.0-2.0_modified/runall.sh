# Graph types 1 2 3.
# gGraphType=1: Fix S, increase N
# gGraphType=2: Fix N, increase S
# gGraphType=3: Fix S,N, increase SNR.

# ./generate_graphs 20   10000 0 1 1 > results/G1_sfft1.txt
# ./generate_graphs 20   10000 0 1 2 > results/G1_sfft2.txt

./generate_graphs 3   4000 0 1 1 results/G1_sff1.out > results/G1_sfft1.txt
./generate_graphs 3   4000 0 1 2 results/G1_sff2.out > results/G1_sfft2.txt

#./generate_graphs 50   100    0 2 1 > data_$FF/G2_sfft1.txt &
#./generate_graphs 50   100    0 2 2 > data_$FF/G2_sfft2.txt &

./generate_graphs 10   100    0 2 1 results/G2_sff1.out > results/G2_sfft1.txt &
./generate_graphs 10   100    0 2 2 results/G2_sff2.out > results/G2_sfft2.txt &

#./generate_graphs 100  1     1 3 1 > data_$FF/G3_sfft1.txt &
#./generate_graphs 100  1     1 3 2 > data_$FF/G3_sfft2.txt &


# gOuterReps = atoi(argv[1]);
# gInnerReps = atoi(argv[2]);
# gComputeError = atoi(argv[3]);
# gGraphType = atoi(argv[4]);
# gVer = atoi(argv[5]);

# std::ofstream fout(argv[6]);