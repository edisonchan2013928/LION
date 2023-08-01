#Compile the program
g++ -c shortest_path.cpp -w -o shortest_path.o -std=c++11
g++ -c KAF.cpp -w -o KAF.o -std=c++11
g++ -c lion.cpp -w -o lion.o -std=c++11
g++ -c alg_NKDV.cpp -w -o alg_NKDV.o -std=c++11
g++ main.cpp -O3 -o main shortest_path.o KAF.o lion.o alg_NKDV.o
exit #Remove it if you want to run our code.

#These are the parameters for calling our code.
#our_model.network_fileName = argv[1]; //The input network file name
#our_model.out_NKDV_fileName = argv[2]; //The output visualization file name
#our_model.method = atoi(argv[3]); //method = 1: RQS, method = 2: SPS, method = 3: ADA, and method = 6: LION 
#our_model.lixel_reg_length = atoi(argv[4]); //The length of the lixel size (\ell in Figure 1).
#our_model.k_type = atoi(argv[5]); //k_type = 1: Triangular kernel, k_type = 2: Epanechnikov kernel, and k_type = 3: Quartic kernel
#our_model.bandwidth = atof(argv[6]); //The bandwidth parameter (in meters)

network_fileName="./Datasets/Detroit"
out_NKDV_fileName="./Results/Detroit_NKDV_output"
method=6
lixel_reg_length=10
k_type=2
bandwidth=1000

./main $network_fileName $out_NKDV_fileName $method $lixel_reg_length $k_type $bandwidth