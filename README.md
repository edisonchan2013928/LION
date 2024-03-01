# LION
This is the implementation of our research paper "LION: Fast and High-Resolution Network Kernel Density Visualization", which is currently accepted in VLDB 2024. Based on our theoretical analysis in this paper, this method can further reduce the time complexity for generating high-resolution network kernel density visualization (NKDV) compared with the state-of-the-art ADA method [a].

[a] Tsz Nam Chan, Zhe Li, Leong Hou U, Jianliang Xu, Reynold Cheng. Fast Augmentation Algorithms for Network Kernel Density Visualization. Proceedings of the VLDB Endowment (PVLDB), 2021.

# Our Code
The code is written in C++, which only depends on the C++ Standard Library (i.e., no additional library is used). Therefore, readers can test our code using the Cygwin software (for Windows OS) and the Terminal (for Linux OS), without any additional software/hardware. We have provided the shell script file "compile_run.sh" (in the "Code" folder) for showing how to compile and run our code. In particular, we have provided detailed descriptions of all parameters and an example for calling our code in "compile_run.sh".

# Datasets
Due to space limitations in this Github link, we do not upload the datasets. However, since all datasets are open to public, readers can simply download these datasets from the source webpages (cf. the "Ref." column in Table 3). To call our code properly, readers need to process these datasets based on our input file structure, which is shown in the file "input_file_structure.txt".

# Output results
After readers call our code, they will obtain an output text file that contains output results. This output text file is based on our output file structure, which is shown in the "output_file_structure.txt".
