#include "alg_NKDV.h"

int main(int argc, char**argv)
{
	alg_NKDV algorithm;
	algorithm.load_network(argv);
	//algorithm.NKDV_compute(argc, argv);
	algorithm.NKDV_compute_text_file(argc, argv);
	algorithm.clear_basic_memory();
}