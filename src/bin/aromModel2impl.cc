/*
 * Converts a models feature weight vector from string notation into a C++
 * implementation form that can be included within the library.
 *
 *  Created on: 29.09.2011
 *      Author: mmann
 */

#include <iostream>

#include <sgd/vectors.h>

int main(int argc, char **argv) {

	if (argc > 1) {
		std::cout <<"\n model2impl : reads a model in string encoding from "
				"STDIN and writes a C++ implementation to STDOUT.\n"
				<<std::endl;
		return 0;
	}

	using namespace sgd;

	SVector vec;

	std::cin >>vec;

	const SVector::Rep * rep = vec.data();
	std::cout <<"\n\n\n";
	std::cout <<"\tint wdatasize = " << rep->npairs <<";"<<std::endl;
	std::cout <<"\tsgd::SVector::Rep wdata[] = {\n ";
	for (int i=0; i<rep->npairs; ++i ) {
		std::cout <<(i>0?", {":" {")<<rep->pairs[i].i <<","<<rep->pairs[i].v<<"}";
		if ((i+1) % 5 == 0) {
			std::cout <<"\n";
		}
	}
	std::cout <<" };"<<std::endl;
	std::cout <<"\n\n\n";


	return 0;
}


