#ifndef VERSION_HH
#define VERSION_HH

#include <iostream>

void
giveVersion()
{
	std::cout	<<"\n " <<PACKAGE_NAME <<" package version " <<PACKAGE_VERSION
				<<"\n"
				<<"\n " <<PACKAGE_URL
				<<"\n"
				<<std::endl;
}

#endif /*VERSION_HH*/
