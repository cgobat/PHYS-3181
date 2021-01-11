#include <stdio.h>
#include "mt64.h"


int main(int argc, char** argv)
{
  int num = atoi(argv[1]);
  for(int j=0; j<num; ++j) printf("%e\n", genrand64_real2());
  return 1;
}
