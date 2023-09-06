/* -*-mode:c; c-basic-offset: 2-*- */
/* compile and run with:
   c++ -std=c++17 -g -fsanitize=address -fno-omit-frame-pointer -Idelaunator-cpp/include/ test-polyspring.cpp && ./a.out
 */

#include "stdio.h"
#include "polyspring.hpp"

int main (int argc, char *argv[])
{
  // create test data
  int   width   = 2;
  int   bufsize = 10;
  float buffer[bufsize * width];
  float *buf[1] = { buffer };

  for (int i = 0; i < bufsize; i++)
  { // fill lower third
    buffer[i * 2    ] = (float) i / bufsize;
    buffer[i * 2 + 1] = (float) (i % 3) / bufsize;
  }
  print_points("init", bufsize, buffer);

  Polyspring<float> poly;

  // set corpus, copies blocks (buffers) into points array
  poly.set_points(bufsize, 1, &bufsize, buf, width, 0, 1);
  print_points("set", bufsize, poly.points_.get_points_interleaved().data());

  while (poly.iterate()  &&  poly.get_count() < 100)
  {
    printf("iter %d, tri %d\n", poly.get_count(), poly.get_triangulation_count());
    print_points("", bufsize, poly.points_.get_points_interleaved().data());
  }
}
