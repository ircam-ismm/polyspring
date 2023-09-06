/* -*-mode:c; c-basic-offset: 2-*- */

#include "stdio.h"
#include "polyspring.hpp"

void print_points(const char * msg, int num, float *buf)
{
  if (msg && msg[0])
    printf("%s\n", msg);
  
  for (int i = 0; i < num; i++)
    printf("%6.3f %6.3f\n", buf[i * 2], buf[i * 2 + 1]);

  printf("\n");
}

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
  poly.set_points(bufsize, 1, &bufsize, (float **) &buffer, width, 0, 1);
  print_points("set", bufsize, poly.points_.get_points_interleaved().data());

  while (poly.iterate())
  {
    printf("iter %d, tri %d\n", poly.get_count(), poly.get_triangulation_count());
    print_points("", bufsize, poly.points_.get_points_interleaved().data());
  }
}
