/* -*-mode:c; c-basic-offset: 2-*- */
/* compile and run with:
   c++ -std=c++17 -g -fsanitize=address -fno-omit-frame-pointer -I /sw/include -I .. -I ../delaunator-cpp/include/ -L/sw/lib -llo test-polyspring.cpp && ./a.out
 */

#include "stdio.h"
#include "unistd.h"	// sleep
#include "polyspring.hpp"
#include "lo/lo.h"

lo_address addr;

bool osc_open (const char *host, const char *port)
{
  addr = lo_address_new(host, port);
  if (addr)
  {
    int err = lo_address_errno(addr);
    if (err != 0)
    {
      printf("osc_open error %d\n", err);
      return false;
    }
    else
    {
      lo_send(addr, "/done_import", "i", 1);
      return true;
    }
  }
  else
  {
    printf("osc_open error\n");
    return false; //addr != NULL; //TODO: other errors?
  }
}

void osc_close ()
{
  lo_address_free(addr);
}

void osc_send_buffer (int bufind, int numrows, float *buffer)
{
  int stat = 0;
  stat += lo_send(addr, "/buffer_index", "i", bufind);
  stat += lo_send(addr, "/clear", "");

  for (int i = 0; i < numrows; i++)
  {
    float x = buffer[i * 2];
    float y = buffer[i * 2 + 1];
    stat += lo_send(addr, "/append", "iffff", i, x, y, x, y);
  }
  //stat += lo_send(addr, "/bounds", "ffff", ....);
  stat += lo_send(addr, "/done_init", "i", 1);
  stat += lo_send(addr, "/update", "i", 1);

  if (stat < 0)
  {
    printf("osc_send_buffer status %d, error %d '%s'\n", stat, lo_address_errno(addr), lo_address_errstr(addr));
  }
}


int main (int argc, char *argv[])
{
  if (!osc_open(NULL /*"127.0.0.1"*/, "8012"))
    return 1;

// create test data
  int   width   = 2;
  int   bufsize = 25;
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
  osc_send_buffer(1, bufsize, poly.points_.get_points_interleaved().data());
  print_points("set", bufsize, poly.points_.get_points_interleaved().data());

  bool keepgoing;
  do
  {
    keepgoing = poly.iterate()  &&  poly.get_count() < 100;
    printf("iter %d  tri %d  go %d\n", poly.get_count(), poly.get_triangulation_count(), keepgoing);
    osc_send_buffer(1, bufsize, poly.points_.get_points_interleaved().data());
    print_points("", bufsize, poly.points_.get_points_interleaved().data());
    usleep(100 * 1000);
  } while (keepgoing);

  osc_close();
}
