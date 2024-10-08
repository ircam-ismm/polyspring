/* -*-mode:c; c-basic-offset: 2-*- */
/* compile and run with:
   c++ -std=c++17 -g -fsanitize=address -fno-omit-frame-pointer -I /sw/include -I .. -I ../delaunator-cpp/include/ -L/sw/lib -llo test-polyspring.cpp && ./a.out
 */

#include "stdio.h"
#include "unistd.h"	// sleep
#include "time.h"	// clock
#include "lo/lo.h"	// osc
#include <fstream>
#include <iterator>
#include <vector>
#include "polyspring.hpp"

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

void osc_send_buffer (int bufind, int numrows, float *buffer, float *orig, int width, int xcol, int ycol)
{
  int stat = 0;
  stat += lo_send(addr, "/buffer_index", "i", bufind);
  stat += lo_send(addr, "/clear", "");

  for (int i = 0; i < numrows; i++)
  {
    float x = buffer[i * 2];
    float y = buffer[i * 2 + 1];
    float origx = x;
    float origy = y;
    if (orig)
    {
      origx = orig[i * width + xcol];
      origy = orig[i * width + ycol];
    }
    stat += lo_send(addr, "/append", "iffff", i, origx, origy, x, y);
  }
  //stat += lo_send(addr, "/bounds", "ffff", ....);
  stat += lo_send(addr, "/done_init", "i", 1);
  stat += lo_send(addr, "/update", "i", 1);

  if (stat < 0)
  {
    printf("osc_send_buffer status %d, error %d '%s'\n", stat, lo_address_errno(addr), lo_address_errstr(addr));
  }
}

void osc_send_tri (std::vector<size_t> &vertices)
{
  int stat = 0;
  int numv = vertices.size();
  lo_message msg = lo_message_new();

    for (int i = 0; i < numv; i++)
  {
    lo_message_add_int32(msg, vertices[i]);
  }

  stat += lo_send_message(addr, "/tri", msg);

  if (stat < 0)
  {
    printf("osc_send_tri status %d, error %d '%s'\n", stat, lo_address_errno(addr), lo_address_errstr(addr));
  }

  lo_message_free(msg);
}

int data_gen_linear (int bufsize, float **bufs, int &width, int &xcol, int &ycol)
{
  float *buffer = (float *) malloc(bufsize * 2 * sizeof(float)); // memleak, but we don't care
  int wrap = bufsize <= 3 ? 2 : 3;

  for (int i = 0; i < bufsize; i++)
  { // fill lower third
    buffer[i * 2    ] = (float) i / bufsize;
    buffer[i * 2 + 1] = (float) (i % wrap) / bufsize;
  }

  bufs[0] = buffer;
  width = 2;
  xcol = 0;
  ycol = 1;
  return bufsize;
}

int data_load (const char *filename, float **bufs, int &width, int &xcol, int &ycol)
{
  std::ifstream file(filename);
  if (file.fail())
  {
    printf("ERROR: can't open %s.\n", filename);
    return -1;
  }
  
  static std::vector<float> vec(std::istream_iterator<float>(file), {});

  // columns in txt files : timetags | original x | original y | final x | final y
  bufs[0] = vec.data();
  width = 5;
  xcol = 1;
  ycol = 2;
  return vec.size() / 5;
}

int main (int argc, char *argv[])
{
  if (!osc_open(NULL /*"127.0.0.1"*/, "8012"))
    return 1;

  // create test data
  float *buf[1] = { NULL };
  int width;	// layout of loaded buffer
  int xcol;
  int ycol;

  //int   bufsize = data_gen_linear(3 , buf, width, xcol, ycol);
  int   bufsize = data_load("/Users/schwarz/src/polyspring/test-data/truth-1.txt", buf, width, xcol, ycol);
  if (bufsize <= 0)  exit(1);

  float *buffer = buf[0];  
  print_points("init", bufsize, buffer, width, xcol, ycol);

  Polyspring<float> poly;

  // set corpus, copies blocks (buffers) into points array
  poly.set_points(bufsize, 1, &bufsize, buf, width, xcol, ycol);
  osc_send_buffer(1, bufsize, poly.points_.get_points_interleaved(true).data(), buffer, width, xcol, ycol);
  print_points("set", bufsize, poly.points_.get_points_interleaved(true).data());

  bool keepgoing;
  const int numiter = 100;
  do
  {
    clock_t start_iter = clock();
    keepgoing = poly.iterate()  &&  poly.get_count() < numiter;
    clock_t stop_iter = clock();
    float dur = (stop_iter - start_iter) / (float) CLOCKS_PER_SEC * 1000.;
  
    printf("iter %d  tri %d  go %d\n", poly.get_count(), poly.get_triangulation_count(), keepgoing);
    osc_send_buffer(1, bufsize, poly.points_.get_points_interleaved(true).data(), buffer, width, xcol, ycol);
    //print_points("", bufsize, poly.points_.get_points_interleaved(true).data());
    osc_send_tri(poly.triangulation_.get_vertices());
    printf("iter %d  tri %d  go %d took %f ms\n", poly.get_count(), poly.get_triangulation_count(), keepgoing, dur);

    usleep(std::max(0., (100 - dur) * 1000.));
  } while (keepgoing);

  osc_close();
}
