/* -*-mode:c; c-basic-offset: 2-*- */

/* required: 

delaunayTriangulation

- Geogram (INRIA)
https://github.com/BrunoLevy/geogram/wiki
https://brunolevy.github.io/geogram/classGEO_1_1Delaunay2d.html

points as interleaved double *


- Delaunator A really fast C++ library for Delaunay triangulation of 2D points.
  MIT License
  delaunator-cpp is a C++ port from https://github.com/mapbox/delaunator a JavaScript implementation of very fast 2D Delaunay algorithm.
  Probably the fastest C++ open source 2D Delaunay implementation

https://github.com/delfrrr/delaunator-cpp/blob/master/include/delaunator.hpp

single header
points as interleaved vector<double> (---> vector<pair<buffer, marker> ???)
returns triangles as vector<int> of point indices √


- bl4ckb0ne / delaunay-triangulation
GPL

https://github.com/bl4ckb0ne/delaunay-triangulation/tree/master

points as vector<vertextype>,
vertextype = vector2D<T>  -->  override with pair<buffer, index>
returns vector<TriangleType>,
TriangleType = vertextype *a, *b, *c

improvement from:
https://github.com/Neverland1026/Delaunay_2D_CPP/commit/30d11bc967f26b2dc9e9567089f6aebc6ad94a34

*/

/*
  force calculation in loop over points or connections????

  max number of triangles = 2 npoints - b - 2  (b = number of borders, min. 3)

  density map from jitter matrix --> mubu
*/


#include <cmath>
#include <vector>
#include <numeric>      // for iota
#include <algorithm>    // for sort


///////////////////////////////////////////////////////////////////////////////
// UTILITIES

// apply lambda on each block buffer of size bufsize, keep corresponding offset into contiguous vector
// https://stackoverflow.com/questions/66224510/simplest-way-to-pass-a-lambda-as-a-method-parameter-in-c17:
// - use type-erasing wrapper:
// void blockwise (int numbuffers, int bufsize[], CoordT *buffers[], std::function< void (int bufsize, CoordT *buffer, int offset)> func)
// - use template to deduce type:
template<typename CoordT, typename F>
void blockwise (int numbuffers, int bufsizes[], CoordT *buffers[], F&& func)
{
  int offset = 0;
  
  for (int i = 0; i < numbuffers; i++)
  {
    func(bufsizes[i], buffers[i], offset);
    offset += bufsizes[i];
  }
}

// vector function wrappers
template<typename CoordT>
void copy_with_strides (int num, CoordT *src, int srcstride, CoordT *dest, int deststride)
{
}

template<typename CoordT>
std::vector<CoordT> vector_add(std::vector<CoordT> &a, std::vector<CoordT> &b)
{
  return a; // a + b
}

template<typename CoordT>
std::vector<CoordT> vector_mul(std::vector<CoordT> &a, std::vector<CoordT> &b)
{
  return a; // a * b
}

template<typename CoordT>
void vector_norm(std::vector<CoordT> &norm, std::vector<CoordT> &vec)
{
  ; // norm = sqrt(vec_x ^ 2 + vec_y ^ 2);
}


///////////////////////////////////////////////////////////////////////////////
// Regions

// interface for region
template<typename CoordT>
class Region
{
public:
  virtual bool point_is_within (CoordT x, CoordT y) = 0;
  virtual void move_point_back (CoordT &x, CoordT &y) = 0;
};

template<typename CoordT>
class SquareRegion : public Region<CoordT>
{
public:
  virtual bool point_is_within (CoordT x, CoordT y) override
  {
    return x >= 0  &&  x <= 1  &&  y >= 0  &&  y <= 1;
  }
  
  virtual void move_point_back (CoordT &x, CoordT &y) override
  {
    if (x < 0)  x = 0;
    if (x > 1)  x = 1;
    if (y < 0)  y = 0;
    if (y > 1)  y = 1;
  }
}; // end class SquareRegion


// macros for accessing x/y in interleaved vectors
#define x(index) ((index) * 2)
#define y(index) ((index) * 2 + 1)


///////////////////////////////////////////////////////////////////////////////
template<typename CoordT>
struct Edges
{ // built from triangulation, needs to be updated after point moved
  int numedges_ = 0;
  std::vector<CoordT> length_;	// length of edge = norm of dist_
  std::vector<CoordT> dist_;	// interleaved x/y vector of edge length
  //float mid_x, mid_y;// only needed for h_dist calculation
  std::vector<CoordT> h_dist_;  // vector of target density evaluated at mid_x/y
  std::vector<int>    a_, b_;	// indices into Points arrays to end points
  const double	edge_correction_ = 2;	// factor taking into account that we store unique edges, while python code visits every ege twice via the near list of point a containing b and reciprocally b containing a
  CoordT (*get_h_)(CoordT x, CoordT y);

  Edges (CoordT (*hfunc)(CoordT x, CoordT y))
  : get_h_(hfunc) {}

  void init (int num)
  {
    numedges_ = num;
    //TDB: .clear()?
    length_.resize(num);
    dist_.resize(num * 2);
    h_dist_.resize(num);
    a_.resize(num);
    b_.resize(num);
  }

  // build edges list from interleaved triangle vertex index list
  void set (std::vector<int> &tri)
  {
    assert(tri.size() == 3 * numedges_);

    // check if edge a, b or b, a already exists
  }

  // update edges after points have moved
  void update (std::vector<CoordT> &points)
  {
/*  for (auto e: edges_)
    {
	e.dist_x = point_[2 * e.b]     - point_[2 * e.a];
	e.dist_y = point_[2 * e.b + 1] - point_[2 * e.a + 1];
	e.length = std::sqrt(e.dist_x * e.dist_x + e.dist_y * e.dist_y);
	CoordT midx = point_[2 * e.a]     + 0.5 * e.dist_x;
	CoordT midy = point_[2 * e.a + 1] + 0.5 * e.dist_y;
	h_dist = get_h(midx, midy);
    }
*/
    for (int i = 0; i < numedges_; i++)
    {
      dist_[x(i)] = points[x(b_[i])] - points[x(a_[i])];
      dist_[y(i)] = points[y(b_[i])] - points[y(a_[i])];
    }
    vector_norm(length_, dist_); // length = sqrt(dist_x ^ 2 + dist_y ^ 2);

    for (int i = 0; i < numedges_; i++)
    { // get middle point and lookup density function
      CoordT midx = points[x(a_[i])] + 0.5 * dist_[x(i)];
      CoordT midy = points[y(a_[i])] + 0.5 * dist_[y(i)];
      h_dist_[i]  = get_h_(midx, midy);
    }
  } // end Edges::update ()

  // called after setting/updating edges
  double scaling_factor ()
  { /* for point in self.points:
	   for near in point.near:
		npair += 1
		midx ,midy = point.midTo(near)
		target_area += 1 / self.h_dist(midx, midy)**2 
      return self.l0_uni * np.sqrt(npair / target_area)
 */
     // loop over edges (do over precalculated density h at edges midpoints)
    //double target_area = sum(1. / (h_dist_ * h_dist_) * edge_correction_);
    double target_area = 1; // sum(1. / vector_mul(h_dist_, h_dist_) * edge_correction_);
    
    return sqrt(numedges_ / target_area);
  }

  // apply spring repulsive force f to edge index i's end points' push vectors
  void apply_force (int i, double f, std::vector<CoordT> &push)
  {
/*  def repulsiveForce(self, f, point):
        angle = np.arctan2(self.y - point.y, self.x - point.x)
        self.push_x += f * np.cos(angle)
        self.push_y += f * np.sin(angle)
        self.shap = ShPoint(
            self.x + self.push_x,
            self.y + self.push_y) # update the shapely point now for outside observation
*/
    double angle = atan2((double) dist_[y(i)], (double) dist_[x(i)]);
    push[x(a_[i])] += f * cos(angle);
    push[y(a_[i])] += f * sin(angle);
    push[x(b_[i])] -= f * cos(angle);
    push[y(b_[i])] -= f * sin(angle);
  }
}; // end struct Edges


///////////////////////////////////////////////////////////////////////////////
template<typename CoordT>
struct Triangulation
{
  std::vector<CoordT> tripoints_;	// interleaved(!) array of coordinates for delaunay triangulation, need to keep for use in dist_since_triangulation
  std::vector<int> vertices_;		// interleaved(!) array of triangle vertex indices (into tripoints array)

  int count_ = 0;

  void init (int num)
  {
    tripoints_.reserve(num * 2); // don't init elements, will be overwritten by triangulate()
    count_ = 0;
  }

  void triangulate (std::vector<CoordT> &points)
  {
    tripoints_.assign(points.begin(), points.end());

    // call delaunay triangulation
    
    
    count_++;
  }

  std::vector<int> &get_vertices() { return vertices_; }
};


///////////////////////////////////////////////////////////////////////////////
template<typename CoordT>
struct Points
{
  int numpoints_ = 0;
  std::vector<CoordT> points_;	// normalised interleaved point coords 0..1, created by pre-uniformisation
  std::vector<CoordT> push_;	// interleaved displacement vector x/y

  std::vector<CoordT> &get_points_interleaved() { return points_; }

  void init (int num)
  {
    numpoints_ = num;
    points_.resize(num * 2);
    push_.resize(num * 2);
  }

  void set (int numtotal, int numbuffers, int bufsizes[], CoordT *buffers[], int bufwidth, int colx, int coly)
  {
    numpoints_ = numtotal;
    points_.resize(2 * numtotal);
    blockwise(numbuffers, bufsizes, buffers,
	      [&](int bufsize, CoordT *buffer, int offset) -> void
	      {
		copy_with_strides(bufsize, buffer + colx, bufwidth, points_.data(),     2);
		copy_with_strides(bufsize, buffer + coly, bufwidth, points_.data() + 1, 2);
	      });
  }
  
  void pre_uniformize ()
  { // replace points (in descriptor coordinates) by normalised sort index
    // create rank array of size numframes
    std::vector<size_t> indices(numpoints_);

    // for each column
    for (int colind = 0; colind < 2; colind++)
    {
      // fill rank array with 0..numframes
      std::iota(begin(indices), end(indices), static_cast<size_t>(0));

      // sort rank array with indirection to data
      std::sort(
        begin(indices), end(indices),
        [&](size_t a, size_t b)
        {
          // ascending order, such that small elements get small indices
          return points_[a * 2 + colind] < points_[b * 2 + colind];
        }
      );
      
      // scale to 0..1 and copy order of elem back to points
      for (size_t i = 0; i < numpoints_; i++)
      {
        size_t order = indices[i];
	points_[order * 2 + colind] = (float) i / (numpoints_ - 1);
      }
    } // end for colind
  } // end Points::pre_uniformize ()
  
  void update ()
  {
    vector_add(points_, push_);    // points_ += push_;
  }

  void end_iteration ()
  { // set push vector to zero
    std::fill(push_.begin(), push_.end(), 0);
  }
  
  bool within_region (int i, Region<CoordT> *region)
  {
    return region->point_is_within(points_[x(i)], points_[y(i)]);
  }

  void move_point_back (int i, Region<CoordT> *region)
  {
    region->move_point_back(points_[x(i)], points_[y(i)]); // changes points coords
  }

  CoordT norm (CoordT x1, CoordT y1)
  {
    return std::sqrt(x1 * x1 + y1 * y1);
  }

  CoordT dist_moved (int i)
  {
    return norm(push_[x(i)], push_[y(i)]);
  }

  CoordT dist_since_triangulation (int i, Triangulation<CoordT> &tri)
  {
    return norm(points_[x(i)] - tri.tripoints_[x(i)], points_[y(i)] - tri.tripoints_[y(i)]);
  }
}; // end struct Points


///////////////////////////////////////////////////////////////////////////////
template<typename CoordT>
class Polyspring
{
public:
  // simulation parameters
  double dt_	   = 0.2;   // simulation step
  double tri_tol_  = 0.1;   // displacement threshold (relative to l0_uni) for retriangularisation
  double int_pres_ = 1.2;   // spring pressure
  double k_	   = 1;     // spring stiffness (supposing mass = 1)
  double stop_tol_ = 0.001;
      
  //h_dist: function to get target distance for point
  Region<CoordT>	*region_ = NULL;

  Points<CoordT>	points_;	// points
  Triangulation<CoordT> triangulation_;	// wrapper around delaunay triangulation
  Edges<CoordT>		edges_;		// keeps list of edges
  double		l0_uni_;	// spring rest length
  int			count_ = 0;
  bool			update_tri_ = false;

public:
  Polyspring ()
  : edges_(get_h) {}

  void set_region (std::string name);
  void set_points (int numtotal, int numbuffers, int bufsizes[], CoordT *buffers[], int bufwidth, int colx, int coly);   // copy points from buffers into vector, do rescaling and pre-uniformisation
  static CoordT get_h (CoordT x, CoordT y) { return 1; } // TODO: evaluate target density function at point
  bool iterate ();
  int get_count() { return count_; };
  int get_triangulation_count() { return triangulation_.count_; };
};


// copy points from buffers into vector, do rescaling
template<typename CoordT>
void Polyspring<CoordT>::set_points (int numtotal, int numbuffers, int bufsizes[], CoordT *buffers[], int bufwidth, int colx, int coly)
{
  points_.set(numtotal, numbuffers, bufsizes, buffers, bufwidth, colx, coly);
  // get min/max x/y

  // set default region 0..1 square, signed distance function, inner box
  double area = 1;

  // pre-uniformization, replaces points by normalised sort index
  // and scale points to fit into region
  points_.pre_uniformize();
    
  // compute the spring rest length
  // self.l0_uni = np.sqrt(2 / (np.sqrt(3) * len(self.points) / self.region.area))
  l0_uni_ = sqrt(2 / (sqrt(3) * numtotal / area));
  
  count_ = 0;
  update_tri_ = true;
}


// main loop
template<typename CoordT>
bool Polyspring<CoordT>::iterate ()
{
  // update triangulation if necessary
/* 	  if update_tri:
                tri_count += 1
                self.delaunayTriangulation()
                update_tri = False
*/
  if (update_tri_)
  {
    triangulation_.triangulate(points_.get_points_interleaved());
    edges_.set(triangulation_.get_vertices());	// construct edges list
    edges_.update(points_.get_points_interleaved());
    update_tri_ = false;
  }
		  
  // compute rest length scaling factor
  double hscale = l0_uni_ * edges_.scaling_factor();

  // sum repulsive actions for each point
  /* for point in self.points: 
  	for near in point.near: // each edge is visited twice!!!!!
              midX ,midY = point.midTo(near)
              f = k * (int_pres * hscale / self.h_dist(midX, midY) - point.distTo(near))
              if f > 0:
                  near.repulsiveForce(dt * f, point) // update push vector with force from near point
  */

  // calculate spring forces 
  for (int i = 0; i < edges_.numedges_; i++) // loop over over precalculated edge length and density h at edges midpoints
  {
    double f = k_ * (int_pres_ * hscale / edges_.h_dist_[i] - edges_.length_[i]); // TODO: vectorise this loop, second loop with apply_force
    if (f > 0)
      edges_.apply_force(i, dt_ * f, points_.push_); // update edge's end points' push vectors with force from spring
  }
		  
  // move point positions, but leave push
  points_.update();

  // second loop after all forces computation
  /* for point in self.points:        // check stop condition if inside region, else move it back inside
                if point.shap.within(self.region): # shap point is already pushed
                    if exit and point.moveDist() / self.l0_uni > stop_tol: 
                        exit = False
                else:
                    point.moveTo(nearest_points(self.region, point.shap)[0].coords[0]) ?????
  */
  bool keep_going = true;
  for (int i = 0; /* keep_going  && */ i < points_.numpoints_; i++)
  { // check stop condition if inside region, else move it back inside
    if (points_.within_region(i, region_))
    {
      if (points_.dist_moved(i) / l0_uni_ <= stop_tol_)
      {
	keep_going = false;
	break;
      }
    }
    else
      points_.move_point_back(i, region_);
  }

  // if live update, rescale to descr. coordinates and write to output buffers
  //points_.get(bounds, outbuffers);
  
  // update edges after points have moved
  edges_.update(points_.points_);
  /* NOT: for (auto e: edges_)     {
		e.dist_x += push_[2 * e.b]     - push_[2 * e.a];
		e.dist_y += push_[2 * e.b + 1] - push_[2 * e.a + 1];
	   }  */

  // loop over points, 
  // check if triangulation needs to be updated: moved from prev_x/y more than tri_tol thresh
  for (int i = 0; i < points_.numpoints_; i++) // todo: vectorise, mostly goes through all points
    if (points_.dist_since_triangulation(i, triangulation_) / l0_uni_ > tri_tol_)
    {
      update_tri_ = true;
      break;
    }

  // set push to 0
  points_.end_iteration();

  return keep_going;
} // end Polyspring::iterate ()
