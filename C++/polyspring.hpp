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


// apply lambdd on each block buffer of size bufsize, keep corresponding offset into contiguous vector
// https://stackoverflow.com/questions/66224510/simplest-way-to-pass-a-lambda-as-a-method-parameter-in-c17:
// - use type-erasing wrapper:
// void blockwise (int numbuffers, int bufsize[], CoordT *buffers[], std::function< void (int bufsize, CoordT *buffer, int offset)> func)
// - use template to deduce type:
template<typename F>
void blockwise (int numbuffers, int bufsizes[], CoordT *buffers[], F&& func)
{
  int offset = 0;
  
  for (int i = 0; i < numbuffers; i++)
  {
    func(bufsizes[i], buffers[i], offset);
    offset += bufsizes[i];
  }
}

// 
void copy_with_strides (int num, CoordT *src, int srcstride, CoordT *dest, int deststride)
{
}


class Point
{
// original descr coords
        self.scaled_og_x = x // unused
        self.scaled_og_y = y
	self.scaled_x = x // unused
        self.scaled_y = y
// scale to 0..1	    
        normalized_x = (x - bounds[0]) / (bounds[1] - bounds[0])
        normalized_y = (y - bounds[2]) / (bounds[3] - bounds[2])
	self.og_x = normalized_x # original position
        self.og_y = normalized_y
        self.x = normalized_x # current position
        self.y = normalized_y
	self.shap = ShPoint(normalized_x, normalized_y) // unused
	self.uni_x = normalized_x # position after uniformisation // unused
        self.uni_y = normalized_y
        self.prev_x = normalized_x # position at previous triangulation
        self.prev_y = normalized_y
        self.push_x = 0.0 # amount of pushing for next step
        self.push_y = 0.0
	self.near = [] # references to indices of points connected by triangles, typically 6
};

class Edge
{   // built from triangulation, needs to be updated after point moved
    float length;
    float dist_x, dist_y;
    //float mid_x, mid_y;// only needed for h_dist
    float h_dist; // evaluated at mid_x/y
    PointIndex a, b;
};

class Polyspring
{
    //std::vector<CoordT> x_, y_;	// normalised coords 0..1, created by pre-uniformisation
    std::vector<CoordT> points_;	// normalised interleaved point coords 0..1, created by pre-uniformisation
    //std::vector<CoordT> push_x_, push_y_; // displacement vector
    std::vector<CoordT> push_;		// interleaved displacement vector
    std::vector<CoordT> tripoints_;	// interleaved(!) array of coordinates for delaunay triangulation, need to keep for use in distFromOrigin
    std::vector<Edge>   edges_;		// list of edges

    const double	edge_correction = 2;	// factor taking into account that we store unique edges, while python code visits every ege twice via the near list of point a containing b and reciprocally b containing a
    double		l0_uni_;	// spring rest length
};


// copy points from buffers into vector, do rescaling
void Polyspring::set_points (int numtotal, int numbuffers, int bufsizes[], CoordT *buffers[], int bufwidth, int colx, int coly)
{
  points_.resize(2 * numtotal);
  blockwise(numbuffers, bufsizes, buffers,
	    [&](int bufsize, CoordT *buffer, int offset) -> void
	    {
	      copy_with_strides(bufsize, buffer + colx, bufwidth, points_.data(),     2);
	      copy_with_strides(bufsize, buffer + coly, bufwidth, points_.data() + 1, 2);
	    });

  // get min/max x/y

  // set default region 0..1 square, signed distance function, inner box
  double area = 1;

  // scale points to fit into region

  // compute the spring rest length
  // self.l0_uni = np.sqrt(2 / (np.sqrt(3) * len(self.points) / self.region.area))
  l0_uni_ = sqrt(2 / (sqrt(3) * numtotal / area));
  

}

void Polyspring::build_edges (triangles)
{
    // check if edge a, b or b, a already exists
}

void Polyspring::update_edges()
{
    for (auto e: edges_)
    {
	e.dist_x = point_[2 * e.b]     - point_[2 * e.a];
	e.dist_y = point_[2 * e.b + 1] - point_[2 * e.a + 1];
	e.length = std::sqrt(e.dist_x * e.dist_x + e.dist_y * e.dist_y);
	CoordT midx = point_[2 * e.a] + 0.5 * e.dist_x;
	CoordT midy = point_[2 * e.a + 1] + 0.5 * e.dist_y;
	h_dist = get_h(midx, midy);
    }
}

// called after updating edges
double Polyspring::get_scaling_factor(self)
{
  double target_area = 0;
  int	 npair = 0;
  /*
  for point in self.points:
    for near in point.near:
      npair += 1
      midx ,midy = point.midTo(near)
      target_area += 1 / self.h_dist(midx, midy)**2 
  return self.l0_uni * np.sqrt(npair / target_area)
  */
  
  for (auto e: edges_) // do over precalculated edges midpoints and density h
    target_area += 1. / (e.h_dist * e.h_dist) * edge_correction;

  return l0_uni * sqrt(edges_.size() / target_area);
}

    def midTo(self, point):
        midx = (self.x + point.x)/2
        midy = (self.y + point.y)/2
        return midx, midy

    def distTo(self, point):
        return np.sqrt((self.x-point.x)**2 + (self.y-point.y)**2)

    def moveDist(self):
        return np.sqrt(self.push_x ** 2 + self.push_y ** 2)


    def repulsiveForce(self, f, point):
        angle = np.arctan2(self.y - point.y, self.x - point.x)
        self.push_x += f * np.cos(angle)
        self.push_y += f * np.sin(angle)
        self.shap = ShPoint(
            self.x + self.push_x,
            self.y + self.push_y) # update the shapely point now for outside observation

    def moveDist(self):
        return np.sqrt(self.push_x ** 2 + self.push_y ** 2)

    def distFromOrigin(self):
        return np.sqrt((self.x-self.prev_x)**2 + (self.y-self.prev_y)**2)
        
    def moveTo(self, coords):
        nextx = coords[0]
        nexty = coords[1]
        self.push_x = nextx - self.x
        self.push_y = nexty - self.y

    def update(self, bounds):
        self.x += self.push_x
        self.y += self.push_y
        self.scaled_x = self.x * (bounds[1] - bounds[0]) + bounds[0]
        self.scaled_y = self.y * (bounds[3] - bounds[2]) + bounds[2]
        self.push_x = 0.0
        self.push_y = 0.0
        self.shap = ShPoint(self.x, self.y)
        
main		    
        # pre-uniformization

        # simulation parameters
        dt = 0.2        # simulation step
        tri_tol = 0.1   # displacement threshold (relative to l0_uni) for retriangularisation
        int_pres = 1.2  # spring pressure
        k = 1           # spring stiffness (supposing mass = 1)

	h_dist: function to get target distance for point

	// main loop
        while not exit:
            exit = True

	// update triangulation if necessary
            if update_tri:
                tri_count += 1
                self.delaunayTriangulation()
		// construct edges list
                update_tri = False
		    

	// compute rest length scaling factor
	// loop over edges
            hscale = self.getScalingFactor()

	// sum repulsive actions for each point
	// loop over edges
            for point in self.points: 
		for near in point.near: // each edge is visited twice!!!!!
                    midX ,midY = point.midTo(near)
                    f = k * (int_pres * hscale / self.h_dist(midX, midY) - point.distTo(near))
                    if f > 0:
                        near.repulsiveForce(dt * f, point) // update push vector with force from near point

	// second loop after all forces computation
	// loop over points, 
            for point in self.points:
		// update point positions
			    points_ += push_;

                	// check stop condition if inside region, else move it back inside
                if point.shap.within(self.region): # shap point is already pushed
                    if exit and point.moveDist() / self.l0_uni > stop_tol: 
                        exit = False
                else:
                    point.moveTo(nearest_points(self.region, point.shap)[0].coords[0]) ?????
			
	// update point positions
	// loop over points, or edges (and set push to 0)
                point.update(self.bounds)
			// if live update, rescale to descr. coordinates and write to output buffers
	// update edges after points have moved
	edges_update();
			/* NOT: for (auto e: edges_)
			{
			    e.dist_x_ += push_[2 * e.b]     - push_[2 * e.a];
			    e.dist_y_ += push_[2 * e.b + 1] - push_[2 * e.a + 1];
			    }
*/

	// reset push to 0

	// loop over points, 
		// check if triangulation needs to be updated: moved from prev_x/y more than tri_tol thresh
                if not update_tri and point.distFromOrigin() / self.l0_uni > tri_tol:
                    update_tri = True
