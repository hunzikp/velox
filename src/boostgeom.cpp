// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]

#include <vector>
#include <Rcpp.h>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/index/rtree.hpp>


using namespace Rcpp;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


typedef bg::model::d2::point_xy<double> Point;
typedef bg::model::multi_point<Point> MultiPoint;
typedef bg::model::polygon<Point> Polygon;
typedef bg::model::ring<Point> Ring;
typedef bg::model::multi_polygon<Polygon> MultiPolygon;
typedef bg::model::box<Point> Box;
typedef bg::model::segment<Point> Segment;
typedef bg::model::linestring<Point> Line;
typedef bg::model::multi_linestring<Line> MultiLine;
typedef std::pair<Box, unsigned> BoxValue;
typedef std::pair<Point, unsigned> PointValue;
typedef std::pair<Polygon, unsigned> PolygonValue;
typedef std::pair<MultiPolygon, unsigned> MultiPolygonValue;
typedef std::pair<MultiLine, unsigned> MultiLineValue;


// GEOMETRY COLLECTION CLASSES & FUNCTIONS

Box makeLookupGeometry(Polygon &p) {
  return bg::return_envelope<Box>(p);
}

Box makeLookupGeometry(Line &l) {
  return bg::return_envelope<Box>(l);
}

Point makeLookupGeometry(Point &p) {
  return p;
}


template <typename LookupContainer, typename LookupGeometry,
          typename MultiGeometry, typename Geometry>
class MultiGeometryCollection : public std::vector<MultiGeometry>{
public:

  MultiGeometryCollection(std::vector<MultiGeometry> &x) {
    for (unsigned i = 0; i < x.size(); ++i) {

      // Store MultiGeometry
      this->push_back(x[i]);

      // Store LookupContainer
      for (unsigned j = 0; j < x[i].size(); ++j) {
        LookupGeometry lg = makeLookupGeometry(x[i][j]);
        LookupContainer luc(lg, std::pair<unsigned, unsigned>(i, j));
        lookupVector.push_back(luc);
      }
    }

    // Index the mofo
    bgi::rtree<LookupContainer, bgi::rstar<16> > rtree_new(lookupVector);
    rtree = rtree_new;
  }

  template <typename TargetCollection, typename TargetMultiGeometry, typename TargetGeometry,
            typename TargetLookupGeometry>
  List intersects(TargetCollection targetCollection) {

    // Iterate over target multigeometries
    std::vector<std::vector<unsigned> > out;
    for (unsigned i = 0; i < targetCollection.size(); ++i) {
      TargetMultiGeometry targetMultiGeometry = targetCollection[i];

      std::vector<unsigned> hitVec;
      for (unsigned j = 0; j < targetMultiGeometry.size(); ++j) {
        TargetGeometry targetGeometry = targetMultiGeometry[j];
        TargetLookupGeometry targetLookupGeometry = makeLookupGeometry(targetGeometry);

        // Lookup using rtree and lookup geometries
        std::vector<LookupContainer> result_s;
        rtree.query(bgi::intersects(targetLookupGeometry), std::back_inserter(result_s));

        // Refine lookup and fill hit vector
        for (unsigned k = 0; k < result_s.size(); ++k) {
          unsigned mGeometryIndex = result_s[k].second.first;

          bool foundBefore = std::find(hitVec.begin(), hitVec.end(), mGeometryIndex + 1) != hitVec.end();
          if (foundBefore) {
            continue;
          }

          unsigned geometryIndex = result_s[k].second.second;
          Geometry sourceGeometry = this->at(mGeometryIndex)[geometryIndex];
          bool intersects = bg::intersects(sourceGeometry, targetGeometry);
          if (intersects) {
            hitVec.push_back(mGeometryIndex + 1);
          }
        }
      }
      out.push_back(hitVec);
    }

    return wrap(out);
  }

  MultiGeometryCollection subset(IntegerVector index) {
    std::vector<MultiGeometry> mgVec;
    for (int i = 0; i < index.size(); ++i) {
      unsigned idx = index[i] - 1;
      if ((idx < this->size()) && (idx >= 0)) {
        MultiGeometry mg = this->at(idx);
        mgVec.push_back(mg);
      }
    }
    MultiGeometryCollection mgc(mgVec);
    return mgc;
  }

  unsigned length() {
    return this->size();
  }

private:
  std::vector<LookupContainer> lookupVector;
  bgi::rtree<LookupContainer, bgi::rstar<16> > rtree;
};

typedef std::pair<Box, std::pair<unsigned, unsigned> > BoxMultiValue;
typedef std::pair<Point, std::pair<unsigned, unsigned> > PointMultiValue;

typedef MultiGeometryCollection<BoxMultiValue, Box, MultiPolygon, Polygon> MultiPolygonCollection;
typedef MultiGeometryCollection<BoxMultiValue, Box, MultiLine, Line> MultiLineCollection;
typedef MultiGeometryCollection<PointMultiValue, Point, MultiPoint, Point> MultiPointCollection;



// GRID COLLECTION CLASSES AND FUNCTIONS

std::vector<Point> convertToPoints(Box &b) {

  Point p1(b.min_corner().x(), b.min_corner().y());
  Point p2(b.min_corner().x(), b.max_corner().y());
  Point p3(b.max_corner().x(), b.min_corner().y());
  Point p4(b.max_corner().x(), b.max_corner().y());

  std::vector<Point> outVec;
  outVec.push_back(p1);
  outVec.push_back(p2);
  outVec.push_back(p3);
  outVec.push_back(p4);

  return outVec;
}

std::vector<Point> convertToPoints(Polygon &p) {

  std::vector<Point> outVec;
  for (unsigned i = 0; i < p.outer().size() - 1; ++i) {
    outVec.push_back(p.outer().at(i));
  }
  for (unsigned i = 0; i < p.inners().size(); ++i) {
    for (unsigned j = 0; j < p.inners().at(i).size() - 1; ++j) {
      outVec.push_back(p.inners().at(i).at(j));
    }
  }

  return outVec;
}

struct ChargedRing {
public:
  ChargedRing () {};
  ChargedRing (Ring &ring) {

    for (unsigned i = 0; i < ring.size() - 1; ++i) {
      polyX.push_back(ring[i].x());
      polyY.push_back(ring[i].y());
    }

    polyCorners = ring.size() - 1;

    constant.resize(polyCorners);
    multiple.resize(polyCorners);

    // Prepare polygon data
    int j = polyCorners-1;
    for(unsigned i = 0; i < polyCorners; i++) {

      // Polygon structure
      if(polyY[j]==polyY[i]) {
        constant[i]=polyX[i];
        multiple[i]=0;
      }
      else {
        constant[i] = polyX[i]-(polyY[i]*polyX[j])/(polyY[j]-polyY[i])+(polyY[i]*polyX[i])/(polyY[j]-polyY[i]);
        multiple[i] = (polyX[j]-polyX[i])/(polyY[j]-polyY[i]);
      }
      j = i;
    }
  }

  unsigned polyCorners;
  std::vector<double> constant;
  std::vector<double> multiple;
  std::vector<double> polyX, polyY;
};

struct ChargedPolygon {
public:
  ChargedPolygon(Polygon &polygon) {
    ChargedRing outerR(polygon.outer());
    outerRing = outerR;
    for (unsigned i = 0; i < polygon.inners().size(); ++i) {
      ChargedRing innerR(polygon.inners().at(i));
      innerRings.push_back(innerR);
    }
  }

  ChargedPolygon(Box &box) {
    Ring outerR;
    bg::convert<Box, Ring>(box, outerR);
    ChargedRing outerCR(outerR);
    outerRing = outerCR;
  }

  ChargedRing outerRing;
  std::vector<ChargedRing> innerRings;
};

bool fastIntersects(Point &p, ChargedRing &cr) {
  double tx = p.x();
  double ty = p.y();

  bool oddNodes = false;
  int j = cr.polyCorners-1;

  for (unsigned i = 0; i < cr.polyCorners; i++) {
    if (((cr.polyY[i]< ty && cr.polyY[j]>=ty)
           ||   (cr.polyY[j]< ty && cr.polyY[i]>=ty))) {
      oddNodes^=(ty*cr.multiple[i]+cr.constant[i]<tx);
    }
    j=i;
  }

  return oddNodes;
}

bool fastIntersects(Point &p, ChargedPolygon &pr) {

  bool intersects = fastIntersects(p, pr.outerRing);
  if (intersects) {
    bool inner = false;
    for (unsigned i = 0; i < pr.innerRings.size(); ++i) {
      inner = fastIntersects(p, pr.innerRings[i]);
      if (inner) {
        intersects = false;
        break;
      }
    }
  }

  return intersects;
}

std::vector<unsigned> fastIntersects(std::vector<PointValue> pointValues, Polygon &p) {

  std::vector<unsigned> out;
  ChargedPolygon cp(p);
  for (unsigned i = 0; i < pointValues.size(); ++i) {
    bool intersects = fastIntersects(pointValues[i].first, cp);
    if (intersects) {
      out.push_back(pointValues[i].second);
    }
  }

  return out;
}

bool fastCovers(ChargedRing &chargedRing, std::vector<Point> pVec) {

  bool covers = true;
  for (unsigned i = 0; i < pVec.size(); ++i) {
    bool intersects = fastIntersects(pVec[i], chargedRing);
    if (!intersects) {
      covers = false;
      break;
    }
  }

  return covers;
}

template <typename Politype>
bool fastCovers(ChargedRing &chargedRing, Politype &pt) {
  std::vector<Point> pVec = convertToPoints(pt);
  bool covers = fastCovers(chargedRing, pVec);
  return covers;
}

template <typename Linetype>
std::vector<Segment> segmentize(Linetype &r) {
  std::vector<Segment> sVec;
  for (unsigned i = 0; i < r.size()-1; ++i) {
    Segment s(r[i], r[i+1]);
    sVec.push_back(s);
  }
  return sVec;
}

std::vector<Segment> segmentize(Polygon &p) {

  Ring outerRing = p.outer();
  std::vector<Segment> sVec = segmentize(outerRing);

  for (unsigned i = 0; i < p.inners().size(); ++i) {
    std::vector<Segment> thisSVec = segmentize(p.inners().at(i));
    sVec.insert(sVec.end(), thisSVec.begin(), thisSVec.end());
  }

  return sVec;
}

std::vector<Segment> segmentize(MultiPolygon &mp) {

  std::vector<Segment> sVec;

  for (unsigned i = 0; i < mp.size(); ++i) {
    std::vector<Segment> thisSVec = segmentize(mp[i]);
    sVec.insert(sVec.end(), thisSVec.begin(), thisSVec.end());
  }

  return sVec;
}

std::vector<unsigned> fastIntersects(std::vector<BoxValue> &boxes, Polygon &polygon) {

  std::vector<unsigned> hitValues;

  // Index segmentized polygon
  std::vector<Segment> pSegments = segmentize(polygon);
  bgi::rtree<Segment, bgi::rstar<16> > rt(pSegments);

  // Create ChargedPolygon
  ChargedPolygon chargedPoly(polygon);

  // Create point vector of polygon
  std::vector<Point> polyPoints = convertToPoints(polygon);

  for (unsigned i = 0; i < boxes.size(); i++) {

    // Check whether boxes intersect with polygon border
    std::vector<Segment> result_s;
    rt.query(bgi::intersects(boxes[i].first), std::back_inserter(result_s));
    if (result_s.size() > 0) {
      hitValues.push_back(boxes[i].second);
      continue;
    }

    // Check whether outer ring of polygon covers box
    bool outerCovered = fastCovers(chargedPoly.outerRing, boxes[i].first);
    if (outerCovered) {
      // Check whether any inner ring of polygon completely covers box
      bool innercovered = false;
      for (unsigned j = 0; j < chargedPoly.innerRings.size(); ++j) {
        ChargedRing innerChargedRing = chargedPoly.innerRings[j];
        innercovered = fastCovers(innerChargedRing, boxes[i].first);
        if (innercovered) {  // One hit is enough info
          break;
        }
      }
      if (!innercovered) {
        hitValues.push_back(boxes[i].second);
        continue;
      }
    }
  }

  return hitValues;
}

template <typename GeomValue>  // For PointValue and BoxValue
std::vector<unsigned> fastIntersects(std::vector<GeomValue> &geomValues, Line &line) {

  std::vector<unsigned> hitValues;

  // Index segmentized line
  std::vector<Segment> pSegments = segmentize(line);
  bgi::rtree<Segment, bgi::rstar<16> > rt(pSegments);

  for (unsigned i = 0; i < geomValues.size(); i++) {

    // Check whether geometries intersect with line segments
    std::vector<Segment> result_s;
    rt.query(bgi::intersects(geomValues[i].first), std::back_inserter(result_s));
    if (result_s.size() > 0) {
      hitValues.push_back(geomValues[i].second);
    }
  }

  return hitValues;
}




// Template for grids of type pair<point, unsigned>, and pair<box, unsigned>
template <typename Gridtype>
class GridCollection : public std::vector<Gridtype> {
public:
  GridCollection(std::vector<Gridtype> &x,
                 unsigned inrow, unsigned incol,
                 double ixres, double iyres,
                 double ixmin, double iymax) {
    for (unsigned i = 0; i < x.size(); ++i) {
      this->push_back(x[i]);
    }
    nrow = inrow;
    ncol = incol;
    xres = ixres;
    yres = iyres;
    xmin = ixmin;
    ymax = iymax;
  }

  GridCollection subset(unsigned startrow, unsigned endrow,
                        unsigned startcol, unsigned endcol) {

    std::vector<Gridtype> gVec;
    for (unsigned j = startcol; j < endcol; ++j) {
      for (unsigned i = startrow; i < endrow; ++i) {
        gVec.push_back(this->at(j*nrow + i));
      }
    }
    unsigned new_nrow = endrow-startrow;
    unsigned new_ncol = endcol-startcol;
    double new_xmin = xmin + new_ncol*xres;
    double new_ymax = ymax - new_nrow*yres;
    GridCollection newGC(gVec,
                         new_nrow, new_ncol,
                         xres, yres,
                         new_xmin, new_ymax);

    return newGC;
  }

  GridCollection crop(double crop_xmin, double crop_xmax,
                      double crop_ymin, double crop_ymax) {

    double xmax = xmin + xres*ncol;
    double ymin = ymax - yres*nrow;

    if (crop_xmin > xmax || crop_xmax < xmin || crop_ymin > ymax || crop_ymax < ymin) {
      std::vector<Gridtype> gVec;
      GridCollection newGC(gVec,
                           0, 0,
                           xres, yres,
                           xmin, ymax);
      return newGC;
    }

    unsigned startcol = 0;
    if (crop_xmin > xmin) {
      startcol = floor((crop_xmin - xmin)/xres);
    }
    unsigned endcol = ncol;
    if (crop_xmax < xmax) {
      endcol = ceil((crop_xmax-xmin)/xres);
    }

    unsigned startrow = 0;
    if (crop_ymax < ymax) {
      startrow = floor((ymax - crop_ymax)/yres);
    }
    unsigned endrow = nrow;
    if (crop_ymin > ymin) {
      endrow = ceil((ymax - crop_ymin)/yres);
    }

    GridCollection newGC = this->subset(startrow, endrow, startcol, endcol);
    return newGC;
  }

  GridCollection crop(Box b) {
    double xmax = b.max_corner().x();
    double ymax = b.max_corner().y();
    double xmin = b.min_corner().x();
    double ymin = b.min_corner().y();
    return this->crop(xmin, xmax, ymin, ymax);
  }

  IntegerVector getValues() {
    std::vector<unsigned> valVec;
    for (unsigned i = 0; i < this->size(); ++i) {
      valVec.push_back((this->at(i)).second);
    }
    return wrap(valVec);
  }

  template <typename Collection, typename MultiGeometry, typename Geometry>
  List intersects(Collection &mgc) {

    std::vector<std::vector<unsigned> > out;

    for (unsigned i = 0; i < mgc.size(); ++i) {
      MultiGeometry* mgPtr = &(mgc[i]);
      std::vector<unsigned> indexVec;

      // iterate over multi
      for (unsigned j = 0; j < (*mgPtr).size(); ++j) {

        Geometry* gPtr = &(mgPtr->at(j));

        // crop GridCollection at this geometries's box
        Box b = bg::return_envelope<Box>(*gPtr);
        GridCollection cropGrid = this->crop(b);

        std::vector<unsigned> hitVec = fastIntersects(cropGrid, *gPtr);
        indexVec.insert(indexVec.end(), hitVec.begin(), hitVec.end());
      }

      // Remove duplicates
      sort( indexVec.begin(), indexVec.end() );
      indexVec.erase( unique( indexVec.begin(), indexVec.end() ), indexVec.end() );

      out.push_back(indexVec);
    }

    return wrap(out);
  }

private:
  unsigned nrow;
  unsigned ncol;
  double xres;
  double yres;
  double xmin;
  double ymax;
};

typedef GridCollection<std::pair<Point, unsigned> > PointGrid;
typedef GridCollection<std::pair<Box, unsigned> > BoxGrid;




// R CONSTRUCTORS AND DECONSTRUCTORS

class BoostFactory {
public:
  BoostFactory() {};

  Box makeBox(NumericVector x, NumericVector res) {
    // Construct from minimum corner and xy resolution (both vectors of length 2)
    Box b;
    b.min_corner().set<0>(x[0]);
    b.min_corner().set<1>(x[1]);
    b.max_corner().set<0>(x[0]+res[0]);
    b.max_corner().set<1>(x[1]+res[1]);
    return b;
  }

  PointGrid makePointGrid(NumericVector origin, IntegerVector dim, NumericVector res) {

    std::vector<std::pair<Point, unsigned> > pointValueVec;
    unsigned counter = 1;
    for (int j = 0; j < dim[1]; ++j) {
      for (int i = 0; i < dim[0]; ++i) {
        double x = origin[0] + (double)j*res[0] + res[0]/2;
        double y = origin[1] - (double)i*res[1] - res[1]/2;
        PointValue pointValue(Point(x, y), counter);
        pointValueVec.push_back(pointValue);
        counter++;
      }
    }

    PointGrid pointGrid(pointValueVec,
                        dim[0], dim[1],
                        res[0], res[1],
                        origin[0], origin[1]);
    return pointGrid;
  }

  BoxGrid makeBoxGrid(NumericVector origin, IntegerVector dim, NumericVector res) {

    std::vector<std::pair<Box, unsigned> > boxValueVec;
    unsigned counter = 1 ;
    for (int j = 0; j < dim[1]; ++j) {
      for (int i = 0; i < dim[0]; ++i) {
        NumericVector minc = NumericVector::create(origin[0] + (double)j*res[0],
                                                   origin[1] - (double)(i+1)*res[1]);
        Box b = makeBox(minc, res);
        BoxValue boxValue(b, counter);
        boxValueVec.push_back(boxValue);
        counter++;
      }
    }

    BoxGrid boxGrid(boxValueVec,
                    dim[0], dim[1],
                    res[0], res[1],
                    origin[0], origin[1]);

    return boxGrid;
  }

  Polygon makePolygon(List x) {
    Polygon poly;

    // Construct polygon from list of matrices (i.e. sf 'POLYGON' format)
    int nInnerRings = x.length() - 1;

    // Add outer ring
    NumericMatrix outerRingMat = x[0];
    for (int j = 0; j < outerRingMat.nrow(); ++j) {
      bg::append(poly.outer(), Point(outerRingMat(j,0), outerRingMat(j,1)));
    }

    // Add inner rings
    if (nInnerRings > 0) {
      poly.inners().resize(nInnerRings);
    }
    for (int j = 0; j < nInnerRings; ++j) {
      NumericMatrix innerRingMat = x[j+1];
      for (int k = 0; k < innerRingMat.nrow(); ++k) {
        bg::append(poly.inners()[j], Point(innerRingMat(k,0), innerRingMat(k,1)));
      }
    }

    return poly;
  }

  MultiPolygon makeMultiPolygon(List x) {
    // Construct multipolygon from list of list of matrices (i.e. sf 'MULTIPOLYGON' format)
    MultiPolygon mpoly;
    int nPolygons = x.length();
    for (int i = 0; i < nPolygons; ++i) {
      List polyList = x[i];
      Polygon poly = makePolygon(polyList);
      mpoly.push_back(poly);
    }

    return mpoly;
  }

  MultiPolygonCollection makeMultiPolygonCollection(List x) {

    std::vector<MultiPolygon> mpVec;

    StringVector classTypes = x.attr("class");
    std::string classType = as<std::string>(classTypes[0]);

    if (!classType.compare("sfc_MULTIPOLYGON")) {
      for (int i = 0; i < x.length(); ++i) {
        List mpolyList = x[i];
        MultiPolygon mpoly = makeMultiPolygon(mpolyList);
        mpVec.push_back(mpoly);
      }
    } else if (!classType.compare("sfc_POLYGON")) {
      for (int i = 0; i < x.length(); ++i) {
        List polyList = x[i];
        Polygon poly = makePolygon(polyList);
        MultiPolygon mpoly;
        mpoly.push_back(poly);
        mpVec.push_back(mpoly);
      }
    }

    MultiPolygonCollection mpc(mpVec);
    return mpc;
  }

  Line makeLine(NumericMatrix x) {
    Line ln;
    for (int i = 0; i < x.nrow(); ++i) {
      bg::append(ln, Point(x(i,0), x(i,1)));
    }
    return ln;
  }

  MultiLine makeMultiLine(List x) {
    MultiLine ml;
    for (int i = 0; i < x.length(); ++i) {
      NumericMatrix mat = x[i];
      Line ln = makeLine(mat);
      ml.push_back(ln);
    }
    return ml;
  }

  MultiLineCollection makeMultiLineCollection(List x) {
    std::vector<MultiLine> mlVec;

    StringVector classTypes = x.attr("class");
    std::string classType = as<std::string>(classTypes[0]);

    if (!classType.compare("sfc_MULTILINESTRING")) {
      for (int i = 0; i < x.length(); ++i) {
        MultiLine ml = makeMultiLine(x[i]);
        mlVec.push_back(ml);
      }
    } else if (!classType.compare("sfc_LINESTRING")) {
      for (int i = 0; i < x.length(); ++i) {
        Line ln = makeLine(x[i]);
        MultiLine ml;
        ml.push_back(ln);
        mlVec.push_back(ml);
      }
    }

    MultiLineCollection mlc(mlVec);

    return mlc;
  }

  Point makePoint(NumericVector x) {
    Point pnt(x[0], x[1]);
    return pnt;
  }

  MultiPoint makeMultiPoint(NumericMatrix x) {
    MultiPoint mp;
    for (int i = 0; i < x.nrow(); ++i) {
      Point pnt(x(i,0), x(i,1));
      mp.push_back(pnt);
    }
    return mp;
  }

  MultiPointCollection makeMultiPointCollection(List x) {
    std::vector<MultiPoint> mpVec;

    StringVector classTypes = x.attr("class");
    std::string classType = as<std::string>(classTypes[0]);

    if (!classType.compare("sfc_MULTIPOINT")) {
      for (int i = 0; i < x.length(); ++i) {
        MultiPoint mp = makeMultiPoint(x[i]);
        mpVec.push_back(mp);
      }
    } else if (!classType.compare("sfc_POINT")) {
      for (int i = 0; i < x.length(); ++i) {
        Point p = makePoint(x[i]);
        MultiPoint mp;
        mp.push_back(p);
        mpVec.push_back(mp);
      }
    }

    MultiPointCollection mpc(mpVec);
    return mpc;
  }

  NumericMatrix makeMatrix(Ring &ring) {
    NumericMatrix mat(ring.size(), 2);
    for (unsigned i = 0; i < ring.size(); ++i) {
      mat(i,0) = ring[i].x();
      mat(i,1) = ring[i].y();
    }
    return mat;
  }

  List makeMultiPolygonList(MultiPolygonCollection &mpc) {

    List out(mpc.size());
    for (unsigned i = 0; i < mpc.size(); ++i) {
      MultiPolygon mp = mpc[i];
      List mpList(mp.size());
      for (unsigned j = 0; j < mp.size(); ++j) {
        Polygon p = mp[j];
        List pList(1 + p.inners().size());
        pList[0] = makeMatrix(p.outer());
        for (unsigned k = 0; k < p.inners().size(); ++k) {
          pList[k+1] = makeMatrix(p.inners()[k]);
        }
        mpList[j] = pList;
      }
      out[i] = mpList;
    }

    return out;
  }

  List makeMultiLineList(MultiLineCollection &mlc) {
    List out(mlc.size());

    for (unsigned i = 0; i < mlc.size(); ++i) {
      MultiLine ml = mlc[i];
      List mlList(ml.size());
      for (unsigned j = 0; j < ml.size(); ++j) {
        Line ln = ml[j];
        NumericMatrix lineMat(ln.size(), 2);
        for (unsigned k = 0; k < ln.size(); ++k) {
          lineMat(k, 0) = ln[k].x();
          lineMat(k, 1) = ln[k].y();
        }
        mlList[j] = lineMat;
      }
      out[i] = mlList;
    }

    return out;
  }

  List makeMultiPointList(MultiPointCollection &mpc) {

    List out(mpc.size());
    for (unsigned i = 0; i < mpc.size(); ++i) {
      MultiPoint mp = mpc[i];
      NumericMatrix mpMat(mp.size(), 2);
      for (unsigned j = 0; j < mp.size(); ++j) {
        mpMat(j, 0) = mp[j].x();
        mpMat(j, 1) = mp[j].y();
      }
      out[i] = mpMat;
    }

    return out;
  }
};


RCPP_EXPOSED_CLASS(BoostFactory)
RCPP_EXPOSED_CLASS_NODECL(MultiPolygonCollection)
RCPP_EXPOSED_CLASS_NODECL(MultiLineCollection)
RCPP_EXPOSED_CLASS_NODECL(MultiPointCollection)
RCPP_EXPOSED_CLASS_NODECL(PointGrid)
RCPP_EXPOSED_CLASS_NODECL(BoxGrid)
RCPP_MODULE(BOOSTGEOM) {
  Rcpp::class_<MultiPolygonCollection>("MultiPolygonCollection")
    .method(  "intersectsMultiPolygon",
              &MultiPolygonCollection::intersects<MultiPolygonCollection, MultiPolygon, Polygon, Box>, "intersects")
    .method(  "intersectsMultiLine",
               &MultiPolygonCollection::intersects<MultiLineCollection, MultiLine, Line, Box>, "intersects")
    .method(  "intersectsMultiPoint",
               &MultiPolygonCollection::intersects<MultiPointCollection, MultiPoint, Point, Point>, "intersects")
    .method(  "subset", &MultiPolygonCollection::subset, "subset")
    .method(  "length", &MultiPolygonCollection::length, "length")
  ;

  Rcpp::class_<MultiLineCollection>("MultiLineCollection")
    .method(  "intersectsMultiLine",
              &MultiLineCollection::intersects<MultiLineCollection, MultiLine, Line, Box>, "intersects")
    .method(  "intersectsMultiPolygon",
              &MultiLineCollection::intersects<MultiPolygonCollection, MultiPolygon, Polygon, Box>, "intersects")
    .method(  "intersectsMultiPoint",
              &MultiLineCollection::intersects<MultiPointCollection, MultiPoint, Point, Point>, "intersects")
    .method(  "subset", &MultiLineCollection::subset, "subset")
    .method(  "length", &MultiLineCollection::length, "length")
  ;

  Rcpp::class_<MultiPointCollection>("MultiPointCollection")
    .method(  "intersectsMultiPoint",
              &MultiPointCollection::intersects<MultiPointCollection, MultiPoint, Point, Point>, "intersects")
  .method(  "intersectsMultiLine",
            &MultiPointCollection::intersects<MultiLineCollection, MultiLine, Line, Box>, "intersects")
  .method(  "intersectsMultiPolygon",
            &MultiPointCollection::intersects<MultiPolygonCollection, MultiPolygon, Polygon, Box>, "intersects")
  .method(  "subset", &MultiPointCollection::subset, "subset")
  .method(  "length", &MultiPointCollection::length, "length")
  ;

  Rcpp::class_<PointGrid>("PointGrid")
    .method("subset",     &PointGrid::subset,     "subset")
    .method("getValues",  &PointGrid::getValues,  "get index values")
    .method("intersectsMultiPolygon", &PointGrid::intersects<MultiPolygonCollection, MultiPolygon, Polygon>,  "intersects")
    .method("intersectsMultiLine", &PointGrid::intersects<MultiLineCollection, MultiLine, Line>,  "intersects")
  ;

  Rcpp::class_<BoxGrid>("BoxGrid")
    .method("subset",     &BoxGrid::subset,     "subset")
    .method("getValues",  &BoxGrid::getValues,  "get index values")
    .method("intersectsMultiPolygon", &BoxGrid::intersects<MultiPolygonCollection, MultiPolygon, Polygon>,  "intersects")
    .method("intersectsMultiLine", &BoxGrid::intersects<MultiLineCollection, MultiLine, Line>,  "intersects")
  ;

  Rcpp::class_<BoostFactory>("BoostFactory")
    .constructor("constructor")
    .method("makeMultiPolygonCollection", &BoostFactory::makeMultiPolygonCollection, "make mpc")
    .method("makeMultiLineCollection", &BoostFactory::makeMultiLineCollection, "make mlc")
    .method("makeMultiPointCollection", &BoostFactory::makeMultiPointCollection, "make mpntc")
    .method("makePointGrid", &BoostFactory::makePointGrid, "make pg")
    .method("makeBoxGrid", &BoostFactory::makeBoxGrid, "make bg")
    .method("makeMultiPolygonList", &BoostFactory::makeMultiPolygonList, "make mpl")
    .method("makeMultiLineList", &BoostFactory::makeMultiLineList, "make mll")
    .method("makeMultiPointList", &BoostFactory::makeMultiPointList, "make mpl")
  ;
}




