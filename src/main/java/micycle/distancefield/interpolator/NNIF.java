package micycle.distancefield.interpolator;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;

import org.tinfour.common.GeometricOperations;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IIncrementalTinNavigator;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.Thresholds;
import org.tinfour.common.Vertex;
import org.tinfour.interpolation.IVertexValuator;
import org.tinfour.interpolation.NaturalNeighborInterpolator;

/**
 * Indifferent to constrained edges etc.
 * 
 * @author MCarleton
 *
 */
public class NNIF extends NaturalNeighborInterpolator {

	IIncrementalTinNavigator navigator;

	final GeometricOperations geoOp;

	private Vertex v0, v1, v2;

	public NNIF(IIncrementalTin tin) {
		super(tin);

		Thresholds thresholds = tin.getThresholds();
		geoOp = new GeometricOperations(thresholds);

		navigator = tin.getNavigator();
	}

	@Override
	public double interpolate(double x, double y, IVertexValuator valuator) {

		List<IQuadEdge> eList = getBowyerWatsonEnvelope(x, y);
		int nEdge = eList.size();
		if (nEdge == 0) {
			// (x,y) is outside defined area
			return Double.NaN;
		} 
//		else if (nEdge == 1) {
//			// (x,y) is an exact match with the one edge in the list
//			IQuadEdge e = eList.get(0);
//			Vertex v = e.getA();
//			return v.getZ();
//		}

//		sumN++;
//		sumSides += eList.size();
		// The eList contains a series of edges definining the cavity
		// containing the polygon.
		double[] w = this.getBarycentricCoordinates(eList, x, y);
		if (w == null) {
			// the coordinate is on the perimeter, no Barycentric coordinates
			// are available.
			return Double.NaN;
		}
		double zSum = 0;
		int k = 0;
		for (IQuadEdge s : eList) {
			double z = s.getA().getZ();
			zSum += w[k++] * z;
		}
		return zSum;

	}

	IQuadEdge e;

	public List<IQuadEdge> getBowyerWatsonEnvelope(double x, double y) {
		// in the logic below, we access the Vertex x and y coordinates directly
		// but we use the getZ() method to get the z value. Some vertices
		// may actually be VertexMergerGroup instances

		ArrayList<IQuadEdge> eList = new ArrayList<>();
//		IQuadEdge e = null;
		if (v2 == null || !ptInTriangle(x, y, v0, v1, v2)) {
			IQuadEdge locatorEdge = navigator.getNeighborEdge(x, y);
			e = locatorEdge;
			v0 = e.getA();
			v1 = e.getB();
			v2 = e.getForward().getB();
			if (v2 == null) {
				return eList; // empty list, NNI undefined.
			}
		}

//		eList.clear();
		double h;

		// by the way the getNeighborEdge() method is defined, if
		// the query is outside the TIN or on the perimeter edge,
		// the edge v0, v1 will be the perimeter edge and v2 will
		// be the ghost vertex (e.g. a null). In either case, v2 will
		// not be defined. So, if v2 is null, the NNI interpolation is not defined.

		// ------------------------------------------------------
		// The fundamental idea of natural neighbor interpolation is
		// based on measuring how the local geometry of a Voronoi
		// Diagram would change if a new vertex were inserted.
		// (recall that the Voronoi is the dual of a Delaunay Triangulation).
		// Thus the NNI interpolation has common element with an
		// insertion into a TIN. In writing the code below, I have attempted
		// to preserve similarities with the IncrementalTIN insertion logic
		// where appropriate.
		//
		// Step 1 -----------------------------------------------------
		// Create an array of edges that would connect to the radials
		// from an inserted vertex if it were added at coordinates (x,y).
		// This array happens to describe a Thiessen Polygon around the
		// inserted vertex.
		ArrayDeque<IQuadEdge> stack = new ArrayDeque<>();
		IQuadEdge c, n0, n1;

		c = e;
		while (true) {
			n0 = c.getDual();
			n1 = n0.getForward();

			if (n1.getB() == null) {
				// the search has reached a perimeter edge
				// just add the edge and continue.
				h = -1;
			} else {
				// TODO use spatial data structure?
				// test for the Delaunay inCircle criterion.
				// see notes about efficiency in the IncrementalTIN class.
				double a11 = n0.getA().x - x;
				double a21 = n1.getA().x - x;
				double a31 = n1.getB().x - x;

				// column 2
				double a12 = n0.getA().y - y;
				double a22 = n1.getA().y - y;
				double a32 = n1.getB().y - y;

				h = (a11 * a11 + a12 * a12) * (a21 * a32 - a31 * a22)
						+ (a21 * a21 + a22 * a22) * (a31 * a12 - a11 * a32)
						+ (a31 * a31 + a32 * a32) * (a11 * a22 - a21 * a12);
			}

			if (h >= 0) {
				// The vertex is within the circumcircle the associated
				// triangle. The Thiessen triangle will extend to include
				// that triangle and, perhaps, its neighbors.
				// So continue the search.
				stack.addFirst(n0);
				c = n1;
			} else {
				eList.add(c);
				c = c.getForward();
				IQuadEdge p = stack.peekFirst();
				while (c.equals(p)) {
					stack.remove();
					c = c.getDual().getForward();
					p = stack.peekFirst();
				}
				if (c.equals(e)) {
					break;
				}
			}
		}

		return eList;
	}

	@Override
	public double[] getBarycentricCoordinates(List<IQuadEdge> polygon, double x, double y) {

		int nEdge = polygon.size();
		if (nEdge < 3) {
			return new double[0];
		}

		// The eList contains a series of edges definining the cavity
		// containing the polygon.
		Vertex a, b, c;

		double[] c0 = new double[2];
		double[] c1 = new double[2];
		double[] c2 = new double[2];
		double[] c3 = new double[2];

		IQuadEdge e0, e1, n, n1;
		double x0, y0, x1, y1, wThiessen, wXY, wDelta;
		double wSum = 0;
		double[] weights = new double[nEdge];
		// TODO cache edge calculation?
		for (int i0 = 0; i0 < nEdge; i0++) {
			int i1 = (i0 + 1);
			if (i1 == nEdge) {
				i1 = 0;
			}
			e0 = polygon.get(i0);
			e1 = polygon.get(i1);
			a = e0.getA();
			b = e1.getA(); // same as e0.getB();
			c = e1.getB();
			double ax = a.getX() - x;
			double ay = a.getY() - y;
			double bx = b.getX() - x;
			double by = b.getY() - y;
			double cx = c.getX() - x;
			double cy = c.getY() - y;

			x0 = (ax + bx) / 2;
			y0 = (ay + by) / 2;
			x1 = (bx + cx) / 2;
			y1 = (by + cy) / 2;

			// for the first edge processed, the code needs to initialize values
			// for c0 and c3. But after that, the code can reuse values from
			// the previous calculation.
			if (i0 == 0) {
				circumcircle(ax, ay, bx, by, 0d, 0d, c0);
				Vertex nb = e0.getForward().getB();
				circumcircle(ax, ay, bx, by, nb.getX() - x, nb.getY() - y, c3);
			} else {
				c0[0] = c1[0];
				c0[1] = c1[1];
			}

			circumcircle(bx, by, cx, cy, 0, 0, c1);

			// compute the reduced "component area" of the Theissen polygon
			// constructed around point B, the second point of edge[i0].
			wXY = (x0 * c0[1] - c0[0] * y0) + (c0[0] * c1[1] - c1[0] * c0[1]) + (c1[0] * y1 - x1 * c1[1]);

			// compute the full "component area" of the Theissen polygon
			// constructed around point B, the second point of edge[i0]
			n = e0.getForward();
			wThiessen = x0 * c3[1] - c3[0] * y0;
			while (!(n.equals(e1))) {
				n1 = n.getDual();
				n = n1.getForward();
				c2[0] = c3[0];
				c2[1] = c3[1];
				a = n1.getA();
				b = n.getA(); // same as n1.getB();
				c = n.getB();
				ax = a.getX() - x;
				ay = a.getY() - y;
				bx = b.getX() - x;
				by = b.getY() - y;
				cx = c.getX() - x;
				cy = c.getY() - y;
				circumcircle(ax, ay, bx, by, cx, cy, c3);
				wThiessen += c2[0] * c3[1] - c3[0] * c2[1];
			}
			wThiessen += c3[0] * y1 - x1 * c3[1];

			// Compute wDelta, the amount of area that the Theissen polygon
			// constructed around vertex B would yield to an insertion at
			// the query point.
			// for convenience, both the wXY and wThiessen weights were
			// computed in a clockwise order, which means they are the
			// negative of what we need for the weight computation, so
			// negate them and -(wTheissen-wXY) becomes wXY-wTheissen
			// Also, there would normally be a divide by 2 factor from the
			// shoelace area formula, but that is ommitted because it will
			// drop out when we unitize the sum of the set of the weights.
			wDelta = wXY - wThiessen;
			wSum += wDelta;
			weights[i1] = wDelta;
		}

		// Normalize the weights
		for (int i = 0; i < weights.length; i++) {
			weights[i] /= wSum;
		}

		return weights;
	}

	private static void circumcircle(double ax, double ay, double bx, double by, double cx, double cy,
			double[] result) {

		double D = (ax - cx) * (by - cy) - (bx - cx) * (ay - cy);
		double px = (((ax - cx) * (ax + cx) + (ay - cy) * (ay + cy)) / 2 * (by - cy)
				- ((bx - cx) * (bx + cx) + (by - cy) * (by + cy)) / 2 * (ay - cy))
				/ D;

		double py = (((bx - cx) * (bx + cx) + (by - cy) * (by + cy)) / 2 * (ax - cx)
				- ((ax - cx) * (ax + cx) + (ay - cy) * (ay + cy)) / 2 * (bx - cx))
				/ D;

		result[0] = px;
		result[1] = py;
	}

	private static boolean ptInTriangle(double x, double y, Vertex p0, Vertex p1, Vertex p2) {
		final double dX = x - p2.x;
		final double dY = y - p2.y;
		final double dX21 = p2.x - p1.x;
		final double dY12 = p1.y - p2.y;
		final double D = dY12 * (p0.x - p2.x) + dX21 * (p0.y - p2.y);
		final double s = dY12 * dX + dX21 * dY;
		final double t = (p2.y - p0.y) * dX + (p0.x - p2.x) * dY;
		if (D < 0) {
			return s <= 0 && t <= 0 && s + t >= D;
		}
		return s >= 0 && t >= 0 && s + t <= D;
	}

}
