package micycle.distancefield.interpolator;

import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IIncrementalTinNavigator;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.Vertex;
import org.tinfour.interpolation.IVertexValuator;
import org.tinfour.interpolation.TriangularFacetInterpolator;

/**
 * 
 * Provides interpolation based on treating the surface as a collection of
 * planar triangular facets.
 * <p>
 * This is a faster version of TriangularFacetInterpolator that is optimised for
 * repeated iteration; it forgoes some robustness checks but that hasn't been a
 * problem (so far...).
 * <p>
 * Most of the speedup comes from first checking whether the coordinate in the
 * current call to interpolate() is contained by the same triangle as the
 * previous coordinate (which is often the case during grid-like iteration).
 * This saves use from walking the triangulation to find the next triangle that
 * contains the point every time.
 * 
 * @author Michael Carleton
 *
 */
public class TFIF extends TriangularFacetInterpolator {

	private Vertex v0, v1, v2;
	private double z0;
	private double ax, ay, az;
	private double nx, ny, nz;
	private IIncrementalTinNavigator navigator;

	public TFIF(IIncrementalTin tin) {
		super(tin);
		navigator = tin.getNavigator();
	}

	/**
	 * Doesn't support VertexMergerGroup.
	 * 
	 * <p>
	 * Most efficient when iterated calls to interpolate() are likely to be within
	 * the same triangle.
	 * 
	 * @param valuator ignored. Uses vertices getZ() values.
	 */
	@Override
	public double interpolate(double x, double y, IVertexValuator valuator) {

		/**
		 * Only call getNeighborEdge when containing triangle changes. Assumes
		 * interpolate is called on closely packed values throughout iteration.
		 */

		if (v2 == null || !ptInTriangle(x, y, v0, v1, v2)) { // (x % 2 == 0 &&
			IQuadEdge e = navigator.getNeighborEdge(x, y);
			v0 = e.getA();
			v1 = e.getB();
			v2 = e.getForward().getB();
			if (v2 == null) {
				return Double.NaN;
			}

			z0 = v0.getZ();
			double z1 = v1.getZ();

			ax = v1.x - v0.x;
			ay = v1.y - v0.y;
			az = z1 - z0;

			double z2 = v2.getZ();

			double bx = v2.x - v0.x;
			double by = v2.y - v0.y;
			double bz = z2 - z0;

			nx = ay * bz - az * by;
			ny = az * bx - ax * bz;
			nz = ax * by - ay * bx;
		}

		double sx = x - v0.x;
		double sy = y - v0.y;

		return z0 - (nx * sx + ny * sy) / nz;
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
