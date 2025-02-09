package com.github.micycle1.distancefield;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Polygonal;
import org.tinfour.common.IConstraint;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.TriangleCollector;

import processing.core.PShape;
import processing.core.PVector;

/**
 * Methods lifted from PGS.
 */
class Triangulator {

	/**
	 * Generates a constrained Delaunay Triangulation from the given shape.
	 * <p>
	 * This method returns the triangulation in its raw form: a Triangulated
	 * Irregular Network (mesh).
	 *
	 * @param shape the shape whose vertices to generate a triangulation from
	 * @return Triangulated Irregular Network object (mesh)
	 * @see #delaunayTriangulationMesh(PShape, Collection, boolean, int, boolean)
	 */
	static IIncrementalTin delaunayTriangulationMesh(PShape shape) {
		return delaunayTriangulationMesh(shape, null, true, 0, true);
	}

	/**
	 * Generates a Delaunay Triangulation from the given shape. The triangulation
	 * can be both constrained (meaning the triangulation is masked by the original
	 * shape) and refined (meaning additional points are inserted, usually leading
	 * to more uniform triangle shapes and sizes).
	 * <p>
	 * This method returns the triangulation in its raw form: a Triangulated
	 * Irregular Network (mesh).
	 *
	 * @param shape         the shape whose vertices to generate a triangulation
	 *                      from. <b>Can be null</b>.
	 * @param steinerPoints A list of additional points to insert into the
	 *                      triangulation in addition to the vertices of the input
	 *                      shape. <b>Can be null</b>.
	 * @param constrain     Constrain the triangulation output using the shape
	 *                      boundary (from point set). With shapes, you'll probably
	 *                      want to this to be true.
	 * @param refinements   The number of triangulation refinement/subdivision
	 *                      passes to perform. Each pass inserts the centroids of
	 *                      every existing triangle into the triangulation. Should
	 *                      be 0 or greater (probably no more than 5).
	 * @param pretty        Whether to maintain the Delaunay nature when
	 *                      constraining the triangulation, and whether to check
	 *                      that centroid locations lie within the shape during
	 *                      refinement. When pretty=true, triangles in the
	 *                      triangulation may be slightly more regular in
	 *                      shape/size. There is a small performance overhead which
	 *                      becomes more considerable at higher refinement levels.
	 *                      When constrain=false and refinements=0, this argument
	 *                      has no effect.
	 * @return Triangulated Irregular Network object (mesh)
	 * @see #delaunayTriangulation(PShape, Collection, boolean, int, boolean)
	 */
	static IIncrementalTin delaunayTriangulationMesh(PShape shape, Collection<PVector> steinerPoints, boolean constrain,
			int refinements, boolean pretty) {
		Geometry g = Conversion.fromVertices(shape);
		final IncrementalTin tin = new IncrementalTin(10);

		final List<Vertex> vertices = new ArrayList<>();
		final Coordinate[] coords = g.getCoordinates();
		for (int i = 0; i < coords.length; i++) {
			vertices.add(new Vertex(coords[i].x, coords[i].y, 0, i));
		}

		tin.add(vertices, null); // initial triangulation

		int vertexIndex = coords.length;
		if (steinerPoints != null) {
			for (PVector v : steinerPoints) { // add steiner points
				tin.add(new Vertex(v.x, v.y, 0, vertexIndex++));
			}
		}

		if (refinements > 0) {

			final IndexedPointInAreaLocator pointLocator = new IndexedPointInAreaLocator(g);
			final ArrayList<Vertex> refinementVertices = new ArrayList<>();

			/*
			 * All vertices must be added to Tin before constraints are added (hence the
			 * need for pointLocator, since visitSimpleTriangles() visits triangles that lie
			 * outside the shape at this stage.
			 */
			for (int i = 0; i < refinements; i++) {
				/*
				 * Must be iterated like this so that triangles to refine are delaunay at each
				 * stage (rather than recurse the original triangles).
				 */
				refinementVertices.clear();
				TriangleCollector.visitSimpleTriangles(tin, t -> {
					if (t.getArea() > 50) { // don't refine small triangles
						final Coordinate center = centroid(t); // use centroid rather than circumcircle center
						if (pretty || pointLocator.locate(center) != Location.EXTERIOR) {
							refinementVertices.add(new Vertex(center.x, center.y, 0));
						}
					}
				});
				tin.add(refinementVertices, null); // add refinement (steiner) points
			}
		}

		if (constrain) {
			List<IConstraint> constraints = new ArrayList<>();
			for (int n = 0; n < g.getNumGeometries(); n++) {
				boolean exterior = true;

				if (g instanceof Polygonal) {
					LinearRingIterator lri = new LinearRingIterator(g.getGeometryN(n));
					for (LinearRing ring : lri) {
						final List<Vertex> points = new ArrayList<>();
						final Coordinate[] c = ring.getCoordinates();
						if (c.length == 0) {
							exterior = false;
							continue;
						}

						for (Coordinate element : c) {
							points.add(new Vertex(element.x, element.y, 0));
						}
						/*
						 * In Tinfour, the shape exterior must be CCW and the holes must be CW. This is
						 * true for most PShapes, but some shapes (like those created from fonts) may
						 * have the rings orientated the other way, which needs to be corrected.
						 */
						if ((exterior && !Orientation.isCCWArea(c)) || (!exterior && Orientation.isCCWArea(c))) {
							Collections.reverse(points);
						}
						constraints.add(new PolygonConstraint(points));
						exterior = false;
					}
				}
			}
			if (!constraints.isEmpty()) {
				tin.addConstraints(constraints, pretty);
			}
		}

		return tin;
	}

	/**
	 * Computes the centroid/barycentre of a triangle.
	 */
	private static Coordinate centroid(final SimpleTriangle t) {
		final Vertex a = t.getVertexA();
		final Vertex b = t.getVertexB();
		final Vertex c = t.getVertexC();
		double x = a.x + b.x + c.x;
		x /= 3;
		double y = a.y + b.y + c.y;
		y /= 3;
		return new Coordinate(x, y);
	}

	/**
	 * Provides convenient iteration of exterior and linear rings (if any) of a
	 * polygonal JTS geometry. Supports MultiGeometries.
	 *
	 * @author Michael Carleton
	 */
	static final class LinearRingIterator implements Iterable<LinearRing> {

		private LinearRing[] array;

		/**
		 * Constructs the iterator for the given geometry. The first ring returned by
		 * the iterator is the exterior ring; all other rings (if any) are interior
		 * rings.
		 *
		 * @param g input geometry
		 */
		public LinearRingIterator(Geometry g) {
			ArrayList<LinearRing> rings = new ArrayList<>(g.getNumGeometries());
			for (int i = 0; i < g.getNumGeometries(); i++) {
				Polygon poly = (Polygon) g.getGeometryN(i);
				rings.add(poly.getExteriorRing());
				for (int j = 0; j < poly.getNumInteriorRing(); j++) {
					rings.add(poly.getInteriorRingN(j));
				}
			}
			array = rings.toArray(new LinearRing[rings.size()]);
		}

		public LinearRing[] getLinearRings() {
			return array;
		}

		@Override
		public Iterator<LinearRing> iterator() {
			return new Iterator<>() {

				private int currentIndex = 0;

				@Override
				public boolean hasNext() {
					return currentIndex < array.length;
				}

				@Override
				public LinearRing next() {
					if (!hasNext()) {
						throw new NoSuchElementException();
					}
					return array[currentIndex++];
				}

				@Override
				public void remove() {
					throw new UnsupportedOperationException();
				}
			};
		}
	}

}
