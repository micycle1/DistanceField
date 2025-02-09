package com.github.micycle1.distancefield;

import static processing.core.PConstants.BEZIER_VERTEX;
import static processing.core.PConstants.CURVE_VERTEX;
import static processing.core.PConstants.QUADRATIC_VERTEX;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateList;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.PrecisionModel;

import processing.core.PConstants;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Methods lifted from PGS.
 *
 */
class Conversion {

	private static final GeometryFactory GEOM_FACTORY = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING_SINGLE));

	/** Approximate distance between successive sample points on bezier curves */
	private static final float BEZIER_SAMPLE_DISTANCE = 2;

	/**
	 * Creates a JTS Polygon from a geometry or path PShape, whose 'kind' is a
	 * polygon or path.
	 */
	static Geometry fromVertices(PShape shape) {

		if (shape.getVertexCount() < 2) { // skip empty / point PShapes
			return GEOM_FACTORY.createPolygon();
		}

		int[] rawVertexCodes = shape.getVertexCodes();
		/*
		 * getVertexCodes() is null for P2D shapes with no contours created via
		 * createShape(), so need to instantiate array here.
		 */
		if (rawVertexCodes == null) {
			rawVertexCodes = new int[shape.getVertexCount()];
			Arrays.fill(rawVertexCodes, PConstants.VERTEX);
		}

		final int[] contourGroups = getContourGroups(rawVertexCodes);
		final int[] vertexCodes = getVertexTypes(rawVertexCodes);

		final List<CoordinateList> contours = new ArrayList<>(); // list of coords representing rings/contours

		int lastGroup = -1;
		for (int i = 0; i < shape.getVertexCount(); i++) {
			if (contourGroups[i] != lastGroup) {
				if (lastGroup == -1 && contourGroups[0] > 0) {
					lastGroup = 0;
				}
				lastGroup = contourGroups[i];
				contours.add(new CoordinateList());
			}

			/**
			 * Sample bezier curves at regular intervals to produce smooth Geometry
			 */
			switch (vertexCodes[i]) { // VERTEX, BEZIER_VERTEX, CURVE_VERTEX, or BREAK
			case QUADRATIC_VERTEX:
				contours.get(lastGroup).addAll(getQuadraticBezierPoints(shape.getVertex(i - 1), shape.getVertex(i),
						shape.getVertex(i + 1), BEZIER_SAMPLE_DISTANCE), false);
				i += 1;
				continue;
			case BEZIER_VERTEX: // aka cubic bezier
				contours.get(lastGroup).addAll(getCubicBezierPoints(shape.getVertex(i - 1), shape.getVertex(i),
						shape.getVertex(i + 1), shape.getVertex(i + 2), BEZIER_SAMPLE_DISTANCE), false);
				i += 2;
				continue;
			default: // VERTEX
				contours.get(lastGroup).add(coordFromPVector(shape.getVertex(i)), false);
				break;
			}
		}

		contours.forEach(contour -> {
			if (shape.isClosed()) {
				contour.closeRing();
			}
		});

		final Coordinate[] outerCoords = contours.get(0).toCoordinateArray();

		if (outerCoords.length == 0) {
			return GEOM_FACTORY.createPolygon(); // empty polygon
		} else if (outerCoords.length == 1) {
			return GEOM_FACTORY.createPoint(outerCoords[0]);
		} else if (outerCoords.length == 2) {
			return GEOM_FACTORY.createLineString(outerCoords);
		} else if (shape.isClosed()) { // closed geometry or path // assume all contours beyond the first represent
										// holes
			LinearRing outer = GEOM_FACTORY.createLinearRing(outerCoords); // should always be valid
			LinearRing[] holes = new LinearRing[contours.size() - 1]; // Create linear ring for each hole in the shape
			for (int j = 1; j < contours.size(); j++) {
				final Coordinate[] innerCoords = contours.get(j).toCoordinateArray();
				holes[j - 1] = GEOM_FACTORY.createLinearRing(innerCoords);
			}
			return GEOM_FACTORY.createPolygon(outer, holes);

		} else { // not closed
			return GEOM_FACTORY.createLineString(outerCoords);
		}
	}

	/**
	 * For every vertexcode, store the group (i.e. hole) it belongs to.
	 *
	 * @param vertexCodes
	 * @return
	 */
	private static int[] getContourGroups(int[] vertexCodes) {

		int group = 0;
		List<Integer> groups = new ArrayList<>(vertexCodes.length * 2);

		for (int i = 0; i < vertexCodes.length; i++) {
			final int vertexCode = vertexCodes[i];
			switch (vertexCode) {
			case PConstants.VERTEX:
				groups.add(group);
				break;

			case QUADRATIC_VERTEX:
				groups.add(group);
				groups.add(group);
				break;

			case BEZIER_VERTEX:
				groups.add(group);
				groups.add(group);
				groups.add(group);
				break;

			case CURVE_VERTEX:
				groups.add(group);
				break;

			case PConstants.BREAK:
				// BREAK marks beginning/end of new contour, and should be proceeded by a VERTEX
				if (i > 0) {
					// In P2D, svg-loaded shapes begin with a break (so we don't want to increment)
					group++;
				}
				break;
			default:
				System.err.println("Unrecognised vertex code: " + vertexCode);
			}
		}

		final int[] vertexGroups = new int[groups.size()];
		Arrays.setAll(vertexGroups, groups::get);
		return vertexGroups;
	}

	/**
	 * Basically getVertexCodes, but returns the vertex type for every vertex
	 *
	 * @param shape
	 * @return
	 */
	private static int[] getVertexTypes(int[] rawVertexCodes) {

		List<Integer> codes = new ArrayList<>(rawVertexCodes.length);

		for (int vertexCode : rawVertexCodes) {
			switch (vertexCode) {
			case PConstants.VERTEX:
				codes.add(PConstants.VERTEX);
				break;

			case QUADRATIC_VERTEX:
				codes.add(QUADRATIC_VERTEX);
				codes.add(QUADRATIC_VERTEX);
				break;

			case BEZIER_VERTEX:
				codes.add(BEZIER_VERTEX);
				codes.add(BEZIER_VERTEX);
				codes.add(BEZIER_VERTEX);
				break;

			case CURVE_VERTEX:
				codes.add(CURVE_VERTEX);
				break;

			case PConstants.BREAK:
				break;

			default:
				System.err.println("Unrecognised vertex code: " + vertexCode);
			}
		}

		final int[] vertexGroups = new int[codes.size()];
		Arrays.setAll(vertexGroups, codes::get);
		return vertexGroups;
	}

	/**
	 * Subdivide/interpolate/discretise along a quadratic bezier curve, given by its
	 * start, end and control points
	 *
	 * @return list of points along curve
	 */
	private static List<Coordinate> getQuadraticBezierPoints(PVector start, PVector controlPoint, PVector end,
			float sampleDistance) {
		final List<Coordinate> coords;

		if (start.dist(end) <= sampleDistance) {
			coords = new ArrayList<>(2);
			coords.add(coordFromPVector(start));
			coords.add(coordFromPVector(end));
			return coords;
		}

		final float length = bezierLengthQuadratic(start, controlPoint, end);
		final int samples = (int) Math.ceil(length / sampleDistance); // sample every x unit length (approximately)
		coords = new ArrayList<>(samples);

		coords.add(coordFromPVector(start));
		for (int j = 1; j < samples; j++) { // start at 1 -- don't sample at t=0
			final PVector bezierPoint = getQuadraticBezierCoordinate(start, controlPoint, end, j / (float) samples);
			coords.add(coordFromPVector(bezierPoint));
		}
		coords.add(coordFromPVector(end));

		return coords;
	}

	/**
	 *
	 * @param start
	 * @param controlPoint
	 * @param end
	 * @param t            0...1
	 * @return
	 */
	private static PVector getQuadraticBezierCoordinate(PVector start, PVector controlPoint, PVector end, float t) {
		float x = (1 - t) * (1 - t) * start.x + 2 * (1 - t) * t * controlPoint.x + t * t * end.x;
		float y = (1 - t) * (1 - t) * start.y + 2 * (1 - t) * t * controlPoint.y + t * t * end.y;
		return new PVector(x, y);
	}

	/**
	 * Approximate bezier length using Gravesen's approach. The insight is that the
	 * actual bezier length is always somewhere between the distance between the
	 * endpoints (the length of the chord) and the perimeter of the control polygon.
	 * For a quadratic Bézier, 2/3 the first + 1/3 the second is a reasonably good
	 * estimate.
	 *
	 * @return
	 */
	private static float bezierLengthQuadratic(PVector start, PVector controlPoint, PVector end) {
		// https://raphlinus.github.io/curves/2018/12/28/bezier-arclength.html
		final float chord = PVector.sub(end, start).mag();
		final float cont_net = PVector.sub(start, controlPoint).mag() + PVector.sub(end, controlPoint).mag();
		return (2 * chord + cont_net) / 3f;

	}

	/**
	 * Generates a list of samples of a cubic bezier curve.
	 *
	 * @param sampleDistance distance between successive samples on the curve
	 * @return
	 */
	private static List<Coordinate> getCubicBezierPoints(PVector start, PVector controlPoint1, PVector controlPoint2, PVector end,
			float sampleDistance) {
		final List<Coordinate> coords;

		if (start.dist(end) <= sampleDistance) {
			coords = new ArrayList<>(2);
			coords.add(coordFromPVector(start));
			coords.add(coordFromPVector(end));
			return coords;
		}

		final float length = bezierLengthCubic(start, controlPoint1, controlPoint2, end);
		final int samples = (int) Math.ceil(length / sampleDistance); // sample every x unit length (approximately)
		coords = new ArrayList<>(samples);

		coords.add(coordFromPVector(start));
		for (int j = 1; j < samples; j++) { // start at 1 -- don't sample at t=0
			final PVector bezierPoint = getCubicBezierCoordinate(start, controlPoint1, controlPoint2, end, j / (float) samples);
			coords.add(coordFromPVector(bezierPoint));
		}
		coords.add(coordFromPVector(end));
		return coords;
	}

	private static PVector getCubicBezierCoordinate(PVector start, PVector controlPoint1, PVector controlPoint2, PVector end,
			float t) {
		final float t1 = 1.0f - t;
		float x = start.x * t1 * t1 * t1 + 3 * controlPoint1.x * t * t1 * t1 + 3 * controlPoint2.x * t * t * t1
				+ end.x * t * t * t;
		float y = start.y * t1 * t1 * t1 + 3 * controlPoint1.y * t * t1 * t1 + 3 * controlPoint2.y * t * t * t1
				+ end.y * t * t * t;
		return new PVector(x, y);
	}

	/**
	 * Approximate bezier length using Gravesen's approach. The insight is that the
	 * actual bezier length is always somewhere between the distance between the
	 * endpoints (the length of the chord) and the perimeter of the control polygon.
	 *
	 * @return
	 */
	private static float bezierLengthCubic(PVector start, PVector controlPoint1, PVector controlPoint2, PVector end) {
		// https://stackoverflow.com/a/37862545/9808792
		final float chord = PVector.sub(end, start).mag();
		final float cont_net = PVector.sub(start, controlPoint1).mag() + PVector.sub(controlPoint2, controlPoint1).mag()
				+ PVector.sub(end, controlPoint2).mag();
		return (cont_net + chord) / 2;

	}

	private static final Coordinate coordFromPVector(PVector p) {
		return new Coordinate(p.x, p.y);
	}
}
