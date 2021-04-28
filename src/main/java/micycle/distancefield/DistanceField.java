package micycle.distancefield;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.function.Consumer;

import org.jgrapht.Graph;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.DefaultUndirectedWeightedGraph;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.tinfour.common.IConstraint;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.Vertex;
import org.tinfour.interpolation.IInterpolatorOverTin;
import org.tinfour.interpolation.NaturalNeighborInterpolator;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.utils.TriangleCollector;

import micycle.distancefield.interpolator.NNIF;
import micycle.distancefield.interpolator.TFIF;
import micycle.pgs.PGS_Conversion;
import micycle.pgs.PGS_Triangulation;

import peasyGradients.gradient.Gradient;
import processing.core.PApplet;
import processing.core.PGraphics;
import processing.core.PImage;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Distance fields for PShapes.
 * <p>
 * TODO: dist between origin and chosen pixel
 * 
 * @author Michael Carleton
 *
 */
public class DistanceField {

	// see https://github.com/openrndr/orx/tree/master/orx-jumpflood

	public enum Interpolator {
		/** faster than NN, but doesn't work on shapes with holes */
		NaturalNeighborQuick,
		// https://github.com/gwlucastrig/Tinfour/wiki/Introduction-to-Natural-Neighbor-Interpolation
		NaturalNeighbor,
		/** much faster than NN, but worse quality */
		Trianglular
	}

	protected final PApplet p;
	protected final PShape shape;
	protected final Gradient gradient;

	protected IncrementalTin tin;
	protected final Graph<Vertex, IQuadEdge> graph;
	protected DijkstraShortestPath<Vertex, IQuadEdge> shortestPaths;
	protected PGraphics shapeRaster;

	protected int[] colorMap;

	private final ExecutorService THREAD_POOL;
	protected final int THREAD_COUNT;
	private double maxD;

	public DistanceField(PApplet p, PShape shape, Gradient gradient) {
		this.p = p;
		p.loadPixels();
		this.shape = shape;
		this.gradient = gradient;

		THREAD_COUNT = Math.min(4, Runtime.getRuntime().availableProcessors() - 1);
		THREAD_POOL = Executors.newFixedThreadPool(THREAD_COUNT);

		final int REFINEMENTS = 0;
		tin = PGS_Triangulation.delaunayTriangulationTin(shape, null, true, REFINEMENTS, true);

		graph = new DefaultUndirectedWeightedGraph<>(IQuadEdge.class);
		tin.edges().forEach(e -> {
			if (e.isConstrainedRegionInterior() || e.isConstrainedRegionBorder()) {
				graph.addVertex(e.getA());
				graph.addVertex(e.getB());
				graph.addEdge(e.getA(), e.getB(), e);
				graph.setEdgeWeight(e.getA(), e.getB(), e.getLength());
			}
		});

		shape.setFill(true);
		shape.setFill(p.color(255));
		shape.setStroke(false);
		shapeRaster = p.createGraphics(p.width, p.height);
		shapeRaster.beginDraw();
		shapeRaster.loadPixels();
		shapeRaster.shape(shape);
		shapeRaster.endDraw();
		shapeRaster.updatePixels();

		p.registerMethod("dispose", this);
	}

	public PImage compute() {
		return null;
	}

	/**
	 * Compute the distance field for the shape using the vertex of the PShape
	 * CLOSEST to given origin point.
	 * 
	 * @param interpolator
	 * @param origin
	 * @return
	 */
	public void computeField(Interpolator interpolator, final PVector origin) {
		IncrementalTin distanceMesh = new IncrementalTin(10);

		shortestPaths = new DijkstraShortestPath<>(graph);
		var originVertex = tin.getNavigator().getNearestVertex(origin.x, origin.y); // vertices only
		double maxDistance = 0;
		if (graph.containsVertex(originVertex)) {
			var paths = shortestPaths.getPaths(originVertex);
			for (Vertex v : graph.vertexSet()) {
				var d = paths.getWeight(v);
				distanceMesh.add(new Vertex(v.x, v.y, d));
				maxDistance = Math.max(maxDistance, d);
			}
		}

		final int l = 1000;
		colorMap = new int[l];
		for (int i = 0; i < l; i++) { // calc LUT
			colorMap[i] = gradient.getColor(i / (float) l);
		}
		maxD = (1 / maxDistance) * colorMap.length;

		p.loadPixels();
		List<Callable<Boolean>> taskList = new ArrayList<>();

		int rowsPerThread = p.height / THREAD_COUNT;
		int startRow = 0;
		for (int i = 0; i < THREAD_COUNT; i++) {
			switch (interpolator) {
				case NaturalNeighbor: {
					taskList.add(new Thread(new NaturalNeighborInterpolator(distanceMesh), startRow, rowsPerThread));
					break;
				}
				case Trianglular: {
					taskList.add(new Thread(new TFIF(distanceMesh), startRow, rowsPerThread));
					break;
				}
				case NaturalNeighborQuick:
					taskList.add(new Thread(new NNIF(distanceMesh), startRow, rowsPerThread));
					break;
				default:
					break;
			}
			startRow += rowsPerThread;
		}

		try {
			THREAD_POOL.invokeAll(taskList);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		p.fill(p.color(255, 0, 0));
		p.noStroke();
		p.ellipse(origin.x, origin.y, 10, 10);
		p.ellipse((float) originVertex.x, (float) originVertex.y, 10, 10);
	}

	/**
	 * Compute the distance line for the shape using the vertex of the PShape
	 * CLOSEST to given origin point.
	 * 
	 * @param interpolator
	 * @param origin
	 * @return
	 */
	public PShape computeIsoline(Interpolator interpolator, final PVector origin, double distance) {
		IncrementalTin distanceMesh = new IncrementalTin(10);
		IncrementalTin constrainedMesh = new IncrementalTin(10);

		shortestPaths = new DijkstraShortestPath<>(graph);
		var originVertex = tin.getNavigator().getNearestVertex(origin.x, origin.y); // vertices only
		double maxDistance = 0;
		if (graph.containsVertex(originVertex)) {
			var paths = shortestPaths.getPaths(originVertex);
			for (Vertex v : graph.vertexSet()) {
				var d = paths.getWeight(v);
				distanceMesh.add(new Vertex(v.x, v.y, d));
				constrainedMesh.add(new Vertex(v.x, v.y, d));
				maxDistance = Math.max(maxDistance, d);
			}
		}

		List<IConstraint> constraints = new ArrayList<>();
		var ring = PGS_Conversion.fromPShape(shape);
		ArrayList<Vertex> points = new ArrayList<>();
		Coordinate[] c = ring.getCoordinates();
		if (Orientation.isCCW(c)) {
			for (int i = 0; i < c.length; i++) {
				points.add(new Vertex(c[i].x, c[i].y, 0));
			}
		} else {
			for (int i = c.length - 1; i >= 0; i--) { // iterate backwards if CW
				points.add(new Vertex(c[i].x, c[i].y, 0));
			}
		}
		constraints.add(new PolygonConstraint(points));

		constrainedMesh.addConstraints(constraints, true); // true/false is negligible?

		PShape isoline = new PShape(PShape.GEOMETRY);
		isoline.setStroke(true);
		isoline.setStrokeCap(PShape.ROUND);
		isoline.setStroke(p.color(255, 0, 0));
		isoline.setStrokeWeight(10);

		isoline.beginShape(PShape.LINES);
		TriangleCollector.visitTrianglesConstrained(constrainedMesh, new Consumer<Vertex[]>() {

			@Override
			public void accept(Vertex[] t) {

				if (t[0].isConstraintMember() || t[1].isConstraintMember() || t[2].isConstraintMember()) {
					return;
				}
				p.beginShape(PShape.LINES);
				p.stroke(PApplet.map((float) t[0].x, 0, p.width, 0, 255), 20,
						PApplet.map((float) t[0].y, 0, p.width, 255, 0));

				int n = 0;
				if (isoVertex(t[0], t[1], isoline, distance)) {
					n++;
				}
				if (isoVertex(t[1], t[2], isoline, distance)) {
					n++;
				}
				if (isoVertex(t[2], t[0], isoline, distance)) {
					n++;
				}
				if (n == 2) { // isoline crosses triangle
//					drawTriangle(t[0], t[1], t[2]);
				}
				p.endShape();
			}
		});

		isoline.endShape();

		return isoline;
	}

	/**
	 * Compute (if possible) the coordinate where the isoline crosses the between
	 * the two vertices.
	 * 
	 * @return
	 */
	private boolean isoVertex(Vertex a, Vertex b, PShape s, double d) {
		Vertex min, max;

		if (a.getZ() > b.getZ()) {
			max = a;
			min = b;
		} else {
			max = b;
			min = a;
		}

		if (d > min.getZ() && d < max.getZ()) {
			double diff = max.getZ() - min.getZ();
			double numerator = d - min.getZ();

			double fract = numerator / diff;
			double xDiff = max.getX() - min.getX();
			double yDiff = max.getY() - min.getY();
			p.vertex((float) (min.getX() + fract * xDiff), (float) (min.getY() + fract * yDiff));
//			return new PVector((float) (min.getX() + fract * xDiff), (float) (min.getY() + fract * yDiff));
			return true;
		}
		return false;
	}

	private void drawTriangle(Vertex a, Vertex b, Vertex c) {
		var v = a;
		if (v != null && b != null && c != null) {
			p.fill(120, PApplet.map((float) a.x, 0, p.width, 0, 255), PApplet.map((float) a.y, 0, p.height, 0, 255));
			p.stroke(120, PApplet.map((float) a.x, 0, p.width, 0, 255), PApplet.map((float) a.y, 0, p.height, 0, 255));
			p.beginShape();
			p.strokeWeight(1);
			p.vertex((float) v.x, (float) v.y);
			v = b;
			p.vertex((float) v.x, (float) v.y);
			v = c;
			p.vertex((float) v.x, (float) v.y);
			p.endShape();
		}
	}

	// get array[][] of normalised distance field values (not colors)
	public void getField() {

	}

	/**
	 * Public to be called by PApplet when window closed.
	 */
	public void dispose() {
		THREAD_POOL.shutdownNow();
	}

	private class Thread implements Callable<Boolean> {

		private final IInterpolatorOverTin interpolator;
		private final int startRow, rows;

		/**
		 * Each thread should get its own interpoaltor due to last-triangle caching.
		 * Threads draw into rows that span the whole width of the PApplet.
		 * 
		 * @param interpolator
		 * @param startRow     index starting row for this thread
		 * @param rows         number of rows this thread should calculate and render
		 *                     for
		 */
		Thread(IInterpolatorOverTin interpolator, int startRow, int rows) {
			this.interpolator = interpolator;
			this.startRow = startRow;
			this.rows = rows;
		}

		@Override
		public Boolean call() throws Exception {
			int pixel = startRow * p.width;
			for (int y = 0; y < rows; y++) {
				for (int x = 0; x < p.width; x++) {
					if (shapeRaster.pixels[pixel] != 0) {
						final double z = interpolator.interpolate(x, y + startRow, null);
						if (!Double.isNaN(z)) {
							int c = (int) (z * maxD);
							p.pixels[pixel] = colorMap[c];
						}
					}
					pixel++;
				}
			}
			p.updatePixels(0, startRow, p.width, startRow + rows);
			return true;
		}
	}
}
