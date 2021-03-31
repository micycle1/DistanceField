package micycle.distancefield;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.jgrapht.Graph;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.DefaultUndirectedWeightedGraph;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.Vertex;
import org.tinfour.interpolation.IInterpolatorOverTin;
import org.tinfour.standard.IncrementalTin;

import micycle.distancefield.interpolator.NNIF;
import micycle.distancefield.interpolator.TFIF;
import micycle.pts.PTSTriangulation;
import peasyGradients.gradient.Gradient;
import processing.core.PApplet;
import processing.core.PGraphics;
import processing.core.PImage;
import processing.core.PShape;
import processing.core.PVector;

/**
 * Distance fields for PShapes.
 * 
 * @author Michael Carleton
 *
 */
public class DistanceField {

	public enum Interpolator {
		/** slow, but great results */
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
	protected ArrayList<PVector> vertices;

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

		vertices = new ArrayList<PVector>();
		for (int i = 0; i < shape.getVertexCount(); i++) {
			vertices.add(shape.getVertex(i));
		}

		tin = PTSTriangulation.delaunayTriangulationTin(shape, null, true, 0, false);

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
	 */
	public void compute(Interpolator interpolator, final PVector origin) {
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
					taskList.add(new Thread(new NNIF(distanceMesh), startRow, rowsPerThread));
					break;
				}
				case Trianglular: {
					taskList.add(new Thread(new TFIF(distanceMesh), startRow, rowsPerThread));
					break;
				}
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

	// get array[][] of normalised distance field values (not colors)
	public void getField() {

	}

	public void dispose() {
		THREAD_POOL.shutdownNow();
	}

	private class Thread implements Callable<Boolean> {

		private final IInterpolatorOverTin interpolator;
		private final int startRow, rows;

		/**
		 * Each thread should get its own interpoaltor due to last-triangle caching
		 * 
		 * @param interpolator
		 * @param startRow
		 * @param rows
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
