package micycle.distancefield;

import micycle.distancefield.DistanceField.Interpolator;
import micycle.pts.PTS;
import micycle.pts.PTSMorphology;
import micycle.pts.utility.PoissonDistribution;
import peasyGradients.gradient.Gradient;
import processing.core.PApplet;
import processing.core.PImage;
import processing.core.PShape;
import processing.core.PVector;

public final class DistanceFieldExample extends PApplet {

	public static void main(String[] args) {
		PApplet.main(DistanceFieldExample.class, args);
	}

	DistanceField df;
	Gradient heatMap;
	PShape shape;
	PImage shapeCache;

	@Override
	public void settings() {
		System.setProperty("prism.allowhidpi", "false");
		size(800, 800, FX2D);
	}

	@Override
	public void setup() {
		heatMap = new Gradient();
		heatMap.pushColor(color(20, 11, 52));
		heatMap.pushColor(color(132, 32, 107));
		heatMap.pushColor(color(229, 92, 48));
		heatMap.pushColor(color(255, 243, 99));
//		heatMap.primeAnimation();
		heatMap.prime();

		PoissonDistribution pd = new PoissonDistribution(0);
		final int buffer = 30;
		var points = pd.generate(buffer, buffer, width - buffer, height - buffer, 25, 10);
		shape = PTSMorphology.concaveHull(points, 5);

		df = new DistanceField(this, shape, heatMap);

		shape.setFill(false);
		shape.setStroke(true);
		shape.setStrokeWeight(2);
		shape.setStroke(color(255));

		textAlign(LEFT, TOP);

		var t = createGraphics(width, height);
		t.smooth(4);
		t.beginDraw();
		t.loadPixels();
		t.shape(shape);
		t.endDraw();
		t.updatePixels();
		shapeCache = t.get();
	}

	@Override
	public void draw() {
		background(0, 0, 40);
		df.compute(Interpolator.Trianglular, new PVector(mouseX, mouseY));
//		heatMap.animate(-0.01f);
		image(shapeCache, 0, 0);
		text(frameRate, 0, 0);
	}

}
