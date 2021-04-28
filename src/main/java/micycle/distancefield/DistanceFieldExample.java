package micycle.distancefield;

import micycle.distancefield.DistanceField.Interpolator;
import micycle.pgs.PGS_Processing;
import micycle.pgs.utility.PoissonDistribution;
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
		size(1000, 1000, FX2D);
		smooth(2);
	}

	int n = 0;

	@Override
	public void setup() {
		frameRate(30);
		heatMap = new Gradient();
		heatMap.pushColor(color(20, 11, 52));
		heatMap.pushColor(color(132, 32, 107));
		heatMap.pushColor(color(229, 92, 48));
		heatMap.pushColor(color(255, 243, 99));
		heatMap.primeAnimation();
		heatMap.prime();

		PoissonDistribution pd = new PoissonDistribution(0);
		final int buffer = 30;
		var points = pd.generate(buffer, buffer, width - buffer, height - buffer, 25, 10);
		points.clear();
		randomSeed(12);
		for (int i = 0; i < 200; i++) {
			points.add(new PVector(random(buffer, width - buffer), random(buffer, height - buffer)));
		}
		shape = PGS_Processing.concaveHull(points, 65);

//		var square = createShape(RECT, 0, 0, 300, 300);
//		square = PGS_Transformation.translateTo(square, width / 2, height / 2);

//		shape = PTS.difference(shape, square);
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

		// ISOLINE
//		shape(shape);
//		image(shapeCache, 0, 0);
//		var d = mouseX;
//		PShape s = df.computeIsoline(Interpolator.Trianglular, new PVector(320, 320), d);
//		shape(s);
		// draw isoline distance circle
//		strokeWeight(8);
//		noFill();
//		stroke(0, 255, 0);
//		ellipse(320, 320, d * 2, d * 2);
//		if (frameCount > 30) {
//			n += 8;
//		}

		// DISTANCE FIELD
		heatMap.animate(-0.01f);
		df.computeField(Interpolator.Trianglular, new PVector(mouseX, mouseY));
		text(frameRate, 0, 0);
		image(shapeCache, 0, 0);

	}

}
