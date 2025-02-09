package com.github.micycle1.distancefield;

import com.github.micycle1.distancefield.DistanceField.Interpolator;

import micycle.peasygradients.gradient.Gradient;
import micycle.pgs.PGS_Conversion;
import micycle.pgs.PGS_Hull;
import micycle.pgs.PGS_Morphology;
import micycle.pgs.PGS_PointSet;
import processing.core.PApplet;
import processing.core.PShape;
import processing.core.PVector;

public class Demo extends PApplet {

	public static void main(String[] args) {
		System.setProperty("prism.allowhidpi", "false");
		System.setProperty("sun.java2d.uiScale", "1");
		PApplet.main(Demo.class);
	}

	DistanceField d;
	PShape shape;

	@Override
	public void settings() {
		size(1000, 1000);
	}

	@Override
	public void setup() {
		Gradient g = new Gradient(color(255, 79, 25), color(21, 8, 77), color(92, 230, 230));

		var points = PGS_PointSet.random(50, 50, width - 50, height - 50, 750);
		shape = PGS_Hull.concaveHullBFS(points, 0.1);
		
//		shape = PGS_PointSet.findShortestTour(points);
		shape = PGS_Morphology.smoothGaussian(shape, 30);
		shape.setFill(false);
		shape.setStroke(color(255));
		
		d = new DistanceField(this, PGS_Conversion.copy(shape), g);
	}

	@Override
	public void draw() {
		background(255);
		d.computeField(Interpolator.Trianglular, new PVector(mouseX, mouseY));
		shape(shape);
//		d.computeIsoline(Interpolator.NaturalNeighbor, new PVector(500,500), mouseX);
	}

}
