package at.ac.unive.utils;

import java.util.ArrayList;

public class Bucket {
	ArrayList<Point> points;
	float[] bucketInterval;
	int bucketHashCode;
	int clusterNumb;
	Centroid centroid;

	public Bucket(float[] bucketInterval, int bucketHashCode) {
		this.bucketInterval = bucketInterval;
		this.bucketHashCode = bucketHashCode;
	}

	public void addPoint(Point point) {
		if (this.points == null) {
			this.points = new ArrayList<Point>();
			point.setBucketHashCode(bucketHashCode);
			this.points.add(point);
		} else {
			this.points.add(point);
			point.setBucketHashCode(bucketHashCode);
		}
	}

	public ArrayList<Point> getPoints() {
		if (this.points == null) {
			this.points = new ArrayList<>();
		}
		return points;
	}

	public void setPoints(ArrayList<Point> points) {
		this.points = points;
	}

	public float[] getBucketInterval() {
		return bucketInterval;
	}

	public void setBucketInterval(float[] bucketInterval) {
		this.bucketInterval = bucketInterval;
	}

	@Override
	public String toString() {
		if (points == null) {
			this.points = new ArrayList<>();
		}
		String string = "Intervall: [" + getBucketInterval()[0] + ";" + getBucketInterval()[1] + "]"
				+ "\nNumber Elements: " + getPoints().size();
		return string;
	}

	public int getBucketHashCode() {
		return bucketHashCode;
	}

	public void setBucketHashCode(int bucketHashCode) {
		this.bucketHashCode = bucketHashCode;
	}

	public int getClusterNumb() {
		return clusterNumb;
	}

	public void setClusterNumb(int clusterNumb) {
		this.clusterNumb = clusterNumb;
	}

	public Centroid getCentroid() {
		return centroid;
	}

	public void setCentroid(Centroid centroid) {
		this.centroid = centroid;
	}
	
	

}
