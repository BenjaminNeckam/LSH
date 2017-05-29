package at.ac.unive.utils;

import java.util.ArrayList;

public class Bucket {
	ArrayList<Point> points;
	float[] bucketInterval;
	
	public Bucket(float[] bucketInterval){
		this.bucketInterval = bucketInterval;
	}
	
	public void addPoint(Point point){
		if(this.points==null){
			this.points = new ArrayList<Point>();
			this.points.add(point);
		}else{
			this.points.add(point);
		}
	}

	public ArrayList<Point> getPoints() {
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
	

}
