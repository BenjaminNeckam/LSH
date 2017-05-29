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
			this.points = new ArrayList<>();
			this.points.add(point);
		}else{
			this.points.add(point);
		}
	}

	public ArrayList<Point> getPoints() {
		if(this.points==null){
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
	public String toString(){
		if(points==null){
			this.points = new ArrayList<>();
		}
		String string = "Intervall: [" + getBucketInterval()[0] + ";" + getBucketInterval()[1]  + "]" 
				+ "\nNumber Elements: " + getPoints().size();
		return string;
	}

}
