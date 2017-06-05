package at.ac.unive.utils;

import java.util.ArrayList;
import java.util.Random;

public class Centroid extends Point {
	ArrayList<Bucket> buckets;
	int bucketHashCode;
	ArrayList<Point> remPoints;
	ArrayList<Point> completePoints;

	public ArrayList<Bucket> getBuckets() {
		return buckets;
	}

	public void setBucket(ArrayList<Bucket> buckets) {
		this.buckets = buckets;
	}

	public int getBucketHashCode() {
		return bucketHashCode;
	}

	public void setBucketHashCode(int bucketHashCode) {
		this.bucketHashCode = bucketHashCode;
	}

	public void updateCentroid() {
		if (super.coordinates == null) {
			super.coordinates = new ArrayList<>();
		}
		super.coordinates.clear();
		ArrayList<Point> points = bucketToList(buckets);
		if(remPoints!=null){
			points.addAll(remPoints);
		}
		
		float sum = 0;
		for (int i = 0; i < points.get(0).getCoordinates().size(); i++) {
			sum = 0;
			for (int j = 0; j < points.size(); j++) {
				sum += points.get(j).getCoordinates().get(i);
			}
			sum = sum / (float)points.size();
			super.coordinates.add(sum);
		}
		
		completePoints = completePointList();
	}
	
	public void addPoint(Point point) {
		if(this.remPoints == null){
			remPoints = new ArrayList<>();
		}
		if(!remPoints.contains(point)){
			point.setClusterNumb(clusterNumb);
			//point.setBucketHashCode1(bucketHashCode1);
			remPoints.add(point);
		}
	}
	
	public void initUpdateCentroid(){
		Random rand = new Random();
		ArrayList<Point> points = bucketToList(buckets);
		Point center = points.get(rand.nextInt(points.size()));
		super.coordinates = new ArrayList<>();
		super.coordinates.clear();
		super.coordinates.add(center.getCoordinates().get(0));
		super.coordinates.add(center.getCoordinates().get(1));
	}

	public ArrayList<Point> bucketToList(ArrayList<Bucket> buckets) {
		ArrayList<Point> points = new ArrayList<>();
		for (Bucket b : buckets) {
			points.addAll(b.getPoints());
		}
		return points;
	}

	
	
	public ArrayList<Point> getRemPoints() {
		return remPoints;
	}

	public void setRemPoints(ArrayList<Point> remPoints) {
		this.remPoints = remPoints;
	}

	public void addBucket(Bucket bucket) {
		if (buckets == null) {
			buckets = new ArrayList<>();
			ArrayList<Point> points = bucket.getPoints();
			for (Point p : points) {
				p.setClusterNumb(getClusterNumb());
			}
			bucket.setClusterNumb(getClusterNumb());
			buckets.add(bucket);
		} else {
			ArrayList<Point> points = bucket.getPoints();
			for (Point p : points) {
				p.setClusterNumb(getClusterNumb());
			}
			bucket.setClusterNumb(getClusterNumb());
			buckets.add(bucket);
		}
	}

	public void removeBucket(Bucket bucket) {
		buckets.remove(bucket);
	}
	
public ArrayList<Point> completePointList(){
		ArrayList<Point> allPoints = new ArrayList<>();
		if(remPoints!=null){
			allPoints.addAll(remPoints);
		}
		if(buckets != null){
			allPoints.addAll(bucketToList(buckets));
		}
		
		return allPoints;
	}



	public ArrayList<Point> getCompletePoints() {
	return completePoints;
}

public void setCompletePoints(ArrayList<Point> completePoints) {
	this.completePoints = completePoints;
}

	@Override
	public String toString() {
		String coord = "BucketHashCode: " + bucketHashCode + "; Clusternumber: " + super.clusterNumb + "-> [ ";

		for (int i = 0; i < super.coordinates.size(); i++) {
			if (i == 0) {
				coord += super.coordinates.get(i) + " ; ";
			} else if (i != super.coordinates.size() - 1) {
				coord += super.coordinates.get(i) + " ; ";
			} else {
				coord += getCoordinates().get(i);
			}
		}
		coord += "]";
		return coord;
	}

}
