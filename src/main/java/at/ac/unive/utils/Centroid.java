package at.ac.unive.utils;

import java.util.ArrayList;

public class Centroid extends Point {
	ArrayList<Bucket> buckets;
	int bucketHashCode;

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
		float sum = 0;
		for (int i = 0; i < points.get(0).getCoordinates().size(); i++) {
			sum = 0;
			for (int j = 0; j < points.size(); j++) {
				sum += points.get(j).getCoordinates().get(i);
			}
			sum = sum / points.size();
			super.coordinates.add(sum);

		}
	}

	public ArrayList<Point> bucketToList(ArrayList<Bucket> buckets) {
		ArrayList<Point> points = new ArrayList<>();
		for (Bucket b : buckets) {
			points.addAll(b.getPoints());
		}
		return points;
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
