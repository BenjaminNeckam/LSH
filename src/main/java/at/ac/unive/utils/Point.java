package at.ac.unive.utils;

import java.util.ArrayList;

public class Point {
	ArrayList<Float> coordinates;
	int clusterNumb;
	float hashValue1;
	float hashValue2;
	int bucketHashCode;
	int compareClass;

	public Point() {
	}

	public Point(ArrayList<Float> coordinates) {
		this.coordinates = coordinates;
	}

	public Point(ArrayList<Float> coordinates, int clusterNumb) {
		this.coordinates = coordinates;
		this.clusterNumb = clusterNumb;
	}

	public ArrayList<Float> getCoordinates() {
		return coordinates;
	}

	
	
	public int getCompareClass() {
		return compareClass;
	}

	public void setCompareClass(int compareClass) {
		this.compareClass = compareClass;
	}

	public void setCoordinates(ArrayList<Float> coordinates) {
		this.coordinates = coordinates;
	}

	public int getClusterNumb() {
		return clusterNumb;
	}

	public void setClusterNumb(int clusterNumb) {
		this.clusterNumb = clusterNumb;
	}

	public void addNewCoordinate(float value) {
		if (this.coordinates == null) {
			this.coordinates = new ArrayList<Float>();
			this.coordinates.add(value);
		} else {
			this.coordinates.add(value);
		}
	}

	/**
	 * Overwrite an exitisting coordinate
	 * 
	 * @param index
	 * @param value
	 */
	public void overwriteCoordinate(int index, float value) {
		if (this.coordinates == null) {
			this.coordinates = new ArrayList<Float>();
			this.coordinates.add(value);
		} else {
			this.coordinates.remove(index);
			this.coordinates.add(value);
		}
	}

	public float getHashValue() {
		return hashValue1;
	}

	public void setHashValue(float hashValue) {
		this.hashValue1 = hashValue;
	}
	

	public float getHashValue2() {
		return hashValue2;
	}

	public void setHashValue2(float hashValue2) {
		this.hashValue2 = hashValue2;
	}

	public int getBucketHashCode() {
		return bucketHashCode;
	}

	public void setBucketHashCode(int bucketHashCode) {
		this.bucketHashCode = bucketHashCode;
	}

	@Override
	public String toString() {
		String coord = "Cluster: " + getClusterNumb() + " -> [ ";

		for (int i = 0; i < coordinates.size(); i++) {
			if (i == 0) {
				coord += coordinates.get(i) + " ; ";
			} else if (i != coordinates.size() - 1) {
				coord += coordinates.get(i) + " ; ";
			} else {
				coord += coordinates.get(i);
			}
		}
		coord += "]";
		return coord;
	}

}
