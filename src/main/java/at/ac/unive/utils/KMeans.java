package at.ac.unive.utils;


import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.jfree.ui.RefineryUtilities;



import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

public class KMeans {
	static float min = Float.MAX_VALUE;
	static float max = -Float.MAX_VALUE;

	public static void lsh(ArrayList<Point> points, Float bucketSize, float minDist) throws Exception {
		long startTime = System.nanoTime();
		long startTimeV = System.nanoTime();
		Point h1 = generateHashVector();
		double estimatedTimeV = (System.nanoTime() - startTimeV) / 1000000000.0;
		System.out.println("\nElapsed Time Generating HashVector: " + estimatedTimeV + " seconds");

		long startTimeB = System.nanoTime();
		// TODO more than one hash
		TreeMap<Object, Bucket> buckets = initBucketsWithValues(points, h1, bucketSize);
		double estimatedTimeB = (System.nanoTime() - startTimeB) / 1000000000.0;
		System.out.println("\nElapsed Time Init Buckets and Hash Points: " + estimatedTimeB + " seconds");

		long startTimeC = System.nanoTime();
		ArrayList<Centroid> centroid = initCentroids(buckets, minDist);
		System.out.println("Centroids!!!: " + centroid.size());
		double estimatedTimeC = (System.nanoTime() - startTimeC) / 1000000000.0;
		System.out.println("\nElapsed Time Init Centroids: " + estimatedTimeC + " seconds");
		// for(Centroid c: centroid){
		// System.out.println("New bucketsize in Cluster " + c.getClusterNumb()
		// + ": " + c.getBuckets().size());
		// }
		// Loop for updating Centroid and collision check


		for (int i = 0; i < 10; i++) {
//			centroidUpdate(centroid, h1, buckets);
			bucketUpdate(centroid);
			for(Centroid c: centroid){
				
				c.updateCentroid();
			}
		}
		plot(centroid);
		
		double estimatedTime = (System.nanoTime() - startTime) / 1000000000.0;
		System.out.println("\nElapsed Time: " + estimatedTime + " seconds");
		// System.out.println(buckets);

	}

	private static void bucketUpdate(ArrayList<Centroid> centroid) {
		for(int i = 0; i < centroid.size(); i++){
			for(Bucket b: centroid.get(i).getBuckets())
			computeDistAndAssignCenter(b, centroid);
		}
	}

	private static TreeMap<Object, Bucket> initBucketsWithValues(ArrayList<Point> points, Point v, Float w)
			throws Exception {
		// Current Bucketsize
		TreeMap<Object, Bucket> buckets = new TreeMap<>();

		float hashValue = 0;
		// First hash and computing Bucketsize

		for (Point point : points) {
			hash(point, v);
			if (point.getHashValue() >= max) {
				max = point.getHashValue();
			}
			if (point.getHashValue() <= min) {
				min = point.getHashValue();
			}
		}

		// Set intervals of buckets
		// TODO Optimize?
		if (w == null) {
			w = (float) (Math.abs(max - min) / (Math.sqrt(points.size()) / Math.log(points.size())));
		}
		System.out.println("Bucketsize: " + w);
		System.out.println("MIN: " + min);
		System.out.println("MAX: " + max);
		int bucketNumbs = (int) Math.ceil(Math.sqrt(points.size()) / Math.log(points.size()));
		System.out.println("Bucketnumbers: " + bucketNumbs);
		float interval = min;
		for (int i = 0; i < bucketNumbs; i++) {
			buckets.put(Float.hashCode(interval),
					new Bucket(new float[] { interval, (float) (interval + (w - 0.0001)) }, Float.hashCode(interval)));
			// System.out.println("IntevalHash: " + Float.hashCode(intervall));
			interval += w;
		}
		System.out.println("Min: " + buckets.firstKey().hashCode() + "Max: " + buckets.lastKey().hashCode());
		// Set points to buckets
		for (Point point : points) {
			hashValue = point.getHashValue();
			if (hashValue < min || hashValue > max) {
				throw new Exception("Hashcode is not in range");
			} else {
				// buckets.get(hashCode).addPoint(point);
				if (Float.hashCode(hashValue) > 0) {
					buckets.floorEntry(Float.hashCode(hashValue)).getValue().addPoint(point);
				} else {
					buckets.ceilingEntry(Float.hashCode(hashValue)).getValue().addPoint(point);
				}

			}
		}
		return buckets;

	}

	public static ArrayList<Point> initCentroids(ArrayList<Point> points, int numOfCtr) {
		ArrayList<Point> pointArray = new ArrayList(points);
		ArrayList<Point> centroids = new ArrayList();
		Random r = new Random();

		Point firstCentroid = pointArray.get(r.nextInt(pointArray.size()));
		centroids.add(firstCentroid);
		pointArray.remove(firstCentroid);

		double[] d2 = new double[pointArray.size()];

		while (centroids.size() < numOfCtr) {
			double sum = 0;
			for (int i = 0; i < centroids.size(); ++i) {
				Point p = pointArray.get(i);
				int index = nearestCentroid(centroids, p);
				double d = dist(p, centroids.get(index));
				sum += d * d;
				d2[i] = sum;
			}

			double d = r.nextDouble() * sum;
			for (int i = 0; i < d2.length; ++i) {
				if (d2[i] >= d) {
					Point p = pointArray.get(i);
					centroids.add(p);
					pointArray.remove(p);
					break;
				}
			}
		}

		return centroids;

	}

	public static int nearestCentroid(ArrayList<Point> centroids, Point p) {
		double minDistance = Double.MAX_VALUE;
		int indexOfCentroid = 0;

		for (int i = 0; i < centroids.size(); ++i) {
			float distance = dist(p, centroids.get(i));
			if (distance < minDistance) {
				minDistance = distance;
				indexOfCentroid = i;
			}
		}

		return indexOfCentroid;
	}

	/**
	 * Function to compare results from Martin Perdacher
	 * 
	 * @param one
	 * @param two
	 * @return
	 */
	public static double NMI(ArrayList<Integer> one, ArrayList<Integer> two) {
		if (one.size() != two.size()) {
			throw new IllegalArgumentException("Sizes don't match!");
		}
		int maxone = Collections.max(one);
		int maxtwo = Collections.max(two);

		double[][] count = new double[maxone + 1][maxtwo + 1];
		// System.out.println(count[1][2]);
		for (int i = 0; i < one.size(); i++) {
			count[one.get(i)][two.get(i)]++;
		}
		// i<maxone=R
		// j<maxtwo=C
		double[] bj = new double[maxtwo + 1];
		double[] ai = new double[maxone + 1];

		for (int m = 0; m < maxtwo + 1; m++) {
			for (int l = 0; l < maxone + 1; l++) {
				bj[m] = bj[m] + count[l][m];
			}
		}
		for (int m = 0; m < maxone + 1; m++) {
			for (int l = 0; l < maxtwo + 1; l++) {
				ai[m] = ai[m] + count[m][l];
			}
		}

		double N = 0;
		for (int i = 0; i < ai.length; i++) {
			N = N + ai[i];
		}
		double HU = 0;
		for (int l = 0; l < ai.length; l++) {
			double c = 0;
			c = ai[l] / N;
			if (c > 0) {
				HU = HU - c * Math.log(c);
			}
		}

		double HV = 0;
		for (int l = 0; l < bj.length; l++) {
			double c = 0;
			c = bj[l] / N;
			if (c > 0) {
				HV = HV - c * Math.log(c);
			}
		}
		double HUstrichV = 0;
		for (int i = 0; i < maxone + 1; i++) {
			for (int j = 0; j < maxtwo + 1; j++) {
				if (count[i][j] > 0) {
					HUstrichV = HUstrichV - count[i][j] / N * Math.log(count[i][j] / bj[j]);
				}
			}
		}
		double IUV = HU - HUstrichV;
		double reto = IUV / Math.max(HU, HV);

		System.out.println("NMI: " + reto);
		return reto;
	}

	/**
	 * Normaldistributed vector for hashfunction
	 * 
	 * @return
	 */
	private static Point generateHashVector() {
		Point v;
		ArrayList<Float> coordinates = new ArrayList<Float>();
		Random random = new Random();
		float temp;
		// 10 weil Dim 10 im Moment
		for (int i = 0; i < 10; i++) {
			temp = (float) (0 + random.nextGaussian());
			coordinates.add(temp);
		}
		v = new Point(coordinates);
		return v;
	}

	private static void hash(Point point, Point hashVector) {
		float sum = 0;
		int size = point.getCoordinates().size();
		for (int i = 0; i < size; i++) {
			sum += point.getCoordinates().get(i) * hashVector.getCoordinates().get(i);
		}
		point.setHashValue(sum);
	}

	@SuppressWarnings("unchecked")
	private static HashMap<Object, Bucket> sortByValues(TreeMap<Object, Bucket> map) {
		List list = new LinkedList(map.entrySet());
		Collections.sort(list, new Comparator<Object>() {

			@Override
			public int compare(Object o1, Object o2) {
				Map.Entry<Object, Bucket> entry1 = (Entry<Object, Bucket>) o1;
				Map.Entry<Object, Bucket> entry2 = (Entry<Object, Bucket>) o2;
				if (entry1.getValue().getPoints().size() > entry2.getValue().getPoints().size()) {
					return -1;
				}
				if (entry1.getValue().getPoints().size() < entry2.getValue().getPoints().size()) {
					return 1;
				}
				return 0;
			}
		});

		HashMap<Object, Bucket> sortedTree = new LinkedHashMap();
		for (Object obj : list) {
			Map.Entry<Object, Bucket> entry = (Entry<Object, Bucket>) obj;
			sortedTree.put(entry.getKey(), entry.getValue());
		}
		return sortedTree;
	}

	public static float dist(Point point1, Point point2) {
		float sum = 0;
		for (int i = 0; i < point1.getCoordinates().size(); i++) {
			sum += Math.pow(point1.getCoordinates().get(i) - point2.getCoordinates().get(i), 2);
		}
		float norm = (float) Math.sqrt(sum);
		return norm;
	}

	// At the moment radical usage: compute similarity of a single or few random
	// point to all center and assign all points of the bucket to the center
	private static void computeDistAndAssignCenter(Bucket bucket, ArrayList<Centroid> centroid) {
		ArrayList<Point> bucketPoints = bucket.getPoints();
		float tmp;
		float tmpDist = 0;
		int size = bucket.getPoints().size();
//		if (bucketPoints.size() < 1000) {
//			size = bucketPoints.size();
//		} else {
//			size = 1000;
//		}
		//or intersection and check intersected points
		int centroidIdx = 0;
		for (int i = 0; i < size; i++) {
			tmp = dist(bucketPoints.get(i), centroid.get(0));
			tmpDist = 0;
			for (int j = 0; j < centroid.size(); j++) {
				tmpDist = dist(bucketPoints.get(i), centroid.get(j));
				if (tmpDist < tmp) {
					tmp = tmpDist;
//					bucket.setClusterNumb(centroid.get(j).getClusterNumb());
//					bucket.setCentroid(centroid.get(j));
					centroidIdx = j;
	
				}
			}
			bucketPoints.get(i).setClusterNumb(centroid.get(centroidIdx).getClusterNumb());
			centroid.get(centroidIdx).addPoint(bucketPoints.get(i));
		}
//		centroid.get(bucket.getClusterNumb()).addBucket(bucket);
	}

	public static void plot(ArrayList<Centroid> centroids) {
		Plot scatterplotdemo4 = new Plot("K-Means Start", centroids);
		scatterplotdemo4.pack();
		RefineryUtilities.centerFrameOnScreen(scatterplotdemo4);
		scatterplotdemo4.setVisible(true);
	}

	private static ArrayList<Centroid> initCentroids(TreeMap<Object, Bucket> buckets, float minDist) {
		// Returns a sorted HashMap by Values (size of pointslist)
		HashMap<Object, Bucket> sortedBuckets = sortByValues(buckets);
		Set<Object> sortedKeySet = sortedBuckets.keySet();
		Object[] keys = sortedKeySet.toArray();
		float distanceGood;
		float maxDistanceGood = Float.MIN_VALUE;
		int maxKey = 0;
		boolean flag = false;
		ArrayList<Bucket> bucketList = new ArrayList<>(buckets.values());
		ArrayList<Bucket> bucketListCopy = new ArrayList<>(buckets.values());
		sortList(bucketList);
		sortList(bucketListCopy);
		// Init Centroids (calcultateCentroidOfBucket just for first centroid
		// calculation)
		ArrayList<Centroid> centroid = new ArrayList<>();

		int counterKeys = 0;
		Centroid tmpCentroid;

		for (int i = 0; i < 15; i++) {
			counterKeys = i;
			flag = false;
			if (i == 0) {
				tmpCentroid = new Centroid();

				buckets.get(keys[i]).setClusterNumb(i);
				buckets.get(keys[i]).setCentroid(tmpCentroid);
				tmpCentroid.setClusterNumb(i);
				tmpCentroid.addBucket(buckets.get(keys[i]));
				tmpCentroid.setBucketHashCode(buckets.get(keys[i]).getBucketHashCode());
//				tmpCentroid.updateCentroid();
				tmpCentroid.initUpdateCentroid();
				centroid.add(tmpCentroid);
//				centroid.get(i).updateCentroid();
				tmpCentroid = null;
				counterKeys++;
			} else {
				maxDistanceGood = Float.MIN_VALUE;
				//TODO maby bug because counterkeys goes up till end
				for(int j=i;j<keys.length;j++){
					if(buckets.get(keys[j]).getCentroid() == null){
						counterKeys = j;
						break;
					}
				}
				
				do {
					Centroid maxCentroid = new Centroid();
					if(counterKeys<keys.length){
						tmpCentroid = new Centroid();
						
						
						tmpCentroid.setClusterNumb(i);
						tmpCentroid.addBucket(buckets.get(keys[counterKeys]));
						tmpCentroid.setBucketHashCode(buckets.get(keys[counterKeys]).getBucketHashCode());
//						tmpCentroid.updateCentroid();
						tmpCentroid.initUpdateCentroid();
//						distanceGood = isDistanceGood(centroid.get(i - 1), tmpCentroid, minDist);
						distanceGood  = dist(centroid.get(i - 1), tmpCentroid);
						if(distanceGood > maxDistanceGood){
							maxDistanceGood = distanceGood;
							maxCentroid = tmpCentroid;
							maxKey = counterKeys;
						}
						if (distanceGood >= minDist) {
							buckets.get(keys[counterKeys]).setClusterNumb(i);
							buckets.get(keys[counterKeys]).setCentroid(tmpCentroid);
							centroid.add(tmpCentroid);
//							centroid.get(i).updateCentroid();
//							tmpCentroid = null;
							flag = true;
						}
						counterKeys++;
					}else{
						buckets.get(keys[maxKey]).setClusterNumb(i);
						buckets.get(keys[maxKey]).setCentroid(maxCentroid);
//						maxCentroid.updateCentroid();
						maxCentroid.initUpdateCentroid();
						centroid.add(maxCentroid);
//						centroid.get(i).updateCentroid();
					}
					
				} while (!flag);
			}

		}
		
		for(Centroid c: centroid){
			c.updateCentroid();
		}
		
		for(int i = 0; i < keys.length; i++){
			if(buckets.get(keys[i]).getCentroid() == null){
				computeDistAndAssignCenter(buckets.get(keys[i]), centroid);
			}
		}

		// Calculate remaining buckets
		// for (int i = 15; i < keys.length; i++) {
		// if (buckets.get(keys[i]).getPoints() != null ||
		// buckets.get(keys[i]).getPoints().size() > 0) {
		// if (buckets.get(keys[i]).getPoints().size() != 0) {
		// computeDistAndAssignCenter(buckets.get(keys[i]), centroid);
		// }
		//
		// }
		// }
		return centroid;
	}

	private static void sortList(ArrayList<Bucket> buckets) {
		Collections.sort(buckets, new Comparator<Bucket>() {

			@Override
			public int compare(Bucket o1, Bucket o2) {

				if (o1.getPoints().size() > o2.getPoints().size()) {
					return -1;
				}
				if (o1.getPoints().size() < o2.getPoints().size()) {
					return 1;
				}
				return 0;
			}
		});
	}

	public static ArrayList<Integer> PointsToIntegerList(ArrayList<Point> points){
		ArrayList<Integer> list = new ArrayList<>();
		Integer tmp; 
		for(Point p: points){
			tmp = new Integer(p.getClusterNumb());
			list.add(tmp);
		}
		return list;
	}
}
