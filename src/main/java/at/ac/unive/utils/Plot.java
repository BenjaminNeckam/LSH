package at.ac.unive.utils;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Shape;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.util.ShapeUtilities;

public class Plot extends ApplicationFrame {

	int dimension;

	/**
	 * Constructor
	 * 
	 * @param title
	 * @param points
	 * @param d
	 */
	public Plot(String title, ArrayList<Point> points, int d) {
		super(title);
		this.dimension = d;
		JPanel jpanel = createDemoPanel(points);
		jpanel.setPreferredSize(new Dimension(640, 480));
		add(jpanel);
	}

	public Plot(String title, ArrayList<Centroid> centroids) {
		super(title);
		this.dimension = 15;
		JPanel jpanel = createDemoPanelCentroids(centroids);
		jpanel.setPreferredSize(new Dimension(640, 480));
		add(jpanel);
	}

	/**
	 * Properties for Panel to visualize dataset
	 * 
	 * @param points
	 * @return
	 */
	private JPanel createDemoPanel(ArrayList<Point> points) {
		final JFreeChart jfreechart = ChartFactory.createScatterPlot("K-Means", "X", "Y", createDataset(points),
				PlotOrientation.VERTICAL, true, true, false);

		Shape cross = ShapeUtilities.createDiagonalCross(3, 1);
		XYPlot xyPlot = (XYPlot) jfreechart.getPlot();
		xyPlot.setDomainCrosshairVisible(true);
		xyPlot.setRangeCrosshairVisible(true);
		XYItemRenderer renderer = xyPlot.getRenderer();
		renderer.setSeriesShape(0, cross);
		renderer.setSeriesPaint(0, Color.red);
		renderer.setSeriesPaint(1, Color.blue);
		renderer.setSeriesPaint(2, Color.green);
		renderer.setSeriesPaint(3, Color.black);

		return new ChartPanel(jfreechart);
	}

	private JPanel createDemoPanelCentroids(ArrayList<Centroid> centroids) {
		final JFreeChart jfreechart = ChartFactory.createScatterPlot("K-Means", "X", "Y",
				createDatasetCentroids(centroids), PlotOrientation.VERTICAL, true, true, false);

		Shape cross = ShapeUtilities.createDiagonalCross(3, 1);
		XYPlot xyPlot = (XYPlot) jfreechart.getPlot();
		xyPlot.setDomainCrosshairVisible(true);
		xyPlot.setRangeCrosshairVisible(true);
		XYItemRenderer renderer = xyPlot.getRenderer();
		renderer.setSeriesShape(0, cross);
		renderer.setSeriesPaint(0, Color.red);
		renderer.setSeriesPaint(1, Color.blue);
		renderer.setSeriesPaint(2, Color.green);
		renderer.setSeriesPaint(3, Color.black);

		return new ChartPanel(jfreechart);
	}

	/**
	 * Creates dataset for plot
	 * 
	 * @param points
	 * @return
	 */
	private XYDataset createDataset(ArrayList<Point> points) {
		XYSeriesCollection collection = new XYSeriesCollection();
		List<XYSeries> series = new ArrayList<XYSeries>();

		for (int i = 0; i < dimension; i++) {
			series.add(new XYSeries("Cluster " + (i + 1)));
		}

		for (Point point : points) {
			series.get(point.getClusterNumb()).add(point.getCoordinates().get(0), point.getCoordinates().get(1));
		}

		for (int i = 0; i < dimension; i++) {
			collection.addSeries(series.get(i));
		}

		return collection;
	}

	private XYDataset createDatasetCentroids(ArrayList<Centroid> centroids) {
		XYSeriesCollection collection = new XYSeriesCollection();
		List<XYSeries> series = new ArrayList<XYSeries>();

		for (int i = 0; i < dimension; i++) {
			series.add(new XYSeries("Cluster " + (i + 1)));
		}

		for (Centroid point : centroids) {
			ArrayList<Point> cPoints = point.completePointList();
			for (Point p : cPoints) {
				series.get(p.getClusterNumb()).add(p.getCoordinates().get(0), p.getCoordinates().get(1));
			}
//			//Centroids
//			XYDataItem cent = new XYDataItem(point.getCoordinates().get(0), point.getCoordinates().get(1));
//			String key = "C";
//			String con = Integer.toString(point.getClusterNumb());
//			key+=con;
//			XYSeries ser = new XYSeries(key);
//			ser.add(cent);
//			collection.addSeries(ser);
		}

		for (int i = 0; i < dimension; i++) {
			collection.addSeries(series.get(i));
		}

		return collection;
	}

}
