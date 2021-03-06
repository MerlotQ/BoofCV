/*
 * Copyright (c) 2011-2018, Peter Abeles. All Rights Reserved.
 *
 * This file is part of BoofCV (http://boofcv.org).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package boofcv.alg.geo.calibration;

import georegression.geometry.ConvertRotation3D_F64;
import georegression.struct.point.Point2D_F64;
import georegression.struct.se.Se3_F64;
import org.ddogleg.optimization.FactoryOptimization;
import org.ddogleg.optimization.UnconstrainedLeastSquares;
import org.ejml.data.DMatrixRMaj;

import java.util.ArrayList;
import java.util.List;

/**
 * <p>
 * Full implementation of the Zhang99 camera calibration algorithm using planar calibration targets.  First
 * linear approximations of camera parameters are computed, which are then refined using non-linear estimation.
 * One difference from the original paper is that tangential distortion can be included. No linear estimate
 * if found for tangential, they are estimated by initializing the non-linear estimate with all zero.
 * </p>
 *
 * <p>
 * When processing the results be sure to take in account the coordinate system being left or right handed.  Calibration
 * works just fine with either coordinate system, but most 3D geometric algorithms assume a right handed coordinate
 * system while most images are left handed.
 * </p>
 *
 * <p>
 * A listener can be provide that will give status updates and allows requests for early termination.  If a request
 * for early termination is made then a RuntimeException will be thrown.
 * </p>
 *
 * <p>
 * [1] Zhengyou Zhang, "Flexible Camera Calibration By Viewing a Plane From Unknown Orientations,",
 * International Conference on Computer Vision (ICCV'99), Corfu, Greece, pages 666-673, September 1999.
 * </p>
 *
 * @author Peter Abeles
 */
public class CalibrationPlanarGridZhang99 {

	// estimation algorithms
	private Zhang99ComputeTargetHomography computeHomography;
	private Zhang99CalibrationMatrixFromHomographies computeK;
	private RadialDistortionEstimateLinear computeRadial;
	private Zhang99DecomposeHomography decomposeH = new Zhang99DecomposeHomography();

	// contains found parameters
	private Zhang99AllParam initial;
	private Zhang99AllParam optimized;

	// optimization algorithm
	private UnconstrainedLeastSquares optimizer;

	// provides information on calibration status
	private Listener listener;

	// where calibration points are layout on the target.
	private List<Point2D_F64> layout;

	/**
	 * Configures calibration process.
	 *
	 * @param layout Layout of calibration points on the target
	 * @param intrinsicParam
	 */
	public CalibrationPlanarGridZhang99(List<Point2D_F64> layout, Zhang99IntrinsicParam intrinsicParam)
	{
		this.layout = layout;
		computeHomography = new Zhang99ComputeTargetHomography(layout);
		computeK = new Zhang99CalibrationMatrixFromHomographies(intrinsicParam.assumeZeroSkew);
		computeRadial = new RadialDistortionEstimateLinear(layout,intrinsicParam.getNumberOfRadial());
		optimized = new Zhang99AllParam(intrinsicParam,0);
		initial = optimized.createLike();
	}

	/**
	 * Used to listen in on progress and request that processing be stopped
	 *
	 * @param listener The listener
	 */
	public void setListener(Listener listener) {
		this.listener = listener;
	}

	/**
	 * Processes observed calibration point coordinates and computes camera intrinsic and extrinsic
	 * parameters.
	 *
	 * @param observations Set of observed grid locations in pixel coordinates.
	 * @return true if successful and false if it failed
	 */
	public boolean process( List<CalibrationObservation> observations ) {

		// compute initial parameter estimates using linear algebra
		if( !linearEstimate(observations,initial) )
			return false;

		status("Non-linear refinement");
		// perform non-linear optimization to improve results
		optimized.setNumberOfViews(observations.size());
		if( !optimizedParam(observations,layout,initial,optimized,optimizer))
			return false;

		return true;
	}

	/**
	 * Find an initial estimate for calibration parameters using linear techniques.
	 */
	protected boolean linearEstimate(List<CalibrationObservation> observations , Zhang99AllParam param )
	{
		status("Estimating Homographies");
		List<DMatrixRMaj> homographies = new ArrayList<>();
		List<Se3_F64> motions = new ArrayList<>();

		for( CalibrationObservation obs : observations ) {
			if( !computeHomography.computeHomography(obs) )
				return false;

			DMatrixRMaj H = computeHomography.getHomography();

			homographies.add(H);
		}

		status("Estimating Calibration Matrix");
		computeK.process(homographies);

		DMatrixRMaj K = computeK.getCalibrationMatrix();

		decomposeH.setCalibrationMatrix(K);
		for( DMatrixRMaj H : homographies ) {
			motions.add(decomposeH.decompose(H));
		}

		status("Estimating Radial Distortion");
		computeRadial.process(K, homographies, observations);

		double distort[] = computeRadial.getParameters();

		convertIntoZhangParam(motions, K,distort, param);
		return true;
	}

	private void status( String message ) {
		if( listener != null ) {
			if( !listener.zhangUpdate(message) )
				throw new RuntimeException("User requested termination of calibration");
		}
	}

	/**
	 * Use non-linear optimization to improve the parameter estimates
	 *
	 * @param observations Observations of calibration points in each image
	 * @param grid Location of calibration points on calibration target
	 * @param initial Initial estimate of calibration parameters.
	 * @param found The refined calibration parameters.
	 * @param optimizer Algorithm used to optimize parameters
	 */
	public boolean optimizedParam( List<CalibrationObservation> observations ,
								   List<Point2D_F64> grid ,
								   Zhang99AllParam initial ,
								   Zhang99AllParam found ,
								   UnconstrainedLeastSquares optimizer )
	{
		if( optimizer == null ) {
			optimizer = FactoryOptimization.levenbergMarquardt(null,true);
		}

		double model[] = new double[ initial.numParameters() ];
		initial.convertToParam(model);

		Zhang99OptimizationFunction func = new Zhang99OptimizationFunction(
				initial.createLike(), grid,observations);

		Zhang99OptimizationJacobian jacobian = initial.getIntrinsic().createJacobian(observations,grid);

		optimizer.setFunction(func,jacobian);
		optimizer.initialize(model,1e-10,1e-25*observations.size());

//		System.out.println("Error before = "+optimizer.getFunctionValue());

		for( int i = 0; i < 500; i++ ) {
//			System.out.println("i = "+i);
			if( optimizer.iterate() ) {
				break;
			} else {
				if( i % 25 == 0 )
					status("Progress "+(100*i/500.0)+"%");
			}
		}

		double param[] = optimizer.getParameters();
		found.setFromParam(param);

//		System.out.println("Error after = "+optimizer.getFunctionValue());

		return true;
	}

	/**
	 * Converts results fond in the linear algorithms into {@link Zhang99AllParam}
	 */
	public void convertIntoZhangParam(List<Se3_F64> motions,
									  DMatrixRMaj K,
									  double[] distort ,
									  Zhang99AllParam param ) {
		param.getIntrinsic().initialize(K,distort);

		param.setNumberOfViews(motions.size());

		for( int i = 0; i < param.views.length; i++ ) {
			Se3_F64 m = motions.get(i);

			Zhang99AllParam.View v = new Zhang99AllParam.View();
			v.T = m.getT();
			ConvertRotation3D_F64.matrixToRodrigues(m.getR(), v.rotation);

			param.views[i] = v;
		}
	}

	/**
	 * Applies radial and tangential distortion to the normalized image coordinate.
	 *
	 * @param normPt point in normalized image coordinates
	 * @param radial radial distortion parameters
	 * @param t1 tangential parameter
	 * @param t2 tangential parameter
	 */
	public static void applyDistortion(Point2D_F64 normPt, double[] radial, double t1 , double t2 )
	{
		final double x = normPt.x;
		final double y = normPt.y;

		double a = 0;
		double r2 = x*x + y*y;
		double r2i = r2;
		for( int i = 0; i < radial.length; i++ ) {
			a += radial[i]*r2i;
			r2i *= r2;
		}

		normPt.x = x + x*a + 2*t1*x*y + t2*(r2 + 2*x*x);
		normPt.y = y + y*a + t1*(r2 + 2*y*y) + 2*t2*x*y;
	}

	/**
	 * Specify which optimization algorithm to use
	 */
	public void setOptimizer(UnconstrainedLeastSquares optimizer) {
		this.optimizer = optimizer;
	}

	public Zhang99AllParam getOptimized() {
		return optimized;
	}

	public static int totalPoints( List<CalibrationObservation> observations ) {
		int total = 0;
		for (int i = 0; i < observations.size(); i++) {
			total += observations.get(i).size();
		}
		return total;
	}

	public interface Listener
	{
		/**
		 * Updated to update the status and request that processing be stopped
		 *
		 * @param taskName Name of the task being performed
		 * @return true to continue and false to request a stop
		 */
		public boolean zhangUpdate( String taskName );
	}
}
