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

package boofcv.alg.geo.selfcalib;

import boofcv.struct.calib.CameraPinhole;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.fail;

/**
 * @author Peter Abeles
 */
public class TestSelfCalibrationFourView extends CommonAutoCalibrationChecks {
	@Test
	public void stuff() {
		CameraPinhole expected = new CameraPinhole(400,410,-3,500,505,0,0);
		List<CameraPinhole> intrinsics = new ArrayList<>();
		for (int i = 0; i < 10; i++) {
			intrinsics.add(expected);
		}

		renderGood(intrinsics);
		SelfCalibrationTwoProjectives alg = new SelfCalibrationTwoProjectives();
		addProjectives(alg);

		alg.process(Q);
		fail("Implement");
	}
}