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

import org.ddogleg.optimization.functions.FunctionNtoS;
import org.ejml.data.DMatrix4x4;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.SingularOps_DDRM;
import org.ejml.equation.Equation;
import org.ejml.equation.Sequence;
import org.ejml.ops.ConvertDMatrixStruct;

/**
 * @author Peter Abeles
 */
public class SelfCalibrationTwoProjectives extends SelfCalibrationBase {

	public void process( DMatrixRMaj Q ) {
		DMatrixRMaj A0 = new DMatrixRMaj(10,10);
		DMatrixRMaj A1 = new DMatrixRMaj(10,10);

		DMatrixRMaj EQ0 = A0.createLike();
		DMatrixRMaj EQ1 = A0.createLike();
		DMatrixRMaj EQ2 = A0.createLike();
		DMatrixRMaj EQ3 = A0.createLike();
		DMatrixRMaj EQ4 = A0.createLike();
		DMatrixRMaj EQ_SUM = EQ0.createLike();

		for (int i = 1; i < projectives.size; i++) {
			Projective P0 = projectives.get(0);
			Projective P1 = projectives.get(i);

			A0.zero();A1.zero();
			I00J01(A0,convert(P0),convert(P1));
			I00J01(A1,convert(P1),convert(P0));
			CommonOps_DDRM.subtract(A0,A1,EQ0);
			CommonOps_DDRM.add(EQ_SUM,EQ0,EQ_SUM);

			A0.zero();A1.zero();
			I01J02(A0,convert(P0),convert(P1));
			I01J02(A1,convert(P1),convert(P0));
			CommonOps_DDRM.subtract(A0,A1,EQ1);
			CommonOps_DDRM.add(EQ_SUM,EQ1,EQ_SUM);

			A0.zero();A1.zero();
			I02J11(A0,convert(P0),convert(P1));
			I02J11(A1,convert(P1),convert(P0));
			CommonOps_DDRM.subtract(A0,A1,EQ2);
			CommonOps_DDRM.add(EQ_SUM,EQ2,EQ_SUM);

			A0.zero();A1.zero();
			I11J12(A0,convert(P0),convert(P1));
			I11J12(A1,convert(P1),convert(P0));
			CommonOps_DDRM.subtract(A0,A1,EQ3);
			CommonOps_DDRM.add(EQ_SUM,EQ3,EQ_SUM);

			A0.zero();A1.zero();
			I12J22(A0,convert(P0),convert(P1));
			I12J22(A1,convert(P1),convert(P0));
			CommonOps_DDRM.subtract(A0,A1,EQ4);
			CommonOps_DDRM.add(EQ_SUM,EQ4,EQ_SUM);
		}

		CommonOps_DDRM.divide(EQ_SUM,EQ_SUM.get(9,9));

		EQ_SUM.print();

		DMatrixRMaj q = new DMatrixRMaj(10,1);
		q.data[0] = Q.get(0,0);
		q.data[1] = Q.get(0,1);
		q.data[2] = Q.get(0,2);
		q.data[3] = Q.get(0,3);
		q.data[4] = Q.get(1,1);
		q.data[5] = Q.get(1,2);
		q.data[6] = Q.get(1,3);
		q.data[7] = Q.get(2,2);
		q.data[8] = Q.get(2,3);
		q.data[9] = Q.get(3,3);

		double sv[] = SingularOps_DDRM.singularValues(EQ_SUM);
		for (int i = 0; i < sv.length; i++) {
			System.out.println("sv["+i+"] = "+sv[i]);
		}

		Equation eq = new Equation(q,"q",EQ_SUM,"EQ");
		eq.print("foo = q'*EQ");
		eq.process("zero = q'*EQ*q/(q'*q)");
		System.out.println("zero = "+eq.lookupDouble("zero"));
		q.print();

//		UnconstrainedMinimization optimizer = FactoryOptimization.unconstrained();
//
//		CommonOps_DDRM.fill(q,1);
//		RandomMatrices_DDRM.fillUniform(q,-0.5,0.5,new Random((234)));
//		optimizer.setVerbose(System.out,0);
//		optimizer.setFunction(new Function(EQ_SUM),null,0);
//		optimizer.initialize(q.data,1e-8,1e-8);
//
//		System.out.println("before = "+optimizer.getFunctionValue());
//		UtilOptimize.process(optimizer,100);
//		System.out.println("after = "+optimizer.getFunctionValue());

//		q.data = optimizer.getParameters();
//		q.print();
		System.out.println("START newton's method");

		eq.alias(EQ_SUM,"Q",q,"x");
		eq.print("value = x'*Q*x");
		for (int i = 0; i < 50; i++) {
			eq.process("f = (4*(x'*Q*x))(0,0)*(Q*x)");
			eq.process("J = 8*(Q*x)*(x'*Q) + 4*(x'*Q*x)(0,0)*Q");
			eq.process("x = x - solve(J,f)");
			q.print();
			eq.print("value = x'*Q*x");
			System.out.println();
		}
	}

	static class Function implements FunctionNtoS {
		DMatrixRMaj Q;
		DMatrixRMaj x = new DMatrixRMaj(10,1);

		Equation eq = new Equation();
		Sequence sequence;
		public Function(DMatrixRMaj q) {
			Q = q;

			eq.alias(Q,"Q",x,"x");
			sequence = eq.compile("value = (x'*Q*x)/(x'*x)");
		}

		public int getNumOfInputsN() {
			return 10;
		}

		@Override
		public double process(double[] input) {
			x.data = input;
			sequence.perform();
			return eq.lookupDouble("value");
		}
	}

	private DMatrix4x4 convert( Projective P ) {
		Equation eq = new Equation(P.A,"A",P.a,"a");

		DMatrixRMaj a = eq.process("a=[A,a;[0,0,0,999]]").lookupDDRM("a");

		DMatrix4x4 M = new DMatrix4x4();
		ConvertDMatrixStruct.convert(a,M);
		return M;
	}

	public void I00J01(DMatrixRMaj A , DMatrix4x4 P , DMatrix4x4 R ) {
		/* 0 0 */A.data[ 0] = P.a11*P.a11*R.a11*R.a21;
		/* 0 1 */A.data[ 1] = 2*P.a11*R.a11*P.a12*R.a21 + P.a11*P.a11*R.a12*R.a21 + P.a11*P.a11*R.a11*R.a22;
		/* 0 2 */A.data[ 2] = 2*P.a11*R.a11*P.a13*R.a21 + P.a11*P.a11*R.a13*R.a21 + P.a11*P.a11*R.a11*R.a23;
		/* 0 3 */A.data[ 3] = 2*P.a11*R.a11*P.a14*R.a21 + P.a11*P.a11*R.a14*R.a21 + P.a11*P.a11*R.a11*R.a24;
		/* 0 4 */A.data[ 4] = R.a11*P.a12*P.a12*R.a21 + P.a11*P.a11*R.a12*R.a22;
		/* 0 5 */A.data[ 5] = 2*R.a11*P.a12*P.a13*R.a21 + P.a11*P.a11*R.a13*R.a22 + P.a11*P.a11*R.a12*R.a23;
		/* 0 6 */A.data[ 6] = 2*R.a11*P.a12*P.a14*R.a21 + P.a11*P.a11*R.a14*R.a22 + P.a11*P.a11*R.a12*R.a24;
		/* 0 7 */A.data[ 7] = R.a11*P.a13*P.a13*R.a21 + P.a11*P.a11*R.a13*R.a23;
		/* 0 8 */A.data[ 8] = 2*R.a11*P.a13*P.a14*R.a21 + P.a11*P.a11*R.a14*R.a23 + P.a11*P.a11*R.a13*R.a24;
		/* 0 9 */A.data[ 9] = R.a11*P.a14*P.a14*R.a21 + P.a11*P.a11*R.a14*R.a24;
		/* 1 1 */A.data[11] = 2*P.a11*P.a12*R.a12*R.a21 + 2*P.a11*R.a11*P.a12*R.a22;
		/* 1 2 */A.data[12] = 2*P.a11*R.a12*P.a13*R.a21 + 2*P.a11*P.a12*R.a13*R.a21 + 2*P.a11*R.a11*P.a13*R.a22 + 2*P.a11*R.a11*P.a12*R.a23;
		/* 1 3 */A.data[13] = 2*P.a11*R.a12*P.a14*R.a21 + 2*P.a11*P.a12*R.a14*R.a21 + 2*P.a11*R.a11*P.a14*R.a22 + 2*P.a11*R.a11*P.a12*R.a24;
		/* 1 4 */A.data[14] = P.a12*P.a12*R.a12*R.a21 + R.a11*P.a12*P.a12*R.a22 + 2*P.a11*P.a12*R.a12*R.a22;
		/* 1 5 */A.data[15] = 2*P.a12*R.a12*P.a13*R.a21 + 2*R.a11*P.a12*P.a13*R.a22 + 2*P.a11*P.a12*R.a13*R.a22 + 2*P.a11*P.a12*R.a12*R.a23;
		/* 1 6 */A.data[16] = 2*P.a12*R.a12*P.a14*R.a21 + 2*R.a11*P.a12*P.a14*R.a22 + 2*P.a11*P.a12*R.a14*R.a22 + 2*P.a11*P.a12*R.a12*R.a24;
		/* 1 7 */A.data[17] = R.a12*P.a13*P.a13*R.a21 + R.a11*P.a13*P.a13*R.a22 + 2*P.a11*P.a12*R.a13*R.a23;
		/* 1 8 */A.data[18] = 2*R.a12*P.a13*P.a14*R.a21 + 2*R.a11*P.a13*P.a14*R.a22 + 2*P.a11*P.a12*R.a14*R.a23 + 2*P.a11*P.a12*R.a13*R.a24;
		/* 1 9 */A.data[19] = R.a12*P.a14*P.a14*R.a21 + R.a11*P.a14*P.a14*R.a22 + 2*P.a11*P.a12*R.a14*R.a24;
		/* 2 2 */A.data[22] = 2*P.a11*P.a13*R.a13*R.a21 + 2*P.a11*R.a11*P.a13*R.a23;
		/* 2 3 */A.data[23] = 2*P.a11*R.a13*P.a14*R.a21 + 2*P.a11*P.a13*R.a14*R.a21 + 2*P.a11*R.a11*P.a14*R.a23 + 2*P.a11*R.a11*P.a13*R.a24;
		/* 2 4 */A.data[24] = P.a12*P.a12*R.a13*R.a21 + 2*P.a11*R.a12*P.a13*R.a22 + R.a11*P.a12*P.a12*R.a23;
		/* 2 5 */A.data[25] = 2*P.a12*P.a13*R.a13*R.a21 + 2*P.a11*P.a13*R.a13*R.a22 + 2*R.a11*P.a12*P.a13*R.a23 + 2*P.a11*R.a12*P.a13*R.a23;
		/* 2 6 */A.data[26] = 2*P.a12*R.a13*P.a14*R.a21 + 2*P.a11*P.a13*R.a14*R.a22 + 2*R.a11*P.a12*P.a14*R.a23 + 2*P.a11*R.a12*P.a13*R.a24;
		/* 2 7 */A.data[27] = P.a13*P.a13*R.a13*R.a21 + R.a11*P.a13*P.a13*R.a23 + 2*P.a11*P.a13*R.a13*R.a23;
		/* 2 8 */A.data[28] = 2*P.a13*R.a13*P.a14*R.a21 + 2*R.a11*P.a13*P.a14*R.a23 + 2*P.a11*P.a13*R.a14*R.a23 + 2*P.a11*P.a13*R.a13*R.a24;
		/* 2 9 */A.data[29] = R.a13*P.a14*P.a14*R.a21 + R.a11*P.a14*P.a14*R.a23 + 2*P.a11*P.a13*R.a14*R.a24;
		/* 3 3 */A.data[33] = 2*P.a11*P.a14*R.a14*R.a21 + 2*P.a11*R.a11*P.a14*R.a24;
		/* 3 4 */A.data[34] = P.a12*P.a12*R.a14*R.a21 + 2*P.a11*R.a12*P.a14*R.a22 + R.a11*P.a12*P.a12*R.a24;
		/* 3 5 */A.data[35] = 2*P.a12*P.a13*R.a14*R.a21 + 2*P.a11*R.a13*P.a14*R.a22 + 2*P.a11*R.a12*P.a14*R.a23 + 2*R.a11*P.a12*P.a13*R.a24;
		/* 3 6 */A.data[36] = 2*P.a12*P.a14*R.a14*R.a21 + 2*P.a11*P.a14*R.a14*R.a22 + 2*R.a11*P.a12*P.a14*R.a24 + 2*P.a11*R.a12*P.a14*R.a24;
		/* 3 7 */A.data[37] = P.a13*P.a13*R.a14*R.a21 + 2*P.a11*R.a13*P.a14*R.a23 + R.a11*P.a13*P.a13*R.a24;
		/* 3 8 */A.data[38] = 2*P.a13*P.a14*R.a14*R.a21 + 2*P.a11*P.a14*R.a14*R.a23 + 2*R.a11*P.a13*P.a14*R.a24 + 2*P.a11*R.a13*P.a14*R.a24;
		/* 3 9 */A.data[39] = P.a14*P.a14*R.a14*R.a21 + R.a11*P.a14*P.a14*R.a24 + 2*P.a11*P.a14*R.a14*R.a24;
		/* 4 4 */A.data[44] = P.a12*P.a12*R.a12*R.a22;
		/* 4 5 */A.data[45] = 2*P.a12*R.a12*P.a13*R.a22 + P.a12*P.a12*R.a13*R.a22 + P.a12*P.a12*R.a12*R.a23;
		/* 4 6 */A.data[46] = 2*P.a12*R.a12*P.a14*R.a22 + P.a12*P.a12*R.a14*R.a22 + P.a12*P.a12*R.a12*R.a24;
		/* 4 7 */A.data[47] = R.a12*P.a13*P.a13*R.a22 + P.a12*P.a12*R.a13*R.a23;
		/* 4 8 */A.data[48] = 2*R.a12*P.a13*P.a14*R.a22 + P.a12*P.a12*R.a14*R.a23 + P.a12*P.a12*R.a13*R.a24;
		/* 4 9 */A.data[49] = R.a12*P.a14*P.a14*R.a22 + P.a12*P.a12*R.a14*R.a24;
		/* 5 5 */A.data[55] = 2*P.a12*P.a13*R.a13*R.a22 + 2*P.a12*R.a12*P.a13*R.a23;
		/* 5 6 */A.data[56] = 2*P.a12*R.a13*P.a14*R.a22 + 2*P.a12*P.a13*R.a14*R.a22 + 2*P.a12*R.a12*P.a14*R.a23 + 2*P.a12*R.a12*P.a13*R.a24;
		/* 5 7 */A.data[57] = P.a13*P.a13*R.a13*R.a22 + R.a12*P.a13*P.a13*R.a23 + 2*P.a12*P.a13*R.a13*R.a23;
		/* 5 8 */A.data[58] = 2*P.a13*R.a13*P.a14*R.a22 + 2*R.a12*P.a13*P.a14*R.a23 + 2*P.a12*P.a13*R.a14*R.a23 + 2*P.a12*P.a13*R.a13*R.a24;
		/* 5 9 */A.data[59] = R.a13*P.a14*P.a14*R.a22 + R.a12*P.a14*P.a14*R.a23 + 2*P.a12*P.a13*R.a14*R.a24;
		/* 6 6 */A.data[66] = 2*P.a12*P.a14*R.a14*R.a22 + 2*P.a12*R.a12*P.a14*R.a24;
		/* 6 7 */A.data[67] = P.a13*P.a13*R.a14*R.a22 + 2*P.a12*R.a13*P.a14*R.a23 + R.a12*P.a13*P.a13*R.a24;
		/* 6 8 */A.data[68] = 2*P.a13*P.a14*R.a14*R.a22 + 2*P.a12*P.a14*R.a14*R.a23 + 2*R.a12*P.a13*P.a14*R.a24 + 2*P.a12*R.a13*P.a14*R.a24;
		/* 6 9 */A.data[69] = P.a14*P.a14*R.a14*R.a22 + R.a12*P.a14*P.a14*R.a24 + 2*P.a12*P.a14*R.a14*R.a24;
		/* 7 7 */A.data[77] = P.a13*P.a13*R.a13*R.a23;
		/* 7 8 */A.data[78] = 2*P.a13*R.a13*P.a14*R.a23 + P.a13*P.a13*R.a14*R.a23 + P.a13*P.a13*R.a13*R.a24;
		/* 7 9 */A.data[79] = R.a13*P.a14*P.a14*R.a23 + P.a13*P.a13*R.a14*R.a24;
		/* 8 8 */A.data[88] = 2*P.a13*P.a14*R.a14*R.a23 + 2*P.a13*R.a13*P.a14*R.a24;
		/* 8 9 */A.data[89] = P.a14*P.a14*R.a14*R.a23 + R.a13*P.a14*P.a14*R.a24 + 2*P.a13*P.a14*R.a14*R.a24;
		/* 9 9 */A.data[99] = P.a14*P.a14*R.a14*R.a24;
	}

	public void I01J02(DMatrixRMaj A , DMatrix4x4 P , DMatrix4x4 R ) {
		/* 0 0 */A.data[ 0] = P.a11*R.a11*P.a21*R.a31;
		/* 0 1 */A.data[ 1] = R.a11*P.a12*P.a21*R.a31 + P.a11*R.a12*P.a21*R.a31 + P.a11*R.a11*P.a22*R.a31 + P.a11*R.a11*P.a21*R.a32;
		/* 0 2 */A.data[ 2] = R.a11*P.a13*P.a21*R.a31 + P.a11*R.a13*P.a21*R.a31 + P.a11*R.a11*P.a23*R.a31 + P.a11*R.a11*P.a21*R.a33;
		/* 0 3 */A.data[ 3] = R.a11*P.a14*P.a21*R.a31 + P.a11*R.a14*P.a21*R.a31 + P.a11*R.a11*P.a24*R.a31 + P.a11*R.a11*P.a21*R.a34;
		/* 0 4 */A.data[ 4] = R.a11*P.a12*P.a22*R.a31 + P.a11*R.a12*P.a21*R.a32;
		/* 0 5 */A.data[ 5] = R.a11*P.a13*P.a22*R.a31 + R.a11*P.a12*P.a23*R.a31 + P.a11*R.a13*P.a21*R.a32 + P.a11*R.a12*P.a21*R.a33;
		/* 0 6 */A.data[ 6] = R.a11*P.a14*P.a22*R.a31 + R.a11*P.a12*P.a24*R.a31 + P.a11*R.a14*P.a21*R.a32 + P.a11*R.a12*P.a21*R.a34;
		/* 0 7 */A.data[ 7] = R.a11*P.a13*P.a23*R.a31 + P.a11*R.a13*P.a21*R.a33;
		/* 0 8 */A.data[ 8] = R.a11*P.a14*P.a23*R.a31 + R.a11*P.a13*P.a24*R.a31 + P.a11*R.a14*P.a21*R.a33 + P.a11*R.a13*P.a21*R.a34;
		/* 0 9 */A.data[ 9] = R.a11*P.a14*P.a24*R.a31 + P.a11*R.a14*P.a21*R.a34;
		/* 1 1 */A.data[11] = P.a12*R.a12*P.a21*R.a31 + P.a11*R.a12*P.a22*R.a31 + R.a11*P.a12*P.a21*R.a32 + P.a11*R.a11*P.a22*R.a32;
		/* 1 2 */A.data[12] = R.a12*P.a13*P.a21*R.a31 + P.a12*R.a13*P.a21*R.a31 + P.a11*R.a13*P.a22*R.a31 + P.a11*R.a12*P.a23*R.a31 + R.a11*P.a13*P.a21*R.a32 + P.a11*R.a11*P.a23*R.a32 + R.a11*P.a12*P.a21*R.a33 + P.a11*R.a11*P.a22*R.a33;
		/* 1 3 */A.data[13] = R.a12*P.a14*P.a21*R.a31 + P.a12*R.a14*P.a21*R.a31 + P.a11*R.a14*P.a22*R.a31 + P.a11*R.a12*P.a24*R.a31 + R.a11*P.a14*P.a21*R.a32 + P.a11*R.a11*P.a24*R.a32 + R.a11*P.a12*P.a21*R.a34 + P.a11*R.a11*P.a22*R.a34;
		/* 1 4 */A.data[14] = P.a12*R.a12*P.a22*R.a31 + P.a12*R.a12*P.a21*R.a32 + R.a11*P.a12*P.a22*R.a32 + P.a11*R.a12*P.a22*R.a32;
		/* 1 5 */A.data[15] = R.a12*P.a13*P.a22*R.a31 + P.a12*R.a12*P.a23*R.a31 + P.a12*R.a13*P.a21*R.a32 + R.a11*P.a13*P.a22*R.a32 + P.a11*R.a13*P.a22*R.a32 + R.a11*P.a12*P.a23*R.a32 + P.a12*R.a12*P.a21*R.a33 + P.a11*R.a12*P.a22*R.a33;
		/* 1 6 */A.data[16] = R.a12*P.a14*P.a22*R.a31 + P.a12*R.a12*P.a24*R.a31 + P.a12*R.a14*P.a21*R.a32 + R.a11*P.a14*P.a22*R.a32 + P.a11*R.a14*P.a22*R.a32 + R.a11*P.a12*P.a24*R.a32 + P.a12*R.a12*P.a21*R.a34 + P.a11*R.a12*P.a22*R.a34;
		/* 1 7 */A.data[17] = R.a12*P.a13*P.a23*R.a31 + R.a11*P.a13*P.a23*R.a32 + P.a12*R.a13*P.a21*R.a33 + P.a11*R.a13*P.a22*R.a33;
		/* 1 8 */A.data[18] = R.a12*P.a14*P.a23*R.a31 + R.a12*P.a13*P.a24*R.a31 + R.a11*P.a14*P.a23*R.a32 + R.a11*P.a13*P.a24*R.a32 + P.a12*R.a14*P.a21*R.a33 + P.a11*R.a14*P.a22*R.a33 + P.a12*R.a13*P.a21*R.a34 + P.a11*R.a13*P.a22*R.a34;
		/* 1 9 */A.data[19] = R.a12*P.a14*P.a24*R.a31 + R.a11*P.a14*P.a24*R.a32 + P.a12*R.a14*P.a21*R.a34 + P.a11*R.a14*P.a22*R.a34;
		/* 2 2 */A.data[22] = P.a13*R.a13*P.a21*R.a31 + P.a11*R.a13*P.a23*R.a31 + R.a11*P.a13*P.a21*R.a33 + P.a11*R.a11*P.a23*R.a33;
		/* 2 3 */A.data[23] = R.a13*P.a14*P.a21*R.a31 + P.a13*R.a14*P.a21*R.a31 + P.a11*R.a14*P.a23*R.a31 + P.a11*R.a13*P.a24*R.a31 + R.a11*P.a14*P.a21*R.a33 + P.a11*R.a11*P.a24*R.a33 + R.a11*P.a13*P.a21*R.a34 + P.a11*R.a11*P.a23*R.a34;
		/* 2 4 */A.data[24] = P.a12*R.a13*P.a22*R.a31 + R.a12*P.a13*P.a21*R.a32 + P.a11*R.a12*P.a23*R.a32 + R.a11*P.a12*P.a22*R.a33;
		/* 2 5 */A.data[25] = P.a13*R.a13*P.a22*R.a31 + P.a12*R.a13*P.a23*R.a31 + P.a13*R.a13*P.a21*R.a32 + P.a11*R.a13*P.a23*R.a32 + R.a12*P.a13*P.a21*R.a33 + R.a11*P.a13*P.a22*R.a33 + R.a11*P.a12*P.a23*R.a33 + P.a11*R.a12*P.a23*R.a33;
		/* 2 6 */A.data[26] = R.a13*P.a14*P.a22*R.a31 + P.a12*R.a13*P.a24*R.a31 + P.a13*R.a14*P.a21*R.a32 + P.a11*R.a14*P.a23*R.a32 + R.a11*P.a14*P.a22*R.a33 + R.a11*P.a12*P.a24*R.a33 + R.a12*P.a13*P.a21*R.a34 + P.a11*R.a12*P.a23*R.a34;
		/* 2 7 */A.data[27] = P.a13*R.a13*P.a23*R.a31 + P.a13*R.a13*P.a21*R.a33 + R.a11*P.a13*P.a23*R.a33 + P.a11*R.a13*P.a23*R.a33;
		/* 2 8 */A.data[28] = R.a13*P.a14*P.a23*R.a31 + P.a13*R.a13*P.a24*R.a31 + P.a13*R.a14*P.a21*R.a33 + R.a11*P.a14*P.a23*R.a33 + P.a11*R.a14*P.a23*R.a33 + R.a11*P.a13*P.a24*R.a33 + P.a13*R.a13*P.a21*R.a34 + P.a11*R.a13*P.a23*R.a34;
		/* 2 9 */A.data[29] = R.a13*P.a14*P.a24*R.a31 + R.a11*P.a14*P.a24*R.a33 + P.a13*R.a14*P.a21*R.a34 + P.a11*R.a14*P.a23*R.a34;
		/* 3 3 */A.data[33] = P.a14*R.a14*P.a21*R.a31 + P.a11*R.a14*P.a24*R.a31 + R.a11*P.a14*P.a21*R.a34 + P.a11*R.a11*P.a24*R.a34;
		/* 3 4 */A.data[34] = P.a12*R.a14*P.a22*R.a31 + R.a12*P.a14*P.a21*R.a32 + P.a11*R.a12*P.a24*R.a32 + R.a11*P.a12*P.a22*R.a34;
		/* 3 5 */A.data[35] = P.a13*R.a14*P.a22*R.a31 + P.a12*R.a14*P.a23*R.a31 + R.a13*P.a14*P.a21*R.a32 + P.a11*R.a13*P.a24*R.a32 + R.a12*P.a14*P.a21*R.a33 + P.a11*R.a12*P.a24*R.a33 + R.a11*P.a13*P.a22*R.a34 + R.a11*P.a12*P.a23*R.a34;
		/* 3 6 */A.data[36] = P.a14*R.a14*P.a22*R.a31 + P.a12*R.a14*P.a24*R.a31 + P.a14*R.a14*P.a21*R.a32 + P.a11*R.a14*P.a24*R.a32 + R.a12*P.a14*P.a21*R.a34 + R.a11*P.a14*P.a22*R.a34 + R.a11*P.a12*P.a24*R.a34 + P.a11*R.a12*P.a24*R.a34;
		/* 3 7 */A.data[37] = P.a13*R.a14*P.a23*R.a31 + R.a13*P.a14*P.a21*R.a33 + P.a11*R.a13*P.a24*R.a33 + R.a11*P.a13*P.a23*R.a34;
		/* 3 8 */A.data[38] = P.a14*R.a14*P.a23*R.a31 + P.a13*R.a14*P.a24*R.a31 + P.a14*R.a14*P.a21*R.a33 + P.a11*R.a14*P.a24*R.a33 + R.a13*P.a14*P.a21*R.a34 + R.a11*P.a14*P.a23*R.a34 + R.a11*P.a13*P.a24*R.a34 + P.a11*R.a13*P.a24*R.a34;
		/* 3 9 */A.data[39] = P.a14*R.a14*P.a24*R.a31 + P.a14*R.a14*P.a21*R.a34 + R.a11*P.a14*P.a24*R.a34 + P.a11*R.a14*P.a24*R.a34;
		/* 4 4 */A.data[44] = P.a12*R.a12*P.a22*R.a32;
		/* 4 5 */A.data[45] = R.a12*P.a13*P.a22*R.a32 + P.a12*R.a13*P.a22*R.a32 + P.a12*R.a12*P.a23*R.a32 + P.a12*R.a12*P.a22*R.a33;
		/* 4 6 */A.data[46] = R.a12*P.a14*P.a22*R.a32 + P.a12*R.a14*P.a22*R.a32 + P.a12*R.a12*P.a24*R.a32 + P.a12*R.a12*P.a22*R.a34;
		/* 4 7 */A.data[47] = R.a12*P.a13*P.a23*R.a32 + P.a12*R.a13*P.a22*R.a33;
		/* 4 8 */A.data[48] = R.a12*P.a14*P.a23*R.a32 + R.a12*P.a13*P.a24*R.a32 + P.a12*R.a14*P.a22*R.a33 + P.a12*R.a13*P.a22*R.a34;
		/* 4 9 */A.data[49] = R.a12*P.a14*P.a24*R.a32 + P.a12*R.a14*P.a22*R.a34;
		/* 5 5 */A.data[55] = P.a13*R.a13*P.a22*R.a32 + P.a12*R.a13*P.a23*R.a32 + R.a12*P.a13*P.a22*R.a33 + P.a12*R.a12*P.a23*R.a33;
		/* 5 6 */A.data[56] = R.a13*P.a14*P.a22*R.a32 + P.a13*R.a14*P.a22*R.a32 + P.a12*R.a14*P.a23*R.a32 + P.a12*R.a13*P.a24*R.a32 + R.a12*P.a14*P.a22*R.a33 + P.a12*R.a12*P.a24*R.a33 + R.a12*P.a13*P.a22*R.a34 + P.a12*R.a12*P.a23*R.a34;
		/* 5 7 */A.data[57] = P.a13*R.a13*P.a23*R.a32 + P.a13*R.a13*P.a22*R.a33 + R.a12*P.a13*P.a23*R.a33 + P.a12*R.a13*P.a23*R.a33;
		/* 5 8 */A.data[58] = R.a13*P.a14*P.a23*R.a32 + P.a13*R.a13*P.a24*R.a32 + P.a13*R.a14*P.a22*R.a33 + R.a12*P.a14*P.a23*R.a33 + P.a12*R.a14*P.a23*R.a33 + R.a12*P.a13*P.a24*R.a33 + P.a13*R.a13*P.a22*R.a34 + P.a12*R.a13*P.a23*R.a34;
		/* 5 9 */A.data[59] = R.a13*P.a14*P.a24*R.a32 + R.a12*P.a14*P.a24*R.a33 + P.a13*R.a14*P.a22*R.a34 + P.a12*R.a14*P.a23*R.a34;
		/* 6 6 */A.data[66] = P.a14*R.a14*P.a22*R.a32 + P.a12*R.a14*P.a24*R.a32 + R.a12*P.a14*P.a22*R.a34 + P.a12*R.a12*P.a24*R.a34;
		/* 6 7 */A.data[67] = P.a13*R.a14*P.a23*R.a32 + R.a13*P.a14*P.a22*R.a33 + P.a12*R.a13*P.a24*R.a33 + R.a12*P.a13*P.a23*R.a34;
		/* 6 8 */A.data[68] = P.a14*R.a14*P.a23*R.a32 + P.a13*R.a14*P.a24*R.a32 + P.a14*R.a14*P.a22*R.a33 + P.a12*R.a14*P.a24*R.a33 + R.a13*P.a14*P.a22*R.a34 + R.a12*P.a14*P.a23*R.a34 + R.a12*P.a13*P.a24*R.a34 + P.a12*R.a13*P.a24*R.a34;
		/* 6 9 */A.data[69] = P.a14*R.a14*P.a24*R.a32 + P.a14*R.a14*P.a22*R.a34 + R.a12*P.a14*P.a24*R.a34 + P.a12*R.a14*P.a24*R.a34;
		/* 7 7 */A.data[77] = P.a13*R.a13*P.a23*R.a33;
		/* 7 8 */A.data[78] = R.a13*P.a14*P.a23*R.a33 + P.a13*R.a14*P.a23*R.a33 + P.a13*R.a13*P.a24*R.a33 + P.a13*R.a13*P.a23*R.a34;
		/* 7 9 */A.data[79] = R.a13*P.a14*P.a24*R.a33 + P.a13*R.a14*P.a23*R.a34;
		/* 8 8 */A.data[88] = P.a14*R.a14*P.a23*R.a33 + P.a13*R.a14*P.a24*R.a33 + R.a13*P.a14*P.a23*R.a34 + P.a13*R.a13*P.a24*R.a34;
		/* 8 9 */A.data[89] = P.a14*R.a14*P.a24*R.a33 + P.a14*R.a14*P.a23*R.a34 + R.a13*P.a14*P.a24*R.a34 + P.a13*R.a14*P.a24*R.a34;
		/* 9 9 */A.data[99] = P.a14*R.a14*P.a24*R.a34;
	}

	public void I02J11(DMatrixRMaj A , DMatrix4x4 P , DMatrix4x4 R ) {
		/* 0 0 */A.data[ 0] = P.a11*R.a21*R.a21*P.a31;
		/* 0 1 */A.data[ 1] = P.a12*R.a21*R.a21*P.a31 + 2*P.a11*R.a21*R.a22*P.a31 + P.a11*R.a21*R.a21*P.a32;
		/* 0 2 */A.data[ 2] = P.a13*R.a21*R.a21*P.a31 + 2*P.a11*R.a21*R.a23*P.a31 + P.a11*R.a21*R.a21*P.a33;
		/* 0 3 */A.data[ 3] = P.a14*R.a21*R.a21*P.a31 + 2*P.a11*R.a21*R.a24*P.a31 + P.a11*R.a21*R.a21*P.a34;
		/* 0 4 */A.data[ 4] = P.a11*R.a22*R.a22*P.a31 + P.a12*R.a21*R.a21*P.a32;
		/* 0 5 */A.data[ 5] = 2*P.a11*R.a22*R.a23*P.a31 + P.a13*R.a21*R.a21*P.a32 + P.a12*R.a21*R.a21*P.a33;
		/* 0 6 */A.data[ 6] = 2*P.a11*R.a22*R.a24*P.a31 + P.a14*R.a21*R.a21*P.a32 + P.a12*R.a21*R.a21*P.a34;
		/* 0 7 */A.data[ 7] = P.a11*R.a23*R.a23*P.a31 + P.a13*R.a21*R.a21*P.a33;
		/* 0 8 */A.data[ 8] = 2*P.a11*R.a23*R.a24*P.a31 + P.a14*R.a21*R.a21*P.a33 + P.a13*R.a21*R.a21*P.a34;
		/* 0 9 */A.data[ 9] = P.a11*R.a24*R.a24*P.a31 + P.a14*R.a21*R.a21*P.a34;
		/* 1 1 */A.data[11] = 2*P.a12*R.a21*R.a22*P.a31 + 2*P.a11*R.a21*R.a22*P.a32;
		/* 1 2 */A.data[12] = 2*P.a13*R.a21*R.a22*P.a31 + 2*P.a12*R.a21*R.a23*P.a31 + 2*P.a11*R.a21*R.a23*P.a32 + 2*P.a11*R.a21*R.a22*P.a33;
		/* 1 3 */A.data[13] = 2*P.a14*R.a21*R.a22*P.a31 + 2*P.a12*R.a21*R.a24*P.a31 + 2*P.a11*R.a21*R.a24*P.a32 + 2*P.a11*R.a21*R.a22*P.a34;
		/* 1 4 */A.data[14] = P.a12*R.a22*R.a22*P.a31 + 2*P.a12*R.a21*R.a22*P.a32 + P.a11*R.a22*R.a22*P.a32;
		/* 1 5 */A.data[15] = 2*P.a12*R.a22*R.a23*P.a31 + 2*P.a13*R.a21*R.a22*P.a32 + 2*P.a11*R.a22*R.a23*P.a32 + 2*P.a12*R.a21*R.a22*P.a33;
		/* 1 6 */A.data[16] = 2*P.a12*R.a22*R.a24*P.a31 + 2*P.a14*R.a21*R.a22*P.a32 + 2*P.a11*R.a22*R.a24*P.a32 + 2*P.a12*R.a21*R.a22*P.a34;
		/* 1 7 */A.data[17] = P.a12*R.a23*R.a23*P.a31 + P.a11*R.a23*R.a23*P.a32 + 2*P.a13*R.a21*R.a22*P.a33;
		/* 1 8 */A.data[18] = 2*P.a12*R.a23*R.a24*P.a31 + 2*P.a11*R.a23*R.a24*P.a32 + 2*P.a14*R.a21*R.a22*P.a33 + 2*P.a13*R.a21*R.a22*P.a34;
		/* 1 9 */A.data[19] = P.a12*R.a24*R.a24*P.a31 + P.a11*R.a24*R.a24*P.a32 + 2*P.a14*R.a21*R.a22*P.a34;
		/* 2 2 */A.data[22] = 2*P.a13*R.a21*R.a23*P.a31 + 2*P.a11*R.a21*R.a23*P.a33;
		/* 2 3 */A.data[23] = 2*P.a14*R.a21*R.a23*P.a31 + 2*P.a13*R.a21*R.a24*P.a31 + 2*P.a11*R.a21*R.a24*P.a33 + 2*P.a11*R.a21*R.a23*P.a34;
		/* 2 4 */A.data[24] = P.a13*R.a22*R.a22*P.a31 + 2*P.a12*R.a21*R.a23*P.a32 + P.a11*R.a22*R.a22*P.a33;
		/* 2 5 */A.data[25] = 2*P.a13*R.a22*R.a23*P.a31 + 2*P.a13*R.a21*R.a23*P.a32 + 2*P.a12*R.a21*R.a23*P.a33 + 2*P.a11*R.a22*R.a23*P.a33;
		/* 2 6 */A.data[26] = 2*P.a13*R.a22*R.a24*P.a31 + 2*P.a14*R.a21*R.a23*P.a32 + 2*P.a11*R.a22*R.a24*P.a33 + 2*P.a12*R.a21*R.a23*P.a34;
		/* 2 7 */A.data[27] = P.a13*R.a23*R.a23*P.a31 + 2*P.a13*R.a21*R.a23*P.a33 + P.a11*R.a23*R.a23*P.a33;
		/* 2 8 */A.data[28] = 2*P.a13*R.a23*R.a24*P.a31 + 2*P.a14*R.a21*R.a23*P.a33 + 2*P.a11*R.a23*R.a24*P.a33 + 2*P.a13*R.a21*R.a23*P.a34;
		/* 2 9 */A.data[29] = P.a13*R.a24*R.a24*P.a31 + P.a11*R.a24*R.a24*P.a33 + 2*P.a14*R.a21*R.a23*P.a34;
		/* 3 3 */A.data[33] = 2*P.a14*R.a21*R.a24*P.a31 + 2*P.a11*R.a21*R.a24*P.a34;
		/* 3 4 */A.data[34] = P.a14*R.a22*R.a22*P.a31 + 2*P.a12*R.a21*R.a24*P.a32 + P.a11*R.a22*R.a22*P.a34;
		/* 3 5 */A.data[35] = 2*P.a14*R.a22*R.a23*P.a31 + 2*P.a13*R.a21*R.a24*P.a32 + 2*P.a12*R.a21*R.a24*P.a33 + 2*P.a11*R.a22*R.a23*P.a34;
		/* 3 6 */A.data[36] = 2*P.a14*R.a22*R.a24*P.a31 + 2*P.a14*R.a21*R.a24*P.a32 + 2*P.a12*R.a21*R.a24*P.a34 + 2*P.a11*R.a22*R.a24*P.a34;
		/* 3 7 */A.data[37] = P.a14*R.a23*R.a23*P.a31 + 2*P.a13*R.a21*R.a24*P.a33 + P.a11*R.a23*R.a23*P.a34;
		/* 3 8 */A.data[38] = 2*P.a14*R.a23*R.a24*P.a31 + 2*P.a14*R.a21*R.a24*P.a33 + 2*P.a13*R.a21*R.a24*P.a34 + 2*P.a11*R.a23*R.a24*P.a34;
		/* 3 9 */A.data[39] = P.a14*R.a24*R.a24*P.a31 + 2*P.a14*R.a21*R.a24*P.a34 + P.a11*R.a24*R.a24*P.a34;
		/* 4 4 */A.data[44] = P.a12*R.a22*R.a22*P.a32;
		/* 4 5 */A.data[45] = P.a13*R.a22*R.a22*P.a32 + 2*P.a12*R.a22*R.a23*P.a32 + P.a12*R.a22*R.a22*P.a33;
		/* 4 6 */A.data[46] = P.a14*R.a22*R.a22*P.a32 + 2*P.a12*R.a22*R.a24*P.a32 + P.a12*R.a22*R.a22*P.a34;
		/* 4 7 */A.data[47] = P.a12*R.a23*R.a23*P.a32 + P.a13*R.a22*R.a22*P.a33;
		/* 4 8 */A.data[48] = 2*P.a12*R.a23*R.a24*P.a32 + P.a14*R.a22*R.a22*P.a33 + P.a13*R.a22*R.a22*P.a34;
		/* 4 9 */A.data[49] = P.a12*R.a24*R.a24*P.a32 + P.a14*R.a22*R.a22*P.a34;
		/* 5 5 */A.data[55] = 2*P.a13*R.a22*R.a23*P.a32 + 2*P.a12*R.a22*R.a23*P.a33;
		/* 5 6 */A.data[56] = 2*P.a14*R.a22*R.a23*P.a32 + 2*P.a13*R.a22*R.a24*P.a32 + 2*P.a12*R.a22*R.a24*P.a33 + 2*P.a12*R.a22*R.a23*P.a34;
		/* 5 7 */A.data[57] = P.a13*R.a23*R.a23*P.a32 + 2*P.a13*R.a22*R.a23*P.a33 + P.a12*R.a23*R.a23*P.a33;
		/* 5 8 */A.data[58] = 2*P.a13*R.a23*R.a24*P.a32 + 2*P.a14*R.a22*R.a23*P.a33 + 2*P.a12*R.a23*R.a24*P.a33 + 2*P.a13*R.a22*R.a23*P.a34;
		/* 5 9 */A.data[59] = P.a13*R.a24*R.a24*P.a32 + P.a12*R.a24*R.a24*P.a33 + 2*P.a14*R.a22*R.a23*P.a34;
		/* 6 6 */A.data[66] = 2*P.a14*R.a22*R.a24*P.a32 + 2*P.a12*R.a22*R.a24*P.a34;
		/* 6 7 */A.data[67] = P.a14*R.a23*R.a23*P.a32 + 2*P.a13*R.a22*R.a24*P.a33 + P.a12*R.a23*R.a23*P.a34;
		/* 6 8 */A.data[68] = 2*P.a14*R.a23*R.a24*P.a32 + 2*P.a14*R.a22*R.a24*P.a33 + 2*P.a13*R.a22*R.a24*P.a34 + 2*P.a12*R.a23*R.a24*P.a34;
		/* 6 9 */A.data[69] = P.a14*R.a24*R.a24*P.a32 + 2*P.a14*R.a22*R.a24*P.a34 + P.a12*R.a24*R.a24*P.a34;
		/* 7 7 */A.data[77] = P.a13*R.a23*R.a23*P.a33;
		/* 7 8 */A.data[78] = P.a14*R.a23*R.a23*P.a33 + 2*P.a13*R.a23*R.a24*P.a33 + P.a13*R.a23*R.a23*P.a34;
		/* 7 9 */A.data[79] = P.a13*R.a24*R.a24*P.a33 + P.a14*R.a23*R.a23*P.a34;
		/* 8 8 */A.data[88] = 2*P.a14*R.a23*R.a24*P.a33 + 2*P.a13*R.a23*R.a24*P.a34;
		/* 8 9 */A.data[89] = P.a14*R.a24*R.a24*P.a33 + 2*P.a14*R.a23*R.a24*P.a34 + P.a13*R.a24*R.a24*P.a34;
		/* 9 9 */A.data[99] = P.a14*R.a24*R.a24*P.a34;
	}

	public void I11J12(DMatrixRMaj A , DMatrix4x4 P , DMatrix4x4 R ) {
		/* 0 0 */A.data[ 0] = P.a21*P.a21*R.a21*R.a31;
		/* 0 1 */A.data[ 1] = 2*P.a21*R.a21*P.a22*R.a31 + P.a21*P.a21*R.a22*R.a31 + P.a21*P.a21*R.a21*R.a32;
		/* 0 2 */A.data[ 2] = 2*P.a21*R.a21*P.a23*R.a31 + P.a21*P.a21*R.a23*R.a31 + P.a21*P.a21*R.a21*R.a33;
		/* 0 3 */A.data[ 3] = 2*P.a21*R.a21*P.a24*R.a31 + P.a21*P.a21*R.a24*R.a31 + P.a21*P.a21*R.a21*R.a34;
		/* 0 4 */A.data[ 4] = R.a21*P.a22*P.a22*R.a31 + P.a21*P.a21*R.a22*R.a32;
		/* 0 5 */A.data[ 5] = 2*R.a21*P.a22*P.a23*R.a31 + P.a21*P.a21*R.a23*R.a32 + P.a21*P.a21*R.a22*R.a33;
		/* 0 6 */A.data[ 6] = 2*R.a21*P.a22*P.a24*R.a31 + P.a21*P.a21*R.a24*R.a32 + P.a21*P.a21*R.a22*R.a34;
		/* 0 7 */A.data[ 7] = R.a21*P.a23*P.a23*R.a31 + P.a21*P.a21*R.a23*R.a33;
		/* 0 8 */A.data[ 8] = 2*R.a21*P.a23*P.a24*R.a31 + P.a21*P.a21*R.a24*R.a33 + P.a21*P.a21*R.a23*R.a34;
		/* 0 9 */A.data[ 9] = R.a21*P.a24*P.a24*R.a31 + P.a21*P.a21*R.a24*R.a34;
		/* 1 1 */A.data[11] = 2*P.a21*P.a22*R.a22*R.a31 + 2*P.a21*R.a21*P.a22*R.a32;
		/* 1 2 */A.data[12] = 2*P.a21*R.a22*P.a23*R.a31 + 2*P.a21*P.a22*R.a23*R.a31 + 2*P.a21*R.a21*P.a23*R.a32 + 2*P.a21*R.a21*P.a22*R.a33;
		/* 1 3 */A.data[13] = 2*P.a21*R.a22*P.a24*R.a31 + 2*P.a21*P.a22*R.a24*R.a31 + 2*P.a21*R.a21*P.a24*R.a32 + 2*P.a21*R.a21*P.a22*R.a34;
		/* 1 4 */A.data[14] = P.a22*P.a22*R.a22*R.a31 + R.a21*P.a22*P.a22*R.a32 + 2*P.a21*P.a22*R.a22*R.a32;
		/* 1 5 */A.data[15] = 2*P.a22*R.a22*P.a23*R.a31 + 2*R.a21*P.a22*P.a23*R.a32 + 2*P.a21*P.a22*R.a23*R.a32 + 2*P.a21*P.a22*R.a22*R.a33;
		/* 1 6 */A.data[16] = 2*P.a22*R.a22*P.a24*R.a31 + 2*R.a21*P.a22*P.a24*R.a32 + 2*P.a21*P.a22*R.a24*R.a32 + 2*P.a21*P.a22*R.a22*R.a34;
		/* 1 7 */A.data[17] = R.a22*P.a23*P.a23*R.a31 + R.a21*P.a23*P.a23*R.a32 + 2*P.a21*P.a22*R.a23*R.a33;
		/* 1 8 */A.data[18] = 2*R.a22*P.a23*P.a24*R.a31 + 2*R.a21*P.a23*P.a24*R.a32 + 2*P.a21*P.a22*R.a24*R.a33 + 2*P.a21*P.a22*R.a23*R.a34;
		/* 1 9 */A.data[19] = R.a22*P.a24*P.a24*R.a31 + R.a21*P.a24*P.a24*R.a32 + 2*P.a21*P.a22*R.a24*R.a34;
		/* 2 2 */A.data[22] = 2*P.a21*P.a23*R.a23*R.a31 + 2*P.a21*R.a21*P.a23*R.a33;
		/* 2 3 */A.data[23] = 2*P.a21*R.a23*P.a24*R.a31 + 2*P.a21*P.a23*R.a24*R.a31 + 2*P.a21*R.a21*P.a24*R.a33 + 2*P.a21*R.a21*P.a23*R.a34;
		/* 2 4 */A.data[24] = P.a22*P.a22*R.a23*R.a31 + 2*P.a21*R.a22*P.a23*R.a32 + R.a21*P.a22*P.a22*R.a33;
		/* 2 5 */A.data[25] = 2*P.a22*P.a23*R.a23*R.a31 + 2*P.a21*P.a23*R.a23*R.a32 + 2*R.a21*P.a22*P.a23*R.a33 + 2*P.a21*R.a22*P.a23*R.a33;
		/* 2 6 */A.data[26] = 2*P.a22*R.a23*P.a24*R.a31 + 2*P.a21*P.a23*R.a24*R.a32 + 2*R.a21*P.a22*P.a24*R.a33 + 2*P.a21*R.a22*P.a23*R.a34;
		/* 2 7 */A.data[27] = P.a23*P.a23*R.a23*R.a31 + R.a21*P.a23*P.a23*R.a33 + 2*P.a21*P.a23*R.a23*R.a33;
		/* 2 8 */A.data[28] = 2*P.a23*R.a23*P.a24*R.a31 + 2*R.a21*P.a23*P.a24*R.a33 + 2*P.a21*P.a23*R.a24*R.a33 + 2*P.a21*P.a23*R.a23*R.a34;
		/* 2 9 */A.data[29] = R.a23*P.a24*P.a24*R.a31 + R.a21*P.a24*P.a24*R.a33 + 2*P.a21*P.a23*R.a24*R.a34;
		/* 3 3 */A.data[33] = 2*P.a21*P.a24*R.a24*R.a31 + 2*P.a21*R.a21*P.a24*R.a34;
		/* 3 4 */A.data[34] = P.a22*P.a22*R.a24*R.a31 + 2*P.a21*R.a22*P.a24*R.a32 + R.a21*P.a22*P.a22*R.a34;
		/* 3 5 */A.data[35] = 2*P.a22*P.a23*R.a24*R.a31 + 2*P.a21*R.a23*P.a24*R.a32 + 2*P.a21*R.a22*P.a24*R.a33 + 2*R.a21*P.a22*P.a23*R.a34;
		/* 3 6 */A.data[36] = 2*P.a22*P.a24*R.a24*R.a31 + 2*P.a21*P.a24*R.a24*R.a32 + 2*R.a21*P.a22*P.a24*R.a34 + 2*P.a21*R.a22*P.a24*R.a34;
		/* 3 7 */A.data[37] = P.a23*P.a23*R.a24*R.a31 + 2*P.a21*R.a23*P.a24*R.a33 + R.a21*P.a23*P.a23*R.a34;
		/* 3 8 */A.data[38] = 2*P.a23*P.a24*R.a24*R.a31 + 2*P.a21*P.a24*R.a24*R.a33 + 2*R.a21*P.a23*P.a24*R.a34 + 2*P.a21*R.a23*P.a24*R.a34;
		/* 3 9 */A.data[39] = P.a24*P.a24*R.a24*R.a31 + R.a21*P.a24*P.a24*R.a34 + 2*P.a21*P.a24*R.a24*R.a34;
		/* 4 4 */A.data[44] = P.a22*P.a22*R.a22*R.a32;
		/* 4 5 */A.data[45] = 2*P.a22*R.a22*P.a23*R.a32 + P.a22*P.a22*R.a23*R.a32 + P.a22*P.a22*R.a22*R.a33;
		/* 4 6 */A.data[46] = 2*P.a22*R.a22*P.a24*R.a32 + P.a22*P.a22*R.a24*R.a32 + P.a22*P.a22*R.a22*R.a34;
		/* 4 7 */A.data[47] = R.a22*P.a23*P.a23*R.a32 + P.a22*P.a22*R.a23*R.a33;
		/* 4 8 */A.data[48] = 2*R.a22*P.a23*P.a24*R.a32 + P.a22*P.a22*R.a24*R.a33 + P.a22*P.a22*R.a23*R.a34;
		/* 4 9 */A.data[49] = R.a22*P.a24*P.a24*R.a32 + P.a22*P.a22*R.a24*R.a34;
		/* 5 5 */A.data[55] = 2*P.a22*P.a23*R.a23*R.a32 + 2*P.a22*R.a22*P.a23*R.a33;
		/* 5 6 */A.data[56] = 2*P.a22*R.a23*P.a24*R.a32 + 2*P.a22*P.a23*R.a24*R.a32 + 2*P.a22*R.a22*P.a24*R.a33 + 2*P.a22*R.a22*P.a23*R.a34;
		/* 5 7 */A.data[57] = P.a23*P.a23*R.a23*R.a32 + R.a22*P.a23*P.a23*R.a33 + 2*P.a22*P.a23*R.a23*R.a33;
		/* 5 8 */A.data[58] = 2*P.a23*R.a23*P.a24*R.a32 + 2*R.a22*P.a23*P.a24*R.a33 + 2*P.a22*P.a23*R.a24*R.a33 + 2*P.a22*P.a23*R.a23*R.a34;
		/* 5 9 */A.data[59] = R.a23*P.a24*P.a24*R.a32 + R.a22*P.a24*P.a24*R.a33 + 2*P.a22*P.a23*R.a24*R.a34;
		/* 6 6 */A.data[66] = 2*P.a22*P.a24*R.a24*R.a32 + 2*P.a22*R.a22*P.a24*R.a34;
		/* 6 7 */A.data[67] = P.a23*P.a23*R.a24*R.a32 + 2*P.a22*R.a23*P.a24*R.a33 + R.a22*P.a23*P.a23*R.a34;
		/* 6 8 */A.data[68] = 2*P.a23*P.a24*R.a24*R.a32 + 2*P.a22*P.a24*R.a24*R.a33 + 2*R.a22*P.a23*P.a24*R.a34 + 2*P.a22*R.a23*P.a24*R.a34;
		/* 6 9 */A.data[69] = P.a24*P.a24*R.a24*R.a32 + R.a22*P.a24*P.a24*R.a34 + 2*P.a22*P.a24*R.a24*R.a34;
		/* 7 7 */A.data[77] = P.a23*P.a23*R.a23*R.a33;
		/* 7 8 */A.data[78] = 2*P.a23*R.a23*P.a24*R.a33 + P.a23*P.a23*R.a24*R.a33 + P.a23*P.a23*R.a23*R.a34;
		/* 7 9 */A.data[79] = R.a23*P.a24*P.a24*R.a33 + P.a23*P.a23*R.a24*R.a34;
		/* 8 8 */A.data[88] = 2*P.a23*P.a24*R.a24*R.a33 + 2*P.a23*R.a23*P.a24*R.a34;
		/* 8 9 */A.data[89] = P.a24*P.a24*R.a24*R.a33 + R.a23*P.a24*P.a24*R.a34 + 2*P.a23*P.a24*R.a24*R.a34;
		/* 9 9 */A.data[99] = P.a24*P.a24*R.a24*R.a34;
	}

	public void I12J22(DMatrixRMaj A , DMatrix4x4 P , DMatrix4x4 R ) {
		/* 0 0 */A.data[ 0] = P.a21*P.a31*R.a31*R.a31;
		/* 0 1 */A.data[ 1] = P.a22*P.a31*R.a31*R.a31 + P.a21*R.a31*R.a31*P.a32 + 2*P.a21*P.a31*R.a31*R.a32;
		/* 0 2 */A.data[ 2] = P.a23*P.a31*R.a31*R.a31 + P.a21*R.a31*R.a31*P.a33 + 2*P.a21*P.a31*R.a31*R.a33;
		/* 0 3 */A.data[ 3] = P.a24*P.a31*R.a31*R.a31 + P.a21*R.a31*R.a31*P.a34 + 2*P.a21*P.a31*R.a31*R.a34;
		/* 0 4 */A.data[ 4] = P.a22*R.a31*R.a31*P.a32 + P.a21*P.a31*R.a32*R.a32;
		/* 0 5 */A.data[ 5] = P.a23*R.a31*R.a31*P.a32 + P.a22*R.a31*R.a31*P.a33 + 2*P.a21*P.a31*R.a32*R.a33;
		/* 0 6 */A.data[ 6] = P.a24*R.a31*R.a31*P.a32 + P.a22*R.a31*R.a31*P.a34 + 2*P.a21*P.a31*R.a32*R.a34;
		/* 0 7 */A.data[ 7] = P.a23*R.a31*R.a31*P.a33 + P.a21*P.a31*R.a33*R.a33;
		/* 0 8 */A.data[ 8] = P.a24*R.a31*R.a31*P.a33 + P.a23*R.a31*R.a31*P.a34 + 2*P.a21*P.a31*R.a33*R.a34;
		/* 0 9 */A.data[ 9] = P.a24*R.a31*R.a31*P.a34 + P.a21*P.a31*R.a34*R.a34;
		/* 1 1 */A.data[11] = 2*P.a22*P.a31*R.a31*R.a32 + 2*P.a21*R.a31*P.a32*R.a32;
		/* 1 2 */A.data[12] = 2*P.a23*P.a31*R.a31*R.a32 + 2*P.a21*R.a31*R.a32*P.a33 + 2*P.a22*P.a31*R.a31*R.a33 + 2*P.a21*R.a31*P.a32*R.a33;
		/* 1 3 */A.data[13] = 2*P.a24*P.a31*R.a31*R.a32 + 2*P.a21*R.a31*R.a32*P.a34 + 2*P.a22*P.a31*R.a31*R.a34 + 2*P.a21*R.a31*P.a32*R.a34;
		/* 1 4 */A.data[14] = 2*P.a22*R.a31*P.a32*R.a32 + P.a22*P.a31*R.a32*R.a32 + P.a21*P.a32*R.a32*R.a32;
		/* 1 5 */A.data[15] = 2*P.a23*R.a31*P.a32*R.a32 + 2*P.a22*R.a31*R.a32*P.a33 + 2*P.a22*P.a31*R.a32*R.a33 + 2*P.a21*P.a32*R.a32*R.a33;
		/* 1 6 */A.data[16] = 2*P.a24*R.a31*P.a32*R.a32 + 2*P.a22*R.a31*R.a32*P.a34 + 2*P.a22*P.a31*R.a32*R.a34 + 2*P.a21*P.a32*R.a32*R.a34;
		/* 1 7 */A.data[17] = 2*P.a23*R.a31*R.a32*P.a33 + P.a22*P.a31*R.a33*R.a33 + P.a21*P.a32*R.a33*R.a33;
		/* 1 8 */A.data[18] = 2*P.a24*R.a31*R.a32*P.a33 + 2*P.a23*R.a31*R.a32*P.a34 + 2*P.a22*P.a31*R.a33*R.a34 + 2*P.a21*P.a32*R.a33*R.a34;
		/* 1 9 */A.data[19] = 2*P.a24*R.a31*R.a32*P.a34 + P.a22*P.a31*R.a34*R.a34 + P.a21*P.a32*R.a34*R.a34;
		/* 2 2 */A.data[22] = 2*P.a23*P.a31*R.a31*R.a33 + 2*P.a21*R.a31*P.a33*R.a33;
		/* 2 3 */A.data[23] = 2*P.a24*P.a31*R.a31*R.a33 + 2*P.a21*R.a31*R.a33*P.a34 + 2*P.a23*P.a31*R.a31*R.a34 + 2*P.a21*R.a31*P.a33*R.a34;
		/* 2 4 */A.data[24] = P.a23*P.a31*R.a32*R.a32 + P.a21*R.a32*R.a32*P.a33 + 2*P.a22*R.a31*P.a32*R.a33;
		/* 2 5 */A.data[25] = 2*P.a23*R.a31*P.a32*R.a33 + 2*P.a23*P.a31*R.a32*R.a33 + 2*P.a22*R.a31*P.a33*R.a33 + 2*P.a21*R.a32*P.a33*R.a33;
		/* 2 6 */A.data[26] = 2*P.a24*R.a31*P.a32*R.a33 + 2*P.a22*R.a31*R.a33*P.a34 + 2*P.a23*P.a31*R.a32*R.a34 + 2*P.a21*R.a32*P.a33*R.a34;
		/* 2 7 */A.data[27] = 2*P.a23*R.a31*P.a33*R.a33 + P.a23*P.a31*R.a33*R.a33 + P.a21*P.a33*R.a33*R.a33;
		/* 2 8 */A.data[28] = 2*P.a24*R.a31*P.a33*R.a33 + 2*P.a23*R.a31*R.a33*P.a34 + 2*P.a23*P.a31*R.a33*R.a34 + 2*P.a21*P.a33*R.a33*R.a34;
		/* 2 9 */A.data[29] = 2*P.a24*R.a31*R.a33*P.a34 + P.a23*P.a31*R.a34*R.a34 + P.a21*P.a33*R.a34*R.a34;
		/* 3 3 */A.data[33] = 2*P.a24*P.a31*R.a31*R.a34 + 2*P.a21*R.a31*P.a34*R.a34;
		/* 3 4 */A.data[34] = P.a24*P.a31*R.a32*R.a32 + P.a21*R.a32*R.a32*P.a34 + 2*P.a22*R.a31*P.a32*R.a34;
		/* 3 5 */A.data[35] = 2*P.a24*P.a31*R.a32*R.a33 + 2*P.a21*R.a32*R.a33*P.a34 + 2*P.a23*R.a31*P.a32*R.a34 + 2*P.a22*R.a31*P.a33*R.a34;
		/* 3 6 */A.data[36] = 2*P.a24*R.a31*P.a32*R.a34 + 2*P.a24*P.a31*R.a32*R.a34 + 2*P.a22*R.a31*P.a34*R.a34 + 2*P.a21*R.a32*P.a34*R.a34;
		/* 3 7 */A.data[37] = P.a24*P.a31*R.a33*R.a33 + P.a21*R.a33*R.a33*P.a34 + 2*P.a23*R.a31*P.a33*R.a34;
		/* 3 8 */A.data[38] = 2*P.a24*R.a31*P.a33*R.a34 + 2*P.a24*P.a31*R.a33*R.a34 + 2*P.a23*R.a31*P.a34*R.a34 + 2*P.a21*R.a33*P.a34*R.a34;
		/* 3 9 */A.data[39] = 2*P.a24*R.a31*P.a34*R.a34 + P.a24*P.a31*R.a34*R.a34 + P.a21*P.a34*R.a34*R.a34;
		/* 4 4 */A.data[44] = P.a22*P.a32*R.a32*R.a32;
		/* 4 5 */A.data[45] = P.a23*P.a32*R.a32*R.a32 + P.a22*R.a32*R.a32*P.a33 + 2*P.a22*P.a32*R.a32*R.a33;
		/* 4 6 */A.data[46] = P.a24*P.a32*R.a32*R.a32 + P.a22*R.a32*R.a32*P.a34 + 2*P.a22*P.a32*R.a32*R.a34;
		/* 4 7 */A.data[47] = P.a23*R.a32*R.a32*P.a33 + P.a22*P.a32*R.a33*R.a33;
		/* 4 8 */A.data[48] = P.a24*R.a32*R.a32*P.a33 + P.a23*R.a32*R.a32*P.a34 + 2*P.a22*P.a32*R.a33*R.a34;
		/* 4 9 */A.data[49] = P.a24*R.a32*R.a32*P.a34 + P.a22*P.a32*R.a34*R.a34;
		/* 5 5 */A.data[55] = 2*P.a23*P.a32*R.a32*R.a33 + 2*P.a22*R.a32*P.a33*R.a33;
		/* 5 6 */A.data[56] = 2*P.a24*P.a32*R.a32*R.a33 + 2*P.a22*R.a32*R.a33*P.a34 + 2*P.a23*P.a32*R.a32*R.a34 + 2*P.a22*R.a32*P.a33*R.a34;
		/* 5 7 */A.data[57] = 2*P.a23*R.a32*P.a33*R.a33 + P.a23*P.a32*R.a33*R.a33 + P.a22*P.a33*R.a33*R.a33;
		/* 5 8 */A.data[58] = 2*P.a24*R.a32*P.a33*R.a33 + 2*P.a23*R.a32*R.a33*P.a34 + 2*P.a23*P.a32*R.a33*R.a34 + 2*P.a22*P.a33*R.a33*R.a34;
		/* 5 9 */A.data[59] = 2*P.a24*R.a32*R.a33*P.a34 + P.a23*P.a32*R.a34*R.a34 + P.a22*P.a33*R.a34*R.a34;
		/* 6 6 */A.data[66] = 2*P.a24*P.a32*R.a32*R.a34 + 2*P.a22*R.a32*P.a34*R.a34;
		/* 6 7 */A.data[67] = P.a24*P.a32*R.a33*R.a33 + P.a22*R.a33*R.a33*P.a34 + 2*P.a23*R.a32*P.a33*R.a34;
		/* 6 8 */A.data[68] = 2*P.a24*R.a32*P.a33*R.a34 + 2*P.a24*P.a32*R.a33*R.a34 + 2*P.a23*R.a32*P.a34*R.a34 + 2*P.a22*R.a33*P.a34*R.a34;
		/* 6 9 */A.data[69] = 2*P.a24*R.a32*P.a34*R.a34 + P.a24*P.a32*R.a34*R.a34 + P.a22*P.a34*R.a34*R.a34;
		/* 7 7 */A.data[77] = P.a23*P.a33*R.a33*R.a33;
		/* 7 8 */A.data[78] = P.a24*P.a33*R.a33*R.a33 + P.a23*R.a33*R.a33*P.a34 + 2*P.a23*P.a33*R.a33*R.a34;
		/* 7 9 */A.data[79] = P.a24*R.a33*R.a33*P.a34 + P.a23*P.a33*R.a34*R.a34;
		/* 8 8 */A.data[88] = 2*P.a24*P.a33*R.a33*R.a34 + 2*P.a23*R.a33*P.a34*R.a34;
		/* 8 9 */A.data[89] = 2*P.a24*R.a33*P.a34*R.a34 + P.a24*P.a33*R.a34*R.a34 + P.a23*P.a34*R.a34*R.a34;
		/* 9 9 */A.data[99] = P.a24*P.a34*R.a34*R.a34;
	}

}
