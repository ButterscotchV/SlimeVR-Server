package dev.slimevr.util

import com.jme3.math.FastMath
import io.github.axisangles.ktmath.Quaternion
import org.apache.commons.math3.linear.EigenDecomposition
import org.apache.commons.math3.linear.MatrixUtils
import org.apache.commons.math3.linear.RealMatrix
import kotlin.math.*

/**
 * This is a stat calculator based on Welford's online algorithm
 * https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford%27s_online_algorithm
 */
class AngularStatsCalculator {
	private var count = 0
	var mean: RealMatrix? = null
		private set
	var meanQuat = Quaternion.IDENTITY
		private set
	private var m2 = 0f

	fun reset() {
		count = 0
		mean = null
		m2 = 0f
	}

	// Using Markley's method (https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070017872.pdf)
	// Alternatively, power iteration can also be used to derive the eigenvector (https://en.wikipedia.org/wiki/Power_iteration)
	private fun toQuaternion(matrix: RealMatrix): Quaternion {
		val eigDecomp = EigenDecomposition(matrix)
		val ev0 = eigDecomp.getEigenvector(0)
		return Quaternion(
			ev0.getEntry(0).toFloat(),
			ev0.getEntry(1).toFloat(),
			ev0.getEntry(2).toFloat(),
			ev0.getEntry(3).toFloat(),
		).unit()
	}

	private fun toMatrix(quaternion: Quaternion): RealMatrix {
		val qVec = MatrixUtils.createRealVector(
			doubleArrayOf(
				quaternion.w.toDouble(),
				quaternion.x.toDouble(),
				quaternion.y.toDouble(),
				quaternion.z.toDouble(),
			),
		)
		return qVec.outerProduct(qVec)
	}

	private fun delta(newVal: Quaternion, mean: Quaternion): Float = mean.angleToR(newVal) / FastMath.TWO_PI

	fun addValue(newValue: Quaternion) {
		count += 1

		val mean = this.mean
		if (mean == null) {
			this.mean = toMatrix(newValue)
			meanQuat = newValue
			return
		}

		val delta = delta(newValue, meanQuat)
		val newMean = mean.add(toMatrix(newValue).scalarMultiply((delta / count).toDouble()))

		val newMeanQuat = toQuaternion(newMean)
		meanQuat = newMeanQuat

		val delta2 = delta(newValue, newMeanQuat)
		m2 += delta * delta2
	}

	val variance: Float
		get() = if (count < 1) {
			Float.NaN
		} else {
			m2 / count
		}
	val sampleVariance: Float
		get() = if (count < 2) {
			Float.NaN
		} else {
			m2 / (count - 1)
		}
	val standardDeviation: Float
		get() = sqrt(variance)
}
