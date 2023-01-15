package dev.slimevr.mango;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;

import com.jme3.math.FastMath;
import com.jme3.math.Quaternion;
import com.jme3.math.Vector3f;

import dev.slimevr.vr.trackers.Tracker;
import io.eiren.util.collections.FastList;


public class MangoGrid {

	Vector3f vecBuffer = new Vector3f();
	Quaternion quatBuffer = new Quaternion();

	// TODO: Make this actually use a grid or whatever
	public final float gridSpacing = 0.05f;
	public final FastList<Triple<Long, Vector3f, Quaternion>> dataPoints = new FastList<Triple<Long, Vector3f, Quaternion>>();

	public MangoGrid() {
	}

	public synchronized void addDataPoint(Vector3f position, Quaternion rotation) {
		// Lol this is so bad
		dataPoints
			.add(
				Triple
					.of(
						System.currentTimeMillis(),
						new Vector3f(position),
						new Quaternion(rotation)
					)
			);
	}

	public synchronized boolean addDataPoint(Tracker viveTracker, Tracker tracker) {
		if (viveTracker.getPosition(vecBuffer) && tracker.getRotation(quatBuffer)) {
			addDataPoint(vecBuffer, quatBuffer);
			return true;
		}

		return false;
	}

	public synchronized Pair<Float, Triple<Long, Vector3f, Quaternion>> nearestNeighbour(
		Vector3f position
	) {
		float closestDist = -1f;
		Triple<Long, Vector3f, Quaternion> closestDataPoint = null;

		for (Triple<Long, Vector3f, Quaternion> dataPoint : dataPoints) {
			float curDist = dataPoint.getMiddle().distanceSquared(position);
			if (closestDist < 0f || curDist < closestDist) {
				closestDist = curDist;
				closestDataPoint = dataPoint;
			}
		}

		return Pair
			.of(closestDist < 0f ? closestDist : FastMath.sqrt(closestDist), closestDataPoint);
	}
}
