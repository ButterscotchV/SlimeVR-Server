package dev.slimevr.mango;

import java.util.Arrays;
import java.util.function.Consumer;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;

import com.jme3.math.FastMath;
import com.jme3.math.Quaternion;
import com.jme3.math.Vector3f;

import dev.slimevr.VRServer;
import dev.slimevr.util.ann.VRServerThread;
import dev.slimevr.vr.processor.TransformNode;
import dev.slimevr.vr.processor.skeleton.BoneType;
import dev.slimevr.vr.processor.skeleton.HumanSkeleton;
import dev.slimevr.vr.processor.skeleton.Skeleton;
import dev.slimevr.vr.trackers.IMUTracker;
import dev.slimevr.vr.trackers.ReferenceAdjustedTracker;
import dev.slimevr.vr.trackers.Tracker;
import dev.slimevr.vr.trackers.TrackerPosition;
import io.eiren.util.collections.FastList;
import io.eiren.util.logging.LogManager;


public class MangoDataCollector {

	protected final VRServer server;
	private final FastList<MangoTrackerState> serverTrackers = new FastList<>();
	private final MangoGrid grid = new MangoGrid();

	public MangoDataCollector(VRServer server) {
		this.server = server;
		server.addNewTrackerConsumer(tracker -> {
			Tracker rawTracker = tracker;

			// Extract ReferenceAdjustedTracker trackers >:c
			while (rawTracker instanceof ReferenceAdjustedTracker) {
				rawTracker = ((ReferenceAdjustedTracker) rawTracker).tracker;
			}

			if (rawTracker instanceof IMUTracker) {
				serverTrackers
					.add(
						new MangoTrackerState(tracker, (IMUTracker) rawTracker, this::onMagUpdate)
					);
			}
		});
		server.addOnTick(this::onTick);
	}

	@VRServerThread
	public void onTick() {
		for (MangoTrackerState mangoTrackerState : serverTrackers) {
			mangoTrackerState.update();
		}
	}

	@VRServerThread
	public void onMagUpdate(MangoTrackerState magTrackerState) {
		// Get all the data needed to map this tracker's current state
		if (server.humanPoseProcessor == null)
			return;
		Skeleton untypedSkeleton = server.humanPoseProcessor.getSkeleton();

		if (!(untypedSkeleton instanceof HumanSkeleton))
			return;
		HumanSkeleton skeleton = (HumanSkeleton) untypedSkeleton;

		Quaternion trackerRotation = new Quaternion();
		if (!magTrackerState.tracker.getRotation(trackerRotation))
			return;

		BoneType trackerBoneType = trackerPositionToBoneType(
			magTrackerState.rawTracker.getBodyPosition()
		);
		if (trackerBoneType == null)
			return;

		// Node for position
		TransformNode trackerBoneTail = skeleton.getTailNodeOfBone(trackerBoneType);
		if (trackerBoneTail == null)
			return;

		// We have everything we need now!
		// Collect data and/or compare to previously collected data
		Vector3f trackerPosition = trackerBoneTail.worldTransform.getTranslation();

		if (
			magTrackerState.magAccuracy > 0f
				&& magTrackerState.magAccuracy <= 0.2f
				&& grid.dataPoints.size() < 1
		) {
			// Collect 1 data point(s)
			Quaternion offsetMag = trackerRotation.inverse().mult(magTrackerState.magRot);

			grid.addDataPoint(trackerPosition, offsetMag);

			float[] angles = offsetMag.toAngles(null);
			for (int i = 0; i < angles.length; i++) {
				angles[i] = angles[i] * FastMath.RAD_TO_DEG;
			}

			LogManager
				.debug(
					"Datapoint collected ("
						+ (grid.dataPoints.size())
						+ "/1): [Position = "
						+ trackerPosition
						+ "] [Rotation = "
						+ Arrays.toString(angles)
						+ " deg] [Magnetometer accuracy = "
						+ magTrackerState.magAccuracy * FastMath.RAD_TO_DEG
						+ " deg]"
				);
		} else {
			// Then continuously report the nearest point data
			Pair<Float, Triple<Long, Vector3f, Quaternion>> nearest = grid
				.nearestNeighbour(trackerPosition);
			if (nearest.getRight() != null) {
				float angleBetween = magTrackerState.magRot
					.angleBetween(nearest.getRight().getRight()) * FastMath.RAD_TO_DEG;

				Quaternion magOffset = magTrackerState.magRot
					.mult(nearest.getRight().getRight().inverse());
				Quaternion gyroOffset = magOffset.mult(trackerRotation);

				float[] angles = gyroOffset.toAngles(null);
				for (int i = 0; i < angles.length; i++) {
					angles[i] = angles[i] * FastMath.RAD_TO_DEG;
				}

				LogManager
					.debug(
						"{NN}: [Dist = "
							+ nearest.getLeft()
							+ " m] [Angle between = "
							+ angleBetween
							+ " deg] [Gyro offset = "
							+ Arrays.toString(angles)
							+ " deg] [Magnetometer accuracy = "
							+ magTrackerState.magAccuracy * FastMath.RAD_TO_DEG
							+ " deg]"
					);
			} else {
				LogManager.debug("No data points found...");
			}
		}
	}

	private BoneType trackerPositionToBoneType(TrackerPosition trackerPosition) {
		switch (trackerPosition) {
			case HMD:
				return BoneType.HEAD;
			case NECK:
				return BoneType.NECK;

			case LEFT_SHOULDER:
				return BoneType.LEFT_SHOULDER;
			case LEFT_UPPER_ARM:
				return BoneType.LEFT_UPPER_ARM;
			case LEFT_LOWER_ARM:
				return BoneType.LEFT_LOWER_ARM;
			case LEFT_HAND:
				return BoneType.LEFT_HAND;
			case LEFT_CONTROLLER:
				return BoneType.LEFT_CONTROLLER;

			case RIGHT_SHOULDER:
				return BoneType.RIGHT_SHOULDER;
			case RIGHT_UPPER_ARM:
				return BoneType.RIGHT_UPPER_ARM;
			case RIGHT_LOWER_ARM:
				return BoneType.RIGHT_LOWER_ARM;
			case RIGHT_HAND:
				return BoneType.RIGHT_HAND;
			case RIGHT_CONTROLLER:
				return BoneType.RIGHT_CONTROLLER;

			case CHEST:
				return BoneType.CHEST;
			case WAIST:
				return BoneType.WAIST;
			case HIP:
				return BoneType.HIP;

			case LEFT_UPPER_LEG:
				return BoneType.LEFT_UPPER_LEG;
			case LEFT_LOWER_LEG:
				return BoneType.LEFT_LOWER_LEG;
			case LEFT_FOOT:
				return BoneType.LEFT_FOOT;

			case RIGHT_UPPER_LEG:
				return BoneType.RIGHT_UPPER_LEG;
			case RIGHT_LOWER_LEG:
				return BoneType.RIGHT_LOWER_LEG;
			case RIGHT_FOOT:
				return BoneType.RIGHT_FOOT;
		}

		return null;
	}

	protected class MangoTrackerState {

		public final Tracker tracker;
		public final IMUTracker rawTracker;
		public final Consumer<MangoTrackerState> onMagUpdate;

		// State storage
		public long lastUpdateCheckTime = 0l;
		public long lastUpdateTime = 0l;

		public final Quaternion magRot = new Quaternion();
		public float magAccuracy = 0f;

		public MangoTrackerState(
			Tracker tracker,
			IMUTracker rawTracker,
			Consumer<MangoTrackerState> onMagUpdate
		) {
			this.tracker = tracker;
			this.rawTracker = rawTracker;
			this.onMagUpdate = onMagUpdate;
			update(true);
		}

		public void update(boolean forceUpdate) {
			lastUpdateCheckTime = System.currentTimeMillis();

			if (
				forceUpdate
					|| Float.compare(magAccuracy, rawTracker.magnetometerAccuracy) != 0
					|| !magRot.equals(rawTracker.rotMagQuaternion)
			) {
				// If forced to update, or any values have changed, then update
				// the state and call relevant code
				lastUpdateTime = System.currentTimeMillis();
				magRot.set(rawTracker.rotMagQuaternion);
				magAccuracy = rawTracker.magnetometerAccuracy;

				onMagUpdate.accept(this);
			}
		}

		public void update() {
			update(false);
		}
	}
}
