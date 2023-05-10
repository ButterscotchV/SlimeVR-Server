package dev.slimevr.autobone.errors

import dev.slimevr.autobone.AutoBoneStep
import dev.slimevr.poseframeformat.trackerdata.TrackerFrames
import dev.slimevr.tracking.processor.skeleton.HumanSkeleton

// The distance of any points to the corresponding absolute position
class PositionError : IAutoBoneError {
	@Throws(AutoBoneException::class)
	override fun getStepError(trainingStep: AutoBoneStep): Float {
		val trackers = trainingStep.trainingFrames.frameHolders
		return (
			(
				getPositionError(
					trackers,
					trainingStep.cursor1,
					trainingStep.skeleton1.skeleton,
					trainingStep.skeletonScale
				) +
					getPositionError(
						trackers,
						trainingStep.cursor2,
						trainingStep.skeleton2.skeleton,
						trainingStep.skeletonScale
					)
				) /
				2f
			)
	}

	companion object {
		fun getPositionError(
			trackers: List<TrackerFrames>,
			cursor: Int,
			skeleton: HumanSkeleton,
			scale: Float = 1f,
		): Float {
			var offset = 0f
			var offsetCount = 0
			for (tracker in trackers) {
				val trackerFrame = tracker.tryGetFrame(cursor) ?: continue
				val position = trackerFrame.tryGetPosition()?.apply { this * scale } ?: continue
				val trackerRole = trackerFrame.tryGetTrackerPosition()?.trackerRole ?: continue

				val computedTracker = skeleton.getComputedTracker(trackerRole) ?: continue

				offset += (position - computedTracker.position).len()
				offsetCount++
			}
			return if (offsetCount > 0) offset / offsetCount else 0f
		}
	}
}
