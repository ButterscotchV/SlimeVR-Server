package dev.slimevr.autobone

import dev.slimevr.config.AutoBoneConfig
import dev.slimevr.config.ConfigManager
import dev.slimevr.poseframeformat.PoseFrames
import dev.slimevr.poseframeformat.player.TrackerFramesPlayer
import dev.slimevr.tracking.processor.BoneType
import dev.slimevr.tracking.processor.HumanPoseManager
import java.util.*
import java.util.function.Consumer

class AutoBoneStep(
	val config: AutoBoneConfig,
	val targetHmdHeight: Float,
	val targetFullHeight: Float,
	val trainingFrames: PoseFrames,
	val intermediateOffsets: EnumMap<BoneType, Float>,
	val epochCallback: Consumer<AutoBone.Epoch?>?,
	serverConfig: ConfigManager,
	var curEpoch: Int = 0,
	var curAdjustRate: Float = 0f,
	var cursor1: Int = 0,
	var cursor2: Int = 0,
	var currentHmdHeight: Float = 0f,
) {

	// The scale to normalize the skeleton with
	val normalSkeletonScale: Float = 1f / targetHmdHeight

	var skeletonScale: Float = 1f
		set(value) {
			field = value
			framePlayer1.setScales(value)
			framePlayer2.setScales(value)
		}

	val maxFrameCount = trainingFrames.maxFrameCount

	val framePlayer1 = TrackerFramesPlayer(trainingFrames)
	val framePlayer2 = TrackerFramesPlayer(trainingFrames)

	val trackers1 = framePlayer1.trackers.toList()
	val trackers2 = framePlayer2.trackers.toList()

	val skeleton1 = HumanPoseManager(trackers1)
	val skeleton2 = HumanPoseManager(trackers2)

	val errorStats = StatsCalculator()

	init {
		// Scale the players to be normalized
		framePlayer1.setScales(skeletonScale)
		framePlayer2.setScales(skeletonScale)
		// Load server configs into the skeleton
		skeleton1.loadFromConfig(serverConfig)
		skeleton2.loadFromConfig(serverConfig)
		// Disable leg tweaks, this will mess with the resulting positions
		skeleton1.setLegTweaksEnabled(false)
		skeleton2.setLegTweaksEnabled(false)
	}

	fun setCursors(cursor1: Int, cursor2: Int, updatePlayerCursors: Boolean) {
		this.cursor1 = cursor1
		this.cursor2 = cursor2

		if (updatePlayerCursors) {
			updatePlayerCursors()
		}
	}

	fun updatePlayerCursors() {
		framePlayer1.setCursors(cursor1)
		framePlayer2.setCursors(cursor2)
		skeleton1.update()
		skeleton2.update()
	}

	val heightOffset: Float
		get() = targetHmdHeight - currentHmdHeight
}
