package dev.slimevr.tracking.processor.config

import dev.slimevr.config.SkeletonConfig
import kotlin.reflect.KParameter

annotation class OffsetEnumVal(val enumVal: SkeletonConfigOffsets)
annotation class ToggleEnumVal(val enumVal: SkeletonConfigToggles)
annotation class ValueEnumVal(val enumVal: SkeletonConfigValues)

data class SkeletonConfigImmutable(
	// Offsets
	@OffsetEnumVal(SkeletonConfigOffsets.HEAD)
	val headShift: Float = 0.1f,
	@OffsetEnumVal(SkeletonConfigOffsets.NECK)
	val neckLength: Float = 0.1f,
	@OffsetEnumVal(SkeletonConfigOffsets.CHEST)
	val chestLength: Float = 0.32f,
	@OffsetEnumVal(SkeletonConfigOffsets.CHEST_OFFSET)
	val chestOffset: Float = 0.0f,
	@OffsetEnumVal(SkeletonConfigOffsets.WAIST)
	val waistLength: Float = 0.20f,
	@OffsetEnumVal(SkeletonConfigOffsets.HIP)
	val hipLength: Float = 0.04f,
	@OffsetEnumVal(SkeletonConfigOffsets.HIP_OFFSET)
	val hipOffset: Float = 0.0f,
	@OffsetEnumVal(SkeletonConfigOffsets.HIPS_WIDTH)
	val hipsWidth: Float = 0.26f,
	@OffsetEnumVal(SkeletonConfigOffsets.UPPER_LEG)
	val upperLegLength: Float = 0.42f,
	@OffsetEnumVal(SkeletonConfigOffsets.LOWER_LEG)
	val lowerLegLength: Float = 0.50f,
	@OffsetEnumVal(SkeletonConfigOffsets.FOOT_LENGTH)
	val footLength: Float = 0.05f,
	@OffsetEnumVal(SkeletonConfigOffsets.FOOT_SHIFT)
	val footShift: Float = -0.05f,
	@OffsetEnumVal(SkeletonConfigOffsets.SKELETON_OFFSET)
	val skeletonOffset: Float = 0.0f,
	@OffsetEnumVal(SkeletonConfigOffsets.SHOULDERS_DISTANCE)
	val shouldersDistance: Float = 0.08f,
	@OffsetEnumVal(SkeletonConfigOffsets.SHOULDERS_WIDTH)
	val shouldersWidth: Float = 0.36f,
	@OffsetEnumVal(SkeletonConfigOffsets.UPPER_ARM)
	val upperArmLength: Float = 0.26f,
	@OffsetEnumVal(SkeletonConfigOffsets.LOWER_ARM)
	val lowerArmLength: Float = 0.26f,
	@OffsetEnumVal(SkeletonConfigOffsets.HAND_Y)
	val handDistanceY: Float = 0.035f,
	@OffsetEnumVal(SkeletonConfigOffsets.HAND_Z)
	val handDistanceZ: Float = 0.13f,
	@OffsetEnumVal(SkeletonConfigOffsets.ELBOW_OFFSET)
	val elbowOffset: Float = 0.0f,
	// Toggles
	@ToggleEnumVal(SkeletonConfigToggles.EXTENDED_SPINE_MODEL)
	val extendedSpine: Boolean = true,
	@ToggleEnumVal(SkeletonConfigToggles.EXTENDED_PELVIS_MODEL)
	val extendedPelvis: Boolean = true,
	@ToggleEnumVal(SkeletonConfigToggles.EXTENDED_KNEE_MODEL)
	val extendedKnee: Boolean = true,
	@ToggleEnumVal(SkeletonConfigToggles.FORCE_ARMS_FROM_HMD)
	val forceArmsFromHMD: Boolean = true,
	@ToggleEnumVal(SkeletonConfigToggles.FLOOR_CLIP)
	val floorClip: Boolean = true,
	@ToggleEnumVal(SkeletonConfigToggles.SKATING_CORRECTION)
	val skatingCorrection: Boolean = true,
	@ToggleEnumVal(SkeletonConfigToggles.VIVE_EMULATION)
	val viveEmulation: Boolean = false,
	@ToggleEnumVal(SkeletonConfigToggles.TOE_SNAP)
	val toeSnap: Boolean = false,
	@ToggleEnumVal(SkeletonConfigToggles.FOOT_PLANT)
	val footPlant: Boolean = true,
	// Values
	@ValueEnumVal(SkeletonConfigValues.WAIST_FROM_CHEST_HIP_AVERAGING)
	val waistFromChestHipAveraging: Float = 0.45f,
	@ValueEnumVal(SkeletonConfigValues.WAIST_FROM_CHEST_LEGS_AVERAGING)
	val waistFromChestLegsAveraging: Float = 0.2f,
	@ValueEnumVal(SkeletonConfigValues.HIP_FROM_CHEST_LEGS_AVERAGING)
	val hipFromChestLegsAveraging: Float = 0.45f,
	@ValueEnumVal(SkeletonConfigValues.HIP_FROM_WAIST_LEGS_AVERAGING)
	val hipFromWaistLegsAveraging: Float = 0.4f,
	@ValueEnumVal(SkeletonConfigValues.HIP_LEGS_AVERAGING)
	val hipLegsAveraging: Float = 0.25f,
	@ValueEnumVal(SkeletonConfigValues.KNEE_TRACKER_ANKLE_AVERAGING)
	val kneeTrackerAnkleAveraging: Float = 0.75f,
) {
	fun getOffset(offset: SkeletonConfigOffsets): Float {
		return when (offset) {
			SkeletonConfigOffsets.HEAD -> headShift
			SkeletonConfigOffsets.NECK -> neckLength
			SkeletonConfigOffsets.CHEST -> chestLength
			SkeletonConfigOffsets.CHEST_OFFSET -> chestOffset
			SkeletonConfigOffsets.WAIST -> waistLength
			SkeletonConfigOffsets.HIP -> hipLength
			SkeletonConfigOffsets.HIP_OFFSET -> hipOffset
			SkeletonConfigOffsets.HIPS_WIDTH -> hipsWidth
			SkeletonConfigOffsets.UPPER_LEG -> upperLegLength
			SkeletonConfigOffsets.LOWER_LEG -> lowerLegLength
			SkeletonConfigOffsets.FOOT_LENGTH -> footLength
			SkeletonConfigOffsets.FOOT_SHIFT -> footShift
			SkeletonConfigOffsets.SKELETON_OFFSET -> skeletonOffset
			SkeletonConfigOffsets.SHOULDERS_DISTANCE -> shouldersDistance
			SkeletonConfigOffsets.SHOULDERS_WIDTH -> shouldersWidth
			SkeletonConfigOffsets.UPPER_ARM -> upperArmLength
			SkeletonConfigOffsets.LOWER_ARM -> lowerArmLength
			SkeletonConfigOffsets.HAND_Y -> handDistanceY
			SkeletonConfigOffsets.HAND_Z -> handDistanceZ
			SkeletonConfigOffsets.ELBOW_OFFSET -> elbowOffset
		}
	}

	fun getToggle(toggle: SkeletonConfigToggles): Boolean {
		return when (toggle) {
			SkeletonConfigToggles.EXTENDED_SPINE_MODEL -> extendedSpine
			SkeletonConfigToggles.EXTENDED_PELVIS_MODEL -> extendedPelvis
			SkeletonConfigToggles.EXTENDED_KNEE_MODEL -> extendedKnee
			SkeletonConfigToggles.FORCE_ARMS_FROM_HMD -> forceArmsFromHMD
			SkeletonConfigToggles.FLOOR_CLIP -> floorClip
			SkeletonConfigToggles.SKATING_CORRECTION -> skatingCorrection
			SkeletonConfigToggles.VIVE_EMULATION -> viveEmulation
			SkeletonConfigToggles.TOE_SNAP -> toeSnap
			SkeletonConfigToggles.FOOT_PLANT -> footPlant
		}
	}

	fun getValue(value: SkeletonConfigValues): Float {
		return when (value) {
			SkeletonConfigValues.WAIST_FROM_CHEST_HIP_AVERAGING -> waistFromChestHipAveraging
			SkeletonConfigValues.WAIST_FROM_CHEST_LEGS_AVERAGING -> waistFromChestLegsAveraging
			SkeletonConfigValues.HIP_FROM_CHEST_LEGS_AVERAGING -> hipFromChestLegsAveraging
			SkeletonConfigValues.HIP_FROM_WAIST_LEGS_AVERAGING -> hipFromWaistLegsAveraging
			SkeletonConfigValues.HIP_LEGS_AVERAGING -> hipLegsAveraging
			SkeletonConfigValues.KNEE_TRACKER_ANKLE_AVERAGING -> kneeTrackerAnkleAveraging
		}
	}

	companion object {
		@JvmStatic
		val default = SkeletonConfigImmutable()

		@JvmStatic
		fun createFromConfigReflected(config: SkeletonConfig): SkeletonConfigImmutable {
			val constructor = SkeletonConfigImmutable::class.constructors.first()
			val paramMap = HashMap<KParameter, Any>()
			for (parameter in constructor.parameters) {
				for (annotation in parameter.annotations) {
					when (annotation) {
						is OffsetEnumVal -> {
							paramMap[parameter] =
								config.offsets[annotation.enumVal.configKey] ?: continue
						}
						is ToggleEnumVal -> {
							paramMap[parameter] =
								config.toggles[annotation.enumVal.configKey] ?: continue
						}
						is ValueEnumVal -> {
							paramMap[parameter] =
								config.values[annotation.enumVal.configKey] ?: continue
						}
					}
				}
			}
			return constructor.callBy(paramMap)
		}

		@JvmStatic
		fun createFromConfigGenerated(config: SkeletonConfig): SkeletonConfigImmutable {
			return SkeletonConfigImmutable(
				// Offsets
				config.offsets[SkeletonConfigOffsets.HEAD.configKey] ?: default.headShift,
				config.offsets[SkeletonConfigOffsets.NECK.configKey] ?: default.neckLength,
				config.offsets[SkeletonConfigOffsets.CHEST.configKey] ?: default.chestLength,
				config.offsets[SkeletonConfigOffsets.CHEST_OFFSET.configKey] ?: default.chestOffset,
				config.offsets[SkeletonConfigOffsets.WAIST.configKey] ?: default.waistLength,
				config.offsets[SkeletonConfigOffsets.HIP.configKey] ?: default.hipLength,
				config.offsets[SkeletonConfigOffsets.HIP_OFFSET.configKey] ?: default.hipOffset,
				config.offsets[SkeletonConfigOffsets.HIPS_WIDTH.configKey] ?: default.hipsWidth,
				config.offsets[SkeletonConfigOffsets.UPPER_LEG.configKey] ?: default.upperLegLength,
				config.offsets[SkeletonConfigOffsets.LOWER_LEG.configKey] ?: default.lowerLegLength,
				config.offsets[SkeletonConfigOffsets.FOOT_LENGTH.configKey] ?: default.footLength,
				config.offsets[SkeletonConfigOffsets.FOOT_SHIFT.configKey] ?: default.footShift,
				config.offsets[SkeletonConfigOffsets.SKELETON_OFFSET.configKey] ?: default.skeletonOffset,
				config.offsets[SkeletonConfigOffsets.SHOULDERS_DISTANCE.configKey] ?: default.shouldersDistance,
				config.offsets[SkeletonConfigOffsets.SHOULDERS_WIDTH.configKey] ?: default.shouldersWidth,
				config.offsets[SkeletonConfigOffsets.UPPER_ARM.configKey] ?: default.upperArmLength,
				config.offsets[SkeletonConfigOffsets.LOWER_ARM.configKey] ?: default.lowerArmLength,
				config.offsets[SkeletonConfigOffsets.HAND_Y.configKey] ?: default.handDistanceY,
				config.offsets[SkeletonConfigOffsets.HAND_Z.configKey] ?: default.handDistanceZ,
				config.offsets[SkeletonConfigOffsets.ELBOW_OFFSET.configKey] ?: default.elbowOffset,
				// Toggles
				config.toggles[SkeletonConfigToggles.EXTENDED_SPINE_MODEL.configKey] ?: default.extendedSpine,
				config.toggles[SkeletonConfigToggles.EXTENDED_PELVIS_MODEL.configKey] ?: default.extendedPelvis,
				config.toggles[SkeletonConfigToggles.EXTENDED_KNEE_MODEL.configKey] ?: default.extendedKnee,
				config.toggles[SkeletonConfigToggles.FORCE_ARMS_FROM_HMD.configKey] ?: default.forceArmsFromHMD,
				config.toggles[SkeletonConfigToggles.FLOOR_CLIP.configKey] ?: default.floorClip,
				config.toggles[SkeletonConfigToggles.SKATING_CORRECTION.configKey] ?: default.skatingCorrection,
				config.toggles[SkeletonConfigToggles.VIVE_EMULATION.configKey] ?: default.viveEmulation,
				config.toggles[SkeletonConfigToggles.TOE_SNAP.configKey] ?: default.toeSnap,
				config.toggles[SkeletonConfigToggles.FOOT_PLANT.configKey] ?: default.footPlant,
				// Values
				config.values[SkeletonConfigValues.WAIST_FROM_CHEST_HIP_AVERAGING.configKey] ?: default.waistFromChestHipAveraging,
				config.values[SkeletonConfigValues.WAIST_FROM_CHEST_LEGS_AVERAGING.configKey] ?: default.waistFromChestLegsAveraging,
				config.values[SkeletonConfigValues.HIP_FROM_CHEST_LEGS_AVERAGING.configKey] ?: default.hipFromChestLegsAveraging,
				config.values[SkeletonConfigValues.HIP_FROM_WAIST_LEGS_AVERAGING.configKey] ?: default.hipFromWaistLegsAveraging,
				config.values[SkeletonConfigValues.HIP_LEGS_AVERAGING.configKey] ?: default.hipLegsAveraging,
				config.values[SkeletonConfigValues.KNEE_TRACKER_ANKLE_AVERAGING.configKey] ?: default.kneeTrackerAnkleAveraging
			)
		}
	}
}
