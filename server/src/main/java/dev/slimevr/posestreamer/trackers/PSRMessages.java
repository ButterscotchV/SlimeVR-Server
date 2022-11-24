package dev.slimevr.posestreamer.trackers;

public enum PSRMessages {
	ADD_TRACKERS(0),
	TRACKERS_DATA(1),
	BODY_OFFSETS(2),
	SET_FRAMERATE(3),
	;

	public final int value;

	private PSRMessages(int value) {
		this.value = value;
	}
}
