package dev.slimevr.posestreamer.trackers;

public enum PSRTrackerData {
	ROTATION(0),
	POSITION(1),
	ACCELERATION(2),
	COMPASS(3),
	COMPASS_ACCURACY(4),
	PING(5),
	;

	public final int flag;

	PSRTrackerData(int id) {
		this.flag = 1 << id;
	}

	public boolean check(int dataFlags) {
		return (dataFlags & this.flag) != 0;
	}
}
