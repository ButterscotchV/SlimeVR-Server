package dev.slimevr.posestreamer.trackers;

public enum PSRTrackerDef {
	NAME(0),
	DESIGNATION(1),
	STATUS(2),
	;

	public final int flag;

	PSRTrackerDef(int id) {
		this.flag = 1 << id;
	}

	public boolean check(int dataFlags) {
		return (dataFlags & this.flag) != 0;
	}
}
