package dev.slimevr.posestreamer.trackers;

import dev.slimevr.vr.trackers.TrackerPosition;
import dev.slimevr.vr.trackers.TrackerStatus;


public class PSRTracker {
	public String name;
	public TrackerStatus trackerStatus;
	public TrackerPosition designation;

	public PSRTracker() {
		// Default constructor, doesn't need anything since it can all be set
		// afterwards
	}

	public boolean serializeDef() {
		return true;
	}
}
