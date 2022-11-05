package dev.slimevr.posestreamer.trackers;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;

import dev.slimevr.posestreamer.PoseDataStream;
import dev.slimevr.posestreamer.PoseStreamer;
import dev.slimevr.vr.processor.skeleton.Skeleton;


public class PSRFileStream extends PoseDataStream {

	protected final BufferedOutputStream bufferedOutputStream;

	public PSRFileStream(OutputStream outputStream) {
		super(outputStream);
		bufferedOutputStream = new BufferedOutputStream(this.outputStream, 4096);
	}

	public PSRFileStream(File file) throws FileNotFoundException {
		super(file);
		bufferedOutputStream = new BufferedOutputStream(outputStream, 4096);
	}

	public PSRFileStream(String file) throws FileNotFoundException {
		super(file);
		bufferedOutputStream = new BufferedOutputStream(outputStream, 4096);
	}

	@Override
	public void writeHeader(Skeleton skeleton, PoseStreamer streamer) throws IOException {
		if (skeleton == null) {
			throw new NullPointerException("skeleton must not be null");
		}
		if (streamer == null) {
			throw new NullPointerException("streamer must not be null");
		}

		// Header >w<
	}

	@Override
	public void writeFrame(Skeleton skeleton) throws IOException {
		// Wow cool frame
	}

	@Override
	public void writeFooter(Skeleton skeleton) throws IOException {
		// No footer :P
	}

	@Override
	public void close() throws IOException {
		bufferedOutputStream.close();
		super.close();
	}
}
