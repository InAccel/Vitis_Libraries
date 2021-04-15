import cv2 as cv
import numpy as np
import math
import time

from ._vision_stereobm import StereoBM


class StereoBM_video:

	def __init__(self, parallel_requests=None):
		self.stereo = []
		self.parallel_requests = 16
		if parallel_requests is not None:
			self.parallel_requests = parallel_requests

		for i in range(0, self.parallel_requests):
			self.stereo.append( StereoBM() )

	def run(self, input_video, output_video, split_dimention=0):
		self.m_runStartTime = int(round(time.time() * 1000000))

		# try to open the video
		vid = cv.VideoCapture(input_video)
		if (vid.isOpened()== False):
			  raise RuntimeError("Error opening video stream or file")

		frame_num = int(vid.get(cv.CAP_PROP_FRAME_COUNT))
		if(split_dimention == 0):
			height = int(vid.get(cv.CAP_PROP_FRAME_HEIGHT))
			width  = int(vid.get(cv.CAP_PROP_FRAME_WIDTH)/2)
		else:
			height = int(vid.get(cv.CAP_PROP_FRAME_HEIGHT)/2)
			width  = int(vid.get(cv.CAP_PROP_FRAME_WIDTH))

		fourcc = int(vid.get(cv.CAP_PROP_FOURCC))
		fps = int(vid.get(cv.CAP_PROP_FPS))
		out = cv.VideoWriter(output_video, fourcc, fps, (height,width))

		frames = []
		for i in range(math.ceil(frame_num/self.parallel_requests)):
			frames.clear()

			for j in range(self.parallel_requests):
				ret, frame = vid.read()
				if ret == False:
					break
				else:
					frames.append(frame)

			for j in range(len(frames)):
				if(split_dimention == 0):
					frame_left = cv.cvtColor(cv.resize(frames[j][:,:width,:], (1280, 720)), cv.COLOR_BGR2GRAY)
					frame_right = cv.cvtColor(cv.resize(frames[j][:,width:,:], (1280, 720)), cv.COLOR_BGR2GRAY)
				else:
					frame_left = cv.cvtColor(cv.resize(frames[j][:height,:,:], (1280, 720)), cv.COLOR_BGR2GRAY)
					frame_right = cv.cvtColor(cv.resize(frames[j][height:,:,:], (1280, 720)), cv.COLOR_BGR2GRAY)

				self.stereo[j].runAsync(frame_left, frame_right)

			for j in range(len(frames)):
				disp_img = self.stereo[j].wait()
				disp_img = cv.resize(cv.cvtColor(disp_img, cv.COLOR_GRAY2BGR), (height,width))
				out.write(disp_img)

		vid.release()
		out.release()
		self.m_runEndTime = int(round(time.time() * 1000000))

	def lastruntime(self):
		duration = self.m_runEndTime - self.m_runStartTime
		return duration
