import inaccel.coral as inaccel
import numpy as np
import time

class StereoBM:
	def __init__(self, cameraMA_l=None, cameraMA_r=None, distC_l=None, distC_r=None, irA_l=None, irA_r=None, bm_state=None ):
		# allocate mem for camera parameters for rectification and bm_state class
		with inaccel.allocator:
			if cameraMA_l is None:
				self.cameraMA_l_fl = np.array([933.173, 0.0, 663.451, 0.0, 933.173, 377.015, 0.0, 0.0, 1.0], dtype=np.float32)
			else:
				self.cameraMA_l_fl = np.array(cameraMA_l, dtype=np.float32)

			if cameraMA_r is None:
				self.cameraMA_r_fl = np.array([933.467, 0.0, 678.297, 0.0, 933.467, 359.623, 0.0, 0.0, 1.0], dtype=np.float32)
			else:
				self.cameraMA_r_fl = np.array(cameraMA_r, dtype=np.float32)

			if distC_l is None:
				self.distC_l_fl = np.array([-0.169398, 0.0227329, 0.0, 0.0, 0.0], dtype=np.float32)
			else:
				self.distC_l_fl = np.array(distC_l, dtype=np.float32)

			if distC_r is None:
				self.distC_r_fl = np.array([-0.170581, 0.0249444, 0.0, 0.0, 0.0], dtype=np.float32)
			else:
				self.distC_r_fl = np.array(distC_r, dtype=np.float32)

			if irA_l is None:
				self.irA_l_fl = np.array([0.0011976323, -0.0000000019, -0.8153011732, 0.0000000007, 0.0011976994, \
	                      								-0.4422348617,  0.0000126839,  0.0000001064, 0.9913820905], dtype=np.float32)
			else:
				self.irA_l_fl = np.array(irA_l, dtype=np.float32)

			if irA_r is None:
				self.irA_r_fl = np.array([0.0011976994,  0.0000000000, -0.8047567905, -0.0000000000, 0.0011976994, \
	                      -0.4420566166, -0.0000000000, -0.0000001064,  1.0000392898], dtype=np.float32)
			else:
				self.irA_r_fl = np.array(irA_r, dtype=np.float32)

			if bm_state is None:
				self.bm_state_arr = np.array([0, 15, 31, 15, 0, 48, 20, 15, 16, 3, 0], dtype=np.int32)
			else:
				self.bm_state_arr = np.array(bm_state, dtype=np.int32)


	def runAsync(self, left_img, right_img):
		self.m_runStartTime = int(round(time.time() * 1000000))

		if left_img is None:
			raise RuntimeError('Invalid left image')
		if right_img is None:
			raise RuntimeError('Invalid right image')
		if left_img.shape[0] != right_img.shape[0] or left_img.shape[1] != right_img.shape[1]:
			raise RuntimeError('Image sizes differ')

		# allocate and initialize buffers
		rows = np.int32(left_img.shape[0]);
		cols = np.int32(left_img.shape[1]);

		with inaccel.allocator:
			self.left_mat = np.array(left_img)
			self.right_mat = np.array(right_img)
			self.disp_mat = np.ndarray((rows, cols), dtype=np.uint16)

		# Create request for stereo accelerator
		req = inaccel.request('com.xilinx.vitis.vision.stereoBM')
		req.arg(self.left_mat).arg(self.right_mat).arg(self.disp_mat)
		req.arg(self.cameraMA_l_fl).arg(self.cameraMA_r_fl)
		req.arg(self.distC_l_fl).arg(self.distC_r_fl)
		req.arg(self.irA_l_fl).arg(self.irA_r_fl)
		req.arg(self.bm_state_arr)
		req.arg(rows).arg(cols)

		self.response = inaccel.submit(req)

	def wait(self):
		# Send request and wait for completion
		self.response.result()

		# Write output image
		disp_mat_scaled = (self.disp_mat.view(np.ndarray)*(256.0 / 48.0) / (16.0)).astype(np.uint8)

		self.m_runEndTime = int(round(time.time() * 1000000))
		return disp_mat_scaled;

	def run(self, left_img, right_img):
		self.runAsync(left_img, right_img)
		return self.wait()

	def lastruntime(self):
		duration = self.m_runEndTime - self.m_runStartTime
		return duration
