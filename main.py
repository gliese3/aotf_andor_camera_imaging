import sys
import os
from collections.abc import Iterable
from PyQt5 import QtCore, QtWidgets, QtGui

import pyqtgraph as pg

import numpy as np
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib_additional import MplCanvas
import matplotlib as mpl
from multiprocessing import Pool
import itertools

from obis_laser import Obis 

import andor_camera
from brimrose_aotf import Aotf
from multi_process import buildSpectra

from interface import Ui_MainWindow

class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, *args, **kwargs):        

        # from compiled file
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.setWindowTitle("Andor real time spectrum")
        
        # set icon
        self.setWindowIcon(QtGui.QIcon("images/microscope.png"))
        
#================================ CONNECTIONS =================================

        self.apply_settings_but.clicked.connect(self.onApplySettingsButton)
        self.start_exp_but.clicked.connect(self.onStartExperimentButton)
        
        self.time_step1_slider.valueChanged.connect(self.onSlider1ValueChanged)
        self.time_step1_spBox.valueChanged.connect(self.onTimeStep1SpBoxValueChanged)
        
        self.time_step2_slider.valueChanged.connect(self.onSlider2ValueChanged)
        self.time_step2_spBox.valueChanged.connect(self.onTimeStep2SpBoxValueChanged)
        
        self.time_step3_slider.valueChanged.connect(self.onSlider3ValueChanged)
        self.time_step3_spBox.valueChanged.connect(self.onTimeStep3SpBoxValueChanged)
        
        self.time_step2_slider.valueChanged.connect(self.plotSpectra)
        self.time_step2_spBox.valueChanged.connect(self.plotSpectra)
        
        self.time_step3_slider.valueChanged.connect(self.plotIntensityRatioSpectra)
        self.time_step3_spBox.valueChanged.connect(self.plotIntensityRatioSpectra)
        
        self.plot_spectra_but.clicked.connect(self.onPlotSpectra)
        self.plot_pixels_but.clicked.connect(self.onPixelPlot)
        
        self.plot_intensity_but.clicked.connect(self.onPlotIntensityRatioSpectra)
        
        self.save_raw_data_action.triggered.connect(self.onSaveRawData)
        self.save_spectra_action.triggered.connect(self.onSaveSpectraData)
        self.save_config_and_times_action.triggered.connect(self.onSaveConfigAndTimes)
        
        self.cooling_chBox.stateChanged.connect(self.onCoolerChackbox)

#================================ /CONNECTIONS ================================

        # TIMERS        
        self.camera_curr_temp_timer = QtCore.QTimer()
        self.camera_curr_temp_timer.timeout.connect(self.updateCameraTemp)

        # VARIABLES
        self.plot1_elements = [
            self.wavelen_chBox,
            self.time_step1_slider,
            self.time_step1_spBox,
        ]
        self.plot2_elements = [
            self.time_step2_slider,
            self.time_step2_spBox,
            self.plot_pixels_but
        ]
         
         
        # CONSTANTS
        #!: should be modified for common case
        #!: now it is OK
        self.IMG_HEIGHT = 512
        self.IMG_WIDTH = 512
        
        self.x_lim = None
        self.y_lim = None
        
        # INIT CAMERA 
        try:
            self.camera = andor_camera.AndorCamera()
            self.camera.coolCamera(-80)
        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.critical(self, "Camera", "Camera connection problem.")
            return
        else:
            
            # read camera temp every 500 millis
            self.camera_curr_temp_timer.start(500)
        
        # INIT AOTF
        try:
            self.aotf = Aotf()
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "AOTF", "AOTF init problem.")
            print(e)
            return
        
        # INIT OBIS LASER
        try:
            self.obis = Obis()
            self.obis.setCwPowerMode()
            self.obis.setPower(0)
            
            # turn on laser, alhough while power is 0, it is off 
            self.obis.laserOn()
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Obis laser", "Obis laser init problem.")
            print(e)
            return
        
        # open window in maximized size
        self.showMaximized()
        
        # use pyqtgraph for speed
        # basic adjustments
        pg.setConfigOptions(imageAxisOrder='row-major')
        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')
        self.graphWidget1 = pg.ImageView()
        
        self.canvas2 = MplCanvas(self) # canvas to plot
        self.toolbar2 = NavigationToolbar(self.canvas2, self) # navigation toolbar
        
        self.canvas3 = MplCanvas(self) # canvas to plot
        self.toolbar3 = NavigationToolbar(self.canvas3, self) # navigation toolbar
        
        self.canvas4 = MplCanvas(self) # canvas to plot
        self.toolbar4 = NavigationToolbar(self.canvas4, self) # navigation toolbar

        
        # arrange matplotlib elements in gui
        self.plot1_layout.addWidget(self.graphWidget1)
        
        self.plot2_layout.addWidget(self.toolbar2)
        self.plot2_layout.addWidget(self.canvas2)
        self.canvas2.fig.subplots(1)
        
        self.plot3_layout.addWidget(self.toolbar3)
        self.plot3_layout.addWidget(self.canvas3)
        self.canvas3.fig.subplots(1)
        
        self.plot4_layout.addWidget(self.toolbar4)
        self.plot4_layout.addWidget(self.canvas4)
        self.canvas4.fig.subplots(1)
        
        # just init plot (not necessary)
        self.canvas2.draw()
        self.canvas3.draw()
        
#=================================== SLOTS ====================================                
    
    def onCoolerChackbox(self):
        
        # check cooling
        if self.cooling_chBox.isChecked():
            self.camera.coolCamera(-80)
        else:
            self.camera.coolerOff()
    
    
    def plotIntensityRatioSpectra(self):
        
        wavelen_array_str = [line.strip() for line in self.wavelengths_for_intensity_lineEdit.text().split(",")]
        wavelen_array_int = [int(wavelen) for wavelen in wavelen_array_str]
        
        iter_norm_time_step = self.iteration_start_spBox.value()        
        current_time_step = self.curr_time_step3
        
        for i, wavelen in enumerate(wavelen_array_int):
            wavelen_ind = self.WAVELEN_IND_MAP[wavelen]
            init_intensity_image = self.images_arr[wavelen_ind, iter_norm_time_step]
            image = self.images_arr[wavelen_ind, current_time_step]
            image = np.divide(image, init_intensity_image)
            
            # check normalization
            if self.normalization_intensity_chBox.isChecked():
                vmin = self.images_norm_dict[wavelen]["vmin"]
                vmax = self.images_norm_dict[wavelen]["vmax"]
            else:
                vmin = np.amin(image)
                vmax = np.amax(image)
                
            # plot
            cmap = self.colormap3_comBox.currentText()
            interpolation = self.interpolation3_comBox.currentText()
            self.im4_arr[i].set(
                data=image,
                clim=(vmin, vmax),
                cmap=cmap,
                interpolation=interpolation)
            
        # draw everything
        self.canvas4.draw_idle()
        
    
    def onPlotIntensityRatioSpectra(self):
        
        # check wavelengths
        wavelen_array_str = [line.strip() for line in self.wavelengths_for_intensity_lineEdit.text().split(",")]
        if len(set(wavelen_array_str).intersection(self.WAVELEN_ARRAY_STR)) != len(wavelen_array_str):
            QtWidgets.QMessageBox.critical(self, "Wavelengths", "Improper wavelengths. Check them again.")
            return
        
        wavelen_array_int = [int(wavelen) for wavelen in wavelen_array_str]
        
        # adjust prog bar
        self.inten_progBar.setRange(0, len(wavelen_array_int) - 1)
        
        # create images norm dict (needed if normalization option is checked)
        self.images_norm_dict = {wavelen : {} for wavelen in wavelen_array_int}
        
        # find colorbar normalization bounds (needed if normalization option is checked)
        iter_norm_time_step = self.iteration_start_spBox.value()
        for i, wavelen in enumerate(wavelen_array_int):
            wavelen_ind = self.WAVELEN_IND_MAP[wavelen]
            init_intensity_image = self.images_arr[wavelen_ind, iter_norm_time_step]            

            vmin = float("inf")
            vmax = float("-inf")
            for iter in range(self.total_num_of_steps):
                loc_min = np.amin( np.divide(self.images_arr[wavelen_ind, iter], init_intensity_image) )
                loc_max = np.amax( np.divide(self.images_arr[wavelen_ind, iter], init_intensity_image) )
                if loc_min < vmin : vmin = loc_min
                if loc_max > vmax : vmax = loc_max

            # update prog bar value
            self.inten_progBar.setValue(i)
            
            # process all gui events
            QtWidgets.QApplication.processEvents()
         
            self.images_norm_dict[wavelen]["vmin"] = vmin
            self.images_norm_dict[wavelen]["vmax"] = vmax
        
        # pre-plot
        wavelen_array_str = [line.strip() for line in self.wavelengths_for_intensity_lineEdit.text().split(",")]
        wavelen_array_int = [int(wavelen) for wavelen in wavelen_array_str]
        
        # prepare canvas for plot
        wavlen_num = len(wavelen_array_int)
        self.canvas4.fig.clear()
        axs = self.canvas4.fig.subplots(1, wavlen_num)    
        
        iter_norm_time_step = self.iteration_start_spBox.value()
        cmap = self.colormap3_comBox.currentText()  
        current_time_step = 0
        interpolation = self.interpolation3_comBox.currentText()
        
        self.im4_arr = []
        
        for wavelen, axis in zip(wavelen_array_int, axs):
            wavelen_ind = self.WAVELEN_IND_MAP[wavelen]
            init_intensity_image = self.images_arr[wavelen_ind, iter_norm_time_step]
            image = self.images_arr[wavelen_ind, current_time_step]
            image = np.divide(image, init_intensity_image)
            
            # check normalization
            if self.normalization_intensity_chBox.isChecked():
                vmin = self.images_norm_dict[wavelen]["vmin"]
                vmax = self.images_norm_dict[wavelen]["vmax"]
            else:
                vmin = np.amin(image)
                vmax = np.amax(image)
                
            # plot
            im = axis.imshow(image,
                        vmin=vmin,
                        vmax=vmax,
                        cmap=cmap,
                        extent=(0, 511, 511, 0),
                        interpolation=interpolation)
            self.im4_arr.append(im)
            divider = make_axes_locatable(axis)
            cax = divider.append_axes("right", size="5%", pad=0.08)
            
            self.canvas4.fig.colorbar(im, cax=cax)
            axis.set_title(f"{wavelen} nm")
        
        # init plot
        self.time_step3_slider.setValue(0)
        
        
    def plotRawData(self):
        weve_len_str = self.wavelen_chBox.currentText()
        time_step = self.curr_time_step1
        
        # to avoid error after multiple applies
        if not weve_len_str:
            return
        
        weve_len_ind = self.WAVELEN_IND_MAP[int(weve_len_str)]
        self.time1_lab.setText(str(self.elapsed_times[time_step][weve_len_ind]))
        
        # PLOT DATA                
        image = self.images_arr[weve_len_ind, time_step]
        
        # find which nornalization should be used
        norm = self.checkNormalization()
        if norm == "global":
            
            # to use full color range
            max_val = self.max_image_val
            vmin = self.min_image_val / max_val
            vmax = 1
            image = image / max_val
            
        elif norm == "wavelength":
            if self.previous_norm != "wavelength":
                max_val = np.amax(self.images_arr[:, time_step])
                vmin = np.amin(self.images_arr[:, time_step]) / max_val
                vmax = 1
                image = image / max_val
            else:
                vmin = self.previous_vmin
                max_val = self.previous_maxval
                vmax = 1
                image = image / max_val
            
            self.previous_vmin = vmin
            self.previous_maxval = max_val
            self.previous_norm = "wavelength"
            
        elif norm == "time":
            if self.previous_norm != "time":
                max_val = np.amax(self.images_arr[weve_len_ind, :])
                vmin = np.amin(self.images_arr[weve_len_ind, :]) / max_val
                vmax = 1
                image = image / max_val
            else:
                vmin = self.previous_vmin
                max_val = self.previous_maxval
                vmax = 1
                image = image / max_val
            
            self.previous_vmin = vmin
            self.previous_maxval = max_val
            self.previous_norm = "time"
            
        elif norm == "none":
            vmin = np.amin(image)
            vmax = np.amax(image)                
        
        # plot image
        self.graphWidget1.setImage(image,
                                   levels=(vmin, vmax),
                                   autoRange=False,
                                   autoHistogramRange=False)


    def plotSpectra(self):
        time_step = self.curr_time_step2
        base_cmap = self.colormap2_comBox.currentText()
        interpolation = self.interpolation2_comBox.currentText()
        
        # our image
        image = self.max_inten_spectra_dict[time_step]
        
        if base_cmap != self.base_cmap_old:
            self.base_cmap_old = base_cmap
            noise_val = 0
            upper_lim = self.WAVELEN_ARRAY_INT[-1] + 1 # only for plotting purposes
            WAVE_LEN_ARRAY_EXT = [noise_val] + self.WAVELEN_ARRAY_INT
            WAVE_LEN_ARRAY_EXT_WITH_LIM = WAVE_LEN_ARRAY_EXT + [upper_lim]
            # add "pre-value" to WAVE_LEN_ARRAY to be considered as noise value
            # it is done because of simplicity implementing it to matplotlib
            # WAVE_LEN_ARRAY_EXT[0] means noise value

            # create proper cmap
            init_cmap = mpl.cm.get_cmap(self.base_cmap_old, len(WAVE_LEN_ARRAY_EXT))

            # extract all colors from cmap
            cmaplist = [init_cmap(i) for i in range(len(WAVE_LEN_ARRAY_EXT))]

            # force the first color entry to be grey (for noise)
            cmaplist[0] = (.5, .5, .5, 1.0)

            # create the new map
            cmap = mpl.colors.LinearSegmentedColormap.from_list("Custom cmap",
                                                                cmaplist,
                                                                len(WAVE_LEN_ARRAY_EXT))
            norm = mpl.colors.BoundaryNorm(WAVE_LEN_ARRAY_EXT_WITH_LIM, cmap.N)
            
            # plot image
            self.im2.set(
                data=image,
                interpolation=interpolation,
                cmap=cmap,
                norm=norm)
        else:
            # plot image
            self.im2.set(
                data=image,
                interpolation=interpolation)        
        
        # draw everything
        self.canvas2.draw_idle()
        
        
    def onApplySettingsButton(self):
        
        # check run settings
        # a list of int at the end
        self.total_num_of_steps_in_run_arr = np.array(list(map(int, map(str.strip, self.steps_in_each_run_lineEdit.text().split(',')))))
        
        # a list of float at the end
        self.time_pause_in_run_arr = np.array(list(map(float, map(str.strip, self.time_pause_in_run_lineEdit.text().split(',')))))
        
        # a list of float at the end
        self.time_pause_between_run_arr = np.array(list(map(float, map(str.strip, self.time_pause_between_runs_lineEdit.text().split(',')))))
        
        # integration times for camera
        # a list of float at the end
        self.integration_times_arr = np.array(list(map(float, map(str.strip, self.integ_times_lineEdit.text().split(',')))))
        
        # a list of float at the end
        self.laser_powers_arr = np.array(list(map(float, map(str.strip, self.laser_powers_lineEdit.text().split(',')))))
        
        check_list = map(len, 
                         [self.total_num_of_steps_in_run_arr, 
                         self.time_pause_in_run_arr,
                         self.time_pause_between_run_arr,
                         self.integration_times_arr,
                         self.laser_powers_arr])
        
        # select unique elements
        check_list = set(check_list)
        
        # all settings arrays should be the same size
        if len(check_list) != 1: #
            QtWidgets.QMessageBox.critical(self, "Experiment settings", "Check experiment settings.")
            return
            
        # disconnect to avoid errors
        try:
            self.time_step1_slider.valueChanged.disconnect(self.plotRawData)
            self.time_step1_spBox.valueChanged.disconnect(self.plotRawData)
            self.wavelen_chBox.currentTextChanged.disconnect(self.plotRawData)
        except Exception as e:
            print(e)
        
        # disable gui elements
        for el in self.plot1_elements + self.plot2_elements:
            el.setEnabled(False)
            
        # some init to avoid errors
        self.curr_time_step1 = 0
        
        # clear combo box
        self.wavelen_chBox.clear()
        
        # CAMERA SETTINGS
        self.camera.setPreAmpGain(self.pre_amp_gain_comBox.currentIndex())
        self.camera.setEMGainMode(1)
        
        if self.emccd_gain_chBox.isChecked():
            self.camera.setEMCCDGain(self.emccd_gain_spBox.value())
        else:
            self.camera.setEMCCDGain(0)
            
        # self.camera.setExposureTime(self.exposure_spBox.value())
        self.camera.setAcquisitionMode(mode=1) # single scan
        self.camera.setReadMode(mode=4) # image
        self.camera.setTriggerMode(mode=0) # internal trigger
        self.camera.setImage(hbin=self.hor_binning_spBox.value(),
                             vbin=self.ver_binning_spBox.value(),
                             hstart=1,
                             hend=512,
                             vstart=1,
                             vend=512)
            
        # INIT EXPERIMENT
        
        # collect of wavelengths we need
        wavelen_start = self.wavelen_start_spBox.value()
        wavelen_finish = self.wavelen_finish_spBox.value()
        wavelen_step = self.wavelen_step_spBox.value()
        init_wavelen_array_str = list(map(str, range(wavelen_start, wavelen_finish + wavelen_step, wavelen_step)))
        additional_wavelen_array_str = [line.strip() for line in self.wavelen_array_lineEdit.text().split(",")]
        
        if additional_wavelen_array_str[0]: # not an empty list 
            self.WAVELEN_ARRAY_STR = sorted(init_wavelen_array_str + additional_wavelen_array_str)
        else:
            self.WAVELEN_ARRAY_STR = sorted(init_wavelen_array_str)
        
        self.WAVELEN_ARRAY_INT = [int(wavelen) for wavelen in self.WAVELEN_ARRAY_STR]
        
        # create index-wavelength map
        self.WAVELEN_IND_MAP = dict( zip(self.WAVELEN_ARRAY_INT,
                                         range(len(self.WAVELEN_ARRAY_INT))) )
        
        # AOTF SETTINGS
        self.aotf.setAmplitude(self.aotf_ampl_spBox.value()) 
        
        # total number of runs
        self.num_of_runs = len(self.total_num_of_steps_in_run_arr)
            
        # total number of steps (for all runs)
        self.total_num_of_steps = np.sum(self.total_num_of_steps_in_run_arr)
        
        # GUI SETTINGS
        self.exp_progBar.setRange(0, self.total_num_of_steps)
        self.exp_progBar.setValue(0)
        
        self.wavelen_chBox.addItems(self.WAVELEN_ARRAY_STR)
        
        self.time_step1_slider.setRange(0, self.total_num_of_steps - 1)
        self.time_step1_spBox.setRange(0, self.total_num_of_steps - 1)
        
        self.time_step3_slider.setRange(0, self.total_num_of_steps - 1)
        self.time_step3_spBox.setRange(0, self.total_num_of_steps - 1)        
        
        # approx time duration, 0.05 sec - overhead
        time_dur = (np.sum(self.total_num_of_steps_in_run_arr * (self.integration_times_arr + 0.05) * len(self.WAVELEN_ARRAY_INT)) +
                    np.sum(self.total_num_of_steps_in_run_arr * self.time_pause_in_run_arr)  + 
                    np.sum(self.time_pause_between_run_arr)) # in sec
        self.exp_duration_lab.setText(f"Aproximate experiment duration is {round(time_dur / 60, 1)} min")
        
        # add connections
        self.time_step1_slider.valueChanged.connect(self.plotRawData)
        self.time_step1_spBox.valueChanged.connect(self.plotRawData)
        self.wavelen_chBox.currentTextChanged.connect(self.plotRawData)
        
        
    def updateCameraTemp(self):
        temp, _ = self.camera.getCameraTemp()
        self.cam_temp_value_lab.setText(str(temp))


    def onSaveRawData(self):
        save_folder = QtWidgets.QFileDialog.getExistingDirectory(self, 
                                                    "Select folder to save data")


        if not save_folder: # no folder chosen
            return
        
        # create "raw_data" folder
        raw_data_path = os.path.normpath(os.path.join(save_folder, "raw_data"))
        os.mkdir(raw_data_path)
        
        for j, wavelen in enumerate(self.WAVELEN_ARRAY_INT):
            
            # make wavelength folders
            save_path = os.path.normpath(os.path.join(save_folder, "raw_data", str(wavelen)))
            os.mkdir(save_path)
            
            for iter in range(self.total_num_of_steps):
                image = self.images_arr[j, iter]
                np.savetxt(f"{save_path}/{iter}.csv", image, delimiter="\t", fmt="%d")
    
    
    def onSaveConfigAndTimes(self):
        save_folder = QtWidgets.QFileDialog.getExistingDirectory(self, 
                                                    "Select folder to save config")

        if not save_folder: # no folder chosen
            return
        
        # save config file
        with open(f"{save_folder}/config.txt", 'w') as f:
            f.write(f"number of runs: {len(self.total_num_of_steps_in_run_arr)}\n")
            f.write(f"number of steps in each run: {self.total_num_of_steps_in_run_arr}\n")
            f.write(f"delay between steps: {self.time_pause_in_run_arr} sec\n")
            f.write(f"pause between runs: {self.time_pause_between_run_arr} sec\n")
            f.write(f"ocean optics integration times: {self.integration_times_arr} millis\n")
            f.write(f"laser power for runs: {self.laser_powers_arr} mW (%)")
        
        
        # save times file
        # save times in a separate file
        np.savetxt(f"{save_folder}/times.csv", self.elapsed_times,
           fmt="%4.3f", 
           delimiter=",",
           header=f"Steps are in a column, different wavelengths are in a row",
           comments="")
        
    def onSaveSpectraData(self):
        save_folder = QtWidgets.QFileDialog.getExistingDirectory(self, 
                                                    "Select folder to save data")

        if not save_folder: # no folder
            return
        
        
        # make wavelength folders
        save_path = os.path.normpath(os.path.join(save_folder, "spectra_data"))
        os.mkdir(save_path)
        
        for iter in range(self.total_num_of_steps):
            image =  self.max_inten_spectra_dict[iter]
            np.savetxt(f"{save_path}/{iter}.csv", image, delimiter="\t", fmt="%d")
        
        
    def onStartExperimentButton(self):
        
        # some init
        self.emergency_stop = False
        
        but_name = self.start_exp_but.text()
        
        if but_name == "Start experiment":
            self.start_exp_but.setText("Stop experiment")
            
            # disable gui elements
            for el in self.plot1_elements + self.plot2_elements:
                el.setEnabled(False)
            self.plot_spectra_but.setEnabled(False)
            self.apply_settings_but.setEnabled(False)
            
            # open camera shutter
            self.camera.setShutter(1, 1, 0, 0)
            
            # for some reason making the first image in series
            # for the camera takes much longer time than subsequent.
            # So, make first image before our real exmeriment series 
            self.camera.getImage()
            
            # init images array
            self.images_arr = np.zeros((len(self.WAVELEN_ARRAY_INT),
                                        self.total_num_of_steps,
                                        self.IMG_HEIGHT,
                                        self.IMG_WIDTH), 
                                        dtype=np.uint16)
            
            self.elapsed_times = np.zeros((self.total_num_of_steps, len(self.WAVELEN_ARRAY_INT)))
            
            # start time for the experiment
            start_time = time.time()
            
            i = 0 # step number
            
            for run in range(self.num_of_runs):
                
                # set laser power
                self.camera.setExposureTime(self.integration_times_arr[run])
                self.obis.setPower(self.laser_powers_arr[run])
                self.obis.laserOn()
            
                for step in range(self.total_num_of_steps_in_run_arr[run]):
                            
                    for j, wave_len in enumerate(self.WAVELEN_ARRAY_INT):
                        
                        if self.emergency_stop:
                            
                            #TODO: add actions here if needed
                            
                            self.camera.abortAcquisition()
                            
                            # close camera shutter
                            self.camera.setShutter(1, 2, 0, 0)
                            return
                        
                        self.statusbar.showMessage(f"RUN: {run + 1} OF {self.num_of_runs} | STEP: {step + 1} OF {self.total_num_of_steps_in_run_arr[run]}", 3000)
                        self.aotf.setWavelength(wave_len)
                        
                        # st = time.time()
                        image = self.camera.getImage()
                        # print(f"image: {time.time() - st} sec")
                        
                        # append elapsed time
                        self.elapsed_times[i, j] = round(time.time() - start_time, 2)
                        
                        # save image
                        self.images_arr[j, i] = image
                        
                        self.exp_progBar.setValue(i + 1)
                        
                        # process all gui events
                        QtWidgets.QApplication.processEvents()
                        
                    i += 1
                    
                    self.statusbar.showMessage("PAUSE...")
                    
                    # process all gui events
                    QtWidgets.QApplication.processEvents()
                    
                    # make pause here
                    time.sleep(self.time_pause_in_run_arr[run])
                    
                if self.time_pause_between_run_arr[run] > 1E-6:
                    
                    # turn off the laser
                    self.obis.laserOff()
                    
                    if run != self.num_of_runs - 1: # last run, no need to wait
            
                        self.statusbar.showMessage("PAUSE BETWEEN RUNS...")
                        QtWidgets.QApplication.processEvents()
                        time.sleep(self.time_pause_between_run_arr[run])
                
            # turn off the laser
            self.obis.laserOff()    
            
            # find global max and min values
            self.max_image_val = np.amax(self.images_arr)
            self.min_image_val = np.amin(self.images_arr)
                
            self.statusbar.showMessage(f"Experiment has finished!", 5000)    
                
            # close camera shutter
            self.camera.setShutter(1, 2, 0, 0)
            
            # enable some gui elements
            for el in self.plot1_elements:
                el.setEnabled(True)
            self.plot_spectra_but.setEnabled(True)
            self.apply_settings_but.setEnabled(True)
            
            # change button name if experiment has finished successfully
            self.start_exp_but.setText("Start experiment")
            
            # some init for plot normalization
            self.previous_norm = None
            self.previous_vmin = None
            self.previous_maxval = None
            
            # init slider and canvas
            self.time_step1_slider.setValue(0)
            self.plotRawData()
        
        elif but_name == "Stop experiment":
            self.start_exp_but.setText("Start experiment")
            self.emergency_stop = True
            self.apply_settings_but.setEnabled(True)
            
        
    def onSlider1ValueChanged(self, value):
        self.time_step1_spBox.setValue(value)
        self.curr_time_step1 = value
        
        
    def onTimeStep1SpBoxValueChanged(self, value):
        self.time_step1_slider.setValue(value)
        self.curr_time_step1 = value
        
        
    def onSlider2ValueChanged(self, value):
        self.time_step2_spBox.setValue(value)
        self.curr_time_step2 = value
        self.time2_lab.setText(str(self.elapsed_times[value][0]))
        
        
    def onTimeStep2SpBoxValueChanged(self, value):
        self.time_step2_slider.setValue(value)
        self.curr_time_step2 = value
        self.time2_lab.setText(str(self.elapsed_times[value][0]))
        
        
    def onSlider3ValueChanged(self, value):
        self.time_step3_spBox.setValue(value)
        self.curr_time_step3 = value
        self.time3_lab.setText(str(self.elapsed_times[value][0]))
        
        
    def onTimeStep3SpBoxValueChanged(self, value):
        self.time_step3_slider.setValue(value)
        self.curr_time_step3 = value
        self.time3_lab.setText(str(self.elapsed_times[value][0]))
        

    def onPlotSpectra(self):
        
        # get threshold from gui
        self.THRESHOLD = self.threshold_spBox.value()
        
        # disable some gui elements
        self.save_spectra_action.setEnabled(False)
        for el in self.plot2_elements:
            el.setEnabled(False)
        
        # multiprocessing part 
        self.statusbar.showMessage("Making spectra...", 5000)
        # process all gui events
        QtWidgets.QApplication.processEvents()
        
        # prepare arguments for function
        steps           = range(self.total_num_of_steps)
        wavelen_arrs    = np.array([self.WAVELEN_ARRAY_INT] * self.total_num_of_steps, dtype=np.uint16)
        
        images_arrs = []
        for step in range(self.total_num_of_steps):
            images_arrs.append(self.images_arr[:, step, :, :])  
        
        thresholds      = [self.THRESHOLD] * self.total_num_of_steps
        img_heights     = [self.IMG_HEIGHT] * self.total_num_of_steps
        img_widths      = [self.IMG_WIDTH] * self.total_num_of_steps
        args = zip(steps, wavelen_arrs, images_arrs, thresholds, img_heights, img_widths)
        
        print(f"Array size is {self.images_arr.nbytes / 1048576 } Mb")
        
        # run multiprocessing
        # looks like processes=4 is the fastest (for small number of steps)
        st = time.time()
        with Pool() as pool:
            results = pool.starmap(buildSpectra, args)
        fn =   time.time()      
            
        self.statusbar.showMessage(f"Making spectra... Done! [{round(fn - st, 1)} secs]")
        # process all gui events
        QtWidgets.QApplication.processEvents()  

        results = list(itertools.chain(*results)) # to flatten list
        self.max_inten_spectra_dict = dict(results[0::2])
        self.spectra_dict = dict(results[1::2])
        
        # clear intentionally memory (may be not needed)
        del steps, wavelen_arrs, images_arrs, thresholds, img_heights, img_widths, args
            
        # gui
        self.time_step2_slider.setRange(0, self.total_num_of_steps - 1)
        self.time_step2_spBox.setRange(0, self.total_num_of_steps - 1)
        
        # pre-plotting
        time_step = 0
        self.base_cmap_old = self.colormap2_comBox.currentText()
        interpolation = self.interpolation2_comBox.currentText()
        
        # create proper cmap
        noise_val = 0
        upper_lim = self.WAVELEN_ARRAY_INT[-1] + 1 # only for plotting purposes
        WAVE_LEN_ARRAY_EXT = [noise_val] + self.WAVELEN_ARRAY_INT
        WAVE_LEN_ARRAY_EXT_WITH_LIM = WAVE_LEN_ARRAY_EXT + [upper_lim]
        # add "pre-value" to WAVE_LEN_ARRAY to be considered as noise value
        # it is done because of simplicity implementing it to matplotlib
        # WAVE_LEN_ARRAY_EXT[0] means noise value

        # create proper cmap
        init_cmap = mpl.cm.get_cmap(self.base_cmap_old, len(WAVE_LEN_ARRAY_EXT))

        # extract all colors from cmap
        cmaplist = [init_cmap(i) for i in range(len(WAVE_LEN_ARRAY_EXT))]

        # force the first color entry to be grey (for noise)
        cmaplist[0] = (.5, .5, .5, 1.0)

        # create the new map
        cmap = mpl.colors.LinearSegmentedColormap.from_list("Custom cmap",
                                                            cmaplist,
                                                            len(WAVE_LEN_ARRAY_EXT))
        norm = mpl.colors.BoundaryNorm(WAVE_LEN_ARRAY_EXT_WITH_LIM, cmap.N)

        # ticks for colorbar
        ticks = [ WAVE_LEN_ARRAY_EXT_WITH_LIM[i] + (WAVE_LEN_ARRAY_EXT_WITH_LIM[i + 1] - WAVE_LEN_ARRAY_EXT_WITH_LIM[i]) / 2 for i in range(len(WAVE_LEN_ARRAY_EXT_WITH_LIM[:-1])) ]

        # tick labels
        tick_labs = ["noise"] + [ f"{i} nm" for i in WAVE_LEN_ARRAY_EXT[1:] ]
        
        # PLOT DATA
        # clear axis
        self.canvas2.fig.clear()
        axs = self.canvas2.fig.subplots(1)
        
        # our image
        image = self.max_inten_spectra_dict[time_step]
        
        y_shape, x_shape = image.shape # reversed x and y            
        
        # plot image
        self.im2 = axs.imshow(image,
                         cmap=cmap,
                         norm=norm,
                         interpolation=interpolation,
                         extent=[0, x_shape, y_shape, 0])
        # previous line is needed for proper functioning onMouseClick() function
        
        # colorbar settings
        cb2 = self.canvas2.fig.colorbar(self.im2, ticks=ticks)
        cb2.ax.set_yticklabels(tick_labs)     
        
        # draw everything
        self.canvas2.draw()
        
        # init slider to plot
        self.time_step2_slider.setValue(0)
        
        # add mouse clicked connection
        self.cid = self.canvas2.fig.canvas.mpl_connect("button_press_event",
                                                       self.onMouseClick)
        
        # enable some gui elements
        self.save_spectra_action.setEnabled(True)
        for el in self.plot2_elements:
            el.setEnabled(True)
        
    
    def onPixelPlot(self):

        # clear fig every plot iteration
        self.canvas3.fig.clear()
        
        # axsinit
        axs_arr = None
                
        rows = self.pixels_tabWidget.rowCount()
        pixels_to_plot = []
        for row in range(rows):
            x_coord = int(self.pixels_tabWidget.item(row, 0).text())
            y_coord = int(self.pixels_tabWidget.item(row, 1).text())
            pixels_to_plot.append((y_coord, x_coord))
        
        if not len(pixels_to_plot): # nothing to plot
            self.canvas3.fig.clear()
            self.canvas3.draw()
            return
        
        # plot spectra from pixels
        axs_arr = self.canvas3.fig.subplots(len(pixels_to_plot))
        wavelengths = self.WAVELEN_ARRAY_INT
        
        time_steps = np.linspace(0,
                                 self.total_num_of_steps - 1,
                                 10,
                                 dtype=int)
        
        if isinstance(axs_arr, Iterable):
            for i, axs in enumerate(axs_arr):
                for step in time_steps: # key is time iteration
                    image = self.spectra_dict[step]
                    spectrum = image[pixels_to_plot[i]][:, 1] # only intensities
                    axs.plot(wavelengths, spectrum, label=f"Iteration: {step}")
                axs.set_title(f"Pixel {tuple(reversed(pixels_to_plot[i]))}")
                axs.grid(True)
                axs.legend() # no space to plot
                
        else: # only one pixel
            for step in time_steps: # key is time iteration
                image = self.spectra_dict[step]
                spectrum = image[pixels_to_plot[0]][:, 1] # only intensities
                axs_arr.plot(wavelengths, spectrum, label=f"Iteration: {step}")
                axs_arr.set_title(f"Pixel {tuple(reversed(pixels_to_plot[0]))}")
                axs_arr.grid(True)
                axs_arr.legend() # no space to plot
            
        self.canvas3.draw()
        
        
    def closeEvent(self, event):
        # close serial port for AOTF and Obis laser
        try:
            self.aotf.stopAotf()
            self.obis.laserOff()
        except Exception as e:
            print(e)
        
        # abort camera possible acquisition
        self.camera.abortAcquisition()
             
        event.accept()
                   
#=================================== /SLOTS ====================================

#================================= FUNCTIONS ===================================  
   
    def checkNormalization(self):
        # check radio buttons
        if self.global_norm_radBut.isChecked():
            return "global"
        elif self.wavelen_norm_radBut.isChecked():
            return "wavelength"
        elif self.time_norm_radBut.isChecked():
            return "time"
        elif self.none_norm_radBut.isChecked():      
            return "none"  
               
   
    def onMouseClick(self, event):
        if not event.inaxes or int(event.button) != 2: return # only middle button is valid
        x_coord, y_coord = int(np.floor(event.xdata)), int(np.floor(event.ydata))
        
        currentRowCount = self.pixels_tabWidget.rowCount()
        self.pixels_tabWidget.insertRow(currentRowCount)
        x_coord_item = QtWidgets.QTableWidgetItem(str(x_coord))
        y_coord_item = QtWidgets.QTableWidgetItem(str(y_coord))
        
        btn = QtWidgets.QPushButton("remove")
        btn.clicked.connect(lambda: self.pixels_tabWidget.removeRow(self.pixels_tabWidget.currentRow()))
        self.pixels_tabWidget.setCellWidget(currentRowCount,2, btn)
        
        self.pixels_tabWidget.setItem(currentRowCount,0, x_coord_item)
        self.pixels_tabWidget.setItem(currentRowCount,1, y_coord_item)               
       
#================================= /FUNCTIONS ==================================
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec_()