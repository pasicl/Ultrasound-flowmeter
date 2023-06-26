# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 17:42:52 2022

@author: scharfetter_admin
"""

#!/usr/bin/env python

# Control program for the Red Pitaya Pulsed NMR system
# Copyright (C) 2015  Pavel Demin
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import struct
import os
import numpy as np
import time
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')
import pandas as pd
from scipy import signal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from PyQt5.uic import loadUiType
from PyQt5.QtCore import QRegExp, QTimer
from PyQt5.QtGui import QRegExpValidator
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QSizePolicy, QMessageBox, QWidget
from PyQt5.QtNetwork import QAbstractSocket, QTcpSocket

Ui_PulsedNMR, QMainWindow = loadUiType('pulsed_ultrasound.ui')

class PulsedNMR(QMainWindow, Ui_PulsedNMR):
  rates = {0:25.0e3, 1:50.0e3, 2:125.0e3, 3:250.0e3, 4:500.0e3, 5:1250.0e3}
  def __init__(self):
    super(PulsedNMR, self).__init__()
    self.setupUi(self)
    self.rateValue.addItems(['25', '50', '125', '250', '500', '1250'])
    # IP address validator
    rx = QRegExp('^(([0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])\.){3}([0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])|rp-[0-9A-Fa-f]{6}\.local$')
    self.addrValue.setValidator(QRegExpValidator(rx, self.addrValue))
    # state variable
    self.idle = True
    # number of samples to show on the plot
    self.size = 50000 #Total size of the signal
    self.counter_switch = 0 # Counter that resets every time the switch has completed the cycle
    self.counter_reps = 0 #Counter of total number of pulses sent 
    self.counter_avg = 0 #Counter of number of averages per flux rate
    # buffer and offset for the incoming samples
    self.buffer = bytearray(16 * self.size)
    self.offset = 0
    self.data = np.frombuffer(self.buffer, np.int32)
    # create figure
    figure = Figure()
    figure.set_facecolor('none')
    self.axes = figure.add_subplot(111)
    self.canvas = FigureCanvas(figure)
    self.plotLayout.addWidget(self.canvas)
    # create navigation toolbar
    self.toolbar = NavigationToolbar(self.canvas, self.plotWidget, False)
    # remove subplots action
    actions = self.toolbar.actions()
    if int(matplotlib.__version__[0]) < 2:
      self.toolbar.removeAction(actions[7])
    else:
      self.toolbar.removeAction(actions[6])
    self.plotLayout.addWidget(self.toolbar)
    # create TCP socket
    self.socket = QTcpSocket(self)
    self.socket.connected.connect(self.connected)
    self.socket.readyRead.connect(self.read_data)
    self.socket.error.connect(self.display_error)
    # connect signals from buttons and boxes
    self.startButton.clicked.connect(self.start)
    self.freqValue.valueChanged.connect(self.set_freq)
    self.deltaValue.valueChanged.connect(self.set_delta)
    self.rateValue.currentIndexChanged.connect(self.set_rate)
    # set rate
    self.rateValue.setCurrentIndex(3)
    # create timer for the repetitions
    self.startTimer = QTimer(self)
    self.startTimer.timeout.connect(self.timeout)
    self.timer = QTimer(self)
    self.timer.timeout.connect(self.start_sequence)

  def start(self):
    """
    Sets the sequence to start
    :param self: An instance of the class containing attributes
    :type self: Object
    """
    if self.idle:
      #Set initial values of the GPIOs
      self.clear_pin(2)
      self.set_pin(3)

      self.startButton.setEnabled(False)
      self.socket.connectToHost(self.addrValue.text(), 1001)
      self.startTimer.start(5000)

    else:
      self.stop()

  def stop(self):
    """
    Sets the sequence to stop
    :param self: An instance of the class containing attributes
    :type self: Object
    """    
    self.idle = True
    self.timer.stop()
    self.socket.abort()
    self.offset = 0
    self.startButton.setText('Start')
    self.startButton.setEnabled(True)
    self.clear_pin(2) #With this command one sets the GPIO value to "low" (0V)
    self.set_pin(3) #With this command one sets the GPIO value to "high" (3V3)

    self.counter_switch = 0

  def timeout(self):
    """
    Timeout: the time to connect to the IP adress has ran out
    :param self: An instance of the class containing attributes
    :type self: Object
    """      
    self.display_error('timeout')

  def connected(self):
    """
    Connected succesfully to the server.
    Sets the frequency, number of pulses, waiting time, sampling rate, starts the sequence and sets the GPIO
    :param self: An instance of the class containing attributes
    :type self: Object
    """      
    self.startTimer.stop()
    self.idle = False
    self.set_freq(self.freqValue.value())
    self.num_pulses = self.numPulses.value()
    self.time = self.wait_time.value()
    self.set_rate(self.rateValue.currentIndex())
    self.start_sequence()
    self.timer.start(self.deltaValue.value())
    self.startButton.setText('Stop')
    self.startButton.setEnabled(True)
    self.phi_M = np.zeros((self.numPulses.value(),1))
    self.clear_pin(2)
    self.set_pin(3)

    self.counter_switch = 0
    self.counter_reps = 0
    self.counter_avg = 0
  def read_data(self):
      
    """
    Read the received data
    Separate into different variables the I and Q data from the two ADCs
    Control the sequence with the counters.
    Call the self.process_pulse function
    When the sequence is finished, call self.results to finish
    :param self: An instance of the class containing attributes
    :type self: object

    """  
    size = self.socket.bytesAvailable()
    if self.offset + size < 16 * self.size:
      self.buffer[self.offset:self.offset + size] = self.socket.read(size)
      self.offset += size
    else:
      self.buffer[self.offset:16 * self.size] = self.socket.read(16 * self.size - self.offset)
      self.offset = 0
      # plot the signal 
     
      self.curve.set_ydata(np.real(self.data.astype(np.float32).view(np.complex64)[0::2] / (1 << 30)))
      self.canvas.draw()
      # self.curve.hold(True)
      
      # self.curve.set_ydata(np.real(self.data.astype(np.float32).view(np.complex64)[1::2] / (1 << 30)))
      # self.canvas.draw()

      data1 = self.data.astype(np.float32).view(np.complex64)[0::2] / (1 << 30)
      data2 = self.data.astype(np.float32).view(np.complex64)[1::2] / (1 << 30)
      self.accumulated_data = np.array(data1)
      data1_corrected_real = np.real(data1)
      data1_corrected_imag = np.imag(data1)
      data2_corrected_real = np.real(data2)
      data2_corrected_imag = np.imag(data2)

      self.my_array = np.array([data1_corrected_real, data1_corrected_imag, data2_corrected_real, data2_corrected_imag])
      
      phi = self.process_pulse(self.my_array)
      
      print("phi_",self.counter_reps, " = ", phi,"\n")
      self.phi_M[self.counter_reps] = phi
      self.dfarray = pd.DataFrame(self.my_array)
      
      
      switch = self.repPulses.value()
      #print(self.counter_switch, "\n")
      self.averages = self.avgValue.value()
      print("Pulse number = ",self.counter_reps,"\n Avg counter = ",self.counter_avg,"\n Switch counter = ",self.counter_switch,"\n")
      
      
         
      if self.counter_switch == switch-1:
          self.set_pin(2)
          self.clear_pin(3)

      if self.counter_switch == switch*2-1:
          self.clear_pin(2)
          self.set_pin(3)

          self.counter_switch = -1
          self.counter_avg = self.counter_avg + 1
      if self.counter_avg == self.averages:
          self.counter_avg = 0
          time.sleep(self.time)

      if self.counter_reps == self.num_pulses-1:
          self.results(self.phi_M)
          self.counter_reps = 0
          
      self.counter_switch = self.counter_switch + 1
      self.counter_reps = self.counter_reps + 1

      #Uncomment the following lines to save the raw signals as a csv file and add path
      #name = str(self.counter_reps)
      #cwd = os.getcwd()
      #path = cwd + "/new"
      #self.dfarray.to_csv(path+"pulse"+name+".csv", index=False)
      #self.df.to_csv('addpath/test'+name+'.csv', index=False)


  def display_error(self, socketError):
    """
    Displays an error
    :param self: An instance of the class containing attributes
    :type self: Object
    :param socketError: Whether there is an error in the socket communication
    :type socketError: string
    """      
    self.startTimer.stop()
    if socketError == 'timeout':
      QMessageBox.information(self, 'PulsedNMR', 'Error: connection timeout.')
    else:
      QMessageBox.information(self, 'PulsedNMR', 'Error: %s.' % self.socket.errorString())
    self.stop()

  def set_freq(self, value):
    """
    Sends the frequency of the pulse to the server application
    :param self: An instance of the class containing attributes
    :type self: Object
    :param value: Frequency value, in MHz
    :type value: double 
    """            
    if self.idle: return
    self.IF = self.IFValue.value()
    self.socket.write(struct.pack('<Q', 0<<60 | int(1.0e6 * value - self.IF)))
    self.socket.write(struct.pack('<Q', 1<<60 | int(1.0e6 * value)))

  def set_rate(self, index):
    """
    Sends the sampling rate to the server application
    Sets up the figure to plot the received pulse
    :param self: An instance of the class containing attributes
    :type self: Object
    :param index: Index of the sampling frequency
    :type index: int 
    """ 
    # time axis
    rate = float(PulsedNMR.rates[index])
    time = np.linspace(0.0, (self.size - 1) * 1000.0 / rate, self.size)
    # reset toolbar
    self.toolbar.home()
    self.toolbar.update()
    # reset plot
    self.axes.clear()
    self.axes.grid()
    # plot zeros and get store the returned Line2D object
    self.curve, = self.axes.plot(time, np.zeros(self.size))
    x1, x2, y1, y2 = self.axes.axis()
    # set y axis limits
    self.axes.axis((x1, x2, -1.1, 1.1))
    self.axes.set_xlabel('time, ms')
    self.canvas.draw()
    if self.idle: return
    self.socket.write(struct.pack('<Q', 2<<60 | int(125.0e6 / rate / 2)))

  def set_delta(self, value):
    """
    Sets the repeatition time of each pulse
    :param self: An instance of the class containing attributes
    :type self: Object
    :param delta: Repetition time in ms
    :type delta: double 
    """       
    if self.idle: return
    self.timer.stop()
    self.timer.start(value)

  def clear_pulses(self):
    """
    Stops sending the pulse
    :param self: An instance of the class containing attributes
    :type self: Object
    """        
    if self.idle: return
    self.socket.write(struct.pack('<Q', 7<<60))

  def add_delay(self, width):
    """
    Adds a delay in the sequence (not used)
    :param self: An instance of the class containing attributes
    :type self: Object
    :param width: Length of the delay
    :type width: double 
    """        
    if self.idle: return
    self.socket.write(struct.pack('<Q', 8<<60 | int(width - 4)))

  def add_pulse(self, level, phase, width):
    """
    Sends the amplitude, phase and width of the pulse to the server application
    :param self: An instance of the class containing attributes
    :type self: Object
    :param level: Amplitude, encoded from 0 to 32766
    :type level: int 
    :param phase: Phase of the sent pulse (is always set to 0 in this application)
    :type: int
    :param width: Width of the pulse, as 125*length, in us
    :type: double
    """   
    if self.idle: return
    self.socket.write(struct.pack('<Q', 8<<60 | int(width)))
    self.socket.write(struct.pack('<Q', 9<<60 | int(phase << 16 | level)))
    
  def set_pin(self, pin):
    """
    Sets the GPIO from extension connector E1 in "high" state
    :param self: An instance of the class containing attributes
    :type self: Object
    :param pin: number of GPIO: 0 -> PIN 3, 1-> PIN 5, every odd PIN of E1 of the GPIO, PIN 3:2:17; the even ones are not working
    :type width: int 
    """         
    if self.idle: return
    self.socket.write(struct.pack('<Q', 4<<60 | int(pin)))

  def clear_pin(self, pin):
    """
    Sets the GPIO from extension connector E2 in "low" state
    :param self: An instance of the class containing attributes
    :type self: Object
    :param pin: number of GPIO: 0 -> PIN 3, 1-> PIN 5, every odd PIN of E1 of the GPIO, PIN 3:2:17; the even ones are not working
    :type width: int 
    """          
    if self.idle: return
    self.socket.write(struct.pack('<Q', 5<<60 | int(pin)))
      

  def start_sequence(self):
    """
    Starts the sequence
    :param self: An instance of the class containing attributes
    :type self: Object
    """          
    if self.idle: return

    awidth = 125 * self.awidthValue.value()
    size = self.size
    self.clear_pulses()
    self.add_pulse(32766, 0, awidth) 

    self.socket.write(struct.pack('<Q', 10<<60 | int(size)))

    
  def process_pulse(self,vector):
      """
      Filter the signals with Spicy library 
      Apply phase estimation algorithm with the filtered signals
      Return phase
      :param self: An instance of the class containing attributes
      :type self: object
      :param vector: variable with all the received data
      :type vector: array
      :return: the value of the phase of the pulse
          - phi_i: float.32
      """
      if self.idle: return
      self.IF = self.IFValue.value()
      #print(float(PulsedNMR.rates[self.rateValue.currentIndex()]))
      l = signal.firwin(49,1.12*self.IF,pass_zero = 'lowpass',fs = float(PulsedNMR.rates[self.rateValue.currentIndex()]))
      h = signal.firwin(49,0.16*self.IF,pass_zero = 'highpass',fs = float(PulsedNMR.rates[self.rateValue.currentIndex()]))
      filtered_1 = np.array([np.convolve(xi, l, mode='valid') for xi in vector])
      filtered_signals = np.array([np.convolve(xi, h, mode='valid') for xi in filtered_1])
      re = np.vdot(filtered_signals[0][3200:],filtered_signals[2][3200:])
      im = np.vdot(filtered_signals[0][3200:],filtered_signals[3][3200:])
      phi_i = np.arctan(im/re);
      if im > 0 and re < 0:
        phi_i = np.pi - abs(phi_i)
      elif im < 0 and re < 0:
        phi_i = phi_i + np.pi;
      elif im < 0 and re > 0:
        phi_i = 2*np.pi - abs(phi_i)
   
      return phi_i
    
  def results(self,phi):
      """
      Separate the pulses in both directions
      Plot the phase all the phases as a function of the pulse number, the phases on both directions and their difference 
      Save the phases in .csv file
      :param self: An instance of the class containing attributes
      :type self: Object
      :param phi: vector with the phases of all the pulses
      :type phi: np.array
      """
      MM = self.numPulses.value()
      sw = self.repPulses.value()
      phase_1 = np.zeros(int(MM/(2*sw)),float)
      phase_2 = np.zeros(int(MM/(2*sw)),float)
      P = len(phase_1)
      for ii in range(0,P):
          phase_1[ii] = np.mean(phi[int(((ii+1)*(2*sw))-(2*sw-1)):int(((ii+1)*(2*sw))-sw)])
          phase_2[ii] = np.mean(phi[int(((ii+1)*(2*sw))-(sw-1)):int(((ii+1)*(2*sw)))])

      dif=phase_2-phase_1;
      
      print(phi,"\n")
      plt.figure(figsize=(20, 12))
      plt.plot(range(0,MM),phi , label="diff")
      plt.figure(figsize=(20, 12))
      plt.plot(range(0,P),dif , label="diff")

      plt.grid()
      plt.title("Phase difference", fontdict={'size': 24})
      plt.figure(figsize=(20, 12))
      plt.plot(range(0,P),phase_1 , label="phase1")     
      plt.plot(range(0,P),phase_2 , label="phase2")     
      plt.title("The two phases", fontdict={'size': 24})

 
      self.dfphi = pd.DataFrame(phi)
      cwd = os.getcwd()
      path = cwd + "/stored_data/"
      self.dfphi.to_csv(path+"data.csv", index=False)

      self.stop()



app = QApplication(sys.argv)
window = PulsedNMR()
window.show()
sys.exit(app.exec_())
