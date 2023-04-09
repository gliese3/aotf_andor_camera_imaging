import serial

class Obis():
    def __init__(self):
        self.BAUDRATE = 9600
        self.PORT = "COM4"
        self.ser_device = serial.Serial(port=self.PORT,
                                        baudrate=self.BAUDRATE,
                                        timeout=0.1) #TODO: should be changed
    
    
    def laserOff(self):
        self.ser_device.write(bytes("SOUR:AM:STATE OFF\r\n", "utf-8")) # turn OFF the laser
        self.ser_device.readline()
        
    
    def laserOn(self):
        self.ser_device.write(bytes("SOUR:AM:STATE ON\r\n", "utf-8")) #turn ON laser
        self.ser_device.readline()
        
    def setCwPowerMode(self):
        self.ser_device.write(bytes("SOUR:AM:INT CWP\r\n", "utf-8")) # set operation mode to internal CW Power
        self.ser_device.readline()
        

    def setPower(self, power):
        power = power / 1000 # in milli watts
        self.ser_device.write(bytes(f"SOUR:POW:LEV:IMM:AMPL {power}\r\n", "utf-8")) # set power level
        self.ser_device.readline()
        
    
    def getCurrentPower(self):
        self.ser_device.write(bytes("SOUR:POW:LEV:IMM:AMPL?\r\n", "utf-8"))
        print(self.ser_device.readline())
        