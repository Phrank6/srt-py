"""motors.py

Module for Controlling Different Motor Types over Serial

"""
import serial

from abc import ABC, abstractmethod
from time import sleep
from math import cos, acos, pi, sqrt, floor
# from pos import *
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, get_body, SkyCoord
from astropy import units as u
from astropy.time import Time
from datetime import datetime
import subprocess
import serial
import time
from scipy.optimize import least_squares



class Motor(ABC):
    """Abstract Class for All Motors Types

    Attributes
    ----------
    port : str
        Serial Port Identifier String for Communicating with the Motor
    baudrate : int
        Baudrate for serial connection
    az_limits : (float, float)
        Tuple of Lower and Upper Azimuth Limits
    el_limits : (float, float)
        Tuple of Lower and Upper Elevation Limits
    serial : serial.Serial
        Serial Object for Communicating with the Motor

    See Also
    --------
    <https://pyserial.readthedocs.io/en/latest/pyserial_api.html>
    """

    def __init__(self, port, baudrate, az_limits, el_limits):
        """Constructor for the Abstract Motor Class

        Parameters
        ----------
        port : str
            Serial Port Identifier String for Communicating with the Motor
        az_limits : (float, float)
            Tuple of Lower and Upper Azimuth Limits
        el_limits : (float, float)
            Tuple of Lower and Upper Elevation Limits
        """
        self.port = port
        self.baudrate = baudrate
        self.az_limits = az_limits
        self.el_limits = el_limits
        self.serial = None

    @abstractmethod
    def point(self, az, el):
        """Abstract Method Prototype for Pointing a Motor at an AzEl Coordinate

        Parameters
        ----------
        az : float
            Azimuth Coordinate Value to Point At
        el : float
            Elevation Coordinate Value to Point At

        Returns
        -------
        (float, float)
            Azimuth and Elevation Coordinate as a Tuple of Floats
        """
        pass

    @abstractmethod
    def status(self):
        """Abstract Method Prototype for Getting a Motor's Current AzEl Position

        Returns
        -------
        (float, float)
            Azimuth and Elevation Coordinate as a Tuple of Floats
        """
        pass

    def __del__(self):
        """Override of Motor Delete Method to Close Serial Port if Necessary

        Returns
        -------
        None
        """
        if self.serial is not None and self.serial.is_open:
            self.serial.close()


class NoMotor(Motor):
    """
    Class for Simulating a Motor or Using a Stationary Telescope
    """

    def __init__(self, port, baudrate, az_limits, el_limits):
        """
        Initializer for Rot2Motor

        Parameters
        ----------
        port : str
            NOT USED - Needed For Abstract Motor Initializer
        baudrate : int
            Baudrate for serial connection
        az_limits : (float, float)
            Tuple of Lower and Upper Azimuth Limits (if Stationary, both should be the same value)
        el_limits : (float, float)
            Tuple of Lower and Upper Elevation Limits (if Stationary, both should be the same value)
        """
        super().__init__(port, baudrate, az_limits, el_limits)
        self.position = (az_limits[0], el_limits[0])

    def point(self, az, el):
        """Changes the Unchanging Position of the Stationary / Simulated Motor

        Parameters
        ----------
        az : float
            Azimuth Coordinate to Point At
        el : float
            Elevation Coordinate to Point At

        Returns
        -------
        None
        """
        self.position = (az, el)

    def status(self):
        """Returns the Unchanging Position of the Stationary / Simulated Motor

        Returns
        -------
        (float, float)
            Current Azimuth and Elevation Coordinate as a Tuple of Floats
        """
        return self.position


class Rot2Motor(Motor):
    """Class for Controlling any ROT2 Protocol-Supporting Motor (e.g. SPID Motors)

    See Also
    --------
    <http://ryeng.name/blog/3>
    <https://github.com/jaidenfe/rot2proG/blob/master/rot2proG.py>
    <https://www.haystack.mit.edu/edu/undergrad/srt/pdf%20files/MD-01%20en.pdf>
    """

    VALID_PULSES_PER_DEGREE = (1, 2, 4)

    def __init__(
        self,
        port,
        baudrate,
        az_limits,
        el_limits,
        pulses_per_degree=2,
        test_pulses_per_degree=True,
    ):
        """Initializer for Rot2Motor

        Parameters
        ----------
        port : str
            Serial Port Identifier String for Communicating with the Motor
        baudrate : int
            Baudrate for serial connection
        az_limits : (float, float)
            Tuple of Lower and Upper Azimuth Limits
        el_limits : (float, float)
            Tuple of Lower and Upper Elevation Limits
        pulses_per_degree : int
            Number of Motor Pulses per Degree of Movement
        test_pulses_per_degree : bool
            Whether to Run A Call to Ask the Motor What its True Pulses per Degree Is (By Calling status)
        """
        Motor.__init__(self, port, baudrate, az_limits, el_limits)
        self.serial = serial.Serial(
            port=self.port,
            baudrate=baudrate,
            bytesize=serial.EIGHTBITS,
            parity="N",
            stopbits=serial.STOPBITS_ONE,
            timeout=None,
        )
        if pulses_per_degree in Rot2Motor.VALID_PULSES_PER_DEGREE:
            self.pulses_per_degree = pulses_per_degree
        else:
            raise ValueError("Invalid Pulse Per Degree Value")
        if test_pulses_per_degree:
            self.status()

    def send_rot2_pkt(self, cmd, az=None, el=None):
        """Builds and Sends a ROT2 Command Packet over Serial

        Parameters
        ----------
        cmd : int
            ROT2 Motor Command Value (0x2F -> Set, 0x1F -> Get, 0x0F -> Stop)
        az : float
            Azimuth Coordinate to Point At (If Applicable)
        el : float
            Elevation Coordinate to Point At (If Applicable)

        Notes
        -----
        All send_rot2_pkt calls should be followed with a receive_rot2_pkt

        Returns
        -------
        None
        """
        if az is not None and el is not None:
            azimuth = int(
                self.pulses_per_degree * (az + 360.0) + 0.5
            )  # Formatted Az Pulse Value
            elevation = int(
                self.pulses_per_degree * (el + 360.0) + 0.5
            )  # Formatted El Pulse Value
        else:
            azimuth = 0
            elevation = 0

        azimuth_ticks = (
            self.pulses_per_degree
        )  # Documentation for Rot2 Says This Is Ignored
        elevation_ticks = (
            self.pulses_per_degree
        )  # Documentation for Rot2 Says This Is Ignored

        cmd_string = "W%04d%c%04d%c%c " % (
            azimuth,
            azimuth_ticks,
            elevation,
            elevation_ticks,
            cmd,
        )
        cmd_bytes = cmd_string.encode("ascii")
        # print("Packet of Size " + str(len(cmd_bytes)))
        # print([hex(val) for val in cmd_bytes])
        self.serial.write(cmd_bytes)

    def receive_rot2_pkt(self):
        """Receives and Parsers an ROT2 Status Packet

        Returns
        -------
        (float, float)
            Azimuth and Elevation Coordinate as a Tuple of Floats
        """
        received_vals = self.serial.read(12)
        az = (
            (received_vals[1] * 100)
            + (received_vals[2] * 10)
            + received_vals[3]
            + (received_vals[4] / 10.0)
            - 360.0
        )
        el = (
            (received_vals[6] * 100)
            + (received_vals[7] * 10)
            + received_vals[8]
            + (received_vals[9] / 10.0)
            - 360.0
        )
        az_pulse_per_deg = received_vals[5]
        el_pulse_per_deg = received_vals[10]
        assert az_pulse_per_deg == el_pulse_per_deg  # Consistency Check
        if az_pulse_per_deg != self.pulses_per_degree:
            print(
                "Motor Pulses Per Degree Incorrect, Changing Value to "
                + str(az_pulse_per_deg)
            )
            self.pulses_per_degree = az_pulse_per_deg
        return az, el

    def point(self, az, el):
        """Point ROT2 Motor at AzEl Coordinate

        Parameters
        ----------
        az : float
            Azimuth Coordinate to Point At
        el : float
            Elevation Coordinate to Point At

        Returns
        -------
        None
        """
        cmd = 0x2F  # Rot2 Set Command
        az_relative = az - self.az_limits[0]
        el_relative = el - self.el_limits[0]
        self.send_rot2_pkt(cmd, az=az_relative, el=el_relative)

    def status(self):
        """Requests the Current Location of the ROT2 Motor

        Returns
        -------
        (float, float)
            Current Azimuth and Elevation Coordinate as a Tuple of Floats
        """
        cmd = 0x1F  # Rot2 Status Command
        self.send_rot2_pkt(cmd)
        az_relative, el_relative = self.receive_rot2_pkt()
        return az_relative + self.az_limits[0], el_relative + self.el_limits[0]

    def stop(self):
        """Stops the ROT2 Motor at its Current Location

        Returns
        -------
        None
        """
        cmd = 0x0F  # Rot2 Stop Command
        self.send_rot2_pkt(cmd)
        # az_relative, el_relative = self.receive_rot2_pkt()
        # return (az_relative + self.az_limits[0], el_relative + self.el_limits[0])


class H180Motor(Motor):  # TODO: Test!
    """
    Class for Controlling any ROT2 Protocol-Supporting Motor (e.g. SPID Motors)
    Based on the h180 function from the C SRT code

    ftp://gemini.haystack.mit.edu/pub/web/src/source_srt_newsrtsource_ver9.tar.gz
    """

    AZCOUNTS_PER_DEG = 52.0 * 27.0 / 120.0
    ELCOUNTS_PER_DEG = 52.0 * 27.0 / 120.0

    def __init__(self, port, baudrate, az_limits, el_limits, counts_per_step=100):
        """Initializer for the H180 Motor, baudrate should be 2400.

        Parameters
        ----------
        port : str
            Serial Port Identifier String for Communicating with the Motor
        baudrate : int
            Baudrate for serial connection
        az_limits : (float, float)
            Tuple of Lower and Upper Azimuth Limits
        el_limits : (float, float)
            Tuple of Lower and Upper Elevation Limits
        counts_per_step : int
            Maximum number of counts to move per call to function
        """
        Motor.__init__(self, port, az_limits, el_limits)
        self.serial = serial.Serial(
            port=port,
            baudrate=baudrate,  # 2400,
            bytesize=serial.EIGHTBITS,
            parity=serial.PARITY_NONE,
            stopbits=serial.STOPBITS_ONE,
            timeout=None,
        )
        self.count_per_step = counts_per_step
        self.az_lower_lim = az_limits[0]
        self.el_lower_lim = el_limits[0]
        self.az_count = 0.0
        self.el_count = 0.0

    def send_h180_cmd(self, az, el, stow):
        """Sends a Command to the H180 Motor

        Parameters
        ----------
        az : float
            Azimuth Coordinate to Point At
        el : float
            Elevation Coordinate to Point At
        stow : bool
            Whether or Not to Stow Antenna (makes az,el irrelevant)

        Returns
        -------
        None
        """
        azz = az - self.az_lower_lim
        ell = el - self.el_lower_lim
        for axis in range(2):
            mm = -1
            count = 0
            if stow:
                if axis == 0:
                    mm = 0
                else:
                    mm = 2
                count = 8000
            else:
                if axis == 0:
                    acount = azz * H180Motor.AZCOUNTS_PER_DEG - self.az_count
                    if self.count_per_step and acount > self.count_per_step:
                        acount = self.count_per_step
                    if self.count_per_step and acount < -self.count_per_step:
                        acount = -self.count_per_step
                    if acount > 0:
                        count = acount + 0.5
                    else:
                        count = acount - 0.5
                    if count > 0:
                        mm = 1
                    if count < 0:
                        mm = 0
                if axis == 1:
                    acount = ell * H180Motor.ELCOUNTS_PER_DEG - self.el_count
                    if self.count_per_step and acount > self.count_per_step:
                        acount = self.count_per_step
                    if self.count_per_step and acount < -self.count_per_step:
                        acount = -self.count_per_step
                    if acount > 0:
                        count = acount + 0.5
                    else:
                        count = acount - 0.5
                    if count > 0:
                        mm = 3
                    if count < 0:
                        mm = 2
                if count < 0:
                    count = -count
            if mm >= 0 and count:
                cmd_string = " move %d %d%1c" % (mm, count, 13)
                self.serial.write(cmd_string.encode("ascii"))
                resp = ""
                sleep(0.01)
                im = 0
                i = 0
                while i < 32:
                    ch = int.from_bytes(self.serial.read(1), byteorder="big")
                    sleep(0.01)
                    if i < 32:
                        resp += chr(ch)
                        i += 1
                    if ch == 13 or ch == 10:
                        break
                status = i
                sleep(0.1)
                for i in range(status):
                    if resp[i] == "M" or resp[i] == "T":
                        im = i
                ccount = int(resp[im:status].split(" ")[-1])
                if resp[im] == "M":
                    if mm == 1:
                        self.az_count += ccount
                    if mm == 0:
                        self.az_count -= ccount
                    if mm == 3:
                        self.el_count += ccount
                    if mm == 2:
                        self.el_count -= ccount
                if resp[im] == "T":
                    if mm == 1:
                        self.az_count += count
                    if mm == 0:
                        self.az_count -= count
                    if mm == 3:
                        self.el_count += count
                    if mm == 2:
                        self.el_count -= count
        if stow:
            self.az_count = 0
            self.el_count = 0

    def point(self, az, el):
        """Points an H180 Motor at a Certain Az, El

        Parameters
        ----------
        az : float
            Azimuth Coordinate to Point At
        el : float
            Elevation Coordinate to Point At

        Returns
        -------
        None
        """
        self.send_h180_cmd(az, el, False)
        return self.status()

    def status(self):
        """Requests the Current Location of the H180 Motor

        Returns
        -------
        (float, float)
            Current Azimuth and Elevation Coordinate as a Tuple of Floats
        """
        azz = self.az_count / H180Motor.AZCOUNTS_PER_DEG
        ell = self.el_count / H180Motor.ELCOUNTS_PER_DEG
        az = azz + self.az_lower_lim
        el = ell + self.el_lower_lim
        return az, el


class PushRodMotor(Motor):  # TODO: Test!
    """
    Controls old SRT PushRod Style Motors. baudrate should be 2000

    WARNING: This is currently a hard port of the azel function in sport.java, so expect some errors
    """

    AZCOUNTS_PER_DEG = (
        8.0 * 32.0 * 60.0 / (360.0 * 9.0)
    )  # Should this be 52.0 * 27.0 / 120.0?
    ELCOUNTS_PER_DEG = 52.0 * 27.0 / 120.0

    def __init__(
        self,
        port,
        baudrate,
        az_limits,
        el_limits,
        rod=(),
        counts_per_step=100,
        count_tol=1,
        count_corr=(0, 0),
    ):
        """

        Parameters
        ----------
        port : str
            Serial Port Identifier String for Communicating with the Motor
        baudrate : int
            Baudrate for serial connection
        az_limits : (float, float)
            Tuple of Lower and Upper Azimuth Limits
        el_limits : (float, float)
            Tuple of Lower and Upper Elevation Limits
        """
        Motor.__init__(self, port, baudrate, az_limits, el_limits)
        self.serial = serial.Serial(
            port=port,
            baudrate=baudrate,
            bytesize=serial.EIGHTBITS,
            parity=serial.PARITY_NONE,
            stopbits=serial.STOPBITS_ONE,
            timeout=0.1,
        )
        self.rod = rod
        self.az_count = 0.0
        self.el_count = 0.0
        self.count_per_step = counts_per_step
        self.count_tol = count_tol
        self.count_corrections = count_corr
        self.az = az_limits[0]
        self.el = el_limits[0]
        self.azatstow = 0
        self.elatstow = 0

    def send_pushrod_cmd(self, az, el, stow):
        """Sends a Command to the Pushrod Motor

        Parameters
        ----------
        az : float
            Azimuth Coordinate to Point At
        el : float
            Elevation Coordinate to Point At
        stow : bool
            Whether or Not to Stow Antenna (makes az, el irrelevant)

        Returns
        -------
        None
        """
        mm = count = 0
        lenzero = 0.0

        az = az % 360  # put az into reasonable range
        az = az + 360.0  # put az in range 180 to 540
        if az > 540.0:
            az -= 360.0
        if az < 180.0:
            az += 360.0

        region1 = region2 = region3 = 0
        if (
            self.az_limits[0] <= az < self.az_limits[1]
            and self.el_limits[0] <= el <= self.el_limits[1]
        ):
            region1 = 1
        if az > self.az_limits[0] + 180.0 and el > (180.0 - self.el_limits[1]):
            region2 = 1
        if az < self.az_limits[1] - 180.0 and el > (180.0 - self.el_limits[1]):
            region3 = 1
        if region1 == 0 and region2 == 0 and region3 == 0:
            raise ValueError("The Azimuth and Elevation Provided are Not Valid")

        flip = 0
        azz = az - self.az_limits[0]
        ell = el - self.el_limits[0]
        azscale = self.AZCOUNTS_PER_DEG
        elscale = self.ELCOUNTS_PER_DEG
        # g.set_slew(0);

        lenzero = (
            self.rod[0] * self.rod[0]
            + self.rod[1] * self.rod[1]
            - 2.0
            * self.rod[0]
            * self.rod[1]
            * cos((self.rod[3] - self.el_limits[0]) * pi / 180.0)
            - self.rod[2] * self.rod[2]
        )
        if lenzero >= 0.0:
            lenzero = sqrt(lenzero)
        else:
            lenzero = 0

        ellcount = (
            self.rod[0] * self.rod[0]
            + self.rod[1] * self.rod[1]
            - 2.0 * self.rod[0] * self.rod[1] * cos((self.rod[3] - el) * pi / 180.0)
            - self.rod[2] * self.rod[2]
        )
        if ellcount >= 0.0:
            ellcount = (-sqrt(ellcount) + lenzero) * self.rod[4]
        else:
            ellcount = 0

        if ellcount > self.el_count * 0.5:
            axis = 1
        else:
            axis = 0

        for ax in range(0, 2):
            if axis == 0:
                if azz * azscale > self.az_count * 0.5 - 0.5:
                    mm = 1
                    count = int(floor(azz * azscale - self.az_count * 0.5 + 0.5))
                if azz * azscale <= self.az_count * 0.5 + 0.5:
                    mm = 0
                    count = int(floor(self.az_count * 0.5 - azz * azscale + 0.5))
            else:
                if ellcount > self.el_count * 0.5 - 0.5:
                    mm = 3
                    count = int(floor(ellcount - self.el_count * 0.5 + 0.5))
                if ellcount <= self.el_count * 0.5 + 0.5:
                    mm = 2
                    count = int(floor(self.el_count * 0.5 - ellcount + 0.5))
            ccount = count
            if stow == 1:  # drive to stow
                count = 5000
                if axis == 0:
                    mm = 0
                    if self.azatstow == 1:
                        count = 0
                if axis == 1:
                    mm = 2  # complete azimuth motion to stow before completely drop in elevation
                    if self.elatstow == 1 or (
                        ccount <= 2.0 * self.count_per_step and self.azatstow == 0
                    ):
                        count = 0
                flip = 0
            if count > self.count_per_step and ccount > self.count_per_step:
                count = self.count_per_step
            if count >= self.count_tol:
                cmd_str = (
                    "  move " + str(mm) + " " + str(count) + "\n"
                )  # need space at start and end
                n = 0
                if count < 5000:
                    str2 = "M " + str(count) + "\n"
                else:
                    str2 = "T " + str(count) + "\n"
                recv = str2
                n = len(str2)
                j = 0
                kk = -1
                try:
                    self.serial.write(cmd_str.encode("ascii"))
                    j = n = rcount = kk = 0
                    resp = ""
                    while 0 <= kk < 3000:
                        result = self.serial.read(1)
                        if len(result) < 1:
                            j = -1
                        else:
                            j = int.from_bytes(result, byteorder="big")
                        kk += 1
                        if j >= 0 and n < 80:
                            resp += chr(j)
                            n += 1
                        if n > 0 and j == -1:
                            kk = -1  # end of message
                        # t.getTsec(g, d, gg)
                    recv = resp
                except Exception as e:
                    print(e)
                if kk != -1 or (recv[0] != "M" and recv[0] != "T"):
                    print("* ERROR comerr")
                    return  # TODO: Should throw error here?
                sleep(0.1)
                cmd_str = recv[0:n]  # String.copyValueOf(recv, 0, n - 1)
                parsed_strings = cmd_str.split(" ")
                try:
                    str2 = parsed_strings[0]
                except IndexError as e:
                    print(e)
                rcount = 0
                try:
                    srt2 = parsed_strings[1]
                    rcount = int(str2)
                except IndexError as e:
                    print(e)
                b2count = 0
                try:
                    str2 = parsed_strings[2]
                    b2count = int(str2)
                except IndexError as e:
                    print(e)
                if count < 5000:
                    fcount = (
                        count * 2 + b2count
                    )  # add extra 1 / 2 count from motor coast
                else:
                    fcount = 0
                if mm == 2 and recv[0] == "T":
                    self.elatstow = 1
                    self.el_count = 0
                if mm == 0 and recv[0] == "T":
                    self.azatstow = 1
                    self.az_count = 0
                if recv[0] == "M":
                    if axis == 0:
                        self.azatstow = 0
                        if mm == 1:
                            self.az_count += fcount
                        else:
                            self.az_count -= fcount
                    if axis == 1:
                        self.elatstow = 0
                        if mm == 3:
                            self.el_count += fcount
                        else:
                            self.el_count -= fcount
                sleep(0.005)
            axis += 1
            if axis > 1:
                axis = 0
        self.az = (
            self.az_limits[0]
            - self.count_corrections[0]
            + self.az_count * 0.5 / azscale
        )
        if self.az > 360.0:
            self.az -= 360.0
        ellnow = -self.el_count * 0.5 / self.rod[4] + lenzero
        ellnow = (
            self.rod[0] * self.rod[0]
            + self.rod[1] * self.rod[1]
            - self.rod[2] * self.rod[2]
            - ellnow * ellnow
        )
        ellnow = ellnow / (2.0 * self.rod[0] * self.rod[1])
        ellnow = -acos(ellnow) * 180.0 / pi + self.rod[3] - self.el_limits[0]
        self.el = self.el_limits[0] - self.count_corrections[1] + ellnow
        if self.el > 90.0:
            if self.az >= 180.0:
                self.az -= 180.0
            else:
                self.az += 180.0
                self.el = 180.0 - self.el

    def point(self, az, el):
        """Points an Pushrod Motor at a Certain Az, El

        Parameters
        ----------
        az : float
            Azimuth Coordinate to Point At
        el : float
            Elevation Coordinate to Point At

        Returns
        -------
        None
        """
        self.send_pushrod_cmd(az, el, 0)

    def status(self):
        """Requests the Current Location of the Pushrod Motor

        Returns
        -------
        (float, float)
            Current Azimuth and Elevation Coordinate as a Tuple of Floats
        """
        return self.az, self.el


class Arduino(Motor):
    
    '''With Arduino control to change the orientation, which input the angles, calculate the leg lengths, 
    and send to the Arduino. Later, also receive the signal from Arduino as leg lengths, and calculate the
    angle of the surface.'''
    
    # code the send function
    # add in necessary parameters [ones defined for calculation]
    # calculate the leg lengths and embed in send function
    # Code the receiving function
    # Calculate teh surface orientation according to the updated signal
    # program it into a rot2prog responsive system, make the return happens after request.
    
    '''now is :
    1. leg lengths communication
    2. Status, setting, stop are sent by computer
    3. return only happens on status and stop
    4. '''
    
    """
    Arduino bridge speaking 18-byte frames:
      b'W' + CMD + L1(5 ASCII) + L2(5) + L3(5) + END
    where L# are absolute encoder counts (0..99999) for the three actuators.
    """

    # ---- protocol constants (change END_BYTE to ord(' ') if you prefer a space) ----
    CMD_SET  = 0x2F
    CMD_STAT = 0x1F
    CMD_STOP = 0x0F
    CMD_RST  = 0x3F
    END_BYTE = ord('K')   # or: ord(' ') for ROT2-like trailing space

    # ---- encoder mapping (match your Arduino firmware) ----
    PULSES_PER_CM = 944.88189
    L0_CM         = 121.5898     # length at 0 counts
    MIN_CM        = 101.6
    MAX_CM        = 223.1898

    def __init__(self, port, baudrate, az_limits, el_limits, timeout=0.2):
        super().__init__(port, baudrate, az_limits, el_limits)
        self.serial = serial.Serial(
            port=self.port,
            baudrate=baudrate,
            bytesize=serial.EIGHTBITS,
            parity=serial.PARITY_NONE,
            stopbits=serial.STOPBITS_ONE,
            timeout=timeout,
            write_timeout=timeout,
        )

    # ------------------------- private helpers -------------------------

    @staticmethod
    def _read_exact(ser: serial.Serial, n: int) -> bytes:
        out = bytearray()
        while len(out) < n:
            chunk = ser.read(n - len(out))
            if not chunk:
                break
            out += chunk
        return bytes(out)

    @classmethod
    def _cm_to_counts(cls, cm: float) -> int:
        cnt = int(round((cm - cls.L0_CM) * cls.PULSES_PER_CM))
        return max(0, min(99999, cnt))

    @classmethod
    def _counts_to_cm(cls, cnt: int) -> float:
        return (int(cnt) / cls.PULSES_PER_CM) + cls.L0_CM

    @classmethod
    def _build_frame(cls, cmd: int, l1: int, l2: int, l3: int) -> bytes:
        assert 0 <= cmd <= 255
        for v in (l1, l2, l3):
            assert 0 <= v <= 99999
        pkt = (
            b'W' +
            bytes([cmd]) +
            f"{l1:05d}{l2:05d}{l3:05d}".encode("ascii") +
            bytes([cls.END_BYTE])
        )
        # safety
        if len(pkt) != 18:
            raise ValueError("frame not 18 bytes")
        return pkt

    @classmethod
    def _parse_frame(cls, buf: bytes):
        if len(buf) != 18 or buf[0] != 0x57 or buf[17] != cls.END_BYTE:
            raise ValueError("bad frame")
        cmd = buf[1]
        l1  = int(buf[2:7].decode("ascii"))
        l2  = int(buf[7:12].decode("ascii"))
        l3  = int(buf[12:17].decode("ascii"))
        return cmd, l1, l2, l3

    def _write_counts3(self, cmd: int, counts3):
        l1, l2, l3 = map(int, counts3)
        frame = self._build_frame(cmd, l1, l2, l3)
        self.serial.write(frame)
        self.serial.flush()

    def _read_counts3(self):
        buf = self._read_exact(self.serial, 18)
        if len(buf) != 18:
            raise TimeoutError("timeout waiting for 18-byte frame")
        return self._parse_frame(buf)

    # ----------------------- geometry glue -----------------------

    def _az_el_to_lengths_cm(self, az_deg: float, el_deg: float):
        # yaw assumed 0; your helper returns radians if out_deg=False
        p, r, _ = az_el_to_pitch_roll(az_deg, el_deg, out_deg=False)
        Lcm = calculate_leg_lengths(float(p), float(r), 0.0)   # ndarray(3,)
        return [float(x) for x in Lcm]

    def _lengths_cm_to_counts3(self, Lcm):
        return [self._cm_to_counts(x) for x in Lcm]

    def _counts3_to_lengths_cm(self, c3):
        return [self._counts_to_cm(x) for x in c3]

    # ----------------------- public API -----------------------

    def send(self, cmd: int, az: float | None = None, el: float | None = None):
        """
        Raw send. For SET, az/el must be provided and are converted to 3 absolute
        counts. For STATUS/STOP/RESET, lengths are ignored on the MCU; we send zeros.
        """
        if cmd == self.CMD_SET:
            if az is None or el is None:
                raise ValueError("SET requires az & el")
            Lcm = self._az_el_to_lengths_cm(az, el)
            counts3 = self._lengths_cm_to_counts3(Lcm)
        else:
            counts3 = (0, 0, 0)
        self._write_counts3(cmd, counts3)

    def point(self, az: float, el: float):
        """Command a move to (az, el)."""
        self.send(self.CMD_SET, az, el)
        # if you want an immediate ACK, uncomment:
        # _ = self._read_counts3()

    def status(self):
        """Request current pose; parse counts -> cm -> (az, el)."""
        self.send(self.CMD_STAT)
        cmd, l1, l2, l3 = self._read_counts3()
        # many firmwares echo CMD_STAT; accept it
        if cmd not in (self.CMD_STAT, self.CMD_STOP, self.CMD_RST):
            raise IOError(f"unexpected reply cmd=0x{cmd:02X}")

        Lcm = self._counts3_to_lengths_cm((l1, l2, l3))
        # soft range check (non-fatal)
        if (min(Lcm) < self.MIN_CM) or (max(Lcm) > self.MAX_CM):
            print("⚠️  lengths outside stroke:", Lcm)

        az, el = leg_lengths_to_az_el(np.array(Lcm, dtype=float))
        return float(az), float(el)

    def stop(self):
        """Tell the MCU to stop; reply is parsed like status."""
        self.send(self.CMD_STOP)
        cmd, l1, l2, l3 = self._read_counts3()
        Lcm = self._counts3_to_lengths_cm((l1, l2, l3))
        az, el = leg_lengths_to_az_el(np.array(Lcm, dtype=float))
        return float(az), float(el)

    def reset(self):
        """Optional homing/reset; expects a status-like reply."""
        self.send(self.CMD_RST)
        cmd, l1, l2, l3 = self._read_counts3()
        Lcm = self._counts3_to_lengths_cm((l1, l2, l3))
        az, el = leg_lengths_to_az_el(np.array(Lcm, dtype=float))
        return float(az), float(el)

# Constants
BASE_RADIUS = 61  # cm; back up, 96.61cm
TOP_RADIUS = 47  # cm
INITIAL_HEIGHT = 171.82  # cm

# -----------------------------------
# **Rotation Matrix Calculation**
# -----------------------------------
def rotation_matrix(pitch, roll, yaw):
    """Create 3D rotation matrix from pitch, roll, and yaw angles (in radians)."""
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(pitch), -np.sin(pitch)],
                   [0, np.sin(pitch), np.cos(pitch)]])
    
    Ry = np.array([[np.cos(roll), 0, np.sin(roll)],
                   [0, 1, 0],
                   [-np.sin(roll), 0, np.cos(roll)]])
    
    Rz = np.array([[np.cos(yaw), -np.sin(yaw), 0],
                   [np.sin(yaw), np.cos(yaw), 0],
                   [0, 0, 1]])
    
    return Rz @ Ry @ Rx

# -----------------------------------
# **Calculate Leg Lengths**
# -----------------------------------
# def calculate_leg_lengths(pitch, roll, yaw):
#     """Calculate leg lengths for given orientation angles."""
    
#     # Base attachment points (120 degrees apart)
#     base_points = np.array([
#         [BASE_RADIUS * np.cos(0), BASE_RADIUS * np.sin(0), 0],
#         [BASE_RADIUS * np.cos(2*np.pi/3), BASE_RADIUS * np.sin(2*np.pi/3), 0],
#         [BASE_RADIUS * np.cos(4*np.pi/3), BASE_RADIUS * np.sin(4*np.pi/3), 0]
#     ])
    
#     # Top attachment points (120 degrees apart)
#     top_points = np.array([
#         [TOP_RADIUS * np.cos(0), TOP_RADIUS * np.sin(0), INITIAL_HEIGHT],
#         [TOP_RADIUS * np.cos(2*np.pi/3), TOP_RADIUS * np.sin(2*np.pi/3), INITIAL_HEIGHT],
#         [TOP_RADIUS * np.cos(4*np.pi/3), TOP_RADIUS * np.sin(4*np.pi/3), INITIAL_HEIGHT]
#     ])
    
#     # Calculate rotation matrix
#     R = rotation_matrix(pitch, roll, yaw)
    
#     # Transform top attachment points
#     transformed_top_points = (R @ (top_points - [0, 0, INITIAL_HEIGHT]).T).T + [0, 0, INITIAL_HEIGHT]
    
#     # Calculate leg vectors
#     leg_vectors = transformed_top_points - base_points
    
#     # Calculate leg lengths
#     leg_lengths = np.linalg.norm(leg_vectors, axis=1)
    
#     if np.all(122 < leg_lengths) and np.all(leg_lengths < 223):
#         return leg_lengths
#     else:
#         raise ValueError("Calculated leg lengths are out of actuator range.")

def calculate_leg_lengths(pitch, roll, yaw):
    base_points = np.array([
        [BASE_RADIUS*np.cos(0),             BASE_RADIUS*np.sin(0),             0],
        [BASE_RADIUS*np.cos(2*np.pi/3),     BASE_RADIUS*np.sin(2*np.pi/3),     0],
        [BASE_RADIUS*np.cos(4*np.pi/3),     BASE_RADIUS*np.sin(4*np.pi/3),     0],
    ])
    top_points = np.array([
        [TOP_RADIUS*np.cos(0),              TOP_RADIUS*np.sin(0),              INITIAL_HEIGHT],
        [TOP_RADIUS*np.cos(2*np.pi/3),      TOP_RADIUS*np.sin(2*np.pi/3),      INITIAL_HEIGHT],
        [TOP_RADIUS*np.cos(4*np.pi/3),      TOP_RADIUS*np.sin(4*np.pi/3),      INITIAL_HEIGHT],
    ])
    R = rotation_matrix(pitch, roll, yaw)
    # rotate about the center at z=INITIAL_HEIGHT
    pts = (R @ (top_points - [0,0,INITIAL_HEIGHT]).T).T + [0,0,INITIAL_HEIGHT]
    legs = pts - base_points
    return np.linalg.norm(legs, axis=1)

def orientation_from_leg_lengths(measured, guess=(0,0,0)):
    """
    Solve for (pitch, roll, yaw) that makes calculate_leg_lengths()
    match the measured array of 3 lengths.
    """
    def resid(angles):
        return calculate_leg_lengths(*angles) - measured

    sol = least_squares(resid, guess)
    return sol.x  # (pitch, roll, yaw)

def leg_lengths_to_az_el(measured_lengths):
    # 1) find orientation
    pitch, roll, yaw = orientation_from_leg_lengths(measured_lengths)

    # 2) compute platform normal (local [0,0,1] rotated into world frame)
    n = rotation_matrix(pitch, roll, yaw) @ np.array([0, 0, 1])

    # 3) azimuth = angle in XY plane from X axis
    az  = np.degrees(np.arctan2(n[1], n[0]))
    # 4) elevation = angle above the XY plane
    el  = np.degrees(np.arctan2(n[2], np.linalg.norm(n[:2])))

    return az, el

def az_el_to_pitch_roll(az_deg, el_deg, out_deg=True):
    """
    Invert az/el to (pitch, roll) with yaw=0.

    Convention:
      - azimuth: bearing clockwise from north (0°=N, 90°=E)
      - elevation: angle above horizon
      - pitch: rotation about +X
      - roll:  rotation about +Z (to match your pr<->az/el pair)

    Returns:
      (pitch, roll) in degrees if out_deg=True, else radians.
    """
    # 1) convert bearing back to CCW-from-north
    a_ccw = np.deg2rad((360.0 - np.asarray(az_deg)) % 360.0)
    e     = np.deg2rad(np.asarray(el_deg))

    # 2) rebuild the platform normal from az/el
    nz = np.sin(e)
    h  = np.cos(e)
    nx = -h * np.sin(a_ccw)
    ny =  h * np.cos(a_ccw)

    # 3) invert n → (pitch, roll)
    #    ny = –sin(pitch)        → pitch = –asin(ny)
    #    nx =  sin(roll)cos(pitch), nz = cos(roll)cos(pitch) → roll = atan2(nx, nz)
    pitch = -np.arcsin(np.clip(ny, -1.0, 1.0))
    roll  =  np.arctan2(nx, nz)

    if out_deg:
        return np.rad2deg(pitch), np.rad2deg(roll),0
    else:
        return pitch, roll, 0
    
def build_rot2_len_packet(cmd_byte: int, l1: int, l2: int, l3: int, end_byte: int = 0x20) -> bytes:
    assert 0 <= cmd_byte <= 255
    for v in (l1, l2, l3):
        assert 0 <= v <= 99999
    pkt = (
        b'W' +
        bytes([cmd_byte]) +
        f"{l1:05d}{l2:05d}{l3:05d}".encode("ascii") +
        bytes([end_byte])
    )
    assert len(pkt) == 18
    return pkt

def parse_rot2_len_packet(buf: bytes):
    if len(buf) != 18 or buf[0] != 0x57:
        raise ValueError("bad frame")
    cmd = buf[1]
    l1  = int(buf[2:7].decode("ascii"))
    l2  = int(buf[7:12].decode("ascii"))
    l3  = int(buf[12:17].decode("ascii"))
    end = buf[17]
    return cmd, l1, l2, l3, end


# -----------------------------------
# **Find Position of Celestial Object**
# -----------------------------------
def celestial_obj_position(object_name):
    '''Calculate the azimuth and altitude of a celestial object.'''
    location = EarthLocation(lat=40.36 * u.deg, lon=-74.67 * u.deg, height=61.87 * u.m)
    time = Time(datetime.now())
    altaz = AltAz(obstime=time, location=location)

    try:
        obj = get_body(object_name, time, location)  # Planets, Sun, Moon
    except:
        obj = SkyCoord.from_name(object_name)  # Other celestial objects
    
    obj_altaz = obj.transform_to(altaz)
    alt, az = obj_altaz.alt.degree, obj_altaz.az.degree
    print(f"{object_name.title()}'s Alt: {alt:.2f}, Az: {az:.2f}")
    return [az, alt]

# -----------------------------------
# **Convert Celestial Coordinates to Pitch/Roll**
# -----------------------------------
def find_p_r(position):
    '''Convert celestial coordinates to pitch, roll, yaw angles.'''
    az, alt = np.radians(position)
    x = np.arctan(np.sin(az) / np.sin(alt))
    y = np.arctan(np.sin(az) / np.cos(alt))
    return [x, y, 0]

# -----------------------------------
# Calculate and genearte signal
# -----------------------------------

def calc_len(obj_name):
    '''provides the calculated lengths in cm'''
    # the obj_name can be either position or name
    if(type(obj_name) == str):
        azalt = celestial_obj_position(obj_name)
    else:
        azalt = obj_name
    
    print(azalt)
    p, r, y = find_p_r(azalt)
    
    lens = calculate_leg_lengths(p, r, y)
    return lens
    

def lconvert(list):
    '''return the string for communication'''
    num = np.int32((list/2.54-47.87)*2376.325)
    return num
    
def generate_signal(goal_azalt):
    '''generate the communication string from star orientation'''
    # calculate the length + conver to correct unit + send
    lens = calc_len(goal_azalt)
    return lconvert(lens)

# -----------------------------------
# Send Serial signal
# -----------------------------------
    
def get_serial():
    '''Get list of /dev/tty.usb* devices (Mac/Linux).'''
    addresses = subprocess.run("ls /dev/tty.usb*", shell=True, capture_output=True, text=True).stdout.splitlines()
    if len(addresses) >= 1:
        return addresses[0]  # only return first one
    else:
        raise ValueError("No Arduino found on /dev/tty.usb*")
    
def send_serial(data, port = None):
    
    if (port==None):
        port = get_serial()
        
    baud_rate = 115200
    message = generate_signal(data)

    try:
        # Open the serial port
        arduino = serial.Serial(port, baud_rate, timeout=1)
        # Give some time for the connection to establish
        time.sleep(2)
        print(f"Connected to Arduino on {port} at {baud_rate} baud.")
    except serial.SerialException as e:
        print(f"Error opening serial port: {e}")
        return

    real_len = None
    try:
        arduino.write((message + "\n").encode('utf-8'))  # Newline at end to indicate message end
        while(True):
            if arduino.in_waiting>0:
                real_len = arduino.read(arduino.in_waiting)
                break
            time.sleep(0.5)
        arduino.close()
        lens = np.int32([real_len[1:6:],real_len[7:12:],real_len[13::]])
        return lens
    except serial.SerialException as e:
        print(f"Error during write: {e}")
        
def pitch_roll_to_az_alt(pitch_deg, roll_deg):
    # Convert degrees to radians
    pitch = np.radians(pitch_deg)
    roll = np.radians(roll_deg)

    # Rotation matrices
    # Roll: rotation around Z axis (body axis)
    R_roll = np.array([
        [np.cos(roll), -np.sin(roll), 0],
        [np.sin(roll),  np.cos(roll), 0],
        [0,             0,            1]
    ])

    # Pitch: rotation around X axis (east-west axis)
    R_pitch = np.array([
        [1, 0,             0],
        [0, np.cos(pitch), -np.sin(pitch)],
        [0, np.sin(pitch),  np.cos(pitch)]
    ])

    # Combined rotation (yaw is zero, so ignored)
    R = R_roll @ R_pitch

    # Reference vector pointing to zenith (Alt=90°, Az=any)
    zenith = np.array([0, 0, 1])

    # Apply rotation
    rotated = R @ zenith

    x, y, z = rotated

    # Compute Altitude
    alt_rad = np.arcsin(z)
    alt_deg = np.degrees(alt_rad)

    # Compute Azimuth
    az_rad = np.arctan2(y, x)
    az_deg = (np.degrees(az_rad) + 360) % 360  # Normalize to 0-360

    return az_deg, alt_deg
