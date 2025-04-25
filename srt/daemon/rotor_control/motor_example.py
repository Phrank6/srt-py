import serial

from srt.daemon.rotor_control.pos import *
from abc import ABC, abstractmethod
from time import sleep
from math import cos, acos, pi, sqrt, floor


class RPS3Motor(Motor):
    """Class for Controlling any ROT2 Protocol-Supporting Motor (e.g. SPID Motors)

    See Also
    --------
    <http://ryeng.name/blog/3>
    <https://github.com/jaidenfe/rot2proG/blob/master/rot2proG.py>
    <https://www.haystack.mit.edu/edu/undergrad/srt/pdf%20files/MD-01%20en.pdf>
    """

    VALID_PULSES_PER_2MM = (188, 189, 190)

    def __init__(
        self,
        port,
        baudrate,
        az_limits,
        el_limits,
        pulses_per_2mm=189,
        test_pulses_per_2mm=True,
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
        if pulses_per_2mm in RPS3Motor.VALID_PULSES_PER_2MM:
            self.pulses_per_2mm = pulses_per_2mm
        else:
            raise ValueError("Invalid Pulse Per Degree Value")
        if test_pulses_per_2mm:
            self.status()

    def send_rps_pkt(self, cmd, az=None, el=None):
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
        
        *** one time writing, therefore, we should write a longer series with the lengths specified

        Returns
        -------
        None
        """
        if az is not None and el is not None:
            pry = generate_signal((az,el))
        else:
            pry = generate_signal((0,0))

        leg1_ticks = (
            self.pulses_per_2mm
        )  # Documentation for Rot2 Says This Is Ignored
        leg2_ticks = (
            self.pulses_per_2mm
        )  # Documentation for Rot2 Says This Is Ignored
        leg3_ticks = (
            self.pulses_per_2mm
        )  # Documentation for Rot2 Says This Is Ignored
        
        cmd_string = "w%c%05d%03d%05d%03d%05d%03d" % (
            cmd,
            pry[0],
            leg1_ticks,
            pry[1],
            leg2_ticks,
            pry[2],
            leg3_ticks,
        ) # make the sending format "aaaaabbbaaaaabbbaaaaabbbk"
        cmd_bytes = cmd_string.encode("ascii")
        print("Packet of Size " + str(len(cmd_bytes)))
        print([hex(val) for val in cmd_bytes])
        self.serial.write(cmd_bytes)

    def receive_rps_pkt(self):
        """Receives and Parsers an ROT2 Status Packet
        *** reads the signal received and convert to az/alt; can embed the az/alt calculation function from 3 linear actuator to this

        Returns
        -------
        (float, float)
            Azimuth and Elevation Coordinate as a Tuple of Floats
        """
        received_vals = self.serial.read(24)
        pry = [received_vals[0]*10000+received_vals[1]*1000+received_vals[2]*100+received_vals[3]*10+received_vals[4],
               received_vals[8]*10000+received_vals[9]*1000+received_vals[10]*100+received_vals[11]*10+received_vals[12],
               received_vals[16]*10000+received_vals[17]*1000+received_vals[18]*100+received_vals[19]*10+received_vals[20]]
        az , el = leg_lengths_to_az_el(pry)
        leg1_pulse_per_mm = received_vals[5]*100+received_vals[6]*10+received_vals[7]
        leg2_pulse_per_mm = received_vals[13]*100+received_vals[14]*10+received_vals[15]
        leg3_pulse_per_mm = received_vals[21]*100+received_vals[22]*10+received_vals[23]

        assert leg1_pulse_per_mm == leg2_pulse_per_mm == leg3_pulse_per_mm  # Consistency Check
        if leg1_pulse_per_mm != self.pulses_per_degree:
            print(
                "Motor Pulses Per Degree Incorrect, Changing Value to "
                + str(leg1_pulse_per_mm)
            )
            self.pulses_per_degree = leg1_pulse_per_mm
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
        cmd = "l"  # Rot2 Set Command
        az_relative = az - self.az_limits[0]
        el_relative = el - self.el_limits[0]
        self.send_rps_pkt(cmd, az=az_relative, el=el_relative)

    def status(self):
        """Requests the Current Location of the ROT2 Motor

        Returns
        -------
        (float, float)
            Current Azimuth and Elevation Coordinate as a Tuple of Floats
        """
        cmd = "s".encode()  # Rot2 Status Command
        self.send_rps_pkt(cmd)
        # az_relative, el_relative = self.receive_rps_pkt()
        # return az_relative + self.az_limits[0], el_relative + self.el_limits[0]

    def stop(self):
        """Stops the ROT2 Motor at its Current Location

        Returns
        -------
        None
        """
        cmd = "t".encode()  # make the current position the same as predicted
        self.send_rps_pkt(cmd) 
        # az_relative, el_relative = self.receive_rot2_pkt()
        # return (az_relative + self.az_limits[0], el_relative + self.el_limits[0])
