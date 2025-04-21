import serial
import time

class RPS3Controller:
    """
    Talk to your Arduino-based 3RPS mock:

      • Send goal pulses with: lddddd,ddddd,ddddd
      • Receive current pulses as: ddddd,ddddd,ddddd
    """

    def __init__(self, port, baudrate=115200, timeout=1.0):
        self.ser = serial.Serial(port, baudrate, timeout=timeout)
        # give the Arduino time to reboot/reset
        time.sleep(2)

    def _format_cmd(self, p1, p2, p3):
        """Format three ints into 'l00000,00000,00000\\n'."""
        return f"l{p1:05d},{p2:05d},{p3:05d}\n".encode("ascii")

    def _parse_reply(self, line):
        """Parse '00000,00000,00000' → (0, 0, 0)."""
        parts = line.strip().split(",")
        if len(parts) != 3:
            raise ValueError(f"Bad reply from Arduino: {line!r}")
        return tuple(int(x) for x in parts)

    def send_goals(self, g1, g2, g3):
        """
        Send three goal‐pulse values to the Arduino.
        """
        cmd = self._format_cmd(g1, g2, g3)
        self.ser.write(cmd)

    def read_positions(self):
        """
        Read one line back from the Arduino and parse it.
        Blocks until timeout if nothing arrives.
        """
        raw = self.ser.readline().decode("ascii")
        if not raw:
            raise TimeoutError("No response from Arduino")
        return self._parse_reply(raw)

    def move_and_wait(self, g1, g2, g3, settle_time=0.1):
        """
        Send goals, wait a bit, then read back current pulses.
        """
        self.send_goals(g1, g2, g3)
        time.sleep(settle_time)
        return self.read_positions()

    def close(self):
        self.ser.close()


if __name__ == "__main__":
    # Example usage:
    ctl = RPS3Controller(port="/dev/ttyACM0")  # adjust to your port

    try:
        # Move to (10000, 20000, 15000) pulses
        current = ctl.move_and_wait(10000, 20000, 15000)
        print("At pulses:", current)

        # Move to (0, 0, 0)
        current = ctl.move_and_wait(0, 0, 0)
        print("Back home:", current)

    finally:
        ctl.close()
