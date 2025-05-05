import numpy as np
from astropy.coordinates import AltAz, EarthLocation, get_body, SkyCoord
from astropy import units as u
from astropy.time import Time
from datetime import datetime
import subprocess
import serial
import time
from scipy.optimize import least_squares


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


# -----------------------------------
# **Test Code**
# -----------------------------------
if __name__ == "__main__":
    
    # port = ac.get_serial()
    # ser = ac.start_serial(port)
    # print("Moving to test position...")
    # move_to_point([pitch, roll, yaw], np.radians([0, 0, 0]), 0, step_time=3, save=True, obs_name="yirenyinjiuzui")
    # array = pd.read_csv('orientaiton_files/test.csv').to_numpy()

    pitch, roll, yaw = np.radians([55, 10, 0])
    print(calculate_leg_lengths(pitch, roll, yaw))

    p, r, y = (find_p_r(celestial_obj_position('capella')))
    print(lconvert(calculate_leg_lengths(p,r,y)))
    
    print(generate_signal('capella'))
    
    print(pitch_roll_to_az_alt(38, 0))