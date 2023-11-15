import os
import unittest
import warnings

import numpy as np
import pyroomacoustics as pra
from parameterized import parameterized, param
from scipy.io import wavfile

from blind_rt60 import BlindRT60


def decaying_chirp(fs, duration: float = 5.0, frequency_start: float = 250.0, frequency_end: float = 1000.0,
                   decay_rate: float = 10) -> np.ndarray:
    """
    Generate a decaying chirp signal.

    Args:
        fs (int): Sampling frequency.
        duration (float): Duration of the signal in seconds.
        frequency_start (float): Start frequency of the chirp.
        frequency_end (float): End frequency of the chirp.
        decay_rate (float): Rate of signal decay.

    Returns:
        numpy.ndarray: Decaying chirp signal.
    """
    t = np.linspace(0, duration, int(duration * fs))

    # Create a chirp signal
    chirp_signal = np.sin(2 * np.pi * np.linspace(frequency_start, frequency_end, len(t)) * t)

    # Apply decay envelope
    decaying_chirp_signal = chirp_signal * np.exp(-decay_rate * t)

    return decaying_chirp_signal


class TestRT60(unittest.TestCase):
    """
    Test cases for the BlindRT60 class.
    """

    @parameterized.expand([
        param(fs_sig=16000, fs_estimator=8000),
        param(fs_sig=4000, fs_estimator=8000),
    ])
    def test_decimate(self, fs_sig: int = 16000, fs_estimator: int = 8000):
        """
        Test the decimate or resample methods.

        Args:
            fs_sig (int): Signal sampling frequency.
            fs_estimator (int): Estimator sampling frequency.
        """
        x = decaying_chirp(fs_sig)
        blind_rt60 = BlindRT60(fs=fs_estimator)
        rt60 = blind_rt60(x, fs_sig)
        self.assertTrue(isinstance(rt60, float))
        self.assertLess(rt60, 1)

    @parameterized.expand([
        param(decay_rate_smaller=1, decay_rate_larger=2),
        param(decay_rate_smaller=5, decay_rate_larger=6),
        param(decay_rate_smaller=10, decay_rate_larger=11),
    ])
    def test_sanity(self, decay_rate_smaller: float, decay_rate_larger: float, fs_sig: int = 8000,
                    fs_estimator: int = 8000):
        """
        Test the basic functionality of BlindRT60.

        Args:
            fs_sig (int): Signal sampling frequency.
            fs_estimator (int): Estimator sampling frequency.
        """
        x1 = decaying_chirp(fs_sig, decay_rate=decay_rate_smaller)
        x2 = decaying_chirp(fs_sig, decay_rate=decay_rate_larger)
        blind_rt60 = BlindRT60(fs=fs_estimator)
        self.assertGreater(blind_rt60(x1, fs_sig), blind_rt60(x2, fs_sig))

    @parameterized.expand([
        param(rt60_tgt=0.3, max_err=0.25),
        param(rt60_tgt=0.8, max_err=0.25),
        param(rt60_tgt=1.2, max_err=0.25),
    ])
    def test_functionality(self, rt60_tgt: float, max_err: float, path: str = r"supplementary_material/data/sp09.wav"):
        """
        Test the functionality of BlindRT60.
        Each parameter set represents a test case with different target RT60 values.

        Args:
        - rt60_tgt (float): Target reverberation time in seconds.
        - path (str, optional): Path to the source wav file. Defaults to "../supplementary_material/data/sp09.wav".
        - max_err (float, optional): Maximum allowable error between BlindRT60 estimation and Schroeder RT60.
                                    Defaults to 0.2.

        Raises:
        - AssertionError: If the absolute error between BlindRT60 estimation and Schroeder RT60
                          is greater than max_err for any microphone.
        """
        # Import a mono wavfile as the source signal
        try:
            fs, audio = wavfile.read(path)
        except Exception as e:
            warnings.warn(str(e))
            fs, audio = wavfile.read(os.path.join('..', path))

        self.assertEqual(np.ndim(audio), 1)

        # Create the room
        room_dim = [10, 7.5, 3.5]  # meters
        e_absorption, max_order = pra.inverse_sabine(rt60_tgt, room_dim)  # Invert Sabine's formula, ISM simulator

        room = pra.ShoeBox(
            room_dim, fs=fs, materials=pra.Material(e_absorption), max_order=max_order
        )

        # Place the source in the room
        room.add_source([2.5, 3.73, 1.76], signal=audio, delay=0.5)

        # Define the locations of the microphones
        mic_locs = np.c_[
            [6.3, 4.87, 1.2], [6.3, 4.93, 1.2],  # mic 1  # mic 2
        ]
        room.add_microphone_array(mic_locs)

        # Run the simulation
        room.simulate()

        # Compute Schroeder RT60
        rt60_schroeder = room.measure_rt60()

        # Compute Blind RT60
        blind_rt60 = BlindRT60()
        rt60_estimations = np.array([blind_rt60(room.mic_array.signals[i, ...], fs)
                                     for i in range(mic_locs.shape[-1])])

        # Calculate absolute error between BlindRT60 estimation and RT60
        err = np.max(np.abs(rt60_estimations - rt60_schroeder))
        self.assertLessEqual(err, max_err)


if __name__ == '__main__':
    unittest.main()
