from typing import Optional
import numpy as np

# Constants
FRAME_LENGTH = 200e-3  # Frame length in seconds


class UpdateMethod:
    NEWTON = 'newton'
    BISECTED = 'bisected'


class BlindRT60:
    def __init__(self, fs: int, framelen: Optional[float] = None, hop: Optional[float] = None,
                 percentile: int = 50., a_init: int = 0.9, sigma2_init: int = 0.5, max_itr: int = 1000,
                 max_err: float = 1e-1, a_range: tuple = (0., 0.99999999), bisected_itr: int = 8,
                 sigma2_range: tuple = (0., 1.)):
        """
        Initialize the BlindRT60 object.

        Parameters:
        - fs: Sample rate
        - framelen: Length of each frame in seconds
        - hop: Hop size in seconds
        - percentile: Pre specified percentile value
        - a_init: Initial value for the parameter 'a'
        - sigma2_init: Initial value for the parameter 'sigma2'
        - max_itr: Maximum number of iterations
        - max_err: Maximum error for convergence
        - a_range: Range of valid values for 'a'
        - bisected_itr: Number of iterations for the bisection method
        - sigma2_range: Range of valid values for 'sigma2'
        """
        self.fs = fs
        self.framelen = int(self.fs * FRAME_LENGTH) if framelen is None else int(self.fs * framelen)
        self.hop = int(self.framelen) // 4 if hop is None else int(self.fs * hop)
        self.percentile = percentile
        self.a_init = a_init
        self.sigma2_init = sigma2_init
        self.max_itr = max_itr
        self.max_err = max_err
        self.a_range = a_range
        self.bisected_itr = bisected_itr
        self.sigma2_range = sigma2_range

        self.a = None
        self.a_bisected = None
        self.sigma2 = None
        self.n = None
        self.framelen_fac = None
        self.converged = None

        self.sanity_check()

    def init_states(self, batch):
        """
        Initialize the internal states of the estimator.

        Parameters:
        - batch: Number of frames in the input signal
        """
        self.a = self.a_init * np.ones((batch, 1))
        self.a_bisected = self.a_range[1] * np.ones((batch, 1))
        self.sigma2 = self.sigma2_init * np.ones((batch, 1))
        self.n = np.expand_dims(np.arange(self.framelen), axis=0)
        self.framelen_fac = self.framelen * (self.framelen - 1) / 2
        self.converged = False

    def sanity_check(self):
        """
        Check the validity of input parameters.
        """
        assert 0. <= self.percentile <= 100., 'gamma should be between 0 to 100'
        assert self.framelen > 0, f'sigma2 should be larger than 0'
        assert 0 < self.hop <= self.framelen, 'hop must be between 0 to framelen'
        assert self.a_range[0] <= self.a_init < self.a_range[
            1], f'a should be between {self.a_range[0]} to {self.a_range[1]}'
        assert self.sigma2_init > 0, f'sigma2 should be larger than 0'

    def a_x_prod(self, x_frames: np.ndarray, bisected: bool = False):
        """
        Calculate the product of 'a' and squared input frames.

        Parameters:
        - x_frames: Input frames
        - bisected: Flag indicating whether to use bisected 'a' values
        """
        _a = self.a_bisected if bisected else self.a
        a_x_prod = _a ** (-2 * self.n) * x_frames ** 2
        return a_x_prod

    def step(self, x_frames: np.ndarray, method: str):
        """
        Update the parameters 'a' and 'sigma2' using the specified method.

        Parameters:
        - x_frames: Input frames
        - method: Update method ('newton' or 'bisected')
        """
        self.sigma2 = np.mean(self.a_x_prod(x_frames), axis=1, keepdims=True)
        self.sigma2 = np.clip(self.sigma2, a_min=self.sigma2_range[0], a_max=self.sigma2_range[1])
        dl_da = 1 / self.a * (1 / self.sigma2 * np.sum(self.n * self.a_x_prod(x_frames), axis=1,
                                                       keepdims=True) - self.framelen_fac)

        if method == UpdateMethod.NEWTON:
            d2l_da2 = self.framelen_fac / (self.a ** 2) + 1 / self.sigma2 * np.sum(
                (1 - 2 * self.n) * self.n * self.a_x_prod(x_frames), axis=1, keepdims=True)
            self.a -= dl_da / d2l_da2
        elif method == UpdateMethod.BISECTED:
            sigma2_bisected = np.mean(self.a_x_prod(x_frames, bisected=True), axis=1, keepdims=True)
            dl_da_bisected = 1 / self.a_bisected * (
                    1 / sigma2_bisected * np.sum(self.n * self.a_x_prod(x_frames, bisected=True), axis=1,
                                                 keepdims=True) - self.framelen_fac)

            changed_sign = np.sign(dl_da) != np.sign(dl_da_bisected)
            not_changed_sign = np.bitwise_not(changed_sign)

            middle_a = 0.5 * (self.a + self.a_bisected)
            if changed_sign.sum() == 0:
                self.a = middle_a
            elif not_changed_sign.sum() == 0:
                self.a_bisected = middle_a
            else:
                self.a_bisected[changed_sign] = middle_a[changed_sign]
                self.a[np.bitwise_not(not_changed_sign)] = middle_a[np.bitwise_not(not_changed_sign)]
        else:
            raise ValueError(f'method {method} should be {UpdateMethod.NEWTON} or {UpdateMethod.BISECTED}')

        self.a = np.clip(self.a, a_min=self.a_range[0], a_max=self.a_range[1])
        updated_dl_da = 1 / self.a * (
                1 / self.sigma2 * np.sum(self.n * self.a_x_prod(x_frames), axis=1, keepdims=True) - self.framelen_fac)
        return updated_dl_da

    def estimate(self, x: np.ndarray):
        """
        Estimate the Room Impulse Response (RT60) from the input signal.

        Parameters:
        - x: Input signal
        """
        self.sanity_check()
        x_frames = np.array([x[i:i + self.framelen] for i in range(0, len(x) - self.framelen + 1, self.hop)])
        self.init_states(x_frames.shape[0])

        itr = 0

        while itr < self.bisected_itr or (itr < self.max_itr and np.any(np.bitwise_not(self.converged))):
            method = UpdateMethod.BISECTED if itr <= self.bisected_itr else UpdateMethod.NEWTON
            dl_da = self.step(x_frames, method=method)
            self.converged = np.abs(dl_da) <= self.max_err
            itr += 1

        taus = -1 / np.log(self.a) / self.fs
        tau = np.percentile(taus[self.converged], q=self.percentile)
        rt60 = -3 * tau / np.log10(np.e ** -1)
        return rt60

    def __call__(self, *args, **kwargs):
        """
        Call the estimate method when the object is called.
        """
        return self.estimate(*args, **kwargs)
