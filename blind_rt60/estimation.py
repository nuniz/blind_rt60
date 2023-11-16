from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig
from matplotlib.figure import Figure

# Constants
FRAME_LENGTH = 200e-3  # Frame length in seconds
EPS = np.finfo('float').eps


class UpdateMethod:
    NEWTON = 'newton'
    BISECTED = 'bisected'


class BlindRT60:
    def __init__(self, fs: int = 8000, framelen: Optional[float] = None, hop: Optional[float] = None,
                 percentile: int = 50., a_init: int = 0.99, sigma2_init: int = 0.5, max_itr: int = 1000,
                 max_err: float = 1e-1, a_range: tuple = (0.99, 0.999999999), bisected_itr: int = 8,
                 sigma2_range: tuple = (0., np.inf), verbose: bool = False):
        """
        Estimate the reverberation time (RT60) from the input signal.


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
        self.verbose = verbose

        self.a = None
        self.rt60 = None
        self.tau = None
        self.taus = None
        self.a_upper = None
        self.a_lower = None
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
        self.a_upper = self.a_range[1] * np.ones((batch, 1))
        self.a_lower = self.a_range[0] * np.ones((batch, 1))
        self.sigma2 = self.sigma2_init * np.ones((batch, 1))
        self.n = np.expand_dims(np.arange(self.framelen, dtype=float), axis=0)
        self.framelen_fac = self.framelen * (self.framelen - 1) / 2
        self.converged = False
        self.taus = None

    def sanity_check(self):
        """
        Check the validity of input parameters.
        """
        assert 0. <= self.percentile <= 100., 'gamma should be between 0 to 100'
        assert self.framelen > 0, f'sigma2 should be larger than 0'
        assert 0. < self.hop <= self.framelen, 'hop must be between 0 to framelen'
        assert self.a_range[0] <= self.a_init < self.a_range[
            1], f'a should be between {self.a_range[0]} to {self.a_range[1]}'
        assert self.sigma2_init > 0., f'sigma2 should be larger than 0'

    def likelihood_derivative(self, a: np.ndarray, x_frames: np.ndarray) -> (np.ndarray, np.ndarray, np.ndarray):
        """
        Calculate the first and second derivatives of the log-likelihood function with respect to 'a'.

        Parameters:
        - a: Parameter 'a' representing the decay rate
        - x_frames: Input frames

        Returns:
        - dl_da: First derivative of the log-likelihood with respect to 'a'
        - d2l_da2: Second derivative of the log-likelihood with respect to 'a'
        - sigma2: Estimated variance of the signal
        """
        a_x_prod = a ** (-2 * self.n) * x_frames ** 2
        sigma2 = np.clip(np.mean(a_x_prod, axis=1, keepdims=True), a_min=self.sigma2_range[0],
                         a_max=self.sigma2_range[1])
        dl_da = 1 / (a + EPS) * (
                1 / (sigma2 + EPS) * np.sum(self.n * a_x_prod, axis=1, keepdims=True) - self.framelen_fac)
        d2l_da2 = self.framelen_fac / (a ** 2 + EPS) + 1 / (sigma2 + EPS) * np.sum(
            (1 - 2 * self.n) * self.n * a_x_prod, axis=1, keepdims=True)
        return dl_da, d2l_da2, sigma2

    def step(self, x_frames: np.ndarray, method: str) -> np.ndarray:
        """
        Update the parameters 'a' and 'sigma2' using the specified method.

        Parameters:
        - x_frames: Input frames
        - method: Update method ('newton' or 'bisected')

        Returns:
        - dl_da: The derivative of the likelihood function with respect to 'a'
        """
        if method == UpdateMethod.NEWTON:
            dl_da, d2l_da2, self.sigma2 = self.likelihood_derivative(self.a, x_frames)
            self.a -= dl_da / (d2l_da2 + EPS)
        elif method == UpdateMethod.BISECTED:
            middle_a = 0.5 * (self.a_lower + self.a_upper)
            dl_da_upper, _, _ = self.likelihood_derivative(self.a_upper, x_frames)
            dl_da_middle, _, _ = self.likelihood_derivative(middle_a, x_frames)

            changed_sign = np.sign(dl_da_upper) != np.sign(dl_da_middle)
            not_changed_sign = np.bitwise_not(changed_sign)
            self.a = middle_a

            # No frames changed sign
            if np.sum(changed_sign) == 0:
                self.a_upper = middle_a
            # All frames changed sign
            elif np.sum(not_changed_sign) == 0:
                self.a_lower = middle_a
            # Some frames changed sign
            else:
                self.a_lower[changed_sign] = middle_a[changed_sign]
                self.a_upper[not_changed_sign] = middle_a[not_changed_sign]
        else:
            raise ValueError(f'method {method} should be {UpdateMethod.NEWTON} or {UpdateMethod.BISECTED}')

        self.a = np.clip(self.a, a_min=self.a_range[0], a_max=self.a_range[1])
        dl_da, d2l_da2, self.sigma2 = self.likelihood_derivative(self.a, x_frames)

        return dl_da

    def visualize(self, x: np.ndarray, fs: int, ylim: Tuple = (0, 1)) -> Figure:
        """
        Visualize the input signal, estimated Room Impulse Response (RT60), and its histogram.

        Parameters:
        - x: Input signal
        - fs: Sampling frequency of the input signal
        - ylim: Tuple specifying the y-axis limits for the plots (default is (0, 1))

        Returns:
        - fig: Matplotlib figure object
        """
        assert self.taus is not None
        assert self.rt60 is not None
        assert self.tau is not None
        assert np.ndim(x) == 1
        x_duration = len(x) / fs

        fig, axs = plt.subplots(nrows=1, ncols=2, width_ratios=(3, 1), sharey=True)

        # Plot the input signal normalized to its maximum value
        axs[0].plot(np.linspace(0.0, x_duration, len(x)), abs(x) / np.max(np.abs(x - np.mean(x))), label='Signal',
                    color='black')
        axs[0].set_xlabel('Time [sec]')
        axs[0].set_ylabel('Samples')
        axs[0].tick_params(axis='y', labelcolor='black')

        # Plot the estimated reverberation times for each frame
        axs0 = axs[0].twinx()
        axs0.plot(np.linspace(0.0, x_duration, len(self.taus)), self.taus, label='Tau [sec]', color='c')
        axs0.tick_params(axis='y', labelcolor='black')
        axs0.set_ylabel('Time Constant [sec]')

        # Create a histogram of taus
        bins = np.arange(ylim[0], ylim[1], 0.05)
        axs[1].hist(self.taus, bins=bins, orientation='horizontal', color='c')
        axs[1].axhline(self.tau, xmin=0.0, color='black')
        axs[1].text(2, self.tau + 0.05, f'Tau {self.tau:.2f} sec', color='black')
        axs[1].set_xlabel('Counts')

        # Set y-axis limits for all subplots
        for ax in [axs[0], axs[1], axs0]:
            ax.set_ylim(ylim)

        # Add a title to the entire visualization
        plt.suptitle(f'Blind RT60 Estimation | RT60 {self.rt60:.2f} sec')

        fig.tight_layout()
        return fig

    def estimate(self, x: np.ndarray, fs: int) -> float:
        """
        Estimate the reverberation time (RT60) from the input signal.

        Parameters:
        - x: Input signal
        - fs: Sampling frequency of the input signal

        Returns:
        - rt60: Estimated RT60
        """
        self.sanity_check()
        assert np.ndim(x) == 1

        x = sig.decimate(x, int(fs // self.fs)) if fs > self.fs else sig.resample(x, int(len(x) * self.fs / fs))
        x_frames = np.array([x[i:i + self.framelen] for i in range(0, len(x) - self.framelen + 1, self.hop)])
        self.init_states(x_frames.shape[0])

        itr = 0
        while itr < self.bisected_itr or (itr < self.max_itr and np.any(np.bitwise_not(self.converged))):
            method = UpdateMethod.BISECTED if itr < self.bisected_itr else UpdateMethod.NEWTON
            dl_da = self.step(x_frames, method=method)
            self.converged = np.abs(dl_da) <= self.max_err
            itr += 1

        self.taus = -1 / np.log(self.a) / self.fs
        self.taus[np.bitwise_not(self.converged)] = np.nan
        self.tau = np.percentile(self.taus[self.converged], q=self.percentile)
        self.rt60 = -3 * self.tau / np.log10(np.e ** -1)

        if self.verbose:
            print(f'Iteration {itr} / {self.max_itr}; rt60 {self.rt60:.2f} sec; tau {self.tau:.2f} sec')

        return self.rt60

    def __call__(self, *args, **kwargs):
        """
        Call the estimate method when the object is called.

        Returns:
        - rt60_output: The estimated RT60
        """
        return self.estimate(*args, **kwargs)
