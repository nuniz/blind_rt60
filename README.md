# BlindRT60
Python implementation of 'Blind estimation of reverberation time' [[1]](https://www.researchgate.net/publication/5923452_Blind_estimation_of_reverberation_time).

[1] Ratnam, Rama & Jones, Douglas & Wheeler, Bruce & O'Brien, William & Lansing, Charissa & Feng, Albert. (2003). Blind estimation of reverberation time. The Journal of the Acoustical Society of America. 114. 2877-92. 10.1121/1.1616578. 

For the evaluation, a speech utterance was taken from the
[NOIZEUS database](https://ecs.utdallas.edu/loizou/speech/noizeus/) [3], a repository of noisy speech corpus.

# Notes

## Model of Sound Decay
We assume that the reverberant tail of a decaying sound
y is the product of a fine structure x that is random process,
and an envelope a that is deterministic. 
$x\left[ n \right]$ is independent and identically random variables drawn from the normal distribution $N\left( {0,\sigma } \right)$.
The model for room decay then suggests that the observations y are specified by $y\left( n \right) = x\left( n \right) \cdot a\left( n \right)$.
Due to the time-varying term $a\left( n \right)$, $y\left( n \right)$ independent but not identically distributed, and their probability density function is $N\left( {0,\sigma  \cdot a\left( n \right)} \right)$.
For each estimation interval the likelihood function of y is,
$$L\left( {y;a,\sigma } \right) = \frac{1}{{\prod\limits_{n = 0}^{N - 1} {a\left( n \right)} }} \cdot {\left( {\frac{1}{{2\pi {\sigma ^2}}}} \right)^{N/2}} \cdot \exp \left( { - \frac{{\sum\limits_{n = 0}^{N - 1} {{{\left( {\frac{{y\left( n \right)}}{{a\left( n \right)}}} \right)}^2}} }}{{2{\sigma ^2}}}} \right)$$
N+1 unknown parameters of the model: $\left\{ {a\left[ {0,..,N - 1} \right],\,\,\sigma } \right\}$.
Describe $a[n]$ by damped free decay $a\left[ n \right] = \exp \left( { - \frac{n}{\tau }} \right) \buildrel \Delta \over = {a^n}$,
$$\[L\left( {y;a,\sigma } \right) = {\left( {\frac{1}{{2\pi {a^{N - 1}}{\sigma ^2}}}} \right)^{N/2}} \cdot \exp \left( { - \frac{{\sum\limits_{n = 0}^{N - 1} {{a^{ - 2n}}y{{\left( n \right)}^2}} }}{{2{\sigma ^2}}}} \right)\]$$

## Maximum Likelihood Estimator
### Equations
Given the likelihood function, the parameters $$a$$ and $$\sigma$$ can be estimated using a maximum-likelihood approach,
$$\frac{{\partial \ln L\left( {y;a,\sigma } \right)}}{{\partial a}} = {a^{ - 1}}\left( {\frac{1}{{{\sigma ^2}}}\sum\limits_{n = 0}^{N - 1} {n \cdot {a^{ - 2n}}y{{\left( n \right)}^2} - \frac{{N\left( {N - 1} \right)}}{2}} } \right)$$
$$\frac{{{\partial ^2}\ln L\left( {y;a,\sigma } \right)}}{{\partial {a^2}}} = \frac{{N\left( {N - 1} \right)}}{2}{a^{ - 2}} + \frac{1}{{{\sigma ^2}}}\sum\limits_{n = 0}^{N - 1} {n\left( {1 - 2n} \right) \cdot {a^{ - 2n}}y{{\left( n \right)}^2}} $$
$$\frac{{\partial \ln L\left( {y;a,\sigma } \right)}}{{\partial \sigma }} =  - \frac{N}{\sigma } + \frac{1}{{{\sigma ^3}}}\sum\limits_{n = 0}^{N - 1} {{a^{ - 2n}}y{{\left( n \right)}^2}} $$

* The geometric ratio is notably compressive, and in actual scenarios, the values of aa are expected to be proximate to 1. Conversely, $$\sigma$$ exhibits a broad range. 
* Examining the gradient of $\frac{{\partial \ln L\left( {y;a,\sigma } \right)}}{{\partial a}}$, initiating the process with an initial value smaller than a requires the root-solving strategy to descend the gradient fast enough.

### Solution
* Solved using numerical and iterative approach $\frac{{\partial \ln L\left( {y;a,\sigma } \right)}}{{\partial a}} = 0$; $\frac{{\partial \ln L\left( {y;a,\sigma } \right)}}{{\partial \sigma }} = 0$.
* Estimating $$a*$$:
	1. The root was bisected until the zero was bracketed.
	2. The Newton–Raphson method was applied to accurate the root, ${a_{n = 1}} = {a_n} - \frac{{\frac{{\partial \ln L\left( {y;{a_n},\sigma } \right)}}{{\partial a}}}}{{\frac{{{\partial ^2}\ln L\left( {y;{a_n},\sigma } \right)}}{{\partial {a^2}}}}}$.

* Estimating $$\sigma$$:
	$${\sigma ^2} = \frac{1}{N}\sum\limits_{n = 0}^{N - 1} {{a^{ - 2n}}y{{\left( n \right)}^2}}$$

## Strategy for Assigning the Correct Decay Rate
The model will fail during (1) estimation Frames Do Not Fall Within a Region of Free Decay, and (2) sound with a gradual rather than rapid offset. 

* In the first case, the damping of sound in a room cannot occur at a rate faster than the free decay. A robust strategy would be to select a threshold value of $$a*$$ such that the left tail of the probability density function of $$a*$$, $a = \arg \left\{ {P\left( x \right) = \gamma ;\,\,\,P\left( x \right) = \int_0^x {p\left( {{a^*}} \right)} d{a^*}} \right\}$.
* In the second case, ${p\left( {{a^*}} \right)}$ is likely to be multimodal. the strategy then is to select the first dominant peak in ${p\left( {{a^*}} \right)}$, $a = \min \arg \left\{ {dp\left( {{a^*}} \right)/d{a^*} = 0} \right\}$.
* For a unimodal symmetric distribution with $\gamma  = 0.5$ the filter will track the peak value, i.e., the median. In connected speech, where peaks cannot be clearly discriminated or the distribution is multi-modal, $$\gamma$$ should peaked based on the statistics of gap durations.


