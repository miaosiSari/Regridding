# Regridding

## Code Structure

Replicate "Regridding reconstruction algorithm for real-time tomographic imaging".

The code starts with $\texttt{gridrec_radon.m}$. The $\texttt{gridrec.m}$ is now depricated.

The program $\texttt{test_fft_slice.m}$ is a verifaction of the [Fourier Slice Theorem](https://en.wikipedia.org/wiki/Projection-slice_theorem)  and is not called in  $\texttt{gridrec_radon.m}$ or $\texttt{gridrec.m}$ .

$\texttt{calc_psnr_ssim.m}$: It calculates the PSNR and SSIM values between two images. It calls $\texttt{metrix_mse.m}$, $\texttt{metrix_psnr.m}$, $\texttt{metrix_ssim.m}$ and $\texttt{ssim_index.m}$.

## Update Logs

Update Aug 23rd:

We successfully implement $\texttt{gridrec_radon.m}$.. Different from "gridrec.m", it uses Matlab function $\texttt{radon}$ instead of $\texttt{imrotate}$ to do projection. Function $\texttt{radon}$ employs more parallel beams to perform projection, therefore could considerably reduce the artifacts.

In addition, we test the gridrec algorithm on natural images from the GOPRO dataset which is famous for i deblurring.

Currently the problems are:

1. The gridrec algorithm generates blurry results when the resolution is low.
2. The algorithm is not well parallelized (TODO).

## Experimental Results

### Verification of the Fourier Slice Theorem

The high frequency components (near the start and the end) suffers from great loss if the rotation angle is not close to 0 or $\frac{\pi}{2}$.

### Verification of the Gridrec Algorithm

Gridrec on the GOPRO image:

PSNR = 32.8999, SSIM = 0.9649.

Gridrec on the Shepp-Logan:

PSNR = 32.2522, SSIM = 0.9646.