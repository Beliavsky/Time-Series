! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for spectral analysis.
module spectral_mod
   !! Re-export Fourier, periodogram, spectrum, and spectral-filter algorithms.
   use fourier_mod, only: real_dft, inverse_real_dft, fft_transform
   use forecast_mod, only: fourier
   use astsa_mod, only: astsa_spectrum_t, astsa_spectrum_ic_t, &
      arma_spectrum, mv_periodogram, spectrum_ic
   use garma_mod, only: garma_periodogram_t, garma_periodogram, &
      garma_semiparametric, garma_spectral_inverse
   use itsmr_mod, only: itsmr_rank_filter_t, spectral_rank_filter, itsmr_fft_smooth
   implicit none
   private

   public :: real_dft, inverse_real_dft, fft_transform, fourier
   public :: astsa_spectrum_t, astsa_spectrum_ic_t
   public :: arma_spectrum, mv_periodogram, spectrum_ic
   public :: garma_periodogram_t, garma_periodogram
   public :: garma_semiparametric, garma_spectral_inverse
   public :: itsmr_rank_filter_t, spectral_rank_filter, itsmr_fft_smooth
end module spectral_mod
