! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Shared special functions adapted from translated GPL time-series packages.
module special_functions_mod
   use kind_mod, only: dp
   use linalg_mod, only: inverse_logdet
   implicit none
   private

   public :: regularized_gamma_q
   public :: regularized_beta
   public :: bessel_k01, log_bessel_k
   public :: normal_log_density
   public :: multivariate_normal_log_density

   interface normal_log_density
      module procedure centered_normal_log_density
      module procedure shifted_normal_log_density
   end interface normal_log_density

contains

   pure elemental subroutine bessel_k01(value, log_k1, ratio)
      !! Approximate log(K1(x)) and K0(x)/K1(x) for positive x.
      real(dp), intent(in) :: value !! Positive Bessel argument.
      real(dp), intent(out) :: log_k1 !! Natural logarithm of K1.
      real(dp), intent(out) :: ratio !! Ratio K0/K1.
      real(dp) :: scaled, polynomial, k0, k1, i0, i1, log_half

      if (value > 600.0_dp) then
         scaled = 1.0_dp/value
         log_k1 = 0.5_dp*(log(acos(-1.0_dp)) - log(2.0_dp) - log(value)) - &
            value + log(1.0_dp + scaled*(0.375_dp + &
            scaled*(-0.1171875_dp + scaled*0.1025390625_dp)))
         ratio = 1.0_dp - 0.5_dp*scaled
      else if (value <= 2.0_dp) then
         scaled = (value/3.75_dp)**2
         i0 = 1.0_dp + scaled*(3.5156229_dp + scaled*(3.0899424_dp + &
            scaled*(1.2067492_dp + scaled*(0.2659732_dp + &
            scaled*(0.0360768_dp + scaled*0.0045813_dp)))))
         i1 = value*(0.5_dp + scaled*(0.87890594_dp + &
            scaled*(0.51498869_dp + scaled*(0.15084934_dp + &
            scaled*(0.02658733_dp + scaled*(0.00301532_dp + &
            scaled*0.00032411_dp))))))
         log_half = log(0.5_dp*value)
         scaled = (0.5_dp*value)**2
         k0 = -log_half*i0 + (-0.57721566_dp + scaled*(0.42278420_dp + &
            scaled*(0.23069756_dp + scaled*(0.03488590_dp + &
            scaled*(0.00262698_dp + scaled*(0.00010750_dp + &
            scaled*0.00000740_dp))))))
         k1 = log_half*i1 + (1.0_dp/value)*(1.0_dp + &
            scaled*(0.15443144_dp + scaled*(-0.67278579_dp + &
            scaled*(-0.18156897_dp + scaled*(-0.01919402_dp + &
            scaled*(-0.00110404_dp + scaled*(-0.00004686_dp)))))))
         log_k1 = log(k1)
         ratio = k0/k1
      else
         scaled = 2.0_dp/value
         polynomial = exp(-value)/sqrt(value)
         k0 = polynomial*(1.25331414_dp + scaled*(-0.07832358_dp + &
            scaled*(0.02189568_dp + scaled*(-0.01062446_dp + &
            scaled*(0.00587872_dp + scaled*(-0.00251540_dp + &
            scaled*0.00053208_dp))))))
         k1 = polynomial*(1.25331414_dp + scaled*(0.23498619_dp + &
            scaled*(-0.03655620_dp + scaled*(0.01504268_dp + &
            scaled*(-0.00780353_dp + scaled*(0.00325614_dp + &
            scaled*(-0.00068245_dp)))))))
         log_k1 = log(k1)
         ratio = k0/k1
      end if
   end subroutine bessel_k01

   pure real(dp) function log_bessel_k(order, value) result(log_value)
      !! Approximate log(K_order(value)) for nonnegative real order.
      real(dp), intent(in) :: order !! Bessel order; its absolute value is used.
      real(dp), intent(in) :: value !! Positive Bessel argument.
      real(dp) :: log_k0, log_k1, ratio, previous, current, following, fraction
      real(dp) :: positive_order
      integer :: integer_order, index

      positive_order = abs(order)
      call bessel_k01(value, log_k1, ratio)
      log_k0 = log_k1 + log(ratio)
      if (positive_order < 1.0_dp) then
         log_value = (1.0_dp - positive_order)*log_k0 + &
            positive_order*log_k1
         return
      end if
      integer_order = int(positive_order)
      fraction = positive_order - real(integer_order, dp)
      previous = exp(log_k0)
      current = exp(log_k1)
      do index = 1, integer_order - 1
         following = 2.0_dp*real(index, dp)*current/value + previous
         previous = current
         current = following
      end do
      if (fraction < 1.0e-14_dp) then
         log_value = log(current)
      else
         following = 2.0_dp*real(integer_order, dp)*current/value + previous
         log_value = (1.0_dp - fraction)*log(current) + fraction*log(following)
      end if
   end function log_bessel_k

   pure real(dp) function regularized_beta(value, first_shape, second_shape) &
      result(probability)
      !! Return the regularized incomplete beta function.
      real(dp), intent(in) :: value !! Evaluation point on the unit interval.
      real(dp), intent(in) :: first_shape !! Positive first beta shape.
      real(dp), intent(in) :: second_shape !! Positive second beta shape.
      real(dp) :: factor

      if (first_shape <= 0.0_dp .or. second_shape <= 0.0_dp) then
         probability = 0.0_dp
         return
      end if
      if (value <= 0.0_dp) then
         probability = 0.0_dp
         return
      end if
      if (value >= 1.0_dp) then
         probability = 1.0_dp
         return
      end if
      factor = exp(log_gamma(first_shape + second_shape) - &
         log_gamma(first_shape) - log_gamma(second_shape) + &
         first_shape*log(value) + second_shape*log(1.0_dp - value))
      if (value < (first_shape + 1.0_dp)/(first_shape + second_shape + &
         2.0_dp)) then
         probability = factor*beta_fraction(value, first_shape, second_shape)/ &
            first_shape
      else
         probability = 1.0_dp - factor*beta_fraction(1.0_dp - value, &
            second_shape, first_shape)/second_shape
      end if
      probability = max(0.0_dp, min(1.0_dp, probability))
   end function regularized_beta

   pure real(dp) function beta_fraction(value, first_shape, second_shape) &
      result(fraction)
      !! Evaluate the incomplete-beta continued fraction.
      real(dp), intent(in) :: value !! Evaluation point on the unit interval.
      real(dp), intent(in) :: first_shape !! Positive first beta shape.
      real(dp), intent(in) :: second_shape !! Positive second beta shape.
      real(dp) :: qab, qap, qam, c, d, h, aa, delta
      integer :: iteration, twice

      qab = first_shape + second_shape
      qap = first_shape + 1.0_dp
      qam = first_shape - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab*value/qap
      if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
      d = 1.0_dp/d
      h = d
      do iteration = 1, 10000
         twice = 2*iteration
         aa = real(iteration, dp)*(second_shape - real(iteration, dp))*value/ &
            ((qam + real(twice, dp))*(first_shape + real(twice, dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
         c = 1.0_dp + aa/c
         if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
         d = 1.0_dp/d
         h = h*d*c
         aa = -(first_shape + real(iteration, dp))* &
            (qab + real(iteration, dp))*value/ &
            ((first_shape + real(twice, dp))*(qap + real(twice, dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
         c = 1.0_dp + aa/c
         if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
         d = 1.0_dp/d
         delta = d*c
         h = h*delta
         if (abs(delta - 1.0_dp) <= 10.0_dp*epsilon(1.0_dp)) exit
      end do
      fraction = h
   end function beta_fraction

   pure real(dp) function multivariate_normal_log_density(value, mean, &
      covariance) result(log_density)
      !! Evaluate a nonsingular multivariate Gaussian log density.
      real(dp), intent(in) :: value(:) !! Observation vector.
      real(dp), intent(in) :: mean(:) !! Mean vector.
      real(dp), intent(in) :: covariance(:, :) !! Covariance matrix.
      real(dp), allocatable :: inverse(:, :), difference(:)
      real(dp) :: log_determinant
      integer :: info

      if (size(value) < 1 .or. size(value) /= size(mean) .or. &
         any(shape(covariance) /= [size(value), size(value)])) then
         log_density = -huge(1.0_dp)
         return
      end if
      allocate(inverse(size(value), size(value)))
      call inverse_logdet(covariance, inverse, log_determinant, info, &
         1.0e-12_dp)
      if (info /= 0) then
         log_density = -huge(1.0_dp)
         return
      end if
      difference = value - mean
      log_density = -0.5_dp*(real(size(value), dp)* &
         log(2.0_dp*acos(-1.0_dp)) + log_determinant + &
         dot_product(difference, matmul(inverse, difference)))
   end function multivariate_normal_log_density

   pure elemental real(dp) function centered_normal_log_density(value, &
      standard_deviation) result(log_density)
      !! Evaluate a centered Gaussian log density.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: standard_deviation !! Standard deviation.

      log_density = shifted_normal_log_density(value, 0.0_dp, standard_deviation)
   end function centered_normal_log_density

   pure elemental real(dp) function shifted_normal_log_density(value, mean, &
      standard_deviation) result(log_density)
      !! Evaluate a Gaussian log density with a specified mean.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: mean !! Mean value or vector.
      real(dp), intent(in) :: standard_deviation !! Standard deviation.

      if (standard_deviation <= 0.0_dp) then
         log_density = -huge(1.0_dp)
      else
         log_density = -0.5_dp*log(2.0_dp*acos(-1.0_dp)) - &
            log(standard_deviation) - &
            0.5_dp*((value - mean)/standard_deviation)**2
      end if
   end function shifted_normal_log_density

   pure real(dp) function regularized_gamma_q(shape, value) result(probability)
      !! Return the complemented regularized incomplete gamma ratio.
      real(dp), intent(in) :: shape !! Shape.
      real(dp), intent(in) :: value !! Input value.
      real(dp) :: term, total, ap, b, c, d, delta, h, log_scale
      integer :: iteration

      if (shape <= 0.0_dp .or. value < 0.0_dp) then
         probability = 0.0_dp
         return
      end if
      if (value == 0.0_dp) then
         probability = 1.0_dp
         return
      end if
      log_scale = -value + shape*log(value) - log_gamma(shape)
      if (value < shape + 1.0_dp) then
         ap = shape
         term = 1.0_dp/shape
         total = term
         do iteration = 1, 10000
            ap = ap + 1.0_dp
            term = term*value/ap
            total = total + term
            if (abs(term) <= epsilon(1.0_dp)*abs(total)) exit
         end do
         probability = 1.0_dp - total*exp(log_scale)
      else
         b = value + 1.0_dp - shape
         c = 1.0_dp/tiny(1.0_dp)
         d = 1.0_dp/b
         h = d
         do iteration = 1, 10000
            term = -real(iteration, dp)*(real(iteration, dp) - shape)
            b = b + 2.0_dp
            d = term*d + b
            if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
            c = b + term/c
            if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
            d = 1.0_dp/d
            delta = d*c
            h = h*delta
            if (abs(delta - 1.0_dp) <= 10.0_dp*epsilon(1.0_dp)) exit
         end do
         probability = h*exp(log_scale)
      end if
      probability = max(0.0_dp, min(1.0_dp, probability))
   end function regularized_gamma_q

end module special_functions_mod
