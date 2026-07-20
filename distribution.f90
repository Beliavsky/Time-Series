! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Reusable distributions translated from rugarch with MIT numerical support.
module distribution_mod
   !! Standardized innovation densities and random generation for time-series models.
   use kind_mod, only: dp
   use random_mod, only: random_uniform, random_standard_normal, random_gamma, &
      random_standard_student, random_standard_johnson_su
   use special_functions_mod, only: log_bessel_k
   implicit none
   private

   integer, parameter, public :: distribution_normal = 1
   integer, parameter, public :: distribution_student = 2
   integer, parameter, public :: distribution_ged = 3
   integer, parameter, public :: distribution_skew_normal = 4
   integer, parameter, public :: distribution_skew_student = 5
   integer, parameter, public :: distribution_skew_ged = 6
   integer, parameter, public :: distribution_johnson_su = 7
   integer, parameter, public :: distribution_nig = 8
   integer, parameter, public :: distribution_ghyp = 9
   integer, parameter, public :: distribution_gh_skew_student = 10

   public :: standardized_log_density, random_standardized
   public :: standardized_cdf, standardized_quantile
   public :: standardized_skewness, standardized_excess_kurtosis
   public :: distribution_log_density, distribution_density
   public :: distribution_cdf, distribution_quantile, random_distribution
   public :: distribution_has_skew, distribution_has_shape
   public :: distribution_has_lambda, distribution_name

contains

   pure elemental real(dp) function distribution_log_density(value, &
      distribution, location, scale, shape, skew, lambda) result(log_density)
      !! Evaluate a location-scale innovation log density.
      real(dp), intent(in) :: value !! Variate at which to evaluate the density.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: location !! Distribution location.
      real(dp), intent(in) :: scale !! Positive distribution scale.
      real(dp), intent(in) :: shape !! Tail or shape parameter.
      real(dp), intent(in) :: skew !! Skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda parameter.

      if (scale <= 0.0_dp) then
         log_density = -huge(1.0_dp)
      else
         log_density = standardized_log_density((value - location)/scale, &
            distribution, shape, skew, lambda) - log(scale)
      end if
   end function distribution_log_density

   pure elemental real(dp) function distribution_density(value, distribution, &
      location, scale, shape, skew, lambda) result(density)
      !! Evaluate a location-scale innovation density.
      real(dp), intent(in) :: value !! Variate at which to evaluate the density.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: location !! Distribution location.
      real(dp), intent(in) :: scale !! Positive distribution scale.
      real(dp), intent(in) :: shape !! Tail or shape parameter.
      real(dp), intent(in) :: skew !! Skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda parameter.

      density = exp(distribution_log_density(value, distribution, location, &
         scale, shape, skew, lambda))
   end function distribution_density

   pure elemental real(dp) function distribution_cdf(value, distribution, &
      location, scale, shape, skew, lambda) result(probability)
      !! Evaluate a location-scale innovation CDF.
      real(dp), intent(in) :: value !! Quantile at which to evaluate the CDF.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: location !! Distribution location.
      real(dp), intent(in) :: scale !! Positive distribution scale.
      real(dp), intent(in) :: shape !! Tail or shape parameter.
      real(dp), intent(in) :: skew !! Skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda parameter.

      if (scale <= 0.0_dp) then
         probability = 0.0_dp
      else
         probability = standardized_cdf((value - location)/scale, &
            distribution, shape, skew, lambda)
      end if
   end function distribution_cdf

   pure elemental real(dp) function distribution_quantile(probability, &
      distribution, location, scale, shape, skew, lambda) result(value)
      !! Evaluate a location-scale innovation quantile.
      real(dp), intent(in) :: probability !! Probability in the unit interval.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: location !! Distribution location.
      real(dp), intent(in) :: scale !! Positive distribution scale.
      real(dp), intent(in) :: shape !! Tail or shape parameter.
      real(dp), intent(in) :: skew !! Skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda parameter.

      value = location + scale*standardized_quantile(probability, distribution, &
         shape, skew, lambda)
   end function distribution_quantile

   real(dp) function random_distribution(distribution, location, scale, shape, &
      skew, lambda) result(value)
      !! Draw one location-scale innovation.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: location !! Distribution location.
      real(dp), intent(in) :: scale !! Positive distribution scale.
      real(dp), intent(in) :: shape !! Tail or shape parameter.
      real(dp), intent(in) :: skew !! Skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda parameter.

      value = location + scale*random_standardized(distribution, shape, skew, &
         lambda)
   end function random_distribution

   pure elemental real(dp) function standardized_log_density(value, &
      distribution, shape, skew, lambda) result(log_density)
      !! Evaluate a zero-mean, unit-variance innovation log density.
      real(dp), intent(in) :: value !! Standardized variate.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: shape !! Tail or shape parameter.
      real(dp), intent(in) :: skew !! Skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda parameter.
      real(dp) :: scale, first_moment, raw_mean, raw_scale
      real(dp) :: transformed, side_scale, inverse_shape, weight, omega
      real(dp) :: normal_value, alpha, beta, delta, location

      select case (distribution)
      case (distribution_normal)
         log_density = -0.5_dp*(log(2.0_dp*acos(-1.0_dp)) + value**2)
      case (distribution_student)
         if (shape <= 2.0_dp) then
            log_density = -huge(1.0_dp)
         else
            log_density = student_log_density(value, shape)
         end if
      case (distribution_ged)
         if (shape <= 0.0_dp) then
            log_density = -huge(1.0_dp)
         else
            log_density = ged_log_density(value, shape)
         end if
      case (distribution_skew_normal)
         first_moment = sqrt(2.0_dp/acos(-1.0_dp))
         call fs_standardization(skew, first_moment, raw_mean, raw_scale)
         transformed = value*raw_scale + raw_mean
         side_scale = merge(skew, 1.0_dp/skew, transformed >= 0.0_dp)
         log_density = log(2.0_dp/(skew + 1.0_dp/skew)) + log(raw_scale) - &
            0.5_dp*log(2.0_dp*acos(-1.0_dp)) - &
            0.5_dp*(transformed/side_scale)**2
      case (distribution_skew_student)
         if (shape <= 2.0_dp .or. skew <= 0.0_dp) then
            log_density = -huge(1.0_dp)
         else
            first_moment = student_absolute_first_moment(shape)
            call fs_standardization(skew, first_moment, raw_mean, raw_scale)
            transformed = value*raw_scale + raw_mean
            side_scale = merge(skew, 1.0_dp/skew, transformed >= 0.0_dp)
            log_density = log(2.0_dp/(skew + 1.0_dp/skew)) + &
               log(raw_scale) + student_log_density(transformed/side_scale, shape)
         end if
      case (distribution_skew_ged)
         if (shape <= 0.0_dp .or. skew <= 0.0_dp) then
            log_density = -huge(1.0_dp)
         else
            scale = ged_scale(shape)
            first_moment = scale*exp(log_gamma(2.0_dp/shape) - &
               log_gamma(1.0_dp/shape))
            call fs_standardization(skew, first_moment, raw_mean, raw_scale)
            transformed = value*raw_scale + raw_mean
            side_scale = merge(skew, 1.0_dp/skew, transformed >= 0.0_dp)
            log_density = log(2.0_dp/(skew + 1.0_dp/skew)) + &
               log(raw_scale) + ged_log_density(transformed/side_scale, shape)
         end if
      case (distribution_johnson_su)
         if (shape <= 0.0_dp) then
            log_density = -huge(1.0_dp)
         else
            inverse_shape = 1.0_dp/shape
            weight = exp(min(inverse_shape**2, 50.0_dp))
            omega = -skew*inverse_shape
            raw_scale = sqrt(1.0_dp/(0.5_dp*(weight - 1.0_dp)* &
               (weight*cosh(2.0_dp*omega) + 1.0_dp)))
            transformed = (value - raw_scale*sqrt(weight)*sinh(omega))/raw_scale
            normal_value = -skew + asinh(transformed)/inverse_shape
            log_density = -log(raw_scale) - log(inverse_shape) - &
               0.5_dp*log(1.0_dp + transformed**2) - &
               0.5_dp*log(2.0_dp*acos(-1.0_dp)) - 0.5_dp*normal_value**2
         end if
      case (distribution_nig, distribution_ghyp)
         if (shape <= 0.0_dp .or. abs(skew) >= 1.0_dp) then
            log_density = -huge(1.0_dp)
         else
            call gh_parameters(skew, shape, merge(-0.5_dp, lambda, &
               distribution == distribution_nig), alpha, beta, delta, location)
            log_density = gh_log_density(value, alpha, beta, delta, location, &
               merge(-0.5_dp, lambda, distribution == distribution_nig))
         end if
      case (distribution_gh_skew_student)
         if (shape <= 4.0_dp) then
            log_density = -huge(1.0_dp)
         else if (abs(skew) <= 1.0e-12_dp) then
            log_density = student_log_density(value, shape)
         else
            call ghst_parameters(skew, shape, beta, delta, location)
            transformed = delta**2 + (value - location)**2
            log_density = 0.5_dp*(1.0_dp - shape)*log(2.0_dp) + &
               shape*log(delta) + 0.5_dp*(shape + 1.0_dp)*log(abs(beta)) + &
               log_bessel_k(0.5_dp*(shape + 1.0_dp), &
               abs(beta)*sqrt(transformed)) + beta*(value - location) - &
               log_gamma(0.5_dp*shape) - 0.5_dp*log(acos(-1.0_dp)) - &
               0.25_dp*(shape + 1.0_dp)*log(transformed)
         end if
      case default
         log_density = -huge(1.0_dp)
      end select
   end function standardized_log_density

   pure elemental real(dp) function standardized_cdf(value, distribution, &
      shape, skew, lambda) result(probability)
      !! Evaluate a standardized innovation CDF by transformed quadrature.
      real(dp), intent(in) :: value !! Standardized quantile.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: shape !! Tail or shape parameter.
      real(dp), intent(in) :: skew !! Skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda parameter.
      integer, parameter :: intervals = 1024
      real(dp) :: upper, width, angle, variate, log_integrand
      integer :: point

      if (value <= -huge(1.0_dp)**0.25_dp) then
         probability = 0.0_dp
         return
      end if
      if (value >= huge(1.0_dp)**0.25_dp) then
         probability = 1.0_dp
         return
      end if
      if (distribution == distribution_normal) then
         probability = 0.5_dp*erfc(-value/sqrt(2.0_dp))
         return
      end if
      upper = 0.5_dp + atan(value)/acos(-1.0_dp)
      width = upper/real(intervals, dp)
      probability = 0.0_dp
      do point = 1, intervals
         angle = acos(-1.0_dp)*((real(point, dp) - 0.5_dp)*width - 0.5_dp)
         variate = tan(angle)
         log_integrand = standardized_log_density(variate, distribution, &
            shape, skew, lambda) + log(acos(-1.0_dp)) - 2.0_dp*log(cos(angle))
         if (log_integrand > log(tiny(1.0_dp))) &
            probability = probability + exp(min(log_integrand, &
            log(huge(1.0_dp))))
      end do
      probability = max(0.0_dp, min(1.0_dp, probability*width))
   end function standardized_cdf

   pure elemental real(dp) function standardized_quantile(probability, &
      distribution, shape, skew, lambda) result(value)
      !! Invert a standardized innovation CDF by safeguarded bisection.
      real(dp), intent(in) :: probability !! Probability in the closed unit interval.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: shape !! Tail or shape parameter.
      real(dp), intent(in) :: skew !! Skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda parameter.
      real(dp) :: lower, upper, midpoint
      integer :: iteration

      if (probability <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      if (probability >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      lower = -2.0_dp
      upper = 2.0_dp
      do while (standardized_cdf(lower, distribution, shape, skew, lambda) > &
         probability)
         lower = 2.0_dp*lower
         if (lower < -1.0e6_dp) exit
      end do
      do while (standardized_cdf(upper, distribution, shape, skew, lambda) < &
         probability)
         upper = 2.0_dp*upper
         if (upper > 1.0e6_dp) exit
      end do
      do iteration = 1, 60
         midpoint = 0.5_dp*(lower + upper)
         if (standardized_cdf(midpoint, distribution, shape, skew, lambda) < &
            probability) then
            lower = midpoint
         else
            upper = midpoint
         end if
      end do
      value = 0.5_dp*(lower + upper)
   end function standardized_quantile

   pure elemental real(dp) function standardized_skewness(distribution, &
      shape, skew, lambda) result(value)
      !! Return standardized-distribution skewness by transformed quadrature.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: shape !! Tail or shape parameter.
      real(dp), intent(in) :: skew !! Skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda parameter.

      if ((distribution == distribution_student .or. &
         distribution == distribution_skew_student .or. &
         distribution == distribution_gh_skew_student) .and. shape <= 3.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      value = standardized_raw_moment(distribution, shape, skew, lambda, 3)
   end function standardized_skewness

   pure elemental real(dp) function standardized_excess_kurtosis(distribution, &
      shape, skew, lambda) result(value)
      !! Return standardized-distribution excess kurtosis by quadrature.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: shape !! Tail or shape parameter.
      real(dp), intent(in) :: skew !! Skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda parameter.

      if ((distribution == distribution_student .or. &
         distribution == distribution_skew_student .or. &
         distribution == distribution_gh_skew_student) .and. shape <= 4.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      value = standardized_raw_moment(distribution, shape, skew, lambda, 4) - &
         3.0_dp
   end function standardized_excess_kurtosis

   pure real(dp) function standardized_raw_moment(distribution, shape, skew, &
      lambda, order) result(moment)
      !! Integrate one raw moment over the tangent-transformed real line.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: shape !! Tail or shape parameter.
      real(dp), intent(in) :: skew !! Skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda parameter.
      integer, intent(in) :: order !! Nonnegative integer moment order.
      integer, parameter :: intervals = 4096
      real(dp) :: width, angle, variate, log_weight, contribution
      integer :: point

      width = 1.0_dp/real(intervals, dp)
      moment = 0.0_dp
      do point = 1, intervals
         angle = acos(-1.0_dp)*((real(point, dp) - 0.5_dp)*width - 0.5_dp)
         variate = tan(angle)
         log_weight = standardized_log_density(variate, distribution, shape, &
            skew, lambda) + log(acos(-1.0_dp)) - 2.0_dp*log(cos(angle))
         contribution = 0.0_dp
         if (log_weight > log(tiny(1.0_dp))) contribution = &
            variate**order*exp(min(log_weight, log(huge(1.0_dp))))
         moment = moment + contribution
      end do
      moment = moment*width
   end function standardized_raw_moment

   real(dp) function random_standardized(distribution, shape, skew, lambda) &
      result(value)
      !! Draw one zero-mean, unit-variance innovation.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: shape !! Tail or shape parameter.
      real(dp), intent(in) :: skew !! Skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda parameter.
      real(dp) :: scale, magnitude, symmetric, first_moment, raw_mean, raw_scale
      real(dp) :: probability, alpha, beta, delta, location, mixing

      select case (distribution)
      case (distribution_normal)
         value = random_standard_normal()
      case (distribution_student)
         value = random_standard_student(shape)
      case (distribution_ged)
         scale = ged_scale(shape)
         magnitude = scale*random_gamma(1.0_dp/shape)**(1.0_dp/shape)
         value = merge(magnitude, -magnitude, random_uniform() >= 0.5_dp)
      case (distribution_skew_normal, distribution_skew_student, &
         distribution_skew_ged)
         select case (distribution)
         case (distribution_skew_normal)
            symmetric = random_standard_normal()
            first_moment = sqrt(2.0_dp/acos(-1.0_dp))
         case (distribution_skew_student)
            symmetric = random_standard_student(shape)
            first_moment = student_absolute_first_moment(shape)
         case default
            scale = ged_scale(shape)
            magnitude = scale*random_gamma(1.0_dp/shape)**(1.0_dp/shape)
            symmetric = merge(magnitude, -magnitude, random_uniform() >= 0.5_dp)
            first_moment = scale*exp(log_gamma(2.0_dp/shape) - &
               log_gamma(1.0_dp/shape))
         end select
         call fs_standardization(skew, first_moment, raw_mean, raw_scale)
         probability = skew/(skew + 1.0_dp/skew)
         value = merge(abs(symmetric)*skew, -abs(symmetric)/skew, &
            random_uniform() < probability)
         value = (value - raw_mean)/raw_scale
      case (distribution_johnson_su)
         value = random_standard_johnson_su(skew, shape)
      case (distribution_nig, distribution_ghyp)
         call gh_parameters(skew, shape, merge(-0.5_dp, lambda, &
            distribution == distribution_nig), alpha, beta, delta, location)
         mixing = random_gig(merge(-0.5_dp, lambda, &
            distribution == distribution_nig), delta**2, alpha**2 - beta**2)
         value = location + beta*mixing + sqrt(mixing)*random_standard_normal()
      case (distribution_gh_skew_student)
         call ghst_parameters(skew, shape, beta, delta, location)
         mixing = 1.0_dp/random_gamma(0.5_dp*shape, 2.0_dp/delta**2)
         value = location + beta*mixing + sqrt(mixing)*random_standard_normal()
      case default
         value = 0.0_dp
      end select
   end function random_standardized

   pure elemental logical function distribution_has_skew(distribution) &
      result(has_skew)
      !! Report whether the distribution has a skew parameter.
      integer, intent(in) :: distribution !! Distribution code.

      has_skew = distribution >= distribution_skew_normal
   end function distribution_has_skew

   pure elemental logical function distribution_has_shape(distribution) &
      result(has_shape)
      !! Report whether the distribution has a shape parameter.
      integer, intent(in) :: distribution !! Distribution code.

      has_shape = distribution /= distribution_normal .and. &
         distribution /= distribution_skew_normal
   end function distribution_has_shape

   pure elemental logical function distribution_has_lambda(distribution) &
      result(has_lambda)
      !! Report whether lambda is estimated for the distribution.
      integer, intent(in) :: distribution !! Distribution code.

      has_lambda = distribution == distribution_ghyp
   end function distribution_has_lambda

   pure function distribution_name(distribution) result(name)
      !! Return the conventional rugarch distribution label.
      integer, intent(in) :: distribution !! Distribution code.
      character(len=12) :: name

      select case (distribution)
      case (distribution_normal)
         name = 'norm'
      case (distribution_student)
         name = 'std'
      case (distribution_ged)
         name = 'ged'
      case (distribution_skew_normal)
         name = 'snorm'
      case (distribution_skew_student)
         name = 'sstd'
      case (distribution_skew_ged)
         name = 'sged'
      case (distribution_johnson_su)
         name = 'jsu'
      case (distribution_nig)
         name = 'nig'
      case (distribution_ghyp)
         name = 'ghyp'
      case (distribution_gh_skew_student)
         name = 'ghst'
      case default
         name = 'unknown'
      end select
   end function distribution_name

   pure elemental real(dp) function student_log_density(value, shape) &
      result(log_density)
      !! Evaluate the standardized Student-t log density.
      real(dp), intent(in) :: value !! Standardized variate.
      real(dp), intent(in) :: shape !! Degrees of freedom.

      log_density = log_gamma(0.5_dp*(shape + 1.0_dp)) - &
         log_gamma(0.5_dp*shape) - &
         0.5_dp*log(acos(-1.0_dp)*(shape - 2.0_dp)) - &
         0.5_dp*(shape + 1.0_dp)*log(1.0_dp + value**2/(shape - 2.0_dp))
   end function student_log_density

   pure elemental real(dp) function ged_log_density(value, shape) &
      result(log_density)
      !! Evaluate the standardized generalized-error log density.
      real(dp), intent(in) :: value !! Standardized variate.
      real(dp), intent(in) :: shape !! Positive GED shape.
      real(dp) :: scale

      scale = ged_scale(shape)
      log_density = log(shape) - log(2.0_dp) - log(scale) - &
         log_gamma(1.0_dp/shape) - (abs(value)/scale)**shape
   end function ged_log_density

   pure elemental real(dp) function ged_scale(shape) result(scale)
      !! Return the unit-variance GED scale.
      real(dp), intent(in) :: shape !! Positive GED shape.

      scale = exp(0.5_dp*(log_gamma(1.0_dp/shape) - &
         log_gamma(3.0_dp/shape)))
   end function ged_scale

   pure elemental real(dp) function student_absolute_first_moment(shape) &
      result(moment)
      !! Return E(abs(X)) for a variance-standardized Student variate.
      real(dp), intent(in) :: shape !! Degrees of freedom greater than two.

      moment = exp(0.5_dp*log(shape - 2.0_dp) + &
         log_gamma(0.5_dp*(shape - 1.0_dp)) - &
         0.5_dp*log(acos(-1.0_dp)) - log_gamma(0.5_dp*shape))
   end function student_absolute_first_moment

   pure elemental subroutine fs_standardization(skew, first_moment, mean, scale)
      !! Return Fernandez-Steel raw mean and standard deviation.
      real(dp), intent(in) :: skew !! Positive skew parameter.
      real(dp), intent(in) :: first_moment !! Base absolute first moment.
      real(dp), intent(out) :: mean !! Raw transformed mean.
      real(dp), intent(out) :: scale !! Raw transformed standard deviation.
      real(dp) :: second_moment

      mean = first_moment*(skew - 1.0_dp/skew)
      second_moment = (skew**3 + skew**(-3))/(skew + 1.0_dp/skew)
      scale = sqrt(max(second_moment - mean**2, tiny(1.0_dp)))
   end subroutine fs_standardization

   pure elemental subroutine gh_parameters(rho, zeta, lambda, alpha, beta, &
      delta, location)
      !! Transform rho-zeta-lambda parameters to standardized GH parameters.
      real(dp), intent(in) :: rho !! Skew ratio in (-1,1).
      real(dp), intent(in) :: zeta !! Positive shape parameter.
      real(dp), intent(in) :: lambda !! GH index parameter.
      real(dp), intent(out) :: alpha !! GH alpha.
      real(dp), intent(out) :: beta !! GH beta.
      real(dp), intent(out) :: delta !! GH delta.
      real(dp), intent(out) :: location !! Centering location.
      real(dp) :: rho_square, kappa, delta_kappa

      rho_square = 1.0_dp - rho**2
      kappa = exp(log_bessel_k(lambda + 1.0_dp, zeta) - &
         log_bessel_k(lambda, zeta))/zeta
      delta_kappa = exp(log_bessel_k(lambda + 2.0_dp, zeta) - &
         log_bessel_k(lambda + 1.0_dp, zeta))/zeta - kappa
      alpha = zeta**2*kappa/rho_square
      alpha = sqrt(alpha*(1.0_dp + rho**2*zeta**2* &
         delta_kappa/rho_square))
      beta = alpha*rho
      delta = zeta/(alpha*sqrt(rho_square))
      location = -beta*delta**2*kappa
   end subroutine gh_parameters

   pure elemental real(dp) function gh_log_density(value, alpha, beta, delta, &
      location, lambda) result(log_density)
      !! Evaluate a generalized-hyperbolic log density.
      real(dp), intent(in) :: value !! Variate.
      real(dp), intent(in) :: alpha !! GH alpha.
      real(dp), intent(in) :: beta !! GH beta.
      real(dp), intent(in) :: delta !! GH delta.
      real(dp), intent(in) :: location !! GH location.
      real(dp), intent(in) :: lambda !! GH index.
      real(dp) :: gamma_square, distance, normalization

      gamma_square = alpha**2 - beta**2
      distance = delta**2 + (value - location)**2
      normalization = 0.5_dp*lambda*log(gamma_square) - &
         0.5_dp*log(2.0_dp*acos(-1.0_dp)) - &
         (lambda - 0.5_dp)*log(alpha) - lambda*log(delta) - &
         log_bessel_k(lambda, delta*sqrt(gamma_square))
      log_density = normalization + 0.5_dp*(lambda - 0.5_dp)*log(distance) + &
         log_bessel_k(lambda - 0.5_dp, alpha*sqrt(distance)) + &
         beta*(value - location)
   end function gh_log_density

   pure elemental subroutine ghst_parameters(skew, shape, beta, delta, location)
      !! Transform GH skew-Student parameters to zero mean and unit variance.
      real(dp), intent(in) :: skew !! GH skew parameter beta-bar.
      real(dp), intent(in) :: shape !! Degrees of freedom greater than four.
      real(dp), intent(out) :: beta !! GH beta.
      real(dp), intent(out) :: delta !! GH delta.
      real(dp), intent(out) :: location !! Centering location.

      delta = sqrt(1.0_dp/(2.0_dp*skew**2/ &
         ((shape - 2.0_dp)**2*(shape - 4.0_dp)) + 1.0_dp/(shape - 2.0_dp)))
      beta = skew/delta
      location = -beta*delta**2/(shape - 2.0_dp)
   end subroutine ghst_parameters

   real(dp) function random_gig(lambda, chi, psi) result(value)
      !! Draw a generalized inverse-Gaussian variate by ratio of uniforms.
      real(dp), intent(in) :: lambda !! GIG index.
      real(dp), intent(in) :: chi !! Nonnegative chi parameter.
      real(dp), intent(in) :: psi !! Nonnegative psi parameter.
      real(dp) :: alpha, beta, mode, reciprocal_sum, upper, lower_root
      real(dp) :: upper_root, bound_a, bound_b, constant, first, second, candidate
      integer :: iteration

      if (chi < sqrt(epsilon(1.0_dp))) then
         value = random_gamma(lambda, 2.0_dp/psi)
         return
      else if (psi < sqrt(epsilon(1.0_dp))) then
         value = 1.0_dp/random_gamma(-lambda, 2.0_dp/chi)
         return
      end if
      alpha = sqrt(psi/chi)
      beta = sqrt(psi*chi)
      mode = (lambda - 1.0_dp + sqrt((lambda - 1.0_dp)**2 + beta**2))/beta
      reciprocal_sum = mode + 1.0_dp/mode
      upper = mode
      do while (gig_root_function(upper, mode, beta, lambda) <= 0.0_dp)
         upper = 2.0_dp*upper
      end do
      lower_root = gig_bisection(sqrt(epsilon(1.0_dp)), mode, mode, beta, lambda)
      upper_root = gig_bisection(mode, upper, mode, beta, lambda)
      bound_a = (upper_root - mode)*(upper_root/mode)**(0.5_dp*(lambda - 1.0_dp))* &
         exp(-0.25_dp*beta*(upper_root + 1.0_dp/upper_root - reciprocal_sum))
      bound_b = (lower_root - mode)*(lower_root/mode)**(0.5_dp*(lambda - 1.0_dp))* &
         exp(-0.25_dp*beta*(lower_root + 1.0_dp/lower_root - reciprocal_sum))
      constant = -0.25_dp*beta*reciprocal_sum + &
         0.5_dp*(lambda - 1.0_dp)*log(mode)
      do iteration = 1, 100000
         first = max(random_uniform(), tiny(1.0_dp))
         second = random_uniform()
         candidate = mode + bound_a*second/first + &
            bound_b*(1.0_dp - second)/first
         if (candidate > 0.0_dp) then
            if (-log(first) >= -0.5_dp*(lambda - 1.0_dp)*log(candidate) + &
               0.25_dp*beta*(candidate + 1.0_dp/candidate) + constant) exit
         end if
      end do
      value = candidate/alpha
   end function random_gig

   pure real(dp) function gig_root_function(value, mode, beta, lambda) &
      result(function_value)
      !! Evaluate the GIG rejection-envelope root equation.
      real(dp), intent(in) :: value !! Positive candidate.
      real(dp), intent(in) :: mode !! Standardized GIG mode.
      real(dp), intent(in) :: beta !! Standardized GIG beta.
      real(dp), intent(in) :: lambda !! GIG index.

      function_value = 0.5_dp*beta*value**3 - value**2* &
         (0.5_dp*beta*mode + lambda + 1.0_dp) + value* &
         ((lambda - 1.0_dp)*mode - 0.5_dp*beta) + 0.5_dp*beta*mode
   end function gig_root_function

   pure real(dp) function gig_bisection(lower, upper, mode, beta, lambda) &
      result(root)
      !! Find a bracketed GIG envelope root by bisection.
      real(dp), intent(in) :: lower !! Lower bracket.
      real(dp), intent(in) :: upper !! Upper bracket.
      real(dp), intent(in) :: mode !! Standardized GIG mode.
      real(dp), intent(in) :: beta !! Standardized GIG beta.
      real(dp), intent(in) :: lambda !! GIG index.
      real(dp) :: left, right, middle, left_value, middle_value
      integer :: iteration

      left = lower
      right = upper
      left_value = gig_root_function(left, mode, beta, lambda)
      do iteration = 1, 100
         middle = 0.5_dp*(left + right)
         middle_value = gig_root_function(middle, mode, beta, lambda)
         if (left_value*middle_value <= 0.0_dp) then
            right = middle
         else
            left = middle
            left_value = middle_value
         end if
      end do
      root = 0.5_dp*(left + right)
   end function gig_bisection

end module distribution_mod
