! SPDX-License-Identifier: MIT
! SPDX-FileComment: Generic draw-driven resampling utilities.
module resampling_mod
   !! Provide deterministic kernels for caller-supplied bootstrap draws.
   use kind_mod, only: dp
   implicit none
   private

   interface resample
      module procedure resample_vector
      module procedure resample_matrix
   end interface resample

   interface block_resample
      module procedure block_resample_vector
      module procedure block_resample_matrix
   end interface block_resample

   interface additive_resample
      module procedure additive_resample_vector
      module procedure additive_resample_matrix
   end interface additive_resample

   interface wild_resample
      module procedure wild_resample_vector
      module procedure wild_resample_matrix
   end interface wild_resample

   public :: resample, block_resample, additive_resample, wild_resample

contains

   pure function resample_vector(values, indices) result(sample)
      !! Select vector entries using caller-supplied indices with replacement.
      real(dp), intent(in) :: values(:) !! Values to resample.
      integer, intent(in) :: indices(:) !! Source index for every output value.
      real(dp), allocatable :: sample(:)

      if (size(values) == 0 .or. any(indices < 1) .or. &
         any(indices > size(values))) then
         allocate(sample(0))
         return
      end if
      sample = values(indices)
   end function resample_vector

   pure function resample_matrix(values, indices) result(sample)
      !! Select matrix rows using caller-supplied indices with replacement.
      real(dp), intent(in) :: values(:, :) !! Rows to resample.
      integer, intent(in) :: indices(:) !! Source row for every output row.
      real(dp), allocatable :: sample(:, :)

      if (size(values, 1) == 0 .or. any(indices < 1) .or. &
         any(indices > size(values, 1))) then
         allocate(sample(0, size(values, 2)))
         return
      end if
      sample = values(indices, :)
   end function resample_matrix

   pure function block_resample_vector(values, block_starts, block_size) &
      result(sample)
      !! Form a circular moving-block sample from supplied block starts.
      real(dp), intent(in) :: values(:) !! Values to resample.
      integer, intent(in) :: block_starts(:) !! One-based start of every sampled block.
      integer, intent(in) :: block_size !! Observations copied from each block.
      real(dp), allocatable :: sample(:)
      integer :: output_index, block, offset, source

      if (size(values) == 0 .or. block_size < 1 .or. &
         size(block_starts)*block_size < size(values) .or. &
         any(block_starts < 1) .or. any(block_starts > size(values))) then
         allocate(sample(0))
         return
      end if
      allocate(sample(size(values)))
      output_index = 0
      do block = 1, size(block_starts)
         do offset = 0, block_size - 1
            if (output_index == size(values)) return
            output_index = output_index + 1
            source = 1 + modulo(block_starts(block) - 1 + offset, size(values))
            sample(output_index) = values(source)
         end do
      end do
   end function block_resample_vector

   pure function block_resample_matrix(values, block_starts, block_size) &
      result(sample)
      !! Form a circular moving-block row sample from supplied block starts.
      real(dp), intent(in) :: values(:, :) !! Rows to resample.
      integer, intent(in) :: block_starts(:) !! One-based start of every sampled block.
      integer, intent(in) :: block_size !! Rows copied from each block.
      real(dp), allocatable :: sample(:, :)
      integer :: output_index, block, offset, source

      if (size(values, 1) == 0 .or. block_size < 1 .or. &
         size(block_starts)*block_size < size(values, 1) .or. &
         any(block_starts < 1) .or. any(block_starts > size(values, 1))) then
         allocate(sample(0, size(values, 2)))
         return
      end if
      allocate(sample(size(values, 1), size(values, 2)))
      output_index = 0
      do block = 1, size(block_starts)
         do offset = 0, block_size - 1
            if (output_index == size(values, 1)) return
            output_index = output_index + 1
            source = 1 + modulo(block_starts(block) - 1 + offset, &
               size(values, 1))
            sample(output_index, :) = values(source, :)
         end do
      end do
   end function block_resample_matrix

   pure function additive_resample_vector(values, perturbations) result(sample)
      !! Add supplied perturbations as in tsDyn's additive wild schemes.
      real(dp), intent(in) :: values(:) !! Values to perturb.
      real(dp), intent(in) :: perturbations(:) !! Additive perturbation for each value.
      real(dp), allocatable :: sample(:)

      if (size(values) /= size(perturbations)) then
         allocate(sample(0))
         return
      end if
      sample = values + perturbations
   end function additive_resample_vector

   pure function additive_resample_matrix(values, perturbations) result(sample)
      !! Add one supplied perturbation to every component of each row.
      real(dp), intent(in) :: values(:, :) !! Rows to perturb.
      real(dp), intent(in) :: perturbations(:) !! Additive perturbation for each row.
      real(dp), allocatable :: sample(:, :)

      if (size(values, 1) /= size(perturbations)) then
         allocate(sample(0, size(values, 2)))
         return
      end if
      sample = values + spread(perturbations, 2, size(values, 2))
   end function additive_resample_matrix

   pure function wild_resample_vector(values, multipliers) result(sample)
      !! Multiply values by caller-supplied wild-bootstrap multipliers.
      real(dp), intent(in) :: values(:) !! Values to perturb.
      real(dp), intent(in) :: multipliers(:) !! Multiplier for each value.
      real(dp), allocatable :: sample(:)

      if (size(values) /= size(multipliers)) then
         allocate(sample(0))
         return
      end if
      sample = values*multipliers
   end function wild_resample_vector

   pure function wild_resample_matrix(values, multipliers) result(sample)
      !! Multiply every row by its caller-supplied wild-bootstrap multiplier.
      real(dp), intent(in) :: values(:, :) !! Rows to perturb.
      real(dp), intent(in) :: multipliers(:) !! Multiplier for each row.
      real(dp), allocatable :: sample(:, :)

      if (size(values, 1) /= size(multipliers)) then
         allocate(sample(0, size(values, 2)))
         return
      end if
      sample = values*spread(multipliers, 2, size(values, 2))
   end function wild_resample_matrix

end module resampling_mod
