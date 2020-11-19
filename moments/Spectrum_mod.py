"""
Contains Spectrum object, which represents frequency spectra.
"""
import logging

logging.basicConfig()
logger = logging.getLogger("Spectrum_mod")
import functools
import operator
import os
import sys
import numpy
from numpy import newaxis as nuax
import scipy.misc as misc
import copy

# Account for difference in scipy installations.
try:
    from scipy.misc import comb
except ImportError:
    try:
        from scipy.special import comb
    except:
        from scipy import comb
from scipy.integrate import trapz
from scipy.special import betainc

import moments.Integration
import moments.Integration_nomig
from . import Numerics

plotting = True
try:
    import moments.ModelPlot
except ImportError:  # if matplotlib is not present, do not import, and do not run plotting
    # functions
    plotting = False


class Spectrum(numpy.ma.masked_array):
    """
    Represents a single-locus biallelic frequency spectrum.

    Spectra are represented by masked arrays. The masking allows us to ignore
    specific entries in the spectrum. When simulating under the standard infinite
    sites model (ISM), the entries we mask are the bins specifying absent or fixed
    variants. When using a reversible mutation model (i.e. the finite genome model),
    we track the density of variants in fixed bins, setting ``mask_corners`` to
    ``False``.

    To do: check that all optional parameters are functional with tests.

    :param data: An array with dimension equal to the number of populations.
        Each dimension has length :math:`n_i+1`, where :math:`n_i` is the
        sample size for the i-th population.
    :type data: array
    :param mask: An optional array of the same size as data. 'True' entries in
        this array are masked in the Spectrum. These represent missing
        data categories. (For example, you may not trust your singleton
        SNP calling.)
    :type maks: array, optional
    :param mask_corners: If True (default), the 'observed in none' and 'observed
        in all' entries of the FS will be masked. Typically these
        entries are masked. In the defaul infinite sites model, moments does
        not reliably calculate the fixed-bin entries, so you will almost always
        want ``mask_corners=True``. The exception is if we are simulating under
        the finite genome model, in which case we track the probability of
        a site to be fixed for either allele.
    :type maks_corners: bool, optional
    :param data_folded: If True, it is assumed that the input data is folded. An
        error will be raised if the input data and mask are not
        consistent with a folded Spectrum.
    :type data_folded: bool, optional
    :param check_folding: If True and data_folded=True, the data and mask will be
        checked to ensure they are consistent with a folded
        Spectrum. If they are not, a warning will be printed.
    :type check_folding: bool, optional
    :param pop_ids: Optional list of strings containing the population labels,
        with length equal to the dimension of ``data``.
    :type pop_ids: list of strings, optional

    :return: A frequency spectrum object, as a masked array.
    """

    def __new__(
        subtype,
        data,
        mask=numpy.ma.nomask,
        mask_corners=True,
        data_folded=None,
        check_folding=True,
        dtype=float,
        copy=True,
        fill_value=numpy.nan,
        keep_mask=True,
        shrink=True,
        pop_ids=None,
    ):
        data = numpy.asanyarray(data)

        if mask is numpy.ma.nomask:
            mask = numpy.ma.make_mask_none(data.shape)

        subarr = numpy.ma.masked_array(
            data,
            mask=mask,
            dtype=dtype,
            copy=copy,
            fill_value=numpy.nan,
            keep_mask=True,
            shrink=True,
        )
        subarr = subarr.view(subtype)

        if hasattr(data, "folded"):
            if data_folded is None or data_folded == data.folded:
                subarr.folded = data.folded
            elif data_folded != data.folded:
                raise ValueError(
                    "Data does not have same folding status as "
                    "was called for in Spectrum constructor."
                )
        elif data_folded is not None:
            subarr.folded = data_folded
        else:
            subarr.folded = False

        # Check that if we're declaring that the input data is folded, it
        # actually is, and the mask reflects this.
        if data_folded:
            total_samples = numpy.sum(subarr.sample_sizes)
            total_per_entry = subarr._total_per_entry()
            # Which entries are nonsense in the folded fs.
            where_folded_out = total_per_entry > int(total_samples / 2)
            if check_folding and not numpy.all(subarr.data[where_folded_out] == 0):
                logger.warn(
                    "Creating Spectrum with data_folded = True, but "
                    "data has non-zero values in entries which are "
                    "nonsensical for a folded Spectrum."
                )
            if check_folding and not numpy.all(subarr.mask[where_folded_out]):
                logger.warn(
                    "Creating Spectrum with data_folded = True, but "
                    "mask is not True for all entries which are "
                    "nonsensical for a folded Spectrum."
                )

        if hasattr(data, "pop_ids"):
            if pop_ids is None or pop_ids == data.pop_ids:
                subarr.pop_ids = data.pop_ids
            elif pop_ids != data.pop_ids:
                logger.warn(
                    "Changing population labels in construction of new " "Spectrum."
                )
                if len(pop_ids) != subarr.ndim:
                    raise ValueError(
                        "pop_ids must be of length equal to "
                        "dimensionality of Spectrum."
                    )
                subarr.pop_ids = pop_ids
        else:
            if pop_ids is not None and len(pop_ids) != subarr.ndim:
                raise ValueError(
                    "pop_ids must be of length equal to " "dimensionality of Spectrum."
                )
            subarr.pop_ids = pop_ids

        if mask_corners:
            subarr.mask_corners()

        return subarr

    # See http://www.scipy.org/Subclasses for information on the
    # __array_finalize__ and __array_wrap__ methods. I had to do some debugging
    # myself to discover that I also needed _update_from.
    # Also, see http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
    # Also, see http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    #
    # We need these methods to ensure extra attributes get copied along when
    # we do arithmetic on the FS.
    def __array_finalize__(self, obj):
        if obj is None:
            return
        numpy.ma.masked_array.__array_finalize__(self, obj)
        self.folded = getattr(obj, "folded", "unspecified")
        self.pop_ids = getattr(obj, "pop_ids", None)

    def __array_wrap__(self, obj, context=None):
        result = obj.view(type(self))
        result = numpy.ma.masked_array.__array_wrap__(self, obj, context=context)
        result.folded = self.folded
        result.pop_ids = self.pop_ids
        return result

    def _update_from(self, obj):
        numpy.ma.masked_array._update_from(self, obj)
        if hasattr(obj, "folded"):
            self.folded = obj.folded
        if hasattr(obj, "pop_ids"):
            self.pop_ids = obj.pop_ids

    # masked_array has priority 20.
    __array_priority__ = 20

    def __repr__(self):
        return "Spectrum(%s, folded=%s, pop_ids=%s)" % (
            str(self),
            str(self.folded),
            str(self.pop_ids),
        )

    # Functions for manipulating frequency spectra.
    def mask_corners(self):
        """
        Mask the 'seen in 0 samples' and 'seen in all samples' entries.
        """
        self.mask.flat[0] = self.mask.flat[-1] = True

    def unmask_all(self):
        """
        Unmask all entires of the frequency spectrum.
        """
        self.mask[tuple([slice(None)] * self.Npop)] = False

    def _get_sample_sizes(self):
        return numpy.asarray(self.shape) - 1

    sample_sizes = property(_get_sample_sizes)

    def _get_Npop(self):
        return self.ndim

    Npop = property(_get_Npop)

    def _ensure_dimension(self, Npop):
        """
        Ensure that fs has Npop dimensions.
        """
        if not self.Npop == Npop:
            raise ValueError("Only compatible with %id spectra." % Npop)

    def project(self, ns):
        """
        Project to smaller sample size.

        ``project`` does *not* act in-place, so that the input frequency
        spectrum is not changed.

        :param ns: Sample sizes for new spectrum.
        :type ns: list of integers
        """
        if len(ns) != self.Npop:
            raise ValueError(
                "Requested sample sizes not of same dimension "
                "as spectrum. Perhaps you need to marginalize "
                "over some populations first?"
            )
        if numpy.any(numpy.asarray(ns) > numpy.asarray(self.sample_sizes)):
            raise ValueError(
                "Cannot project to a sample size greater than "
                "original. Original size is %s and requested size "
                "is %s." % (self.sample_sizes, ns)
            )

        original_folded = self.folded
        # If we started with an folded Spectrum, we need to unfold before
        # projecting.
        if original_folded:
            output = self.unfold()
        else:
            output = self.copy()

        # Iterate over each axis, applying the projection.
        for axis, proj in enumerate(ns):
            if proj != self.sample_sizes[axis]:
                output = output._project_one_axis(proj, axis)

        output.pop_ids = self.pop_ids

        # Return folded or unfolded as original.
        if original_folded:
            return output.fold()
        else:
            return output

    def _project_one_axis(self, n, axis=0):
        """
        Project along a single axis.
        """
        # This gets a little tricky with fancy indexing to make it work
        # for fs with arbitrary number of dimensions.
        if n > self.sample_sizes[axis]:
            raise ValueError(
                "Cannot project to a sample size greater than "
                "original. Called sizes were from %s to %s."
                % (self.sample_sizes[axis], n)
            )

        newshape = list(self.shape)
        newshape[axis] = n + 1
        # Create a new empty fs that we'll fill in below.
        pfs = Spectrum(numpy.zeros(newshape), mask_corners=False)

        # Set up for our fancy indexes. These slices are currently like
        # [:,:,...]
        from_slice = [slice(None) for ii in range(self.Npop)]
        to_slice = [slice(None) for ii in range(self.Npop)]
        proj_slice = [nuax for ii in range(self.Npop)]

        proj_from = self.sample_sizes[axis]
        # For each possible number of hits.
        for hits in range(proj_from + 1):
            # Adjust the slice in the array we're projecting from.
            from_slice[axis] = slice(hits, hits + 1)
            # These are the least and most possible hits we could have in the
            #  projected fs.
            least, most = max(n - (proj_from - hits), 0), min(hits, n)
            to_slice[axis] = slice(least, most + 1)
            # The projection weights.
            proj = Numerics._cached_projection(n, proj_from, hits)
            proj_slice[axis] = slice(least, most + 1)
            # Do the multiplications
            pfs.data[tuple(to_slice)] += (
                self.data[tuple(from_slice)] * proj[tuple(proj_slice)]
            )
            pfs.mask[tuple(to_slice)] = numpy.logical_or(
                pfs.mask[tuple(to_slice)], self.mask[tuple(from_slice)]
            )

        return pfs

    def marginalize(self, over, mask_corners=None):
        """
        Reduced dimensionality spectrum summing over the set of populations
        given by ``over``.

        ``marginalize`` does not act in-place, so the input frequency spectrum
        will not be altered.

        :param over: List of axes to sum over. For example (0,2) will marginalize
            populations 0 and 2.
        :type over: list of integers
        :param mask_corners: If True, the fixed bins of the resulting spectrum will be
            masked. The default behavior is to mask the corners only if at least one
            of the corners of the input frequency spectrum is masked. If either
            corner is masked, the output frequency spectrum masks the fixed bins.
        :type mask_corners: bool, optional
        """
        if plotting:
            # Update ModelPlot
            model = moments.ModelPlot._get_model()
            if model is not None:
                model.extinction(over)

        original_folded = self.folded
        # If we started with an folded Spectrum, we need to unfold before
        # marginalizing.
        if original_folded:
            output = self.unfold()
        else:
            output = self.copy()

        orig_mask = output.mask.copy()
        orig_mask.flat[0] = orig_mask.flat[-1] = False
        if numpy.any(orig_mask):
            logger.warn(
                "Marginalizing a Spectrum with internal masked values. "
                "This may not be a well-defined operation."
            )

        # Do the marginalization
        for axis in sorted(over)[::-1]:
            output = output.sum(axis=axis)
        pop_ids = None
        if self.pop_ids is not None:
            pop_ids = list(self.pop_ids)
            for axis in sorted(over)[::-1]:
                del pop_ids[axis]
        output.folded = False
        output.pop_ids = pop_ids

        if mask_corners is None:
            if self.mask.flat[0] == True or self.mask.flat[-1] == True:
                mask_corners = True
            else:
                mask_corners = False
        if mask_corners:
            output.mask_corners()

        # Return folded or unfolded as original.
        if original_folded:
            return output.fold()
        else:
            return output

    def swap_axes(self, ax1, ax2):
        """
        Uses numpy's swapaxes function, but also swaps pop_ids as appropriate
        if pop_ids are given.

        Note that `fs.swapaxes(ax1, ax2)` will still work, but if population
        ids are given, it won't swap the pop_ids entries as expected.

        :param ax1: The index of the first population to swap.
        :type ax1: int
        :param ax2: The index of the second population to swap.
        :type ax2: int
        """
        output = numpy.swapaxes(self, ax1, ax2)
        if output.pop_ids is not None:
            pop1, pop2 = output.pop_ids[ax1], output.pop_ids[ax2]
            output.pop_ids[ax1], output.pop_ids[ax2] = pop2, pop1
        return output

    def _counts_per_entry(self):
        """
        Counts per population for each entry in the fs.
        """
        ind = numpy.indices(self.shape)
        # Transpose the first access to the last, so ind[ii,jj,kk] = [ii,jj,kk]
        ind = ind.transpose(list(range(1, self.Npop + 1)) + [0])
        return ind

    def _total_per_entry(self):
        """
        Total derived alleles for each entry in the fs.
        """
        return numpy.sum(self._counts_per_entry(), axis=-1)

    def log(self):
        """
        Returns the natural logarithm of the entries of the frequency spectrum.

        Only necessary because numpy.ma.log now fails to propagate extra
        attributes after numpy 1.10.
        """
        logfs = numpy.ma.log(self)
        logfs.folded = self.folded
        logfs.pop_ids = self.pop_ids
        return logfs

    def fold(self):
        """
        Returns a folded frequency spectrum.

        The folded fs assumes that information on which allele is ancestral or
        derived is unavailable. Thus the fs is in terms of minor allele
        frequency.  Note that this makes the fs into a "triangular" array.

        Note that if a masked cell is folded into non-masked cell, the
        destination cell is masked as well.

        Note also that folding is not done in-place. The return value is a new
        Spectrum object.
        """
        if self.folded:
            raise ValueError("Input Spectrum is already folded.")

        # How many samples total do we have? The folded fs can only contain
        # entries up to total_samples/2 (rounded down).
        total_samples = numpy.sum(self.sample_sizes)

        total_per_entry = self._total_per_entry()

        # Here's where we calculate which entries are nonsense in the folded fs.
        where_folded_out = total_per_entry > int(total_samples / 2)

        original_mask = self.mask
        # Here we create a mask that masks any values that were masked in
        # the original fs (or folded onto by a masked value).
        final_mask = numpy.logical_or(
            original_mask, Numerics.reverse_array(original_mask)
        )

        # To do the actual folding, we take those entries that would be folded
        # out, reverse the array along all axes, and add them back to the
        # original fs.
        reversed = Numerics.reverse_array(numpy.where(where_folded_out, self, 0))
        folded = numpy.ma.masked_array(self.data + reversed)
        folded.data[where_folded_out] = 0

        # Deal with those entries where assignment of the minor allele is
        # ambiguous.
        where_ambiguous = total_per_entry == total_samples / 2.0
        ambiguous = numpy.where(where_ambiguous, self, 0)
        folded += -0.5 * ambiguous + 0.5 * Numerics.reverse_array(ambiguous)

        # Mask out the remains of the folding operation.
        final_mask = numpy.logical_or(final_mask, where_folded_out)

        outfs = Spectrum(
            folded, mask=final_mask, data_folded=True, pop_ids=self.pop_ids
        )
        return outfs

    def unfold(self):
        """
        Returns an unfolded frequency spectrum.

        It is assumed that each state of a SNP is equally likely to be
        ancestral.

        Note also that unfolding is not done in-place. The return value is a new
        Spectrum object.
        """
        if not self.folded:
            raise ValueError("Input Spectrum is not folded.")

        # Unfolding the data is easy.
        reversed_data = Numerics.reverse_array(self.data)
        newdata = (self.data + reversed_data) / 2.0

        # Unfolding the mask is trickier. We want to preserve masking of entries
        # that were masked in the original Spectrum.
        # Which entries in the original Spectrum were masked solely because
        # they are incompatible with a folded Spectrum?
        total_samples = numpy.sum(self.sample_sizes)
        total_per_entry = self._total_per_entry()
        where_folded_out = total_per_entry > int(total_samples / 2)

        newmask = numpy.logical_xor(self.mask, where_folded_out)
        newmask = numpy.logical_or(newmask, Numerics.reverse_array(newmask))

        outfs = Spectrum(newdata, mask=newmask, data_folded=False, pop_ids=self.pop_ids)
        return outfs

    # Functions that apply demographic events, including integration.
    def split(self, idx, n0, n1, new_ids=None):
        """
        Splits a population in the SFS into two populations, with the extra
        population placed at the end. Returns a new frequency spectrum.

        :param idx: The index of the population to split.
        :type idx: int
        :param n0: The sample size of the first split population.
        :type n0: int
        :param n1: The sample size of the second split population.
        :type n1: int
        :param new_ids: The population IDs of the split populations. Can only be
            used if pop_ids are given for the input spectrum.
        :type new_ids: list of strings, optional
        """
        if self.folded:
            raise ValueError("Cannot perform split on folded spectrum.")
        if self.pop_ids is None and new_ids is not None:
            raise ValueError("Trying to assign ids to a SFS with no pop_ids.")
        if new_ids is not None and len(new_ids) != 2:
            raise ValueError("new_ids must be a list of two population id strings")
        fs_split = moments.Manips.split_by_index(self, idx, n0, n1)
        if new_ids is not None:
            fs_split.pop_ids = self.pop_ids
            fs_split.pop_ids[idx] = new_ids[0]
            fs_split.pop_ids.append(new_ids[1])
        return fs_split

    def admix(self, idx0, idx1, num_lineages, proportion, new_id=None):
        """
        Returns a new frequency spectrum with an admixed population that arose through
        admixture from indexed populations with given number of lineages and
        proportions from parental populations. This serves as a wrapper for
        ``Manips.admix_into_new``, with the added feature of handling pop_ids.

        Note that if the number of lineages that move are equal to the number
        of lineages previously present in a source population, that source
        population is marginalized.

        :param idx0: Index of first source population.
        :type idx0: int
        :param idx1: Index of second source population.
        :type idx1: int
        :param num_lineages: Number of lineages in the new population. Cannot be
            greater than the number of existing lineages in either source
            populations.
        :type num_lineages: int
        :param proportion: The proportion of lineages that come from the first
            source population (1-proportion acestry comes from the second source
            population). Must be a number between 0 and 1.
        :type proportion: float
        :param new_id: The ID of the new population. Can only be used if the
            population IDs are specified in the input SFS.
        :type new_id: str, optional
        """
        if new_id is not None and self.pop_ids is None:
            raise ValueError("Cannot specify new pop ids if input SFS has no pop_ids")
        if proportion < 0 or proportion > 1:
            raise ValueError("proportion must be between 0 and 1")
        if idx0 == idx1:
            raise ValueError("Cannot admix population with itself")
        if idx0 < 0 or idx0 >= self.Npop or idx1 < 0 or idx1 >= self.Npop:
            raise ValueError(f"Population indexes must be between 0 and {self.Npop-1}")
        fs_admix = moments.Manips.admix_into_new(
            self, idx0, idx1, num_lineages, proportion
        )
        if new_id is not None:
            new_pop_ids = copy.copy(self.pop_ids)
            # remove pop ids for marginalized pops
            for idx in sorted([idx0, idx1])[::-1]:
                if self.sample_sizes[idx] == num_lineages:
                    del new_pop_ids[idx]
            new_pop_ids.append(new_id)
            fs_admix.pop_ids = new_pop_ids
        return fs_admix

    def pulse_migrate(self, idx_from, idx_to, keep_from, proportion):
        """
        Mass migration (pulse admixture) between two existing populations. The
        target (destination) population has the same number of lineages in the
        output SFS, and the source population has ``keep_from`` number of lineages
        after the pulse event. The proportion is the expected ancestry proportion
        in the target population that comes from the source population.

        This serves as a wrapper for ``Manips.admix_inplace``.

        Note that depending on the proportion and number of lineages, because this
        is an approximate operation, we often need a large number of lineages from
        the source population to maintain accuracy.

        :param idx_from: Index of source population.
        :type idx_from: int
        :param idx_to: Index of targeet population.
        :type idx_to: int
        :param keep_from: Number of lineages to keep in source population.
        :type keep_from: int
        :param proportion: Ancestry proportion of source population that migrates
            to target population.
        :type proportion: float
        """
        if idx_from < 0 or idx_from >= self.Npop or idx_to < 0 or idx_to >= self.Npop:
            raise ValueError(f"Invalid population index for {self.Npop}D SFS.")
        if proportion < 0 or proportion > 1:
            raise ValueError("proportion must be between 0 and 1")
        if idx_from == idx_to:
            raise ValueError("Cannot admix population into itself")
        fs_pulse = moments.Manips.admix_inplace(
            self, idx_from, idx_to, keep_from, proportion
        )
        if self.pop_ids is not None:
            fs_pulse.pop_ids = self.pop_ids
        return fs_pulse

    # Integrate the SFS in-place
    def integrate(
        self,
        Npop,
        tf,
        dt_fac=0.02,
        gamma=None,
        h=None,
        m=None,
        theta=1.0,
        adapt_dt=False,
        finite_genome=False,
        theta_fd=None,
        theta_bd=None,
        frozen=[False],
    ):
        """
        Method to simulate the spectrum's evolution for a given set of demographic
        parameters. The SFS is integrated forward-in-time, and the integration
        occurs in-place, meaning you need only call ``fs.integrate( )``, and the
        ``fs`` is updated.

        :param Npop: List of populations' relative effective sizes. Can be given
            as a list of positive values for constant sizes, or as a function that
            returns a list of sizes at a given time.
        :type Npop: list or function that returns a list
        :param tf: The total integration time in genetic units.
        :type tf: float
        :param dt_fac: The timestep factor, default is 0.02
        :type dt_fac: float, optional
        :param gamma: The selection coefficient (:math:`2 N_e s`), or list of selection
            coefficients if more than one population.
        :type gamma: float or list of floats, optional
        :param h: The dominance coefficient, or list of dominance coefficients within
            each population, if more than one population.
        :type h: float or list of floats, optional
        :param m: The migration rates matrix as an N-D array, where m[i,j] is the
            migration rate from pop j to pop i, normalized by :math:`2N_e`.
        :type m: array-like, optional
        :param theta: The scaled mutation rate :math:`4 N_e u`, which defaults to 1.
            ``theta`` can be used in the reversible model in the case of symmetric
            mutation rates. In this case, ``theta`` must be set to << 1.
        :type theta: float, optional
        :param adapt_dt: flag to allow dt correction avoiding negative entries.
        :type adapt_dt: bool, optional
        :param finite_genome: If True, simulate under the finite-genome model with
            reversible mutations. If using this model, we specify the forward
            and backward mutation rates, which are per-base rates that are not
            scaled by number of mutable loci (different from the standard ISM
            model). Defaults to False.
        :type finite_genome: bool, optional
        :param theta_fd: The forward mutation rate :math:`4 Ne u`.
        :type theta_fd: float, optional
        :param theta_bd: The backward mutation rate :math:`4 Ne v`.
        :type theta_bd: float, optional
        :param frozen: list of same length as number of pops, with True for frozen
            populations at the corresponding index.
        :type frozen: list of bools
        """
        n = numpy.array(self.shape) - 1

        if m is not None:
            m = numpy.array(m)

        if finite_genome == True and (theta_fd == None or theta_bd == None):
            if theta >= 1:
                raise ValueError(
                    "In the finite genome model, theta must be much less than 1. "
                    "If symmetric mutation rates, can use theta << 1. Otherwise, "
                    "theta_fd and theta_bd must be specified."
                )
            else:
                theta_fd = theta_bd = theta

        if hasattr(Npop, "__len__"):
            if numpy.any(frozen) and len(Npop) != len(frozen):
                raise ValueError(
                    "If one or more populations are frozen, length "
                    "of frozen must match number of simulated pops."
                )
        else:
            if numpy.any(frozen) and len(Npop(0)) != len(frozen):
                raise ValueError(
                    "If one or more populations are frozen, length "
                    "of frozen must match number of simulated pops."
                )

        if plotting:
            model = moments.ModelPlot._get_model()
            if model is not None:
                model.evolve(tf, Npop, m)

        if len(n) == 1:
            if gamma is None:
                gamma = 0.0
            if h is None:
                h = 0.5
            if gamma == 0:
                self.data[:] = moments.Integration_nomig.integrate_neutral(
                    self.data,
                    Npop,
                    tf,
                    dt_fac,
                    theta,
                    finite_genome=finite_genome,
                    theta_fd=theta_fd,
                    theta_bd=theta_bd,
                    frozen=frozen,
                )
            else:
                # self.data[:] = integrate_1D(self.data, Npop, n, tf, dt_fac, dt_max, gamma, h, theta)
                self.data[:] = moments.Integration_nomig.integrate_nomig(
                    self.data,
                    Npop,
                    tf,
                    dt_fac,
                    gamma,
                    h,
                    theta,
                    finite_genome=finite_genome,
                    theta_fd=theta_fd,
                    theta_bd=theta_bd,
                    frozen=frozen,
                )
        else:
            if gamma is None:
                gamma = numpy.zeros(len(n))
            if h is None:
                h = 0.5 * numpy.ones(len(n))
            if m is None:
                m = numpy.zeros([len(n), len(n)])
            if (m == 0).all():
                # for more than 2 populations, the sparse solver seems to be faster than the tridiag...
                if (numpy.array(gamma) == 0).all() and len(n) < 3:
                    self.data[:] = moments.Integration_nomig.integrate_neutral(
                        self.data,
                        Npop,
                        tf,
                        dt_fac,
                        theta,
                        finite_genome=finite_genome,
                        theta_fd=theta_fd,
                        theta_bd=theta_bd,
                        frozen=frozen,
                    )
                else:
                    self.data[:] = moments.Integration_nomig.integrate_nomig(
                        self.data,
                        Npop,
                        tf,
                        dt_fac,
                        gamma,
                        h,
                        theta,
                        finite_genome=finite_genome,
                        theta_fd=theta_fd,
                        theta_bd=theta_bd,
                        frozen=frozen,
                    )
            else:
                self.data[:] = moments.Integration.integrate_nD(
                    self.data,
                    Npop,
                    tf,
                    dt_fac,
                    gamma,
                    h,
                    m,
                    theta,
                    adapt_dt,
                    finite_genome=finite_genome,
                    theta_fd=theta_fd,
                    theta_bd=theta_bd,
                    frozen=frozen,
                )

    # Functions for computing statistics from frequency spetra.
    def Fst(self):
        """
        Wright's Fst between the populations represented in the fs.

        This estimate of Fst assumes random mating, because we don't have
        heterozygote frequencies in the fs.

        Calculation is by the method of Weir and Cockerham _Evolution_ 38:1358
        (1984).  For a single SNP, the relevant formula is at the top of page
        1363. To combine results between SNPs, we use the weighted average
        indicated by equation 10.
        """
        # This gets a little obscure because we want to be able to work with
        # spectra of arbitrary dimension.

        # First quantities from page 1360
        r = self.Npop
        ns = self.sample_sizes
        nbar = numpy.mean(ns)
        nsum = numpy.sum(ns)
        nc = (nsum - numpy.sum(ns ** 2) / nsum) / (r - 1)

        # counts_per_pop is an r+1 dimensional array, where the last axis simply
        # records the indices of the entry.
        # For example, counts_per_pop[4,19,8] = [4,19,8]
        counts_per_pop = numpy.indices(self.shape)
        counts_per_pop = numpy.transpose(
            counts_per_pop, axes=list(range(1, r + 1)) + [0]
        )

        # The last axis of ptwiddle is now the relative frequency of SNPs in
        # that bin in each of the populations.
        ptwiddle = 1.0 * counts_per_pop / ns

        # Note that pbar is of the same shape as fs...
        pbar = numpy.sum(ns * ptwiddle, axis=-1) / nsum

        # We need to use 'this_slice' to get the proper aligment between
        # ptwiddle and pbar.
        this_slice = [slice(None)] * r + [numpy.newaxis]
        s2 = numpy.sum(ns * (ptwiddle - pbar[tuple(this_slice)]) ** 2, axis=-1) / (
            (r - 1) * nbar
        )

        # Note that this 'a' differs from equation 2, because we've used
        # equation 3 and b = 0 to solve for hbar.
        a = (
            nbar
            / nc
            * (s2 - 1 / (2 * nbar - 1) * (pbar * (1 - pbar) - (r - 1) / r * s2))
        )
        d = 2 * nbar / (2 * nbar - 1) * (pbar * (1 - pbar) - (r - 1) / r * s2)

        # The weighted sum over loci.
        asum = (self * a).sum()
        dsum = (self * d).sum()

        return asum / (asum + dsum)

    def S(self):
        """
        Returns the number of segregating sites in the frequency spectrum.
        """
        oldmask = self.mask.copy()
        self.mask_corners()
        S = self.sum()
        self.mask = oldmask
        return S

    def Watterson_theta(self):
        """
        Returns Watterson's estimator of theta.

        Note that is only sensible for 1-dimensional spectra.
        """
        if self.Npop != 1:
            raise ValueError("Only defined on a one-dimensional fs.")

        n = self.sample_sizes[0]
        S = self.S()
        an = numpy.sum(1.0 / numpy.arange(1, n))

        return S / an

    def theta_L(self):
        """
        Returns theta_L as defined by Zeng et al. "Statistical Tests for Detecting
        Positive Selection by Utilizing High-Frequency Variants" (2006)
        Genetics

        Note that is only sensible for 1-dimensional spectra.
        """
        if self.Npop != 1:
            raise ValueError("Only defined on a one-dimensional fs.")

        n = self.sample_sizes[0]
        return numpy.sum(numpy.arange(1, n) * self[1:n]) / (n - 1)

    def Zengs_E(self):
        """
        Returns Zeng et al.'s E statistic.

        From Zeng et al. "Statistical Tests for Detecting Positive Selection by
        Utilizing High-Frequency Variants" (2006) Genetics
        """
        num = self.theta_L() - self.Watterson_theta()

        n = self.sample_sizes[0]

        # See after Eq. 3
        an = numpy.sum(1.0 / numpy.arange(1, n))
        # See after Eq. 9
        bn = numpy.sum(1.0 / numpy.arange(1, n) ** 2)
        s = self.S()

        # See immediately after Eq. 12
        theta = self.Watterson_theta()
        theta_sq = s * (s - 1.0) / (an ** 2 + bn)

        # Eq. 14
        var = (n / (2.0 * (n - 1.0)) - 1.0 / an) * theta + (
            bn / an ** 2
            + 2.0 * (n / (n - 1.0)) ** 2 * bn
            - 2 * (n * bn - n + 1.0) / ((n - 1.0) * an)
            - (3.0 * n + 1.0) / (n - 1.0)
        ) * theta_sq

        return num / numpy.sqrt(var)

    def pi(self):
        """
        Returns the estimated expected number of pairwise differences between two
        chromosomes in the population.

        Note that this estimate includes a factor of
        sample_size / (sample_size-1) to make :math:`\\mathbb{E}[\\pi] = \\theta`.
        """
        if self.ndim != 1:
            raise ValueError("Only defined for a one-dimensional SFS.")

        n = self.sample_sizes[0]
        # sample frequencies p
        p = numpy.arange(0, n + 1, dtype=float) / n
        # This expression derives from Gillespie's _Population_Genetics:_A
        # _Concise_Guide_, 2nd edition, section 2.6.
        return n / (n - 1.0) * 2 * numpy.ma.sum(self * p * (1 - p))

    def Tajima_D(self):
        """
        Returns Tajima's D.

        Following Gillespie "Population Genetics: A Concise Guide" pg. 45

        """
        if not self.Npop == 1:
            raise ValueError("Only defined on a one-dimensional SFS.")

        S = self.S()

        n = 1.0 * self.sample_sizes[0]
        pihat = self.pi()
        theta = self.Watterson_theta()

        a1 = numpy.sum(1.0 / numpy.arange(1, n))
        a2 = numpy.sum(1.0 / numpy.arange(1, n) ** 2)
        b1 = (n + 1) / (3 * (n - 1))
        b2 = 2 * (n ** 2 + n + 3) / (9 * n * (n - 1))
        c1 = b1 - 1.0 / a1
        c2 = b2 - (n + 2) / (a1 * n) + a2 / a1 ** 2

        C = numpy.sqrt((c1 / a1) * S + c2 / (a1 ** 2 + a2) * S * (S - 1))

        return (pihat - theta) / C

    # Functions for saving and loading frequency spectra.
    # Make from_file a static method, so we can use it without an instance.
    @staticmethod
    def from_file(fid, mask_corners=True, return_comments=False):
        """
        Read frequency spectrum from file.

        See ``to_file`` for details on the file format.

        :param fid: string with file name to read from or an open file object.
        :type fid: string
        :param mask_corners: If True, mask the 'absent in all samples' and 'fixed in
            all samples' entries.
        :type mask_corners: bool, optional
        :param return_comments: If true, the return value is (fs, comments), where
            comments is a list of strings containing the comments
            from the file (without #'s).
        :type return_comments: bool, optional
        """
        newfile = False
        # Try to read from fid. If we can't, assume it's something that we can
        # use to open a file.
        if not hasattr(fid, "read"):
            newfile = True
            fid = open(fid, "r")

        line = fid.readline()
        # Strip out the comments
        comments = []
        while line.startswith("#"):
            comments.append(line[1:].strip())
            line = fid.readline()

        # Read the shape of the data
        shape_spl = line.split()
        if "folded" not in shape_spl and "unfolded" not in shape_spl:
            # This case handles the old file format
            shape = tuple([int(d) for d in shape_spl])
            folded = False
            pop_ids = None
        else:
            # This case handles the new file format
            shape, next_ii = [int(shape_spl[0])], 1
            while shape_spl[next_ii] not in ["folded", "unfolded"]:
                shape.append(int(shape_spl[next_ii]))
                next_ii += 1
            folded = shape_spl[next_ii] == "folded"
            # Are there population labels in the file?
            if len(shape_spl) > next_ii + 1:
                pop_ids = line.split('"')[1::2]
            else:
                pop_ids = None

        data = numpy.fromstring(
            fid.readline().strip(), count=numpy.product(shape), sep=" "
        )
        # fromfile returns a 1-d array. Reshape it to the proper form.
        data = data.reshape(*shape)

        maskline = fid.readline().strip()
        if not maskline:
            # The old file format didn't have a line for the mask
            mask = None
        else:
            # This case handles the new file format
            mask = numpy.fromstring(maskline, count=numpy.product(shape), sep=" ")
            mask = mask.reshape(*shape)

        # If we opened a new file, clean it up.
        if newfile:
            fid.close()

        fs = Spectrum(data, mask, mask_corners, data_folded=folded, pop_ids=pop_ids)

        if not return_comments:
            return fs
        else:
            return fs, comments

    fromfile = from_file

    def to_file(self, fid, precision=16, comment_lines=[], foldmaskinfo=True):
        """
        Write frequency spectrum to file.

        The file format is:

        - Any number of comment lines beginning with a '#'
        - A single line containing N integers giving the dimensions of the fs
          array. So this line would be '5 5 3' for an SFS that was 5x5x3.
          (That would be 4x4x2 *samples*.)
        - On the *same line*, the string 'folded' or 'unfolded' denoting the
          folding status of the array
        - On the *same line*, optional strings each containing the population
          labels in quotes separated by spaces, e.g. "pop 1" "pop 2"
        - A single line giving the array elements. The order of elements is
          e.g.: fs[0,0,0] fs[0,0,1] fs[0,0,2] ... fs[0,1,0] fs[0,1,1] ...
        - A single line giving the elements of the mask in the same order as
          the data line. '1' indicates masked, '0' indicates unmasked.

        :param fid: string with file name to write to or an open file object.
        :type fid: string
        :param precision: precision with which to write out entries of the SFS. (They
            are formated via %.<p>g, where <p> is the precision.) Defaults to 16.
        :type precision: int, optional
        :param comment_lines: list of strings to be used as comment lines in the header
            of the output file.
        :type comment_lines: list of strings, optional
        :param foldmaskinfo: If False, folding and mask and population label
            information will not be saved.
        :type foldmaskinfo: bool, optional
        """
        # Open the file object.
        newfile = False
        if not hasattr(fid, "write"):
            newfile = True
            fid = open(fid, "w")

        # Write comments
        for line in comment_lines:
            fid.write("# ")
            fid.write(line.strip())
            fid.write(os.linesep)

        # Write out the shape of the fs
        for elem in self.data.shape:
            fid.write("%i " % elem)

        if foldmaskinfo:
            if not self.folded:
                fid.write("unfolded")
            else:
                fid.write("folded")
            if self.pop_ids is not None:
                for label in self.pop_ids:
                    fid.write(' "%s"' % label)

        fid.write(os.linesep)

        # Write the data to the file
        self.data.tofile(fid, " ", "%%.%ig" % precision)
        fid.write(os.linesep)

        if foldmaskinfo:
            # Write the mask to the file
            numpy.asarray(self.mask, int).tofile(fid, " ")
            fid.write(os.linesep)

        # Close file
        if newfile:
            fid.close()

    ## Overide the (perhaps confusing) original numpy tofile method.
    tofile = to_file

    def fixed_size_sample(self, nsamples, include_masked=False):
        """
        Generate a resampled SFS from the current one. Thus, the resampled SFS
        follows a multinomial distribution given by the proportion of sites
        in each bin in the original SFS.

        :param nsamples: Number of samples to include in the new SFS.
        :type nsamples: int
        :param include_masked: If True, use all bins from the SFS. Otherwise,
            use only non-masked bins. Defaults to False.
        :type include_masked: bool, optional
        """
        flat = self.flatten()
        pvals = flat.data
        if include_masked is False:
            pvals[flat.mask] = 0
        pvals /= pvals.sum()

        sample = numpy.random.multinomial(int(nsamples), pvals)
        sample = sample.reshape(self.shape)

        return moments.Spectrum(sample, mask=self.mask, pop_ids=self.pop_ids)

    def sample(self):
        """
        Generate a Poisson-sampled fs from the current one.

        Note: Entries where the current fs is masked or 0 will be masked in the
        output sampled fs.
        """
        import scipy.stats

        # These are entries where the sampling has no meaning. Either the fs is
        # 0 there or masked.
        bad_entries = numpy.logical_or(self == 0, self.mask)
        # We convert to a 1-d array for passing into the sampler
        means = self.ravel().copy()
        # Filter out those bad entries.
        means[bad_entries.ravel()] = 1
        # Sample
        samp = scipy.stats.distributions.poisson.rvs(means, size=len(means))
        # Replace bad entries...
        samp[bad_entries.ravel()] = 0
        # Convert back to a properly shaped array
        samp = samp.reshape(self.shape)
        # Convert to a fs and mask the bad entries
        samp = Spectrum(
            samp, mask=bad_entries, data_folded=self.folded, pop_ids=self.pop_ids
        )
        return samp

    @staticmethod
    def from_ms_file(
        fid,
        average=True,
        mask_corners=True,
        return_header=False,
        pop_assignments=None,
        pop_ids=None,
        bootstrap_segments=1,
    ):
        """
        Read frequency spectrum from file of ms output.

        :param fid: string with file name to read from or an open file object.
        :param average: If True, the returned fs is the average over the runs in the ms
            file. If False, the returned fs is the sum.
        :param mask_corners: If True, mask the 'absent in all samples' and 'fixed in
            all samples' entries.
        :param return_header: If True, the return value is (fs, (command,seeds), where
            command and seeds are strings containing the ms
            commandline and the seeds used.
        :param pop_assignments: If None, the assignments of samples to populations is
            done automatically, using the assignment in the ms
            command line. To manually assign populations, pass a
            list of the from [6,8]. This example places
            the first 6 samples into population 1, and the next 8
            into population 2.
        :param pop_ids: Optional list of strings containing the population labels.
            If pop_ids is None, labels will be "pop0", "pop1", ...
        :param bootstrap_segments: If bootstrap_segments is an integer greater than 1,
            the data will be broken up into that many segments
            based on SNP position. Instead of single FS, a list
            of spectra will be returned, one for each segment.
        """
        newfile = False
        # Try to read from fid. If we can't, assume it's something that we can
        # use to open a file.
        if not hasattr(fid, "read"):
            newfile = True
            fid = open(fid, "r")

        # Parse the commandline
        command = line = fid.readline()
        command_terms = line.split()

        if command_terms[0].count("ms"):
            runs = int(command_terms[2])
            try:
                pop_flag = command_terms.index("-I")
                num_pops = int(command_terms[pop_flag + 1])
                pop_samples = [
                    int(command_terms[pop_flag + ii]) for ii in range(2, 2 + num_pops)
                ]
            except ValueError:
                num_pops = 1
                pop_samples = [int(command_terms[1])]
        else:
            raise ValueError("Unrecognized command string: %s." % command)

        total_samples = numpy.sum(pop_samples)
        if pop_assignments:
            num_pops = len(pop_assignments)
            pop_samples = pop_assignments

        sample_indices = numpy.cumsum([0] + pop_samples)
        bottom_l = sample_indices[:-1]
        top_l = sample_indices[1:]

        seeds = line = fid.readline()
        while not line.startswith("//"):
            line = fid.readline()

        counts = numpy.zeros(len(pop_samples), numpy.int_)
        fs_shape = numpy.asarray(pop_samples) + 1
        dimension = len(counts)

        if dimension > 1:
            bottom0 = bottom_l[0]
            top0 = top_l[0]
            bottom1 = bottom_l[1]
            top1 = top_l[1]
        if dimension > 2:
            bottom2 = bottom_l[2]
            top2 = top_l[2]
        if dimension > 3:
            bottom3 = bottom_l[3]
            top3 = top_l[3]
        if dimension > 4:
            bottom4 = bottom_l[4]
            top4 = top_l[4]
        if dimension > 5:
            bottom5 = bottom_l[5]
            top5 = top_l[5]

        all_data = [
            numpy.zeros(fs_shape, numpy.int_) for boot_ii in range(bootstrap_segments)
        ]
        for run_ii in range(runs):
            line = fid.readline()
            segsites = int(line.split()[-1])

            if segsites == 0:
                # Special case, need to read 3 lines to stay synced.
                for _ in range(3):
                    line = fid.readline()
                continue
            line = fid.readline()
            while not line.startswith("positions"):
                line = fid.readline()

            # Read SNP positions for creating bootstrap segments
            positions = [float(_) for _ in line.split()[1:]]
            # Where we should break our interval to create our bootstraps
            breakpts = numpy.linspace(0, 1, bootstrap_segments + 1)
            # The indices that correspond to those breakpoints
            break_iis = numpy.searchsorted(positions, breakpts)
            # Correct for searchsorted behavior if last position is 1,
            # to ensure all SNPs are captured
            break_iis[-1] = len(positions)

            # Read the chromosomes in
            chromos = fid.read((segsites + 1) * total_samples)

            # For each bootstrap segment, relevant SNPs run from start_ii:end_ii
            for boot_ii, (start_ii, end_ii) in enumerate(
                zip(break_iis[:-1], break_iis[1:])
            ):
                # Use the data array corresponding to this bootstrap segment
                data = all_data[boot_ii]
                for snp in range(start_ii, end_ii):
                    # Slice to get all the entries that refer to a given SNP
                    this_snp = chromos[snp :: segsites + 1]
                    # Count SNPs per population, and record them.
                    if dimension == 1:
                        data[this_snumpy.count("1")] += 1
                    elif dimension == 2:
                        data[
                            this_snp[bottom0:top0].count("1"),
                            this_snp[bottom1:top1].count("1"),
                        ] += 1
                    elif dimension == 3:
                        data[
                            this_snp[bottom0:top0].count("1"),
                            this_snp[bottom1:top1].count("1"),
                            this_snp[bottom2:top2].count("1"),
                        ] += 1
                    elif dimension == 4:
                        data[
                            this_snp[bottom0:top0].count("1"),
                            this_snp[bottom1:top1].count("1"),
                            this_snp[bottom2:top2].count("1"),
                            this_snp[bottom3:top3].count("1"),
                        ] += 1
                    elif dimension == 5:
                        data[
                            this_snp[bottom0:top0].count("1"),
                            this_snp[bottom1:top1].count("1"),
                            this_snp[bottom2:top2].count("1"),
                            this_snp[bottom3:top3].count("1"),
                            this_snp[bottom4:top4].count("1"),
                        ] += 1
                    elif dimension == 6:
                        data[
                            this_snp[bottom0:top0].count("1"),
                            this_snp[bottom1:top1].count("1"),
                            this_snp[bottom2:top2].count("1"),
                            this_snp[bottom3:top3].count("1"),
                            this_snp[bottom4:top4].count("1"),
                            this_snp[bottom5:top5].count("1"),
                        ] += 1
                    else:
                        # This is noticably slower, so we special case the cases
                        # above.
                        for dim_ii in range(dimension):
                            bottom = bottom_l[dim_ii]
                            top = top_l[dim_ii]
                            counts[dim_ii] = this_snp[bottom:top].count("1")
                        data[tuple(counts)] += 1

            # Read to the next iteration
            line = fid.readline()
            line = fid.readline()

        if newfile:
            fid.close()

        all_fs = [
            Spectrum(data, mask_corners=mask_corners, pop_ids=pop_ids)
            for data in all_data
        ]
        if average:
            all_fs = [fs / runs for fs in all_fs]

        # If we aren't setting up for bootstrapping, return fs, rather than a
        # list of length 1. (This ensures backward compatibility.)
        if bootstrap_segments == 1:
            all_fs = all_fs[0]

        if not return_header:
            return all_fs
        else:
            return all_fs, (command, seeds)

    @staticmethod
    def from_sfscode_file(
        fid,
        sites="all",
        average=True,
        mask_corners=True,
        return_header=False,
        pop_ids=None,
    ):
        """
        Read frequency spectrum from file of sfs_code output.

        :param fid: string with file name to read from or an open file object.
        :param sites: If sites=='all', return the fs of all sites. If sites == 'syn',
            use only synonymous mutations. If sites == 'nonsyn', use
            only non-synonymous mutations.
        :param average: If True, the returned fs is the average over the runs in the
            file. If False, the returned fs is the sum.
        :param mask_corners: If True, mask the 'absent in all samples' and 'fixed in
            all samples' entries.
        :param return_header: If true, the return value is (fs, (command,seeds), where
            command and seeds are strings containing the ms
            commandline and the seeds used.
        :param pop_ids: Optional list of strings containing the population labels.
            If pop_ids is None, labels will be "pop0", "pop1", ...
        """
        newfile = False
        # Try to read from fid. If we can't, assume it's something that we can
        # use to open a file.
        if not hasattr(fid, "read"):
            newfile = True
            fid = open(fid, "r")

        if sites == "all":
            only_nonsyn, only_syn = False, False
        elif sites == "syn":
            only_nonsyn, only_syn = False, True
        elif sites == "nonsyn":
            only_nonsyn, only_syn = True, False
        else:
            raise ValueError(
                "'sites' argument must be one of ('all', 'syn', " "'nonsyn')."
            )

        command = fid.readline()
        command_terms = command.split()

        runs = int(command_terms[2])
        num_pops = int(command_terms[1])

        # sfs_code default is 6 individuals, and I assume diploid pop
        pop_samples = [12] * num_pops
        if "--sampSize" in command_terms or "-n" in command_terms:
            try:
                pop_flag = command_terms.index("--sampSize")
                pop_flag = command_terms.index("-n")
            except ValueError:
                pass
            pop_samples = [
                2 * int(command_terms[pop_flag + ii]) for ii in range(1, 1 + num_pops)
            ]

        pop_samples = numpy.asarray(pop_samples)
        pop_digits = [str(i) for i in range(num_pops)]
        pop_fixed_str = [",%s.-1" % i for i in range(num_pops)]
        pop_count_str = [",%s." % i for i in range(num_pops)]

        seeds = fid.readline()
        line = fid.readline()

        data = numpy.zeros(numpy.asarray(pop_samples) + 1, numpy.int_)

        # line = //iteration...
        line = fid.readline()
        for iter_ii in range(runs):
            for ii in range(5):
                line = fid.readline()

            # It is possible for a mutation to be listed several times in the
            # output.  To accomodate this, I keep a dictionary of identities
            # for those mutations, and hold off processing them until I've seen
            # all mutations listed for the iteration.
            mut_dict = {}

            # Loop until this iteration ends.
            while not line.startswith("//") and line != "":
                split_line = line.split(";")
                if split_line[-1] == "\n":
                    split_line = split_line[:-1]

                # Loop over mutations on this line.
                for mut_ii, mutation in enumerate(split_line):
                    counts_this_mut = numpy.zeros(num_pops, numpy.int_)

                    split_mut = mutation.split(",")

                    # Exclude synonymous mutations
                    if only_nonsyn and split_mut[7] == "0":
                        continue
                    # Exclude nonsynonymous mutations
                    if only_syn and split_mut[7] == "1":
                        continue

                    ind_start = len(",".join(split_mut[:12]))
                    by_individual = mutation[ind_start:]

                    mut_id = ",".join(split_mut[:4] + split_mut[5:11])

                    # Count mutations in each population
                    for pop_ii, fixed_str, count_str in zip(
                        range(num_pops), pop_fixed_str, pop_count_str
                    ):
                        if fixed_str in by_individual:
                            counts_this_mut[pop_ii] = pop_samples[pop_ii]
                        else:
                            counts_this_mut[pop_ii] = by_individual.count(count_str)

                    # Initialize the list that will track the counts for this
                    # mutation. Using setdefault means that it won't overwrite
                    # if there's already a list stored there.
                    mut_dict.setdefault(mut_id, [0] * num_pops)
                    for ii in range(num_pops):
                        if counts_this_mut[ii] > 0 and mut_dict[mut_id][ii] > 0:
                            sys.stderr.write(
                                "Contradicting counts between "
                                "listings for mutation %s in "
                                "population %i." % (mut_id, ii)
                            )
                        mut_dict[mut_id][ii] = max(
                            counts_this_mut[ii], mut_dict[mut_id][ii]
                        )

                line = fid.readline()

            # Now apply all the mutations with fixations that we deffered.
            for mut_id, counts in mut_dict.items():
                if numpy.any(numpy.asarray(counts) > pop_samples):
                    sys.stderr.write(
                        "counts_this_mut > pop_samples: %s > "
                        "%s\n%s\n" % (counts, pop_samples, mut_id)
                    )
                    counts = numpy.minimum(counts, pop_samples)
                data[tuple(counts)] += 1

        if newfile:
            fid.close()

        fs = Spectrum(data, mask_corners=mask_corners, pop_ids=pop_ids)
        if average:
            fs /= runs

        if not return_header:
            return fs
        else:
            return fs, (command, seeds)

    def scramble_pop_ids(self, mask_corners=True):
        """
        Spectrum corresponding to scrambling individuals among populations.

        This is useful for assessing how diverged populations are.
        Essentially, it pools all the individuals represented in the fs and
        generates new populations of random individuals (without replacement)
        from that pool. If this fs is significantly different from the
        original, that implies population structure.
        """
        original_folded = self.folded
        # If we started with an folded Spectrum, we need to unfold before
        # projecting.
        if original_folded:
            self = self.unfold()

        total_samp = numpy.sum(self.sample_sizes)

        # First generate a 1d sfs for the pooled population.
        combined = numpy.zeros(total_samp + 1)
        # For each entry in the fs, this is the total number of derived alleles
        total_per_entry = self._total_per_entry()
        # Sum up to generate the equivalent 1-d spectrum.
        for derived, counts in zip(total_per_entry.ravel(), self.ravel()):
            combined[derived] += counts

        # Now resample back into a n-d spectrum
        # For each entry, this is the counts per popuation.
        #  e.g. counts_per_entry[3,4,5] = [3,4,5]
        counts_per_entry = self._counts_per_entry()
        # Reshape it to be 1-d, so we can iterate over it easily.
        counts_per_entry = counts_per_entry.reshape(numpy.prod(self.shape), self.ndim)
        resamp = numpy.zeros(self.shape)
        for counts, derived in zip(counts_per_entry, total_per_entry.ravel()):
            # The probability here is
            # (t1 choose d1)*(t2 choose d2)/(ntot choose derived)
            lnprob = sum(
                Numerics._lncomb(t, d) for t, d in zip(self.sample_sizes, counts)
            )
            lnprob -= Numerics._lncomb(total_samp, derived)
            prob = numpy.exp(lnprob)
            # Assign result using the appropriate weighting
            resamp[tuple(counts)] += prob * combined[derived]

        resamp = Spectrum(resamp, mask_corners=mask_corners)
        if not original_folded:
            return resamp
        else:
            return resamp.fold()

    @staticmethod
    def from_data_dict(
        data_dict, pop_ids, projections, mask_corners=True, polarized=True
    ):
        """
        Spectrum from a dictionary of polymorphisms.

        The data dictionary should be organized as:
        
        .. code-block::

            {snp_id: {
                'segregating': ['A','T'],
                'calls': {
                    'YRI': (23,3),
                    'CEU': (7,3)
                },
                'outgroup_allele': 'T'
            }}

        The 'calls' entry gives the successful calls in each population, in the
        order that the alleles are specified in 'segregating'.
        Non-diallelic polymorphisms are skipped.

        :param pop_ids: list of which populations to make fs for.
        :param projections: list of sample sizes to project down to for each
            population.
        :param polarized: If True, the data are assumed to be correctly polarized by
            'outgroup_allele'. SNPs in which the 'outgroup_allele'
            information is missing or '-' or not concordant with the
            segregating alleles will be ignored.
            If False, any 'outgroup_allele' info present is ignored,
            and the returned spectrum is folded.
        """
        Npops = len(pop_ids)
        fs = numpy.zeros(numpy.asarray(projections) + 1)
        for snp, snp_info in data_dict.items():
            # Skip SNPs that aren't biallelic.
            if len(snp_info["segregating"]) != 2:
                continue

            allele1, allele2 = snp_info["segregating"]
            if not polarized:
                # If we don't want to polarize, we can choose which allele is
                # derived arbitrarily since we'll fold anyways.
                outgroup_allele = allele1
            elif (
                "outgroup_allele" in snp_info
                and snp_info["outgroup_allele"] != "-"
                and snp_info["outgroup_allele"] in snp_info["segregating"]
            ):
                # Otherwise we need to check that it's a useful outgroup
                outgroup_allele = snp_info["outgroup_allele"]
            else:
                # If we're polarized and we didn't have good outgroup info, skip
                # this SNP.
                continue

            # Extract the allele calls for each population.
            allele1_calls = [snp_info["calls"][pop][0] for pop in pop_ids]
            allele2_calls = [snp_info["calls"][pop][1] for pop in pop_ids]
            # How many chromosomes did we call successfully in each population?
            successful_calls = [
                a1 + a2 for (a1, a2) in zip(allele1_calls, allele2_calls)
            ]

            # Which allele is derived (different from outgroup)?
            if allele1 == outgroup_allele:
                derived_calls = allele2_calls
            elif allele2 == outgroup_allele:
                derived_calls = allele1_calls

            # To handle arbitrary numbers of populations in the fs, we need
            # to do some tricky slicing.
            slices = [[numpy.newaxis] * len(pop_ids) for ii in range(Npops)]
            for ii in range(len(pop_ids)):
                slices[ii][ii] = slice(None, None, None)

            # Do the projection for this SNP.
            pop_contribs = []
            iter = zip(projections, successful_calls, derived_calls)
            for pop_ii, (p_to, p_from, hits) in enumerate(iter):
                contrib = Numerics._cached_projection(p_to, p_from, hits)[
                    slices[pop_ii]
                ]
                pop_contribs.append(contrib)
            fs += functools.reduce(operator.mul, pop_contribs)
        fsout = Spectrum(fs, mask_corners=mask_corners, pop_ids=pop_ids)
        if polarized:
            return fsout
        else:
            return fsout.fold()

    @staticmethod
    def _from_count_dict(count_dict, projections, polarized=True, pop_ids=None):
        """
        Frequency spectrum from data mapping SNP configurations to counts.

        :param count_dict: Result of Misc.count_data_dict
        :param projections: List of sample sizes to project down to for each
            population.
        :param polarized: If True, only include SNPs that count_dict marks as polarized.
            If False, include all SNPs and fold resulting Spectrum.
        :param pop_ids: Optional list of strings containing the population labels.
        """

        # create slices for projection calculation
        slices = [[numpy.newaxis] * len(projections) for ii in range(len(projections))]
        for ii in range(len(projections)):
            slices[ii][ii] = slice(None, None, None)

        fs_total = moments.Spectrum(
            numpy.zeros(numpy.array(projections) + 1), pop_ids=pop_ids
        )
        for (
            (called_by_pop, derived_by_pop, this_snp_polarized),
            count,
        ) in count_dict.items():
            if polarized and not this_snp_polarized:
                continue
            pop_contribs = []
            iter = zip(projections, called_by_pop, derived_by_pop)
            for pop_ii, (p_to, p_from, hits) in enumerate(iter):
                contrib = Numerics._cached_projection(p_to, p_from, hits)[
                    slices[pop_ii]
                ]
                pop_contribs.append(contrib)
            fs_proj = functools.reduce(operator.mul, pop_contribs)

            # create slices for adding projected fs to overall fs
            fs_total += count * fs_proj
        if polarized:
            return fs_total
        else:
            return fs_total.fold()

    @staticmethod
    def _data_by_tri(data_dict):
        """
        Nest the data by derived context and outgroup base.

        The resulting dictionary contains only SNPs which are appropriate for
        use of Hernandez's ancestral misidentification correction. It is
        organized as {(derived_tri, outgroup_base): {snp_id: data,...}}
        """
        result = {}
        genetic_bases = "ACTG"
        for snp, snp_info in data_dict.items():
            # Skip non-diallelic polymorphisms
            if len(snp_info["segregating"]) != 2:
                continue
            allele1, allele2 = snp_info["segregating"]
            # Filter out SNPs where we either non-constant ingroup or outgroup
            # context.
            try:
                ingroup_tri = snp_info["context"]
                outgroup_tri = snp_info["outgroup_context"]
            except KeyError:
                continue
            if not outgroup_tri[1] == snp_info["outgroup_allele"]:
                raise ValueError(
                    "Outgroup context and allele are inconsistent "
                    "for polymorphism: %s." % snp
                )
            outgroup_allele = outgroup_tri[1]

            # These are all the requirements to apply the ancestral correction.
            # First 2 are constant context.
            # Next 2 are sensible context.
            # Next 1 is that outgroup allele is one of the segregating.
            # Next 2 are that segregating alleles are sensible.
            if (
                outgroup_tri[0] != ingroup_tri[0]
                or outgroup_tri[2] != ingroup_tri[2]
                or ingroup_tri[0] not in genetic_bases
                or ingroup_tri[2] not in genetic_bases
                or outgroup_allele not in [allele1, allele2]
                or allele1 not in genetic_bases
                or allele2 not in genetic_bases
            ):
                continue

            if allele1 == outgroup_allele:
                derived_allele = allele2
            elif allele2 == outgroup_allele:
                # In this case, the second allele is non_outgroup
                derived_allele = allele1
            derived_tri = ingroup_tri[0] + derived_allele + ingroup_tri[2]
            result.setdefault((derived_tri, outgroup_allele), {})
            result[derived_tri, outgroup_allele][snp] = snp_info
        return result

    @staticmethod
    def from_data_dict_corrected(
        data_dict, pop_ids, projections, fux_filename, force_pos=True, mask_corners=True
    ):
        """
        Spectrum from a dictionary of polymorphisms, corrected for ancestral
        misidentification.

        The correction is based upon Hernandez, Williamson & Bustamante _Mol_Biol_Evol_
            24:1792 (2007)

        :param force_pos: If the correction is too agressive, it may leave some small
            entries in the fs less than zero. If force_pos is true,
            these entries will be set to zero, in such a way that the
            total number of segregating SNPs is conserved.
        :param fux_filename: The name of the file containing the
            misidentification probabilities.
            The file is of the form:
            - # Any number of comments lines beginning with #
            - AAA T 0.001
            - AAA G 0.02
            - ...
            Where every combination of three + one bases is considered
            (order is not important).  The triplet is the context and
            putatively derived allele (x) in the reference species. The
            single base is the base (u) in the outgroup. The numerical
            value is 1-f_{ux} in the notation of the paper.

        The data dictionary should be organized as:
            {snp_id:{'segregating': ['A','T'],
                     'calls': {'YRI': (23,3),
                                'CEU': (7,3)
                                },
                     'outgroup_allele': 'T',
                     'context': 'CAT',
                     'outgroup_context': 'CAT'
                    }
            }
        The additional entries are 'context', which includes the two flanking
        bases in the species of interest, and 'outgroup_context', which
        includes the aligned bases in the outgroup.

        This method skips entries for which the correction cannot be applied.
        Most commonly this is because of missing or non-constant context.
        """
        # Read the fux file into a dictionary.
        fux_dict = {}
        f = open(fux_filename)
        for line in f.readlines():
            if line.startswith("#"):
                continue
            sp = line.split()
            fux_dict[(sp[0], sp[1])] = 1 - float(sp[2])
        f.close()

        # Divide the data into classes based on ('context', 'outgroup_allele')
        by_context = Spectrum._data_by_tri(data_dict)

        fs = numpy.zeros(numpy.asarray(projections) + 1)
        while by_context:
            # Each time through this loop, we eliminate two entries from the
            # data dictionary. These correspond to one class and its
            # corresponding misidentified class.
            (derived_tri, out_base), nomis_data = by_context.popitem()

            # The corresponding bases if the ancestral state had been
            # misidentifed.
            mis_out_base = derived_tri[1]
            mis_derived_tri = derived_tri[0] + out_base + derived_tri[2]
            # Get the data for that case. Note that we default to an empty
            # dictionary if we don't have data for that class.
            mis_data = by_context.pop((mis_derived_tri, mis_out_base), {})

            fux = fux_dict[(derived_tri, out_base)]
            fxu = fux_dict[(mis_derived_tri, mis_out_base)]

            # Get the spectra for these two cases
            Nux = Spectrum.from_data_dict(nomis_data, pop_ids, projections)
            Nxu = Spectrum.from_data_dict(mis_data, pop_ids, projections)

            # Equations 5 & 6 from the paper.
            Nxu_rev = Numerics.reverse_array(Nxu)
            Rux = (fxu * Nux - (1 - fxu) * Nxu_rev) / (fux + fxu - 1)
            Rxu = Numerics.reverse_array(
                (fux * Nxu_rev - (1 - fux) * Nux) / (fux + fxu - 1)
            )

            fs += Rux + Rxu

        # Here we take the negative entries, and flip them back, so they end up
        # zero and the total number of SNPs is conserved.
        if force_pos:
            negative_entries = numpy.minimum(0, fs)
            fs -= negative_entries
            fs += Numerics.reverse_array(negative_entries)

        return Spectrum(fs, mask_corners=mask_corners, pop_ids=pop_ids)

    @staticmethod
    def from_demes(
        g, sampled_demes, sample_sizes, sample_times=None, Ne=None, unsampled_Ne=4
    ):
        """
        Takes a deme graph and computes the SFS. ``demes`` is a package for
        specifying demographic models in a user-friendly, human-readable YAML
        format. This function automatically parses the demographic description
        and returns a SFS for the specified populations and sample sizes.

        This function is new in version 1.1.0. Future developments will allow for
        inference using ``demes``-based demographic descriptions.

        :param g: A ``demes`` DemeGraph from which to compute the SFS. The DemeGraph
            can either be specified as a YAML file, in which case `g` is a string,
            or as a ``DemeGraph`` object.
        :type g: str or :class:`demes.DemeGraph`
        :param sampled_demes: A list of deme IDs to take samples from. We can repeat
            demes, as long as the sampling of repeated deme IDs occurs at distinct
            times.
        :type sampled_demes: list of strings
        :param sample_sizes: A list of the same length as ``sampled_demes``,
            giving the sample sizes for each sampled deme.
        :type sample_sizes: list of ints
        :param sample_times: If None, assumes all sampling occurs at the end of the
            existence of the sampled deme. If there are
            ancient samples, ``sample_times`` must be a list of same length as
            ``sampled_demes``, giving the sampling times for each sampled
            deme. Sampling times are given in time units of the original deme graph,
            so might not necessarily be generations (e.g. if ``g.time_units`` is years)
        :type sapmle_times: list of floats, optional
        :param Ne: reference population size. If none is given, we use the initial
            size of the root deme.
        :type Ne: float, optional
        :param unsampled_n: The default sample size of unsampled demes, which must be
            greater than or equal to 4.
        :type unsampled_n: int, optional
        :return: A ``moments`` site frequency spectrum, with dimension equal to the
            length of ``sampled_demes``, and shape equal to ``sample_sizes`` plus one
            in each dimension, indexing the allele frequency in each deme from 0
            to n[i], where i is the deme index.
        :rtype: :class:`moments.Spectrum`
        """
        raise ValueError(
            "The method ``from_demes`` is not yet implemented - coming ver 1.1"
        )

    # The code below ensures that when I do arithmetic with Spectrum objects,
    # it is not done between a folded and an unfolded array. If it is, I raise
    # a ValueError.

    # While I'm at it, I'm also fixing the annoying behavior that if a1 and a2
    # are masked arrays, and a3 = a1 + a2. Then wherever a1 or a2 was masked,
    # a3.data ends up with the a1.data values, rather than a1.data + a2.data.
    # Note that this fix doesn't work for operation by numpy.ma.exp and
    # numpy.ma.log. Guess I can't have everything.

    # I'm using exec here to avoid copy-pasting a dozen boiler-plate functions.
    # The calls to check_folding_equal ensure that we don't try to combine
    # folded and unfolded Spectrum objects.

    # I set check_folding = False in the constructor because it raises useless
    # warnings when, for example, I do (model + 1).

    # These functions also ensure that the pop_ids
    # get properly copied over.

    # This is pretty advanced Python voodoo, so don't fret if you don't
    # understand it at first glance. :-)
    for method in [
        "__add__",
        "__radd__",
        "__sub__",
        "__rsub__",
        "__mul__",
        "__rmul__",
        "__div__",
        "__rdiv__",
        "__truediv__",
        "__rtruediv__",
        "__floordiv__",
        "__rfloordiv__",
        "__rpow__",
        "__pow__",
    ]:
        exec(
            """
def %(method)s(self, other):
    self._check_other_folding(other)
    if isinstance(other, numpy.ma.masked_array):
        newdata = self.data.%(method)s (other.data)
        newmask = numpy.ma.mask_or(self.mask, other.mask)
    else:
        newdata = self.data.%(method)s (other)
        newmask = self.mask
    newpop_ids = self.pop_ids
    if hasattr(other, 'pop_ids'):
        if other.pop_ids is None:
            newpop_ids = self.pop_ids
        elif self.pop_ids is None:
            newpop_ids = other.pop_ids
        elif other.pop_ids != self.pop_ids:
            logger.warn('Arithmetic between Spectra with different pop_ids. '
                        'Resulting pop_id may not be correct.')
    outfs = self.__class__.__new__(self.__class__, newdata, newmask, 
                                   mask_corners=False, data_folded=self.folded,
                                   check_folding=False, pop_ids=newpop_ids)
    return outfs
"""
            % {"method": method}
        )

    # Methods that modify the Spectrum in-place.
    for method in [
        "__iadd__",
        "__isub__",
        "__imul__",
        "__idiv__",
        "__itruediv__",
        "__ifloordiv__",
        "__ipow__",
    ]:
        exec(
            """
def %(method)s(self, other):
    self._check_other_folding(other)
    if isinstance(other, numpy.ma.masked_array):
        self.data.%(method)s (other.data)
        self.mask = numpy.ma.mask_or(self.mask, other.mask)
    else:
        self.data.%(method)s (other)
    if hasattr(other, 'pop_ids') and other.pop_ids is not None\
             and other.pop_ids != self.pop_ids:
        logger.warn('Arithmetic between Spectra with different pop_ids. '
                    'Resulting pop_id may not be correct.')
    return self
"""
            % {"method": method}
        )

    def _check_other_folding(self, other):
        """
        Ensure other Spectrum has same .folded status
        """
        if isinstance(other, self.__class__) and other.folded != self.folded:
            raise ValueError(
                "Cannot operate with a folded Spectrum and an " "unfolded one."
            )


# Allow spectrum objects to be pickled.
# See http://effbot.org/librarybook/copy-reg.htm
try:
    import copy_reg

    def Spectrum_unpickler(data, mask, data_folded, pop_ids):
        return moments.Spectrum(
            data,
            mask,
            mask_corners=False,
            data_folded=data_folded,
            check_folding=False,
            pop_ids=pop_ids,
        )

    def Spectrum_pickler(fs):
        return (
            Spectrum_unpickler,
            (fs.data, fs.mask, fs.folded, fs.pop_ids),
        )

    copy_reg.pickle(Spectrum, Spectrum_pickler, Spectrum_unpickler)
except:
    import copyreg

    def Spectrum_unpickler(data, mask, data_folded, pop_ids):
        return moments.Spectrum(
            data,
            mask,
            mask_corners=False,
            data_folded=data_folded,
            check_folding=False,
            pop_ids=pop_ids,
        )

    def Spectrum_pickler(fs):
        return (
            Spectrum_unpickler,
            (fs.data, fs.mask, fs.folded, fs.pop_ids),
        )

    copyreg.pickle(Spectrum, Spectrum_pickler, Spectrum_unpickler)
