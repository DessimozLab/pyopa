import bisect
from . import backend


class Aligner(object):
    def __init__(self, env=None):
        if env is None:
            envs = backend.pyopa.load_default_environments()
            self.environments = envs['environments']
            self.pam1 = envs['log_pam1']
        else:
            raise NotImplementedError("non-default environments are not yet implemented")
        self._pam_distances = [env.pam for env in self.environments]

    def environment_at_distance(self, dist):
        """returns the :class:`AlignmentEnvironment` for a given distance.

        The returned AlignmentEnvironment is optimized for the given
        evolutionary distance, i.e. sequences that diverged for this
        evolutionary distance will reach the highest alignment score
        with this scoring matrix.

        :param numeric dist: the evolutionary distance"""
        index = bisect.bisect_right(self._pam_distances, dist)

        if index == len(self._pam_distances):
            return self.environments[index - 1]
        return min(self.environments[index - 1], self.environments[index],
                   key=lambda x: abs(x.pam - dist))

    def score(self, s1, s2, env=None, **kwargs):
        """Compute the alignment score using Farrer's algorithm.

        Compute the optimal pairwise alignment score of the two sequences
        `s1` and `s2` for either a fixed scoring environment, or compute
        the overall best score
        """
        pass

    def align(self, s1, s2, env=None, **kwargs):
        pass

    def estimate_distance(self, s1, s2):
        """Estimates the evolutionary distance for a fixed alignment.

        This method requires already aligned sequences `s1` and `s2`. """
        pass