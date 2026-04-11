"""
Local environment descriptor and hashing.

Encodes the local chemical environment around a site into a hashable
descriptor. Two sites with the same descriptor have the same local
environment and therefore the same event set and rates.

Descriptor components:
  - center atom type
  - sorted nearest-neighbor composition tuple
  - adsorbate state at center
  - adsorbate count among neighbors
  - coordination number
  - (optional) motif class, strain
"""

from collections import defaultdict
from hashlib import sha256
import struct


class LocalEnvironment:
    """
    Descriptor of the local chemical environment around a site.

    This is the key object for environment caching: two sites with
    identical LocalEnvironments share the same event set and rates.

    Parameters
    ----------
    center_type : str
        Atom type at the center site.
    center_adsorbate : int
        Adsorbate species ID at center (0=empty).
    nn_composition : tuple of (str, int)
        Sorted tuple of (atom_type, count) for nearest neighbors.
        Example: (('Au', 2), ('Pd', 4))
    nn_adsorbate_count : int
        Number of occupied nearest-neighbor sites.
    coordination : int
        Total number of nearest neighbors.
    layer : int
        Layer index of the center site.
    """

    __slots__ = ('center_type', 'center_adsorbate', 'nn_composition',
                 'nn_adsorbate_count', 'coordination', 'layer', '_hash')

    def __init__(self, center_type, center_adsorbate, nn_composition,
                 nn_adsorbate_count, coordination, layer=0):
        self.center_type = center_type
        self.center_adsorbate = center_adsorbate
        self.nn_composition = tuple(sorted(nn_composition))
        self.nn_adsorbate_count = nn_adsorbate_count
        self.coordination = coordination
        self.layer = layer
        self._hash = None

    def env_hash(self):
        """
        Compute a deterministic hash for this environment.

        Two LocalEnvironments with the same hash are considered identical
        for event caching purposes.

        Returns
        -------
        int
            64-bit hash value.
        """
        if self._hash is not None:
            return self._hash

        key = (
            self.center_type,
            self.center_adsorbate,
            self.nn_composition,
            self.nn_adsorbate_count,
            self.coordination,
            self.layer,
        )
        h = sha256(repr(key).encode()).digest()
        self._hash = struct.unpack('<Q', h[:8])[0]
        return self._hash

    def __eq__(self, other):
        if not isinstance(other, LocalEnvironment):
            return False
        return (self.center_type == other.center_type
                and self.center_adsorbate == other.center_adsorbate
                and self.nn_composition == other.nn_composition
                and self.nn_adsorbate_count == other.nn_adsorbate_count
                and self.coordination == other.coordination
                and self.layer == other.layer)

    def __hash__(self):
        return self.env_hash()

    def __repr__(self):
        comp_str = ','.join(f"{t}:{n}" for t, n in self.nn_composition)
        ads = f"+sp{self.center_adsorbate}" if self.center_adsorbate else ""
        return (f"Env({self.center_type}{ads}|NN=[{comp_str}]"
                f"|ads_nn={self.nn_adsorbate_count}|CN={self.coordination})")


class EnvHash:
    """
    Utility class for building LocalEnvironments from a DynamicSurface.
    """

    @staticmethod
    def from_surface(surface, site_idx):
        """
        Extract the LocalEnvironment for a site on the surface.

        Parameters
        ----------
        surface : DynamicSurface
        site_idx : int

        Returns
        -------
        LocalEnvironment
        """
        site = surface.sites[site_idx]

        # Check cache
        if site._env_hash is not None:
            return site._env_hash

        # NN composition
        comp = defaultdict(int)
        nn_ads = 0
        for j in surface.neighbors[site_idx]:
            nn = surface.sites[j]
            comp[nn.atom_type] += 1
            if nn.is_occupied:
                nn_ads += 1

        env = LocalEnvironment(
            center_type=site.atom_type,
            center_adsorbate=site.adsorbate,
            nn_composition=tuple(comp.items()),
            nn_adsorbate_count=nn_ads,
            coordination=len(surface.neighbors[site_idx]),
            layer=site.layer,
        )

        # Cache on the site
        site._env_hash = env
        return env

    @staticmethod
    def from_surface_all(surface):
        """
        Extract LocalEnvironments for all sites.

        Returns
        -------
        list of LocalEnvironment
        """
        return [EnvHash.from_surface(surface, i)
                for i in range(surface.n_sites)]

    @staticmethod
    def count_unique(surface):
        """
        Count the number of unique environments on the surface.

        Returns
        -------
        dict
            {LocalEnvironment: count}
        """
        envs = EnvHash.from_surface_all(surface)
        counts = defaultdict(int)
        for env in envs:
            counts[env] += 1
        return dict(counts)
