"""
Event cache — reuse event sets for identical local environments.

Central to making the dynamic KMC practical: "reuse, don't recompute."
When a local environment has been seen before, its entire event set
(with rates) is retrieved from cache instead of regenerated.
"""

from collections import defaultdict


class CachedEventSet:
    """
    A cached set of events for a specific local environment.

    Stores event templates (without site-specific indices) and their rates.
    When reused, the templates are instantiated with the actual site index.
    """

    __slots__ = ('env', 'event_templates', 'hit_count')

    def __init__(self, env, event_templates):
        self.env = env
        self.event_templates = event_templates  # list of (event_class, kwargs, rate)
        self.hit_count = 0


class EventCache:
    """
    Cache for environment → event set mapping.

    When the local environment at a site matches a previously seen
    environment (same hash), the cached event set is reused directly.
    This avoids regenerating events and re-querying the rate estimator.

    Parameters
    ----------
    max_size : int
        Maximum number of cached environments. 0 = unlimited.
    debug : bool
    """

    def __init__(self, max_size=0, debug=False):
        self.max_size = max_size
        self.debug = debug
        self._cache = {}  # env_hash → CachedEventSet
        self._hits = 0
        self._misses = 0

    @property
    def size(self):
        return len(self._cache)

    @property
    def hit_rate(self):
        total = self._hits + self._misses
        return self._hits / total if total > 0 else 0.0

    def lookup(self, env):
        """
        Look up a cached event set for an environment.

        Parameters
        ----------
        env : LocalEnvironment

        Returns
        -------
        CachedEventSet or None
        """
        h = env.env_hash()
        if h in self._cache:
            cached = self._cache[h]
            # Double-check full equality (hash collision guard)
            if cached.env == env:
                cached.hit_count += 1
                self._hits += 1
                return cached
        self._misses += 1
        return None

    def store(self, env, event_templates):
        """
        Store an event set in the cache.

        Parameters
        ----------
        env : LocalEnvironment
        event_templates : list of (event_class, kwargs, rate)
        """
        h = env.env_hash()

        # Evict if at capacity (LRU-like: remove least hit)
        if self.max_size > 0 and len(self._cache) >= self.max_size:
            min_key = min(self._cache, key=lambda k: self._cache[k].hit_count)
            del self._cache[min_key]

        self._cache[h] = CachedEventSet(env, event_templates)

        if self.debug:
            print(f"CACHE: Stored env {env}, cache size={self.size}")

    def invalidate(self, env):
        """Remove a specific environment from the cache."""
        h = env.env_hash()
        if h in self._cache:
            del self._cache[h]

    def clear(self):
        """Clear entire cache."""
        self._cache.clear()
        self._hits = 0
        self._misses = 0

    def summary(self):
        """Print cache statistics."""
        print(f"EventCache: {self.size} entries, "
              f"hit_rate={self.hit_rate:.1%} "
              f"({self._hits} hits / {self._misses} misses)")
        if self.debug:
            for h, cs in self._cache.items():
                print(f"  {cs.env}: {len(cs.event_templates)} events, "
                      f"hits={cs.hit_count}")

    def __repr__(self):
        return (f"EventCache(size={self.size}, "
                f"hit_rate={self.hit_rate:.1%})")
