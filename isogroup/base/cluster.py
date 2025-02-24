from typing import List, Union, Iterator, Self
from isogroup.base.feature import Feature


class Cluster:

    def __init__(self, features: List[Feature], cluster_id = None):
        self.features = features
        self.cluster_id = cluster_id

    def __repr__(self) -> str:
        return f"Cluster({self.features})"
    
    
    def __len__(self) -> int:
        return len(self.features)

    def __add__(self, other: Union[Feature, Self]) -> Self:
        """
        Extend the cluster with either an extra feature or another cluster.
        :param other: Feature or Cluster
        :return: self
        """
        pass

    def __iadd__(self, other: Union[Feature, Self]) -> Self:
        """
        Extend the cluster with either an extra feature or another cluster.
        :param other: Feature or Cluster
        :return: self
        """
        pass

    def __contains__(self, item):
        pass

    def __iter__(self) -> Iterator[Feature]:
        return iter(self.features)
    
    
    def get_cluster(metabolite: str, sample: str) -> Self:
        """
        Returns the cluster of the metabolite in the given sample
        """
        pass

    @property
    def lowest_rt(self) -> float:
        return min([f.rt for f in self.features])

    @property
    def highest_rt(self) -> float:
        return max([f.rt for f in self.features])

    @property
    def lowest_mz(self) -> float:
        return min([f.mz for f in self.features])

    @property
    def highest_mz(self) -> float:
        return max([f.mz for f in self.features])

    @property
    def metabolite(self) -> str:
        """
        Returns the metabolite of the annotated features in the cluster
        """
        return self.features[0].metabolite
    
    @property
    def isotopologues(self) -> List[int]:
        """
        Returns the isotopologues of the annotated features in the cluster
        """
        return [f.isotopologue for f in self.features]
 
    @property
    def is_complete(self) -> bool:
        """
        Returns True if the cluster is complete
        """
        delta_mz_tracer = 1.00335 # Should be added to the attributes of the class for untargeted analysis ?
        # Suppose that the lowest mz is M0 and the highest mz is Mn
        # Check only for annotated clusters ?
        excepted_length = int(round((self.highest_mz - self.lowest_mz) / delta_mz_tracer, 1)) + 1
        return len(self.features) == excepted_length

    @property
    def missing_isopologue(self) -> List[int]:
        """
        Returns a list of missing isotopologues in the cluster
        """
        pass    
    
    @property
    def is_adduct(self) -> tuple[bool, str]:
        pass
