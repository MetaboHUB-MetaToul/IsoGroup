from typing import List, Union, Iterator, Self

from isogroup.base.feature import Feature


class Cluster:

    def __init__(self, features: List[Feature]):
        self.features = features

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
        # Returns the metabolite of the annotated features in the cluster
        return self.features[0].metabolite
    
    @property
    def isotopologues(self) -> List[int]:
        # Returns the isotopologues of the annotated features in the cluster
        return [f.isotopologue for f in self.features]

    @property
    def missing_isopologue(self) -> List[int]:
        # Returns a list of missing isotopologues in the cluster
        pass

    @property
    def is_adduct(self) -> tuple[bool, str]:
        # Retourne True si le cluster est un adduct from another cluster
        pass
