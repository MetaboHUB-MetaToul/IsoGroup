from __future__ import annotations
import numpy as np
from typing import List, Iterator
from isogroup.base.feature import Feature

class Cluster:
    """
    Represents a cluster of mass spectrometry features.
    A cluster is a group of mass features originating from the same molecule, sharing the same elemental composition but different isotopic compositions.
    Clusters are used to group features related to the same metabolite or chemical compound.

    """

    def __init__(self, features: list, cluster_id :str, name :str=None):
        """ 
        :param features: List of features in the cluster.
        :param cluster_id: Unique identifier for the cluster.
        :param name: Name of the cluster, usually corresponding to the annotated name of the metabolite or compound.
        
        """
        self.features = features
        self.cluster_id = cluster_id
        self.tracer_element = features[0]._tracer_element if features is not None else None
        self.tracer = features[0].tracer if features is not None else None
        self.name = name
        self._formula = None

    def __repr__(self) -> str:
        return f"Cluster({self.cluster_id}, {self.features})"
    
    def __len__(self) -> int:
        """
        Returns the number of features in the cluster
        """
        return len(self.features)

    def __iter__(self) -> Iterator[Feature]:
        """
        Return a list iterator of the features in the cluster
        """
        return iter(self.features)
    
    # def __contains__(self, item):
    #     pass

    
    @property
    def lowest_rt(self) -> float:
        """
        Returns the lowest retention time (RT) of the features in the cluster.
        """
        return min([f.rt for f in self.features])

    @property
    def highest_rt(self) -> float:
        """
        Returns the highest retention time (RT) of the features in the cluster.
        """
        return max([f.rt for f in self.features])

    @property
    def lowest_mz(self) -> float:
        """
        Returns the lowest mass-to-charge ratio (m/z) of the features in the cluster.
        """
        return min([f.mz for f in self.features])

    @property
    def highest_mz(self) -> float:
        """
        Returns the highest mass-to-charge ratio (m/z) of the features in the cluster.
        """
        return max([f.mz for f in self.features])
    
    @property
    def mean_rt(self) -> float:
        """
        Returns the mean retention time (RT) of the features in the cluster.
        """
        return np.mean([f.rt for f in self.features])
    
    @property
    def mean_mz(self) -> float:
        """
        Returns the mean mass-to-charge ratio (m/z) of the features in the cluster.
        """
        return np.mean([f.mz for f in self.features])
    
    @property
    def metabolite(self):
        """
        Returns the list of metabolite annotations for features in the cluster.
        """
        return [f.metabolite for f in self.features]
    
    @property
    def chemical(self):
        """
        Returns the list of chemical objects for the cluster.
        Based on the metabolite name matching to the cluster name.
        """
        for feature in self.features:
            if self.name in feature.metabolite:
                idx = feature.metabolite.index(self.name)
                return feature.chemical[idx]

    @property
    def element_number(self) -> int:
        """
        Returns the number of tracer elements in the cluster.
        """
        self.formula
        return self._formula[self.tracer_element]

    
    @property
    def isotopologues(self) -> List[int]:
        """
        Returns the list of isotopologues in the cluster.
        Based on the metabolite name matching to the cluster name.
        """
        isotopologues = []
        for feature in self.features:
            if self.name in feature.metabolite:
                isotopologues.append(feature.cluster_isotopologue[self.name])
                # idx = feature.metabolite.index(self.name)
                # isotopologues.append(feature.isotopologue[idx])
        return isotopologues


    @property
    def formula(self) -> str:
        """
        Returns the formula of the cluster.
        Based on the metabolite name matching to the cluster name.
        """
        if self._formula is None:
            # Check if the cluster is annotated
            if self.name is None:
                raise ValueError("No cluster found. Run clusterize() first")
            
            if self.name in self.features[0].metabolite:
                idx = self.features[0].metabolite.index(self.name)
                self._formula = self.features[0].formula[idx]
            else:  
                raise ValueError(f"No metabolite matching '{self.name}' found in the cluster features")
            
        return self._formula


    @property
    def expected_isotopologues_in_cluster(self) -> List[int]:
        """
        Returns the list of expected isotopologues in the cluster.
        Based on the number of tracer element in its formula.
        """
        return list(range(self.element_number + 1))
                           

    @property
    def is_complete(self) -> bool:
        """
        Returns True if the cluster is complete (i.e contains all isotopologues expected).
        """   
        return len(self) == self.element_number + 1 and self.expected_isotopologues_in_cluster == self.isotopologues
    
    @property
    def is_incomplete(self) -> bool:
        """
        Returns True if the cluster is incomplete (i.e contains less isotopologues than expected).
        """
        return len(self) < self.element_number + 1 or len(set(self.isotopologues)) != len(self.expected_isotopologues_in_cluster)

    @property
    def is_duplicated(self) -> bool:
        """
        Returns True if the cluster contains duplicated isotopologues.
        """
        return len(set(self.isotopologues)) != len(self.isotopologues)

    @property
    def is_corrupted(self) -> bool:
        """
        Returns True if the cluster is corrupted (overfilled ?) (i.e contains isotopologues not expected)
        """ 
        pass


    @property 
    def status(self) -> str:
        """
        Returns the status of the cluster based on its completeness, incompleteness, and duplication.
        """
        status = []

        if self.is_complete:
            status.append("Complete")
        
        # if self.is_incomplete and self.is_duplicated:
        #     return "Incomplete, duplicated isotopologues"

        if self.is_incomplete:
            status.append("Incomplete")

        if self.is_duplicated:
            status.append("Duplicated isotopologues")
        
        # if self.is_corrupted:
        #     return "Corrupted"
        return str(", ".join(status))

    @property
    def missing_isotopologues(self) -> List[int]:
        """
        Returns a list of missing isotopologues in the annotated cluster.
        Based on the expected isotopologues in the cluster.
        """
        if self.is_incomplete:
            return [i for i in self.expected_isotopologues_in_cluster if i not in self.isotopologues]
        

    @property
    def duplicated_isotopologues(self) -> List[int]:
        """
        Returns a list of duplicated isotopologues in the cluster.
        """
        if self.is_duplicated:
            return [i for i in set(self.isotopologues) if self.isotopologues.count(i) > 1]
        
    
    @property
    def is_adduct(self) -> tuple[bool, str]:
        pass

    @property
    def summary(self) -> dict:
        """
        Returns a summary of the cluster.

        """
        return {
            "ClusterID": self.cluster_id,
            "Name": self.name,
            "Number_of_features": len(self),
            "Isotopologues": self.isotopologues,
            "Status": self.status,
            "Missing_isotopologues": self.missing_isotopologues,
            "Duplicated_isotopologues": self.duplicated_isotopologues
        }
    

    # def __add__(self, other: Union[Feature, Self]) -> Self:
    #     """
    #     Extend the cluster with either an extra feature or another cluster.
    #     :param other: Feature or Cluster
    #     :return: self
    #     """
    #     pass

    # def __iadd__(self, other: Union[Feature, Self]) -> Self:
    #     """
    #     Extend the cluster with either an extra feature or another cluster.
    #     :param other: Feature or Cluster
    #     :return: self
    #     """
    #     pass

    # def __remove__(self, other: Union[Feature, Self]) -> Self:
    #     """
    #     Remove a feature from the cluster
    #     :param other: Feature
    #     :return: self
    #     """
    #     pass
